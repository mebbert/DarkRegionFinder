/**
 * 
 */
package ebbert.cgf;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset;


/**
 * @author markebbert
 *
 */
public class CamoGeneFinder {
	
	private static Logger logger = Logger.getLogger(CamoGeneFinder.class);
	private static int MAPQ_THRESHOLD, MIN_CAMO_DEPTH, DARK_DEPTH,
						DARK_LOWER = 1, DARK_UPPER = 9, MIN_REGION_SIZE,
						MIN_MAPQ_MASS;
	private static ValidationStringency SAM_VALIDATION_STRINGENCY;
	
	BufferedWriter camoWriter, darkWriter, incWriter;

	private SAMFileHeader header;
	private SamReader reader;
	
	private IndexedFastaSequenceFile hgRefReader;
	

	/**
	 * Instantiate a CamoGeneFinder object.
	 * 
	 * @param samFile
	 * @param outCamoBed
	 * @param hgRef
	 * @param mapqThreshold
	 * @param minMapQMass
	 * @param minRegionSize
	 * @param minCamoDepth
	 * @param vs
	 * @throws IOException
	 */
	public CamoGeneFinder(final File samFile, final File outCamoBed,
			File outDarkBed, File outIncBed,
			File hgRef, final int mapQThreshold, final int minMapQMass,
			final int minRegionSize, final int minCamoDepth, int darkDepth,
			final ValidationStringency vs/*, int startWalking,
			int endWalking*/) throws IOException {
		
		// Would be nice to be able to specify start/end locations
//		this.startWalking = startWalking;
//		this.endWalking = endWalking;
		
		camoWriter = new BufferedWriter(new OutputStreamWriter(
	              new FileOutputStream(outCamoBed), "utf-8"));
		camoWriter.write("chrom\tstart\tend\tnMapQBelowThreshold\tdepth\tpercMapQ0\n");

		darkWriter = new BufferedWriter(new OutputStreamWriter(
	              new FileOutputStream(outDarkBed), "utf-8"));
		darkWriter.write("chrom\tstart\tend\tdepth\n");

		incWriter = new BufferedWriter(new OutputStreamWriter(
	              new FileOutputStream(outIncBed), "utf-8"));
		incWriter.write("chrom\tstart\tend\n");

		CamoGeneFinder.MAPQ_THRESHOLD = mapQThreshold;
		CamoGeneFinder.MIN_REGION_SIZE = minRegionSize;
		CamoGeneFinder.MIN_CAMO_DEPTH = minCamoDepth;
		CamoGeneFinder.DARK_DEPTH = darkDepth;
		CamoGeneFinder.MIN_MAPQ_MASS = minMapQMass;
		CamoGeneFinder.SAM_VALIDATION_STRINGENCY = vs;
		
		this.hgRefReader = new IndexedFastaSequenceFile(hgRef);

		this.reader = CamoGeneFinder.openSam(samFile,
				CamoGeneFinder.SAM_VALIDATION_STRINGENCY);
		this.header = reader.getFileHeader();
		
		/* Get sample name(s) from the sam/bam file */
        TreeSet<String> samples = new TreeSet<String>();
        for(SAMReadGroupRecord group : header.getReadGroups()){
        	samples.add(group.getSample());
        }
	}


	/**
	 * @throws Exception 
	 */
	public void startWalkingByLocus() throws Exception{

		SamLocusIterator sli = new SamLocusIterator(reader);
		
		/* Set the base quality score cutoff to 0 
		 * so the iterator won't try to validate base qualities. Any
		 * secondary alignment won't have the qualities, and they're
		 * all made up anyway.
		 */
		int baseQualityScoreCutoff = 0;
		sli.setQualityScoreCutoff(baseQualityScoreCutoff);

		
		
		/* Walk along genome identifying 'dark' and 'camouflaged' regions */

		LocusInfo locus;
		int consecCamo = 0, consecDark = 0, consecInc = 0, pos, depth,
				nMapQBelowThreshold, nMapQBetween1And9;
		ArrayList<String> camoRegion = new ArrayList<String>(),
				darkRegion = new ArrayList<String>(),
				incRegion = new ArrayList<String>(),
				ignore = new ArrayList<String>();
		String contig; byte[] bases; byte base;
		double percMapQ0, percMapQBetween1And9;
		while(sli.hasNext()){

			locus = sli.next();
			
			contig = locus.getSequenceName();
			
			/* If this contig is not in the ref, then skip */
			if(ignore.contains(contig)) {
				continue;
			}

			/* Returns 1-based position */
			pos = locus.getPosition();  
			
			/* Ensure sequence is present in provided reference */
			if(hgRefReader.getSequenceDictionary().getSequence(contig) == null){
				logger.warn("BAM file contains alignments for " + contig
						+ " but this sequence was not found in the provided"
						+ " reference. Skipping.");
				ignore.add(contig);
				continue;
			}

			/* Expects 1-based position (inclusive to inclusive) */
			bases = hgRefReader.getSubsequenceAt(contig, pos, pos).getBases();

			/* Track progress */ 
        	if(pos % 1000000 == 0){
        		logger.debug("Assessed " + pos + " loci on " + contig);
        	}		

			
			/* Sanity check that we're only getting one base back */
			if(bases.length > 1) {
				logger.error("Expected only one base, but received " + bases.length);
				camoWriter.close();
				darkWriter.close();
				incWriter.close();
				sli.close();
				throw new Exception("Expected only one base, but received " + bases.length);
			}

			base = bases[0];

			/* Record incomplete genomic regions (i.e., 'N') */
			if(base == 'N' || base == 'n'){
				incRegion.add(incompleteRegionToString(contig, pos));
				consecInc++;

				/* Write dark and camo regions if large enough */
				if(consecDark >= CamoGeneFinder.MIN_REGION_SIZE){
					writeRegion(darkRegion, darkWriter);
				}
				if(consecCamo >= CamoGeneFinder.MIN_REGION_SIZE){
					writeRegion(camoRegion, camoWriter);
				}

				/* Clear regardless (i.e., even if the region wasn't large enough) */
				camoRegion.clear();
				consecCamo = 0;
				darkRegion.clear();
				consecDark = 0;

				continue;
			}
	

			/* Write incomplete regions if large enough. Clear in either case. */
			if(consecInc >= CamoGeneFinder.MIN_REGION_SIZE) {
				writeRegion(incRegion, incWriter);
			}

			/* Clear regardless because we know we're outside an incomplete
			 * region
			 */
			incRegion.clear();
			consecInc = 0; 		

        	
			/* Get depth at this position. */
			depth = locus.getRecordAndPositions().size();
			nMapQBelowThreshold = 0;
			nMapQBetween1And9 = 0;

			
			/* Get number of reads with MAPQ ≤ threshold */
			List<RecordAndOffset> recs = locus.getRecordAndPositions();
			int mapq;
			for(RecordAndOffset rec : recs){
				mapq = rec.getRecord().getMappingQuality();
				if(mapq <= CamoGeneFinder.MAPQ_THRESHOLD){
					nMapQBelowThreshold++;
				}
				else if(mapq >= CamoGeneFinder.DARK_LOWER &&
						mapq <= CamoGeneFinder.DARK_UPPER) {
					/* Yes, I just used magic numbers */
					nMapQBetween1And9++;
				}
			}

			/* If depth ≥ minimum required depth to trust mass AND we're above
			 * required MAPQ mass
			 */
			percMapQ0 = depth > 0 ? nMapQBelowThreshold / depth * 100 : -1;
			if(depth >= CamoGeneFinder.MIN_CAMO_DEPTH &&
					percMapQ0 >= CamoGeneFinder.MIN_MAPQ_MASS){

				camoRegion.add(camoRegionToString(contig, pos,
						nMapQBelowThreshold, depth, percMapQ0));
				consecCamo++;
			}
			else if(consecCamo >= CamoGeneFinder.MIN_REGION_SIZE){
				writeRegion(camoRegion, camoWriter);
				camoRegion.clear();
				consecCamo = 0;
			}
			else {
				/* Also clear if we didn't reach the region size */
				camoRegion.clear();
				consecCamo = 0;
			}

			/* Check if we're in a dark region, regardless of whether we're
			 * in a camo region. A region is 'dark' if depth is < DARK_DEPTH
			 * OR if the percentage of reads with a 0 < MAPQ < 10 is ≥
			 * MIN_MAPQ_MASS
			 */
			percMapQBetween1And9 = depth > 0 ? nMapQBetween1And9 / depth * 100 : -1;
			if(depth <= CamoGeneFinder.DARK_DEPTH ||
					percMapQBetween1And9 >= CamoGeneFinder.MIN_MAPQ_MASS) {

				/* Save 'dark' regions with low coverage. */
				darkRegion.add(darkRegionToString(contig, pos, nMapQBetween1And9,
						depth, percMapQBetween1And9));
				consecDark++;
			}
			else if(consecDark >= CamoGeneFinder.MIN_REGION_SIZE){
				/* Write dark regions if large enough */
				writeRegion(darkRegion, darkWriter);
				darkRegion.clear();
				consecDark = 0;
			}
			else {
				/* Also clear if we didn't reach the region size */
				darkRegion.clear();
				consecDark = 0;
			}

		}
		
		/* Write regions if large enough */
		if(consecCamo >= CamoGeneFinder.MIN_REGION_SIZE) {
			writeRegion(camoRegion, camoWriter);
		}
		if(consecInc >= CamoGeneFinder.MIN_REGION_SIZE) {
			writeRegion(incRegion, incWriter);
		}
		if(consecDark >= CamoGeneFinder.MIN_REGION_SIZE){
			writeRegion(darkRegion, darkWriter);
		}

		camoWriter.close();
		darkWriter.close();
		incWriter.close();
		sli.close();
	}
	
	/**
	 * @param contigName
	 * @param position
	 * @return
	 */
	private String incompleteRegionToString(String contigName, int position) {

		/* Bed files are 0-based. locus.getPosition() returns 1-based. #Annoying */
		return contigName + "\t" + (position - 1) + "\t" + position + "\n";
	}

	/**
	 * @param seqName
	 * @param position
	 * @param depth
	 * @return
	 */
	private String darkRegionToString(String contigName, int position,
			int nMapQBetween1And9, int depth, double percentMapQBetween1And9) {

		/* Bed files are 0-based. locus.getPosition() returns 1-based. #Annoying */
		return contigName + "\t" + (position - 1) + "\t" + position +
				"\t" + nMapQBetween1And9 + "\t" + depth +
				"\t" + percentMapQBetween1And9 + "\n";
	}

	/**
	 * @param seqName
	 * @param position
	 * @param nMapQBelowThreshold
	 * @param depth
	 * @param percentMapQBelowThreshold
	 * @return
	 */
	private String camoRegionToString(String contigName, int position,
			int nMapQBelowThreshold, int depth,
			double percentMapQBelowThreshold) {

		/* Bed files are 0-based. locus.getPosition() returns 1-based. #Annoying */
		return contigName + "\t" + (position - 1) + "\t" + position +
				"\t" + nMapQBelowThreshold + "\t" + depth +
				"\t" + percentMapQBelowThreshold + "\n";
	}
	
	/**
	 * @param camoRegions
	 * @param writer
	 * @throws IOException
	 */
	private void writeRegion(ArrayList<String> camoRegions,
			BufferedWriter writer) throws IOException {
		for(String s : camoRegions){
			writer.write(s);
		}
	}
    
	/**
	 * 
	 * Open a SAM/BAM file for reading and return the SamReader obj
	 * 
	 * @param samFile
	 * @param vs
	 * @return SamReader
	 */
	 private static SamReader openSam(final File samFile, ValidationStringency vs) {
		 
		final SamReaderFactory factory =
				  SamReaderFactory.makeDefault()
					  .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS,
							  SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
					  .validationStringency(vs);

		final SamReader reader = factory.open(samFile);
		
		 return reader;
	 }
    
}
