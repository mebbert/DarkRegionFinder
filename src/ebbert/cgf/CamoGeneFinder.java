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
	private static int MAPQ_THRESHOLD, MIN_DEPTH, MIN_REGION_SIZE, MIN_MAPQ_0_MASS;
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
	 * @param minMapQ0Mass
	 * @param minRegionSize
	 * @param minDepth
	 * @param vs
	 * @throws IOException
	 */
	public CamoGeneFinder(final File samFile, final File outCamoBed,
			File outDarkBed, File outIncBed,
			File hgRef, final int mapQThreshold, final int minMapQ0Mass,
			final int minRegionSize, final int minDepth,
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
		CamoGeneFinder.MIN_DEPTH = minDepth;
		CamoGeneFinder.MIN_MAPQ_0_MASS = minMapQ0Mass;
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
				nMapQBelowThreshold;
		ArrayList<String> camoRegion = new ArrayList<String>(),
				darkRegion = new ArrayList<String>(),
				incRegion = new ArrayList<String>();
		String contig; byte[] bases; byte base;
		while(sli.hasNext()){

			locus = sli.next();
			
			contig = locus.getSequenceName();

			/* Returns 1-based position */
			pos = locus.getPosition();  
			
			/* Expects 1-based position (inclusive to inclusive) */
			bases = hgRefReader.getSubsequenceAt(contig, pos, pos).getBases();

        	if(pos % 1000000 == 0){
        		logger.debug("Assessed " + pos + " loci.");
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

        	
			/* Get depth at this position. */
			depth = locus.getRecordAndPositions().size();
			nMapQBelowThreshold = 0;
			
			/* If depth ≥ minimum required depth to trust mass */
			if(depth >= CamoGeneFinder.MIN_DEPTH){

				/* Write dark and incomplete regions if large enough */
				if(consecDark >= CamoGeneFinder.MIN_REGION_SIZE){
					writeRegion(darkRegion, darkWriter);
				}
				if(consecInc >= CamoGeneFinder.MIN_REGION_SIZE) {
					writeRegion(incRegion, incWriter);
				}

				/* Clear regardless (i.e., even if the region wasn't large enough) */
				darkRegion.clear();
				consecDark = 0;
				incRegion.clear();
				consecInc = 0; 

				List<RecordAndOffset> recs = locus.getRecordAndPositions();
				int mapq;
				for(RecordAndOffset rec : recs){
					mapq = rec.getRecord().getMappingQuality();
					if(mapq <= CamoGeneFinder.MAPQ_THRESHOLD){
						nMapQBelowThreshold++;
					}
				}
				
				double percMapQ0 = nMapQBelowThreshold / depth * 100;
				if(percMapQ0 >= CamoGeneFinder.MIN_MAPQ_0_MASS){
					consecCamo++;
					camoRegion.add(camoRegionToString(contig, pos,
							nMapQBelowThreshold, depth, percMapQ0));
				}
				else if(consecCamo >= CamoGeneFinder.MIN_REGION_SIZE){
					writeRegion(camoRegion, camoWriter);
					camoRegion.clear();
					consecCamo = 0;
				}
			}
			else {
				/* Write camo region if large enough */
				if(consecCamo >= CamoGeneFinder.MIN_REGION_SIZE) {
					writeRegion(camoRegion, camoWriter);
				}
				if(consecInc >= CamoGeneFinder.MIN_REGION_SIZE) {
					writeRegion(incRegion, incWriter);
				}
				
				/* Clear regardless (i.e., even if the region wasn't large enough) */
				camoRegion.clear();
				consecCamo = 0;
				incRegion.clear();
				consecInc = 0;
				
				/* Save 'dark' regions with low coverage. */
				darkRegion.add(darkRegionToString(contig, pos, depth));
				consecDark++;
				
			}
			
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
			int depth) {

		/* Bed files are 0-based. locus.getPosition() returns 1-based. #Annoying */
		return contigName + "\t" + (position - 1) + "\t" + position +
				"\t" + depth + "\n";
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
