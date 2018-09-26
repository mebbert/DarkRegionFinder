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
import java.util.HashSet;
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
	private static int MAPQ_CAMO_THRESHOLD, DARK_DEPTH,
						MAPQ_DARK_THRESHOLD = 9, MIN_REGION_SIZE,
						MIN_DARK_MAPQ_MASS, MIN_CAMO_MAPQ_MASS,
                        MAX_ARRAY_SIZE = 10000;
	private static ValidationStringency SAM_VALIDATION_STRINGENCY;
	
	BufferedWriter camoWriter, darkWriter, incWriter;

	private SAMFileHeader header;
	private SamReader reader;
	
	private IndexedFastaSequenceFile hgRefReader;


    /**
     *
     * @param samFile
     * @param outCamoBed
     * @param outDarkBed
     * @param outIncBed
     * @param hgRef
     * @param mapQThreshold
     * @param minCamoMapQMass
     * @param minDarkMapQMass
     * @param minRegionSize
     * @param darkDepth
     * @param vs
     * @throws IOException
     */
	public CamoGeneFinder(final File samFile, final File outCamoBed,
			File outDarkBed, File outIncBed, File hgRef,
            final int mapQThreshold, final int minCamoMapQMass,
			final int minDarkMapQMass, final int minRegionSize,
            int darkDepth, final ValidationStringency vs/*, int startWalking,
			int endWalking*/) throws IOException {
		
		// Would be nice to be able to specify start/end locations
//		this.startWalking = startWalking;
//		this.endWalking = endWalking;
		
		camoWriter = new BufferedWriter(new OutputStreamWriter(
	              new FileOutputStream(outCamoBed), "utf-8"));
		camoWriter.write("chrom\tstart\tend\tnMapQBelowThreshold\tdepth\tpercMapQBelowThreshold\n");

		darkWriter = new BufferedWriter(new OutputStreamWriter(
	              new FileOutputStream(outDarkBed), "utf-8"));
		darkWriter.write("chrom\tstart\tend\tnMapQBelowThreshold\tdepth\tpercMapQBelowThreshold\n");

		incWriter = new BufferedWriter(new OutputStreamWriter(
	              new FileOutputStream(outIncBed), "utf-8"));
		incWriter.write("chrom\tstart\tend\n");

		CamoGeneFinder.MAPQ_CAMO_THRESHOLD = mapQThreshold;
		CamoGeneFinder.MIN_REGION_SIZE = minRegionSize;
		CamoGeneFinder.DARK_DEPTH = darkDepth;
		CamoGeneFinder.MIN_CAMO_MAPQ_MASS = minCamoMapQMass;
		CamoGeneFinder.MIN_DARK_MAPQ_MASS = minDarkMapQMass;
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
		int consecCamo = 0, consecDark = 0, consecInc = 0, pos,
				nMapQBelowCamoThreshold, mapq, nMapQBelowDarkThreshold;
		ArrayList<String> camoRegion = new ArrayList<String>(),
				darkRegion = new ArrayList<String>(),
				incRegion = new ArrayList<String>();
		HashSet<String> ignore = new HashSet<String>();
		String contig; byte[] bases; byte base;
		List<RecordAndOffset> recs = null;
		double percMapQBelowDarkThreshold, percMapQBelowCamoThreshold, depth;

		while(sli.hasNext()){

		    /* write out and clear regions if the arrays are getting too big (in order to save memory)
		       assuming MAX_ARRAY_SIZE > MIN_REGION_SIZE so we don't need to check consecCamo before printing
		     */
		    if ( consecInc > CamoGeneFinder.MIN_REGION_SIZE && incRegion.size() > CamoGeneFinder.MAX_ARRAY_SIZE) {
		        writeRegion(incRegion, incWriter);
		        incRegion.clear();
            }
            if ( consecDark > CamoGeneFinder.MIN_REGION_SIZE && darkRegion.size() > CamoGeneFinder.MAX_ARRAY_SIZE) {
                writeRegion(darkRegion, darkWriter);
                darkRegion.clear();
            }
            if ( consecCamo > CamoGeneFinder.MIN_REGION_SIZE && camoRegion.size() > CamoGeneFinder.MAX_ARRAY_SIZE) {
                writeRegion(camoRegion, camoWriter);
                camoRegion.clear();
            }

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


        	/* bases array contains only one element, so extract this base */
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

			
			/* Get number of reads with MAPQ ≤ threshold */
			recs = locus.getRecordAndPositions();
			nMapQBelowCamoThreshold = 0;
			nMapQBelowDarkThreshold = 0;
			for(RecordAndOffset rec : recs){
				mapq = rec.getRecord().getMappingQuality();
				if(mapq <= CamoGeneFinder.MAPQ_CAMO_THRESHOLD){
					nMapQBelowCamoThreshold++;
				}
				if(mapq <= CamoGeneFinder.MAPQ_DARK_THRESHOLD){
					nMapQBelowDarkThreshold++;
				}
			}


			percMapQBelowDarkThreshold = depth > 0 ? Math.round(nMapQBelowDarkThreshold / depth * 100) : -1;
			percMapQBelowCamoThreshold = depth > 0 ? Math.round(nMapQBelowCamoThreshold / depth * 100) : -1;

            /* Check if we're in a dark region
             * A region is 'dark' if depth is < DARK_DEPTH OR if the
             * percentage of reads with a MAPQ < 10 is ≥ MIN_DARK_MAPQ_MASS
             */
            if(depth <= CamoGeneFinder.DARK_DEPTH ||
                    percMapQBelowDarkThreshold >= CamoGeneFinder.MIN_DARK_MAPQ_MASS) {

                /* Save 'dark' regions with low coverage or high mass MapQ < 10. */
                darkRegion.add(darkRegionToString(contig, pos, nMapQBelowDarkThreshold,
                        depth, percMapQBelowDarkThreshold));
                consecDark++;

                /*  With in a dark region check if we're also in a camo region
                 *  camo regions are a subset of dark regions where depth is sufficient (> DARK_DEPTH) AND
                 *  the percentage of reads with MAPQ = 0 is ≥ MIN_CAMO_MAPQ_MASS
                 */
                if(depth > CamoGeneFinder.DARK_DEPTH &&
                        percMapQBelowCamoThreshold >= CamoGeneFinder.MIN_CAMO_MAPQ_MASS){

                    camoRegion.add(camoRegionToString(contig, pos,
                            nMapQBelowCamoThreshold, depth, percMapQBelowCamoThreshold));
                    consecCamo++;
                }
                else if(consecCamo >= CamoGeneFinder.MIN_REGION_SIZE){
                    /* just exited a Camo region, write out region if big enough */
                    writeRegion(camoRegion, camoWriter);
                    camoRegion.clear();
                    consecCamo = 0;
                }
                else {
                    /* Also clear if we didn't reach the region size */
                    camoRegion.clear();
                    consecCamo = 0;
                }

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
        /* Use StringBuilder to save memory */
        StringBuilder sb = new StringBuilder();
        sb.append(contigName).append("\t")
                .append(position - 1).append("\t")
                .append(position).append("\n");
        return sb.toString();
	}

    /**
     *
     * @param contigName
     * @param position
     * @param nMapQBelowDarkThreshold
     * @param depth
     * @param percentMapQBelowDarkThreshold
     * @return
     */
	private String darkRegionToString(String contigName, int position,
			int nMapQBelowDarkThreshold, double depth, double percentMapQBelowDarkThreshold) {

		/* Bed files are 0-based. locus.getPosition() returns 1-based. #Annoying */
        StringBuilder sb = new StringBuilder();
        sb.append(contigName).append("\t")
                .append(position - 1).append("\t")
                .append(position).append("\t")
                .append(nMapQBelowDarkThreshold).append("\t")
                .append((int) depth).append("\t")
                .append(percentMapQBelowDarkThreshold).append("\n");
        return sb.toString();
	}

    /**
     *
     * @param contigName
     * @param position
     * @param nMapQBelowCamoThreshold
     * @param depth
     * @param percentMapQBelowCamoThreshold
     * @return
     */
	private String camoRegionToString(String contigName, int position,
			int nMapQBelowCamoThreshold, double depth,
			double percentMapQBelowCamoThreshold) {

		/* Bed files are 0-based. locus.getPosition() returns 1-based. #Annoying */
        StringBuilder sb = new StringBuilder();
        sb.append(contigName).append("\t")
                .append(position - 1).append("\t")
                .append(position).append("\t")
                .append(nMapQBelowCamoThreshold).append("\t")
                .append((int) depth).append("\t")
                .append(percentMapQBelowCamoThreshold).append("\n");
        return sb.toString();
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
