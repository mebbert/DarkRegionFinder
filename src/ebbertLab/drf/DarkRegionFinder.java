/**
 * 
 */
package ebbertLab.drf;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.TreeSet;
import java.util.zip.GZIPOutputStream;

import org.apache.log4j.Logger;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset;


/**
 * @author markebbert
 *
 */
public class DarkRegionFinder {
	
	private static Logger logger = Logger.getLogger(DarkRegionFinder.class);
	private static int MAPQ_THRESHOLD, MIN_DEPTH,
						MIN_REGION_SIZE, MIN_MAPQ_MASS,
                        MAX_ARRAY_SIZE = 10000;
	private static boolean EXCLUSIVE_REGIONS;
	private static IntervalList intervalList;
	
	Writer lowMapQWriter, lowDepthWriter, incWriter;

	private SAMFileHeader header;
	private SamReader samReader;
	
	private IndexedFastaSequenceFile hgRefReader;
	private SAMSequenceDictionary hgRefDictionary;

    /**
     *
     * @param samFile
     * @param outDepthBed
     * @param outMapQBed
     * @param outIncBed
     * @param hgRef
     * @param mapQThreshold
     * @param minMapQMass
     * @param minRegionSize
     * @param minDepth
     * @param vs
     * @throws IOException
     * @throws ReferenceDictionaryNotFoundException 
     */
	public DarkRegionFinder(final File samFile, final File outDepthBed,
                        File outMapQBed, File outIncBed, File hgRef,
                        final int mapQThreshold, final int minMapQMass, final int minRegionSize,
                        int minDepth, final boolean exclusiveRegions, final ValidationStringency vs
                        , List<String> intervalStringList) throws IOException {
		

		lowDepthWriter = new OutputStreamWriter(new GZIPOutputStream(
	              new FileOutputStream(outDepthBed)), "utf-8");
		// lowDepthWriter.write("chrom\tstart\tend\tnMapQBelowThreshold\tdepth\tpercMapQBelowThreshold\n");

		lowMapQWriter = new OutputStreamWriter(new GZIPOutputStream(
	              new FileOutputStream(outMapQBed)), "utf-8");
		// lowMapQWriter.write("chrom\tstart\tend\tnMapQBelowThreshold\tdepth\tpercMapQBelowThreshold\n");

		incWriter = new OutputStreamWriter(new GZIPOutputStream(
	              new FileOutputStream(outIncBed)), "utf-8");
		// incWriter.write("chrom\tstart\tend\n");

		DarkRegionFinder.MAPQ_THRESHOLD = mapQThreshold;
		DarkRegionFinder.MIN_REGION_SIZE = minRegionSize;
		DarkRegionFinder.MIN_DEPTH = minDepth;
		DarkRegionFinder.MIN_MAPQ_MASS = minMapQMass;
		DarkRegionFinder.EXCLUSIVE_REGIONS = exclusiveRegions;
		
		this.hgRefReader = new IndexedFastaSequenceFile(hgRef);
		this.hgRefDictionary = hgRefReader.getSequenceDictionary();
		
		if(null == hgRefDictionary) {
			String err = "The reference provided does not have a dictionary (.dict) file. Generate"
					+ " the dictionary before proceeding.";
			logger.error(err);
			throw new IOException(err);
		}

		this.samReader = DarkRegionFinder.openSam(samFile, vs, hgRef);
		this.header = samReader.getFileHeader();

		DarkRegionFinder.intervalList = generateIntervalList(intervalStringList, header);

		
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

		SamLocusIterator sli;
		
		/*
		 * If the user specified an interval list, only read those intervals.
		 * Otherwise, start from the beginning.
		 */
		if(null != DarkRegionFinder.intervalList && DarkRegionFinder.intervalList.size() > 0) {
			sli = new SamLocusIterator(samReader, DarkRegionFinder.intervalList, true);
			
			/* Tell SLI to include indels. This is necessary to report
			 * reads that overlap deletions.
			 */
			sli.setIncludeIndels(true);
			
			for(Interval interval : intervalList) {
				System.out.println("Interval: " + interval.toString());
			}
		}
		else {
			sli = new SamLocusIterator(samReader);
		}
		
		/*
		 * Set the base quality score cutoff to 0 
		 * so the iterator won't try to validate base qualities. Any
		 * secondary alignment won't have the qualities, and they're
		 * all made up anyway.
		 */
		int baseQualityScoreCutoff = 0;
		sli.setQualityScoreCutoff(baseQualityScoreCutoff);
		
		
		/* Walk along genome identifying 'dark' and 'camouflaged' regions */

		LocusInfo locus;
		int consecLowDepth = 0, consecLowMapQ = 0, consecInc = 0, pos,
				nMapQBelowThreshold, mapq, nLociAssessed = 0;
		ArrayList<String> lowDepthRegion = new ArrayList<String>(),
				lowMapQRegion = new ArrayList<String>(),
				incRegion = new ArrayList<String>();
		HashSet<String> ignore = new HashSet<String>();
		String contig; byte[] bases; byte base;
		List<RecordAndOffset> recs = null;
		double percMapQBelowThreshold, depthExcludingIndels, totalDepthIncludingIndels, nDeletedInRecord;
		boolean low_depth;

		while(sli.hasNext()){

		    /* write out and clear regions if the arrays are getting too big (in order to save memory)
		     */
		    if ( consecInc > DarkRegionFinder.MIN_REGION_SIZE && incRegion.size() > DarkRegionFinder.MAX_ARRAY_SIZE) {
		        writeRegion(incRegion, incWriter);
		        incRegion.clear();
            }
            if ( consecLowDepth > DarkRegionFinder.MIN_REGION_SIZE && lowDepthRegion.size() > DarkRegionFinder.MAX_ARRAY_SIZE) {
                writeRegion(lowDepthRegion, lowDepthWriter);
                lowDepthRegion.clear();
            }
            if ( consecLowMapQ > DarkRegionFinder.MIN_REGION_SIZE && lowMapQRegion.size() > DarkRegionFinder.MAX_ARRAY_SIZE) {
                writeRegion(lowMapQRegion, lowMapQWriter);
                lowMapQRegion.clear();
            }

			locus = sli.next();
			
			contig = locus.getSequenceName();
			
			/*
			 * If this contig is not in the ref, then skip. There's probably
			 * a way to skip an entire contig. A good to-do.
			 */
			if(ignore.contains(contig)) {
				continue;
			}

			/* Returns 1-based position */
			pos = locus.getPosition();  

			
			/* Ensure sequence is present in provided reference */
			if(null == this.hgRefDictionary.getSequence(contig)){
				logger.warn("BAM file contains alignments for " + contig
						+ " but this sequence was not found in the provided"
						+ " reference. Skipping.");
				ignore.add(contig);
				continue;
			}

			/* Expects 1-based position (inclusive to inclusive) */
			bases = hgRefReader.getSubsequenceAt(contig, pos, pos).getBases();

			/* Track progress */ 
        	if(nLociAssessed > 0 && nLociAssessed % 1000000 == 0){
        		logger.debug("Assessed " + nLociAssessed + " loci on " + contig);
        	}		


        	/* bases array contains only one element, so extract this base */
			base = bases[0];

			/* Record incomplete genomic regions (i.e., 'N') */
			if(base == 'N' || base == 'n'){
				incRegion.add(incompleteRegionToString(contig, pos));
				consecInc++;

				/* Write dark regions if large enough */
				if(consecLowDepth >= DarkRegionFinder.MIN_REGION_SIZE){
					writeRegion(lowDepthRegion, lowDepthWriter);
				}
				if(consecLowMapQ >= DarkRegionFinder.MIN_REGION_SIZE){
					writeRegion(lowMapQRegion, lowMapQWriter);
				}

				/* Clear regardless (i.e., even if the region wasn't large enough) */
				lowMapQRegion.clear();
				consecLowMapQ = 0;
				lowDepthRegion.clear();
				consecLowDepth = 0;

				continue;
			}
	

			/* Write incomplete regions if large enough. Clear in either case. */
			if(consecInc >= DarkRegionFinder.MIN_REGION_SIZE) {
				writeRegion(incRegion, incWriter);
			}

			/* Clear regardless because we know we're outside an incomplete
			 * region
			 */
			incRegion.clear();
			consecInc = 0; 		


			low_depth = false;

			/* Get depth, nOverlappingReads, and nDeletedInRecord at this position. */
			depthExcludingIndels = locus.getRecordAndOffsets().size();
			totalDepthIncludingIndels = locus.size();
			nDeletedInRecord = locus.getDeletedInRecord().size();
			

			/* 
			 * Get number of reads with MAPQ ≤ threshold 
			 */
			
			/*
			 * Get all records at this position, including those
			 * that have a deletion at this position. Records/reads that have
			 * a deletion are not included in the regular list of records, so
			 * have to add them separately (locus.getDeletedInRecord()).
			 */
			recs = new ArrayList<RecordAndOffset>(locus.getRecordAndOffsets());
			recs.addAll(locus.getDeletedInRecord()); 
			nMapQBelowThreshold = 0;
			for(RecordAndOffset rec : recs){
				mapq = rec.getRecord().getMappingQuality();
				if(mapq <= DarkRegionFinder.MAPQ_THRESHOLD){
					nMapQBelowThreshold++;
				}
			}


			percMapQBelowThreshold = totalDepthIncludingIndels > 0 ? Math.round(nMapQBelowThreshold / totalDepthIncludingIndels * 100) : -1;

            /* Check if we're in a low depth Dark Region
             * A region is 'dark' by low_depth if depth is < MIN_DEPTH
             */
            if(depthExcludingIndels <= DarkRegionFinder.MIN_DEPTH ) {

                /* Save low-depth 'dark' regions with low coverage */
                low_depth = true;
                lowDepthRegion.add(lowDepthRegionToString(contig, pos, nMapQBelowThreshold,
                        depthExcludingIndels, nDeletedInRecord, totalDepthIncludingIndels, percMapQBelowThreshold));
                consecLowDepth++;
            }
            else if ( consecLowDepth > DarkRegionFinder.MIN_REGION_SIZE ) {
                /* write dark region then clear */
                writeRegion(lowDepthRegion, lowDepthWriter);

                lowDepthRegion.clear();
                consecLowDepth = 0;
            }
            else {
                lowDepthRegion.clear();
                consecLowDepth = 0;
            }


            /* check if Exclusive and already in dark:
             * if Exclusive is true and locus was already in low_depth, cannot be low mapQ so write out low MapQ and clear
             * else if not exclusive or not low_depth check if it is a low MapQ region
             */
            if (DarkRegionFinder.EXCLUSIVE_REGIONS && low_depth ) {

                /* print out lowMapQ Region if long enough */
                if ( consecLowMapQ > DarkRegionFinder.MIN_REGION_SIZE) {
                    writeRegion(lowMapQRegion, lowMapQWriter);
                }

                /* clear lowMapQ Region buffer regardless of length */
                lowMapQRegion.clear();
                consecLowMapQ = 0;
            }
            else if (percMapQBelowThreshold >= DarkRegionFinder.MIN_MAPQ_MASS) {

                /* Save lowMapQ 'dark' region which has at mass > MIN_MAPQ_MASS of reads with mapq < MAPQ_THRESHOLD */
                lowMapQRegion.add(lowMapQRegionToString(contig, pos,
                        nMapQBelowThreshold, depthExcludingIndels, nDeletedInRecord,
                        totalDepthIncludingIndels, percMapQBelowThreshold));
                consecLowMapQ++;

            }
            else if ( consecLowMapQ > DarkRegionFinder.MIN_REGION_SIZE ) {
                /* write out and clear lowMapQ region since it is long enough */
                writeRegion(lowMapQRegion, lowMapQWriter);
                lowMapQRegion.clear();
                consecLowMapQ = 0;
            }
            else {
                lowMapQRegion.clear();
                consecLowMapQ = 0;
            }
            
            nLociAssessed++;
		}
		
		/* Write regions if large enough */
        if(consecLowDepth >= DarkRegionFinder.MIN_REGION_SIZE){
            writeRegion(lowDepthRegion, lowDepthWriter);
        }
		if(consecLowMapQ >= DarkRegionFinder.MIN_REGION_SIZE) {
			writeRegion(lowMapQRegion, lowMapQWriter);
		}
		if(consecInc >= DarkRegionFinder.MIN_REGION_SIZE) {
			writeRegion(incRegion, incWriter);
		}

		lowDepthWriter.close();
		lowMapQWriter.close();
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
     * Create a string for this low-depth region for output in a .bed file. Will
     * convert to 0-based positions since .bed files are 0-based.
     * 
     * @param contigName
     * @param position
     * @param nMapQBelowThreshold
     * @param depthExcludingIndels
     * @param percentMapQBelowThreshold
     * @return
     */
	private String lowDepthRegionToString(String contigName, int position,
			int nMapQBelowThreshold, double depthExcludingIndels, double nDeletedInRecord,
			double totalDepthIncludingIndels, double percentMapQBelowThreshold) {

		/* Bed files are 0-based. locus.getPosition() returns 1-based. #Annoying */
		int posZeroBased = position - 1;

        StringBuilder sb = new StringBuilder();
        sb.append(contigName).append("\t")
                .append(posZeroBased).append("\t")
                .append(posZeroBased).append("\t")
                .append(nMapQBelowThreshold).append("\t")
                .append(percentMapQBelowThreshold).append("\t")
                .append((int) depthExcludingIndels).append("\t")
                .append(nDeletedInRecord).append("\t")
                .append(totalDepthIncludingIndels).append("\n");
        return sb.toString();
	}

    /**
     *
     * @param contigName
     * @param position
     * @param nMapQBelowThreshold
     * @param depthExcludingIndels
     * @param percentMapQBelowThreshold
     * @return
     */
	private String lowMapQRegionToString(String contigName, int position,
			int nMapQBelowThreshold, double depthExcludingIndels, double nDeletedInRecord,
			double totalDepthIncludingIndels, double percentMapQBelowThreshold) {

		/* Bed files are 0-based. locus.getPosition() returns 1-based. #Annoying */
        StringBuilder sb = new StringBuilder();
        sb.append(contigName).append("\t")
                .append(position - 1).append("\t")
                .append(position).append("\t")
                .append(nMapQBelowThreshold).append("\t")
                .append(percentMapQBelowThreshold).append("\t")
                .append((int) depthExcludingIndels).append("\t")
                .append(nDeletedInRecord).append("\t")
                .append(totalDepthIncludingIndels).append("\n");
        return sb.toString();
	}
	
	/**
	 * @param lowMapQRegions
	 * @param writer
	 * @throws IOException
	 */
	private void writeRegion(ArrayList<String> lowMapQRegions,
			Writer writer) throws IOException {
		for(String s : lowMapQRegions){
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
	 private static SamReader openSam(final File samFile, ValidationStringency vs, File hgRef) {
		 
		final SamReaderFactory factory =
				  SamReaderFactory.makeDefault()
					  .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS,
							  SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
					  .validationStringency(vs)
					  .referenceSource(new ReferenceSource(hgRef));
		

		final SamReader samReader = factory.open(samFile);
		
		return samReader;
	 }
	 
		
	/**
	 * 
	 * Convert the List<String> of intervals to an IntervalList of Interval objects
	 * 
	 * @param intervalStringList
	 * @param header
	 * @return
	 */
	private static IntervalList generateIntervalList(List<String> intervalStringList, SAMFileHeader header) {
		
		if(null == intervalStringList) {
			return null;
		}

		IntervalList intervalList = new IntervalList(header);
		String contig, range;
		int start, end, lastColon;
		String[] rangeToks;
		for(String interval : intervalStringList) {
			
			/*
			 * Some contigs contain colon ':' in the contig name
			 * (e.g., 'HLA-A*01:01:01:02N:1-3291'), so we will get the last
			 * index of ';'.
			 */
			lastColon = interval.lastIndexOf(':');
			contig = interval.substring(0, lastColon);
			range = interval.substring(lastColon + 1);
			
			rangeToks = range.split("-");
			
			start = Integer.parseInt(rangeToks[0]);
			end= Integer.parseInt(rangeToks[1]);
			
			intervalList.add(new Interval(contig, start, end));
		}
		
		return intervalList;
	}
    
}
