package ebbert.cgf;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;

import org.apache.log4j.Logger;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

public class ThreadedDarkRegionFinder implements IDarkRegionFinder {

    private static int MAPQ_THRESHOLD, MIN_DEPTH,
            MIN_REGION_SIZE, MIN_MAPQ_MASS,
            MAX_ARRAY_SIZE = 10000;
    private static boolean EXCLUSIVE_REGIONS;
    private static ValidationStringency SAM_VALIDATION_STRINGENCY;
    private static String TMPDIR;
    private static int threads;

    private File lowMapQBed, lowDepthBed, incBed, SamFile, hgRef;

    private class IntervalLocusWalker implements Runnable {

        private Logger threadLogger = Logger.getLogger(IntervalLocusWalker.class);
        BufferedWriter tmpMapQWriter, tmpDepthWriter, tmpIncWriter;

        private SamReader reader;
        private SAMFileHeader header;

        private IndexedFastaSequenceFile hgRefReader;

        private IntervalList intervalList;

        IntervalLocusWalker(File tmpDepthBed, File tmpMapQBed, File tmpIncBed,
                                   List<Interval> intervalList) throws IOException {

            tmpDepthWriter = new BufferedWriter(new OutputStreamWriter(
                    new FileOutputStream(tmpDepthBed), StandardCharsets.UTF_8));

            tmpMapQWriter = new BufferedWriter(new OutputStreamWriter(
                    new FileOutputStream(tmpMapQBed), StandardCharsets.UTF_8));

            tmpIncWriter = new BufferedWriter(new OutputStreamWriter(
                    new FileOutputStream(tmpIncBed), StandardCharsets.UTF_8));

            this.reader = ThreadedDarkRegionFinder.openSam(SamFile,
                    ThreadedDarkRegionFinder.SAM_VALIDATION_STRINGENCY);
            this.header = reader.getFileHeader();
            this.hgRefReader = new IndexedFastaSequenceFile(hgRef);


            /*
            For some reason Samtools IntervalList object cannot be initialized with a List<Interval> ??
            To get around this we are initializing an IntervalList from SamHeader object,
            then subtracting that IntervalList from itself to get an empty IntervalList
            Then adding the List<Interval> to get an IntervalList object with the intervals we want
             */
            this.intervalList = new IntervalList(header);
            this.intervalList = IntervalList.subtract(this.intervalList, this.intervalList);
            this.intervalList.addall(intervalList);
        }

        @Override
        public void run() {
            try {
                    startWalkingByLocus();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        /**
         * @throws Exception throws possible exception
         */
        void startWalkingByLocus() throws Exception{

            SamLocusIterator sli;
            sli = new SamLocusIterator(reader, intervalList);

            /* Set the base quality score cutoff to 0
             * so the iterator won't try to validate base qualities. Any
             * secondary alignment won't have the qualities, and they're
             * all made up anyway.
             */
            int baseQualityScoreCutoff = 0;
            sli.setQualityScoreCutoff(baseQualityScoreCutoff);

            /* Walk along genome identifying 'dark' and 'camouflaged' regions */

            SamLocusIterator.LocusInfo locus;
            int consecLowDepth = 0, consecLowMapQ = 0, consecInc = 0, pos,
                    nMapQBelowThreshold, mapq;
            ArrayList<String> lowDepthRegion = new ArrayList<>(),
                    lowMapQRegion = new ArrayList<>(),
                    incRegion = new ArrayList<>();
            String contig; byte[] bases; byte base;
            List<SamLocusIterator.RecordAndOffset> recs;
            double percMapQBelowThreshold, depth;
            boolean low_depth;

            while(sli.hasNext()){

                /* write out and clear regions if the arrays are getting too big (in order to save memory)
                 */
                if ( consecInc > ThreadedDarkRegionFinder.MIN_REGION_SIZE && incRegion.size() > ThreadedDarkRegionFinder.MAX_ARRAY_SIZE) {
                    writeRegion(incRegion, tmpIncWriter);
                    incRegion.clear();
                }
                if ( consecLowDepth > ThreadedDarkRegionFinder.MIN_REGION_SIZE && lowDepthRegion.size() > ThreadedDarkRegionFinder.MAX_ARRAY_SIZE) {
                    writeRegion(lowDepthRegion, tmpDepthWriter);
                    lowDepthRegion.clear();
                }
                if ( consecLowMapQ > ThreadedDarkRegionFinder.MIN_REGION_SIZE && lowMapQRegion.size() > ThreadedDarkRegionFinder.MAX_ARRAY_SIZE) {
                    writeRegion(lowMapQRegion, tmpMapQWriter);
                    lowMapQRegion.clear();
                }


                locus = sli.next();

                contig = locus.getSequenceName();


                /* Returns 1-based position */
                pos = locus.getPosition();

                /* Expects 1-based position (inclusive to inclusive) */
                bases = hgRefReader.getSubsequenceAt(contig, pos, pos).getBases();

                /* Track progress */
                if(pos % 1000000 == 0){
                    threadLogger.debug("Assessed " + pos + " loci on " + contig);
                }


                /* bases array contains only one element, so extract this base */
                base = bases[0];

                /* Record incomplete genomic regions (i.e., 'N') */
                if(base == 'N' || base == 'n'){
                    incRegion.add(incompleteRegionToString(contig, pos));
                    consecInc++;

                    /* Write dark regions if large enough */
                    if(consecLowDepth >= ThreadedDarkRegionFinder.MIN_REGION_SIZE){
                        writeRegion(lowDepthRegion, tmpDepthWriter);
                    }
                    if(consecLowMapQ >= ThreadedDarkRegionFinder.MIN_REGION_SIZE){
                        writeRegion(lowMapQRegion, tmpMapQWriter);
                    }

                    /* Clear regardless (i.e., even if the region wasn't large enough) */
                    lowMapQRegion.clear();
                    consecLowMapQ = 0;
                    lowDepthRegion.clear();
                    consecLowDepth = 0;

                    continue;
                }


                /* Write incomplete regions if large enough. Clear in either case. */
                if(consecInc >= ThreadedDarkRegionFinder.MIN_REGION_SIZE) {
                    writeRegion(incRegion, tmpIncWriter);
                }

                /* Clear regardless because we know we're outside an incomplete
                 * region
                 */
                incRegion.clear();
                consecInc = 0;


                low_depth = false;

                /* Get depth at this position. */
                depth = locus.getRecordAndPositions().size();


                /* Get number of reads with MAPQ ≤ threshold */
                recs = locus.getRecordAndPositions();
                nMapQBelowThreshold = 0;
                for(SamLocusIterator.RecordAndOffset rec : recs){
                    mapq = rec.getRecord().getMappingQuality();
                    if(mapq <= ThreadedDarkRegionFinder.MAPQ_THRESHOLD){
                        nMapQBelowThreshold++;
                    }
                }


                percMapQBelowThreshold = depth > 0 ? Math.round(nMapQBelowThreshold / depth * 100) : -1;

                /* Check if we're in a low depth Dark Region
                 * A region is 'dark' by low_depth if depth is < MIN_DEPTH
                 */
                if(depth <= ThreadedDarkRegionFinder.MIN_DEPTH ) {

                    /* Save low-depth 'dark' regions with low coverage */
                    low_depth = true;
                    lowDepthRegion.add(lowDepthRegionToString(contig, pos, nMapQBelowThreshold,
                            depth, percMapQBelowThreshold));
                    consecLowDepth++;
                }
                else if ( consecLowDepth > ThreadedDarkRegionFinder.MIN_REGION_SIZE ) {
                    /* write dark region then clear */
                    writeRegion(lowDepthRegion, tmpDepthWriter);

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
                if (ThreadedDarkRegionFinder.EXCLUSIVE_REGIONS && low_depth ) {

                    /* print out lowMapQ Region if long enough */
                    if ( consecLowMapQ > ThreadedDarkRegionFinder.MIN_REGION_SIZE) {
                        writeRegion(lowMapQRegion, tmpMapQWriter);
                    }

                    /* clear lowMapQ Region buffer regardless of length */
                    lowMapQRegion.clear();
                    consecLowMapQ = 0;
                }
                else if ( percMapQBelowThreshold >= ThreadedDarkRegionFinder.MIN_MAPQ_MASS) {

                    /* Save lowMapQ 'dark' region which has at mass > MIN_MAPQ_MASS of reads with mapq < MAPQ_THRESHOLD */
                    lowMapQRegion.add(lowMapQRegionToString(contig, pos,
                            nMapQBelowThreshold, depth, percMapQBelowThreshold));
                    consecLowMapQ++;

                }
                else if ( consecLowMapQ > ThreadedDarkRegionFinder.MIN_REGION_SIZE ) {
                    /* write out and clear lowMapQ region since it is long enough */
                    writeRegion(lowMapQRegion, tmpMapQWriter);
                    lowMapQRegion.clear();
                    consecLowMapQ = 0;
                }
                else {
                    lowMapQRegion.clear();
                    consecLowMapQ = 0;
                }

            }

            /* Write regions if large enough */
            if(consecLowDepth >= ThreadedDarkRegionFinder.MIN_REGION_SIZE){
                writeRegion(lowDepthRegion, tmpDepthWriter);
            }
            if(consecLowMapQ >= ThreadedDarkRegionFinder.MIN_REGION_SIZE) {
                writeRegion(lowMapQRegion, tmpMapQWriter);
            }
            if(consecInc >= ThreadedDarkRegionFinder.MIN_REGION_SIZE) {
                writeRegion(incRegion, tmpIncWriter);
            }

            tmpDepthWriter.close();
            tmpMapQWriter.close();
            tmpIncWriter.close();

            sli.close();
        }

        /**
         * @param contigName name of Chrom
         * @param position position of incomplete base
         * @return String to be written in Incomplete Bed
         */
        private String incompleteRegionToString(String contigName, int position) {

            /* Bed files are 0-based. locus.getPosition() returns 1-based. #Annoying */
            /* Use StringBuilder to save memory */
            return contigName + "\t" +
                    (position - 1) + "\t" +
                    position + "\n";
        }

        /**
         *
         * @param contigName chromosome name
         * @param position start position
         * @param nMapQBelowThreshold number of reads at position with MAPQ less than the MAPQ threshold
         * @param depth number of reads at position
         * @param percentMapQBelowThreshold percentage of reads with MAPQ less than the MAPQ threshold
         * @return string to be written to lowDepth Bed
         */
        private String lowDepthRegionToString(String contigName, int position,
                                              int nMapQBelowThreshold, double depth, double percentMapQBelowThreshold) {

            /* Bed files are 0-based. locus.getPosition() returns 1-based. #Annoying */
            return contigName + "\t" +
                    (position - 1) + "\t" +
                    position + "\t" +
                    nMapQBelowThreshold + "\t" +
                    (int) depth + "\t" +
                    percentMapQBelowThreshold + "\n";
        }

        /**
         *
         * @param contigName chromosome name
         * @param position position of dark base
         * @param nMapQBelowThreshold number of reads with MAPQ lower than MAPQ threshold
         * @param depth number of reads at base
         * @param percentMapQBelowThreshold percent of reads with MAPQ lower than MAPQ threshold at this base
         * @return string ready to be written in lowMAPQ Region bed
         */
        private String lowMapQRegionToString(String contigName, int position,
                                             int nMapQBelowThreshold, double depth,
                                             double percentMapQBelowThreshold) {

            /* Bed files are 0-based. locus.getPosition() returns 1-based. #Annoying */
            return contigName + "\t" +
                    (position - 1) + "\t" +
                    position + "\t" +
                    nMapQBelowThreshold + "\t" +
                    (int) depth + "\t" +
                    percentMapQBelowThreshold + "\n";
        }

        /**
         * @param lowMapQRegions Array of Low MAPQ region strings
         * @param writer BufferedWriter used to write those strings to the out file
         * @throws IOException thrown if IOException occurs while writing file
         */
        private void writeRegion(ArrayList<String> lowMapQRegions,
                                 BufferedWriter writer) throws IOException {
            for(String s : lowMapQRegions){
                writer.write(s);
            }
        }
    }

    /**
     *
     * @param samFile File object for the sample Bam
     * @param outDepthBed File object where the Low Depth Bed will be written
     * @param outMapQBed File object where the Low MAPQ file will be written
     * @param outIncBed File object where the Incomplete Bed will be written
     * @param hgRef File containing the human reference
     * @param mapQThreshold MAP threshold
     * @param minMapQMass The minimum percentage of reads below threshold to be considered dark
     * @param minRegionSize minimum number of contiguous dark bases needed to report region as dark
     * @param minDepth minimum depth under which a region is considered dark by depth
     * @param vs validation stringency for parsing bam file
     * @throws IOException thrown if there is a problem opening / writing to any of the files supplied
     */
    ThreadedDarkRegionFinder(final File samFile, final File outDepthBed,
                                    File outMapQBed, File outIncBed, File hgRef,
                                    final int mapQThreshold, final int minMapQMass, final int minRegionSize,
                                    int minDepth, final boolean exclusiveRegions, final ValidationStringency vs,
                                    String TMPDIR, int threads) throws IOException {

        // Would be nice to be able to specify start/end locations
//		this.startWalking = startWalking;
//		this.endWalking = endWalking;

        lowDepthBed = outDepthBed;
        BufferedWriter lowDepthWriter = new BufferedWriter(new OutputStreamWriter(
                new FileOutputStream(outDepthBed), StandardCharsets.UTF_8));
        lowDepthWriter.write("chrom\tstart\tend\tnMapQBelowThreshold\tdepth\tpercMapQBelowThreshold\n");
        lowDepthWriter.close();

        lowMapQBed = outMapQBed;
        BufferedWriter lowMapQWriter = new BufferedWriter(new OutputStreamWriter(
                new FileOutputStream(outMapQBed), StandardCharsets.UTF_8));
        lowMapQWriter.write("chrom\tstart\tend\tnMapQBelowThreshold\tdepth\tpercMapQBelowThreshold\n");
        lowMapQWriter.close();

        incBed = outIncBed;
        BufferedWriter incWriter = new BufferedWriter(new OutputStreamWriter(
                new FileOutputStream(outIncBed), StandardCharsets.UTF_8));
        incWriter.write("chrom\tstart\tend\n");
        incWriter.close();


        ThreadedDarkRegionFinder.MAPQ_THRESHOLD = mapQThreshold;
        ThreadedDarkRegionFinder.MIN_REGION_SIZE = minRegionSize;
        ThreadedDarkRegionFinder.MIN_DEPTH = minDepth;
        ThreadedDarkRegionFinder.MIN_MAPQ_MASS = minMapQMass;
        ThreadedDarkRegionFinder.EXCLUSIVE_REGIONS = exclusiveRegions;
        ThreadedDarkRegionFinder.SAM_VALIDATION_STRINGENCY = vs;
        ThreadedDarkRegionFinder.TMPDIR = TMPDIR;

        File tmpDir = new File(TMPDIR);
        if (!tmpDir.exists()) {
            if (!tmpDir.mkdir()) throw new IOException("ERROR: could not create TMP directory.");
        }

        ThreadedDarkRegionFinder.threads = threads;

        this.hgRef = hgRef;

        this.SamFile = samFile;

    }

    private List<List<Interval>> getChunckedRefIntervals() throws FileNotFoundException {
        List<List<Interval>> chunkedRefIntervals = new ArrayList<>();
        SAMSequenceDictionary sequenceDictionary = new IndexedFastaSequenceFile(hgRef).getSequenceDictionary();


        double total = 0;
        for(SAMSequenceRecord contig : sequenceDictionary.getSequences()) {
            total += contig.getSequenceLength();
        }

        //Number of bases per thread
        long numBases = (long) ((total + ThreadedDarkRegionFinder.threads - 1) / ThreadedDarkRegionFinder.threads);

        List<Interval> threadIntervals = new ArrayList<>();
        long targetBases = numBases;
        for(SAMSequenceRecord contig : sequenceDictionary.getSequences()) {
            int start = 1;
            while(contig.getSequenceLength() - start + 1 > targetBases) {
                threadIntervals.add(new Interval(contig.getSequenceName(), start, (int) (start + targetBases - 1)));
                chunkedRefIntervals.add(threadIntervals);
                threadIntervals = new ArrayList<>();
                start += targetBases;
                targetBases = numBases;
            }
            if (contig.getSequenceLength() - start + 1 == targetBases) {
                threadIntervals.add(new Interval(contig.getSequenceName(), start, contig.getSequenceLength()));
                chunkedRefIntervals.add(threadIntervals);
                targetBases = numBases;
            }
            else {
                threadIntervals.add(new Interval(contig.getSequenceName(), start, contig.getSequenceLength()));
                targetBases -= (contig.getSequenceLength() - start + 1);
            }
        }
        if (!threadIntervals.isEmpty()) chunkedRefIntervals.add(threadIntervals);

        return chunkedRefIntervals;
    }

    public void startWalkingByLocus() throws Exception {
        List<List<Interval>> chunkedThreadIntervals = getChunckedRefIntervals();
        List<File> tmpDepthBeds = new ArrayList<>();
        List<File> tmpMapQBeds = new ArrayList<>();
        List<File> tmpIncBeds = new ArrayList<>();
        String depthPrefix = String.format("%s/tmp.depth.", ThreadedDarkRegionFinder.TMPDIR);
        String mapqPrefix = String.format("%s/tmp.mapq.", ThreadedDarkRegionFinder.TMPDIR);
        String incPrefix = String.format("%s/tmp.inc.", ThreadedDarkRegionFinder.TMPDIR);

        List<Thread> threads = new ArrayList<>();
        for (int i = 0; i < ThreadedDarkRegionFinder.threads; i++) {
            File tmpDepth = new File(depthPrefix + i);
            tmpDepthBeds.add(tmpDepth);
            File tmpMapQ = new File(mapqPrefix + i);
            tmpMapQBeds.add(tmpMapQ);
            File tmpInc = new File(incPrefix + i);
            tmpIncBeds.add(tmpInc);

           threads.add(new Thread( new IntervalLocusWalker(tmpDepth, tmpMapQ, tmpInc, chunkedThreadIntervals.get(i))));
           threads.get(i).start();
        }

        for (int i = 0; i < ThreadedDarkRegionFinder.threads; i++) {
            threads.get(i).join();
        }

        cat(tmpDepthBeds, lowDepthBed);
        cat(tmpMapQBeds, lowMapQBed);
        cat(tmpIncBeds, incBed);

    }

    private void cat(List<File> files, File outFile) throws IOException {
        OutputStream out = new FileOutputStream(outFile, true);
        byte[] buf = new byte[1000]; //read 1000 bytes at a time so as not to load too much of the file in memory at once
        for (File f : files) {
            InputStream in = new FileInputStream(f);
            int b;
            while ( (b = in.read(buf)) >= 0) {
                out.write(buf, 0, b);
                out.flush();
            }
            f.delete();
        }
        out.close();
    }


    /**
     *
     * Open a SAM/BAM file for reading and return the SamReader obj
     *
     * @param samFile SAM/BAM file for reading
     * @param vs validation stringency
     * @return SamReader
     */
    private static SamReader openSam(final File samFile, ValidationStringency vs) {

        final SamReaderFactory factory =
                SamReaderFactory.makeDefault()
                        .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS,
                                SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
                        .validationStringency(vs);

        return factory.open(samFile);
    }
}
