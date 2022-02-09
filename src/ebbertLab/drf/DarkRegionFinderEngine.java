/**
 * 
 */
package ebbertLab.drf;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.Random;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;

import htsjdk.samtools.ValidationStringency;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentGroup;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

/**
 * @author markebbert
 *
 */
public class DarkRegionFinderEngine {

	private static Logger logger = Logger.getLogger(DarkRegionFinderEngine.class);
	
	public DarkRegionFinderEngine() {
		return;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		BasicConfigurator.configure();
		DarkRegionFinderEngine cgfe = new DarkRegionFinderEngine();
		ArgumentParser parser = cgfe.init(args);

		cgfe.findDarkRegions(parser, args);
	}
	
	/**
	 * Init the parser options
	 * 
	 * @param args
	 */
	private ArgumentParser init(String[] args){
		ArgumentParser parser = ArgumentParsers.newArgumentParser("DarkRegionFinder");
		parser.description("The Dark Region Finder will identify"
				+ " regions where the genome is 'dark' or"
				+ " 'incomplete'. Dark regions are the regions of the genome where it is"
                + " generally possible to call variants using standard methods."
                + " We define dark regions in two categories:"
				+ " (1) Low-coverage regions where depth is less than"
				+ " the defined '--min-depth' threshold, and (2) low-MapQ"
				+ " regions where 0 < MAPQ < 10 (default) for a majority of reads in the"
				+ " pileup. Incomplete regions are simply those where"
				+ " nucleotides are unknown (i.e., 'N' or 'n'). Incomplete"
				+ " and dark regions are mutually exclusive (i.e.,"
				+ " incomplete are not included as a dark region, and vice-versa).");
		parser.defaultHelp(true);
		
		ArgumentGroup drfOptions = parser.addArgumentGroup("Dark Region Finder arguments");
		ArgumentGroup ioOptions = parser.addArgumentGroup("input/output arguments");
			
		/* Setup SVC options */
		drfOptions
				.addArgument("-s", "--min-region-size")
				.dest("MIN_SIZE")
				.metavar("SIZE")
				.setDefault(1)
				.type(Integer.class)
				.help("The minimum size dark region to consider. Details will"
						+ " still be written at the individual base level, but"
						+ " only for regions that meet this size requirement."
						+ " Regions can be merged using bedtools, if desired.");

		drfOptions
				.addArgument("-t", "--mapq-threshold")
				.dest("MAPQ_THRESHOLD")
				.metavar("THRESH")
				.setDefault(9)
				.type(Integer.class)
				.help("The MAPQ threshold (≤) at which a read is \'inadequately\'"
						+ " aligned and considered \'dark\'. Generally"
						+ " recommended to use MAPQ ≤ 9"
						+ " as this is the default MAPQ cut off used by GATK"
						+ " when filtering reads");

		drfOptions
				.addArgument("-m", "--min-mapq-mass")
				.dest("MIN_MAPQ_MASS")
				.metavar("MAPQ_MASS")
				.setDefault(90)
				.type(Integer.class)
				.help("The minimum percentage (≥) of reads below the"
						+ " --mapq-threshold for the locus to be considered dark."
						+ " Dark loci where the percentage is below this threshold"
						+ " will still be reported if the depth ≤ --min-depth");

		drfOptions
				.addArgument("-d", "--min-depth")
				.dest("MIN_DEPTH")
				.metavar("MIN_DEPTH")
				.setDefault(5)
				.type(Integer.class)
				.help("The depth (≤) at which a region is always considered"
						+ " 'dark'. The default value is meant to be a stringent"
						+ " cutoff to identify regions where were very few reads"
						+ " align (i.e., they're simply missing). Regions where"
						+ " depth ≤ --min-depth will be reported regardless of MAPQ mass");

		drfOptions
                .addArgument("-e","--region-exclusivity")
                .dest("EXCLUSIVE")
                .metavar("EXCLUSIVE")
                .action(Arguments.storeTrue())
                .type(Boolean.class)
                .help("Whether to treat regions as exclusive or not. If present, a locus"
                        + " will either be dark by low coverage or dark by low MAPQ"
                        + " (set by --mapq-threshold), but NOT both. In this case,"
                        + " the locus will be considered dark by low coverage."
                        + " Otherwise, a locus can be in both categories if"
                        + " depth ≤ --min-depth and MAPQ mass ≥ --min-mapq-mass"
                        + " thresholds");

		drfOptions
				.addArgument("-v", "--validation-stringency")
				.dest("STRINGENCY")
				.setDefault("STRICT")
				.choices("STRICT", "LENIENT", "SILENT")
				.type(String.class)
				.help("The validation stringency when parsing a SAM/BAM"
						+ " file. 'STRICT' will throw errors if something "
						+ " is amiss, 'LENIENT' will give warnings but continue,"
						+ " and 'SILENT' will continue AND keep our mouth shut.");
			
		/* Setup IO options */
		ioOptions
				.addArgument("-i", "--input")
				.dest("SAM")
				.metavar("SAM/BAM")
				.type(String.class)
				.required(true)
				.help("The input file. This can be a SAM or BAM file.");
		
		ioOptions
				.addArgument("-g", "--human-ref")
				.dest("HG_REF")
				.type(String.class)
				.required(true)
				.help("The human genome reference file. Must also be indexed "
						+ "by 'samtools faidx' and have a GATK/Picard sequence"
						+ " dictionary (e.g., gatk CreateSequenceDictionary -R <ref.fa>).");

		ioOptions
				.addArgument("-c", "--low-coverage-bed-output")
				.dest("LOW_COV_BED")
				.type(String.class)
				.setDefault("low_coverage.dark.bed")
				.help("The output BED file for low-coverage dark regions. Low coverage dark"
						+ " regions are those where depth ≤ --min_depth. Columns for this file"
						+ " are: chromosome, start, end, MapQBelowThreshold, depth, percMapQBelowThreshold");

		ioOptions
				.addArgument("-a", "--low-mapq-bed-output")
				.dest("LOW_MAPQ_BED")
				.type(String.class)
				.setDefault("low_mapq.dark.bed")
				.help("The output BED file for low MAPQ dark regions. Low MAPQ regions are"
						+ " those with low mapping quality (loci with"
						+ " MAPQ ≤ --mapq_threshold and MAPQ mass ≥ --min_mapq_mass)."
						+ " Columns for this file are: chromosome, start, end, MapQBelowThreshold,"
						+ "depth, percMapQBelowThreshold");

		ioOptions
				.addArgument("-n", "--incomplete-bed-output")
				.dest("INC_BED")
				.type(String.class)
				.setDefault("incomplete.bed")
				.help("The output BED file for incomplete regions. Incomplete"
						+ " regions are those where the bases are unknown"
						+ " (i.e., 'N' or 'n'). Columns for this file are: chromosome, start, end.");
		
		ioOptions
				.addArgument("-L", "--interval-list")
				.dest("INTERVAL_LIST")
				.type(String.class)
				.nargs("+")
				.help("Specific intervals to include. Intervals are 0-based and should"
						+ " be formatted the same as samtools intervals (i.e.,"
						+ " <contig_name>:<start>-<end>; e.g., chr1:207496157-207641765),"
						+ " where <start> is inclusive while <end> is exclusive."
						+ " The contig name should reflect what is used in the reference"
						+ " genome the sample was aligned to (i.e., with or without 'chr',"
						+ " etc.)."
						+ "\n\nIf this argument is not provided, DRF will start at the"
						+ " beginning of the genome (based on how the reference and bam are"
						+ " sorted. This argument makes it possible to parallelize DRF by"
						+ " submitting multiple jobs with different regions. We chose this"
						+ " approach because it makes it easier to split jobs across nodes"
						+ " in a computer cluster, rather than having a single job manage"
						+ " all treads. The results will need to be combined for the sample."
						+ " \n\nIf this argument is provided, DRF will add a random string to"
						+ " the user-specified filenames to avoid multiple runs writing to the"
						+ " same file.");
		

		return parser;

	}
	
	
	public void findDarkRegions(ArgumentParser parser, String[] args){

		Namespace parsedArgs = null;
		try{
			parsedArgs = parser.parseArgs(args);
		} catch (ArgumentParserException e){
			parser.handleError(e);
			System.exit(1);
		}
		
		/*
		 * Input files
		 */
		String sam = parsedArgs.getString("SAM");
		String hgRef = parsedArgs.getString("HG_REF");

		/*
		 * Output files
		 */
		String lowDepthBed = parsedArgs.getString("LOW_COV_BED");
		String lowMapQBed = parsedArgs.getString("LOW_MAPQ_BED");
		String incBed = parsedArgs.getString("INC_BED");

		int minMapQMass = parsedArgs.getInt("MIN_MAPQ_MASS");
		int minRegionSize = parsedArgs.getInt("MIN_SIZE");
		int minDepth = parsedArgs.getInt("MIN_DEPTH");
		int mapQThresh = parsedArgs.getInt("MAPQ_THRESHOLD");
		boolean exclusive = parsedArgs.getBoolean("EXCLUSIVE");
		String stringency = parsedArgs.getString("STRINGENCY");
		
		List<String> intervalList = parsedArgs.getList("INTERVAL_LIST");
		ValidationStringency vs = null;
		
		if("strict".equalsIgnoreCase(stringency)){
			vs = ValidationStringency.STRICT;
		}
		else if("lenient".equalsIgnoreCase(stringency)){
			vs = ValidationStringency.LENIENT;
		}
		else if("silent".equalsIgnoreCase(stringency)){
			vs = ValidationStringency.SILENT;
		}
		
		try {
			
			/*
			 * DRF will write to .gz file. Add .gz to file names if not present.
			 */
			if(!lowDepthBed.startsWith("/dev/null")) {
				lowDepthBed = lowDepthBed.endsWith(".gz") ? lowDepthBed : lowDepthBed + ".gz";
			}
			if(!lowMapQBed.startsWith("/dev/null")) {
				lowMapQBed = lowMapQBed.endsWith(".gz") ? lowMapQBed : lowMapQBed + ".gz";
			}
			if(!incBed.startsWith("/dev/null")) {
				incBed = incBed.endsWith(".gz") ? incBed : incBed + ".gz";
			}
			
			File lowDepthBedFile = new File(lowDepthBed);
			File lowMapQBedFile = new File (lowMapQBed);
			File incBedFile = new File(incBed);
			
			/*
			 * If an interval list is specified, append random string to output
			 * file names.
			 */
			if(null != intervalList) {
				
				File[] newOutputFiles = DarkRegionFinderEngine.createUniqueOutputFileNames(lowDepthBed, lowMapQBed, incBed);
				
				lowDepthBedFile = newOutputFiles[0];
				lowMapQBedFile = newOutputFiles[1];
				incBedFile = newOutputFiles[2];
			}
			
			// Do your thing.
			DarkRegionFinder cgf = new DarkRegionFinder(new File(sam),
					lowDepthBedFile, lowMapQBedFile, incBedFile,
					new File(hgRef), mapQThresh, minMapQMass, minRegionSize, minDepth,
                    exclusive, vs, intervalList);

			cgf.startWalkingByLocus();

		} catch (FileNotFoundException e) {
			DarkRegionFinderEngine.printErrorUsageHelpAndExit(parser, logger, e);
		} catch (IOException e) {
			DarkRegionFinderEngine.printErrorAndExit(e);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.exit(1);
		}
	}

	
	/**
	 * Print the error. Then print the usage and help
	 * information and exit
	 * @param e
	 */
	private static void printErrorUsageHelpAndExit(ArgumentParser parser, Logger logger, Exception e){
		System.err.println("\nERROR: " + e.getMessage() + "\n");
//		logger.error(e.getMessage());
		DarkRegionFinderEngine.printUsageHelpAndExit(parser);
		System.exit(1);		
	}
	
	
	/**
	 * Print only the usage and help information and exit.
	 */
	private static void printUsageHelpAndExit(ArgumentParser parser){
		parser.printUsage();
		parser.printHelp();
		System.exit(1);		
	}
	
	private static void printErrorAndExit(Exception e){
		System.err.println("\nERROR: " + e.getMessage() + "\n");
		System.exit(1);		
	}
	
	/**
	 * Generate a random String to append to user-defined file names
	 * when -L (--interval-list) is used to prevent multiple DRF instances
	 * from writing to the same file.
	 * 
	 * Obtained from: https://stackoverflow.com/questions/20536566/creating-a-random-string-with-a-z-and-0-9-in-java
	 * 
	 * @param length: how long to make the random string
	 * @return a salt string of length 'length'
	 */
	private static String getSaltString(int length) {
		
		/*
		 * Define possible salt characters. Included numbers twice to increase
		 * their presence.
		 */
        String SALTCHARS = "a1bc2de3fg4hi5jk6lm7no8pq9rs0tuvwxyz1234567890";
        StringBuilder salt = new StringBuilder();
        Random rnd = new Random();
        while (salt.length() < length) { // length of the random string.
            int index = rnd.nextInt(SALTCHARS.length());
            salt.append(SALTCHARS.charAt(index));
        }
        String saltStr = salt.toString();
        return saltStr;
    }
	
	private static File[] createUniqueOutputFileNames(String lowDepthBed, String lowMapQBed, String incBed) throws IOException {
		String[] lowDepthFileAndExtension = DarkRegionFinderEngine.getFileAndExtension(lowDepthBed);
		String[] lowMAPQFileAndExtension = DarkRegionFinderEngine.getFileAndExtension(lowMapQBed);
		String[] incFileAndExtension = DarkRegionFinderEngine.getFileAndExtension(incBed);
			
		int saltLength = 5;
		String saltString, devNull = "/dev/null";
		File lowDepthBedFile, lowMapQBedFile, incBedFile;
		while(true) {
			saltString = DarkRegionFinderEngine.getSaltString(saltLength);
			

			lowDepthBedFile = lowDepthBed.startsWith(devNull) ?
					new File(lowDepthBed) :
						new File(lowDepthFileAndExtension[0] + ".salt_"
								+ saltString + lowDepthFileAndExtension[1]);

			lowMapQBedFile = lowMapQBed.startsWith(devNull) ?
					new File(lowMapQBed) :
						new File(lowMAPQFileAndExtension[0] + ".salt_"
								+ saltString + lowMAPQFileAndExtension[1]);

			incBedFile = incBed.startsWith(devNull) ? 
					new File(incBed) :
						new File(incFileAndExtension[0] + ".salt_"
								+ saltString + incFileAndExtension[1]);
			
			/*
			 * Verify we haven't already created a file with this exact name,
			 * including the salt string. If any of them exist, loop again.
			 */
			if( (lowDepthBedFile.isFile() && !lowDepthBedFile.getAbsolutePath().startsWith(devNull)) ||
					(lowMapQBedFile.isFile() && !lowMapQBedFile.getAbsolutePath().startsWith(devNull)) ||
					(incBedFile.isFile() && !incBedFile.getAbsolutePath().startsWith(devNull)) ) {
				continue;
			}
			
			/*
			 * We created a unique saltString for these files.
			 */
			break;
		}
		return new File[] {lowDepthBedFile, lowMapQBedFile, incBedFile};
	}
	
	/**
	 * Remove file extension from provided file string.
	 * 
	 * Obtained from (and modified): https://stackoverflow.com/questions/941272/how-do-i-trim-a-file-extension-from-a-string-in-java
	 * 
	 * @param fileString
	 * @return
	 * @throws IOException 
	 */
	private static String[] getFileAndExtension(String fileString) throws IOException {

	    String separator = System.getProperty("file.separator");
	    String path, filename;

	    // Remove the path upto the filename.
	    int lastSeparatorIndex = fileString.lastIndexOf(separator);
	    if (lastSeparatorIndex == -1) {
	    	path = "";
	        filename = fileString;
	    } else {
	        path = fileString.substring(0, lastSeparatorIndex + 1);
	        filename = fileString.substring(lastSeparatorIndex + 1);
	    }
	    
	    /*
	     * Catch case where user specifies a hidden file naem (beginning with '.'.
	     */
	    if(filename.startsWith(".")) {
	    	throw new IOException("The provided file name is a hidden file"
	    			+ " (starts with '.'). Please specify a non-hidden file"
	    			+ " name.");
	    }

	    // Return filename without extension and extension in a String[]
	    int extensionIndex = filename.lastIndexOf(".");
	    if (extensionIndex == -1) {
	        return new String[] {filename, ""};
	    }

	    return new String[] {path + filename.substring(0, extensionIndex),
	    		filename.substring(extensionIndex)};
	}

}
