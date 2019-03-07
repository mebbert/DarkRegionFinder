/**
 * 
 */
package ebbert.cgf;

import java.io.File;
import java.io.FileNotFoundException;

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
		DarkRegionFinderEngine drfe = new DarkRegionFinderEngine();
		ArgumentParser parser = drfe.init(args);

		drfe.findDarkGenes(parser, args);
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
		
		ArgumentGroup svcOptions = parser.addArgumentGroup("Dark Region Finder arguments");
		ArgumentGroup ioOptions = parser.addArgumentGroup("input/output arguments");
			
		/* Setup SVC options */
		svcOptions
				.addArgument("-s", "--min-region-size")
				.dest("MIN_SIZE")
				.metavar("SIZE")
				.setDefault(1)
				.type(Integer.class)
				.help("The minimum size dark region to consider. Details will"
						+ " still be written at the individual base level, but"
						+ " only for regions that meet this size requirement."
						+ " Regions can be merged using bedtools, if desired.");

		svcOptions
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

		svcOptions
				.addArgument("-m", "--min-mapq-mass")
				.dest("MIN_MAPQ_MASS")
				.metavar("MAPQ_MASS")
				.setDefault(90)
				.type(Integer.class)
				.help("The minimum percentage (≥) of reads below the"
						+ " --mapq-threshold for the locus to be considered dark."
						+ " Dark loci where the percentage is below this threshold"
						+ " will still be reported if the depth ≤ --min-depth");

		svcOptions
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

		svcOptions
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

		svcOptions
				.addArgument("-v", "--validation-stringency")
				.dest("STRINGENCY")
				.setDefault("STRICT")
				.choices("STRICT", "LENIENT", "SILENT")
				.type(String.class)
				.help("The validation stringency when parsing a SAM/BAM"
						+ " file. 'STRICT' will throw errors if something "
						+ " is amiss, 'LENIENT' will give warnings but continue,"
						+ " and 'SILENT' will continue AND keep our mouth shut.");

		svcOptions
				.addArgument("-j", "--threads")
				.dest("THREADS")
				.type(Integer.class)
				.setDefault(1)
				.help("Number of threads to use for multithreaded parallelism. Default is 1 thread, i.e. no parallelism.");

		svcOptions
				.addArgument("-k", "--tmp-dir")
				.dest("TMPDIR")
				.type(String.class)
				.setDefault("./")
				.help("Where to write temporary files created during computation."
						+ " Temporary files are only created if parallelism is enabled (i.e. --threads ≥ 2)"
						+ " Defaults to current working directory.");

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
						+ "by 'samtools faidx' and have a Picard sequence"
						+ " dictionary.");

		ioOptions
				.addArgument("-c", "--low-coverage-bed-output")
				.dest("LOW_COV_BED")
				.type(String.class)
				.setDefault("low_coverage.dark.bed")
				.help("The output BED file for low-coverage dark regions. Low coverage dark"
						+ " regions are those where depth ≤ --min_depth");

		ioOptions
				.addArgument("-a", "--low-mapq-bed-output")
				.dest("LOW_MAPQ_BED")
				.type(String.class)
				.setDefault("low_mapq.dark.bed")
				.help("The output BED file for low MAPQ dark regions. Low MAPQ regions are"
						+ " those with low mapping quality (loci with"
						+ " MAPQ ≤ --mapq_threshold and MAPQ mass ≥ --min_mapq_mass).");

		ioOptions
				.addArgument("-n", "--incomplete-bed-output")
				.dest("INC_BED")
				.type(String.class)
				.setDefault("incomplete.bed")
				.help("The output BED file for incomplete regions. Incomplete"
						+ " regions are those where the bases are unknown"
						+ " (i.e., 'N' or 'n').");


		return parser;

	}
	public void findDarkGenes(ArgumentParser parser, String[] args){

		Namespace parsedArgs = null;
		try{
			parsedArgs = parser.parseArgs(args);
		} catch (ArgumentParserException e){
			parser.handleError(e);
			System.exit(1);
		}
		
		String sam = parsedArgs.getString("SAM");
		String lowDepthBed = parsedArgs.getString("LOW_COV_BED");
		String lowMapQBed = parsedArgs.getString("LOW_MAPQ_BED");
		String incBed = parsedArgs.getString("INC_BED");
		String hgRef = parsedArgs.getString("HG_REF");

		int minMapQMass = parsedArgs.getInt("MIN_MAPQ_MASS");
		int minRegionSize = parsedArgs.getInt("MIN_SIZE");
		int minDepth = parsedArgs.getInt("MIN_DEPTH");
		int mapQThresh = parsedArgs.getInt("MAPQ_THRESHOLD");
		boolean exclusive = parsedArgs.getBoolean("EXCLUSIVE");
		String stringency = parsedArgs.getString("STRINGENCY");
		String tmpdir = parsedArgs.getString("TMPDIR");
		int threads = parsedArgs.getInt("THREADS");
//		boolean ignoreLowCovRegions = parsedArgs.getBoolean("IGNORE_LOW_COV");
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
			// Do your thing.

			IDarkRegionFinder drf = null;
			if (threads == 1) {
				drf = new DarkRegionFinder(new File(sam),
						new File(lowDepthBed), new File(lowMapQBed), new File(incBed),
						new File(hgRef), mapQThresh, minMapQMass, minRegionSize, minDepth,
						exclusive, vs);
			} else if (threads > 1) {
				drf = new ThreadedDarkRegionFinder(new File(sam),
						new File(lowDepthBed), new File(lowMapQBed), new File(incBed),
						new File(hgRef), mapQThresh, minMapQMass, minRegionSize, minDepth,
						exclusive, vs, tmpdir, threads);
			}

			drf.startWalkingByLocus();

		} catch (FileNotFoundException e) {
			DarkRegionFinderEngine.printErrorUsageHelpAndExit(parser, logger, e);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
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
	}
	
	
	/**
	 * Print only the usage and help information and exit.
	 */
	private static void printUsageHelpAndExit(ArgumentParser parser){
		parser.printUsage();
		parser.printHelp();
		System.exit(1);		
	}

}
