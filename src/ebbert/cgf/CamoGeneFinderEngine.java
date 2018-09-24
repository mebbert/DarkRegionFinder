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
import net.sourceforge.argparse4j.inf.ArgumentGroup;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

/**
 * @author markebbert
 *
 */
public class CamoGeneFinderEngine {

	private static Logger logger = Logger.getLogger(CamoGeneFinderEngine.class);
	
	public CamoGeneFinderEngine() {
		return;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		BasicConfigurator.configure();
		CamoGeneFinderEngine cgfe = new CamoGeneFinderEngine();
		ArgumentParser parser = cgfe.init(args);

		cgfe.findCamoGenes(parser, args);
	}
	
	/**
	 * Init the parser options
	 * 
	 * @param args
	 */
	private ArgumentParser init(String[] args){
		ArgumentParser parser = ArgumentParsers.newArgumentParser("CamoGeneFinder");
		parser.description("The Camo Gene Finder will identify"
				+ " regions where the genome is either 'camouflaged', 'dark', or"
				+ " 'incomplete'. Camouflaged regions are those where the aligner"
				+ " assigns a low MAPQ (e.g., MAPQ == 0 in BWA) because the read"
				+ " aligns to multiple locations equally well. Dark regions are"
				+ " those where depth is less than the defined '--min-depth'"
				+ " threshold. Incomplete regions are simply those where"
				+ " nucleotides are unknown (i.e., 'N' or 'n'). Incomplete"
				+ " regions are considered separate from dark regions (i.e.,"
				+ " they are not included as a dark region, and vice-versa).");
		parser.defaultHelp(true);
		
		ArgumentGroup svcOptions = parser.addArgumentGroup("Camo gene finder arguments");
		ArgumentGroup ioOptions = parser.addArgumentGroup("input/output arguments");
			
		/* Setup SVC options */
		svcOptions
				.addArgument("-s", "--min-region-size")
				.dest("MIN_SIZE")
				.metavar("SIZE")
				.setDefault(1)
				.type(Integer.class)
				.help("The minimum size camo region to consider. Details will"
						+ " still be written at the individual base level, but"
						+ " only for regions that meet this size requirement."
						+ " Regions can be merged using bedtools, if desired.");

		svcOptions
				.addArgument("-t", "--mapq-threshold")
				.dest("MAPQ_THRESHOLD")
				.metavar("THRESH")
				.setDefault(0)
				.type(Integer.class)
				.help("The MAPQ threshold (≤) at which a read is \'inadequately\'"
						+ " aligned and considered 'camouflaged'. Generally"
						+ " recommended to use MAPQ == 0"
						+ " when looking for camo-genes in standard illumina"
						+ " data, when aligned by BWA. Can also be used to"
						+ " identify regions where"
						+ " mapping quality is too low to trust for variant"
						+ " calling (e.g., MAPQ ≤ 19).");
	
		svcOptions
				.addArgument("-m", "--min-mapq-mass")
				.dest("MIN_MASS")
				.metavar("MASS")
				.setDefault(90)
				.type(Integer.class)
				.help("The minimum percentage (≥) of reads below the MAPQ threshold"
						+ " for the locus to"
						+ " be considered camoflauged. Loci where the percentage"
						+ " is below this threshold will not be printed.");
		
		svcOptions
				.addArgument("-d", "--min-camo-depth")
				.dest("MIN_CAMO_DEPTH")
				.metavar("MIN_CAMO_DEPTH")
				.setDefault(5)
				.type(Integer.class)
				.help("The minimum depth (≥) required to identify a camo"
						+ " region.");

		svcOptions
				.addArgument("-r", "--dark-depth")
				.dest("DARK_DEPTH")
				.metavar("DARK_DEPTH")
				.setDefault(20)
				.type(Integer.class)
				.help("The depth (≤) at which a region is considered 'dark' (i.e.,"
						+ " it's not deep enough to reliably call heterozygous"
						+ " mutations in the organism). Regions where"
						+ " MIN_CAMO_DEPTH < DEPTH < DARK_DEPTH, and that are"
						+ " camouflaged, will be reported"
						+ " in both the camo and dark bed files.");

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
						+ "by 'samtools faidx'.");

		ioOptions
				.addArgument("-c", "--camo-bed-output")
				.dest("CAMO_BED")
				.type(String.class)
				.required(true)
				.setDefault("camo.bed")
				.help("The output BED file for camouflaged regions. Camouflaged"
						+ " regions are those where reads map equally well to"
						+ " multiple regions, but not so many the aligner gives"
						+ " up (e.g., 10000 for BWA MEM)");
		
		ioOptions
				.addArgument("-a", "--dark-bed-output")
				.dest("DARK_BED")
				.type(String.class)
				.required(true)
				.setDefault("dark.bed")
				.help("The output BED file for dark regions. Dark regions are"
						+ " those with low coverage (depth < --min-depth).");

		ioOptions
				.addArgument("-n", "--incomplete-bed-output")
				.dest("INC_BED")
				.type(String.class)
				.required(true)
				.setDefault("incomplete.bed")
				.help("The output BED file for incomplete regions. Incomplete"
						+ " regions are those where the bases are unknown"
						+ " (i.e., 'N' or 'n').");
		

		return parser;

	}
	public void findCamoGenes(ArgumentParser parser, String[] args){

		Namespace parsedArgs = null;
		try{
			parsedArgs = parser.parseArgs(args);
		} catch (ArgumentParserException e){
			parser.handleError(e);
			System.exit(1);
		}
		
		String sam = parsedArgs.getString("SAM");
		String camoBed = parsedArgs.getString("CAMO_BED");
		String darkBed = parsedArgs.getString("DARK_BED");
		String incBed = parsedArgs.getString("INC_BED");
		String hgRef = parsedArgs.getString("HG_REF");

		int minMass = parsedArgs.getInt("MIN_MASS");
		int minRegionSize = parsedArgs.getInt("MIN_SIZE");
		int minCamoDepth = parsedArgs.getInt("MIN_CAMO_DEPTH");
		int darkDepth = parsedArgs.getInt("DARK_DEPTH");
		int mapQThresh = parsedArgs.getInt("MAPQ_THRESHOLD");
		String stringency = parsedArgs.getString("STRINGENCY");
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
			CamoGeneFinder cgf = new CamoGeneFinder(new File(sam),
					new File(camoBed), new File(darkBed), new File(incBed),
					new File(hgRef), mapQThresh, minMass, minRegionSize,
					minCamoDepth, darkDepth, vs);

			cgf.startWalkingByLocus();

		} catch (FileNotFoundException e) {
			CamoGeneFinderEngine.printErrorUsageHelpAndExit(parser, logger, e);
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
		CamoGeneFinderEngine.printUsageHelpAndExit(parser);
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
