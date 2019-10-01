package model;
import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

enum Strand{NO, YES, REVERSE};

public class Parameters 
{
	public static final String currentVersion = "1.0";
	
	// Input parameters
	public static String outputFolder = null;
	public static String sample1 = null;
	public static String sample2 = null;
	public static String crossTag = null;
	public static File inputVCFFile = null;
	public static File inputBAMFile = null;
	public static File inputGTFFile = null;
	public static long chunkSize = 1000000;
	public static Strand stranded = Strand.YES;
	
	// Computed variables
	public static long readMatchingCross = 0;
	public static long readOverlapSNPs = 0;
	public static long notIndels = 0;
	public static long unknownSNP = 0;
	public static long nbReads = 0;
	public static long noFeature = 0;
	public static long notUnique = 0;
	public static long ambiguous = 0;
	public static long mapped = 0;
	public static long unmapped = 0;
	public static long toolowAqual = 0;
	public static HashSet<String> uniqueGeneId =  null;
	public static DecimalFormat pcFormatter = new DecimalFormat("##.##");
	
	// Temporary
	public static HashMap<String, String> eQTLNAIVE = null;
	public static HashMap<String, String> eQTLTREATED = null;
	public static HashMap<String, ArrayList<SNP>> eQTLNAIVEPerGENE = null;
	public static HashMap<String, ArrayList<SNP>> eQTLTREATEDPerGENE = null;
	
	public static void load(String[] args) throws Exception
	{
		if(args.length == 0)
		{
			printHelp();
			System.exit(0);
		}
		for(int i = 0; i < args.length; i++) 
		{
			if(args[i].startsWith("-"))
			{
				switch(args[i])
				{
					case "-o":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						File f = new File(outputFolder);
						if(!f.exists()) 
						{
							System.out.println("Output folder does not exist. Creating it");
							f.mkdirs();
						}
						else if(!f.isDirectory()) throw new Exception(outputFolder + " is not a folder.");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						break;
					case "-s1":
						i++;
						try
						{
							if(args[i].trim().equals("")) throw new Exception("The '-s1' option should be followed by name of sample 1 from the cross (from VCF file). You entered " + args[i]);
							sample1 = args[i];
						}
						catch(Exception e)
						{
							System.err.println("The '-s1' option should be followed by name of sample 1 from the cross (from VCF file) " + e.getMessage() + ". You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-s2":
						i++;
						try
						{
							if(args[i].trim().equals("")) throw new Exception("The '-s2' option should be followed by name of sample 2 from the cross (from VCF file). You entered " + args[i]);
							sample2 = args[i];
						}
						catch(Exception e)
						{
							System.err.println("The '-s2' option should be followed by name of sample 2 from the cross (from VCF file) " + e.getMessage() + ". You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-cross":
						i++;
						try
						{
							if(args[i].trim().equals("")) throw new Exception("The '-cross' option should be followed by name of the cross (from BAM file @SM tag in Read Group Header). You entered " + args[i]);
							crossTag = args[i];
						}
						catch(Exception e)
						{
							System.err.println("The '-cross' option should be followed by name of the cross (from BAM file @SM tag in Read Group Header) " + e.getMessage() + ". You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-s":
						i++;
						switch(args[i])
						{
						case "no":
							Parameters.stranded = Strand.NO;
							break;
						case "yes":
							Parameters.stranded = Strand.YES;
							break;
						case "reverse":
							Parameters.stranded = Strand.REVERSE;
							break;
						default:
							System.err.println("The '-s' option should be followed by one of the following parameters: [no, yes, reverse].");
							System.exit(-1);
						}
						break;
					case "-gtf":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) throw new Exception("No file at path " + args[i]);
							if(!c.isFile()) throw new Exception(args[i] + " is not a file");
							inputGTFFile = c;
						}
						catch(Exception e)
						{
							System.err.println("The '-gtf' option should be followed by GTF file path. " + e.getMessage() + ". You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-vcf":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) throw new Exception("No file at path " + args[i]);
							if(!c.isFile()) throw new Exception(args[i] + " is not a file");
							inputVCFFile = c;
						}
						catch(Exception e)
						{
							System.err.println("The '-vcf' option should be followed by VCF file path. " + e.getMessage() + ". You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-bam":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) throw new Exception("No file at path " + args[i]);
							if(!c.isFile()) throw new Exception(args[i] + " is not a file");
							inputBAMFile = c;
						}
						catch(Exception e)
						{
							System.err.println("The '-bam' option should be followed by BAM file path. " + e.getMessage() + ". You entered " + args[i]);
							System.exit(-1);
						}
						break;
					default:
						System.err.println("Unused argument: " + args[i]);
				}
			}
		}
		if(inputVCFFile == null)
		{
			System.err.println("Please use '-vcf' option to specify VCF file");
			System.exit(-1);
		}
		if(inputBAMFile == null)
		{
			System.err.println("Please use '-bam' option to specify BAM file");
			System.exit(-1);
		}
		if(sample1 == null)
		{
			System.err.println("Please use '-s1' option to specify name of sample1 from the cross.");
			System.exit(-1);
		}
		if(sample2 == null)
		{
			System.err.println("Please use '-s2' option to specify name of sample1 from the cross.");
			System.exit(-1);
		}
		if(crossTag == null)
		{
			System.err.println("Please use '-cross' option to specify tag name of the cross in BAM file (from Read Group @SM tag in header).");
			System.exit(-1);
		}
		if(inputGTFFile == null)
		{
			System.err.println("Please use '-gtf' option to specify the path of the GTF file");
			System.exit(-1);
		}
		if(outputFolder == null)
		{
			String path = inputVCFFile.getAbsolutePath();
			path = path.replaceAll("\\\\", "/");
			path = path.substring(0, path.lastIndexOf("/"));
			System.out.println("No output Folder is specified, using default one: \""+path+"\". You can specify an output path by using the '-o' option.");
			outputFolder = path;
		}
		System.out.println("Stranded = " + Parameters.stranded);
		outputFolder = outputFolder.replaceAll("\\\\", "/");
		if(!outputFolder.endsWith("/")) outputFolder+= "/";
		new File(outputFolder).mkdirs();
	}
	
	private static void printHelp()
	{
		System.out.println("-vcf %s \t\t[Required] Path of VCF file.");
		System.out.println("-bam %s \t\t[Required] Path of aligned BAM file [do not need to be sorted or indexed].");
		System.out.println("-gtf %s \t[Required] Path of GTF file.");
		System.out.println("-s1 %s \t\t[Required] Name of sample 1 from the cross.");
		System.out.println("-s2 %s \t\t[Required] Name of sample 2 from the cross.");
		System.out.println("-cross %s\t\t[Required] Tag name of the cross in BAM file (from Read Group @SM tag in header).");
		System.out.println("-o %s \t\tOutput folder [default = folder of VCF file]");
		System.out.println("-s %s \t\tDo you want to count only reads falling on same strand than gene? [no, yes, reverse] [default = yes, since BRB-seq is stranded protocol].");
	}
}
