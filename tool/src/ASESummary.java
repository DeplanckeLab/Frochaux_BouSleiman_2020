import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.zip.GZIPInputStream;

import model.EfficientHashMap;
import model.SNP;
import tools.Utils;

/**
 * @author Vincent Gardeux
 * @see vincent.gardeux@epfl.ch
 *
 */
public class ASESummary
{	
	public static HashMap<String, Integer> sampleIndex = new HashMap<String, Integer>();
	public static EfficientHashMap variantIndex = new EfficientHashMap();
	public static SNP[] variantDesc = new SNP[1278439];
	public static int[][] variants = new int[1278439][33]; // I know the number of SNPs/Samples in my VCF. CHEAT!
	public static int[][] variants_ref = new int[1278439][33]; // I know the number of SNPs/Samples in my VCF. CHEAT!
	public static int[][] variants_alt = new int[1278439][33]; // I know the number of SNPs/Samples in my VCF. CHEAT!
	public static int classicCol = 0; // How much to shift the table
	
	public static void main(String[] args) throws Exception
	{
		System.out.println("ASESummary\n");
		
		//Read VCF file first
		System.out.println("Reading VCF...");
		//readVCF("M:/vgardeux/ChIPseq_2018_10_24_Riccardo_VCM/Genotype/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.filtered.recode.vcf");
		readVCF("M:/vgardeux/ChIPseq_2018_11_23_Sebastian_VCM/Genotype/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.filtered.NA12878.vcf");
		for(int i = 0; i < variants.length; i++) for(int j = 0; j < variants[i].length; j++) {
			variants_alt[i][j] = variants[i][j];
			variants_ref[i][j] = variants[i][j];
		}
		
		// Read all ASE files and summarize
		System.out.println("Reading ASE files...");
		 // All ASEs
		ArrayList<File> aseFiles = listASEFiles("M:/vgardeux/ChIPseq_2018_11_23_Sebastian_VCM/ASE/ENCODE_MEF2A/");
		// Only NA11993
		//ArrayList<File> aseFiles = new ArrayList<File>(); //listASEFiles("M:/vgardeux/ChIPseq_2018_10_24_Riccardo_VCM/ASE/");
		//aseFiles.add(new File("M:/vgardeux/ChIPseq_2018_11_23_Sebastian_VCM/ASE/ENCODE_MEF2A/NA12878.2.ASE.table"));
		//aseFiles.add(new File("M:/vgardeux/ChIPseq_2018_11_23_Sebastian_VCM/ASE/NA11993.ASE.table"));
		// Only NA12489
		//aseFiles.add(new File("M:/vgardeux/ChIPseq_2018_10_24_Riccardo_VCM/ASE/NA12489.ASE.table"));
		// Only NA12286
		//aseFiles.add(new File("M:/vgardeux/ChIPseq_2018_10_24_Riccardo_VCM/ASE/NA12286.ASE.table"));
		// Only NA11830
		//aseFiles.add(new File("M:/vgardeux/ChIPseq_2018_10_24_Riccardo_VCM/ASE/NA11830.ASE.table"));
		// Only NA12873
		//aseFiles.add(new File("M:/vgardeux/ChIPseq_2018_10_24_Riccardo_VCM/ASE/NA12873.ASE.table"));
		// Only NA12761
		//aseFiles.add(new File("M:/vgardeux/ChIPseq_2018_10_24_Riccardo_VCM/ASE/NA12761.ASE.table"));
		// Only NA12154
		//aseFiles.add(new File("M:/vgardeux/ChIPseq_2018_10_24_Riccardo_VCM/ASE/NA12154.ASE.table"));
		// Only NA12283
		//aseFiles.add(new File("M:/vgardeux/ChIPseq_2018_10_24_Riccardo_VCM/ASE/NA12283.ASE.table"));
		System.out.println(aseFiles.size() + " ASE files were found in folder.");
		for(File f:aseFiles) addToVCF(f);
		
		// Writing results
		System.out.println("Writing results...");
		BufferedWriter bw = new BufferedWriter(new FileWriter("M:/vgardeux/ChIPseq_2018_11_23_Sebastian_VCM/ASE/ENCODE_MEF2A/data.ase.all.txt"));
		bw.write("CHROM\tPOS\tID\tREF\tALT\tref_count\talt_count\tp_val\n");
		for(int i = 0; i < variants.length; i++) 
		{
			SNP snp = variantDesc[i];
			if(snp != null && variantIndex.get(snp.chr, snp.rsId) != null)
			{
				bw.write(snp.chr + "\t" + snp.loc + "\t" + snp.rsId + "\t" + snp.refAllele + "\t" + snp.altAllele);
				int alt = 0;
				int ref = 0;
				for(int j = 0; j < variants[i].length; j++) 
				{
					if(variants_alt[i][j] != -1) alt += variants_alt[i][j];
					if(variants_ref[i][j] != -1) ref += variants_ref[i][j];
				}
				bw.write("\t" + ref + "\t" + alt + "\t" + Utils.binomialTest(ref, alt) + "\n");
			}
		}
		bw.close();
		System.out.println("DONE");
	}
	
	public static void addToVCF(File f) throws Exception
	{
		String name = f.getName().substring(0, 7);
		if(sampleIndex.get(name) == null) {System.out.println("Sample " + name + " is not in the VCF..."); return;}
		System.out.println("Reading " + name + " ASE file and adding to VCF...");
		long start = System.currentTimeMillis();
		int count = 0;
		BufferedReader br = new BufferedReader(new FileReader(f));
		HashMap<String, Integer> index = new HashMap<>();
		String[] header = br.readLine().split("\t");
		for(int i = 0; i < header.length; i++) index.put(header[i], i);
		String line = br.readLine();
        while(line != null)
        {
            String[] tokens = line.split("\t");
            String variantID = tokens[index.get("variantID")];
            int refCount = Integer.parseInt(tokens[index.get("refCount")]);
            int altCount = Integer.parseInt(tokens[index.get("altCount")]);
            String chr = tokens[index.get("contig")];
            
            Integer n = variantIndex.get(chr, variantID);
            Integer p = sampleIndex.get(name);
            if(n != null) // Can be null if it was deleted because duplicated
            {            
	            SNP snp = variantDesc[n];
	            if(snp.loc != Integer.parseInt(tokens[index.get("position")]))
	            {
	        		System.out.println("WEIRD POS");
	        		System.exit(-1);
	            }
	            if(variants[n][p - classicCol] != -1)
	            {
	            	if(variants[n][p - classicCol] != 0) // Only possible if adding twice the same guy
	            	{
	            		variants_ref[n][p - classicCol] += refCount;
	            		variants_alt[n][p - classicCol] += altCount;
	            		variants[n][p - classicCol] = variants_ref[n][p - classicCol] + variants_alt[n][p - classicCol];
	            	}
	            	else
	            	{
	            		variants_ref[n][p - classicCol] = refCount;
	            		variants_alt[n][p - classicCol] = altCount;
	            		variants[n][p - classicCol] = variants_ref[n][p - classicCol] + variants_alt[n][p - classicCol];
	            	}
	            }
            }
            
            line = br.readLine(); count++;
        }
		br.close();
        System.out.println(count + " ASE variants analyzed [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
	}
		
	public static ArrayList<File> listASEFiles(String path)
	{
		ArrayList<File> results = new ArrayList<File>();
		File[] files = new File(path).listFiles();
		for(File file:files) if(file.isFile() && file.getName().endsWith(".table")) results.add(file);
		return results;
	}
	
	public static void readVCF(String filenameVCF) throws Exception
	{
		long start = System.currentTimeMillis();
		BufferedReader br;
		if(filenameVCF.endsWith("vcf.gz")) br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(filenameVCF))));
		else br = new BufferedReader(new FileReader(filenameVCF));
		
        String[] knownCols = {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"};

        String line = br.readLine();
        while(line.startsWith("##")) line = br.readLine();
        if(line.startsWith("#CHROM")) // Header
        {
    		String[] header = line.substring(1).split("\t"); // remove the # and split
    		for (int i = 0; i < header.length; i++) 
    		{
    			if(Utils.contains(knownCols,header[i])) classicCol++;
    			else sampleIndex.put(header[i], i);
    		}
        }
        line = br.readLine();
        System.out.println(classicCol + " classic columns found.");
        System.out.print(sampleIndex.size() + " samples found: ");
        for(String s:sampleIndex.keySet()) System.out.print(s + " ");
        System.out.println();
        	
        int count = 0;
        HashSet<SNP> toDelete = new HashSet<SNP>(); // Handling duplicates
        
        while(line != null)
        {
            String[] tokens = line.split("\t");

    		SNP s = new SNP();
    		s.chr = tokens[0]; // CHROM is 0
    		s.loc = Long.parseLong(tokens[1]); // POS is 1
    		s.rsId = tokens[2]; // ID is 2
        	s.refAllele = tokens[3]; // REF is 3
        	s.altAllele = tokens[4]; // ALT is 4
        	
        	if(variantIndex.get(s.chr, s.rsId) != null) toDelete.add(s); // Duplicate, marked for deletion
        	else variantIndex.put(s.chr, s.rsId, count);
        	variantDesc[count] = s;

        	for(String nameSampleToExtract:sampleIndex.keySet())
        	{
        		int col = sampleIndex.get(nameSampleToExtract);
            	String geno = tokens[col];
            	if(!geno.equals("./.") && !geno.equals(".|."))
            	{
            		// Which genotype in Ind file
	            	if(geno.equals("1/0") || geno.equals("0/1") || geno.equals("1|0") || geno.equals("0|1")) variants[count][col - classicCol] = 0; // Heterozygous, I want that
	            	else variants[count][col - classicCol] = -1;// Homozygous or else... I don't want that
            	}
        	}
            
            line = br.readLine(); count++;
            if(count % 100000 == 0) System.out.println(count + " variants analyzed [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
        }
        
        for(SNP s:toDelete) { System.out.print(s.rsId + " "); variantIndex.remove(s.chr, s.rsId); }
        System.out.println(" were deleted.");
        
        br.close();
        System.out.println(count + " variants analyzed [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
        System.out.println("The VCF file contained " + toDelete.size() + " duplicated variants that were removed.");
        System.out.println("The VCF file contains " + variantIndex.size() + " unique rsIDs");
	}
}
