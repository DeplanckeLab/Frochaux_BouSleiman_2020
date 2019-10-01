package model;

import java.io.BufferedReader;
import java.util.HashMap;
import java.util.HashSet;

import tools.Utils;

public class VCF
{
	private static String[] knownCols = {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"};
	private static HashSet<Character> authorized = new HashSet<Character>(); // Check if the variants are correctly formed
	
	static {
        authorized.add('a'); authorized.add('c'); authorized.add('g'); authorized.add('t');
        authorized.add('A'); authorized.add('C'); authorized.add('G'); authorized.add('T');
        authorized.add('N'); authorized.add('.');
	}

	public static void readVCF() throws Exception
	{
		System.out.println("\nReading VCF file provided: " + Parameters.inputVCFFile.getAbsolutePath());
		BufferedReader br = Utils.readVCF(Parameters.inputVCFFile);
        
        // Reading Header
		int l = 0;
		HashMap<String, Integer> indexes = new HashMap<>();
        String line = br.readLine(); l++;
        while(!line.startsWith("#CHROM")) {line = br.readLine(); l++;}
    	String[] header = line.substring(1).split("\t"); // remove the # and split
    	for (int i = 0; i < header.length; i++) indexes.put(header[i], i);
    	       
        System.out.println((indexes.size() - knownCols.length) + " samples found in VCF");
        if(indexes.get(Parameters.sample1) == null) System.err.println(Parameters.sample1 + " was not found in the VCF...");
        if(indexes.get(Parameters.sample2) == null) System.err.println(Parameters.sample2 + " was not found in the VCF...");
        if(indexes.get(Parameters.sample1) == null || indexes.get(Parameters.sample2) == null) System.exit(-1);
        System.out.println(Parameters.sample1 + " is at position " + indexes.get(Parameters.sample1));
        System.out.println(Parameters.sample2 + " is at position " + indexes.get(Parameters.sample2));
        
        line = br.readLine(); l++;
        int count = 0, het = 0, inGenes = 0;
        while(line != null)
        {
        	count++;
            String[] tokens = line.split("\t");
        	SNP s = new SNP();
        	s.chr = tokens[indexes.get("CHROM")];
        	s.loc = Long.parseLong(tokens[indexes.get("POS")]);
            s.refAllele = tokens[indexes.get("REF")]; // May be different than in REF
            s.altAllele = tokens[indexes.get("ALT")]; // May be different than in REF
            s.rsId = tokens[indexes.get("ID")];

            if(checkIntegrity(s.refAllele, l) && checkIntegrity(s.altAllele, l))
            {
	           	String geno = tokens[indexes.get(Parameters.sample1)];
	           	int pos = geno.indexOf(":");
	           	if(pos != - 1) geno = geno.substring(0, geno.indexOf(":"));
	           	int genotype1 = Utils.getGeno(geno);
	           	if(genotype1 == 0 || genotype1 == 2) // If Sample1 is homozygous
	           	{
	           		geno = tokens[indexes.get(Parameters.sample2)];
	           		pos = geno.indexOf(":");
	           		if(pos != - 1) geno = geno.substring(0, geno.indexOf(":"));
	           		int genotype2 = Utils.getGeno(geno);
	           		if((genotype1 == 0 && genotype2 == 2) || (genotype1 == 2 && genotype2 == 0)) 
	           		{
	           			het++;
	           			HashSet<IntervalLabelled> genes = GTF.findOverlappingGenes(s.chr, (int)s.loc, (int)(s.loc + 1));
	           			for(IntervalLabelled i:genes) 
	           			{
	           				if(s.refAllele.length() == 1 && s.altAllele.length() == 1) //TODO consider indels
	           				{
	           					Parameters.notIndels++;
	           					if(genotype1 == 0) s.isRef = 1;
	           					else s.isRef = 2;
	           					s.isSNP = true;
	           					i.variants.add(s);
	           				}
	           				inGenes++;
	           			}
	           		}
	           	}
            }
            line = br.readLine(); l++;
        }
        
        br.close();
        System.out.println("The VCF File contained " + count + " variants (" + het + " of them are heterozygous in the cross " + Parameters.sample1 + " x "+ Parameters.sample2 + " and " + inGenes + " of them overlap known gene exons" + ")");
        System.out.println(Parameters.notIndels + " variants are SNPs (not INDELs) and thus will be used in the rest of the pipeline");
	}
		
	private static boolean checkIntegrity(String allele, int l)
	{
		String[] tokens = allele.split(",");
		HashSet<String> set = new HashSet<String>();
		for(String all:tokens)
		{
			if(set.contains(all))
			{
				System.out.println("[IGNORED] Malformed VCF at l. "+l+": Duplicate allele added to VariantContext: " + all); 
				return false; 
			}
			set.add(all);
			if(all.equals("")) 
			{ 
				System.out.println("[IGNORED] Malformed VCF at l. "+l+": empty alleles are not permitted in VCF records"); 
				return false; 
			}
	        for(int i = 0; i < all.length(); i++) 
	        {
	        	if(!authorized.contains(all.charAt(i))) 
	        	{ 
	        		System.out.println("[IGNORED] Malformed VCF at l. "+l+": unparsable vcf record with allele "+all.charAt(i)); 
	        		return false;
	        	}
	        }
		}
		return true;
	}
}
