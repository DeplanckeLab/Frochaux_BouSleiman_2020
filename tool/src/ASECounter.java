import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import htsjdk.tribble.index.interval.Interval;
import htsjdk.tribble.index.interval.IntervalTree;
import model.BAM;
import model.GTF;
import model.IntervalLabelled;
import model.Parameters;
import model.SNP;
import model.VCF;
import tools.Utils;

/**
 * @author Vincent Gardeux
 * @see vincent.gardeux@epfl.ch
 *
 */
public class ASECounter
{	
	public static void main(String[] args) throws Exception
	{
		System.out.println("ASECounter " + Parameters.currentVersion + "\n");
		Parameters.load(args);
		
		System.out.println("\nReading eQTL data...");
		Parameters.eQTLNAIVEPerGENE = new HashMap<String, ArrayList<SNP>>();
		Parameters.eQTLNAIVE = new HashMap<String, String>();
		BufferedReader br = new BufferedReader(new FileReader(Parameters.outputFolder + "MatrixeQTL.Naive.output.cis"));
		br.readLine(); // Skip first line
		int nbEQTLs = 0, nbLines = 0;
		String line = br.readLine();
		while(line != null)
		{
			nbLines++;
			String[] tokens = line.split("\t");
			double fdr = Double.parseDouble(tokens[5]);
			if(fdr <= 0.05)
			{
				nbEQTLs++;
				String[] snpData = tokens[0].split("\\.");
				SNP snp = new SNP();
				snp.chr = snpData[0];
				snp.loc = Long.parseLong(snpData[1]);
				ArrayList<SNP> list = Parameters.eQTLNAIVEPerGENE.get(tokens[1]);
				if(list == null) list = new ArrayList<>();
				list.add(snp);
				Parameters.eQTLNAIVEPerGENE.put(tokens[1], list);
				Parameters.eQTLNAIVE.put(snp.chr+":"+snp.loc, tokens[1]);
			}
			line = br.readLine();
		}
		br.close();
		System.out.println(Parameters.eQTLNAIVE.size());
		System.out.println("Read " + nbLines + " lines\n" + nbEQTLs + " NAIVE eQTLs (FDR<5%, " + Parameters.eQTLNAIVEPerGENE.size() + " unique genes)");
		
		Parameters.eQTLTREATEDPerGENE = new HashMap<String, ArrayList<SNP>>();
		Parameters.eQTLTREATED = new HashMap<String, String>();
		br = new BufferedReader(new FileReader(Parameters.outputFolder + "MatrixeQTL.Treated.output.cis"));
		br.readLine(); // Skip first line
		line = br.readLine();
		nbEQTLs = 0; nbLines = 0;
		while(line != null)
		{
			nbLines++;
			String[] tokens = line.split("\t");
			double fdr = Double.parseDouble(tokens[5]);
			if(fdr <= 0.05)
			{
				nbEQTLs++;
				String[] snpData = tokens[0].split("\\.");
				SNP snp = new SNP();
				snp.chr = snpData[0];
				snp.loc = Long.parseLong(snpData[1]);
				ArrayList<SNP> list = Parameters.eQTLTREATEDPerGENE.get(tokens[1]);
				if(list == null) list = new ArrayList<>();
				list.add(snp);
				Parameters.eQTLTREATEDPerGENE.put(tokens[1], list);
				Parameters.eQTLTREATED.put(snp.chr+":"+snp.loc, tokens[1]);
			}
			line = br.readLine();
		}
		br.close();
		System.out.println(Parameters.eQTLTREATED.size());
		System.out.println("Read " + nbLines + " lines\n" + nbEQTLs + " TREATED eQTLs (FDR<5%, " + Parameters.eQTLTREATEDPerGENE.size() + " unique genes)");
		
		GTF.readGTF();
		VCF.readVCF();
		BAM.readBAM();
		GTF.writeOutputPerGene();
		GTF.writeOutputPerSNP();
		
		HashSet<String> processed = new HashSet<String>();
		BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + Parameters.crossTag + "_" + Parameters.sample1 + "_"+ Parameters.sample2 + "_per_eqtl_treated.txt"));
		bw.write("gene\tsnp_chr\tsnp_pos\tsnp_ref\tsnp_alt\t"+Parameters.sample1+"\t"+Parameters.sample2+"\tbinomial_p\t"+Parameters.sample1+"_reads\t"+Parameters.sample2+"_reads\n");
		for(String chr:GTF.forest.keySet())
		{
			IntervalTree chrom = (IntervalTree)GTF.forest.get(chr);
			for(Interval i:chrom.getIntervals())
			{
				IntervalLabelled exon = (IntervalLabelled)i;
				ArrayList<SNP> eqtls = Parameters.eQTLTREATEDPerGENE.get(exon.gene);
				if(eqtls != null)
				{
					for(SNP eqtl:eqtls)
					{
						// Check if this eQTL is in our data
						for(SNP snp:exon.variants)
						{
							if(snp.chr.equals(eqtl.chr) && snp.loc == eqtl.loc && !processed.contains(exon.gene+":"+snp.chr+":"+snp.loc)) // It IS!
							{
								processed.add(exon.gene+":"+snp.chr+":"+snp.loc);
								bw.write(exon.gene + "\t" + snp.chr + "\t" + snp.loc + "\t" + snp.refAllele + "\t" + snp.altAllele + "\t" + snp.refReads.size() + "\t" + snp.altReads.size() + "\t" + Utils.binomialTest(snp.refReads.size(), snp.altReads.size()) + "\t" + snp.refReads + "\t" + snp.altReads + "\n");
								break;
							}
						}
					}
				}
			}
		}
		bw.close();
		
		processed = new HashSet<String>();
		bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + Parameters.crossTag + "_" + Parameters.sample1 + "_"+ Parameters.sample2 + "_per_eqtl_naive.txt"));
		bw.write("gene\tsnp_chr\tsnp_pos\tsnp_ref\tsnp_alt\t"+Parameters.sample1+"\t"+Parameters.sample2+"\tbinomial_p\t"+Parameters.sample1+"_reads\t"+Parameters.sample2+"_reads\n");
		for(String chr:GTF.forest.keySet())
		{
			IntervalTree chrom = (IntervalTree)GTF.forest.get(chr);
			for(Interval i:chrom.getIntervals())
			{
				IntervalLabelled exon = (IntervalLabelled)i;
				ArrayList<SNP> eqtls = Parameters.eQTLNAIVEPerGENE.get(exon.gene);
				if(eqtls != null)
				{
					for(SNP eqtl:eqtls)
					{
						// Check if this eQTL is in our data
						for(SNP snp:exon.variants)
						{
							if(snp.chr.equals(eqtl.chr) && snp.loc == eqtl.loc && !processed.contains(exon.gene+":"+snp.chr+":"+snp.loc)) // It IS!
							{
								processed.add(exon.gene+":"+snp.chr+":"+snp.loc);
								bw.write(exon.gene + "\t" + snp.chr + "\t" + snp.loc + "\t" + snp.refAllele + "\t" + snp.altAllele + "\t" + snp.refReads.size() + "\t" + snp.altReads.size() + "\t" + Utils.binomialTest(snp.refReads.size(), snp.altReads.size()) + "\t" + snp.refReads + "\t" + snp.altReads + "\n");
								break;
							}
						}
					}
				}
			}
		}
		bw.close();
	}
}
