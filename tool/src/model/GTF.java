package model;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import htsjdk.tribble.index.interval.Interval;
import htsjdk.tribble.index.interval.IntervalTree;
import tools.Utils;

public class GTF 
{
	public static HashMap<String, IntervalTree> forest = null;
	
	public static void readGTF() throws Exception
	{
		System.out.println("\nReading GTF file provided: " + Parameters.inputGTFFile.getAbsolutePath());
		BufferedReader br = Utils.readGTF(Parameters.inputGTFFile);
		String line = br.readLine();
		forest = new HashMap<>();
		int nbExons = 0;
		int nbGenes = 0;
		Parameters.uniqueGeneId = new HashSet<String>();
		ArrayList<String> uniqueGeneName = new ArrayList<String>(); // It's actually not unique and should not be since two geneID can have the same gene Name => ArrayList
		while(line != null)
		{
	    	if(!line.startsWith("#"))
    		{
				String[] tokens = line.split("\t");
				// Parse Line
				String[] params = tokens[8].split(";");
				long start = Long.parseLong(tokens[3]);
				long end = Long.parseLong(tokens[4]);
				String gene_name = null;
				String gene_id = null;
				String chr = tokens[0];
				String type = tokens[2];
				boolean strand = tokens[6].equals("+");
				for(String param:params) 
				{
					String[] values = param.trim().split("\\s+");
					values[1] = values[1].replaceAll("\"", "");
					if(values[0].equals("gene_name")) gene_name = values[1];
					else if(values[0].equals("gene_id")) gene_id = values[1];
				}
				if(gene_name == null) gene_name = gene_id;
				// Which type is it?
				if(type.equals("exon")) 
				{
					nbExons++;
					IntervalTree tree = forest.get(chr);
					if(tree == null) tree = new IntervalTree();
					if(Parameters.uniqueGeneId.add(gene_id)) uniqueGeneName.add(gene_name);
					tree.insert(new IntervalLabelled((int)start, (int)end, gene_id, strand));
					forest.put(chr, tree);
				}
				else if(type.equals("gene")) nbGenes++;
			}
			line = br.readLine();
		}
		br.close();
		
		if(nbGenes == 0) {
			System.out.println("No Genes were detected in the GTF file. Probably the \"gene\" annotation is missing from the GTF file 3rd column?");
			System.out.println("Collapsing exons to their annotated gene_id...");
			nbGenes = Parameters.uniqueGeneId.size();
		}

		System.out.println(nbExons + " 'exons' are annotating " + Parameters.uniqueGeneId.size() + " unique genes in the provided GTF file. In total " + nbGenes + " 'gene' annotations are found in the GTF file.");
		
		if(nbGenes == 0) {
			System.err.println("We couldn't parse the GTF file. Please report this problem if the GTF is in standard format. Or use another GTF from another source.");
			System.exit(-1);
		}
	}
	
	public static HashSet<IntervalLabelled> findOverlappingGenes(String chr, int start, int end, boolean readNegativeStrandFlag)
	{
		HashSet<IntervalLabelled> result = new HashSet<>();
		IntervalTree itree = forest.get(chr);
		if(itree != null)
		{
			List<Interval> found = itree.findOverlapping(new Interval(start, end));
			for(Interval i:found) 
			{
				IntervalLabelled g = (IntervalLabelled)i;
				if(Parameters.stranded == Strand.NO) result.add(g);
				else if(Parameters.stranded == Strand.REVERSE && g.readNegativeStrandFlag == readNegativeStrandFlag) result.add(g);
				else if(Parameters.stranded == Strand.YES && g.readNegativeStrandFlag != readNegativeStrandFlag) result.add(g);
			}
		}
		return result;
	}
	
	public static HashSet<IntervalLabelled> findOverlappingGenes(String chr, int start, int end)
	{
		HashSet<IntervalLabelled> result = new HashSet<>();
		IntervalTree itree = forest.get(chr);
		if(itree != null)
		{
			List<Interval> found = itree.findOverlapping(new Interval(start, end));
			for(Interval i:found) 
			{
				IntervalLabelled g = (IntervalLabelled)i;
				result.add(g);
			}
		}
		return result;
	}
	
	public static void writeOutputPerGene() throws IOException
	{
		HashMap<String, Result> counts = new HashMap<String, Result>();
		for(String geneId:Parameters.uniqueGeneId) counts.put(geneId, new Result());
		
		for(String chr:forest.keySet())
		{
			IntervalTree chrom = (IntervalTree)forest.get(chr);
			for(Interval i:chrom.getIntervals())
			{
				IntervalLabelled exon = (IntervalLabelled)i;
				counts.get(exon.gene).total += exon.count;
				counts.get(exon.gene).allSnps += exon.variants.size();
				for(SNP snp:exon.variants)
				{
					if(snp.refReads.size() != 0 || snp.altReads.size() != 0) counts.get(exon.gene).coveredSnps++;
					if(snp.isRef == 1)
					{
						counts.get(exon.gene).s1.addAll(snp.refReads);
						counts.get(exon.gene).s2.addAll(snp.altReads);
					}
					else if(snp.isRef == 2)
					{
						counts.get(exon.gene).s2.addAll(snp.refReads);
						counts.get(exon.gene).s1.addAll(snp.altReads);
					}
				}
			}
		}
		
		// Compute stats
		String[] keys = Utils.sortKeys(counts);
		double[] pvalue = new double[keys.length];
		for(int i = 0; i < keys.length; i++) pvalue[i] = Utils.binomialTest(counts.get(keys[i]).s1.size(), counts.get(keys[i]).s2.size());
		double[] fdr = Utils.p_adjust(pvalue, "fdr");
		
		// Write output file
		BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + Parameters.crossTag + "_" + Parameters.sample1 + "_"+ Parameters.sample2 + "_per_gene.txt"));
		bw.write("gene\tall_snps\tcovered_snps\t"+Parameters.sample1+"\t"+Parameters.sample2+"\tbinomial_p\tbinomial_FDR\tcovering_reads\n");
		for(int i = 0; i < keys.length; i++) bw.write(keys[i] + "\t" + counts.get(keys[i]).allSnps + "\t" + counts.get(keys[i]).coveredSnps + "\t" + counts.get(keys[i]).s1.size() + "\t" + counts.get(keys[i]).s2.size() + "\t" + pvalue[i] + "\t" + fdr[i] + "\t" + counts.get(keys[i]).total + "\n");
		bw.close();
	}
	
	public static void writeOutputPerSNP() throws IOException
	{
		BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + Parameters.crossTag + "_" + Parameters.sample1 + "_"+ Parameters.sample2 + "_per_snp.txt"));
		bw.write("snp_chr\tsnp_pos\tsnp_ref\tsnp_alt\t"+Parameters.sample1+"\t"+Parameters.sample2+"\tbinomial_p\t"+Parameters.sample1+"_reads\t"+Parameters.sample2+"_reads\n");
		for(String chr:forest.keySet())
		{
			IntervalTree chrom = (IntervalTree)forest.get(chr);
			for(Interval i:chrom.getIntervals())
			{
				IntervalLabelled exon = (IntervalLabelled)i;
				for(SNP snp:exon.variants)
				{
					if(snp.isRef == 1)
					{
						bw.write(snp.chr + "\t" + snp.loc + "\t" + snp.refAllele + "\t" + snp.altAllele + "\t" + snp.refReads.size() + "\t" + snp.altReads.size() + "\t" + Utils.binomialTest(snp.refReads.size(), snp.altReads.size()) + "\t" + snp.refReads + "\t" + snp.altReads + "\n");
					}
					else if(snp.isRef == 2)
					{
						bw.write(snp.chr + "\t" + snp.loc + "\t" + snp.refAllele + "\t" + snp.altAllele + "\t" + snp.altReads.size() + "\t" + snp.refReads.size() + "\t" + Utils.binomialTest(snp.altReads.size(), snp.refReads.size()) + "\t" + snp.altReads + "\t" + snp.refReads + "\n");
					}
				}
			}
		}
		bw.close();
	}
}

