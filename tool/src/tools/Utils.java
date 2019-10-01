package tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;

import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import model.GTF;
import model.IntervalLabelled;

public class Utils 
{
	private static BinomialTest bt = new BinomialTest();

	/**
	 * Read VCF different types of VCF file
	 * @param VCF input file
	 * @return BufferedReader handle
	 * @throws Exception
	 */
	public static BufferedReader readVCF(File vcf) throws Exception
	{
		if(vcf.getAbsolutePath().endsWith(".vcf"))
		{
			return new BufferedReader(new FileReader(vcf));
		}
		else if(vcf.getAbsolutePath().endsWith(".vcf.gz"))
		{
			return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(vcf))));
		}
		else
		{
			System.err.println("The extension of the VCF file is not recognized : " + vcf.getAbsolutePath());
			System.err.println("It should be '.vcf', or '.vcf.gz'");
			System.exit(-1);
		}
		return null;
	}
	
	public static double binomialTest(int v1, int v2)
	{
		double p = bt.binomialTest(v1 + v2, v1, 0.5, AlternativeHypothesis.TWO_SIDED);
		if(p > 1) p = 1;
		return p;
	}	
	
    public static double[] p_adjust(final double[] pvalues, String adjMethod)
    {
        if (!adjMethod.equals("BH") && !adjMethod.equals("fdr") && !adjMethod.equals("bonferroni") && !adjMethod.equals("none"))
        {
            System.out.println("This adjustment method is not implemented.");
            System.exit(0);
        }
        if (adjMethod.equals("fdr")) adjMethod = "BH";
        if (adjMethod.equals("none")) return pvalues;

        Integer[] idx = new Integer[pvalues.length];
        for (int i = 0; i < idx.length; i++) idx[i] = i;
        double[] adj_p_value = new double[pvalues.length];

        Arrays.sort(idx, new Comparator<Integer>()
        { // I sort the indexes with respect to the double values
            @Override
            public int compare(Integer o1, Integer o2)
            {
                return Double.compare(pvalues[o2], pvalues[o1]);
            }
        });
        int[] index = new int[idx.length];
        for (int i = 0; i < index.length; i++)
        {
            index[idx[i]] = i;
        }

        for (int i = 0; i < pvalues.length; i++)
        {
            if (adjMethod.equals("BH"))
            {
                adj_p_value[i] = pvalues.length / (pvalues.length - (double) index[i]) * pvalues[i];
            }
            if (adjMethod.equals("bonferroni"))
            {
                adj_p_value[i] = pvalues.length * pvalues[i];
            }
        }

        double min = Double.MAX_VALUE;
        for (int i = 0; i < index.length; i++) // cummin
        {
            double adjP = adj_p_value[idx[i]];
            if (adjP < min)
            {
                min = adjP;
            } else
            {
                adjP = min;
            }
            adj_p_value[idx[i]] = Math.min(1, adjP);
        }

        return adj_p_value;
    }
	
	public static BufferedReader readGTF(File gtf) throws Exception
	{
		if(gtf.getAbsolutePath().endsWith(".gtf"))
		{
			return new BufferedReader(new FileReader(gtf));
		}
		else if(gtf.getAbsolutePath().endsWith(".gtf.gz"))
		{
			return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(gtf))));
		}
		else
		{
			System.err.println("The extension of the GTF file is not recognized : " + gtf.getAbsolutePath());
			System.err.println("It should be '.gtf', or '.gtf.gz'");
			System.exit(-1);
		}
		return null;
	}
	
	public static int getGeno(String geno)
	{		            	
		if(geno.equals("0/0") || geno.equals("0|0")) return 0;
		if(geno.equals("1/0") || geno.equals("0/1") || geno.equals("1|0") || geno.equals("0|1")) return 1;
        if(geno.equals("1/1") || geno.equals("1|1")) return 2;
        return -1; // Handle ./. or 1/2 complex phenotypes => Will be ignored
	}
	
	public static boolean contains(String[] array, String value)
	{
		for(String val:array) if(val.equals(value)) return true;
		return false;
	}
	
	public static HashSet<IntervalLabelled> getOverlappingGenes(String chr, int start, int end, Cigar c, boolean readNegativeStrandFlag) throws Exception
	{
		HashSet<IntervalLabelled> res = new HashSet<>();
		List<CigarElement> l = c.getCigarElements();
		int s = start;
		for(CigarElement cigar:l)
		{
			switch(cigar.getOperator())
			{
				case M:
					res.addAll(GTF.findOverlappingGenes(chr, s, s + cigar.getLength() - 1, readNegativeStrandFlag)); // -1 Because the last letter is at the index before
					s += cigar.getLength();
					break;
				case N:
					s += cigar.getLength();
					break;
				case D:
					s += cigar.getLength();
					break;
				case EQ:
					System.out.println("CIGAR = " + c);
					return new HashSet<>(); // TODO
				case H:
					System.out.println("CIGAR = " + c);
					return new HashSet<>(); // TODO
				case I:
					// Do nothing
					break;
				case P:
					System.out.println("CIGAR = " + c);
					return new HashSet<>(); // TODO
				case S:
					// Do nothing (alignment start is after any S & alignment end before any S)
					break;
				case X:
					System.out.println("CIGAR = " + c);
					return new HashSet<>(); // TODO
			}
		}
		
		if(s != end + 1) throw new Exception("Error while reading CIGAR");
		return res;
	}
	
	public static String[] sortKeys(Map<String, ?> map)
	{
		String[] keys = map.keySet().toArray(new String[map.keySet().size()]);
		Arrays.sort(keys);
		return keys;
	}
	
	public static String toReadableTime(long ms)
	{
		if(ms < 1000) return ""+ms+" ms";
		long s = ms / 1000;
		ms = ms % 1000;
		if(s < 60) return s+" s "+ms+" ms";
		long mn = s / 60;
		s = s % 60;
		if(mn < 60) return mn +" mn "+s+" s "+ms+" ms";
		long h = mn / 60;
		mn = mn % 60;
		if(h < 24) return h +" h "+ mn +" mn "+s+" s "+ms+" ms";
		long d = h / 24;
		h = h % 24;
		return d+ " d " + h +" h "+ mn +" mn "+s+" s "+ms+" ms";
	}
}
