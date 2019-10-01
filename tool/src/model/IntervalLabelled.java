package model;

import java.util.HashSet;

import htsjdk.tribble.index.interval.Interval;

public class IntervalLabelled extends Interval
{
	public HashSet<SNP> variants; // One gene contains many variants
	public String gene;
	public int count;
	public boolean readNegativeStrandFlag;
	
	public IntervalLabelled(int start, int end, String gene, boolean readNegativeStrandFlag) 
	{
		super(start, end);
		this.gene = gene;
		this.readNegativeStrandFlag = readNegativeStrandFlag;
		this.variants = new HashSet<>();
		this.count = 0;
	}
	
	@Override
	public int compareTo(Object o) {
		int res = super.compareTo(o);
		if(res == 0)
		{
			IntervalLabelled i = (IntervalLabelled)o;
			return this.gene.compareTo(i.gene);
		}
		return res;
	}
	
	@Override
	public boolean equals(Object o) {
		boolean res = super.equals(o);
		if(res)
		{
			IntervalLabelled i = (IntervalLabelled)o;
			return this.gene.equals(i.gene);
		}
		return res;
	}
	
	@Override
	public String toString() {
		return super.toString() + " : " + gene + " : " + variants.size();
	}
}
