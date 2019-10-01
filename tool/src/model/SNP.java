package model;

import java.util.HashSet;
import java.util.List;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class SNP
{
	public String chr;
	public String rsId;
	public long loc;
	public boolean isSNP = false;
	public String refAllele;
	public String altAllele;
	public int isRef = -1;
	public HashSet<String> altReads = new HashSet<>();
	public HashSet<String> refReads = new HashSet<>();
	
	public String getKey()
	{
		return chr+":"+loc+":"+refAllele+":"+altAllele;
	}
	
	public void assign(SAMRecord samRecord)
	{
		if(isSNP) // TODO consider INDELS
		{
			String seq = getAlignedStringAtPos(samRecord);
			if(seq != null)
			{
				if(seq.equals(refAllele)) refReads.add(samRecord.getReadName());
				else if(seq.equals(altAllele)) altReads.add(samRecord.getReadName());
				else Parameters.unknownSNP++;
			}
		}
	}
	
	/**
	 * ONLY for SNPs (optimized in time)
	 * @param samRecord
	 * @return
	 */
	private String getAlignedStringAtPos(SAMRecord samRecord) // WARNING: Only for SNPs
	{
		int readStart = samRecord.getAlignmentStart();
		if(readStart > loc) return null; // Restrict to reads where the SNP fall inside
		if(samRecord.getAlignmentEnd() < loc) return null; // Restrict to reads where the SNP fall inside
		String res = samRecord.getReadString();
		List<CigarElement> l = samRecord.getCigar().getCigarElements();
		int readEnd = readStart;
		int posInRead = 0;
		boolean firstOcc = true;
		for(CigarElement cigar:l)
		{
			switch(cigar.getOperator())
			{
				case M:
					if(readEnd + cigar.getLength() >= loc) 
					{
						if(loc >= readEnd) return ""+res.charAt((int)loc - readEnd + posInRead); // The SNP is within the read
						else return null; // Missed it (prob. because of N or D CIGAR)
					}
					readEnd += cigar.getLength();
					posInRead += cigar.getLength();
					break;
				case N:
					readEnd += cigar.getLength();
					if(readEnd >= loc) return null; // The snp is within a "not mapped" area
					break;
				case D:
					readEnd += cigar.getLength();
					if(readEnd >= loc) return null;
					break;
				case EQ:
					System.err.println("CIGAR = " + samRecord.getCigar());
					System.exit(-1);
				case H:
					System.err.println("CIGAR = " + samRecord.getCigar());
					System.exit(-1);
				case I:
					posInRead += cigar.getLength(); // For an insertion I don't need to test the pos again
					break;
				case P:
					System.err.println("CIGAR = " + samRecord.getCigar());
					System.exit(-1);
				case S:
					posInRead += cigar.getLength();
					if(firstOcc) readStart -= cigar.getLength(); // S is at the beginning of the read, i.e. the beginning of the read string starts before
					else readEnd += cigar.getLength(); // If this is at the end of the read
					break;
				case X:
					System.err.println("CIGAR = " + samRecord.getCigar());
					System.exit(-1);
			}
			firstOcc = false;
		}
		return null;
	}
}