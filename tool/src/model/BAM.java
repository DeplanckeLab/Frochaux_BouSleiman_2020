package model;

import java.util.HashSet;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import tools.Utils;

public class BAM 
{
	/**
	 * Using Picard to read the reads from the BAM file created by the alignment tool
	 * @throws Exception Yes I know...
	 */
	public static void readBAM() throws Exception
	{
		HashSet<String> correctRG = new HashSet<>();
		Long start = System.currentTimeMillis();
		System.out.println("\nReading the reads from the BAM file...");
		SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		SamReader samReader = samReaderFactory.open(Parameters.inputBAMFile);
		SAMRecordIterator it = samReader.iterator();
		SAMFileHeader header = samReader.getFileHeader();
		for(SAMReadGroupRecord readGroup:header.getReadGroups()) if(readGroup.getSample().equals(Parameters.crossTag)) correctRG.add(readGroup.getId());
		System.out.print(correctRG.size() + " read group id(s) are matching your sample name ("+Parameters.crossTag+"):");
		for(String rg:correctRG) System.out.print("\t" + rg);
		System.out.println();
		if(correctRG.size() == 0)
		{
			System.err.println("Error: No @RG tag is matching your sample name ('-cross' option) in BAM header");
			System.exit(-1);
		}
		
		// Start reading the BAM file
		HashSet<String> uniqueGenes = new HashSet<String>();
		Parameters.nbReads = 0;
		while(it.hasNext())
		{
			SAMRecord samRecord = it.next();
			Parameters.nbReads++;
			String rgtag = samRecord.getStringAttribute("RG");
			if(rgtag == null)
			{
				System.err.println("Error: Cannot find @RG tag for read " + samRecord.getReadName() + " (l." + Parameters.nbReads + ")");
				System.exit(-1);
			}
			if(correctRG.contains(rgtag)) Parameters.readMatchingCross++;
			if(samRecord.getReadUnmappedFlag()) Parameters.unmapped++;
			else if(samRecord.getMappingQuality() < 10) Parameters.toolowAqual++;
			else
			{
				if(samRecord.getSupplementaryAlignmentFlag()) Parameters.notUnique++; // Still consider multiple mapped reads
				HashSet<IntervalLabelled> overlappingGenes = Utils.getOverlappingGenes(samRecord.getReferenceName(), samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(), samRecord.getCigar(), samRecord.getReadNegativeStrandFlag());
				HashSet<String> genes = new HashSet<>();
				for(IntervalLabelled i:overlappingGenes) genes.add(i.gene);
				if(genes.size() == 0) Parameters.noFeature++;
				else if(genes.size() == 1) // Unique gene is mapped
				{
					IntervalLabelled gene = overlappingGenes.iterator().next();
					uniqueGenes.add(gene.gene);
					gene.count++;
					Parameters.mapped++;
					if(!gene.variants.isEmpty() && correctRG.contains(rgtag))
					{
						Parameters.readOverlapSNPs++;
						for(SNP snp:gene.variants) snp.assign(samRecord);
					}
				}
				else Parameters.ambiguous++;
			}
			if(Parameters.nbReads %Parameters.chunkSize == 0) System.out.println(Parameters.nbReads + " reads were processed from BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
		}
		samReader.close();
		
		System.out.println(Parameters.nbReads + " reads were processed from BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");

		System.out.println();
		System.out.println(uniqueGenes.size() + " Unique genes were detected (at least one read).");
		System.out.println(Parameters.mapped + " 'Mapped' reads (" + Parameters.pcFormatter.format(((float)Parameters.mapped / Parameters.nbReads) * 100) + "%) from which " + Parameters.notUnique + " are 'not Unique' reads (" + Parameters.pcFormatter.format(((float)Parameters.notUnique / Parameters.nbReads) * 100) + "%)");
		System.out.println(Parameters.ambiguous + " 'Ambiguous' reads (" + Parameters.pcFormatter.format(((float)Parameters.ambiguous / Parameters.nbReads) * 100) + "%)");
		System.out.println(Parameters.noFeature + " 'No Features' reads (" + Parameters.pcFormatter.format(((float)Parameters.noFeature / Parameters.nbReads) * 100) + "%)");
		System.out.println(Parameters.unmapped + " 'Not Aligned' reads (" + Parameters.pcFormatter.format(((float)Parameters.unmapped / Parameters.nbReads) * 100) + "%)");
		System.out.println(Parameters.toolowAqual + " 'Too Low aQual' reads (" + Parameters.pcFormatter.format(((float)Parameters.toolowAqual / Parameters.nbReads) * 100) + "%)");

		System.out.println();
		System.out.println(Parameters.unknownSNP +  " reads have SNPs that does not map to REF or ALT as defined in VCF file");
		
		System.out.println();
		System.out.println(Parameters.readMatchingCross + " reads are matching the cross of interest (" + Parameters.pcFormatter.format(((float)Parameters.readMatchingCross / Parameters.nbReads) * 100) + "%)");
		System.out.println(Parameters.readOverlapSNPs + " reads are overlapping SNPs of interest (" + Parameters.pcFormatter.format(((float)Parameters.readOverlapSNPs / Parameters.nbReads) * 100) + "%)");
	}
}
