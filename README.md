# public

cisVar.py is my pipeline written in python3 and uses regression_qtls.R. Put both of these scripts in the directory you want all the outputs to be in. regression_qtls.R is called inside of cisVar.py so you will not run it directly. Right now I have a MAF>1% in the R script but if you’d like to change it you can easily do that in the R script.

cisVar.py mpileup -F <SampleName> -f <fastaFile> -p <mpileupBEDfile> -B <sortedBam>
cisVar.py post -F <SampleName> -r <readDepth> -a <allelesFile>
cisVar.py geno -F <SampleName> -r <readDepth> -i <individualsFile> -g <genotypesFile>
cisVar.py qtls -F <SampleName> -r <readDepth> -n <numberIndividuals>

 
These are all the files you’ll need. Also use samtoolsv1.9< and bedtools.

-h	help

-g	The genotypes file is a tab separated file, sorted by individuals(columns) and by snps (rows). 

chr   position    ref alt NA18486 NA18489 ...
chr1  10583   g   a   0   1 ...
chr1  10611   C   G   0   0 ...
chr1  13302   c   t   1   0 … 

-i	The individualsFile is just a list of the individuals (one per line) in your pool that are to be extracted from a larger list of individuals.
		EX:   NA18486
			    NA18489
			    ….

-a	The alleles file is a tab-sep file with no header (chr#, position, Ref_allele, Alt_allele):
chr1	10505	A	T
chr1	10506	C	G
chr1	10511	G	A
	* this file should match the length of the pileup file.

-f	fasta file with all chromosomes concatenated and sorted.

-p	The mpileup.bed file is a tab-sep standard bed file with no header (where the SNP location is the third column and the second is SNP-1.
    EX:

    chr1	10504	10505
    chr1	10505	10506
    chr1	10510	10511
		
-B	sortedBam is the sorted without the -n flag i.e; samtools sort [in.bam]



