## Variant Identification
### Germline Callers
1. **[IMPUTE2] (http://mathgen.stats.ox.ac.uk/impute/impute_v2.html)**
  * Description: phasing observed genotypes and imputing missing genotypes uses reference panels to provide all available halotypes, does not use population labels or genome-wide measures; designed to represent variation in one population; Fairly popular
  * Input: 
   * Reference Haplotypes: Links to 1000 Genomes and HapMap downloads
  * Output: 
1. **[FreeBayes] (https://github.com/ekg/freebayes)**
  * Description: finds SNPs, Indels, MNPs; reports variants based on alignment; haplotype based
  * Input: BAM- uses BAMtools API to parse
   * Reference genome: FASTA
  * Output: VCF
1. **[SOAPindel] (http://soap.genomics.org.cn/soapindel.html)**
  * Description: detects indels from NGS paired-end sequencing
  * Input: files with read alignment can be SOAP or SAM formats, users must also give raw reads in Fasta or Fastq
   * Reference Sequence used to align reads: FASTA
  * Output:
  <hr/>
1. **[2Kplus2] (https://github.com/danmaclean/2kplus2)**
  * Description: algorithm searches graphs produced by de novo assembler Cortex; c++ source code for SNP detection
“2kplus2.cpp is a c++ source code for the detection and the classification of single nucleotide polymorphisms in transformed De Bruijn graphs using Cortex assembler.”
  * Input:
  * Output:
1. **[Atlas 2] ( https://www.hgsc.bcm.edu/software/atlas-2)**
  * Description: specializes in separation of true SNPs and indels from sequencing and mapping errors, last update January 2013
  * Input: takes BAM file,
   * Reference Genome: FASTA
  * Output: produces VCF
1. **[CRISP] (https://sites.google.com/site/vibansal/software/crisp)**
  * Description: identifies SNPs and INDELs from pooled high-throughput NGS, not used for analysis of single samples; implemented in C and uses SAMtools API; latest version should work with diploid genomes
  * Input: requires BAM files (aligned with GATK)
   * Reference Genome: indexed FASTA file
  * Output: VCF files
1. **[Dindel] (http://www.sanger.ac.uk/resources/software/dindel/)**
  * Description: (Wellcome Trust Sanger) calls small indels from short-read sequences, only can handle Illumina data; cannot test candidate indels; written in C++, used on Linux based and Mac computers (not tested in windows)
  * Input: BAM files
  * Output: VCF
1. **[discoSnp++] (http://colibread.inria.fr/software/discosnp/)**
  * Description: detects homozygous and heterozygous SNPs and Indels; software composed of 2 modules (kissnp2 and kissreads)
  * Input: raw NGS datasets; fasta, fastq, gzipped or not;
   * no reference genome required; read pairs can be given
  * Output: FASTA
1. **[FamSeq] (http://odin.mdacc.tmc.edu/~wwang7/FamSeqIndex.html)**
  * Description: family-based sequencing studies- provides probability of an individual carrying variant based on family’s raw measurements; accommodates de novo mutations, can perform variant calling at chrX;
  * Input: VCF
  * Output: VCF
1. **[GeneticThesaurus] (http://sourceforge.net/p/ geneticthesaurus/wiki/Example/)**
  * Description: “Annotation of genetic variants in repetitive regions”
  * Input: Initial variant calling from bam ? vcf output
   * Reference Genome: need to provide own fasta file for hg19 genome,  
  * Output: vcf.gz, vtf.gz, and baf.tsv.gz output
1. **[glfMultiples] (http://genome.sph.umich.edu/wiki/GlfMultiples)**
  * Description: command-line, variant caller
  * Input: GLF
  * Output: VCF
1. **[glfSingle] (http://genome.sph.umich.edu/wiki/GlfSingle)**
  * Description: uses likelihood-based model for variant calling, starts from genotype likelihoods that have been computed from other tools (ex. Samtools BAQ), the likelihoods combine with individual-based prior p(genotype) to generate posterior probabilities
  * Input: GLF
  * Output: VCF
1. **[Halvade] (https://github.com/ddcap/halvade)**
  * Description: command-line; written in Java, “to run halvade a reference is needed for both GATK and BWA and a SNP (dbSNP!) database is required
  * Input: FASTQ
  * Output: VCF
1. **[indelMINER] (https://github.com/aakrosh/indelMINER)**
  * Description: identifies indels from paired-end reads
  * Input: BAM (aligned in SAMtools API)
  * Output: VCF
1. **[Indelocator] (https://www.broadinstitute.org/cancer/cga/indelocator)**
  * Description: (Broad Institute): does not perform realignment, relies on alignments in BAM files (BAM files need aligned before put into indelocator); recommended to use GATK prior;
  * Input: 2 BAM files(tumor & normal), annotated as germline or somatic; also has single sample mode
  * Output: “Output of Indelocator is a high-sensitivity list of putative indel events containing large numbers of false positives. The statistics reported for each event have to be used to custom-filter the list in order to lower false positive rate”
1. **[Isaac Variant Caller] (https://github.com/sequencing/isaac_variant_caller)**
  * Description: detects SNPs and small indels from diploid sample; designed to run on “nux-like platforms”
  * Input: BAM
  * Output: VCF
1. **[KvarQ] (http://www.swisstph.ch/kvarq)**
  * Description: in silico genotyping for selected loci in bacterial genome, written in Python and C
  * Input: FASTQ
   * reference genome or de novo assembly not needed
  * Output:
1. **[LoFreq] (http://sourceforge.net/projects/lofreq/files/)**
  * Description: SNV caller, Python language, standalone program, uncovers cell-population heterogeneity from high-throughput sequencing datasets; calls variants found in <.05% of the population
  * Input: BAM file input? suggest running through GATK
  * Output:
1. **[Manta] (https://github.com/Illumina/manta)**
  * Description: Calls indels and SVs from paired end reads; standalone, command line program; Written in C++ and Python
  * Input: BAM (can tolerate non-paired-end reads); a matched tumor sample may be provided as well
  * Output: VCF
1. **[MarginAlign] (https://github.com/benedictpaten/marginAlign)**
  * Description: SNV caller, specifically tailored to Oxford Nanopore Reads, written in Python; Package comes with 3 programs, marginAlign, marginCaller (calls SNVs), marginStats (computes qc stats on sam files)
  * Input: SAM
  * Output: SAM
1. **[MendelScan] (http://gmt.genome.wustl.edu/packages/mendelscan/)**
  * Description: Last release March 2014; for analyzing sequencing data in family studies of inherited diseases; variant calls for a family in VCF file; still in alpha-testing
on github, example data uses 1000 genomes dataset
  * Input:
  * Output:
1. **[nanopore] (https://github.com/mitenjain/nanopore)**
  * Description: UCSC Nanopore group (group at UCSC studying using ion channels for analysis of single RNA/DNA structures) software pipeline; tailored to Oxford Nanopore Reads; command line program
  * Input: FASTQ
   * Reference files: FASTA
  * Output: “For each possible pair of read file, reference genome and mapping algorithm an experiment directory will be created in the nanopore/output directory.”
1. **[Platypus] (http://omictools.com/platypus-s1989.html)**
  * Description: Package program, written in C, Python, Cython;  Can identify SNPs, MNPs, short indels, and larger variants; has been tested on very large datasets (1000 genomes)
  * Input: BAM
   * Reference Genome: FASTA (files must be indexed using Samtools or similar program
  * Output: VCF
1. **[QualitySNPng] (http://www.bioinformatics.nl/QualitySNPng/)**
  * Description: detection of SNPs; “can be used as a standalone application with graphical user interface as part of pipeline system”; does not require fully sequenced reference genome; haplotype strategy
  * Input:SAM, ACE
  * Output: GUI
1. **[ReviSTER] (http://revister.sourceforge.net/)**
  * Description: command line program; automated pipeline; utilizes BWA, BLAT, and SAMTools;  utilizes BWA mapping program;
  * Input: FASTQ,
   * Reference sequence file and list file containing STR locations as inputs
  * Output: SAM
1. **[RVD] (http://dna-discovery.stanford.edu/software/rvd/)**
  * Description: command-line program, detection of rare SNVs, relies upon Samtools, can be run in MATLAB
  * Input: BAM
   * Reference Genome: FASTA
  * Output: “The algorithm output is a call table -- a comma-separated file with one line for each base position and each line in the following format:
   * AlginmentReferencePosition, AlignmentBase, Call ,SecondBase, CenteredErrorPrc, ReferenceErrorPrc, SecondBasePrc”
1. **[SNVer] (http://snver.sourceforge.net/)**
  * Description: calls common and rare variants in pool or individual NGS data, reports overall p-value, operating system independent statistical tool, identifies SNPs and INDELs, written in Java, no dependencies, straightforward command-line    
   * (SNVerGUI=GUI version) --SNVerGUI: desktop tool for variant detection
  * Input: chrX annotation, sam.zip, bam.zip
   * reference file must be aligned to the data file
  * Output:
1. **[SNVMix] (http://compbio.bccrc.ca/software/snvmix/)**
  * Description: detects SNVs from NGS, post-alignment tool
  * Input: pileupformat (Maq or Samtools)
  * Output:
1. **[SV-M] (http://www.bsse.ethz.ch/mlcb/research/bioinformatics-and-computational-biology/structural-variant-machine--sv-m-.html)**
  * Description: Structural Variant Machine - predicts indels, uses split read alignment profiles, validated by Sanger Sequencng
  * Input:paired-end Illumina reads from 1001 genomes project (uses ref plant- 1001genomes.org)
  * Ouptut:
1. **[SNPest] (https://github.com/slindgreen/SNPest)**
  * Description: Standalone program, language C++, Perl
  * Input: mpileup (SAMtools)
  * Output: VCF
1. **[TrioCaller] (http://genome.sph.umich.edu/wiki/TrioCaller)**
  * Description:Command line program, relies on BWA and samtools; genotype calling for unrelated individuals and parent-offspring trios
  * Input: BAM (that has been aligned in BWA and Samtools
  * Output: BCF that can be formatted to VCF using bcftools
1. **[Snippy] (http://www.vicbioinformatics.com/software.snippy.shtml)**
  * Description: finds indels between haploid reference genome and NGS sequence reads
  * Input:read files- FASTQ or FASTA (can be .gz compressed), output- .aln, .tab, .txt
   * Reference genome in FASTA or GENBANK
  * Output: 
1. **[VntrSeek] (http://orca.bu.edu/vntrseek/)**
  * Description: pipeline for discovering microsatellite tandem repeats with high-throughput sequencing data
  * Input: gzip-compressed FASTA or FASTQ
  * Output: VCF files; one for TRs and observed alleles, another file contains link to viewer
