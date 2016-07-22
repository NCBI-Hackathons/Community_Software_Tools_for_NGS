## CNV (copy number variations) Identification  
1. **[cn.MOPS] (http://bioconductor.org/packages/2.12/bioc/html/cn.mops.html)**
  * Description: Bioconductor; “Mixture of Poissons for CNV detection of NGS data”; models depth of coverage at each genomic position to reduce read count bias; Bayesian approach, written in C++
  * Input: Read count matrices, genomic range objects (has ability to convert BAM files to these)
  * Output:
1. **[CNA-Seq] (http://chipster.csc.fi/manual/cna-define-experiment.html)**
  * Description: Tool in Chipster analytical pipeline, count reads in bins, segment, and call copy number aberrations, plots CNV profiles, identifies common regions, tests DNA copy number induced differential expression, plots combined profiles and copy numbers, adds cytogenic bands, and counts overlapping CNVs by comparing database of genomic variants (http://dgv.tcag.ca/dgv/app/home)
  * Input:
  * Output:
1. **[CNValidator] (https://code.google.com/p/cnvalidator/)**
  * Description: Package program, Command line interface; “Identifies CNVs based on homozygous SNPs and the ratio of heterozygous SNP reads”;
  * Input: text based tab-delimited files
  * Output: txt output
1. **[CNVer] (http://compbio.cs.toronto.edu/cnver/)**
  * Description: Package program, command line interface; “supplements the depth-of-coverage with paired-end mapping”; Donor graph- unified framework- reduces sequence biases; only works with human build hg19 and hg18
  * Input: BAM, tab delimited text
  * Output:
1. **[CNVnator] (https://github.com/abyzovlab/CNVnator)**
  * Description:Tool for CNV discovery and genotyping 
  * Input: Sam/Bam
  * Output: CNV calls or genotype results
1. **[CNVrd2] (http://www.bioconductor.org/packages/devel/bioc/html/CNVrd2.html)**
  * Description: Part of Bioconductor; package program with command line interface; “measure human gene copy number for multiple samples”; written in R
  * Input: multiple; imports from DNAcopy, IRange, Rsamtools
  * Output:
1. **[CONTRA] (http://sourceforge.net/projects/contra-cnv/)**
  * Description: detects CNVs in whole-exome data, requires Python and R, dependent on BEDtools and SAMtools
  * Input: BAM/SAM
   * Reference genome: Ensembl or 1000genomes
  * Output: VCF
1. **[CopyCat] (http://tvap.genome.wustl.edu/tools/copycat/)**
  * Description: Currently under “active development” according to website; utilizes multi-core architecture; “detects copy number aberrations by measuring the depth of coverage
  * Input:
  * Output:
1. **[GASV] (https://code.google.com/p/gasv/)**
  * Description: Geometric Analysis of Structural Variants
  * Input: BAM
  * Output:
1. **[GENSENG/AS-GENSENG] (http://sourceforge.net/projects/asgenseng/?source=directory)**
  * Description: uses Markhov model; read-depth-based method; identifies regions of discrete copy-number changes while accounting for confounders
  * Input: BAM
  * Output: two txt files
1. **[GROM-RD] (http://grigoriev.rutgers.edu/software/grom-rd/index.html)**
  * Description: Standalone program, command line interface; analyzes multiple biases to detect CNVs. More sensitive than CNVnator or RDXplorer
  * Input: BAM
  * Output: “union set of predicted CNV’s from two pipelines”
1. **[m-HMM] (https://www.stt.msu.edu/users/hengwang/mHMM.html)**
  * Description: written in R, the “mixture hidden Markov model” to detect CNVs, package program, command-line interface
  * Input:
  * Output: “data frame of the resulting CNV detection”
1. **[Magnolya] (http://sourceforge.net/projects/magnolya/)**
  * Description: “de novo CNV detection by co-assembly”; No reference genome needed, compared two NGS datasets; package program, command-line interface, written in Python
  * Input: FASTA
  * Output:
1. **[modSaRa] (http://c2s2.yale.edu/software/modSaRa/)**
  * Description: package program, command-line interface; written in R; “modified Screening and Ranking algorithm”; performs quantile normalization, searches for change-point canidates, eliminates unlikely change points, outputs potential CNV segments
  * Input:
  * Output:
1. **[PEMer] (http://sv.gersteinlab.org/pemer/)**
  * Description: Package program, command-line interface; restricted to academic users only, written in Perl and Python; simulation-based error models; three modules- PEMer workflow, SV-Simulation,   BreakDB (own database);
  * Input: FASTA
  * Output:
1. **[QDNAseq] (http://www.bioconductor.org/packages/release/bioc/html/QDNAseq.html)**
  * Description: (Bioconductor) package program, written in R; WGS method for CNV analysis; dependent on R; “the genome is divided into non-overlapping fixed-sized bins, number of sequence reads in each counted, adjusted with a simultaneous two-dimensional loess correction for sequence mappability and GC content, and filtered to remove spurious regions in the genome”
  * Input: BAM
  * Output: “read counts per bin, which have been corrected, filtered, normalized, and optionally log2-transformed”
1. **[readDepth] (https://github.com/chrisamiller/readDepth)**
  * Description: package program, command-line; written in R; measures depth of coverage obtained by massively parallel sequencing; “niche uses...obsolete on human data”
  * Input: BED
  * Output: segs.dat or alts.dat
1. **[SLOPE] (http://www-genepi.med.utah.edu/suppl/SLOPE/slope_guide.txt)**
  * Description: package program, command-line; written in C++; “detects structural variants from targeted short DNA reads”; quickly detects insertions and deletions
  * Input:SAM/FASTQ/MAQ
   * Reference file: FASTA
  * Output: html file
1. **[TIGRA] (http://bioinformatics.mdanderson.org/main/TIGRA)**
  * Description: package program, command-line, for academic users only; written in C++; “conducts targeted local assembly of SV”; Uses data from 1000 Genomes
  * Input: “SV calls and a set of bam files that contain reads mapped to a reference genome such as NCBI build36.”
  * Output:
1. **[WISECONDOR] (http://omictools.com/wisecondor-s1816.html)**
  * Description: “WIthin-SamplE COpy Number aberration”; package program, command-line; written in Python; “Detect fetal trisomies and smaller CNV's in a maternal plasma sample using whole-genome data.”
  * Input: SAM, BAM
   * Reference file: FASTA
  * Output: STDOUT plot in PDF

##SV Identification
1. **[BreakDancer] (http://gmt.genome.wustl.edu/packages/breakdancer/)**
  * Description: Genome wide detection of structural variants; Includes two complementary programs, BreakDancerMini (small indel detection) and BreakDancerMax (5 types of structural variants)
  * Input: Set of map files produced by a front-end aligner (MAQ, BWA, NovoAlign, and Bfast) and a tab-delimited file containing locations of the map files, detection parameters and sample information
  * Output:File containing 14 columns
1. **[Delly2] (https://github.com/tobiasrausch/delly)**
  * Description: Can discover and genotype deletions, tandem duplications, inversions and translocations; uses paired-ends and split-reads; SVs can be visualized using Delly-maze and Delly-sauve
  * Input: Bam
  * Output: VCF
1. **[Genome STRiP] (http://www.broadinstitute.org/software/genomestrip/)**
  * Description: Set of tools for structural variant identification and genotyping; requires genomes from multiple individuals (20-30 minimum); Uses GATK and consists of multiple modules
  * Input: BAM
  * Output: VCF
1. **[mrFast] (http://mrfast.sourceforge.net/)**
  * Description: micro-read Fast Alignment Search Tool; designed to map short reads generate with Illumina
  * Input:
  * Output: SAM
1. **[mrsFast] (http://sfu-compbio.github.io/mrsfast/)**
  * Description: Micro-read substitution-only Fast Alignment Search Tool; designed to map short reads to reference genome assemblies
  * Input: FastQ
   *Reference Genome: FASTA
  * Output: SAM
1. **[NovelSeq] (http://novelseq.sourceforge.net)**
  * Description: Designed to detect novel sequence insertions using high-throughput paired-end whole genome sequencing data
  * Input: FASTA or FASTQ
  * Output: FASTA file containing paired-end sequences
1. **[PEMer] (http://sv.gersteinlab.org/pemer/)**
  * Description: Command line package composed of three modules (PEMer workflow, SV-Simulation, and BreakDB) for SV analysis
  * Input: FASTA
  * Output:
1. **[Pindel] (http://gmt.genome.wustl.edu/packages/pindel/)**
  * Description: Command line package program for detection of “large deletions, medium sized insertions, inversions, tandem duplications, and other structural variations”
  * Input: SAM/BAM
  * Output: pindel raw output format (can be converted to VCF)
1. **[SVMerge] (http://svmerge.sourceforge.net/)**
  * Description: Pipeline for detecting structural variants by integrating calls from different SV callers
  * Input: BAM
  * Output: multiple (depends on SV callers used)
1. **[TIGRA] (http://bioinformatics.mdanderson.org/main/TIGRA)**
  * Description: Command line package program that “conducts targeted local assembly of structural variants using targeted iterative graph routing assembly algorithm” Uses data from 1000 Genomes project
  * Input: SV file; group of BAM files
  * Output: Fasta
