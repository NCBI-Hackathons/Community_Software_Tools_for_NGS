## Downstream Analysis of Variants
1. **[PrediXcan] (https://github.com/hakyimlab/PrediXcan https://github.com/hriordan/PrediXcan/)**
  * Description: command-line, standalone package program; available in Perl, Python, and R versions; predicts liklihood of a gene being related to a certain phenotype- “that directly tests the molecular mechanisms through which genetic variation affects phenotype.”; no actual expression data used, only in silico expression; “PrediXcan can detect known and novel genes associated with disease traits and provide insights into the mechanism of these associations.”
  * Input: genotype and phenotype file (doesn’t specify file type)
  * Output:default values: genelist, dosages (file format:  snpid rsid) , dosage_prefix, weights, output
  <hr/>
1. **[ATHENA] (http://ritchielab.psu.edu/software/athena-downloads)**
  * Description: Analysis Tool for Heritable and Environmental Network Associations; software package, combines machine learning model with biology and statistics to predict non-linear interactions
  * Input: Configuration file, Data file, Map file (includes rsID)
  * Output: Summary file, Best model file, dot file, individual score file, cross-validation file
1. **[CCRaVAT and QuTie] (http://www.sanger.ac.uk/resources/software/rarevariant/#t_2)**
  * Description: (Wellcome Trust Sanger) Case-Control Rare Variant Analysis Tool and Quantitative Trait; software packages for large-scale analysis of rare variants
  * Input: PED file and MAP file
  * Output: Five tab-delimited txt files
1. **[GCTA] (http://cnsgenomics.com/software/gcta/)**
  * Description: Genome Wide Complex Trait Analysis; package program, command line interface; estimates variance by all SNPs; 5 main functions: “data management, estimation of the genetic relationships from SNPs, mixed linear model analysis of variance explained by the SNPs, estimation of the linkage disequilibrium structure, and GWAS simulation”
  * Input: PLINK binary PED files, MACH output format
  * Output:
1. **[GenomeComb] (http://genomecomb.sourceforge.net/)**
  * Description: package for analysis of complete genome data; annotation using public data or custom tracks, automated primer desing for Sanger or Sequenom validation; “The cg process_illumina command can be used to generate annotated multisample data starting from fastq files, using tools such as bwa for alignment and GATK and samtools for variant calling. Sequencing data can also be imported from Complete Genomics (cg_process_sample command), Real Time Genomics (cg_process_rtgsample command) and VariantCallFormat (VCF) variant files (vcf2sft command).”
  * Input: Sequencing data from Complete Genomics, Illumina, SOLiD and VCF;
  * Output: standard file format used is a simple tab delimited file (.sft, .tsv)
1. **[Genome Track Analyzer] (http://ancorr.eimb.ru/)**
  * Description: compares genome tracks; allows user to compare DNA expression/binding;
  * Input: multiple: SGR/TXT, BED, BED6, GFF;  if using prealigned sequence data- use MACS peak caller: BAM, BED, SAM, ELAND
  * Output:
1. **[GVCBLUP] (http://animalgene.umn.edu/gvcblub)**
  * Description: animal gene mapping; “genomic prediction and variance component estimation of additive and dominance effects”; standalone program, command line interface, writting in C++ and Java
  * Input:
  * Output:
1. **[HOMOG] (http://www.jurgott.org/linkage/homog.htm)**
  * Description: Analyzes heterogeneity with respect to single marker loci or known maps of markers; Carries out homogeneity test for alternative hypothesis “Two family types, one with linkage betweeen a trait to a marker or map of markers, the other without linkage”
  * Input: HOMOG.DAT - described on website
  * Output: HOMOG.OUT
1. **[INTERSNP] (http://intersnp.meb.uni-bonn.de/)**
  * Description: GWIA for case-control SNP and quantitative traits; selected for joint analysis using priori information; Provides linear regression framework, Pathway Association Analysis, Genome-wide Haplotype Analysis,
  * Input: PLINK input formats (ped/map, tped/tfam, bed/bim/fam) Compatible with SetID files
   * Gene reference file: Ensembl Release 75
  * Output: covariance matrix for regression models
1. **[mtSet] (https://github.com/PMBio/mtSet)**
  * Description: Currently only the standalone version available, but moving to LIMIX software suite; offers set tests- allows for testing between variants and traits; accounts for confounding factors ex. relatedness
  * Input: sample-to-sample genetic covariance matrix needs to be computed; multiple types of input; simulator requires input genotype and relatedness component;
  * Output: resdir (result file of analysis), outfile (test statistics and p-values), manhattan_plot (flag)
1. **[MultiBLUP] (http://dougspeed.com/multiblup/)**
  * Description: Package program, command line interface; constructs linear prediction models; Best Linear Unbiased Prediction; improves upon BLUP involving kinship matrices; options: pre-specified kinships, regional kinships, adaptive multiblups, LD weightings
  * Input: PLINK format
  * Output:.reml, .indi.blp
