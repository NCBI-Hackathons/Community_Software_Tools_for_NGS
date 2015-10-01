# List of Community Software Tools for Next Generation Sequencing
## Commercial tools
1. **[Strand NGS] (http://www.strand-ngs.com/)**
  * offers many different tools including alignment, RNA-Seq, DNA-Seq, ChIP-Seq, Small RNA-Seq, Genome Browser, visualizations, Biological Interpretation, etc. Supports workflows
“one can import the sample data in FASTA, FASTQ or tag-count format. In addition, prealigned data in SAM, BAM or Illumina-specific ELAND format can be directly imported for analysis.”
  * Alignment feature: Supports alignment from Illumina, Ion Torrent, 454 (Roche), and Pac Bio
  * DNA-Seq Feature, can annotate with dbSNP
1. **[CLC Genomics Workbench] (http://www.clcbio.com/desktop-applications/top-features/)**<br>
  * (QIAGEN). Features include: resequencing, workflow, read mapping, de novo assembly, variant detection, RNA-Seq, ChIP-Seq, Genome Browser, etc (entire list on website); Main Workbench offers database search (Genbank, Blast, Pubmed); 2000 organizations have invested in CLC
  * Accepts VCF files from 1000 Genomes Project
  * Accepts downloaded tracks from dbSNP
  * Also accepts: FASTA, GFF/GTF/GVF, BED, Wiggle, Cosmic, UCSC variant database, complete genomics master var file
  * Read mapping: “In addition to Sanger sequence data, reads from these high-throughput sequencing machines are supported: The 454 FLX System and the 454 GS Junior System from Roche, Illumina Genome Analyzer, Illumina HiSeq, Illumina HiScan, and Illumina MiSeq sequencing systems, SOLiD system from Life Technologies, Ion Torrent system from Life Technologies, Helicos from Helicos BioSciences”
  * De novo assembly: “In addition to Sanger sequence data, reads from these high-throughput sequencing machines are supported The 454 FLX System and the 454 GS Junior System from Roche, Illumina Genome Analyzer, Illumina HiSeq, Illumina HiScan, and Illumina MiSeq sequencing systems, SOLiD system from Life Technologies, Ion Torrent system from Life Technologies”
  * Annotation tracks from Ensembl
1. **[DNAnexus] (https://www.dnanexus.com/product-overview)**
  * Private cloud repository -- formerly a redistributor of SRA and other NCBI resources; command-line or via web, can fetch data from a URL,  build custom pipeline/ workflow
has sra.dnanexus.com site: data downloads come directly from NCBI
1. **[Ingenuity Variant Analysis] (http://www.ingenuity.com/products/variant-analysis)**
  * (QIAGEN) allows for variant identification and analysis, uses NCI-60 data set for cancer, Supported third part informatin: Entrez Gene, RefSeq, ClinVar; gives contextual details of results instead of just A to B relationship
  * Has own database-- “knowledge base” based on COSMIC, OMIM, and TCGA databases
1. **[NextGENe] (http://www.softgenetics.com/NextGENe.html)**
  * “perfect analytical partner for the analysis of desktop sequencing data produced by the ION PGM™, Roche Junior, Illumina MiSeq as well as high throughput systems as the Ion Torrent Proton, Roche FLX, Applied BioSystems SOLiD™ and Illumina® platforms.” runs on Windows, free-standing multi-application package-- SNP/Indel analysis, CNV prediction and disease discovery, whole genome alignment, etc.  
  * Data can be imported from Clinvar, dbSNP, Genbank: http://www.softgenetics.com/PDF/NextGene_UsersManual_web.pdf
1. **[Partek Genomics Suite (Partek GS)] (http://www.partek.com/pgs)**
  * cited in over 800 articles in 2013, “ability to support all next generation sequencing.. platforms”, analysis and visualization capabilities
  * Can input GEO SOFT files→ http://www.partek.com/Tutorials/microarray/User_Guides/GEO_Import_Guide.pdf
  * RNAseq- Sequenced on Illumina, output- BAM that can supported by Genomics Suite, ELAND, Bowtie, BWA, TopHat-- http://www.partek.com/Tutorials/microarray/NextGen/RNASeqTutorial.pdf
1. **[Golden Helix: SNP and Variation Suite] (http://goldenhelix.com/SNP_Variation/)**
  * used for managing, analyzing and visualizing genotypic and phenotypic data; Features: Genome-wide association studies, genomic prediction, copy number analysis, small sample DNA-Seq workflows, large sample DNA-seq analysis, RNA-seq analysis.
Supported files: .txt, excel XLS & XLSX, CEL, CHP, CNT, Illumina, Plink PED, TPED, BED, Agilent files, NimbleGen data summary files, VCF files, Impute2 GWAS files, HapMap format, MACH output, + 50 other formats
consumes NCBI data directly
1. **[Genomatix] (https://www.genomatix.de/)**
  * Applications: ChIP-Seq, DNA-Seq, RNA-Seq, DNA methylation; enable personalized medicine,
  * Mining Stations: Supports all established NGS sequencing platforms- SOLiD, 454 Life Sciences, Genome Analyzer, HiSeq, MiSeq, IonTorrent
  * Software Suite: can upload sequence of BED files
  * Genome browser: BED and BAM files, Public data- 1500 BED files available for every user
1. **[Biodatomics] (http://www.biodatomics.com/)**
  * Open source platform (SaaS), analysis and genome sequencing tools, integrates over 400 genomic analysis open source tools and pipelines, have a private and public cloud version. Features: genomic data visualization, drag and drop interface, accelerated analysis, real-time collaboration
  *They have a couple modules to do so, and have enabled parts of the sra toolkit
1. **[SolveBio] (https://www.solvebio.com/)**
  * Software product, for clinical genomics professionals, manage, curate, report genomic variation
  * Has own data library -- data from NCBI

## Variant Identification:
### Germline Callers
1. **[IMPUTE2] (http://mathgen.stats.ox.ac.uk/impute/impute_v2.html)**
  * Description: phasing observed genotypes and imputing missing genotypes uses reference panels to provide all available halotypes, does not use population labels or genome-wide measures; designed to represent variation in one population; Fairly popular
  * Input: 
   * Reference Haplotypes: Links to 1000 Genomes and HapMap downloads
  * Output: 
1. **[FreeBayes] (https://github.com/ekg/freebayes)**
  * Description: finds SNPs, Indels, MNPs; reports variants based on alignment; haplotype based
  * Input: BAM- uses BAMtools API to parse
   *Reference genome: FASTA
  *Output: VCF
1. **[SOAPinde] (http://soap.genomics.org.cn/soapindel.html)**
  * Description: detects indels from NGS paired-end sequencing
  * Input: files with read alignment can be SOAP or SAM formats, users must also give raw reads in Fasta or Fastq
   * Reference Sequence used to align reads: FASTA
  * Output:
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
  * Input: Initial variant calling from bam → vcf output
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
  * Input: BAM file input→ suggest running through GATK
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

### Somatic Callers
1. **[Cake] (http://cakesomatic.sourceforge.net/)**
  * Description: standalone program, “pipeline for the integrated analysis of somatic variants in cancer genomes”; integrates four algorithms; written in Perl; required tools: samtools, tabix, vcftools, VarScan2, bambino, cmake, somaticsniper (User guide; workflow page)
  * Input: tumor and normal reads in BAM files, run through variant calling programs to generate intermediate VCF
  * Output: VCF
1. **[MuTect] (http://www.broadinstitute.org/cancer/cga/mutect)**
  * Description: Broad Institute, identification of somatic point mutations in cancer genomes; requires preprocessing of reads (GATK)
  * Input: same as GATK (FASTA reference genome, SAM read files)
  * Output: call-stats, VCF, wiggle files
1. **[Polymutt] (http://genome.sph.umich.edu/wiki/Polymutt)**
  * Description: calls SNVs and detects de novo point mutations in families
  * Input: GLF or BAM or VCF (must have identical chromosome orders)
  * Output: VCF
1. **[Bassovac] (http://tvap.genome.wustl.edu/tools/bassovac/)**
  * Description: Improved Bayesian inversion somatic caller; unlike other software packages, treats effects fully probabilisticallys instead of using ad-hoc modeling; effects are integrated at the atomic level and standard probability theory integrates read tallies to the sample level and to the tumor-normal pair level; "pending public release"
  * Input:
  * Output:
1. **[CLImAT] (http://bioinformatics.ustc.edu.cn/CLImAT/)**
  * Description: standalone program; “accurate detection of copy number alteration and loss of heterozygosity in impure and aneuploid tumor samples using whole genome sequencing data”
  * Input: depth file generated by DFExtract and a config file
  * Output: .results file, .Gtype, LOG.txt, also generates visualization
1. **[DeNovoGear] (http://denovogear.sourceforge.net/)**
  * Description: de-novo variant calling and interpretation; standalone program; dependencies C++ compiler, CMake, HTSlib, Eigen, Boost
  * Input: PED and BCF
  * Output: “The output format is a single row for each putative de novo mutation (DNM), with the following fields”
1. **[EBCall] (https://github.com/friend1ws/EBCall)**
  * Description: Empirical Baysian Mutation Calling; standalone program; uses tumor/normal paired reads and non-paired normal reference samples; dependent on samtools, R and VGAM pack for R
  * Input: BAM
  * Output: not sure what exact type of file- “The format of the result is suitable for adding annotation by annovar.”
1. **[HapMuc] (https://github.com/usuyama/hapmuc)**
  * Description: standalone program; “utilizes the information of heterozygous germline variants near candidate mutations”; Dependent upon- Boost, SAMtools, BEDtools; 3 step workflow
  * Input: BAM
  * Output: BED
1. **[MultiGeMS] (https://github.com/cui-lab/multigems)**
  * Description: Multi-sample Genotype Model Selection
  * Input: .txt, pileup (SAM/BAM converted to pileup format)
  * Output: VCF
1. **[MultiSNV] (https://bitbucket.org/joseph07/multisnv/wiki/Home)**
  * Description: command-line program; calls SNVs from NGS data from multiple samples from the same patient; dependent on R, Git, cmake, Boost and compile libraries
  * Input: BAM or pileup
  * Output: VCF
1. **[MutationSeq] (http://compbio.bccrc.ca/software/mutationseq/)**
  * Description: standalone program, somatic SNV detection in tumor/normal samples; dependent on python, bamtools, boost, and LAPACK
  * Input: BAM
  * Output: VCF4.1 consisting of two parts (meta information & data lines)
1. **[qSNP] (http://www.qcmg.org/bioinformatics/tiki-index.php)**
  * Description: standalone program; SNV caller for somatic variants in “low cellularity cancer samples”
  * Input: BAM, dbSNP data, Illumina data, chrConv
  * Output: “qSNP output files are named using a 4-element pattern: <subject-id>.<snp-type>.<description>.<extension>”
1. **[RADIA] (https://github.com/aradenbaugh/radia/)**
  * Description: RNA and DNA Integrated Analysis for Somatic Mutation Detection; DNA only Method(tumor/normal pair, ignores RNA) or Triple BAM Method (uses all three datasets from same patient); dependent upon python, samtoools, pysam API, BLAT, SnpEff
  * Input: BAM
   * Reference Genome: FASTA indexed with SAMtools faidx
  * Output: VCF
1. **[RVD2] (http://genomics.wpi.edu/rvd2/)**
  * Description: sensitive, variant detection for low-depth targeted NGS data; python module or command- line program;
  * Input: tab- deliminted depth chart format (converted from pileup files)
  * Output: three hdf5 files and a vcf file
1. **[Shimmer] (https://github.com/nhansen/Shimmer)**
  * Description: standalone program; detects somatic SNVs with multiple testing correction, uses Fisher’s exact test; dependent on git, samtools, R, R statmod package; for tumor/normal matched samples
  * Input: BAM
  * Output: VCF
1. **[SNV-PPILP] (http://www.cs.helsinki.fi/en/gsa/snv-ppilp/)**
  * Description: Refines GATK’s Unified Genotyper SNV calls for “multiple samples assumed to form a phylogeny”
  * Input:
  * Output:
1. **[SomaticSniper] (http://gmt.genome.wustl.edu/packages/somatic-sniper/)**
  * Description: command-line application to identify SNPs between tumor/normal pairs- predicts probability of difference between two
  * Input: BAM
   *Reference Genome in FASTA
  * Output: VCF
1. **[Strelka] (https://sites.google.com/site/strelkasomaticvariantcaller/)**
  * Description: somatic variant calling workflow for matched tumor-normal samples; detects indels; runs on *nux-like platform
  * Input: BAM (must be sorted and indexed)- Strelka does own realignment around indels-- don’t need to do this type of pre-processing
  * Output: pair of VCF files
1. **[Triodenovo] (http://www.pitt.edu/~wec47/triodenovo.html)**
  * Description: Bayesian framework for calling de novo mutations in trios
  * Input: VCF file with PL or GL fields (recommend using GATK or samtools to generate)
  * Output: out_vcf
1. **[UNCeqr] (http://lbg.med.unc.edu/~mwilkers/unceqr_dist/)**
  * Description: finds somatic mutations using integration of DNA and RNA seq data-- boosts sensitivity for low purity tumors and rare mutations;
  * Input:”can accept a variety of sequencing inputs and configurations”
  * Output: “table of somatically mutated sites and associated information. These somatic mutations can be annotated with predicted transcript and protein effects using third party tools, such as Annovar”
1. **[Virmid] (http://sourceforge.net/projects/virmid/)**
  * Description: Virtual Microdissection for SNP calling; Java based; for disease-control matched samples; uncovers SNPs with low allele frequency by considering alpha contamination
  * Input: BAM (must be sorted and indexed- samtools sort)
  * Output: VCF and report file

### Germline + Somatic  Callers:
1. **[VarScan 2] (http://massgenomics.org/varscan)**
  * Description: identify germline variants, private and shared variants, somatic mutations, and somatic CNVs; detects indels
  * Input: SAMtools pileup
  * Output: VCF
1. **[BAYSIC] (http://genformatic.com/baysic/)**
  * Description: Bayesian method; combines variant calls from different methods (GATK, FreeBayes, Atlas, Samtools, etc)
  * Input: VCF format from one or more variant calling programs
  * Output: VCF file containing integrated set of variant calls
1. **[MSIsensor] (https://github.com/ding-lab/msisensor)**
  * Description: Microsatellite instability detection; C++ program, detects somatic and germline variants in tumor-normal paired data
  * Input: BAM index files (normal and tumor)
  * Output:
1. **[QuadGT] (http://www.iro.umontreal.ca/~csuros/quadgt/)**
  * Description: software package, SNV calling from normal-tumor pair and two parent genomes; quantifies descent-by-modification relationships; Written in Java
  * Input: BAM files (parsed by Picard/Samtools API)
   * Reference Genome; FASTA
  * Output: VCF
1. **[RAREVATOR] (http://sourceforge.net/projects/rarevator/)**
  * Description: RAre REference VAriant annotaTOR; command line;  “identification and annotation of germline and somatic variants in rare reference allele loci from second generation sequencing data”; Bayesian genotype likelihood model
  * Input: BED or VCF files from GATK
  * Output: two VCF files (one for SNVs, one for Indels)
1. **[Scalpel] (http://scalpel.sourceforge.net/)**
  * Description: Used for detecting indels in a reference genome; performs localized micro-assembly of specific regions of interest; can do single, de novo, somatic reads; requires that raw reads are aligned with BWA
  * Input: BAM
  * Output: either VCF or ANNOVAR
1. **[SOAPsnp] (http://soap.genomics.org.cn/soapsnp.html)**
  * Description: based on Baye’s theorem; calls consensus genotype
  * Input:SOAP short read alignment results
  * Output: GLF, option of flat tabular format
1. **[VariantMaster] (http://sourceforge.net/projects/variantmaster/)**
  * Description: “extract causative variants for monogenic and sporadic genetic diseases”; uses ANNOVAR;
  * Input: BAM or VCF files (from SAMtools, GATK)
  * Output:

## Downstream Analysis of Variants:
1. **[PrediXcan] (https://github.com/hakyimlab/PrediXcan https://github.com/hriordan/PrediXcan/)**
  * Description: command-line, standalone package program; available in Perl, Python, and R versions; predicts liklihood of a gene being related to a certain phenotype- “that directly tests the molecular mechanisms through which genetic variation affects phenotype.”; no actual expression data used, only in silico expression; “PrediXcan can detect known and novel genes associated with disease traits and provide insights into the mechanism of these associations.”
  * Input: genotype and phenotype file (doesn’t specify file type)
  * Output:default values: genelist, dosages (file format:  snpid rsid) , dosage_prefix, weights, output
1. **[ATHENA] (http://ritchielab.psu.edu/software/athena-downloads)**
  * Description: Analysis Tool for Heritable and Environmental Network Associations; software package, combines machine learning model with biology and statistics to predict non-linear interactions
  * Input: Configuration file, Data file, Map file (includes rsID)
  * Output: Summary file, Best model file, dot file, individual score file, cross-validation file
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
  *Output:
1. **[GVCBLUP] (http://animalgene.umn.edu/gvcblub)**
  * Description: animal gene mapping; “genomic prediction and variance component estimation of additive and dominance effects”; standalone program, command line interface, writting in C++ and Java
  * Input:
  * Output:
1. **[INTERSNP] (http://intersnp.meb.uni-bonn.de/)**
  * Description: GWIA for case-control SNP and quantitative traits; selected for joint analysis using priori information; Provides linear regression framework, Pathway Association Analysis, Genome-wide Haplotype Analysis,
  * Input: PLINK input formats (ped/map, tped/tfam, bed/bim/fam) Compatible with SetID files
   *Gene reference file: Ensembl Release 75
  * Output: covariance matrix for regression models
1. **[mtSet] (https://github.com/PMBio/mtSet)**
  * Description: Currently only the standalone version available, but moving to LIMIX software suite; offers set tests- allows for testing between variants and traits; accounts for confounding factors ex. relatedness
  * Input: sample-to-sample genetic covariance matrix needs to be computed; multiple types of input; simulator requires input genotype and relatedness component;
  * Output: resdir (result file of analysis), outfile (test statistics and p-values), manhattan_plot (flag)
1. **[MultiBLUP] (http://dougspeed.com/multiblup/)**
  * Description: Package program, command line interface; constructs linear prediction models; Best Linear Unbiased Prediction; improves upon BLUP involving kinship matrices; options: pre-specified kinships, regional kinships, adaptive multiblups, LD weightings
  * Input: PLINK format
  * Output:.reml, .indi.blp

### CNV (copy number variations) Identification:   
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
  *Output:
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

### Variant Annotation
1. **[ANNOVAR] (http://annovar.openbioinformatics.org/en/latest/)**
  * Description: command-line tool, supports SNPs, INDELs, CNVs and block substitutions, provides wide variety of annotation techniques, depends upon multiple databases (each needing to be downloaded); annotates genetic variants; utilizes RefSeq, UCSC Genes, and the Ensembl gene annotation systems; can compare mutations detected in dpSNP or 1000 Genomes Project; Very popular
   *“The final command run TABLE_ANNOVAR, using dbSNP version 138, 1000 Genomes Project 2014 Oct version, NIH-NHLBI 6500 exome database version 2 (referred to as esp6400siv2), dbNFSP version 2.6 (referred to as ljb26), dbSNP version 138 (referred to as snp138) databases and remove all temporary files, and generates the output file called myanno.hg19_multianno.txt”
  * Input: VCF, ANNOVAR input format (simple text-based format); can convert other formats into ANNOVAR input format
  * Output: VCF (if input VCF), output file with multiple columns, tab-delimited output file
1. **[wANNOVAR] (http://wannovar.usc.edu/)**
  * provides web-based access to ANNOVAR software
1. **[PolyPhen-2] (http://genetics.bwh.harvard.edu/pph2/)**
  * Description: Very popular; Polymorphism Phenotyping; Web application; predicts impact of amino acid substitution on protein; Calculates Bayes posterior probability (Last update July 2015)
  * Input: FASTA
  * Output:
1. **[SIFT] (http://sift.jcvi.org/)**
  * Description: predicts how an amino acid substitution will affect protein function; Based on degree of conservation of amino acid residues- collected though PSI-BLAST; can be applied to nonsynonymous polymorphisms or laboratory-induced missense mutations; links to dbSNP 132, GRCh37; Standalone or web app program; Very popular
  * Input: Uniprot ID or Accession, Go term ID, Function name, Species Name or ID, etc
  * Output:
1. **[snpEff] (http://snpeff.sourceforge.net/)**
  * Description: Genetic variant annotation and effect prediction toolbox; integrated with Galaxy, GATK, and GNKO; can annotate SNPs, INDELs, and multiple-nucleotide polymorphisms; categorizes effects into classes by functionality; Very popular; Standalone or Web app; Claims to calculate all SNPs in 1000 genomes (EMBI) in less than 15 minutes; can annotate SNPs, MNPs, and insertions and deletions; Provides assessment of impact of the variant ( low, medium or high)
  * Input: VCF, BED
  * Output: VCF (with new ANN field, also used in ANNOVAR and VEP), HTML summary files
1. **[SnpSIFT] (http://snpeff.sourceforge.net/SnpSift.html)**
  * Description: Filter and manipulate annotated files; Part of SnpEff main distribution; one variants have been annotated, this can be used to filter your data to find relevant variants
  * Input:
  * Output:
1. **[VAAST 2] (http://www.yandell-lab.org/software/vaast.html)**
  * Description: Variant Annotation, Analysis, and Search Tool; probabilistic search tool for identifying damage genes and the disease causing variants; can score both coding and non-coding variants; Four tools: VAT (Variant annotation tool), VST (Variant Selection Tool), VAAST, pVAAST (for pedigree data); updated April 2015
  * Input: FASTA, GFF3, GVF
  * Output: CDR (condenser file), VAAST file (both unique to VAAST)
1. **[VEP] (http://useast.ensembl.org/info/docs/tools/vep/index.html?redirect=no)**
  * Description: (Ensembl) Variant Effect Predictor; determines effect of variants on genes, transcripts, and protein sequence; uses SIFT and PolyPhen
  * Input: Coordinates of variants and nucleotide changes; whitespace- separated format, VCF, pileup, HGVS
  * Output: VCF, JSON, Statistics
1. **[ABSOLUTE] (http://www.broadinstitute.org/cancer/cga/absolute)**
  * Description: (Broad Institute); can estimate purity and ploidy to compute absolute copy number and mutation multiplicitie; reextracts data from the mixed DNA population
  * Input: HAPSEQ segdat or segmentation file
  * Output: per-sample output directory and subdirectory providing per-sample text files containing standard out being emitted from R
1. **[Alamut Batch] (http://www.interactive-biosoftware.com/alamut-batch/)**
  * Description: high-throughput annotation software for NGS analysis; for “intensive variant analysis workflows”; “enriches raw NGS variants with dozens of attributes”; based on clinically oriented Alamut database; Supports human genes; easy to integrate into pipeline (Latest Release- July 2015)
  * Input:VCF, tab-delimted file
  * Output: tab-separated file of annotations
1. **[AVIA] (http://avia.abcc.ncifcrf.gov/apps/site/index)**
  * Description: Annotation, Visualization, and Impact Analysis; “The tool is based on coupling a comprehensive annotation pipeline with a flexible visualization method. We leveraged the ANNOVAR (Wang et. al, 2010) framework for assigning functional impact to genomic variations by extending its list of reference annotation databases (RefSeq, UCSC, SIFT, Polyphen etc.) with additional in-house developed sources (Non-B DB, PolyBrowse).”
  * Input: BED
  * Output: Table of annotations with gene annotation features
1. **[BioR] (http://bioinformaticstools.mayo.edu/research/bior/)**
  * Description: (Mayo Clinic) (Page last updated June 2015) Biological Reference Repository; “data integration tool that enables coordinate based searches and joins based on strings”; “BioR consists of two parts 1) the BioR toolkit which depends on Java…. 2) the BioR catalogs which are the data files used by the system”
  * Input: VCF
   * BioR-Supported Catalogs (tar-gzip files): dbSNP, 1000 genomes, HapMap, OMIM, NCBIGene
  * Output: VCF + JSON
1. **[CADD] (http://cadd.gs.washington.edu/)**
  * Description: Combined Annotation Dependent Depletion; tool for scoring SNV deletions/insertions; “integrates multiple annotations into one metric”; Score strongly correlates with allelic diversity and pathogenicity; links to 1000 Genome variants; uses Ensembl Variant Effect Predictor
  * Input: VCF
  * Output: CADD score
1. **[CandiSNPer] (http://www2.hu-berlin.de/wikizbnutztier/software/CandiSNPer/)**
  * Description: web application, characterizes SNPs located in vicinity of SNP of interest;
  * Input: enter SNP ID (rsID), choose population, region, measure for LD, threshold plot format, color of SNPs, and chose to show genes
  * Output: Imagefile
1. **[CanvasDB] (https://github.com/UppsalaGenomeCenter/CanvasDB)**
  * Description: “local database infrastructure for analysis of targeted- and whole genome re-sequencing projects”; dependent on MySQL, R, and ANNOVAR
  * Input:
  * Output:
1. **[CHASM] (http://wiki.chasmsoftware.org/index.php/Main_Page)**
  * Description: Cancer-specific High-throughput Annotation of Somatic Mutations; Last updated May 2014; uses Random Forest Method to “distinguish between driver and passenger somatic mutations”; Positive driver class curated from COSMIC database; packed together with  SNVBox (database)
  * Input:Passenger mutation rates, Transcript and amino acid change, Genomic coordinates
  * Output: CHASM score, p-value, FDR
1. **[CRAVAT] (http://www.cravat.us/)**
  * Description: Cancer-Related Analysis of Variants Toolkit; Web application; Uses CHASM, VEST, SNVGet; “CRAVAT provides predictive scores for germline variants, somatic mutations and relative gene importance, as well as annotations from published literature and databases” Latest Release May 2015;
  * Input: VCF, CRAVAT format
  * Output: CRAVAT report- MS Excel spreadsheet or tab-separated file (emailed)
1. **[CUPSAT] (http://cupsat.tu-bs.de/)**
  * Description: Cologne University Protein Stability Analysis Tool; “tool to predict changes in protein stability upon point mutations”; web service program; Can predict mutant stability from existing PDB structures or custom protein structures
  * Input:for PDB- provide PDB ID and Amino Acid Residue Number; for custom- PDB file format
  * Output:
1. **[DANN] (https://cbcl.ics.uci.edu/public_data/DANN/)**
  * Description: Deleterious Annotation of genetic variants; standalone program, uses “the same feature set and training data as CADD to train a deep neural network”; can catch nonlinear relationships; “There are four different datasets: training, validation, testing, and ClinVar_ESP...The ClinVar_ESP dataset is also a testing set containing a set of “gold standard” pathogenic and benign variants”
  * Input: 
  * Output: 
1. **[ESEfinder] (http://rulai.cshl.edu/cgi-bin/tools/ESE3/esefinder.cgi?process=matrices)**
  * Description: Exonic Splicing Enhancer; useful for interpretation of point mutations/polymorphisms that are disease-associated; GUI interface; web app program
  * Input: FASTA
  * Output: html or plain text format, graphical display of results
1. **[Exomiser] (http://www.sanger.ac.uk/resources/software/exomiser/)**
  * Description: Wellcome Trust Sanger; functionally annotates variants from whole-exome sequencing data; Based on Jannovar and uses UCSC KnownGene; Java program; web app program (Page last modified Feb 2015)
  * Input: VCF
  * Output: TSV, VCF
1. **[FamAn] (https://sites.google.com/site/famannotation/home)**
  * Description: Automated variant annotation pipeline for family-based sequencing studies; Annotaties SNVs and INDELs; 4 models- autosomal dominant, autosomal recessive, de novo mutations and a general model; “A variety of annotations are provided for each segregating variant: number of family (and family ID) each variant hits, variant genomic location and coding effect (based on snpEff), loss-of-function mutation annotation, selected ENCODE annotation, allele frequency in the 1000 Genomes Project, allele frequency in the Exome Variant Server (ESP6500), segmental duplication annotation,  SIFT, PolyPhen2, LRT, MutationTaster, GERP++, PhyloP, SiPhy, etc.” (Last updated May 2014)
  * Input: VCF
  * Output: two excel compatible outputs
1. **[GeneTalk] (http://www.gene-talk.de/)**
  * Description: Combines tool for filtering and data analysis with an online network for genetic professionals; Different degrees- basic license, premium license, in-house solution (the last ones are paid for- Commercial tool?)
  * Input: VCF
  * Output: GeneTalk Annotation- includes clinical data, medical relevance, scientific relevance (http://www.gene-talk.de/public/GeneTalk_Whitepaper_Annotations.pdf)
1. **[GeneVetter] (http://genevetter.kidneyomics.org/)**
  * Description: “GeneVetter is a tool designed for investigation of the background prevalence of exonic variation in the Phase 3 1000 Genomes data under user defined filtering criteria”; web app program; GeneVetter uses GRch37p4 (hs37d5.fa.gz), dbSNP build 138, 1000G Phase 3, clinvar_2014072
  * Input: VCF
  * Output: TIMS score, summary table, PCA plot
1. **[GSITIC] (http://www.broadinstitute.org/software/cprg/?q=node/31)**
  * Description: (Broad Institute) Last update- July 2014; Identifies genomic regions that are significantly “amplified or deleted”; Each is given a G score; gives genomic locations and q-values from aberrant regions
  * Input: segmentation file -seg, markers file -mk (required); -array file list -alf, CNV file -cnv
   * Reference genome: -refgene (created in MATLAB, GISITIC provides four reference genomes: hg16.mat, hg17.mat, hg18.mat, hg19.mat
  * Output: All lesions file (text file), amplifications file (text file), deletion genes file (text file), Gistic Scores file, Segmented copy number (pdf file), amplification score GISTIC plot (pdf file), Deletion score/q-vale GISTIC plot (pdf file)
1. **[HOPE] (http://www.cmbi.ru.nl/hope/about)**
  * Description: Have yOur Protein Explained; Web app program; Automatic mutant analysis server that provides structural effects of a mutation; Uses BLAST against UniProt and PDB along with homology modeling
  * Input: FASTA protein sequence, or accession code of protein of interest
  * Output: a report containing information from a “decision tree” and illustrated figures and animations
1. **[Human Splicing Finder] (http://umd.be/HSF/)**
  * Description: Last update: May 2013; aimed to help study pre-mRNA splicing; combines 12 algorithms to identify mutations’ effect on splicing motifs; uses ensembl database 70
  * Input: Gene Name, Ensembl transcript ID, Ensembl Gene ID, Consensus CDS, RefSeq Peptide ID, or own sequence (looks like you can enter FASTA)
  * Output: Chart with columns for predicted signal, predicted algorithm, cDNA position and interpretation
1. **[LARVA] (http://larva.gersteinlab.org/)**
  * Description: Large-scale Analysis of Variants in noncoding Annotations; New version released July 2015; Command-line program; used for studying noncoding variants; integrates comprehensive set of noncoding elements, modeling their mutation count; Dependent on C++ and BEDtools
  * Input: multiple
  * Output:
1. **[MAC] (http://sourceforge.net/projects/mnvannotationcorrector/)**
  * Description: MNV Annotation Corrector; Ad hoc software, fixes incorrect amino acid predictions that are caused by multiple nucleotide variations; Uses existing annotators ANNOVAR, SnpEff, VEP (last update April 2015) (only 1 download this week → not popular)
  * Input: List of called SNVs and corresponding BAM
  * Output: Report identifying block of mutation within codon (BMCs)
1. **[mit-o-matic] (http://genome.igib.res.in/mitomatic/)**
  * Description: focuses on mtDNA, provides clinically relevant information from different resources; two component pipeline: command link for alignment of NGS reads and online version that provides genetic report on mitocondrial variants
  * Input:FASTQ, pileup
   * Reference sequence: rCRSm
  * Output: Online version gives comprehensive genetic report
1. **[Mutadelic] (http://krauthammerlab.med.yale.edu/mutadelic/index.html)**
  * Description: Web App program; “This application generates reports on inherited mutations in five genes (ANK1, SLC4A1, SPTA1, SPTB and EPB42) associated with the following rare Mendelian blood disorders: Hereditary Spherocytosis (HS), Hereditary Elliptocytosis (HE) and Hereditary Pyropoikilocytosis”; Newer program- recently validated on omictools
  * Input: Can upload coordinates of DNA variants or VEP
  * Output: Displayed on web or can be downloaded in Excel or RDF format
1. **[MutationTaster] (http://www.mutationtaster.org/)**
  * Description: (Last post on site 2014) Web app program; Rapid evaluation of disease causing alterations; uses NCBI 37 and Ensembl 69
  * Input: HGNC symbol, NCBI GeneID, or Ensembl ID,
  * Output: Report containing prediction, summary, name of alteration, etc
1. **[MutPred] (http://mutpred.mutdb.org/)**
  * Description: web app tool; Classifies amino acids substituation as disease associated or neutral in humans; Last modified Feb. 2014; Based on SIFT, trained using Human Gene Mutation Database
  * Input:
  * Output: “The output of MutPred contains a general score (g), i.e., the probability that the amino acid substitution is deleterious/disease-associated, and top 5 property scores (p), where p is the P-value that certain structural and functional properties are impacted.”
1. **[MutSigCV] (http://www.broadinstitute.org/cancer/cga/mutsig)**
  * Description: (Broad Institute) Mutation Significance (CV= covariates); Analyzes mutations discovered in DNA sequencing to identify genes that were mutated more often than expected
  * Input: mutations.maf, coverage.txt, covariates.txt
  * Output: output.txt
1. **[NGS-SNP] (http://stothard.afns.ualberta.ca/downloads/NGS-SNP/)**
  * Description: Collection of command-line scripts for providing rich SNP annotations; “NCBI, Ensembl, and Uniprot IDs are provided for genes, transcripts and proteins when applicable”;
  * Input: Samtools consensus pileup, Maq, diBayes, Genetic format, VCF
  * Output: File containing annotated SNPs is copied from SNP list and some classes are added
1. **[Oncotator] (http://www.broadinstitute.org/oncotator)**
  * Description: (Broad Institute) “Tool for annotating human genomic point mutations and data relevant to cancer researchers”; Web app; Supports annotation of data from ClinVar, dbSNP, 1000 genomes (plus many other external sites); Only GRCh27 coordinates supported; Last update: April 2015
  * Input: tal-delimited file
  * Output: tab-delimited MAF
1. **[PANTHER] (http://omictools.com/panther-s649.html)**
  * Description: Protein ANalysis THrough Evolutionary Relationships; Web app program, also has its own database; Classification system used to classify proteins and their genes; Also, “Estimates the likelihood of a particular nonsynonymous (amino-acid changing) coding SNP to cause a functional impact on the protein”; Updated in 2015
  * Input: Data from PANTHER, IDs from Ensembl, EntrezGene, NCBI GI numbers, NCBI UniGene IDs HUGO, UniProt; if ID type is not one of the above, can input txt file or excel format
  * Output: Analysis results displayed online
1. **[PESX] (http://cubio.biology.columbia.edu/pesx/pesx/)**
  * Description: Putative Exonic Splicing Enhancers/Silencers; (Can’t tell if this is outdated or not)
  * Input: FASTA or plain text
  * Output: Excel spread sheet
1. **[Phen-Gen] (http://phen-gen.org/index.html)**
  * Description: Combines patient's’ disease symptoms with sequencing data; Standalone or Web app version; Only excepts 1 family per run, in order to evaluate unrelated individuals, each sample needs to be run individually
  * Input: Variant- VCF; Pheotype- HPO; Pedigree- PED
  * Output: Combined scores file, variants for top genes file
1. **[PMUT] (http://mmb.pcb.ub.es/PMut/)**
  * Description: Aimed at annotation and prediction of pathological mutations; based on different kinds of sequence info and neural networks to process information
  * Input: FASTA
  * Output; Simple yes/no and reliability index
1. **[PROVEAN] (http://provean.jcvi.org/index.php)**
  * Description: Protein Variation Effect Analyzer; predicts whether an amino acid substitution or indel has impact on biological function of the protein; “comparable to SIFT or Polyphen-2”; Standalone, Web app, Command line or GUI; Last update May 2014
  * Input: FASTA, list of variants;
  * Output: tab-separated columns including Variant, Provean Score and prediciton
1. **[Rescue-ESE] (http://genes.mit.edu/burgelab/rescue-ese/)**
  * Description: “An online tool for identifying candidate ESEs in vertebrate exons”; Web application; For human, mouse, zebrafish, pufferfish
  * Input: multi-FASTA or plain text
  * Output:
1. **[SCAN] (http://scandb.org/newinterface/index_v1.html)**
  * Description: Web application program, includes a database as well; Database contains physical-based SNP annotations and functional annotations; “Information on physical, functional, and LD annotation served on the SCAN database comes directly from public resources, including the HapMap (release 23a), NCBI (dbSNP 129), or is information created by us using data downloaded from these public resources”; “SCAN can be utilized in several ways including: (i) queries of the SNP and gene databases; (ii) analysis using the attached tools and algorithms; (iii) downloading files with SNP annotation for various GWA platforms”
  * Input:
  * Output: HTML, comma-delimited, tab-delimited
1. **[SeattleSeq Annotation] (http://snp.gs.washington.edu/SeattleSeqAnnotation137/)**
  * Description: “SeattleSeqAnnotation137 was most recently updated October 13, 2013. The current version is 8.08. The most recent site, based on dbSNP build 141, and hg38/NCBI 38”; Provides annotations for SNVs and Indels- includes dbSNP rsID, gene names and accession numbers, variation functions, protein positions and amino acid changes, conservation scores, HapMap frequencies, PolyPhen predictions and clinical association.
  * Input: Maq, gff, CASAVA, VCF, GATK bed, custom
  * Output: “default output file format is a header line (starting with "#") followed by tab-separated annotations”; VCF
1. **[seqminer 3.7] (https://cran.r-project.org/web/packages/seqminer/)**
  * Description: “Efficiently Read Sequence Data (VCF Format, BCF Format and METAL Format) into R”; Command line package program; Published August 2015
  * Input: VCF, BCF
  * Output: VCF
1. **[SG Adviser] (https://genomics.scripps.edu/ADVISER/Home.jsp)**
  * Description: Scripps Genome Annotation and Distributed Variant Interpretation Server, web developed applications for variant annotation, “Downstream applications of variant annotation include: Clinical sequencing applications including: carrier testing, or identification of causal variants in molecular diagnosis, tumor sequencing, or diagnostic odyssey. Prioritization of variants prior to statistical analysis of sequence based disease association studies, especially for automated set-generation and enrichment of likely functional variants within sets. Identification of causal variants in post-GWAS/linkage sequencing studies. Identification of causal variants in forward genetic screens (stay tuned for non-human annotation)”
  * Input: SNV- VCF, BED, and a few others; CNV- BED, CNVator, plus others
  * Output: tab-delimited file
1. **[SNAP-2] (https://rostlab.org/services/snap/)**
  * Description: “SNAP2 is a trained classifier that is based on a machine learning device called "neural network". It distinguishes between effect and neutral variants/non-synonymous SNPs by taking a variety of sequence and variant features into account”; predicts impact of amino acid substitution on protein
  * Input:Protein Sequence in FASTA
  * Output: “Each model outputs one score for each output class (neutral/effect). These scores of 10 models are averaged in a jury decision. The final score is calculated as the difference of the average score for effect and the average score for neutral’
1. **[SNiPA] (http://snipa.helmholtz-muenchen.de/snipa/)**
  * Description: “SNiPA offers both functional annotations and linkage disequilibrium information for bi-allelic genomic variants (SNPs and SNVs)”; Based on EMBI 1000 Genomes Project; uses GRCh37 and Ensembl
  * Input: rs identifiers
  * Output: list of SNiPA cards
1. **[SNPAAMapper] (http://www.ccmb.med.umich.edu/ccdu/SNPAAMapper)**
  * Description: “SNPAAMapper is a downstream variant annotation program that can effectively classify variants by region (e.g. exon, intron, etc), predict amino acid change type (e.g. synonymous, non-synonymous mutation, etc), and prioritize mutation effects (e.g. CDS versus 5’UTR, etc).”
  * Input: VCF in tab-delimited format
  * Output: spreadsheet result file
1. **[SNPMeta] (http://www.tc.umn.edu/~konox006/Code/SNPMeta/)**
  * Description: Generates metadata for SNPs, including gene name, whether the SNP is coding or noncoding, and whether the SNP is synonymous or nonsynonymous; last update October 2013
  * Input: FASTA
  * Output: dbSNP submission report format or tab-delimited file
1. **[SNPnexus] (http://www.snp-nexus.org/)**
  * Description: SNPnexus allows single queries using dbSNP identifiers or chromosomal regions for annotating known variants. The users are also allowed to provide novel in-house SNPs/indels using genomic coordinates on clones, contigs and chromosomes. For practical purposes, SNPnexus allows batch queries comprising SNP data using dbSNP identifiers or genomic coordinates; Web app, also has its own database; Last update May 2014;
  * Input: genomic position, chromosomal region, dbSNP id, tab-delimited, VCF
  * Output: Tables for genomic mapping, gene/protein consequences, effect on protein function, hapmap population, regulatory elements, conservation, phenotype and disease association, structural variants
1. **[SNPs&Go] (http://snps-and-go.biocomp.unibo.it/snps-and-go/)**
  * Description: Predicts human disease-related mutations in proteins with functional annotations; web application; collects framework information;
  * Input: UNIPROT accession number, mutation postition, wild-type residue, substituting residue
  * Output:
1. **[TRAMS] (http://figshare.com/articles/Tool_for_rapid_annotation_of_microbial_SNPs_TRAMS_a_simple_program_for_rapid_annotation_of_genomic_variation_in_prokaryotes_/782261)**
  * Description: Tool for Rapid Annotation of Microbial SNPs; functional annotation of genomic SNPs, annotates as synonymous, nonsynonymous or nonsense; standalone program; have developed workflow in Galaxy
  * Input: tab-delimited file containing SNP locations, reference nucleotide and SNPs in different strands
   * Reference genome: Genbank or EMBL format
  * Output:  
1. **[Trinotate] (http://trinotate.github.io/)**
  * Description: Transcriptome functional annotation and analysis; “Trinotate makes use of a number of different well referenced methods for functional annotation including homology search to known sequence data (BLAST+/SwissProt/Uniref90), protein domain identification (HMMER/PFAM), protein signal peptide and transmembrane domain prediction (singalP/tmHMM), and comparison to currently curated annotation databases (EMBL Uniprot eggNOG/GO Pathways databases)”
  * Input: trinity.fasta, BLAST homologies (confusing workflow, don’t really understand, has a ton of dependencies)
  * Output: Annotation report in .xls
1. **[VAGrENT] (http://www.sanger.ac.uk/resources/software/vagrent/)**
  * Description: (Wellcome Trust Sanger) Variation Annotation GENeraTor; Suite of perl modules, compares variants with reference genome annotations and generates effects of each variant; Last modified Feb 2014
  * Input:
  * Output:
1. **[VARIANT] (http://variant.bioinfo.cipf.es/)**
  * Description: VARIant ANalysis Tool; reports properties of any variant in human, mouse or rat genomes; Requires log in
  * Input: VCF, GFF, BED
  * Output: Filtered Variants (same as input); Genes with Variants (txt), Consequence type histogram (txt)
1. **[VariantAnnotation] (http://www.bioconductor.org/packages/2.12/bioc/html/VariantAnnotation.html)**
  * Description: (Bioconductor) Annotates and filters genetic variants; data is read in VCF and variants are identified according to region (conding, intron, intergenic, spice site, etc); Uses SIFT and PolyPhen to provide protein function predictions; uses 1000 Genomes Chromosome 22 for sample data
  * Input; VCF
  * Output:
1. **[Variobox] (http://bioinformatics.ua.pt/software/variobox/)**
  * Description: Desktop tool for annotation, analysis and comparison of human genes; Standalone program; Annotation data obtained from wave, PDB, or UniProt, LRG, or RefSeq; Researched genes compared to sequences for LRG and RefSeq;
  * Input: DNA Sequence chromatogram file (.scf, .abi); DNA electropherogram file (.ab1); FASTA
  * Output: “Gene Panel” (looks like variation viewer)
1. **[VarioWatch] (http://genepipe.ncgm.sinica.edu.tw/variowatch/main.do)**
  * Description: annotations of genomic variants; Last update February 2015; uses dbSNP 141, annotation release 105;
  * Input: gene name, SNP ID, position
  * Output: specially designed views
1. **[VAT] (http://vat.gersteinlab.org/)**
  * Description: Variant Annotation Lab; Updated May 2014,
  * Input:VCF, Interval
  * Output: Gene and sample summaries
1. **[WhopGenome] (https://cran.r-project.org/web/packages/WhopGenome/index.html)**
  * Description: fast access to whole genome, population scale variation data, Links to UCSC Bioconductor, and AmiGO (for pedigree formats)
  * Input: VCF, FASTA, Phylip, MAF, plus others
   * Reference sequence: FASTA
  * Output:
