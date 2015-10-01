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
1. **[HapMuc] (https://github.com/usuyama/hapmuc)
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
1. **[MutationSeq] (http://compbio.bccrc.ca/software/mutationseq/)
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
1. **[RVD2] (http://genomics.wpi.edu/rvd2/)
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
1. **[Triodenovo] (http://www.pitt.edu/~wec47/triodenovo.html)
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
1. **[MSIsensor] (https://github.com/ding-lab/msisensor)
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
1. **[VariantMaster] (http://sourceforge.net/projects/variantmaster/)
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
