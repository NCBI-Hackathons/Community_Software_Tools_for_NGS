# List of Community Software Tools for Next Generation Sequencing
## Commercial tools (not in any particular order):
1. **[Strand NGS] (http://www.strand-ngs.com/)**<br>
  * offers many different tools including alignment, RNA-Seq, DNA-Seq, ChIP-Seq, Small RNA-Seq, Genome Browser, visualizations, Biological Interpretation, etc. Supports workflows
“one can import the sample data in FASTA, FASTQ or tag-count format. In addition, prealigned data in SAM, BAM or Illumina-specific ELAND format can be directly imported for analysis.” <br>
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
  *“perfect analytical partner for the analysis of desktop sequencing data produced by the ION PGM™, Roche Junior, Illumina MiSeq as well as high throughput systems as the Ion Torrent Proton, Roche FLX, Applied BioSystems SOLiD™ and Illumina® platforms.” runs on Windows, free-standing multi-application package-- SNP/Indel analysis, CNV prediction and disease discovery, whole genome alignment, etc.  
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

