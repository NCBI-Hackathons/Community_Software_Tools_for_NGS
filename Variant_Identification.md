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
