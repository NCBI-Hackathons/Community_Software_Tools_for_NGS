## Tool Repositories, Analytical Pipelines, & Workflow Systems
### ToolKits
1. **[bamtools] (https://github.com/pezmaster31/bamtools)**
  * Description: command line tool kit, also API version; for reading, writing and manipulating BAM (genome alignment) files
  * Input: BAM
  * Output: BAM
1. **[Bioconductor] (http://www.bioconductor.org/)**
  * Description: Open source software project, Popular R module library with over 1000 modules, can analyze diverse high-throughput genomic data; options for multiple different workflows; also available on AMI
  * Input: Multiple
  * Output: Multiple
1. **[Galaxy Tool Shed] (https://toolshed.g2.bx.psu.edu/)**
  * Description: Open, web-based platform; Has 3374 valid tools for many different purposes. Looks like this site has a different tool repository for each different data type, They have agreed to support NCBI genomes, and SRA on a limited basis
  * Input: Different tools for different types of data
  * Output: Multiple
1. **[GATK] (https://www.broadinstitute.org/gatk/)**
  * Description: (Broad Institute) Genome Analysis ToolKit; software library including tools (depth of coverage analyzers, quality score recalibrator, local realigner, SNP/INDEL caller), requires Ant and Java, extensive documentation and wiki system, active community, can also be used to identify somatic mutations
  * Input: Multiple
  * Output: Multiple
1. **[picard] (http://broadinstitute.github.io/picard/)**
  * Description: (Broad Institute) Java based command line tools for manipulating high throughput sequencing data in the BAM format; supported through GATK
  * Input: BAM, SAM, VCF
  * Output: BAM, SAM
1. **[SAMtools] (http://samtools.sourceforge.net/)**
  * Description: provides utilities for manipulating alignments in the SAM format, including sorting, merging, indexing, and generating alignments; can be used for calling SNPs and Indels
  * Input: Conforms to the specifications produced by GA4GH File Formats working group
  * Output:
1. **[VCFtools] (https://vcftools.github.io/specs.html)**
  * Description: set of tools written in Perl and C++ for working with VCF files
  * Input: VCF, GZVCF, BCF
  * Output: Multiple, can be piped into another program
