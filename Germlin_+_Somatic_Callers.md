## Variant Identification
### Germline + Somatic  Callers
1. **[VarScan 2] (http://massgenomics.org/varscan)**
  * Description: identify germline variants, private and shared variants, somatic mutations, and somatic CNVs; detects indels
  * Input: SAMtools pileup
  * Output: VCF
  <hr/>
1. **[BAYSIC] (http://genformatic.com/baysic/)**
  * Description: Bayesian method; combines variant calls from different methods (GATK, FreeBayes, Atlas, Samtools, etc)
  * Input: VCF format from one or more variant calling programs
  * Output: VCF file containing integrated set of variant calls
1. **[MSIsensor] (https://github.com/ding-lab/msisensor)**
  * Description: Microsatellite instability detection; C++ program, detects somatic and germline variants in tumor-normal paired data
  * Input: BAM index files (normal and tumor)
  * Output:
1. **[Beagle version 4] (http://faculty.washington.edu/browning/beagle/beagle.html)**
  * Description: software package: genotype calling, phasing, imputation of ungenotyped markers, and identity-by-descent segment detection:unsure if this one is in the right category; genotype calling, phasing, imputation of ungenotyped markers, and identity-by-descent segment detection; 
  * Input: VCF
  * Output: VCF
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
