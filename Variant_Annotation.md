## Variant Annotation
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
  <hr/>
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
1. **[CAROL] (http://www.sanger.ac.uk/resources/software/carol/)**
  * Description: (Wellcome Trust Sanger); Combined Annotation scoRing toOL; Combined functional annotation score of nonsynonymous coding variants; Combines information from PolyPhen-2 and SIFT
  * Input: tab-delimited with columns obtained from PolyPhen-2 and SIFT output
  * Output: tab-delimited file
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
1. **[LINKAGE] (http://www.jurgott.org/linkage/LinkagePC.html)**
  * Description:three main programs: mlink (calculates lod scores at fixed values for the recombination fraction in one interval of a genetic map), linkmap (calculates location scores for positions of a disease locus along a marker), and ilink (estimates parameters including recombination fractions, allele frequencies, penetrances, etc)
  * Input: pedfile (processed by MAKEPED) and datafile (reflects loci for each individual; set in PREPLINK)
  * Output:
1. **[MAC] (http://sourceforge.net/projects/mnvannotationcorrector/)**
  * Description: MNV Annotation Corrector; Ad hoc software, fixes incorrect amino acid predictions that are caused by multiple nucleotide variations; Uses existing annotators ANNOVAR, SnpEff, VEP (last update April 2015) (only 1 download this week ? not popular)
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
