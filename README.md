# CAGLogger

<img width="662" alt="image" src="https://github.com/hchetia/CAGLogger/assets/64626735/76b39313-8ffd-4e03-bf76-8020ca1f3a8f">

**Broad Goal:** 
Detecting CAG/CTG repeats that are extremely long (> 110) from targeted amplicon sequencing. The tested loci are CAG repeat loci in HTT exon1 and CTG repeat loci in ATXN3 genes . We use an in-house alignment free detection script with BWA and tally CAG/CTG repeats from targeted MiSeq data from genomic DNA.  
**A. Define Type A Error**: 
1 base substitution per one CAG repeat (concept of mismatches or gap not applied here).  
**B. Define Type B Error:**  
Two consecutive CAG repeats with Type A Error.  
**C. Error Cut-offs**:   
Maximum Type A Error allowed per CAG stretch is 5. Maximum threshold for Type B Error per CAG stretch is 1.  
**D. Inputs to CAGlogger**  
FastQ files (generated using the protocol in References Section).  
**E. Outputs of CAGlogger**  
A frequency CSV file all C[A/T]G repeats > 35; A metric file with percentage of extremely long CAG repeats.  
**F. Percentage of reads over 110**  
(Total reads with CAG over 110/Total reads with CAG over 35) * 100.  

**Note**  
[Our in-house and published alignment-based CAG genotyping studies have revealed that the most common right flank structure of mutCAG allele is ....(CAG)nCAACAGCCG_CCA_CCGCCG. The algo halts counting the structure at CCA as it has two base substitutions within one triplet codon, thereby calculating 3 CAGs less in each allele structure reported. This issue has been appropiately addressed in our code.]  

**References:**  
[Ciosi et al 2021](https://content.iospress.com/articles/journal-of-huntingtons-disease/jhd200433)  
[Matlik et al 2023](https://www.nature.com/articles/s41588-024-01653-6).
