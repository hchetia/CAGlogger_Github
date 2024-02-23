# CAGLogger

<img width="662" alt="image" src="https://github.com/hchetia/CAGLogger/assets/64626735/76b39313-8ffd-4e03-bf76-8020ca1f3a8f">

**Broad Goal:**
Taking the HTT genotyping based on the Ciosi Assay one step further by detecting the detecting CAG repeats that are extremely long (> 113). We use an in-house alignment free detection program and tally CAG repeats from targeted MiSeq data from genomic DNA. 
**A. Define Type A Error**: 
1 base substitution per one CAG repeat (concept of mismatches or gap not applied here).  
**B. Define Type B Error:** 
Two consecutive CAG repeats with Type A Error.  
**C. Error Cut-off**:   
Maximum Type A Error allowed per CAG stretch is 5. DONE   
Maximum threshold for Type B Error per CAG stretch is 1.   CHECK 
% of reads over 110= Total reads with CAG over 110/Total reads with CAG over 35 * 100.  
**D. What are the inputs?**
FastQ files (generated using the protocol in References Section).  
**E. What are the outputs?**
A frequency CSV file all CAG repeats > 15 and their frequencies [*Adjust for 3 repeats*]. DONE
A frequency distribution plot for > 35. DONE     
A metric file with percentage of extremely long CAG repeats. DONE.   
**F. References:**
Targeted Miseq data from HTT exon1 CAG repeat (Strategy- [Ciosi et al 2021](https://content.iospress.com/articles/journal-of-huntingtons-disease/jhd200433), [Matlik et al 2023](https://www.biorxiv.org/content/10.1101/2023.04.24.538082v2.abstract)).
~~Next steps will be reducing the program's runtime per sample.~~  
_Catch- _Our in-house and published alignment-based CAG genotyping studies have revealed that the most common right flank structure of mutCAG allele is ....(CAG)nCAACAGCCG_CCA_CCGCCG. The algo halts counting the structure at CCA as it has two base substitutions within one triplet codon, thereby calculating 3 CAGs less in each allele structure reported. This issue has been appropiately addressed in our code.   

**NEXT in 2024-**   
~~1. Get the conda env and slurm+R script together.~~  
2. Share with Nick.  
~~3. Get the trgt database path repeats (31 repeats atm)~~  
4. Select the repeats that have an identical assay (with open source geo datasets) (FEB END)  
5. Dig into literature if there are any new path repeats that aren't in these 31 ref.  
6. Create functions (MARCH)  
7. Download GEO raw data and replicate the results. (MARCH)  
8. Create the shiny app OR a viewer (something like REVIEWER Maybe)?? (APR)  
9. Users should be able download some format of visualisation to test with alignment data. (APR)  
10. Write the package. Collab with BRC!?   
11. Draft of paper (Bioinformatics journal).  
12. Make it squeaky clean, bug-free and PUBLISHHHHH in a free bioinfo journal. (Ivy got a headstart on this. Draft ready for submission by June)
