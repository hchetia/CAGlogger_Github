# CAGLogger

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
A frequency CSV file all CAG repeats > 15 and their frequencies [*Adjust for 3 repeats*]. DONE (gotta minus 3)     
A frequency distribution plot for > 15. DONE     
A frequency distribution plot for > 35. DONE     
A metric file with percentage of extremely long CAG repeats. NOT DONE.   

**F. References:**
Targeted Miseq data from HTT exon1 CAG repeat (Strategy- [Ciosi et al 2021](https://content.iospress.com/articles/journal-of-huntingtons-disease/jhd200433), [Matlik et al 2023](https://www.biorxiv.org/content/10.1101/2023.04.24.538082v2.abstract)).


Next steps will be improving the look of the output graph, the metric report in csv and accelerating the program's runtime (Ivy).
