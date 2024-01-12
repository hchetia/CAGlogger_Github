# CAGLogger

**Broad Goal:**
Taking the HTT genotyping based on the Ciosi Assay one step further by detecting the detecting CAG repeats that are extremely long (> 113). We use an in-house alignment free detection program and tally CAG repeats from targeted MiSeq data from genomic DNA. 

**A. Define Type A Error**: 
1 base substitution per one CAG repeat (concept of mismatches or gap not applied here).  
**B. Define Type B Error:** 
Two consecutive CAG repeats with Type A Error.  
**C. Error Cut-off**:   
Maximum Type A Error allowed per CAG stretch is 3. Maximum threshold for Type B Error per CAG stretch is 1.  DONE  
Try 5 and 7. Compare % of reads over 100= Total reads with CAG over 100/Total reads with CAG over 35 * 100.  
**D. What are the inputs?**
FastQ files (generated using the protocol in References Section).  
**E. What are the outputs?**
A frequency file of CAG repeats and their frequencies [*Adjust for 3 repeats*]. DONE (gotta minus 3)
A frequency distribution plot. DONE
A metric file with percentage of extremely long CAG repeats. DONE  
**F. References:**
Targeted Miseq data from HTT exon1 CAG repeat (Strategy- [Ciosi et al 2021](https://content.iospress.com/articles/journal-of-huntingtons-disease/jhd200433), [Matlik et al 2023](https://www.biorxiv.org/content/10.1101/2023.04.24.538082v2.abstract)).


DONE!

Testing with the 7094 sample (which has known no. of CAGn) and using Kert's suggestion about percentage calculation:  
% of reads over 100 is 0.09% for 3 errors, 0.4 % for 5 errors and 0.2% for 7 errors.  
The code is working well. Will settle for one of these error-cutoffs after today's discussion.

Moving ahead with 5 and 7 type A errors with selected samples, going to compare them with the alignments seen in those samples and settle on the cut-off of this error type.

Next steps will be improving the look of the output graph, the metric report in csv and accelerating the program somehow (Ivy?).

