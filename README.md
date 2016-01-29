# petcat-rt
###Codes
**cevt.c**  
Usage: ./cevt.c or ./cevt -c   
Synthesize two waveforms of different data file.  
-c option enables CFD based timing adjustment  

**cevt_list.c**  
Usage: ./cevt_list.c -l runlist.txt  
Do same procedure as cevt.c according to runlist.txt.  
Synthesize the waveforms with all combinations.  
-c option enables CFD based timing adjustment  
-> Output will be wrote in ./cevtout/  


**vegcat.c**  
Usage: ./vegcat -l runlist.txt  
Reconstruct the interaction position by xi square fitting.  
-> Output will be wrote in ./vegcatout/  
Comment: On Jan26, the mapfile is modified. Now detmap_Q4pos4_CC09_Ecal.txt is used in vegcat.h to apply the appropriate Ecal. But it seems not working well.  
  

###Root macros
**tr2hist.C**  
Convert tr file createc by cevet.c to waveform graph, make an output PDF file.

**cfdoffsetdist.C**
Make time histogram of the CFD offset result.
