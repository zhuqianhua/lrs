# lrs
**Long Reads Segmentation**
## Introduction
lrs is a software for Long Reads segmentation. The HTSlib  
was neccessary, which was used to deal the bam files.  
The primers against the reads was performed by using  
Smith-Waterman.
## Getting started
	git clone https://github.com/zhuqianhua/lrs.git  
	modify to the correct path of htslib in Makefile
	make      
## Parameters
Options | Type| Description
--- | --- | ---
-i,--input|*<s>|input file *.bam or *.fa  
-f, --format|*<s>|input file format, b: bam and f: fasta  
-o, --out|<s>|prefix of output, default ./lrs  
-u, --umi|<i>|length of umi, default 8  
-a, --primer5|<s>|5' primer sequence, default AAGCAGTGGTATCAACGCAGAGTACGGGGG  
-b, --primer3|<s>|3' primer sequence, default AAGCAGTGGTATCAACGCAGAGTAC  
-G, --pm5g|<b>|5' primer with GGGGG and 3' primer without GGGGG  
-L, --cln|<i>|minimum continuously align length of primer, default 5  
-l, --len|<i>|minimum align length for primer, default 16  
-p, --process|<i>|number of threads, default 6  
-h, --help|<b>|print this information

*note * indicates required parameters*
## Output


## Notes
It relies on multithreaded execution and your system needs  
to support -lpthread. If you have questions about lrs, you  
may send the questions to zhuqianhua@bgi.com.

## Reference
1)	http://bioinformatics.oxfordjournals.org/content/early/2016/12/21/bioinformatics.btw741.full.pdf+html