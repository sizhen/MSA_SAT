# MSA_SAT
solving multiple sequence alignment instances using pySMT

## Dependencies
python = 3.7.4
pySMT = 0.9.0

## RUN
mas01.py can ran on the multiple sequences over binary alphabet (0 and 1) with exactly two ones in each sequence.  
``python msa01.py 0100000100 000010001 0001000100``

masRNA.py can ran on the multiple sequences over four letters (A, G, U, C).  
``python msaRNA.py "GAUCAGUC" "AUCAGC" "GACGAUAAAA"``