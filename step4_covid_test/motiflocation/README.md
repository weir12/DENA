# motiflocation
A small program to find motif locations in FASTA files. 

Usage
-----
To compile, you should just be able to type 'make'.

Usage is as follows:

`./motiflocation -file1 <fasta_file> -startstr <motif> -out <out_file>`

+ Parameters
  - -file1 : Specify the input fasta file.
  - -startstr : Specific the motif to find in the fasta file.
  - -out : Specify the output file
  
+ Optional parameters
  - -chr : Specify the chromosome on which to run the algorithm.
  - -mutations : Specify the number of mutations (mismatches) to allow when finding motifs. 
  - -uppercase_only : Only report uppercase matches.
  - -lowercase_only : Only report lowercase matches.

Output
------

The program produces an output file with 5 columns.

1. Chromosome
2. Motif start position
3. Motif end position
4. Strand
5. Motif sequence
