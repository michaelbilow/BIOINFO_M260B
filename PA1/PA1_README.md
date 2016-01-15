# Programming Assignment 1


## Getting Started
If you've set everything up correctly, you should have no trouble running
*basic_aligner.py* and *basic_pileup.py* to get some properly formatted output
which you can submit on the [course site](https://cm124.herokuapp.com)

## Instructions
Don't worry about insertions and deletions for this project. All you need to do is figure out how to tell the true SNPs from the false positives in the consensus sequence.

## Outline of provided scripts
The reference and reads are the inputs to *basic_aligner.py*, which are converted to the aligned file using a trivial alignment algorithm.

The aligned data is fed into the *basic_pileup.py* script, which generates a consensus sequence by picking the most common base at each position.  That "consensus" file is the heart of this assignment; it has the reference, the aligned reads, and then an asterisk at every position where the consensus sequence differs from the reference. If you can understand what is going on there, you'll understand what you can change to improve mapping true SNPs, and then go on to mapping structural variations.

The *basic_pileup.py* script also makes a file that starts snps that is used to format the output properly.  It notes all of the differences between the consensus and the reference and notes the position where they occur. It also zips that file so you can sumbit it directly to the course site.

## What You Should Change


## If you need more help

### Thinking about the biology

.1% of your genome is SNPs;