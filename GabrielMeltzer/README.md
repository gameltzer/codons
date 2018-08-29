#### Name: Gabriel Meltzer
#### Date: Originally written Novemeber 2016, small revisions August 29, 2018
#### Python version Python 2.7.14 (may have issues when running in Python 3)

##Description

Asks the user for the path of a FASTA or FASTQ file containing nucleotide sequecnes. It then returns a count of the individual amino acids found in each record. It also includes a record of the amnino acids found after converting the DNA sequence to RNA. (This should be identical.) It returns an error if the file entered is not in the proper format.

##Additional Notes
Python 2.7 is recommended. In addition, Biopython will need to be installed to get the code to work.

To install biopython,

`pip install biopython`