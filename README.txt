KFERQ_searcher

Version history

0.4	2018-02-18
Starting version with fully implemented search from data base
choose between organisms

0.5	2018-02-20
search in an unknown amino acid sequence

0.6	2018-02-21
released version
enables displaying identified motifs in color if a sequence is entered

0.61	2018-02-22
convert all entered sequences to upper case to find motifs even if lower case characters are entered
also implemented for UniProt identifiers

0.64	2018-02-23
changed color of canonical motifs from "yellow" to "gold" for better readability
Included handling of FASTA formats (headers) and possible whitespace in the sequence
CSS was added to the UI function to allow for line breaks in long words (i.e. protein sequences)

Starting from this point Git version control is enabled

0.65	2018-04-15
To make launching the website easier I added a "launcher.R" script that will take care of loading the required packages. This could also be used to check for files and such
