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

0.8	2018-07-30
Major update of the ui and server code. Data is now loaded on demand to increase performance. Motifs from precomputed tables and new sequences can now be downloaded as tables. Multiple sequences can be pasted or uploaded from file. Colored text output is moved to a separate tab in the from sequence window. Complete redesign of the interface to hopefully improve the usability and make the option of batch processing easier to see. The "more" tab is replaced with an "about" tab containing a brief rundown of the available functions. Styling information is moved to a separate CSS file ("styles.css").
This version is intended to be loaded on its own and not through the launcher (additional packages that are required now are "readr" and "shinyjs"). 