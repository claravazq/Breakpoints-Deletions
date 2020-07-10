# Breakpoints and deletions
Learn how to make breakpoints plot using BED files.
BED files can come to you in two ways:

  1.- With a column containing something like: chrX:zzzz-kkkk
  
  2.- Treating the coordinates like blocks and having cells with the begining and the end of the blocks separated by commas
  
In this script, both ways are considered and resolved to use in R. 

Plots resembling the IgH locus and the breakpoints as histogram are then written:
- For the whole IgH locus 
- For a specific region, used to zoom into different switch regions
