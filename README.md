# G4PromFinder

G4PromFinder is an algorithm for the promoter prediction in intergenic regions of a bacterial genome. It is recommended for GC rich genomes. G4PromFinder predicts putative promoters based on AT-rich elements and G-quadruplex DNA motifs.

We created 2 versions of G4PromFinder: 
- G4PromFinder (Method 1) that consider all the elements identified in a determined intergenic region;
- G4PromFinder (Method 2) that consider only the predict elements closest the coding sequence.


Input: Genome file and annotation file.
Genome file must be in fasta format.
Annotation file must be in txt format. It must be organized in 3 columns. The first column must contain CDS start, the second CDS end and the third the strand (+1 for strand+ and -1 for strand-).


Output: a text file containing promoter coordinates and a file containing informations about them.

