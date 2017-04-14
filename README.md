# SiteOut: a tool to design motif-free DNA sequences

SiteOut allows to remove specific nucelotide motifs, such as trancription factor binding sites, from a DNA sequence. It can be used to design a new sequence from scratch, to refine a predefined sequence, or to create neutral spacers between functional sequences. The code will either look for explicit motifs or will use Patser to check for predicted motifs based on the provided frequency matrices (FMs) -collections of FMs ready to be used with SiteOut for D. melanogaster and yeast transcription factors are available in this repository-.
