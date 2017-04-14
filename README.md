# SiteOut: a tool to design motif-free DNA sequences

SiteOut allows to remove specific nucelotide motifs, such as trancription factor binding sites, from a DNA sequence. It can be used to design a new sequence from scratch, to refine a predefined sequence, or to create neutral spacers between functional sequences. The code will either look for explicit motifs or will use Patser to check for predicted motifs based on the provided frequency matrices (FMs) -collections of FMs ready to be used with SiteOut for D. melanogaster and yeast transcription factors are available in this repository inside the Dmel_pwms and yeast_pwms zipfiles respectively-.

To run SiteOut locally you need to compile Patser (http://stormo.wustl.edu/resources.html) and give it execution permits. Once you've done that, you will have to edit line 147 in SiteOut.py to indicate your path to Patser

os.system("your_path_to_patser/patser-v3b -w -v -p -f fileName -c -l %.4f > %s " % (cutoff, tempFile))
