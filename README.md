primeScaff
==========

primeScaff is intended to automate the sometimes tedious process of manually designing specific primer pairs around gaps of genomic scaffolds and speedup the genome finishing stage of a genome sequencing project. It incorporates de-novo repeat finding using RECON to avoid as much as possible designing primers in repetitive regions and offers the possibility to easily fine-tune the primer design options using primer3. It outputs the repeat, gap and primer annotations in gff2 and gff3 to easily load and visualize the results using tools like UGENE (http://ugene.unipro.ru/) (recommended) or artemis (http://www.sanger.ac.uk/resources/software/artemis/).

## Features
* Iterative de novo repeat masking using iteratively RECON
* Automatic design of primers around hundreds of gaps in scaffolds in seconds
* Can be used in any assembly pipeline, just need a standard FASTA file
* Solution for de novo assembly projects with no close reference
