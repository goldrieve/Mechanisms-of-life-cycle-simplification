# Mechanisms of life cycle simplification in field-derived and laboratory-selected African trypanosomes

## Abstract

African trypanosomes undergo development to transmissible stumpy forms in their mammalian host to favor uptake by their tsetse fly vector. However, Trypanosoma brucei evansi and Trypanosoma brucei equiperdum have simplified their lifecycle by escaping dependence on tsetse allowing an expanded geographical range, with direct transmission achieved via biting flies or through sexual transmission between animals. Concomitantly, stumpy formation is lost, and the isolates are described as monomorphic. Through genomic analysis of distinct field isolates we identified and functionally confirmed molecular changes that reduce stumpy formation. Further, by laboratory selection for reduced stumpy formation, we identified reversible steps in the initial development to monomorphism. This identifies a trajectory of events that simplify the trypanosome life cycle with impact on disease spread, geographical range and virulence.

## Scripts and data

Snakemake was used to analyse data
The following script was used to summarise data and produce the target databases. Python scripts are included in this repository which summarise output files from major tools, these are called from within the Snakefiles. Finally, the targets are identified with the final_db.py script. 

```
snakemake Snakefile.mapping
snakemake Snakefile.calling
snakemake Snakefile.delly
snakemake Snakefile.paml
python final_db.py
```

All figures in the manuscript were produced in R, using the four R scripts in the R_scripts directory.

You can download the aligned fasta sequences, one per isolate, for each CDS.

```
wget https://github.com/goldrieve/Mechanisms-of-life-cycle-simplification/raw/master/alignment_fasta/alignment.fasta.tar.gz
tar -xvf alignment.fasta.tar.gz
```
