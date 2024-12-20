# off_target_reads
Python script for counting off-target reads in targeted sequencing data

## Reference
Peter Micallef, Manuel Luna Santamaría, Daniel Andersson, Tobias Österlund, Stefan Filges, Gustav Johansson, Anders Ståhlberg. Digital sequencing using structured unique molecular identifiers. Manuscript

## Usage
`python3 off_target.py -i filename.sorted.bam -b bedfile.bed`

Optionally, the argument `-o output_directory` can be used to define a different output directory.

## Usage example
This scripts takes as input sorted BAM (*.sorted.bam*) alignment files:
`python3 off_target.py -i Sample_1.sorted.bam -b Sample_1.bed`

This command will index the bam file in case it was not previously indexed. Then, it will analyze the alignment file using the bed file as reference for the genomic position of the targeted amplicons. The script produces a *.tsv* file containing the number of reads detected for each target and the total number of off-target reads in the alignment, as well as the percentage of the total reads aligned to each target and percentage of off-target reads:

|Assay|Reads|Reads_percent|
|---|---|---|
|Assay 1|485558|12.47|
|Assay 2|452108|11.61|
|Assay 3|654047|16.81|
|Off-target|2301184|59.11|

## Requirements
Python 3 and the python package pysam are needed.
