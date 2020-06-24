# Demultiplexing FASTQ Reads by Inline Barcode
Demultiplexing using inline barcode is the process of extracting reads into separated file based on the barcode at the beginning of the read sequence.

This package contains the command line program for demultiplexing FASTQ reads by inline barcodes. The program is designed to utilize the power of multi-core or multi-CPU computer, like the VMs on GCP or AWS with 16 or 32 CPUs. See [Project Readme](/Readme.md) for more information about multi-processing.

## Command
```
python -m Cancer.run demux_inline [-h] --r1 R1 [R1 ...] --r2 R2 [R2 ...]
                           [--barcode BARCODE [BARCODE ...]]
                           [--unmatched UNMATCHED] [--error_rate ERROR_RATE]
                           [--score SCORE] [--penalty PENALTY] [--stats STATS]
                           [--name NAME] [--header HEADER]

optional arguments:
  -h, --help            show this help message and exit
  --r1 R1 [R1 ...]      FASTQ R1 files
  --r2 R2 [R2 ...]      FASTQ R2 files
  --barcode BARCODE [BARCODE ...]
                        Barcode and the output file prefix in the format of
                        BARCODE=PREFIX
  --unmatched UNMATCHED
                        File path for saving the unmatched reads.
  --error_rate ERROR_RATE
                        Max Error Allowed, defaults to 20%
  --score SCORE         Score for each bp matched, defaults to 1
  --penalty PENALTY     Penalty for each bp mis-matched or gap, defaults to 10
  --stats STATS         Specify a CSV file path to save the statistics.
  --name NAME           For statistics, Sample name
  --header HEADER       For statistics, the column header for matched read
                        counts
```

## Example

For example, assume we have following reads with 10bp inline barcode:
```
@MN01031:94:000H32JLF:1:11102:17417:1081 1:N:0:GTTGTTCG+CTACAAGG
AGTCGACATGCCAGAGATGAAATCATAGGGAGTTGAAGCTGACCTCTTGCACTGAGTCAGTTCCTCAGTCAGGT
+
FFFFFFFFFFFFFFFFFFFFFFFFFFAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF6FFF
@MN01031:94:000H32JLF:1:11102:18381:1081 1:N:0:GTTGTTCG+CTACAAGG
GACAAACACAAATCTGTCCCTTTTCACAATGTTTCCAGGAAAATACTCTCGTTGGCTGTAATAACAACTTTGAC
+
FAFFFFFFFAFFFFFFFAAAFFA/F=AFF/FFFFF/FFFFFFFFFF/FAFAFA/FFAF=FFFFFFFAFA/AAFA
```

### Extract reads matching different barcodes into different files
The barcode for the first read will be `AGTCGACATG` and the barcode for the second read will be `GACAAACACA`. The demultiplxing process can extract the two reads into two different files. The barcode sequences (along with the corresponding quality scores) are trimmed before the reads are saved into new files. The following command extracts: 
* Reads with barcode `AGTCGACATG` into FILE_PREFIX_1.R1.fastq.gz and FILE_PREFIX_1.R2.fastq.gz.
* Reads with barcode `GACAAACACA` into FILE_PREFIX_2.R1.fastq.gz and FILE_PREFIX_2.R2.fastq.gz.
```
python -m Cancer.run demux_inline --r1 R1_FILE_PATH.fastq.gz --r2 R2_FILE_PATH.fastq.gz --barcode AGTCGACATG=FILE_PREFIX_1 GACAAACACA=FILE_PREFIX_2
```

All commands in the `Cancer` package can be called using `python -m Cancer.run SUB_COMMAND`. In this case, `SUB_COMMAND` is `demux_inline`.

The barcodes are specified in the format of `BARCODE=PREFIX`. `BARCODE` is the barcode we are look for at the beginning of the read sequence. `PREFIX` is the output file path prefix for saving the corresponding extracted reads. For paired-end data, the `R1.fastq.gz` and `R2.fastq.gz` are appended the prefix to form the full file paths. This program will extract the read pairs if the barcode is matching the beginning of either one of the paired-end reads.

This program requires the user to specify the inline barcodes. Not every read will matched with a barcode. In the above example, if only one barcode, e.g. `AGTCGACATG`, is specified, only the read matching the barcode will be trimmed and extracted.
```
python -m Cancer.run demux_inline --r1 R1_FILE_PATH.fastq.gz --r2 R2_FILE_PATH.fastq.gz --barcode AGTCGACATG=FILE_PREFIX_1
```

### Extract reads matching multiple barcodes into the same file
We can extract reads matching either barcode into the same set of files by specifying the same `PREFIX` for both barcodes:
```
python -m Cancer.run demux_inline --r1 R1_FILE_PATH.fastq.gz --r2 R2_FILE_PATH.fastq.gz --barcode AGTCGACATG=FILE_PREFIX_1 GACAAACACA=FILE_PREFIX_1
```
Alternatively, we can also specify multiple barcodes as space-separated strings in quotes:
```
python -m Cancer.run demux_inline --r1 R1_FILE_PATH.fastq.gz --r2 R2_FILE_PATH.fastq.gz --barcode "AGTCGACATG GACAAACACA"=FILE_PREFIX_1
```

### Extract reads not matching any barcode
The reads not matching any barcode can also be extracted by specifying another prefix using the `--unmatched` parameter. These reads will not be trimmed. The following command extracts the unmatched reads into UNMATCHED_FILE_PREFIX.R1.fastq.gz and UNMATCHED_FILE_PREFIX.R2.fastq.gz
```
python -m Cancer.run demux_inline --r1 R1_FILE_PATH.fastq.gz --r2 R2_FILE_PATH.fastq.gz --barcode "AGTCGACATG GACAAACACA"=FILE_PREFIX_1 --unmatched UNMATCHED_FILE_PREFIX
```

## Barcode Matching Algorithm
When barcode and the read sequence may not match exactly, they are matched using semi-global alignment (sg_de) algorithm implemented in the [parasail](https://github.com/jeffdaily/parasail-python) library. Gaps at the end of the read sequence are not penalized. 

### Mismatch and Error Rate
Mismatch is defined as the edit distance between the barcode and the beginning of the read sequence. One mismatch can be a substitution, an insertion or a deletion of a bp. Mismatches at the end of the read sequence are ignored. The maximum number of mismatches allowed is calculated as the error rate (specified by `--error_rate`) times the length of the barcode (and rounding down to the nearest integer). For example:
* A 10bp barcode with 0.2 error rate will allow a maximum of 2 mismatches
* An 8bp barcode with 0.2 error rate will allow a maximum of 1 mismatch.

### Score and Penalty
While our goal is to identify the alignment with minimum mismatch, the best alignment may not be unique. To rank the possible alignments more precisely, the semi-global algorithm calculates alignment scores using additional parameters:
* `--score`: The score for one bp match. Defaults to 1.
* `--penalty`: The penalty for one bp substitution, insertion or deletion (gap). Defaults to 10.

The alignment with the highest score will be considered as the best alignment. Semi-global alignment algorithm allows more complicated scoring by using different penalties for substitution, opening a gap or extending a gap. This program uses the same penalty score to simplify the parameters. The default score and penalty setting guarantees that the best alignment will have the minimum number of mismatches. However, this may not be the case for some other score and penalty values.

In general, a barcode should only contain "A", "C", "G" and "T". The read sequence may also contain "N". This program considers having an "N" in the read sequence as a mismatch. For example, when score is 1 and penalty is 10. The following penalty matrix is used:
```
     A    C    G    T    N
A    1  -10  -10  -10  -10
C  -10    1  -10  -10  -10
G  -10  -10    1  -10  -10
T  -10  -10  -10    1  -10
N  -10  -10  -10  -10    1
```

See also: [Parasail: SIMD C library for global, semi-global, and local pairwise sequence alignments](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0930-z)

## Output Statistics
User can also save the statistics of the demultiplexing process into a CSV file by specifying a file path using the `--stats` parameter. Each row in the CSV file contain the statistics for one barcode, including the following columns (reads means read pairs in the following description):
* `sample`, The sample name specified by the `--name` parameter.
* `barcode`, The barcode.
* `read1_percent`, The percentage of reads with R1 matching the barcode.
* `read2_percent`, The percentage of reads with R2 matching the barcode.
* `total_percent_HEADER`, The percentage of reads matching the barcode.
* `total`, Total number of reads from FASTQ files. This number should be the same for all rows.
* `HEADER_reads`: Total number of reads matching this barcode.
* `nonHEADER_reads`: Total number of unmatched reads. This number should be the same for all rows.
In the above column names, `HEADER` is the specified by the `--header` parameter.
