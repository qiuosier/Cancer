# Demultiplexing FASTQ Reads by Inline Barcode
Demultiplexing using inline barcode is the process of extracting reads into separated file based on the barcode at the beginning of the read sequence.

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
The barcode for the first read will be `AGTCGACATG` and the barcode for the second read will be `GACAAACACA`. The demultiplxing process can extract the two reads into two different files. The barcode sequences (along with the corresponding quality scores) are trimmed before the read is saved into the new files. The following command extracts: 
* Reads with barcode `AGTCGACATG` into FILE_PREFIX_1.R1.fastq.gz and FILE_PREFIX_1.R2.fastq.gz.
* Reads with barcode `GACAAACACA` into FILE_PREFIX_2.R1.fastq.gz and FILE_PREFIX_2.R2.fastq.gz
```
python -m Cancer.run demux_inline --r1 R1_FILE_PATH.fastq.gz --r2 R2_FILE_PATH.fastq.gz --barcode AGTCGACATG=FILE_PREFIX_1 GACAAACACA=FILE_PREFIX_2
```
The barcodes are specified in the format of `BARCODE=PREFIX`.


This program requires the user to specify the inline barcodes. Not every read will matched with a barcode. In the above example, if only one barcode, e.g. `AGTCGACATG`, is specified, only the read matching the barcode will be trimmed and extracted.

The reads not matching any barcode can also be extracted using the `--unmatched` parameter.

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
  --penalty PENALTY     Penalty for each bp mis-matched, defaults to 10
  --stats STATS         Specify a CSV file path to save the statistics.
  --name NAME           For statistics, Sample name
  --header HEADER       For statistics, the column header for matched read
                        counts
```