# Usage Count of Inline Barcodes in FASTQ Files
 Here, Inline Barcode is defined as a fixed length string starting from a fixed position (fixed substring window) in a read sequence.
 
 This command obtains the usage count of inline barcodes by checking the fixed substring window in each read sequence of the FASTQ files. Mismatches and shifts are not considered.

 The results are stored in a dictionary, where the keys are the barcodes and the values are the corresponding read counts. 
 
 For paired-end sequencing files, this program counts R1 and R2 separately. In the results, `_1` and `_2` are appended to the barcode to indicate if the counts are from R1 or R2. Note that for paired-end files, the total number of reads is twice the number of read pairs in the FASTQ files, the percentages in the outputs are with respect to the total number of reads processed.

## Command:
```
python -m Cancer.run count_inline_barcode [-h] --r1 R1 [R1 ...] [--r2 R2 [R2 ...]]
                                   [-s START] -l LENGTH

optional arguments:
  -h, --help            show this help message and exit
  --r1 R1 [R1 ...]      FASTQ R1 files
  --r2 R2 [R2 ...]      FASTQ R2 files
  -s START, --start START
                        Starting position of the barcode (0-based), defaults to 0
  -l LENGTH, --length LENGTH
                        Length of the barcode
```

For example, if the first 10bp of a read are considered as the barcode, the following command will count the barcode usage for a pair of FASTQ files.
```
python -m Cancer.run count_inline_barcode -l 10 --r1 R1_FILENAME --r2 R2_FILENAME
```

## API
```
from Cancer.fastq.barcode import BarcodeCounter

# Parameters:
# start: Starting position of the barcode (0-based)
# length: Length of the barcode
# fastq_files: A list of 2-tuples, each is a pair of file paths for paired-end fastq files
barcode_counter = BarcodeCounter(start, length).start(fastq_files)

# The counts attribute contains the usage counts
print(barcode_counter.counts)
```
