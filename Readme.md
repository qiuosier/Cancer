# The Cancer Toolbox
Python programs for processing FASTQ files on multiprocessor computer.
* [Accessing FASTQ files on Illumina BaseSpace](docs/BaseSpace.md)
* Demultiplexing undetermined reads from Illumina paired-end sequencing.
* [Demultiplexing FASTQ reads by inline barcode.](docs/DemultiplexInline.md)
* Compare reads in two pairs of paired-end FASTQ files.
* [Usage Count of Inline Barcodes in FASTQ Files](docs/CountInlineBarcode.md)

## Multiprocessing
The programs are designed to be executed on a server with multiple CPUs, such as the VM instances on GCP/AWS/Azure, or a Kubernetes node with multiple CPUs. A 16 CPU machine is recommended for processing a single pair of compressed FASTQ (fastq.gz) files. Higher number of CPUs may not be fully utilized due to the bottleneck of decompressing a single file. Programs like demultiplexing are capable of processing multiple pairs of FASTQ files for the same sample at the same time. The bottleneck generally shifts to the speed of the hard disk as the reads are coming from more pairs of FASTQ files. [Read more...](docs/MultiProcessing.md)

The [dnaio](https://github.com/marcelm/dnaio/) packages is used to read FASTQ files. Behind the scene it uses [xopen](https://github.com/marcelm/xopen/) to open compressed files. The `xopen` module uses [pigz](https://zlib.net/pigz/) to exploit multi-threading for compressing and decompressing data for a single file.

## FASTQ Processing Framework
A `FASTQProcessor` class (in `fastq.processor`) is designed to provide a framework for processing FASTQ files. It implements the logic for reading the file and processing the reads. This framework is capable of processing multiple pairs of FASTQ files for the same sample (for example, Illumina MiniSeq produces 4 pairs of FASTQ files for each sample). A reader process will be used to read and decompress each pair of files. The reads are put into a queue for processing. A number of workers are also started at the same time to process the reads from the queue. The workers are also responsible to compress and write reads into new FASTQ files if needed.

## Requirements
This toolbox requires Python 3.6 and all packages listed in the `requirements.txt`. The [pigz](https://zlib.net/pigz/) tool is required in order to use multiple threads for compression/decompression.

## Docker Image
This toolbox is available as a docker image: [qiuosier/caner](https://hub.docker.com/repository/docker/qiuosier/cancer).

## Usage
The command line programs uses `python -m Cancer.run` as entry point. Running `python -m Cancer.run --help` will display all available programs. Running `python -m Cancer.run SUB_PROGRAM --help` will display all options for the `SUB_PROGRAM`.
