# The Cancer Toolbox
Python programs for processing FASTQ files on multiprocessor computer.
* [Accessing FASTQ files on Illumina BaseSpace](docs/BaseSpace.md)
* Demultiplexing undetermined reads from Illumina paired-end sequencing.
* [Demultiplexing FASTQ reads by inline barcode.](docs/DemultiplexInline.md)
* Compare reads in two pairs of paired-end FASTQ files.
* [Usage Count of Inline Barcodes in FASTQ Files](docs/CountInlineBarcode.md)

## Multiprocessing
The programs are designed to be executed on a server with multiple CPUs, such as the VM instances on GCP/AWS/Azure, or a Kubernetes node with multiple CPUs. A 32 CPU machine is recommended for processing a single pair of compressed FASTQ (fastq.gz) files. Higher number of CPUs may not be fully utilized due to the bottleneck of decompressing a single file. Programs like demultiplexing are capable of processing multiple pairs of FASTQ files for the same sample at the same time. The bottleneck generally shifts to the speed of the hard disk as the reads are coming from more pairs of FASTQ files.

## Reading FASTQ files
The [dnaio](https://github.com/marcelm/dnaio/) packages is used to read FASTQ files. Behind the scene it uses [xopen](https://github.com/marcelm/xopen/) to open compressed files. The `xopen` module uses [pigz](https://zlib.net/pigz/) to exploit multiple CPUs for compressing and decompressing data for a single file.

## Requirements
This toolbox requires Python 3.6 and all packages listed in the `requirements.txt`.

## Docker Image
This toolbox is available as a docker image: [qiuosier/caner](https://hub.docker.com/repository/docker/qiuosier/cancer).

## Usage
The command line programs uses `python -m Cancer.run` as entry point. Running `python -m Cancer.run --help` will display all available programs.
