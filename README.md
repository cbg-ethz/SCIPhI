# SCIPhI

## Description

**Single-cell mutation identification via phylogenetic inference** (SCIPhI) is a new approach to mutation detection in individual tumor cells by leveraging the evolutionary relationship among cells. SCIPhI, jointly calls mutations in each individual cell and estimates the tumor phylogeny on these cells. Employing a Markov Chain Monte Carlo scheme we robustly account for the various sources of noise in single-cell sequencing data. Our approach enables us to reliably call mutations in each single cell even in experiments with high dropout rates and missing data

## Availability

**SCIPhI** is freely available under a GNU General Public License v3.0 at https://github/cbg-ETHZ/SCIPhI

## How to install **SCIPhI**

SCIPhI has the **following dependencies** which need to be installed:

* `Boost >= 1.6.x`,
* `SeqAn >= 2.3.2`,
* `DLIB >= 19.9`
* `zlib`

Clone SCIPhI with submodules from the soucre code:

`git clone --recurse-submodules https://github.com/cbg-ethz/SCIPhI.git`

In order to install SCIPhI issue the following commands in the github directory:

`autoreconf -vif`

`./configure --with-boost=BOOST_PATH SEQAN_INCLUDEDIR=SEQAN_PATH DLIB_INCLUDEDIR=DLIB_PATH`

`make`

An alternative is to use CMake (tested with version >= 3.11.1). Simply issue the commands:

`mkdir build`

`cd build`

`cmake ..`

`cmake --build . --config Release`

Afterwards you will find SCIPhI in the build directory.

## Run SCIPhI

SCIPhI expects the sequencing information to be passed in form of the well known mpileup format (http://www.htslib.org/doc/samtools.html). In order to generate such a file you need to align your fastq files to a reference and post process the result (e.g., following the instuctions here: https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165). 

In order to see all available options type

`sciphi -h`

Executing

`sciphi -o result --in cellNames.txt --seed 42 example.mpileup`

will run SCIPhI using the cell names provided in *cellNames.txt* (same order as in the mpileup file). Note that *cellNames.txt* is a tab delimited file with the cell name in the first column and a cell type identifier in the second column. The cell type can be either *CT* (tumor cell), *CN* (control normal cell), or *BN* (control bulk normal). Note that SCIPhI assumes a pileup against a reference and ignores positions with 'N' as reference.
