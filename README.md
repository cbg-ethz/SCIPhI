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
* 'zlib'

In order to install SCIPhI issue the following commands:

autoreconf -vif
./configure --with-boost=BOOST_PATH SEQAN_INCLUDEDIR=SEQAN_PATH DLIB_INCLUDEDIR=DLIB_PATH
make

## Run SCIPhI

In order to see all available options

./sciphi -h
