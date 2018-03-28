#!/bin/bash

module load eth_proxy;

snakemake -s SRP044380.snake --configfile config.json --use-conda -p -n
