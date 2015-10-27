#!/bin/bash
model=$1
echo "python limitPlotter.py -M Asymptotic  -v -p ../${model}"
python limitPlotter.py -M Asymptotic -v -p ${model}  -e -r
