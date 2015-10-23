#!/bin/bash
model=$1
echo "python scripts/limitPlotter.py -M Asymptotic  -v -p ../${model}"
python scripts/limitPlotter.py -M Asymptotic -v -p ${model}  -e -r
