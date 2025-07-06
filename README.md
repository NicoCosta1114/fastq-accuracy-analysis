# FASTQ Accuracy Analysis

This project includes a Python script to evaluate the **base call accuracy** of a FASTQ file using **Phred quality scores**, commonly produced by sequencing technologies such as Illumina.

## What does it do?

It reads a `.fastq` file, calculates the **mean base call accuracy** at each position across all reads, and generates a **plot with confidence intervals** (95% and 99%) to assess sequencing quality.

> Accuracy is estimated from Phred scores using the formula:  
> **accuracy = 1 - 10^(-Q/10)**


## Files included

- `fastq_accuracy.py` – Main script for accuracy analysis.
- `unknown_illumina_2024.fastq` – Example FASTQ file (Illumina reads).


## How to run

Make sure you have **Python 3** and the required libraries installed:

bash
pip install biopython matplotlib numpy

Then run:

python fastq_accuracy.py unknown_illumina_2024.fastq

A graph will be displayed showing the mean accuracy per base position with confidence intervals.

## Technologies used

- Python 3
- Biopython
- NumPy
- Matplotlib

