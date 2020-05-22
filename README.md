# fqc
![github version](https://img.shields.io/badge/Version-0.0.1-informational)
[![pypi version](https://img.shields.io/pypi/v/fqc)](https://pypi.org/project/fqc/0.0.1/)
![python versions](https://img.shields.io/pypi/pyversions/fqc)
![status](https://github.com/pachterlab/fqc/workflows/CI/badge.svg)
[![pypi downloads](https://img.shields.io/pypi/dm/fqc)](https://pypi.org/project/fqc/)
[![license](https://img.shields.io/pypi/l/fqc)](LICENSE)

Detect single-cell technology that was used to generate a set of FASTQ files
or a single BAM file. If a BAM file is provided, it is split into FASTQs.

## Installation

```
pip install fqc
```

## Usage

### Detect the technology of a set of FASTQ files
```
fqc [FASTQ1] [FASTQ2] ...
```
where `[FASTQ1]` and `[FASTQ2]` are FASTQ files.

### Detect the technology of a single BAM file and split it into FASTQs
```
fqc [BAM]
```
where `[BAM]` is a BAM file.
