import argparse
import glob
import logging
import os

import fqc.fqc as fqc

logging.getLogger('fqc').setLevel(logging.ERROR)


def read_runs(path):
    with open(path, 'r') as f:
        lines = f.readlines()
    return [
        line.strip().split(',')
        for line in lines
        if not line.isspace() and not line.startswith('#')
    ]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('path')
    parser.add_argument('-s', default=1000, type=int)
    parser.add_argument('-n', default=100000, type=int)
    args = parser.parse_args()

    technology = os.path.splitext(os.path.basename(args.path))[0]
    runs = read_runs(args.path)

    try:
        for i, run in enumerate(runs):
            print(f'{technology} {i+1}/{len(runs)}')
            print(' '.join(run))
            if run[0].endswith('.bam'):
                assert technology == fqc.fqc_bam(run[0]).name
            else:
                fastqs, technologies = fqc.fqc_fastq(run, args.s, args.n)
                assert fastqs == run
                assert len(technologies) == 1
                assert technologies[0].technology.name == technology
    finally:
        # Remove any bam indexes that were created
        for index in glob.glob('*.bai'):
            os.remove(index)
