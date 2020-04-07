from .config import N_READS, SKIP_READS
from .technologies import TECHNOLOGIES
from .utils import fastq_reads


def filter_files(fastqs, technologies=None):
    """Filter for possible technologies using the number of files (`n_files`).

    :param fastqs: list of paths to FASTQ files
    :type fastqs: list
    :param technologies: list of possible technologies, defaults to `None`
    :type technologies: list, optional

    :return: list of filtered technologies. All technologies in the returned
             list are possible technologies the FASTQs were derived from.
    :rtype: list
    """
    technologies = technologies or TECHNOLOGIES

    # Filter
    possible = []
    for technology in technologies:
        if len(fastqs) == technology.n_files:
            possible.append(technology)

    return possible


def filter_barcodes_umis(fastqs, technologies=None):
    """Filter for possible technologies using barcodes and UMI positions.

    :param fastqs: list of paths to FASTQ files
    :type fastqs: list
    :param technologies: list of possible technologies, defaults to `None`
    :type technologies: list, optional

    :return: list of filtered technologies. All technologies in the returned
             list are possible technologies the FASTQs were derived from.
    :rtype: list
    """
    technologies = technologies or TECHNOLOGIES

    for fastq in fastqs:
        for i, read in enumerate(fastq_reads(fastq)):
            # Ignore first SKIP_READS, and consider next N_READS only
            if i < SKIP_READS:
                continue
            if i >= SKIP_READS + N_READS:
                break

            # Extract all barcode and UMI positions.


def fqc(fastqs):
    pass
