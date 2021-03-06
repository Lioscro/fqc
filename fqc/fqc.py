import logging
from collections import OrderedDict
from itertools import permutations

import numpy as np
import scipy.stats as stats

from .bam import BAM
from .fastq import Fastq
from .technologies import OrderedTechnology, TECHNOLOGIES
from .utils import open_as_text

logger = logging.getLogger(__name__)


def all_ordered_technologies(technologies=None, n=1):
    """Given a list of technologies, return all possible OrderedTechnology objects.

    :param technologies: list of Technology objects, defaults to `None`
    :type technologies: list, optional
    :param n: number of FASTQs, defaults to `1`
    :type n: int, optional

    :return: list of OrderedTechnology objects
    :rtype: list
    """
    technologies = technologies or TECHNOLOGIES
    ordered = []
    for technology in technologies:
        ordered.extend([
            OrderedTechnology(technology, permutation)
            for permutation in permutations(range(n))
        ])
    return ordered


def extract_barcodes_umis(reads, technologies=None):
    """Extract all sequences in barcode and UMI positions for each given
    technology for all possible orderings of the FASTQs.

    This function returns a 3-tuple of `barcodes`, `umis`, `invalids`.
    `barcodes`: Has technology names as outer keys, and FASTQ orderings as inner
                keys. For example, `barcodes['10xv2'][(1, 0)]` contains all
                barcode sequences extracted according to the 10xv2 technology,
                where read 0 comes from FASTQ 1 and read 1 comes from FASTQ 0.
                The list of barcodes is actually a list of list of barcodes,
                since there can be multiple barcode sections for a given technology.
    `umis`: Same structure as `barcodes`, but instead contains UMI sequences.
    `invalids`: A dictionary with technology names as keys and a set
                of FASTQ orderings as values. Any FASTQ ordering in the `invalids`
                dictionary, along with the specific technology, is invalid.

    :param reads: an ordered dictionary with the path to fastqs as keys and
                  a list of reads as values
    :type reads: OrderedDict
    :param technologies: list of OrderedTechnology objects to consider, defaults to `None`
    :type technologies: list, optional

    :return: 3-tuple (`barcodes`, `umis`, `invalids`) of dictionaries
    :rtype: list
    """
    logger.debug('Extracting barcode and UMI sequences')
    technologies = technologies or all_ordered_technologies(
        TECHNOLOGIES, len(reads)
    )

    # Dictionaries that contain list of barcodes/umis for each possible
    # permutation.
    barcodes = {}
    umis = {}
    invalids = {}
    for i, reads in enumerate(zip(*list(reads.values()))):
        for ordered in technologies:
            technology = ordered.technology
            permutation = ordered.permutation

            t_barcodes = barcodes.setdefault(technology.name, {})
            t_umis = umis.setdefault(technology.name, {})
            t_invalid = invalids.setdefault(technology.name, set())

            # If the permutation is [1, 0, 2], then read 0 is from fastq 1,
            # read 1 is from fastq 0, and read 2 is from fastq 2
            if permutation not in t_invalid:
                p_barcodes = t_barcodes.setdefault(permutation, [])
                p_umis = t_umis.setdefault(permutation, [])

                bc_list = []
                for substring in technology.barcode_positions:
                    read = reads[permutation[substring.file]]

                    if len(read) <= substring.start or len(read
                                                           ) < substring.stop:
                        logger.debug((
                            f'Technology {ordered} is '
                            'invalid due to barcode sequence length.'
                        ))
                        t_invalid.add(permutation)
                        break

                    bc_list.append(read[substring.start:substring.stop])
                p_barcodes.append(bc_list)

                umi_list = []
                for substring in technology.umi_positions:
                    read = reads[permutation[substring.file]]

                    if len(read) <= substring.start or len(read
                                                           ) < substring.stop:
                        logger.debug((
                            f'Technology {ordered} is '
                            'invalid due to UMI sequence length.'
                        ))
                        t_invalid.add(permutation)
                        break

                    umi_list.append(read[substring.start:substring.stop])
                p_umis.append(umi_list)

                if permutation in t_invalid:
                    del barcodes[technology.name][permutation]
                    del umis[technology.name][permutation]

    return barcodes, umis, invalids


def filter_files(reads, technologies=None):
    """Filter for possible technologies using the number of files (`n_files`).

    :param reads: an ordered dictionary with the path to fastqs as keys and
                  a list of reads as values
    :type reads: OrderedDict
    :param technologies: list of possible OrderedTechnology objects, defaults to `None`
    :type technologies: list, optional

    :return: list of OrderedTechnology objects. All technologies in the returned
             list are possible technologies the FASTQs were derived from.
    :rtype: list
    """
    technologies = technologies or all_ordered_technologies(
        TECHNOLOGIES, len(reads)
    )

    # Filter
    possible = []
    for ordered in technologies:
        technology = ordered.technology
        if len(reads) == technology.n_files:
            possible.append(ordered)

    return possible


def filter_barcodes_umis(reads, technologies=None):
    """Filter for possible technologies using barcodes and UMI positions.

    :param reads: an ordered dictionary with the path to fastqs as keys and
                  a list of reads as values
    :type reads: OrderedDict
    :param technologies: list of possible OrderedTechnology objects, defaults to `None`
    :type technologies: list, optional

    :return: list of OrderedTechnology objects. All technologies in the returned
             list are possible technologies the FASTQs were derived from.
    :rtype: list
    """
    technologies = technologies or all_ordered_technologies(
        TECHNOLOGIES, len(reads)
    )
    barcodes, umis, invalids = extract_barcodes_umis(reads, technologies)
    possible = []

    # Filter with barcodes.
    # For all technologies with available whitelist, count the number of
    # sequences that match the barcodes.
    barcode_technologies = [
        technology for technology in TECHNOLOGIES
        if technology.whitelist_path and technology.name in
        [ordered.technology.name for ordered in technologies]
    ]

    logger.debug(
        f'Checking technologies with whitelists: {", ".join(str(technology) for technology in barcode_technologies)}'
    )
    max_ordered = None
    max_count = 0
    for technology in barcode_technologies:
        with open_as_text(technology.whitelist_path, 'r') as f:
            whitelist = set(f.read().splitlines())

        for p in barcodes[technology.name]:
            ordered = OrderedTechnology(technology, p)
            n = len(barcodes[technology.name][p])
            count = 0
            for bc_list in barcodes[technology.name][p]:
                if ''.join(bc_list) in whitelist:
                    count += 1

            logger.debug(
                f'Technology {ordered} has {count}/{n} matching barcodes.'
            )
            if count / n > 0.5 and count > max_count:
                max_ordered = ordered
                max_count = count
    if max_ordered is not None:
        possible.append(max_ordered)
        logger.debug(
            f'Technology {max_ordered} passed whitelist filter with {max_count}/{n} matching barcodes'
        )
    return possible

    # # Check technologies without whitelist
    # other_technologies = [
    #     technology for technology in TECHNOLOGIES
    #     if not technology.whitelist_path and technology.name in
    #     [ordered.technology.name for ordered in technologies]
    # ]
    # logger.debug(
    #     f'Checking technologies without whitelists: {", ".join(str(ordered) for ordered in other_technologies)}'
    # )
    #
    # def passes_stats_threshold(mean, std, kurtosis, skew):
    #     """Helper function to check if the barcode statistics passes a
    #     certain threshold.
    #
    #     TODO: is there a way to make these numbers not arbitrary/manually selected?
    #     """
    #     if mean > 2 and std > 10 and kurtosis < 2000 and skew < 50:
    #         return True
    #     return False
    #
    # max_ordered = None
    # max_mean = 0
    # for technology in other_technologies:
    #     for p in barcodes[technology.name]:
    #         bcs = [''.join(bc_list) for bc_list in barcodes[technology.name][p]]
    #         counts = list(Counter(bcs).values())
    #         mean = np.mean(counts)
    #         std = np.std(counts)
    #         kurtosis = stats.kurtosis(counts)
    #         skew = stats.skew(counts)
    #         logger.debug((
    #             f'Barcodes for technology {technology} has {len(set(bcs))}/{n} unique barcodes, '
    #             f'mean={mean:.2f}, std={std:.2f}, kurtosis={kurtosis:.2f}, skew={skew:.2f}'
    #         ))
    #         if mean > max_mean and passes_stats_threshold(mean, std, kurtosis,
    #                                                       skew):
    #             max_ordered = OrderedTechnology(technology, p)
    #             max_mean = mean
    # if max_ordered is not None:
    #     possible.append(max_ordered)
    #     logger.debug(f'Technology {max_ordered} passed barcode filter')
    #     return possible
    #
    # return possible


def is_single_cell(reads):
    """Given a list of list of reads as values, determine if they are from a
    single-cell experiment.

    Only supports paired-end reads for now.
    1) If the read lengths differ, is single-cell.
    2) If one of the reads has more A's at the end of the read, is single-cell.
    3) Take a highly-exposed housekeeping gene and align the reads to that.
       For bulk paired-end reads, the pairs will map to both ends of the gene.

    :param reads: list of lists, with the inner list containing reads
    :type reads: list

    :return: whether or not the reads are single-cell
    :rtype: bool
    """
    if len(reads) != 2:
        return False

    # Check read lengths.
    lengths = [[len(r) for r in rs] for rs in reads]
    _, p_value = stats.ttest_ind(*lengths, equal_var=False)
    logger.debug(f'Lengths p-value={p_value}')
    if not np.isnan(p_value) and p_value < 0.05:
        return True

    # Check for poly-A, allowing N's to be considered as A's.
    # TODO: should we consider 'islands' of A's as well?
    polys = [[l - len(r.upper().rstrip('AN'))
              for r, l in zip(rs, ls)]
             for rs, ls in zip(reads, lengths)]
    means = [np.mean(poly) for poly in polys]
    logger.debug(f'Poly-A means={means}')
    # TODO: is there a better way to decide whether to run the T-test on
    # the number of trailing A's?
    if any(mean > 1 for mean in means):
        _, p_value = stats.ttest_ind(*polys, equal_var=False)
        logger.debug(f'Poly-A p-value={p_value}')
        if p_value < 0.05:
            return True

    # TODO: implement aligning to housekeeping gene
    return False


def fqc_bam(path, split=False, prefix='', threads=4):
    bam = BAM(path)
    if split:
        return bam.to_fastq(prefix=prefix, threads=threads)
    return bam.technology


def fqc_fastq(fastqs, skip, n, technologies=None):
    """Detect single-cell technology and file ordering.

    :param fastqs: paths to FASTQs
    :type fastqs: list
    :param skip: number of reads to skip at the beginning
    :type skip: int
    :param n: number of reads to consider
    :type n: int
    :param technologies: list of possible OrderedTechnology objects, defaults to `None`
    :type technologies: list, optional

    :return: tuple of a list of paths to FASTQs and a list of TechnologyOrdering objects
    :rtype: tuple
    """
    # Read n reads, starting from read skip
    reads = OrderedDict()
    for path in fastqs:
        fastq = Fastq(path)
        rs = fastq[skip:skip + n]
        logger.info(
            f'Read {len(rs)} reads after skipping the first {skip} reads from {path}'
        )

        # Check if index fastq, which will have very low variation.
        if len(set(rs)) / len(rs) < 0.05:
            logger.warning((
                f'FASTQ {path} has {len(set(rs))}/{len(rs)} unique sequences. '
                'This file will be considered an index read and will be ignored.'
            ))
            continue
        reads[path] = rs
    logger.info('Only the following FASTQs will be considered:')
    for path in reads.keys():
        logger.info(f'\t{path}')

    # if not is_single_cell(list(reads.values())):
    #     raise Exception(
    #         'The provided FASTQs are not from a single-cell experiment.'
    #     )

    logger.info(f'Filtering based on number of files: {len(reads)}')
    technologies = filter_files(reads)
    logger.debug(
        f'{len(technologies)} passed the filter: {", ".join(str(technology) for technology in technologies)}'
    )

    logger.info('Filtering based on barcode and UMI sequences')
    technologies = filter_barcodes_umis(reads, technologies)
    logger.debug(
        f'{len(technologies)} passed the filter: {", ".join(str(technology) for technology in technologies)}'
    )

    return list(reads.keys()), technologies
