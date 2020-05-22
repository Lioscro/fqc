import logging
from itertools import permutations

from .config import BARCODE_THRESHOLD
from .technologies import OrderedTechnology, TECHNOLOGIES
from .utils import fastq_reads, open_as_text

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


def extract_barcodes_umis(fastqs, skip, n, technologies=None):
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

    :param fastqs: list of paths to FASTQ files
    :type fastqs: list
    :param skip: number of reads to skip at the beginning
    :type skip: int
    :param n: number of reads to consider
    :type n: int
    :param technologies: list of OrderedTechnology objects to consider, defaults to `None`
    :type technologies: list, optional

    :return: 3-tuple (`barcodes`, `umis`, `invalids`) of dictionaries
    :rtype: list
    """
    logger.debug('Extracting barcode and UMI sequences')
    technologies = technologies or all_ordered_technologies(
        TECHNOLOGIES, len(fastqs)
    )

    # Dictionaries that contain list of barcodes/umis for each possible
    # permutation.
    barcodes = {}
    umis = {}
    invalids = {}
    for i, reads in enumerate(zip(*[fastq_reads(fastq) for fastq in fastqs])):
        # Ignore first SKIP_READS, and consider next N_READS only
        if i < skip:
            continue
        if i >= skip + n:
            break

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


def filter_files(fastqs, technologies=None):
    """Filter for possible technologies using the number of files (`n_files`).

    :param fastqs: list of paths to FASTQ files
    :type fastqs: list
    :param technologies: list of possible OrderedTechnology objects, defaults to `None`
    :type technologies: list, optional

    :return: list of OrderedTechnology objects. All technologies in the returned
             list are possible technologies the FASTQs were derived from.
    :rtype: list
    """
    technologies = technologies or all_ordered_technologies(
        TECHNOLOGIES, len(fastqs)
    )

    # Filter
    possible = []
    for ordered in technologies:
        technology = ordered.technology
        if len(fastqs) == technology.n_files:
            possible.append(ordered)

    return possible


def filter_barcodes_umis(fastqs, skip, n, technologies=None):
    """Filter for possible technologies using barcodes and UMI positions.

    :param fastqs: list of paths to FASTQ files
    :type fastqs: list
    :param skip: number of reads to skip at the beginning
    :type skip: int
    :param n: number of reads to consider
    :type n: int
    :param technologies: list of possible OrderedTechnology objects, defaults to `None`
    :type technologies: list, optional

    :return: list of OrderedTechnology objects. All technologies in the returned
             list are possible technologies the FASTQs were derived from.
    :rtype: list
    """
    technologies = technologies or all_ordered_technologies(
        TECHNOLOGIES, len(fastqs)
    )
    barcodes, umis, invalids = extract_barcodes_umis(
        fastqs, skip, n, technologies
    )
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
            count = 0
            for bc_list in barcodes[technology.name][p]:
                if ''.join(bc_list) in whitelist:
                    count += 1

            if count > n * BARCODE_THRESHOLD and count > max_count:
                max_ordered = OrderedTechnology(technology, p)
                max_count = count
    if max_ordered is not None:
        possible.append(max_ordered)
        logger.debug(
            f'Technology {max_ordered} passed whitelist filter with {max_count}/{n} matching barcodes'
        )
        return possible

    # Check technologies without whitelist
    other_technologies = [
        technology for technology in TECHNOLOGIES
        if not technology.whitelist_path and technology.name in
        [ordered.technology.name for ordered in technologies]
    ]
    logger.debug(
        f'Checking technologies without whitelists: {", ".join(str(ordered) for ordered in other_technologies)}'
    )
    max_ordered = None
    max_diff = 0
    max_unique = 0
    for technology in other_technologies:
        # Expected number of unique barcodes if nucleotide probability is uniform
        n_barcodes = 4**sum(
            pos.stop - pos.start for pos in technology.barcode_positions
        )
        expected = n_barcodes * (1 - ((n_barcodes - 1) / n_barcodes)**n)

        for p in barcodes[technology.name]:
            unique = set([
                ''.join(bc_list) for bc_list in barcodes[technology.name][p]
            ])
            diff = expected - len(unique)

            if diff > max_diff:
                max_ordered = OrderedTechnology(technology, p)
                max_diff = diff
                max_unique = len(unique)
    if max_ordered is not None:
        possible.append(max_ordered)
        logger.debug(
            f'Technology {max_ordered} passed barcode filter with {max_unique}/{n} unique barcodes'
        )
        return possible

    return possible


def fqc(fastqs, skip, n, technologies=None):
    """Detect single-cell technology and file ordering.

    :param fastqs: paths to FASTQs
    :type fastqs: list
    :param skip: number of reads to skip at the beginning
    :type skip: int
    :param n: number of reads to consider
    :type n: int
    :param technologies: list of possible OrderedTechnology objects, defaults to `None`
    :type technologies: list, optional

    :return: TechnologyOrdering object
    :rtype: TechnologyOrdering
    """
    technologies = technologies or all_ordered_technologies(
        TECHNOLOGIES, len(fastqs)
    )

    logger.info('Filtering based on number of files')
    technologies = filter_files(fastqs, technologies)
    logger.info(
        f'{len(technologies)} passed the filter: {", ".join(str(technology) for technology in technologies)}'
    )

    logger.info('Filtering based on barcode and UMI sequences')
    technologies = filter_barcodes_umis(fastqs, skip, n, technologies)
    logger.info(
        f'{len(technologies)} passed the filter: {", ".join(str(technology) for technology in technologies)}'
    )

    return technologies
