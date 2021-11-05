# Copyright (c) 2021 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
import argparse
import functools
import sys
from typing import Callable, Generator, Iterable, List

import dnaio

import xopen  # type: ignore

from .optimized_algorithms import qualmean, qualmedian

DEFAULT_COMPRESSION_LEVEL = 2


def file_to_fastq_records(filepath: str) -> Generator[dnaio.Sequence,
                                                      None, None]:
    """Parse a FASTQ file into a generator of Sequence objects"""
    opener = functools.partial(xopen.xopen, threads=0)
    with dnaio.open(filepath, opener=opener) as record_h:  # type: ignore
        yield from record_h


def fastq_records_to_file(records: Iterable[dnaio.Sequence], filepath: str,
                          compression_level: int = DEFAULT_COMPRESSION_LEVEL):
    with xopen.xopen(filepath, mode='wb', threads=0,
                     compresslevel=compression_level) as output_h:
        for record in records:
            output_h.write(record.fastq_bytes())


def mean_quality_filter(quality: float, record: dnaio.Sequence) -> bool:
    """The mean quality of the FASTQ record is equal or above the given
    quality value."""
    return qualmean(record.qualities_as_bytes()) >= quality  # type: ignore


def median_quality_filter(quality: float, record: dnaio.Sequence) -> bool:
    """The median quality of the FASTQ record is equal or above the given
    quality value."""
    return qualmedian(record.qualities_as_bytes()) >= quality  # type: ignore


def min_length_filter(min_length: int, record: dnaio.Sequence) -> bool:
    """The length of the sequence in the FASTQ record is at least min_length"""
    return len(record.sequence) >= min_length


def max_length_filter(max_length: int, record: dnaio.Sequence) -> bool:
    """The length of the sequence in the FASTQ record is at most max_length"""
    return len(record.sequence) <= max_length


# Store filter names for use on the command line interface. Also store a
# tuple of types so the command line arguments (strings) can be converted
# in the appropiate types.
FILTERS = {"mean_quality": (mean_quality_filter, (float,), ("quality",)),
           "median_quality": (median_quality_filter, (float,), ("quality",)),
           "min_length": (min_length_filter, (int,), ("length",)),
           "max_length": (max_length_filter, (int,), ("length",))}


def print_filter_help():
    for filter_name, filter_tuple in FILTERS.items():
        filter_func, _, arg_names = filter_tuple
        # Reuse the docstring for the filter explanation.
        # Convert all whitespace in docstring to space.
        filter_usage = f"{filter_name}:<{'>,<'.join(arg_names)}>"
        filter_explanation = " ".join(filter_func.__doc__.split())
        print(f"{filter_usage:<30} {filter_explanation}")


def filter_string_to_filters(filter_string: str
                             ) -> List[Callable[[dnaio.Sequence], bool]]:
    """Convert a filter string such as 'min_length:50|mean_quality:20 into
    a list of filter functions that can be used by Python's builtin filter
    function."""
    filters: List[Callable[[dnaio.Sequence], bool]] = []
    for single_filter_string in filter_string.split('|'):
        filter_name, filter_argstring = single_filter_string.split(':')
        try:
            filter_function, filter_argtypes, _ = FILTERS[filter_name]
        except KeyError:
            raise ValueError(f"Unknown filter: {filter_name}. Choose one of:"
                             f" {' '.join(FILTERS.keys())}")
        # Convert the strings from the command line in the appropriate types
        filter_args = [filter_argtypes[pos](arg) for pos, arg
                       in enumerate(filter_argstring.split(','))]

        filters.append(functools.partial(filter_function, *filter_args))
    return filters


def filter_fastq(filter_string: str,
                 input_file: str,
                 output_file: str,
                 compression_level: int = DEFAULT_COMPRESSION_LEVEL):
    """
    Filter a FASTQ input file with the filters in filter_string and write
    the results to the output file.

    :param filter_string: A string representing one or multiple filters. For
    more information see the documentation.
    :param input_file: A FASTQ input filename. Compressed files are handled
    automatically.
    :param output_file: A FASTQ output filename. Compressed files are handled
    automatically.
    :param compression_level: Compression level for the output files (if
    applicable)
    """
    fastq_records = file_to_fastq_records(input_file)
    filtered_fastq_records: Iterable[dnaio.Sequence] = fastq_records
    for filter_func in filter_string_to_filters(filter_string):
        filtered_fastq_records = filter(filter_func, filtered_fastq_records)
    fastq_records_to_file(filtered_fastq_records, output_file,
                          compression_level=compression_level)


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "filters",
        help="Filters and arguments. For example: mean_quality:20, for "
             "filtering all reads with an average quality below 20. Multiple "
             "filters can be applied by separating with the | symbol. For "
             "example: min_length:100|mean_quality:20.  Make sure to use "
             "faster filters (length) before slower ones (quality) for "
             "optimal performance. Use --help-filters to print all the "
             "available filters.")
    parser.add_argument("input",
                        help="Input FASTQ file. Compression format "
                             "automatically detected. ")
    parser.add_argument("--help-filters", action="store_true",
                        help="Print all the available filters.")
    parser.add_argument("-o", "--output",
                        default="-",
                        help="Output FASTQ file. Compression format "
                             "automatically determined by file extension. "
                             "Default: stdout.")
    parser.add_argument("-l", "--compression-level", type=int,
                        default=DEFAULT_COMPRESSION_LEVEL,
                        help=f"Compression level for the output files. "
                             f"Relevant when output files have a .gz "
                             f"extension. Default: {DEFAULT_COMPRESSION_LEVEL}"
                        )
    return parser


def main():
    if "--help-filters" in sys.argv[1:]:
        print_filter_help()
        sys.exit(0)
    args = argument_parser().parse_args()
    filter_fastq(filter_string=args.filters,
                 input_file=args.input,
                 output_file=args.output,
                 compression_level=args.compression_level)


if __name__ == "__main__":  # pragma: no cover
    main()
