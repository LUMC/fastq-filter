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
import warnings
from typing import Callable, Generator, Iterable, List

import dnaio

import xopen  # type: ignore

from ._filters import (
    AverageErrorRateFilter,
    DEFAULT_PHRED_SCORE_OFFSET,
    MaximumLengthFilter,
    MedianQualityFilter,
    MinimumLengthFilter,
    qualmean,
    qualmedian
)

__all__ = [
    "file_to_fastq_records",
    "fastq_records_to_file",
    "filter_fastq",
    "AverageErrorRateFilter",
    "MaximumLengthFilter",
    "MedianQualityFilter",
    "MinimumLengthFilter",
    "qualmean",
    "qualmedian",
    "DEFAULT_PHRED_SCORE_OFFSET"
]

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


def filter_fastq(input_file: str,
                 output_file: str,
                 filters: List[Callable[[dnaio.SequenceRecord], bool]],
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
    for filter_func in filters:
        filtered_fastq_records = filter(filter_func, filtered_fastq_records)
    fastq_records_to_file(filtered_fastq_records, output_file,
                          compression_level=compression_level)


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.description = "Filter FASTQ files on various metrics."
    parser.add_argument("input",
                        help="Input FASTQ file. Compression format "
                             "automatically detected. Use - for stdin.")
    parser.add_argument("-o", "--output",
                        default="-",
                        help="Output FASTQ file. Compression format "
                             "automatically determined by file extension. "
                             "Default: stdout.")
    parser.add_argument("-l", "--min-length", type=int,
                        help="The minimum length for a read.")
    parser.add_argument("-L", "--max-length", type=int,
                        help="The maximum length for a read.")
    parser.add_argument("-e", "--average-error-rate", type=float,
                        help="The minimum average per base error rate.")
    parser.add_argument("-q", "--mean-quality", type=int,
                        help="Average quality. Same as the "
                             "'--average-error-rate' option but specified "
                             "with a phred score. I.e '-q 30' is equivalent "
                             "to '-e 0.001'.")
    parser.add_argument("-Q", "--median-quality", type=int,
                        help="The minimum median phred score.")
    parser.add_argument("-c", "--compression-level", type=int,
                        default=DEFAULT_COMPRESSION_LEVEL,
                        help=f"Compression level for the output files. "
                             f"Relevant when output files have a .gz "
                             f"extension. Default: {DEFAULT_COMPRESSION_LEVEL}"
                        )
    return parser


def main():
    args = argument_parser().parse_args()
    filters = []
    # Filters are ordered from low cost to high cost.
    if args.min_length:
        filters.append(MinimumLengthFilter(args.min_length))
    if args.max_length:
        filters.append(MaximumLengthFilter(args.max_length))
    if args.average_error_rate:
        filters.append(AverageErrorRateFilter(args.average_error_rate))
    if args.mean_quality:
        average_error_rate = 10 ** -(args.mean_quality / 10)
        filters.append(AverageErrorRateFilter(average_error_rate))
    if args.median_quality:
        filters.append(MedianQualityFilter(args.median_quality))
    filter_fastq(filters=filters,
                 input_file=args.input,
                 output_file=args.output,
                 compression_level=args.compression_level)


if __name__ == "__main__":  # pragma: no cover
    main()
