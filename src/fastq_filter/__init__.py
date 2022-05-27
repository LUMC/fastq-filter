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
import contextlib
import functools
from typing import Callable, Iterable, Iterator, List, Tuple

import dnaio

import xopen  # type: ignore

from ._filters import (
    DEFAULT_PHRED_SCORE_OFFSET,
    average_error_rate,
    qualmean,
    qualmedian
)

__version__ = "0.2.0"

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


def file_to_fastq_records(filepath: str) -> Iterator[dnaio.Sequence]:
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


def multiple_files_to_records(input_files: List[str],
                              ) -> Iterator[Tuple[dnaio.SequenceRecord, ...]]:
    readers = [file_to_fastq_records(f) for f in input_files]
    iterators = [iter(reader) for reader in readers]
    for records in zip(*iterators):
        if len(records) > 1 and not dnaio.records_are_mates(*records):
            raise dnaio.FastqFormatError(
                f"Records are out of sync, names "
                f"{', '.join(r.name for r in records)} do not match.",
                line=None
            )
        yield records
    # Check if all iterators are exhausted.
    for iterator in iterators:
        try:
            _ = next(iterator)
            raise dnaio.FastqFormatError("Input files have an unequal number"
                                         " of FASTQ records.", line=None)
        except StopIteration:
            pass


def filter_fastq(input_files: List[str], output_files: List[str],
                 filters: List[Callable[[Tuple[dnaio.SequenceRecord, ...]], bool]],
                 compression_level: int = DEFAULT_COMPRESSION_LEVEL):
    """
    Filter FASTQ input files with the filters in filters and write
    the results to the output file.

    :param filters: Functions that filter a tuple of dnaio.sequence records
    :param input_files: FASTQ input filenames. Compressed files are handled
    automatically.
    :param output_files: FASTQ output filenames. Compressed files are handled
    automatically.
    :param compression_level: Compression level for the output files (if
    applicable)
    """
    if len(input_files) != len(output_files):
        raise ValueError("Number of inputs and outputs should be equal.")
    filtered_fastq_records = multiple_files_to_records(input_files)
    for filter_func in filters:
        filtered_fastq_records = filter(filter_func, filtered_fastq_records)
    with contextlib.ExitStack() as output_stack:
        outputs = [output_stack.enter_context(
                   xopen.xopen(output_file, threads=0, mode="wb",
                               compresslevel=compression_level))
                   for output_file in output_files]
        for records in filtered_fastq_records:
            for record, output in zip(records, outputs):
                output.write(record.fastq_bytes())


def MinimumLengthFilter(threshold: int):
    def filterfunc(record: dnaio.SequenceRecord):
        return len(record) >= threshold

    def combined_filter(records: Tuple[dnaio.SequenceRecord, ...]):
        return any(map(filterfunc, records))

    return combined_filter


def MaximumLengthFilter(threshold: int):
    def filterfunc(record: dnaio.SequenceRecord):
        return len(record) <= threshold

    def combined_filter(records: Tuple[dnaio.SequenceRecord, ...]):
        return all(map(filterfunc, records))

    return combined_filter


def AverageErrorRateFilter(threshold: float, phred_offset=DEFAULT_PHRED_SCORE_OFFSET):
    def combined_filter(records: Tuple[dnaio.SequenceRecord, ...]):
        phred_scores = "".join([record.qualities for record in records
                                if record.qualities is not None])
        return average_error_rate(phred_scores, phred_offset) <= threshold
    return combined_filter


def MedianQualityFilter(threshold: float, phred_offset=DEFAULT_PHRED_SCORE_OFFSET):
    def combined_filter(records: Tuple[dnaio.SequenceRecord, ...]):
        phred_scores = "".join([record.qualities for record in records
                                if record.qualities is not None])
        return qualmedian(phred_scores, phred_offset) >= threshold
    return combined_filter


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.description = "Filter FASTQ files on various metrics."
    parser.add_argument("input",
                        help="Input FASTQ files. Compression format "
                             "automatically detected. Use - for stdin.",
                        nargs='+')
    parser.add_argument("-o", "--output",
                        help="Output FASTQ files. Compression format "
                             "automatically determined by file extension. "
                             "Flag can be used multiple times. An output must "
                             "be given for each input. Default: stdout.",
                        action='append')
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
    output = args.output if args.output else ["-"]
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
                 input_files=args.input,
                 output_files=output,
                 compression_level=args.compression_level)


if __name__ == "__main__":  # pragma: no cover
    main()
