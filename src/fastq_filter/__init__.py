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
import math
import typing
from typing import Generator

import numpy as np
from numpy import typing as npt

import xopen


PHRED_SCORE_OFFSET = 33


class FastqRecord(typing.NamedTuple):
    name: bytes
    sequence: bytes
    plus: bytes
    qualities: bytes

    def phred_scores(self) -> npt.NDArray:
        return np.frombuffer(self.qualities, dtype=np.uint8) - \
               PHRED_SCORE_OFFSET


def qualmean(phred_scores: npt.NDArray) -> float:
    # 10 ** (phred_scores / -10) -> 10 ** (phred_scores * -0.1)
    # a ** (p * q) == (a**p)**q according to the math rules.
    # So we can do (10 ** -0.1) ** phred_scores.
    # That way we only do one calculation on the phred_score array instead of
    # two, which improves performance.
    qualities = np.power((10 ** -0.1), phred_scores)
    average = np.average(qualities)
    return -10 * math.log10(average)


def file_to_fastq_records(filepath) -> Generator[FastqRecord, None, None]:
    with xopen.xopen(filepath, "rb", threads=0) as file_h:
        while True:
            name = file_h.readline()
            if name == b"":
                return
            if not name.startswith(b"@"):
                raise ValueError("Record header should start with '@'.")
            try:
                sequence = next(file_h)
                plus = next(file_h)
                qualities = next(file_h)
            except StopIteration:
                raise ValueError("Truncated fastq record at EOF")
            if len(sequence) != len(qualities):
                raise ValueError(f"Fastq record with sequence and qualities "
                                 f"of unequal length at record: "
                                 f"{name.rstrip()}")
            yield FastqRecord(name.rstrip(),
                              sequence.rstrip(),
                              plus.rstrip(),
                              qualities.rstrip())


def mean_quality_filter(quality: float, record: FastqRecord) -> bool:
    return qualmean(record.phred_scores()) >= quality


def median_quality_filter(quality: int, record: FastqRecord) -> bool:
    return np.median(record.phred_scores()) >= quality


def min_length_filter(min_length: int, record: FastqRecord) -> bool:
    return len(record.sequence) >= min_length


def max_length_filter(max_length: int, record: FastqRecord) -> bool:
    return len(record.sequence) <= max_length


# Store filter names for use on the command line interface. Also store a
# tuple of types so the command line arguments (strings) can be converted
# in the appropiate types.
FILTERS = {"mean_quality": (mean_quality_filter, (float,)),
           "median_quality": (median_quality_filter, (int,)),
           "min_length": (min_length_filter, (int,)),
           "max_length": (max_length_filter, (int,))}


def argument_parser() -> argparse.ArgumentParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("filter", nargs="+",
                        help="Filter and arguments. For example: "
                             "mean_quality:20, for filtering all reads with "
                             "an average quality below 20. Available filters: "
                             "mean_quality:<quality>, "
                             "median_quality:<quality>, "
                             "min_length:<length>, "
                             "max_length:<length>. "
                             "Multiple filters can be used. Make sure to use "
                             "faster filters (length) before slower ones "
                             "(quality) for optimal performance.")
    parser.add_argument("input",
                        help="Input FASTQ file. Compression format "
                             "automatically detected. ")
    parser.add_argument("output",
                        help="Output FASTQ file. Compression format "
                             "automatically determined by file extension.")
    return parser


def main():
    args = argument_parser().parse_args()
    fastq_records = file_to_fastq_records(args.input)
    filtered_fastq_records = fastq_records
    for filter_string in args.filter:
        filter_name, filter_argstring = filter_string.split(':')
        try:
            filter_function, filter_argtypes = FILTERS[filter_name]
        except KeyError:
            raise ValueError(f"Unknown filter: {filter_name}. Choose one of:"
                             f" {' '.join(FILTERS.keys())}")

        # Convert the strings from the command line in the appropiate types
        filter_args = [filter_argtypes[pos](arg) for pos, arg
                       in enumerate(filter_argstring.split(','))]

        filter_function_with_args = functools.partial(
            filter_function, *filter_args)

        # Add the filter. Due to this setup we can add multiple filters
        filtered_fastq_records = filter(filter_function_with_args,
                                        filtered_fastq_records)

    with xopen.xopen(args.output, "wb", threads=0) as output_h:
        for record in filtered_fastq_records:
            output_h.write(b"\n".join(record) + b"\n")


if __name__ == "__main__":
    main()
