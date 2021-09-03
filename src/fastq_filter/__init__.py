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
import io
import math
import typing
from typing import Generator

import numpy as np


class FastqRecord(typing.NamedTuple):
    name: bytes
    sequence: bytes
    plus: bytes
    phred_scores: bytes


def file_to_fastq_records(filepath) -> Generator[FastqRecord, None, None]:
    with open(filepath, "rb") as file_h:  # type: io.BufferedReader
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
                                 f"of unequal length at byte: "
                                 f"{file_h.tell()}")
            yield FastqRecord(name.rstrip(),
                              sequence.rstrip(),
                              plus.rstrip(),
                              qualities.rstrip())


def mean_quality_filter(record: FastqRecord, quality: float) -> bool:
    phred_scores = np.frombuffer(record.phred_scores, dtype=np.uint8)
    qualities = 10 ^ (phred_scores / -10)
    average = np.average(qualities)
    return -10 * math.log10(average) >= quality


def main():
    pass


if __name__ == "__main__":
    main()
