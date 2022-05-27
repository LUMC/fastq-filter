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
import array
import itertools
import math
import statistics
import sys
from typing import List

from dnaio import Sequence

import fastq_filter
from fastq_filter import DEFAULT_PHRED_SCORE_OFFSET, qualmean, qualmedian

import pytest  # type: ignore


def quallist_to_string(quallist: List[int]):
    return array.array(
        "B", [qual + DEFAULT_PHRED_SCORE_OFFSET for qual in quallist]
    ).tobytes().decode("ascii")


QUAL_STRINGS = [
    "I?>DC:>@?IDC9??G?>EH9E@66=9<?@E?DC:@<@BBFG>=FIC@F9>7CG?IC?I;CD9>>>A@C7>>"
    "8>>D9GCB<;?DD>C;9?>5G>?H?=6@>:G6B<?==A7?@???8IF<75C=@A:BEA@A;C89D:=1?=<A"
    ">D=>B66C",
    "C:@?;8@=DC???>E>E;98BBB?9D=?@B;D?I:??FD8CH?A7?<H>ABD@C@C?>;;B<><;9@8BAFD"
    "?;:>I3DB<?<B=?A??CI>2E>><BD?A??FCBCE?DAI><B:8D>?C>@BA=F<>7=E=?DC=@9GG=>?"
    "C@><CA;>",
]


@pytest.mark.parametrize("qualstring", QUAL_STRINGS)
def test_qualmean(qualstring):
    offset = DEFAULT_PHRED_SCORE_OFFSET
    qualities = [qual - offset for qual in
                 array.array("b", qualstring.encode('ascii'))]
    probabilities = [10 ** (qual / -10) for qual in qualities]
    average_prob = statistics.mean(probabilities)
    phred = - 10 * math.log10(average_prob)
    assert phred == pytest.approx(qualmean(qualstring))


@pytest.mark.parametrize("qualstring", QUAL_STRINGS)
def test_qualmedian(qualstring):
    offset = DEFAULT_PHRED_SCORE_OFFSET
    qualities = [qual - offset for qual in
                 array.array("b", qualstring.encode('ascii'))]
    median_quality = statistics.median(qualities)
    assert median_quality == qualmedian(qualstring)


def test_qualmedian_correct():
    # Make sure qualmedian also returns averages.
    qualities = "AACEGG"  # Median value should be D. ord("D") == 68
    result = qualmedian(qualities, 0)
    assert result == 68.0
    assert type(result) == float


TOO_LOW_PHREDS = [chr(x) for x in range(33)]
TOO_HIGH_PHREDS = [chr(127)]
OUTSIDE_RANGE_PHREDS = TOO_LOW_PHREDS + TOO_HIGH_PHREDS
NON_ASCII_PHREDS = [chr(x) for x in range(128, 256)]


@pytest.mark.parametrize(["func", "quals"],
                         itertools.product([qualmean, qualmedian],
                                           OUTSIDE_RANGE_PHREDS))
def test_outside_range_phreds(func, quals):
    with pytest.raises(ValueError) as error:
        func(quals)
    assert error.match("outside of valid phred range")


@pytest.mark.parametrize(["func", "quals"],
                         itertools.product([qualmean, qualmedian],
                                           NON_ASCII_PHREDS))
def test_non_ascii_phreds(func, quals):
    with pytest.raises(ValueError) as error:
        func(quals)
    assert error.match("phred_scores must be ASCII encoded.")


def test_fastq_records_to_file(tmp_path):
    records = [Sequence("TEST", "A", "A")] * 3
    out = tmp_path / "test.fq"
    fastq_filter.fastq_records_to_file(records, str(out))
    assert out.read_bytes() == b"@TEST\nA\n+\nA\n" \
                               b"@TEST\nA\n+\nA\n" \
                               b"@TEST\nA\n+\nA\n"


def test_file_to_fastq_records(tmp_path):
    out = tmp_path / "test.fq"
    out.write_bytes(b"@TEST\nA\n+\nA\n@TEST\nA\n+\nA\n@TEST\nA\n+\nA\n")
    assert list(fastq_filter.file_to_fastq_records(str(out))) == [
        Sequence("TEST", "A", "A")] * 3


def test_main(tmp_path):
    in_f = tmp_path / "in.fq"
    out_f = tmp_path / "out.fq"
    in_f.write_bytes(b"@TEST\nAA\n+\nAA\n@TEST\nA\n+\n-\n@TEST\nA\n+\nA\n")
    sys.argv = ["", "-q", "20", "-l", "2", str(in_f), "-o", str(out_f)]
    fastq_filter.main()
    assert out_f.read_bytes() == b"@TEST\nAA\n+\nAA\n"


@pytest.mark.parametrize("func", [qualmean, qualmedian])
def test_empty_quals_returns_nan(func):
    assert math.isnan(func(""))
