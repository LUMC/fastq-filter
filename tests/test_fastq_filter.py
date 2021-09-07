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
import math
import statistics

from fastq_filter import DEFAULT_PHRED_SCORE_OFFSET, qualmean, qualmedian

import pytest  # type: ignore

QUAL_STRINGS = [
    b"I?>DC:>@?IDC9??G?>EH9E@66=9<?@E?DC:@<@BBFG>=FIC@F9>7CG?IC?I;CD9>>>A@C7>>"
    b"8>>D9GCB<;?DD>C;9?>5G>?H?=6@>:G6B<?==A7?@???8IF<75C=@A:BEA@A;C89D:=1?=<A"
    b">D=>B66C",
    b"C:@?;8@=DC???>E>E;98BBB?9D=?@B;D?I:??FD8CH?A7?<H>ABD@C@C?>;;B<><;9@8BAFD"
    b"?;:>I3DB<?<B=?A??CI>2E>><BD?A??FCBCE?DAI><B:8D>?C>@BA=F<>7=E=?DC=@9GG=>?"
    b"C@><CA;>",
]


@pytest.mark.parametrize("qualstring", QUAL_STRINGS)
def test_qualmean(qualstring):
    offset = DEFAULT_PHRED_SCORE_OFFSET
    qualities = [qual - offset for qual in array.array("b", qualstring)]
    probabilities = [10 ** (qual / -10) for qual in qualities]
    average_prob = statistics.mean(probabilities)
    phred = - 10 * math.log10(average_prob)
    assert phred == pytest.approx(qualmean(qualstring))


@pytest.mark.parametrize("qualstring", QUAL_STRINGS)
def test_qualmedian(qualstring):
    offset = DEFAULT_PHRED_SCORE_OFFSET
    qualities = [qual - offset for qual in array.array("b", qualstring)]
    median_quality = statistics.median(qualities)
    assert median_quality == pytest.approx(qualmedian(qualstring))
