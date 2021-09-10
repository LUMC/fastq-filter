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

import math
import statistics

from .constants import DEFAULT_PHRED_SCORE_OFFSET


def qualmean(qualities: bytes, phred_offset: int = DEFAULT_PHRED_SCORE_OFFSET
             ) -> float:
    phred_constant = 10 ** -0.1
    sum_probabilities = 0.0
    for score in qualities:
        sum_probabilities += phred_constant ** score
    average = sum_probabilities / len(qualities)
    return -10 * math.log10(average) - phred_offset


def qualmedian(qualities: bytes, phred_offset: int = DEFAULT_PHRED_SCORE_OFFSET
               ) -> float:
    """Calculate the median phred score from a raw FASTQ quality string."""
    return statistics.median(qualities) - phred_offset