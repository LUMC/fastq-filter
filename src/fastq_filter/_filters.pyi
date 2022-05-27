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

from typing import Union

from dnaio import SequenceRecord

DEFAULT_PHRED_SCORE_OFFSET: int = ...

class _Filter:
    threshold: Union[int, float]
    passed: int
    total: int

    def __init__(self): ...

    def __call__(self, __record: SequenceRecord) -> bool: ...

class _QualityFilter(_Filter):
    phred_offset: int

    def __init__(self, threshold: float,
                 phred_offset: int = DEFAULT_PHRED_SCORE_OFFSET): ...


class _LengthFilter(_Filter):
    def __init__(self, threshold: int): ...


class AverageErrorRateFilter(_QualityFilter): ...
class MedianQualityFilter(_QualityFilter): ...
class MinimumLengthFilter(_LengthFilter): ...
class MaximumLengthFilter(_LengthFilter): ...

def qualmean(phred_scores: str, phred_offset: int = DEFAULT_PHRED_SCORE_OFFSET): ...

def qualmedian(phred_scores: str, phred_offset: int = DEFAULT_PHRED_SCORE_OFFSET): ...

def average_error_rate(phred_scores: str,
                       phred_offset: int = DEFAULT_PHRED_SCORE_OFFSET): ...
