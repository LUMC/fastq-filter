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
from typing import List

from dnaio import SequenceRecord

from fastq_filter import (
    AverageErrorRateFilter,
    DEFAULT_PHRED_SCORE_OFFSET,
    MaximumLengthFilter,
    MedianQualityFilter,
    MinimumLengthFilter
)

import pytest


def quallist_to_string(quallist: List[int]):
    return array.array(
        "B", [qual + DEFAULT_PHRED_SCORE_OFFSET for qual in quallist]
    ).tobytes().decode("ascii")


# def test_average_error_rate_filter_new():
#     filter = AverageErrorRateFilter(0.001, phred_offset=20)
#     assert filter.threshold == 0.001
#     assert filter.phred_offset == 20
#
# def test_median_quality_filter_new():
#     filter = MedianQualityFilter(0.001, phred_offset=20)
#     assert filter.threshold == 0.001
#     assert filter.phred_offset == 20
#
#
# def test_minimum_length_filter_new():
#     filter = MinimumLengthFilter(20)
#     assert filter.threshold == 20
#
#
# def test_maximum_length_filter_new():
#     filter = MaximumLengthFilter(20)
#     assert filter.threshold == 20


@pytest.mark.parametrize(
    ["threshold", "qualities", "result"], (
        (0.001, [chr(63)], True),
        (0.001, [chr(64)], True),
        (0.001, [chr(62)], False),
        (10 ** -(10 / 10), [quallist_to_string([9, 9, 9])], False),
        (10 ** -(8 / 10), [quallist_to_string([9, 9, 9])], True),
        (10 ** -(8 / 10), [quallist_to_string([8, 8, 8]),
                           quallist_to_string([8, 7, 8])], False),
    ))
def test_average_error_rate_filter(threshold, qualities, result):
    filter = AverageErrorRateFilter(threshold)
    records = [SequenceRecord("name", len(qual) * 'A', qual)
               for qual in qualities]
    assert filter(tuple(records)) is result
    # assert filter.total == 1
    # if result:
    #     assert filter.passed == 1
    # else:
    #     assert filter.passed == 0


@pytest.mark.parametrize(
    ["threshold", "qualities", "result"], (
        (30, [chr(63)], True),
        (30, [chr(64)], True),
        (30, [chr(62)], False),
        (10, [quallist_to_string([9, 9, 9, 10, 10])], False),
        (8, [quallist_to_string([9, 9, 9])], True),
        (8, [quallist_to_string([1, 1, 1, 8, 9, 9, 9])], True),
        (8, [quallist_to_string([1, 1, 1, 8, 9, 9, 9]),
             quallist_to_string([1, 1, 1, 7, 9, 9, 9])], False),
        (8, [quallist_to_string([1, 1, 1, 8, 9, 9, 9]),
             quallist_to_string([1, 1, 1, 8, 9, 9, 9])], True)
    ))
def test_median_quality_filter(threshold, qualities, result):
    filter = MedianQualityFilter(threshold)
    records = [SequenceRecord("name", len(qual) * 'A', qual)
               for qual in qualities]
    assert filter(tuple(records)) is result
    # assert filter.total == 1
    # if result:
    #     assert filter.passed == 1
    # else:
    #     assert filter.passed == 0


TOO_LOW_PHREDS = [chr(x) for x in range(33)]
TOO_HIGH_PHREDS = [chr(127)]
OUTSIDE_RANGE_PHREDS = TOO_LOW_PHREDS + TOO_HIGH_PHREDS


@pytest.mark.parametrize(
    ["filter_class", "quals"],
    itertools.product(
        [AverageErrorRateFilter, MedianQualityFilter], OUTSIDE_RANGE_PHREDS)
)
def test_outside_range(filter_class, quals):
    record = SequenceRecord("name", "A", quals)
    filter = filter_class(1)
    with pytest.raises(ValueError) as error:
        filter((record,))
    error.match("outside of valid phred range")


@pytest.mark.parametrize(
    ["threshold", "lengths", "result"], (
        (10, [10], True),
        (11, [10], True),
        (10, [11], False),
        (10, [11, 10], False),
        (10, [10, 9], True)
    ))
def test_maximum_length_filter(threshold, lengths, result):
    filter = MaximumLengthFilter(threshold)
    records = [SequenceRecord("name", length * 'A', length * 'H')
               for length in lengths]
    assert filter(tuple(records)) is result
    # assert filter.total == 1
    # if result:
    #     assert filter.passed == 1
    # else:
    #     assert filter.passed == 0


@pytest.mark.parametrize(
    ["threshold", "lengths", "result"], (
        (10, [10], True),
        (11, [10], False),
        (10, [11], True),
        (10, [11, 10], True),
        (10, [10, 9], True),
        (10, [10, 8, 8], True),
        (10, [9, 9], False)
    ))
def test_minimum_length_filter(threshold, lengths, result):
    filter = MinimumLengthFilter(threshold)
    records = [SequenceRecord("name", length * 'A', length * 'H')
               for length in lengths]
    assert filter(tuple(records)) is result
    # assert filter.total == 1
    # if result:
    #     assert filter.passed == 1
    # else:
    #     assert filter.passed == 0
