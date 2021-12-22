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

from cpython.buffer cimport PyObject_GetBuffer, PyBuffer_Release, PyBUF_SIMPLE

from libc.stdint cimport uint8_t
from libc.string cimport memset, memcmp
from libc.math cimport log10

DEFAULT_PHRED_SCORE_OFFSET = 33
cdef double[128] QUAL_LOOKUP = [10 ** (-0.1 * x) for x in range(128)]

def qualmean(qualities, uint8_t phred_offset = DEFAULT_PHRED_SCORE_OFFSET):
    """
    Returns the mean of the quality scores in an ASCII quality string as stored
    in a FASTQ file.
    """
    cdef Py_buffer buffer
    # Cython makes sure error is handled when acquiring buffer fails.
    PyObject_GetBuffer(qualities, &buffer, PyBUF_SIMPLE)
    cdef uint8_t *scores = <uint8_t *>buffer.buf
    cdef uint8_t score
    cdef double sum_probabilities = 0.0
    cdef double average
    try:
        if buffer.len == 0:
            raise ValueError("Empty quality string")

        for i in range(buffer.len):
            score = scores[i]
            if not (phred_offset <= score < 127):
                raise ValueError(f"Value outside phred range "
                                 f"({phred_offset}-127) detected in qualities: "
                                 f"{repr(qualities)}.")
            sum_probabilities += QUAL_LOOKUP[score - phred_offset]
        average = sum_probabilities / <double>buffer.len
        return -10 * log10(average)
    finally:
        PyBuffer_Release(&buffer)


def qualmean_precise(qualities, uint8_t phred_offset = DEFAULT_PHRED_SCORE_OFFSET):
    """
    Returns the mean of the quality scores in an ASCII quality string as stored
    in a FASTQ file.

    This version minimizes floating point addition by using a histogram. This
    comes at the cost of speed.
    """
    cdef PhredHistogram histogram
    cdef Py_buffer buffer
    cdef double average_probability
    # Cython makes sure error is handled when acquiring buffer fails.
    PyObject_GetBuffer(qualities, &buffer, PyBUF_SIMPLE)

    try:
        if buffer.len == 0:
            raise ValueError("Empty quality string")

        create_phred_histogram(histogram, <uint8_t *>buffer.buf, buffer.len)
        if not histogram_scores_in_phred_range(histogram, phred_offset):
            raise ValueError(f"Value outside phred range "
                             f"({phred_offset}-127) detected in qualities: "
                             f"{repr(qualities)}.")
        average_probability = average_probability_from_phred_histogram(histogram, buffer.len, phred_offset)
        return -10 * log10(average_probability)
    finally:
        PyBuffer_Release(&buffer)
 

ctypedef Py_ssize_t HistogramInt
ctypedef HistogramInt[256] PhredHistogram
cdef const PhredHistogram EMPTY_HIST

def qualmedian(qualities, uint8_t phred_offset = DEFAULT_PHRED_SCORE_OFFSET):
    """
    Returns the median of the quality scores in an ASCII quality string as
    stored in a FASTQ file.
    """

    # To calculate the median this implements half a counting sort:
    # https://en.wikipedia.org/wiki/Counting_sort
    # This is the fastest sort if there is a limited number of discrete values.
    # The complexity is O(n+k) where n is the length of the quality array and
    # k the number of values. Which equals 126 - phred_offset.
    # In counting sort usually a key converter function is needed, but the
    # phred scores themselves can be used as indexes for the count array here.

    # Unlike a counting sort, we are not interested in the actual sorted array
    # only the median. We find the median by finding the point in the counts
    # array where half or over half of the values have been counted.

    cdef PhredHistogram histogram
    cdef Py_buffer buffer
    # Cython makes sure error is handled when acquiring buffer fails.
    PyObject_GetBuffer(qualities, &buffer, PyBUF_SIMPLE)

    try:
        if buffer.len == 0:
            raise ValueError("Empty quality string")

        create_phred_histogram(histogram, <uint8_t *>buffer.buf, buffer.len)
        if not histogram_scores_in_phred_range(histogram, phred_offset):
            raise ValueError(f"Value outside phred range "
                             f"({phred_offset}-127) detected in qualities: "
                             f"{repr(qualities)}.")
        return median_from_phred_histogram(histogram, buffer.len, phred_offset)

    finally:
        PyBuffer_Release(&buffer)
 
cdef object median_from_phred_histogram(HistogramInt *histogram,
                                        Py_ssize_t no_items, 
                                        uint8_t phred_offset):
    cdef bint odd_number_of_items = no_items % 2
    cdef Py_ssize_t half_of_items = no_items // 2  # First middle value of 50 = 25
    if odd_number_of_items:
        half_of_items += 1  # Middle value of 49 = 25 
    cdef Py_ssize_t counted_items = 0
    cdef int i, j
    for i in range(phred_offset, 127):
        counted_items += histogram[i]
        if counted_items >= half_of_items:
            if odd_number_of_items:  # Only one median value
                return i - phred_offset
            if counted_items > half_of_items:
                return i - phred_offset  # The two middle values where the same
            for j in range(i+1, 127):  # Look for the other median value
                if histogram[j] > 0:
                    return <double>(i + j - 2 * phred_offset) / 2.0
    raise RuntimeError("Unable to find median. This is an error in the "
                       "code. Please contact the developers.")
    
cdef void create_phred_histogram(HistogramInt *histogram, uint8_t * scores,
                                 Py_ssize_t length):
    cdef Py_ssize_t i
    # Set all counts to zero
    memset(histogram, 0, sizeof(PhredHistogram))
    for i in range(length):
        histogram[scores[i]] += 1

cdef bint histogram_scores_in_phred_range(HistogramInt *histogram, uint8_t phred_offset):
    """Check if values in the histogram are outside the phred score range."""
    # Memcmp to prepared empty hist seems to be the fastest method.
    # Faster than anything involving array iteration.
    cdef int below_range = memcmp(histogram, EMPTY_HIST, phred_offset * sizeof(HistogramInt))
    cdef int above_range = memcmp(histogram + 127, EMPTY_HIST, (256-127) * sizeof(HistogramInt))
    if below_range != 0 or above_range !=0:
        return 0
    return 1

cdef double average_probability_from_phred_histogram(HistogramInt *histogram,
                                                     Py_ssize_t number_of_items,
                                                     uint8_t phred_offset):
    cdef double sum_probabilities = 0
    cdef int phred_score
    for phred_score in range(phred_offset, 127):
        sum_probabilities += histogram[phred_score] * QUAL_LOOKUP[phred_score - phred_offset]
    return sum_probabilities / <double>number_of_items
