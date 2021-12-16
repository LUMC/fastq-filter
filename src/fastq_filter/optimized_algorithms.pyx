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
cdef const Py_ssize_t[256] EMPTY_HIST
cdef double[256] QUAL_LOOKUP = [10 ** (-0.1 * x) for x in range(256)]

cdef void create_histogram(Py_ssize_t * histogram, uint8_t * scores,
                           Py_ssize_t length):
    cdef Py_ssize_t i
    # Set all counts to zero
    memset(histogram, 0, sizeof(Py_ssize_t) * 256)
    for i in range(length):
        histogram[scores[i]] += 1

cdef bint histogram_scores_in_phred_range(Py_ssize_t * histogram, uint8_t phred_offset):
    """Check if values in the histogram are outside the phred score range."""
    # Memcmp to prepared empty hist seems to be the fastest method.
    # Faster than anything involving array iteration.
    cdef int below_range = memcmp(histogram, EMPTY_HIST, phred_offset * sizeof(Py_ssize_t))
    cdef int above_range = memcmp(histogram + 127, EMPTY_HIST, (256-127) * sizeof(Py_ssize_t))
    if below_range != 0 or above_range !=0:
        return 0
    return 1

def qualmean(qualities, double phred_offset = DEFAULT_PHRED_SCORE_OFFSET):
    """
    Returns the mean of the quality scores in an ASCII quality string as stored
    in a FASTQ file.
    """
    # All intermediate calculations are calculated with single precision as
    # this saves a double power (**) and double addition calculation.
    # This is 4 times faster on the PC of this developer.
    # For the average_phred, double values are used since Python uses doubles
    # internally and this prevents casting.
    cdef double sum_probabilities = 0.0
    cdef double average
    cdef double average_phred
    cdef Py_buffer buffer_data
    cdef Py_buffer* buffer = &buffer_data
    # Cython makes sure error is handled when acquiring buffer fails.
    PyObject_GetBuffer(qualities, buffer, PyBUF_SIMPLE)
    cdef uint8_t *scores = <uint8_t *>buffer.buf
    cdef uint8_t score
    try:
        if buffer.len == 0:
            raise ValueError("Empty quality string")
        for i in range(buffer.len):
            score = scores[i]
            if not (phred_offset <= score < 127):
                raise ValueError(f"Value outside phred range "
                                 f"({phred_offset}-127) detected in qualities: "
                                 f"{repr(qualities)}.")
            sum_probabilities += QUAL_LOOKUP[score]
        average = sum_probabilities / <double>buffer.len
        average_phred = -10 * log10(average) - phred_offset
        return average_phred
    finally:
        PyBuffer_Release(buffer)


def qualmedian(qualities, int phred_offset = DEFAULT_PHRED_SCORE_OFFSET):
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

    cdef Py_ssize_t[256] counts
    cdef Py_ssize_t total = 0
    cdef Py_buffer buffer_data
    cdef Py_ssize_t half
    cdef bint odd
    cdef Py_ssize_t i
    cdef int j
    cdef int k
    cdef Py_buffer* buffer = &buffer_data
    # Cython makes sure error is handled when acquiring buffer fails.
    PyObject_GetBuffer(qualities, buffer, PyBUF_SIMPLE)
    cdef uint8_t *scores = <uint8_t *>buffer.buf

    try:
        if buffer.len == 0:
            raise ValueError("Empty quality string")

        odd = buffer.len % 2
        if odd:
            # If len = 49, we have reached the middle value at count 25
            half = buffer.len // 2 + 1
        else:
            # If len = 50, we have reached the first of the two middle values
            # at count 25
            half = buffer.len // 2

        create_histogram(counts, scores, buffer.len)
        if not histogram_scores_in_phred_range(counts, phred_offset):
            raise ValueError(f"Value outside phred range "
                             f"({phred_offset}-127) detected in qualities: "
                             f"{repr(qualities)}.")
        # Check at which value of j, we have counted half of the values.
        for j in range(phred_offset, 127):
            total += counts[j]
            if total >= half:
                if odd:  # There is only one value
                    return j - phred_offset
                if total > half:  # The two middle values were the same.
                    return j - phred_offset
                # The highest middle value is higher than j.
                for k in range(j + 1, 127):
                    if counts[j] > 0:
                        # Cast to double to prevent integer scores here.
                        return (<double>(j + k) / 2.0 - <double>phred_offset)
        raise RuntimeError("Unable to find median. This is an error in the "
                           "code. Please contact the developers.")
    finally:
        PyBuffer_Release(buffer)
