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

from cpython.buffer cimport PyObject_GetBuffer, PyBuffer_Release, PyBUF_C_CONTIGUOUS

from libc.stdint cimport uint8_t, uint_fast32_t
from libc.math cimport log10

from .constants import DEFAULT_PHRED_SCORE_OFFSET


def qualmean(bytes qualities, double phred_offset = DEFAULT_PHRED_SCORE_OFFSET):
    """
    Returns the mean of the quality scores in an ASCII quality string as stored
    in a FASTQ file. Returns a float value.
    """
    # All intermediate calculations are calculated with single precision as
    # this saves a double power (**) and double addition calculation.
    # This is 4 times faster on the PC of this developer.
    # For the average_phred, double values are used since Python uses doubles
    # internally and this prevents casting.
    cdef float phred_constant = 10 ** -0.1
    cdef float sum_probabilities = 0.0
    cdef float average
    cdef double average_phred
    cdef Py_buffer buffer_data
    cdef Py_buffer* buffer = &buffer_data
    # Cython makes sure error is handled when acquiring buffer fails.
    PyObject_GetBuffer(qualities, buffer, PyBUF_C_CONTIGUOUS)
    cdef uint8_t *scores = <uint8_t *>buffer.buf
    try:
        if buffer.len == 0:
            raise ValueError("Empty quality string")
        for i in range(buffer.len):
            sum_probabilities += phred_constant ** scores[i]
        average = sum_probabilities / buffer.len
        average_phred = -10 * log10(<double>average) - phred_offset
        return average_phred
    finally:
        PyBuffer_Release(buffer)
