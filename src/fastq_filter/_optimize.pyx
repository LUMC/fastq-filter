# cython: binding=True
# cython: language_level=3
# cython: cdivision=True
import math
from cpython.buffer cimport PyObject_GetBuffer, PyBuffer_Release, PyBUF_C_CONTIGUOUS

from libc.stdint cimport uint8_t, uint_fast32_t
from libc.math cimport log10

def qualmean(bytes qualities, double phred_offset = 33):
    """
    Returns the mean of the quality scores in an ASCII quality string as stored
    in a FASTQ file. Returns a float value.
    """
    cdef double phred_constant = 10 ** -0.1
    cdef double sum_probabilities = 0.0
    cdef double average
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
        average_phred = -10 * log10(average) - phred_offset
        return average_phred
    finally:
        PyBuffer_Release(buffer)
