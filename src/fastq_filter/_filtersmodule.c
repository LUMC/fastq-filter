// Copyright (c) 2021 Leiden University Medical Center
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "structmember.h"

#include "score_to_error_rate.h"
#define MAXIMUM_PHRED_SCORE 126
#define DEFAULT_PHRED_OFFSET 33

/**
 * @brief Returns the average error rate based on an array of phred scores. 
 * 
 * @param phred_scores The array of phred scores.
 * @param phred_length The length of the phred scores array.
 * @param phred_offset The offset for the phred scores
 * @return double The average error rate or -1.0L on error.
 */
static double 
average_error_rate(const uint8_t *phred_scores, size_t phred_length, uint8_t phred_offset)
{
    double total_error_rate = 0.0;
    uint8_t score;
    uint8_t max_score = MAXIMUM_PHRED_SCORE - phred_offset;
    for (size_t i=0; i<phred_length; i+=1) {
        score = phred_scores[i] - phred_offset;
        if (score > max_score) {
            PyErr_Format(
                PyExc_ValueError,
                "Character %c outside of valid phred range ('%c' to '%c')",
                phred_scores[i], phred_offset, MAXIMUM_PHRED_SCORE);
            return -1.0L;
        }
        total_error_rate += SCORE_TO_ERROR_RATE[score];
    }
    return total_error_rate / (double)phred_length;
}

/**
 * @brief Returns a rounded up median phred_score. (i.e. 26.5 -> 27)
 * 
 * @param phred_scores Array of phred scores to calculate the median off.
 * @param phred_length The length of the phred_scores
 * @param phred_offset The phred offset
 * @return ssize_t The median, or -1 on error.
 */
static double
qualmedian(const uint8_t *phred_scores, size_t phred_length, uint8_t phred_offset) {
    size_t histogram[128];
    uint8_t score;
    uint8_t max_score = MAXIMUM_PHRED_SCORE - phred_offset;
    memset(histogram, 0, 128 * sizeof(size_t));
    for (size_t i=0; i < phred_length; i+= 1) {
        score = phred_scores[i] - phred_offset;
        if (score > max_score) {
            PyErr_Format(
                PyExc_ValueError,
                "Character %c outside of valid phred range ('%c' to '%c')",
                phred_scores[i], phred_offset, MAXIMUM_PHRED_SCORE);
            return -1.0L;
        }
        histogram[phred_scores[i]] += 1;
    }
    int odd_number_of_items = phred_length % 2;
    size_t half_of_items = phred_length / 2;  // First middle value of 50 = 25
    if (odd_number_of_items) {
        half_of_items += 1;  // Middle value of 49 = 25
    }
    size_t counted_items = 0;
    for (uint8_t i=0; i <= max_score; i+=1) {
        counted_items += histogram[i];
        if (counted_items >= half_of_items) {
            if (odd_number_of_items) { 
                // Only one median value
                return (double)i;
            }
            if (counted_items > half_of_items) {
                // The two middle values were the same
                return (double)i;
            } 
            for (uint8_t j=i+1; j<max_score; j+=1) {
                if (histogram[j] > 0) {
                    return (double)(i + j) / 2.0L;
                }
            } 
        }
    } 
    PyErr_SetString(PyExc_RuntimeError, 
                    "Unable to find median. This is an error in the code. "
                    "Please contact the developers.");
    return -1.0L;
}

typedef struct {
    PyObject_HEAD
    unsigned long long total; 
    unsigned long long pass;
    double threshold_d;
    Py_ssize_t threshold_i;
    uint8_t phred_offset;
} FastqFilter;

#define GENERIC_FILTER_MEMBERS \
    {"total", T_ULONGLONG, offsetof(FastqFilter, total), READONLY, \
     "the total number of reads checked by this filter"}, \
    {"passed", T_ULONGLONG, offsetof(FastqFilter, pass), READONLY, \
     "the total number of reads to pass this filter"},

static PyMemberDef GenericQualityFilterMembers[] = {
    GENERIC_FILTER_MEMBERS
    {"threshold", T_DOUBLE, offsetof(FastqFilter, threshold_d), READONLY, 
     "The threshold for this filter."},
    {"phred_offset", T_UBYTE, offsetof(FastqFilter, phred_offset), READONLY,
     "The phred offset used for this filter."},
    {NULL}
};

static PyMemberDef GenericLengthFilterMembers[] = {
    GENERIC_FILTER_MEMBERS
    {"threshold", T_PYSSIZET, offsetof(FastqFilter, threshold_i), READONLY, 
     "The threshold for this filter."},
    {NULL}
};

static PyObject *
GenericQualityFilter__new__(PyTypeObject *type, PyObject *args, PyObject *kwargs) 
{
    uint8_t phred_offset = DEFAULT_PHRED_OFFSET;
    double threshold_d = 0.0L;
    static char *kwarg_names[] = {"threshold", "phred_offset", NULL};
    static const char *format = "d|$b:";
    if (!PyArg_ParseTupleAndKeywords(
        args, kwargs, format, kwarg_names,
        &threshold_d,
        &phred_offset)) {
            return NULL;
    }
    FastqFilter *self = PyObject_New(FastqFilter, type);
    self->phred_offset = phred_offset;
    self->threshold_d = threshold_d;
    self->threshold_i = 0;
    self->total = 0;
    self->pass = 0;
    return (PyObject *)self;
}

static PyObject *
GenericLengthFilter__new__(PyTypeObject *type, PyObject *args, PyObject *kwargs) 
{
    uint8_t phred_offset = DEFAULT_PHRED_OFFSET;
    Py_ssize_t threshold_i = 0L;
    static char *kwarg_names[] = {"threshold", NULL};
    static const char *format = "n|:";
    if (!PyArg_ParseTupleAndKeywords(
        args, kwargs, format, kwarg_names,
        &threshold_i)) {
            return NULL;
    }
    FastqFilter *self = PyObject_New(FastqFilter, type);
    self->phred_offset = phred_offset;
    self->threshold_i = threshold_i;
    self->threshold_d = 0.0L;
    self-> total = 0;
    self->pass = 0;
    return (PyObject *)self;
}

static int 
CheckASCIIString(const char *argname, PyObject *arg) {
    if (!PyUnicode_CheckExact(arg)) {
        PyErr_Format(PyExc_TypeError, 
                     "%s must be of type str, got %s",
                     argname, Py_TYPE(arg)->tp_name);
        return 0;
    }
    if (!PyUnicode_IS_COMPACT_ASCII(arg)) {
        PyErr_Format(PyExc_ValueError,
                     "%s must contain only ASCII characters", 
                     argname);
        return 0;
    }
    return 1;
}

PyDoc_STRVAR(GenericQualityFilter_doc,
"passes_filter($self, phred_scores, /)\n"
"--\n"
"\n"
"Check if the phred scores pass the filter.\n"
"\n"
"  phred_scores\n"
"    phred_scores as an ASCII-encoded string.\n"
"\n");

static PyObject *
AverageErrorRateFilter_passes_filter(FastqFilter *self, PyObject *phred_scores) 
{
    if (!CheckASCIIString("phred_scores", phred_scores)) {
        return NULL;
    }
    uint8_t *phreds = PyUnicode_DATA(phred_scores);
    Py_ssize_t phred_length = PyUnicode_GET_LENGTH(phred_scores);
    double error_rate = average_error_rate(phreds, phred_length, self->phred_offset);
    if (error_rate < 0) {
        return NULL; 
    }
    self->total += 1;
    int pass = error_rate >= self->threshold_d;
    if (pass) {
        self->pass += 1;
    }
    return PyBool_FromLong(pass);
}

static PyObject *
MedianQualityFilter_passes_filter(FastqFilter *self, PyObject *phred_scores)
{
    if (!CheckASCIIString("phred_scores", phred_scores)) {
        return NULL;
    }
    uint8_t *phreds = PyUnicode_DATA(phred_scores);
    Py_ssize_t phred_length = PyUnicode_GetLength(phred_scores);
    double median = qualmedian(phreds, phred_length, self->phred_offset);
    if (median < 0) {
        return NULL;
    }
    self->total += 1;
    int pass = median >= self->threshold_d;
    if (pass) {
        self->pass += 1;
    }
    return PyBool_FromLong(pass);
}

PyDoc_STRVAR(GenericSequenceFilter_doc,
"passes_filter($self, sequence, /)\n"
"--\n"
"\n"
"Check if the sequence passes the filter.\n"
"\n"
"  sequence\n"
"    The nucleotide sequence as an ASCII-encoded string.\n"
"\n");
static PyObject * 
MinLengthFilter_passes_filter(FastqFilter *self, PyObject *sequence)
{
    if (!CheckASCIIString("sequence", sequence)) {
        return NULL;
    }
    Py_ssize_t length = PyUnicode_GET_LENGTH(sequence);
    int pass = length >= self->threshold_i;
    self->total += 1;
    if (pass) {
        self->pass += 1;
    }
    return PyBool_FromLong(pass);
}

static PyObject *
MaxLengthFilter_passes_filter(FastqFilter *self, PyObject *sequence)
{
    if (!CheckASCIIString("sequence", sequence)) {
        return NULL;
    }
    Py_ssize_t length = PyUnicode_GET_LENGTH(sequence);
    int pass = length <= self->threshold_i;
    self->total += 1;
    if (pass) {
        self->pass += 1;
    }
    return PyBool_FromLong(pass);   
}

static PyMethodDef AverageErrorRateFilter_methods[] = {
    {"passes_filter", (PyCFunction)AverageErrorRateFilter_passes_filter, METH_O,
     GenericQualityFilter_doc},
    {NULL}, 
};

static PyMethodDef MedianQualityFilter_methods[] = {
    {"passes_filter", (PyCFunction)MedianQualityFilter_passes_filter, METH_O,
     GenericQualityFilter_doc},
    {NULL}, 
};

static PyMethodDef MinLengthFilter_methods[] = {
    {"passes_filter", (PyCFunction)MinLengthFilter_passes_filter, METH_O,
     GenericSequenceFilter_doc},
    {NULL}, 
};

static PyMethodDef MaxLengthFilter_methods[] = {
    {"passes_filter", (PyCFunction)MaxLengthFilter_passes_filter, METH_O,
     GenericSequenceFilter_doc},
    {NULL}, 
};


static PyTypeObject AverageErrorRateFilter_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_filter.AverageErrorRateFilter",
    .tp_basicsize = sizeof(FastqFilter),
    .tp_new = GenericQualityFilter__new__,
    .tp_methods = AverageErrorRateFilter_methods,
    .tp_members = GenericQualityFilterMembers,
};

static struct PyModuleDef _filters_module = {
    PyModuleDef_HEAD_INIT,
    "_filters",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,
    NULL  /* module methods */
};


PyMODINIT_FUNC
PyInit__filters(void)
{
    PyObject *m;

    m = PyModule_Create(&_filters_module);
    if (m == NULL) {
        return NULL;
    }

    PyTypeObject *AverageErrorRateFilter = &AverageErrorRateFilter_Type;
    if (PyType_Ready(AverageErrorRateFilter) == -1) {
        return NULL;
    } 
    PyModule_AddObject(m, "AverageErrorRateFilter", 
                       (PyObject *)AverageErrorRateFilter);

    PyModule_AddIntMacro(m, DEFAULT_PHRED_OFFSET);
    return m;
}