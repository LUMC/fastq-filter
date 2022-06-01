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
#include <math.h>

#include "structmember.h"

#include "score_to_error_rate.h"
#define MAXIMUM_PHRED_SCORE 126
#define DEFAULT_PHRED_SCORE_OFFSET 33

static PyTypeObject *SequenceRecord = NULL;
static PyObject *QualitiesAttrString = NULL;


static inline double 
sum_error_rate(const uint8_t *phred_scores, size_t phred_length, uint8_t phred_offset) {
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
    return total_error_rate;
}

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
    double total_error_rate = sum_error_rate(phred_scores, phred_length, phred_offset);
    if (total_error_rate < 0.0) {
        return -1.0L;
    }
    return total_error_rate / (double)phred_length;
}

static inline int 
make_histogram(size_t *histogram, const uint8_t *phred_scores, size_t phred_length, uint8_t phred_offset) 
{
    uint8_t score;
    uint8_t max_score = MAXIMUM_PHRED_SCORE - phred_offset;

    for (size_t i=0; i < phred_length; i+= 1) {
        score = phred_scores[i] - phred_offset;
        if (score > max_score) {
            PyErr_Format(
                PyExc_ValueError,
                "Character %c outside of valid phred range ('%c' to '%c')",
                phred_scores[i], phred_offset, MAXIMUM_PHRED_SCORE);
            return -1;
        }
        histogram[score] += 1;
    }
    return 0;
}

static inline double 
median_from_histogram(size_t *histogram, size_t phred_length, uint8_t phred_offset) 
{
    uint8_t max_score = MAXIMUM_PHRED_SCORE - phred_offset;

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

/**
 * @brief Returns a rounded up median phred_score. (i.e. 26.5 -> 27)
 * 
 * @param phred_scores Array of phred scores to calculate the median off.
 * @param phred_length The length of the phred_scores
 * @param phred_offset The phred offset
 * @return double The median, or -1.0 on error.
 */
static double
qualmedian(const uint8_t *phred_scores, size_t phred_length, uint8_t phred_offset)
{
    if (phred_length == 0) {
        return NAN;
    }
    size_t histogram[128];
    memset(histogram, 0, 128 * sizeof(size_t));
    int ret = make_histogram(histogram, phred_scores, phred_length, phred_offset);
    if (ret != 0) {
        return -1.0L;
    }
    return median_from_histogram(histogram, phred_length, phred_offset);
}

PyDoc_STRVAR(qualmean__doc__,
"qualmean($self, phred_scores, /, phred_offset=DEFAULT_PHRED_SCORE_OFFSET)\n"
"--\n"
"\n"
"Returns the mean quality score. \n"
"\n"
"  phred_scores\n"
"    ASCII string with the phred scores.\n"
);

#define QUALMEAN_METHODDEF    \
    {"qualmean", (PyCFunction)(void(*)(void))qualmean, \
     METH_VARARGS | METH_KEYWORDS, qualmean__doc__}

static PyObject *
qualmean(PyObject *module, PyObject *args, PyObject *kwargs)
{
    PyObject *phred_scores = NULL;
    uint8_t phred_offset = DEFAULT_PHRED_SCORE_OFFSET;
    char *kwarg_names[] = {"", "phred_offset", NULL};
    const char *format = "O!|b:qualmean";
    if (!PyArg_ParseTupleAndKeywords(
        args, kwargs, format, kwarg_names,
        &PyUnicode_Type,
        &phred_scores,
        &phred_offset)) {
            return NULL;
    }

    if (!PyUnicode_IS_COMPACT_ASCII(phred_scores)) {
        PyErr_SetString(PyExc_ValueError,
                        "phred_scores must be ASCII encoded.");
        return NULL;
    }
    uint8_t *phreds = PyUnicode_DATA(phred_scores);
    size_t phred_length = PyUnicode_GET_LENGTH(phred_scores);
    double error_rate = average_error_rate(phreds, phred_length, phred_offset);
    if (error_rate < 0.0L) {
        return NULL;
    }
    double phred_score = -10 * log10(error_rate);
    return PyFloat_FromDouble(phred_score);
}


PyDoc_STRVAR(average_error_rate__doc__,
"average_error_rate($self, phred_scores, /, phred_offset=DEFAULT_PHRED_SCORE_OFFSET)\n"
"--\n"
"\n"
"Returns the average_error_rate. \n"
"\n"
"  phred_scores\n"
"    ASCII string with the phred scores.\n"
);

#define AVERAGE_ERROR_RATE_METHODDEF    \
    {"average_error_rate", (PyCFunction)(void(*)(void))average_error_rate_py, \
     METH_VARARGS | METH_KEYWORDS, average_error_rate__doc__}

static PyObject *
average_error_rate_py(PyObject *module, PyObject *args, PyObject *kwargs)
{
    PyObject *phred_scores = NULL;
    uint8_t phred_offset = DEFAULT_PHRED_SCORE_OFFSET;
    char *kwarg_names[] = {"", "phred_offset", NULL};
    const char *format = "O!|b:average_error_rate";
    if (!PyArg_ParseTupleAndKeywords(
        args, kwargs, format, kwarg_names,
        &PyUnicode_Type,
        &phred_scores,
        &phred_offset)) {
            return NULL;
    }

    if (!PyUnicode_IS_COMPACT_ASCII(phred_scores)) {
        PyErr_SetString(PyExc_ValueError,
                        "phred_scores must be ASCII encoded.");
        return NULL;
    }
    uint8_t *phreds = PyUnicode_DATA(phred_scores);
    size_t phred_length = PyUnicode_GET_LENGTH(phred_scores);
    double error_rate = average_error_rate(phreds, phred_length, phred_offset);
    if (error_rate < 0.0L) {
        return NULL;
    }
    return PyFloat_FromDouble(error_rate);
}

PyDoc_STRVAR(qualmedian__doc__,
"qualmedian($self, phred_scores, /, phred_offset=DEFAULT_PHRED_SCORE_OFFSET)\n"
"--\n"
"\n"
"Returns the median quality score. \n"
"\n"
"  phred_scores\n"
"    ASCII string with the phred scores.\n"
);

#define QUALMEDIAN_METHODDEF    \
    {"qualmedian", (PyCFunction)(void(*)(void))qualmedian_py, \
     METH_VARARGS | METH_KEYWORDS, qualmedian__doc__}

static PyObject *
qualmedian_py(PyObject *module, PyObject *args, PyObject *kwargs)
{
    PyObject *phred_scores = NULL;
    uint8_t phred_offset = DEFAULT_PHRED_SCORE_OFFSET;
    char *kwarg_names[] = {"", "phred_offset", NULL};
    const char *format = "O!|b:qualmedian";
    if (!PyArg_ParseTupleAndKeywords(
        args, kwargs, format, kwarg_names,
        &PyUnicode_Type,
        &phred_scores,
        &phred_offset)) {
            return NULL;
    }

    if (!PyUnicode_IS_COMPACT_ASCII(phred_scores)) {
        PyErr_SetString(PyExc_ValueError,
                        "phred_scores must be ASCII encoded.");
        return NULL;
    }
    uint8_t *phreds = PyUnicode_DATA(phred_scores);
    size_t phred_length = PyUnicode_GET_LENGTH(phred_scores);
    double median = qualmedian(phreds, phred_length, phred_offset);
    if (median < 0.0) {
        return NULL;
    }
    return PyFloat_FromDouble(median);
}

static PyMethodDef _filters_functions[] = {
    AVERAGE_ERROR_RATE_METHODDEF,
    QUALMEAN_METHODDEF,
    QUALMEDIAN_METHODDEF,
    {NULL}
};

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
    uint8_t phred_offset = DEFAULT_PHRED_SCORE_OFFSET;
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
    uint8_t phred_offset = DEFAULT_PHRED_SCORE_OFFSET;
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

static PyObject *
GenericFilter_ParseArgsToRecordTuple(PyObject *args, PyObject *kwargs) 
{
    if (kwargs != NULL) {
        PyErr_Format(PyExc_TypeError, 
                     "filter takes exactly 0 keyword arguments, got %d",
                     PyDict_GET_SIZE(kwargs));
        return NULL;
    }
    if (PyTuple_GET_SIZE(args) != 1) {
        PyErr_Format(PyExc_TypeError, 
                     "filter takes exactly 1 positional argument, got %d",
                     PyTuple_GET_SIZE(args));
        return NULL;
    }
    PyObject *arg = PyTuple_GET_ITEM(args, 0);
    if (!PyTuple_CheckExact(arg)) {
        PyErr_Format(PyExc_TypeError, 
                     "filter argument must be a tuple, got %s",
                     Py_TYPE(arg)->tp_name);
        return NULL;
    }
    Py_ssize_t record_tuple_length = PyTuple_GET_SIZE(arg);
    PyObject *record;
    for (Py_ssize_t i=0; i < record_tuple_length; i++) {
        record = PyTuple_GET_ITEM(arg, i);
        if (!(Py_TYPE(record) == SequenceRecord)) {
            PyErr_Format(
                PyExc_TypeError, 
                "All records must be of type dnaio.SequenceRecord, "
                "got %s at index %zd", 
                Py_TYPE(arg)->tp_name, i);
        return NULL;
        }
    }
    return arg;
}

static PyObject *
AverageErrorRateFilter__call__(FastqFilter *self, PyObject *args, PyObject *kwargs) 
{
    PyObject *record_tuple = GenericFilter_ParseArgsToRecordTuple(args, kwargs);
    if (record_tuple == NULL) {
        return NULL;
    }
    PyObject *record;
    PyObject *phred_scores;
    uint8_t *phreds;
    Py_ssize_t phred_length;
    uint8_t phred_offset = self->phred_offset;
    double total_error_sum = 0.0;
    size_t length_sum = 0;
    Py_ssize_t record_tuple_length = PyTuple_GET_SIZE(record_tuple);
    for (Py_ssize_t i=0; i < record_tuple_length; i++) {
        record = PyTuple_GET_ITEM(record_tuple, i);
        phred_scores = PyObject_GetAttr(record, QualitiesAttrString);
        if (phred_scores == NULL) {
            return NULL;
        }
        if (phred_scores == Py_None) {
            PyErr_Format(
                PyExc_ValueError,
                "SequenceRecord object with name %R does not have quality scores "
                "(FASTA record)", PyObject_GetAttrString(record, "name")
            );
            Py_DECREF(phred_scores);
            return NULL;
        }
        phreds = PyUnicode_DATA(phred_scores);
        phred_length = PyUnicode_GET_LENGTH(phred_scores);
        double error_sum = sum_error_rate(phreds, phred_length, phred_offset);
        Py_DECREF(phred_scores);
        if (error_sum < 0) {
            return NULL; 
        }
        total_error_sum += error_sum;
        length_sum += phred_length;
    }
    double error_rate = total_error_sum / (double)length_sum;
    self->total += 1;
    int pass = error_rate <= self->threshold_d;
    if (pass) {
        self->pass += 1;
    }
    return PyBool_FromLong(pass);
}

static PyObject *
MedianQualityFilter__call__(FastqFilter *self, PyObject *args, PyObject *kwargs) 
{
    PyObject *record_tuple = GenericFilter_ParseArgsToRecordTuple(args, kwargs);
    if (record_tuple == NULL) {
        return NULL;
    }
    Py_ssize_t record_tuple_length = PyTuple_GET_SIZE(record_tuple);
    uint8_t phred_offset = self->phred_offset;
    size_t total_phred_length = 0;
    int ret;
    PyObject *record;
    size_t histogram[128];
    memset(histogram, 0, sizeof(size_t) * 128);
    for (Py_ssize_t i=0; i < record_tuple_length; i++) {
        record = PyTuple_GET_ITEM(record_tuple, i);
        PyObject *phred_scores = PyObject_GetAttr(record, QualitiesAttrString);
        if (phred_scores == NULL) {
            return NULL;
        }
        if (phred_scores == Py_None) {
            PyErr_Format(
                PyExc_ValueError,
                "SequenceRecord object with name %R does not have quality scores "
                "(FASTA record)", PyObject_GetAttrString(record, "name")
            );
            Py_DECREF(phred_scores);
            return NULL;
        }
        uint8_t *phreds = PyUnicode_DATA(phred_scores);
        Py_ssize_t phred_length = PyUnicode_GetLength(phred_scores);
        ret = make_histogram(histogram, phreds, phred_length, phred_offset);
        Py_DECREF(phred_scores);
        if (ret != 0) {
            return NULL;
        }
        total_phred_length += phred_length;
    }
    double median = median_from_histogram(histogram, total_phred_length, phred_offset);
    if (median < 0.0) {
        return NULL;
    }
    self->total += 1;
    int pass = median >= self->threshold_d;
    if (pass) {
        self->pass += 1;
    }
    return PyBool_FromLong(pass);
}


static PyObject * 
MinLengthFilter__call__(FastqFilter *self, PyObject *args, PyObject *kwargs)
{
    PyObject *record_tuple = GenericFilter_ParseArgsToRecordTuple(args, kwargs);
    if (record_tuple == NULL) {
        return NULL;
    }
    PyObject *record;
    Py_ssize_t record_tuple_length = PyTuple_GET_SIZE(record_tuple);
    for (Py_ssize_t i=0; i < record_tuple_length; i++) {
        record = PyTuple_GET_ITEM(record_tuple, i);
        Py_ssize_t length = PyObject_Length(record);
        if (length < 0) {
            return NULL;
        }
        // If any of the records passes the minimum length we pass.
        // R1 and R2 sequence the same molecule so this is valid.
        if (length >= self->threshold_i) {
            self->pass += 1;
            self->total += 1;
            Py_RETURN_TRUE;
        }
    }
    self->total += 1;
    Py_RETURN_FALSE;
}

static PyObject *
MaxLengthFilter__call__(FastqFilter *self, PyObject *args, PyObject *kwargs) 
{
    PyObject *record_tuple = GenericFilter_ParseArgsToRecordTuple(args, kwargs);
    if (record_tuple == NULL) {
        return NULL;
    }
    PyObject *record;
    Py_ssize_t record_tuple_length = PyTuple_GET_SIZE(record_tuple);
    for (Py_ssize_t i=0; i < record_tuple_length; i++) {
        record = PyTuple_GET_ITEM(record_tuple, i);
        Py_ssize_t length = PyObject_Length(record);
        if (length < 0) {
            return NULL;
        }
        // If any of the records exceeds the maximum length we fail.
        // R1 and R2 sequence the same molecule so this is valid.
        if (length > self->threshold_i) {
            self->total += 1;
            Py_RETURN_FALSE;
        }
    }
    self->pass += 1;
    self->total += 1;
    Py_RETURN_TRUE;
}

static PyObject *
AverageErrorRateFilter_get_name(PyObject *self, void *closure)
{
    return PyUnicode_FromString("average error rate");
}

static PyObject *
MedianQualityFilter_get_name(PyObject *self, void *closure)
{
    return PyUnicode_FromString("median quality");
}

static PyObject *
MinimumLengthFilter_get_name(PyObject *self, void *closure)
{
    return PyUnicode_FromString("minimum length");
}

static PyObject *
MaximumLengthFilter_get_name(PyObject *self, void *closure)
{
    return PyUnicode_FromString("maximum length");
}

static PyGetSetDef AverageErrorRateFilter_properties[] = {
    {"name", AverageErrorRateFilter_get_name, NULL, NULL, NULL}, {NULL}};

static PyGetSetDef MedianQualityFilter_properties[] = {
    {"name", MedianQualityFilter_get_name, NULL, NULL, NULL}, {NULL}};

static PyGetSetDef MinimumLengthFilter_properties[] = {
    {"name", MinimumLengthFilter_get_name, NULL, NULL, NULL}, {NULL}};

static PyGetSetDef MaximumLengthFilter_properties[] = {
    {"name", MaximumLengthFilter_get_name, NULL, NULL, NULL}, {NULL}};

static PyTypeObject AverageErrorRateFilter_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_filter.AverageErrorRateFilter",
    .tp_basicsize = sizeof(FastqFilter),
    .tp_new = GenericQualityFilter__new__,
    .tp_call = (ternaryfunc)AverageErrorRateFilter__call__,
    .tp_members = GenericQualityFilterMembers,
    .tp_getset = AverageErrorRateFilter_properties,
};

static PyTypeObject MedianQualityFilter_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_filter.MedianQualityFilter",
    .tp_basicsize = sizeof(FastqFilter),
    .tp_new = GenericQualityFilter__new__,
    .tp_call = (ternaryfunc)MedianQualityFilter__call__,
    .tp_members = GenericQualityFilterMembers,
    .tp_getset = MedianQualityFilter_properties,
};

static PyTypeObject MinimumLengthFilter_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_filter.MinimumLengthFilter",
    .tp_basicsize = sizeof(FastqFilter),
    .tp_new = GenericLengthFilter__new__,
    .tp_call = (ternaryfunc)MinLengthFilter__call__,
    .tp_members = GenericLengthFilterMembers,
    .tp_getset = MinimumLengthFilter_properties
};

static PyTypeObject MaximumLengthFilter_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_filter.MaximumLengthFilter",
    .tp_basicsize = sizeof(FastqFilter),
    .tp_new = GenericLengthFilter__new__,
    .tp_call = (ternaryfunc)MaxLengthFilter__call__,
    .tp_members = GenericLengthFilterMembers,
    .tp_getset = MaximumLengthFilter_properties,
};

static struct PyModuleDef _filters_module = {
    PyModuleDef_HEAD_INIT,
    "_filters",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,
    _filters_functions  /* module methods */
};

#define MODULE_ADD_TYPE(module, typename, type) \
    PyTypeObject *typename = &type; \
        if (PyType_Ready(typename) == -1) { \
        return NULL; \
    } \
    PyModule_AddObject(module, #typename, \
                       (PyObject *)typename); 

PyMODINIT_FUNC
PyInit__filters(void)
{
    PyObject *m;

    m = PyModule_Create(&_filters_module);
    if (m == NULL) {
        return NULL;
    }
    PyObject *dnaio = PyImport_ImportModule("dnaio");
    if (dnaio == NULL) {
        return NULL;
    }
    SequenceRecord = (PyTypeObject *)PyObject_GetAttrString(dnaio, "SequenceRecord");
    if (SequenceRecord == NULL) {
        return NULL;
    }

    QualitiesAttrString = PyUnicode_FromString("qualities");
    if (QualitiesAttrString == NULL) {
        return NULL; 
    }

    MODULE_ADD_TYPE(m, AverageErrorRateFilter, AverageErrorRateFilter_Type)
    MODULE_ADD_TYPE(m, MedianQualityFilter, MedianQualityFilter_Type)
    MODULE_ADD_TYPE(m, MinimumLengthFilter, MinimumLengthFilter_Type)
    MODULE_ADD_TYPE(m, MaximumLengthFilter, MaximumLengthFilter_Type)

    PyModule_AddIntMacro(m, DEFAULT_PHRED_SCORE_OFFSET);
    return m;
}