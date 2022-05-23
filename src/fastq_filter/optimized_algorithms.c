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

#include "score_to_error_rate.h"
#define MAXIMUM_PHRED_SCORE 126
#define DEFAULT_PHRED_OFFSET 33

PyDoc_STRVAR(qualmean__doc__,
"qualmean($self, phred_scores, /, phred_offset=DEFAULT_PHRED_OFFSET)\n"
"--\n"
"\n"
"Returns the average error rate as a float. \n"
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
    uint8_t phred_offset = DEFAULT_PHRED_OFFSET;
    char *kwarg_names[] = {"", "phred_offset", NULL};
    const char *format = "O!|$b:qualmean";
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
    double total_error_rate = 0.0;
    uint8_t *scores = PyUnicode_DATA(phred_scores);
    uint8_t score;
    uint8_t max_score = MAXIMUM_PHRED_SCORE - phred_offset;
    Py_ssize_t length = PyUnicode_GET_LENGTH(phred_scores);
    for (Py_ssize_t i=0; i<length; i+=1) {
        score = scores[i] - phred_offset;
        if (score > max_score) {
            PyErr_Format(
                PyExc_ValueError,
                "Character %c outside of valid phred range ('%c' to '%c')",
                scores[i], phred_offset, MAXIMUM_PHRED_SCORE);
            return NULL;
        }
        total_error_rate += SCORE_TO_ERROR_RATE[score];
    }
    double average_error = total_error_rate / (double)length;
    return PyFloat_FromDouble(average_error);
}

static PyMethodDef optimized_algorithms_functions[] = {
    QUALMEAN_METHODDEF,
    {NULL}
};

static struct PyModuleDef optimized_algorithms_module = {
    PyModuleDef_HEAD_INIT,
    "optimized_algorithms",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,
    optimized_algorithms_functions  /* module methods */
};

PyMODINIT_FUNC
PyInit_optimized_algorithms(void)
{
    PyObject *m;

    m = PyModule_Create(&optimized_algorithms_module);
    if (m == NULL) {
        return NULL;
    }
    PyModule_AddIntMacro(m, DEFAULT_PHRED_OFFSET);
    return m;
}