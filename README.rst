=============
fastq-filter
=============

A fast FASTQ filter program.

Fastq-filter correctly takes into account that quality scores are log scores
when calculating the mean.

Installation
============

For the latest development version

.. code-block::

    pip install git+https://github.com/LUMC/fastq-filter


Quickstart
==========
.. code-block::

    fastq-filter mean_quality:20 my.fastq

This will filter out all fastq files that have a mean quality below 20.

Other filters are ``median_quality``, ``min_length`` and ``max_length``.
For more information use: ``fastq-filter --help-filters`` or see the filters
chapter below.

Fastq-filter can also chain filters together:

.. code-block::

    fastq-filter 'min_length:100|mean_quality:20' my.fastq

It is advisible to put the fastest filters (length) before the slower ones
(quality) to optimize performance.

Usage
=====

.. code-block::

    usage: fastq-filter [-h] [--help-filters] [-o OUTPUT] filters input

    positional arguments:
      filters               Filters and arguments. For example: mean_quality:20,
                            for filtering all reads with an average quality below
                            20. Multiple filters can be applied by separating with
                            the | symbol. For example:
                            min_length:100|mean_quality:20. Make sure to use
                            faster filters (length) before slower ones (quality)
                            for optimal performance. Use --help-filters to print
                            all the available filters.
      input                 Input FASTQ file. Compression format automatically
                            detected.

    optional arguments:
      -h, --help            show this help message and exit
      --help-filters        Print all the available filters.
      -o OUTPUT, --output OUTPUT
                            Output FASTQ file. Compression format automatically
                            determined by file extension. Default: stdout.

Filters
=======

============================== ===================================================================================
mean_quality:<quality>         The mean quality of the FASTQ record is equal or above the given quality value.
median_quality:<quality>       The median quality of the FASTQ record is equal or above the given quality value.
min_length:<length>            The length of the sequence in the FASTQ record is at least min_length
max_length:<length>            The length of the sequence in the FASTQ record is at most max_length
============================== ===================================================================================

Optimizations
=============

fastq-filter has used the following optimizations to be fast:

- Filters can be chained together to minimize IO.
- The python filter function is used. Which is a a shorthand for python code
  that would otherwise need to be interpreted.
- The mean and median quality algorithms are implemented in Cython.
- The mean quality algorithm has been optimized to use the minimum number of
  calculations.
- The mean quality algorithm calculates the per-base quality from the phred
  score as a single-precision (32-bit) floating point number instead of a
  double-precision (64-bit) floating point number. This is 4 times faster.
  Given that phred scores have at most 2 significant digits, using
  double-precision does not add any meaningful accuracy.
- The median quality algorithm implements a counting sort, which is really
  fast but not applicable for generic data. Since FASTQ qualities are uniquely
  suited for a counting sort, median calculation can be performed very quickly.
- `dnaio <https://github.com/marcelm/dnaio>`_ is used as FASTQ parser.  This
  parses the FASTQ files with a parser written in Cython.
- `xopen <https://github.com/pycompression/xopen>`_ is used to read and write
  files. This allows for support of gzip compressed files which are opened
  using `python-isal <https://github.com/pycompression/python-isal>`_ which
  reads gzip files 2 times faster and writes gzip files 5 times faster than
  the python ``gzip`` module implementation.
