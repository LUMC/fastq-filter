.. image:: https://img.shields.io/pypi/v/fastq-filter.svg
  :target: https://pypi.org/project/isal/
  :alt:

.. image:: https://img.shields.io/conda/v/bioconda/fastq-filter.svg
  :target: https://bioconda.github.io/recipes/fastq-filter/README.html
  :alt:

.. image:: https://img.shields.io/pypi/pyversions/fastq-filter.svg
  :target: https://pypi.org/project/isal/
  :alt:

.. image:: https://img.shields.io/pypi/l/fastq-filter.svg
  :target: https://github.com/LUMC/fastq-filter/blob/main/LICENSE
  :alt:

.. image:: https://codecov.io/gh/LUMC/fastq-filter/branch/main/graph/badge.svg?token=E85BEYDQ45
  :target: https://codecov.io/gh/LUMC/fastq-filter
  :alt:

=============
fastq-filter
=============

A fast FASTQ filter program.

Fastq-filter correctly takes into account that quality scores are log scores
when calculating the mean. It also provides an option to filter on average
error rate directly.

FASTQ Q=30 stands for an average error rate of 0.001, Q=20 for 0.01 and Q=10
for 0.1. This is not very intuitive. Q=20 has 10 times more errors than Q=30
though the numbers (20 and 30) do little to convey this difference. Using
0.01 and 0.001 correctly conveys that these error rates are an order of
magnitude apart. This also means that the phred scores cannot be naively
averaged. Q=10 and Q=30 do not average Q=20. The actual average error rate
is (0.001 + 0.1) / 2 = 0.0505. Roughly 1 in 20. Q=20 means 0.01: 1 in 100.
By naively averaging the quality is overestimated by a factor of 5! This
means any tool that averages naively is unusable in practice.

Unfortunately many tools do this. fastq-filter was written to provide a very
fast filtering solution so the correct filtering can be applied at a very low
cost.

Installation
============

+ With pip: ``pip install fastq-filter``
+ For the latest development version: ``pip install git+https://github.com/LUMC/fastq-filter``
+ With conda ``conda install -c conda-forge -c bioconda fastq-filter``


Usage
=====

Single fastq files can be filtered with::

    fastq-filter -e 0.001 -o output.fastq input.fastq

Multiple fastq files can be filtered with::

    fastq-filter -e 0.001 -o r1_filtered.fastq.gz -o r2_filtered.fastq.gz r1.fastq.gz r2.fastq.gz

Fastq-filter ensures the output is in sync. It is not limited to two inputs
so also ``R1.fq``, ``R2.fq`` and ``R3.fq`` can be filtered together.

In the following section 'pair' is used to note when 2 or more FASTQ records are
evaluated. When multiple FASTQ files are given the filters behave as follows:

+ average error rate: The average of the combined phred scores is used.
+ median quality: The median of the combined phred scores is used.
+ Minimum length: at least one of the records of the pair must meet the minimum length.
+ Maximum length: None of the records in the pair must exceed the maximum length.

The rationale for the length filters is that R1 and R2 both sequence the same
molecule and the canonical length is the longest of both.

.. code-block::

    usage: fastq-filter [-h] [-o OUTPUT] [-l MIN_LENGTH] [-L MAX_LENGTH]
                        [-e AVERAGE_ERROR_RATE] [-q MEAN_QUALITY]
                        [-Q MEDIAN_QUALITY] [-c COMPRESSION_LEVEL]
                        input [input ...]

    Filter FASTQ files on various metrics.

    positional arguments:
      input                 Input FASTQ files. Compression format automatically
                            detected. Use - for stdin.

    optional arguments:
      -h, --help            show this help message and exit
      -o OUTPUT, --output OUTPUT
                            Output FASTQ files. Compression format automatically
                            determined by file extension. Flag can be used
                            multiple times. An output must be given for each
                            input. Default: stdout.
      -l MIN_LENGTH, --min-length MIN_LENGTH
                            The minimum length for a read.
      -L MAX_LENGTH, --max-length MAX_LENGTH
                            The maximum length for a read.
      -e AVERAGE_ERROR_RATE, --average-error-rate AVERAGE_ERROR_RATE
                            The minimum average per base error rate.
      -q MEAN_QUALITY, --mean-quality MEAN_QUALITY
                            Average quality. Same as the '--average-error-rate'
                            option but specified with a phred score. I.e '-q 30'
                            is equivalent to '-e 0.001'.
      -Q MEDIAN_QUALITY, --median-quality MEDIAN_QUALITY
                            The minimum median phred score.
      -c COMPRESSION_LEVEL, --compression-level COMPRESSION_LEVEL
                            Compression level for the output files. Relevant when
                            output files have a .gz extension. Default: 2


Optimizations
=============

fastq-filter has used the following optimizations to be fast:

- Multiple filters can applied simultaneously to minimize IO.
- fastq-filter can be used in pipes to minimize IO
- The python filter function is used. Which is a a shorthand for python code
  that would otherwise need to be interpreted.
- The mean and median quality algorithms are implemented in C with bindings to
  Python.
- The mean quality algorithm uses a lookup table since there are only 93
  possible phred scores encoded in FASTQ. That saves a lot of power
  calculations to calculate the probabilities.
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
