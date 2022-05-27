==========
Changelog
==========

.. Newest changes should be on top.

.. NOTE: This document is user facing. Please word the changes in such a way
.. that users understand how the changes affect the new version.

0.2.0-dev
--------------------
+ Drop Python 3.6 support as Python 3.6 is end of life.
+ Filters are now applied by command flags. This is easier to understand than
  the previous method (which selected methods using strings). It is also
  easier to program and document.
+ Improve performance by refactoring the filtering pipeline. Filters are now
  instantiated as a Python C class. Any thresholds
  and phred offsets are converted from Python variables to C variables only
  once during filter instantiation, instead of every time a filter is called.
+ Improve performance by not double checking guarantees given by dnaio. The
  minimum dnaio version requirement was bumped to 0.8.0.

0.1.1
--------------------
+ Various documentation fixes
+ Include a galaxy wrapper

0.1.0
--------------------
+ Use a lookup table for faster quality lookup
+ Add checking to ensure all quality values are correct.
+ Create a fastq-filter program with functions for length, mean quality and
  median quality filtering.