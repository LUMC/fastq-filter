==========
Changelog
==========

.. Newest changes should be on top.

.. NOTE: This document is user facing. Please word the changes in such a way
.. that users understand how the changes affect the new version.

1.0.0-dev
--------------------
+ Refactored the _filters module code slightly to remove global variables.

0.3.0
--------------------
+ Improved performance by:

  + Slight rewrites in the multiple FASTQ pipeline. Cases with one, two, three
    and multiple FASTQ files are now handled separately. This vastly improves
    the speed of especially the single and paired end instances.
  + Rewriting the filters in C.
+ Added logging with basic stats. When ``--verbose`` is set, the counts for
  individual filters are reported. With ``--quiet`` logging can be turned off.

0.2.0
--------------------
+ Add ability for filtering paired sequences.
+ The minimum dnaio version requirement was bumped to 0.9.0.
+ Drop Python 3.6 support as Python 3.6 is end of life.
+ Filters are now applied by command flags. This is easier to understand than
  the previous method (which selected methods using strings). It is also
  easier to program and document.

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