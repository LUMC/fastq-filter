import cProfile
import sys

from fastq_filter import filter_fastq

cProfile.run(f"filter_fastq('{sys.argv[1]}', '{sys.argv[2]}', '/dev/null')")