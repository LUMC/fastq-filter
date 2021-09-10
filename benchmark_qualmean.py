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

import textwrap
QUALMEAN_CODE = {
    "Numpy: Naive setup with math.log10": textwrap.dedent(
        """
        import math
        import numpy as np
    
        def qualmean(qualities: bytes, phred_offset: int = 33) -> float:
            phred_scores = np.frombuffer(qualities, dtype=np.int8)
            probabilities = np.power(10, ((phred_scores - phred_offset) / -10))
            average = np.average(probabilities)
            return -10 * math.log10(average)
        """
    ),
    "Numpy: Factor out the phred offset and the /-10": textwrap.dedent(
        """
        import math
        import numpy as np
    
        def qualmean(qualities: bytes, phred_offset: int = 33) -> float:
            phred_scores = np.frombuffer(qualities, dtype=np.int8)
            probabilities = np.power((10 ** -0.1), phred_scores)
            average = np.average(probabilities)
            return -10 * math.log10(average) - phred_offset
        """
    ),
    "Numpy: Factored and float32": textwrap.dedent(
        """
        import math
        import numpy as np
    
        def qualmean(qualities: bytes, phred_offset: int = 33) -> float:
            phred_scores = np.frombuffer(qualities, dtype=np.int8)
            probabilities = np.power((10 ** -0.1), phred_scores, dtype=np.float32)
            average = np.average(probabilities)
            return -10 * math.log10(average) - phred_offset
        """
    ),
    "Fast pure python implementation": textwrap.dedent(
        """
        import math
        import array
    
        def qualmean(qualities: bytes, phred_offset: int = 33) -> float:
            phred_constant = 10 ** -0.1
            sum_probabilities = 0.0
            for score in qualities:
                sum_probabilities += phred_constant ** score
            average = sum_probabilities / len(qualities)
            return -10 * math.log10(average) - phred_offset  
        """
    ),
    "Cython implementation": textwrap.dedent(
        """
        from fastq_filter.optimized_algorithms import qualmean
        """
    ),
}
QUALMEDIAN_CODE = {
    "Median numpy implementation with float": textwrap.dedent(
        """
        import numpy as np
    
        def qualmedian(qualities: bytes, phred_offset: int = 33) -> float:
            phred_scores = np.frombuffer(qualities, dtype=np.int8)
            return float(np.median(phred_scores)) - phred_offset
        """),
    "Median simple python implementation": textwrap.dedent(
        """
        import statistics
        import array
    
        def qualmedian(qualities: bytes, phred_offset: int = 33) -> float:
            return statistics.median(qualities) - phred_offset
        """),
    "Adaptive implementation": textwrap.dedent(
        """
        import statistics
        import array
        import numpy as np
    
        def qualmedian(qualities: bytes, phred_offset: int = 33) -> float:
            if len(qualities) < 500:
                return statistics.median(qualities) - phred_offset
            else:
                phred_scores = np.frombuffer(qualities, dtype=np.int8)
                return float(np.median(phred_scores)) - phred_offset
        """),
    "Cython implementation": textwrap.dedent(
        """
        from fastq_filter.optimized_algorithms import qualmedian
        """),
}

import timeit
import statistics

def benchmark(name, line, setup, loops=10000, runs=10):
    print(f"{name}")
    results = [timeit.timeit(line, setup, number=loops) for _ in range(runs)]
    # Calculate microseconds
    results = [(result / loops) * 1_000_000 for result in results]
    print(f"average: {round(statistics.mean(results), 2)}, "
          f"range: {round(min(results), 2)}-{round(max(results),2)} "
          f"stdev: {round(statistics.stdev(results),2)}")


if __name__ == "__main__":
    for name, setup in QUALMEAN_CODE.items():

        benchmark(name,
                  "qualmean(b'GFFAFFFFDEDDEFFFC?FFFFF8FF9FDFEFCDFFFGEDDFCFEFFD"
                  "FFF/FEED?FFF7FFF;FFFEDFEFC8>FCF>EDDEBDA2F@;FFCFFDDEFF@E9FFD"
                  "FFEFFFFFFFFFFFFEBFF?A1GCFFFF=FCFFDGFECEEEDG')",
                  setup)

    for name, setup in QUALMEDIAN_CODE.items():

        benchmark(name,
                  "qualmedian(b'GFFAFFFFDEDDEFFFC?FFFFF8FF9FDFEFCDFFFGEDDFCFEF"
                  "FDFFF/FEED?FFF7FFF;FFFEDFEFC8>FCF>EDDEBDA2F@;FFCFFDDEFF@E9F"
                  "FDFFEFFFFFFFFFFFFEBFF?A1GCFFFF=FCFFDGFECEEEDG')",
                  setup)

