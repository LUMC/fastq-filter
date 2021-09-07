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
    "Naive setup with np.log10": textwrap.dedent(
        """
        import math
        import numpy as np
        
        def qualmean(qualities: bytes, phred_offset: int = 33) -> float:
            phred_scores = np.frombuffer(qualities, dtype=np.int8)
            probabilities = np.power(10, ((phred_scores - phred_offset) / -10))
            average = np.average(probabilities)
            return -10 * np.log10(average)
        """
    ),
    "Naive setup with math.log10": textwrap.dedent(
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
    "Use * -0.1 instead of division": textwrap.dedent(
        """
        import math
        import numpy as np
    
        def qualmean(qualities: bytes, phred_offset: int = 33) -> float:
            phred_scores = np.frombuffer(qualities, dtype=np.int8)
            probabilities = np.power(10, ((phred_scores - phred_offset) * -0.1))
            average = np.average(probabilities)
            return -10 * math.log10(average)
        """
    ),
    "Factor out the -0.1": textwrap.dedent(
        """
        import math
        import numpy as np
    
        def qualmean(qualities: bytes, phred_offset: int = 33) -> float:
            phred_scores = np.frombuffer(qualities, dtype=np.int8)
            probabilities = np.power((10 ** -0.1), ((phred_scores - phred_offset)))
            average = np.average(probabilities)
            return -10 * math.log10(average)
        """
    ),
    "Factor out the phred offset": textwrap.dedent(
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
    "Use constant": textwrap.dedent(
        """
        import math
        import numpy as np
        
        PHRED_CONSTANT = 10 ** -0.1
    
        def qualmean(qualities: bytes, phred_offset: int = 33) -> float:
            phred_scores = np.frombuffer(qualities, dtype=np.int8)
            probabilities = np.power(PHRED_CONSTANT, phred_scores)
            average = np.average(probabilities)
            return -10 * math.log10(average) - phred_offset
        """
    ),
    "Use np.mean": textwrap.dedent(
        """
        import math
        import numpy as np

        PHRED_CONSTANT = 10 ** -0.1

        def qualmean(qualities: bytes, phred_offset: int = 33) -> float:
            phred_scores = np.frombuffer(qualities, dtype=np.int8)
            probabilities = np.power(PHRED_CONSTANT, phred_scores)
            average = np.mean(probabilities)
            return -10 * math.log10(average) - phred_offset
        """
    ),
    "Fast python implementation": textwrap.dedent(
        """
        import math
        import array
    
        def qualmean(qualities: bytes, phred_offset: int = 33) -> float:
            phred_scores = array.array('b', qualities)
            phred_constant = 10 ** -0.1
            probabilities = [phred_constant ** score for score in phred_scores]
            average = sum(probabilities) / len(qualities)
            return -10 * math.log10(average) - phred_offset
        """
    ),
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
