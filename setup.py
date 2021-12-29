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

from pathlib import Path

from Cython.Build import cythonize

from setuptools import Extension, find_packages, setup

EXT_MODULES = cythonize(
    Extension("fastq_filter.optimized_algorithms",
              ["src/fastq_filter/optimized_algorithms.pyx"]),
    compiler_directives=dict(language_level="3", binding=True,
                             cdivision=True, profile=True))

LONG_DESCRIPTION = Path("README.rst").read_text()

setup(
    name="fastq-filter",
    version="0.1.1",
    description="A fast FASTQ filter progam.",
    author="Leiden University Medical Center",
    author_email="r.h.p.vorderman@lumc.nl",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/x-rst",
    license="MIT",
    keywords="FASTQ filter",
    zip_safe=False,
    packages=find_packages('src'),
    package_dir={'': 'src'},
    package_data={'fastq_filter': ['*.pyx', '*.pyi']},
    ext_modules=EXT_MODULES,
    url="https://github.com/lumc/fastq-filter",
    classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires=">=3.6",
    install_requires=[
        "xopen>=1.2.1",
        "dnaio>=0.6.0"
    ],
    entry_points={"console_scripts": [
        "fastq-filter = fastq_filter:main"]}
)
