from setuptools import setup
from glob import glob
import re
import io
import os


# Get version strip
def read(*names, **kwargs):
    with io.open(os.path.join(os.path.dirname(__file__), *names),
                 encoding=kwargs.get("encoding", "utf8")) as fp:
        return fp.read()

def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


setup(
    name="MutTui",
    version=find_version("MutTui/__init__.py"),
    author="Christopher Ruis",
    description=
    "A pipeline to reconstruct mutational spectra for bacterial and viral datasets",
    long_description_content_type="text/markdown",
    url="https://github.com/chrisruis/MutTui",
    install_requires=[
        'gffutils', 'biopython', 'phylo-treetime', 'PyQt5', 'matplotlib', 'pandas', 'sklearn'
    ],
    python_requires='>=3.8.0',
    packages=['MutTui'],
    keywords='transmission clustering metagenomics',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    entry_points={
        'console_scripts':
        ['MutTui = MutTui.__main__:main',],
    },
    extras_require={"test": "pytest"},
)



