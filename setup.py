
import io
from setuptools import setup
import sys
from distutils.core import Extension
from Cython.Build import cythonize

EXTRA_COMPILE_ARGS = ['-std=c++11']
EXTRA_LINK_ARGS = []
if sys.platform == "darwin":
    EXTRA_COMPILE_ARGS += ["-stdlib=libc++", "-mmacosx-version-min=10.9"]
    EXTRA_LINK_ARGS += ["-stdlib=libc++", "-mmacosx-version-min=10.9"]

reader = cythonize([
    Extension('bgen.reader',
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=['bgen/bgen.pyx',
            'src/bgen.cpp',
            'src/genotypes.cpp',
            'src/header.cpp',
            'src/samples.cpp',
            'src/utils.cpp',
            'src/variant.cpp'],
        include_dirs=['src/', 'src/zstd/lib'],
        libraries=['z'],
        language='c++'),
    ])

setup(name='bgen',
    description='Package for loading data from bgen files',
    long_description=io.open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    version='1.0.0',
    author='Jeremy McRae',
    author_email='jmcrae@illumina.com',
    license="MIT",
    url='https://github.com/jeremymcrae/bgen',
    packages=['bgen'],
    install_requires=['cython',
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    ext_modules=reader,
    test_suite='tests',
    )
