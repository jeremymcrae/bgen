
import io
from pathlib import Path
from setuptools import setup
import sys

from distutils.core import Extension
from distutils.ccompiler import new_compiler
from Cython.Build import cythonize

EXTRA_COMPILE_ARGS = ['-std=c++11', '-I/usr/include']
EXTRA_LINK_ARGS = []
if sys.platform == "darwin":
    EXTRA_COMPILE_ARGS += ["-stdlib=libc++",
        "-mmacosx-version-min=10.9",
        "-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include",
        ]
    EXTRA_LINK_ARGS += ["-stdlib=libc++", "-mmacosx-version-min=10.9"]

def flatten(*lists):
    return [str(x) for sublist in lists for x in sublist]

def build_zstd():
    ''' compile zstd code to object files for linking with bgen c++ code
    
    This needs to be compiles independently of the bgen c++ code, as zstd is in
    c, so cannot be compiled directly with the bgen code. zstd is compiled to
    object files, and these are staticly linked in with the bgen code.
    
    I tried to link to a shared/dynamic zstd library, but couldn't specify a
    location that would be found by the compiled bgen library after relocating
    files into to a python package.
    
    Returns:
        list of paths to compiled object code
    '''
    folder = Path('src/zstd/lib')
    include_dirs = ['src/zstd/lib/', 'src/zstd/lib/common']
    sources = flatten(
        (folder / 'common').glob('*.c'),
        (folder / 'compress').glob('*.c'),
        (folder / 'decompress').glob('*.c'),
        (folder / 'dictBuilder').glob('*.c'),
        (folder / 'deprecated').glob('*.c'),
        (folder / 'legacy').glob('*.c'),  # TODO: drop some legacy versions
    )
    extra_compile_args = ['-std=gnu11', '-fPIC']
    
    cc = new_compiler()
    return cc.compile(sources, include_dirs=include_dirs,
        extra_preargs=extra_compile_args)

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
        extra_objects=build_zstd(),
        include_dirs=['src/', 'src/zstd/lib'],
        libraries=['z'],
        library_dirs=['bgen'],
        language='c++'),
    ])

setup(name='bgen',
    description='Package for loading data from bgen files',
    long_description=io.open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    version='1.2.0',
    author='Jeremy McRae',
    author_email='jmcrae@illumina.com',
    license="MIT",
    url='https://github.com/jeremymcrae/bgen',
    packages=['bgen'],
    install_requires=[
        'cython',
        'numpy',
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    ext_modules=reader,
    test_suite='tests',
    )
