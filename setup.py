
import glob
import io
from pathlib import Path
from setuptools import setup
import sys
import platform

from distutils.core import Extension
from distutils.ccompiler import new_compiler
from Cython.Build import cythonize

EXTRA_COMPILE_ARGS = ['-std=c++11', '-I/usr/include', '-O2']
EXTRA_LINK_ARGS = []
if sys.platform == "darwin":
    EXTRA_COMPILE_ARGS += [
        "-stdlib=libc++",
        "-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1",
        "-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include",
        ]
    EXTRA_LINK_ARGS += [
        "-L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib",
        ]
elif sys.platform == "win32":
    EXTRA_COMPILE_ARGS += ['/std:c++14', '/O2']

if platform.machine() == 'x86_64':
    EXTRA_COMPILE_ARGS += ['-mavx', '-mavx2']

def flatten(*lists):
    return [str(x) for sublist in lists for x in sublist]

def build_zlib():
    ''' compile zlib code to object files for linking with bgen on windows
    
    Returns:
        list of paths to compiled object code
    '''
    include_dirs = ['src/zlib/']
    sources = list(glob.glob('src/zlib/*.c'))
    extra_compile_args = ['/O2']
    
    cc = new_compiler()
    return cc.compile(sources, include_dirs=include_dirs,
        extra_preargs=extra_compile_args)

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
    print('building zstd')
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
    extra_compile_args = ['-std=gnu11', '-fPIC', '-O2']
    
    # newer zstd versions have an asm file, which needs to be to compiled and
    # used as a library for full zstd compilation.
    cc = new_compiler()
    # the unix compiler needs to allow files with .S extension
    cc.src_extensions += ['.S']
    compiled = cc.compile(sources=flatten((folder / 'decompress').glob('*.S')),)
    
    if len(compiled) > 0:
        extra_compile_args += [f'-L{" ".join(compiled)}']
    
    return cc.compile(sources, include_dirs=include_dirs,
        extra_preargs=extra_compile_args) + compiled

if sys.platform == 'win32':
    zlib, libs = build_zlib(), []
else:
    zlib, libs = [], ['z']
zstd = build_zstd()

ext = cythonize([
    Extension('bgen.reader',
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=['src/bgen/reader.pyx',
            'src/reader.cpp',
            'src/genotypes.cpp',
            'src/header.cpp',
            'src/samples.cpp',
            'src/utils.cpp',
            'src/variant.cpp'],
        extra_objects=zstd + zlib,
        include_dirs=['src', 'src/zstd/lib', 'src/zlib'],
        libraries=libs,
        language='c++'),
    Extension('bgen.writer',
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=['src/bgen/writer.pyx',
            'src/writer.cpp',
            'src/genotypes.cpp',
            'src/utils.cpp',
            ],
        extra_objects=zstd + zlib,
        include_dirs=['src', 'src/zstd/lib', 'src/zlib'],
        libraries=libs,
        language='c++'),
    ])

setup(name='bgen',
    description='Package for loading data from bgen files',
    long_description=io.open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    version='1.7.0',
    author='Jeremy McRae',
    author_email='jmcrae@illumina.com',
    license="MIT",
    url='https://github.com/jeremymcrae/bgen',
    packages=['bgen'],
    package_dir={'': 'src'},
    install_requires=[
        'numpy',
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    ext_modules=ext,
    test_loader='unittest:TestLoader',
    )
