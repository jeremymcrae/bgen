
import glob
import io
import os
from pathlib import Path
from setuptools import setup
import subprocess
import sys
import platform

from distutils.core import Extension
from distutils.ccompiler import new_compiler
from Cython.Build import cythonize

EXTRA_COMPILE_ARGS = []
EXTRA_LINK_ARGS = []
if sys.platform == 'linux':
    EXTRA_COMPILE_ARGS += ['-std=c++11', '-I/usr/include', '-O2']
elif sys.platform == "darwin":
    EXTRA_COMPILE_ARGS += [
        "-stdlib=libc++",
        "-std=c++11",
        "-O2",
        "-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1",
        "-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include",
        ]
    EXTRA_LINK_ARGS += [
        "-L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib",
        ]
elif sys.platform == "win32":
    EXTRA_COMPILE_ARGS += ['/std:c++14', '/O2']

if platform.machine() == 'x86_64' and sys.platform != "darwin":
    EXTRA_COMPILE_ARGS += ['-mavx', '-mavx2']

def flatten(*lists):
    return [str(x) for sublist in lists for x in sublist]

def build_zlib():
    ''' compile zlib code to object files for linking
    
    Returns:
        list of paths to compiled object code
    '''
    cur_dir = Path.cwd()
    source_dir = cur_dir / 'src' / 'zlib-ng'
    build_dir = cur_dir / 'zlib_build'
    build_dir.mkdir(exist_ok=True)
    os.chdir(build_dir)
    
    cmd = ['cmake', '-S', source_dir, '-B', build_dir,
        '-DZLIB_COMPAT=ON',
        '-DZLIB_ENABLE_TESTS=OFF',
        '-DBUILD_SHARED_LIBS=OFF',
        '-DCMAKE_C_FLAGS="-fPIC"',
    ]
    subprocess.run(cmd)
    subprocess.run(['cmake', '--build', build_dir, '-v', '--config', 'Release'])
    os.chdir(cur_dir)
    
    objs = [str(build_dir / 'libz.a')]
    if sys.platform == 'win32':
        objs = [str(build_dir / 'Release' / 'zlibstatic.lib'),
                ]
    
    return str(build_dir), objs

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

zlib_dir, zlib = build_zlib()
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
        include_dirs=['src', 'src/zstd/lib', zlib_dir],
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
        include_dirs=['src', 'src/zstd/lib', zlib_dir],
        language='c++'),
    ])

setup(name='bgen',
    description='Package for loading data from bgen files',
    long_description=io.open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    version='1.9.0',
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
