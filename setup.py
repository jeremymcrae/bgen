
import glob
import io
import os
from pathlib import Path
from setuptools import setup, Extension
from setuptools.command.build_clib import build_clib
from setuptools.command.build_ext import build_ext
import subprocess
import sys

from Cython.Build import cythonize

ROOT = Path(__file__).resolve().parent
ZLIB_DIR = str(ROOT / 'zlib_build')

EXTRA_COMPILE_ARGS = []
EXTRA_LINK_ARGS = []
if sys.platform == 'linux':
    EXTRA_COMPILE_ARGS += ['-std=c++11', '-O2']
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

def flatten(*lists):
    return [str(x) for sublist in lists for x in sublist]

def build_zlib():
    ''' compile zlib code to object files for linking
    
    Returns:
        list of paths to compiled object code
    '''
    prev_dir = Path.cwd()
    source_dir = ROOT / 'src' / 'zlib-ng'
    build_dir = ROOT / 'zlib_build'
    build_dir.mkdir(exist_ok=True)
    os.chdir(build_dir)
    
    cmd = ['cmake', '-S', source_dir, '-B', build_dir,
        '-DZLIB_COMPAT=ON',
        '-DZLIB_ENABLE_TESTS=OFF',
        '-DBUILD_SHARED_LIBS=OFF',
        '-DCMAKE_C_FLAGS="-fPIC"',
    ]
    try:
        subprocess.run(cmd, check=True)
        subprocess.run(['cmake', '--build', build_dir, '-v', '--config', 'Release'],
                       check=True)
    finally:
        os.chdir(prev_dir)
    
    objs = [str(build_dir / 'libz.a')]
    if sys.platform == 'win32':
        objs = [str(build_dir / 'Release' / 'zlibstatic.lib'),
                ]
    
    return str(build_dir), objs

def zstd_sources():
    ''' list the zstd sources to compile into a static library
    
    zstd is in c, whereas the bgen code is c++, but setuptools selects the
    compiler per source file, so these can be built as a static library and
    linked in with the bgen code.
    
    Newer zstd versions also include an x86-64 assembly implementation of huffman
    decoding. That is compiled from a .S file, which needs the c preprocessor to
    resolve its includes, and which compiles to an empty object file on other
    architectures (it is wrapped in ZSTD_ENABLE_ASM_X86_64_BMI2 checks). MSVC
    cannot assemble .S files, but zstd only defines that macro for gnu-style
    compilers, so the c code does not expect the asm symbols there.
    
    Returns:
        list of paths to the zstd source files
    '''
    folder = Path('src/zstd/lib')
    sources = flatten(
        (folder / 'common').glob('*.c'),
        (folder / 'compress').glob('*.c'),
        (folder / 'decompress').glob('*.c'),
        (folder / 'dictBuilder').glob('*.c'),
        (folder / 'deprecated').glob('*.c'),
        (folder / 'legacy').glob('*.c'),  # TODO: drop some legacy versions
    )
    if sys.platform != 'win32':
        sources += flatten((folder / 'decompress').glob('*.S'))
    
    return sources

ZSTD_LIB = ('zstd', {
    'sources': zstd_sources(),
    'include_dirs': ['src/zstd/lib', 'src/zstd/lib/common'],
    })

class build_clib_subclass(build_clib):
    ''' allow zstd's .S assembly file to be compiled alongside its c code
    '''
    def build_libraries(self, libraries):
        # src_extensions is inherited from the compiler class, and += would add .S
        # to that shared list, so replace it with a copy that only this command
        # uses, to avoid affecting the compiler used for the extensions
        self.compiler.src_extensions = list(self.compiler.src_extensions) + ['.S']
        super().build_libraries(libraries)

class build_ext_subclass(build_ext):
    ''' compile the zlib dependency when the extensions are built
    
    zlib is only needed to link the extension modules, so compiling it here
    (rather than at import) keeps metadata-only commands such as egg_info and
    sdist from invoking cmake.
    '''
    def run(self):
        # build_ext links against the static libraries built by build_clib, but
        # does not run that command itself, so build zstd before the extensions
        self.run_command('build_clib')
        _, zlib = build_zlib()
        for ext in self.extensions:
            ext.extra_objects = zlib
        super().run()

extensions = [
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
        include_dirs=['src', 'src/zstd/lib', ZLIB_DIR],
        language='c++'),
    Extension('bgen.writer',
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=['src/bgen/writer.pyx',
            'src/writer.cpp',
            'src/genotypes.cpp',
            'src/utils.cpp',
            ],
        include_dirs=['src', 'src/zstd/lib', ZLIB_DIR],
        language='c++'),
    ]

setup(
    package_dir={'': 'src'},
    libraries=[ZSTD_LIB],
    ext_modules=cythonize(extensions),
    cmdclass={'build_clib': build_clib_subclass,
              'build_ext': build_ext_subclass},
    )
