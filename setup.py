
from pathlib import Path
from setuptools import setup, Extension
from setuptools.command.build_clib import build_clib
from setuptools.command.build_ext import build_ext
import subprocess
import sys

from Cython.Build import cythonize

ROOT = Path(__file__).resolve().parent

EXTRA_COMPILE_ARGS = []
# -ffp-contract=off avoids fused multiply-adds (e.g. on aarch64), so decoded floats are
# identical across architectures
if sys.platform == 'linux':
    EXTRA_COMPILE_ARGS += ['-std=c++11', '-O2', '-ffp-contract=off']
elif sys.platform == "darwin":
    EXTRA_COMPILE_ARGS += ["-stdlib=libc++", "-std=c++11", "-O2", '-ffp-contract=off']
elif sys.platform == "win32":
    EXTRA_COMPILE_ARGS += ['/std:c++14', '/O2']

def flatten(*lists):
    return [str(x) for sublist in lists for x in sublist]

def build_zlib(build_dir):
    ''' compile zlib-ng into a static library to link with
    
    Args:
        build_dir: folder to build in, which also gets the generated zlib.h
    
    Returns:
        list holding the path to the static library
    '''
    source_dir = ROOT / 'src' / 'zlib-ng'
    cmd = ['cmake', '-S', source_dir, '-B', build_dir,
        '-DZLIB_COMPAT=ON',
        '-DZLIB_ENABLE_TESTS=OFF',
        '-DBUILD_SHARED_LIBS=OFF',
        '-DCMAKE_POSITION_INDEPENDENT_CODE=ON',
    ]
    subprocess.run(cmd, check=True)
    subprocess.run(['cmake', '--build', build_dir, '-v', '--config', 'Release'],
                   check=True)
    
    if sys.platform == 'win32':
        return [str(build_dir / 'Release' / 'zlibstatic.lib')]
    return [str(build_dir / 'libz.a')]

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
    ''' link the vendored zlib and zstd into the extension modules
    
    zlib is only needed to link the extension modules, so compiling it here
    (rather than at import) keeps metadata-only commands such as egg_info and
    sdist from invoking cmake.
    
    '''
    def run(self):
        # build_ext links against the static libraries built by build_clib, but
        # does not run that command itself, so build zstd before the extensions
        self.run_command('build_clib')
        # hide shared libs during build to avoid linking to conda-based libzstd.so
        libraries = self.distribution.libraries
        self.distribution.libraries = []
        try:
            super().run()
        finally:
            self.distribution.libraries = libraries
    
    def build_extensions(self):
        # unlike run, this is only reached once self.compiler exists, so the
        # platform's static library naming can be asked for rather than guessed
        clib = self.get_finalized_command('build_clib')
        zstd = Path(clib.build_clib) / self.compiler.library_filename(
            'zstd', lib_type='static')
        if not zstd.exists():
            raise RuntimeError(f'the zstd static library is missing from {zstd}')
        # built alongside the other temporary files, rather than in the source tree
        zlib_dir = Path(self.build_temp) / 'zlib-ng'
        zlib = build_zlib(zlib_dir)
        for ext in self.extensions:
            ext.include_dirs.append(str(zlib_dir))
            ext.extra_objects = [str(zstd)] + zlib
        super().build_extensions()

extensions = [
    Extension('bgen.reader',
        extra_compile_args=EXTRA_COMPILE_ARGS,
        sources=['src/bgen/reader.pyx',
            'src/reader.cpp',
            'src/genotypes.cpp',
            'src/header.cpp',
            'src/samples.cpp',
            'src/streaming.cpp',
            'src/utils.cpp',
            'src/variant.cpp'],
        include_dirs=['src', 'src/zstd/lib'],
        language='c++'),
    Extension('bgen.writer',
        extra_compile_args=EXTRA_COMPILE_ARGS,
        sources=['src/bgen/writer.pyx',
            'src/writer.cpp',
            'src/genotypes.cpp',
            'src/utils.cpp',
            ],
        include_dirs=['src', 'src/zstd/lib'],
        language='c++'),
    ]

setup(
    package_dir={'': 'src'},
    libraries=[ZSTD_LIB],
    ext_modules=cythonize(extensions),
    cmdclass={'build_clib': build_clib_subclass,
              'build_ext': build_ext_subclass},
    )
