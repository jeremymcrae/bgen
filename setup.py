
import io
from pathlib import Path
from setuptools import setup
import sys

from distutils.core import Extension
from Cython.Build import cythonize, build_ext

EXTRA_COMPILE_ARGS = ['-std=c++11']
EXTRA_LINK_ARGS = []
if sys.platform == "darwin":
    EXTRA_COMPILE_ARGS += ["-stdlib=libc++", "-mmacosx-version-min=10.9"]
    EXTRA_LINK_ARGS += ["-stdlib=libc++", "-mmacosx-version-min=10.9"]

def flatten(*lists):
    return [str(x) for sublist in lists for x in sublist]

folder = Path('src/zstd/lib')
zstd = Extension('bgen.libzstd',
    extra_compile_args=["-std=gnu11"],
    extra_link_args=['-shared',],
    sources=flatten(
        (folder / 'common').glob('*.c'),
        (folder / 'compress').glob('*.c'),
        (folder / 'dictBuilder').glob('*.c'),
        (folder / 'deprecated').glob('*.c'),
        (folder / 'legacy').glob('*.c'),  # TODO: drop some legacy versions
    ),
    include_dirs=['src/zstd/lib/',
        'src/zstd/lib/common',
    ],
    language='c'),
builder = build_ext('temp')
print(dir(zstd[0]))
builder.build_extension(zstd[0])


# ext_libraries = [['zstd', {
#     'extra_compile_args': ["-std=gnu11"],
#     # 'extra_link_args': ['-shared'],
#     'sources':  flatten(
#         (folder / 'common').glob('*.c'),
#         (folder / 'compress').glob('*.c'),
#         (folder / 'dictBuilder').glob('*.c'),
#         (folder / 'deprecated').glob('*.c'),
#         (folder / 'legacy').glob('*.c'),  # TODO: drop some legacy versions
#     ),
#     'include_dirs': ['src/zstd/lib/',
#         'src/zstd/lib/common',
#     ],
#     'language': 'c'
#     }
# ]]

reader = cythonize([
    # Extension('bgen.libzstd',
    #     extra_compile_args=["-std=gnu11"],
    #     extra_link_args=['-shared',],
    #     sources=flatten(
    #         (folder / 'common').glob('*.c'),
    #         (folder / 'compress').glob('*.c'),
    #         (folder / 'dictBuilder').glob('*.c'),
    #         (folder / 'deprecated').glob('*.c'),
    #         (folder / 'legacy').glob('*.c'),  # TODO: drop some legacy versions
    #     ),
    #     include_dirs=['src/zstd/lib/',
    #         'src/zstd/lib/common',
    #     ],
    #     language='c'),
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
        libraries=['z', 'zstd.cpython-37m-darwin'],
        # library_dirs=[str(x) for x in Path('build').glob('**') if x.is_dir() and 'lib' in str(x)],
        library_dirs=['bgen'],
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
    # ext_modules=[zstd()],
    ext_modules=reader,
    # libraries=ext_libraries,
    test_suite='tests',
    )
