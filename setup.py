from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize

extensions = [
    Extension(
        name="fractionalstd2inner.fractionals",
        sources=["src/fractionalstd2inner/fractionals.pyx"],  # relative path!
        language="c++",
        libraries=["gmp"],
        extra_compile_args=["-std=c++17"],
        include_dirs=["src/fractionalstd2inner"],  # relative path
    )
]

setup(
    name="fractionalstd2inner",
    version="1.0.5",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    ext_modules=cythonize(extensions, compiler_directives={"language_level": "3"}),
    include_package_data=True,
)
