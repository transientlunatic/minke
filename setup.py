#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools.command.build_ext import build_ext as _build_ext

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

# see https://stackoverflow.com/a/21621689/1862861 for why this is here
class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())

setup_requirements = [
    'numpy',
    'setuptools_scm'
]

with open("requirements.txt") as requires_file:
    requirements = requires_file.read().split("\n")

test_requirements = [
    "py",
    "pytest",
    "coverage"
]

setup(
    name='minke',
    #version='1.0.1',
    use_scm_version=True,
    description="Minke is a Python package to produce Mock Data Challenge data sets for LIGO interferometers.",
    long_description=readme + '\n\n' + history,
    author="Daniel Williams",
    author_email='daniel.williams@ligo.org',
    url='https://github.com/transientlunatic/minke',
    packages=[
        'minke',
    ],
    package_dir={'minke':
                 'minke'},
    include_package_data=True,
    setup_requires = setup_requirements,
    install_requires=requirements,
    license="ISCL",
    zip_safe=False,
    keywords='minke',
    classifiers=[
        'Development Status :: 5 - Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: ISC License (ISCL)',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
