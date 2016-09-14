#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    
#    'numpy',
#    'matplotlib',
#    'pandas',
#    'scipy',
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='minke',
    version='1.0.0',
    description="Minke is a Python package to produce Mock Data Challenge data sets for LIGO interferometers.",
    long_description=readme + '\n\n' + history,
    author="Daniel Williams",
    author_email='d.williams.2@research.gla.ac.uk',
    url='https://github.com/transientlunatic/minke',
    packages=[
        'minke',
    ],
    package_dir={'minke':
                 'minke'},
    include_package_data=True,
    install_requires=requirements,
    license="ISCL",
    zip_safe=False,
    keywords='minke',
    classifiers=[
        'Development Status :: 4 - Beta',
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
