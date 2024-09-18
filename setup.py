#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [ ]

test_requirements = [ ]

setup(
    author="Christian Fellinger",
    author_email='christian.fellinger@univie.ac.at',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="PHarmacophore bAsed scoriNg funcTiOn iMplementations using DRug interAction Grail scOre calculatioNs",
    entry_points={
        'console_scripts': [
            'phantomdragon=phantomdragon.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='phantomdragon',
    name='phantomdragon',
    packages=find_packages(include=['phantomdragon', 'phantomdragon.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/dragon3221/phantomdragon',
    version='0.1.0',
    zip_safe=False,
)
