#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ["scikit-learn","xgboost","numpy==1.21.5","pandas","matplotlib","seaborn","scipy", "umap-learn", "rdkit","joblib","CDPKit","ipykernel","umap-learn[plot]"]

test_requirements = [ ]

setup(
    author="Christian Fellinger, Thomas Seidel",
    author_email='christian.fellinger@univie.ac.at, thomas.seidel@univie.ac.at',
    python_requires='>3.7,<3.11',
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
    description="This reposatory includes all scripts nececarry to calculate the GRAIL-based Descriptor (GRADE) and the extended GRAIL-based Descriptor (X-GRADE) and the scripts that are nececarry to reproduce the results of the paper 'GRADE and X-GRADE: Unveiling novel Protein-Ligand Interaction Fingerprints based on GRAIL-Scores'. This includes the testing suite PHANTOMDRAGON (PHarmacophore bAsed scoriNg funcTiOn iMplementations using DRug interAction Grail scOre calculatioNs).",
    entry_points={
        'console_scripts': [
            'phantomdragon=phantomdragon.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='GRADE, PHANTOMDRAGON, GRAIL, X-GRADE, GRAIL-Scores, Protein-Ligand Interaction Fingerprints',
    name='PHANTOMDRAGON',
    packages=find_packages(include=['phantomdragon', 'phantomdragon.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/molinfo-vienna/GRADE',
    version='1.0.0',
    zip_safe=False,
)
