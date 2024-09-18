=============
GRADE
=============

This reposatory includes all scripts nececarry to calculate the GRAIL-based Descriptor (GRADE) 
and the extended GRAIL-based Descriptor (X-GRADE) and the scripts that are nececarry to reproduce 
the results of the paper "GRADE and X-GRADE: Unveiling novel Protein-Ligand Interaction Fingerprints based on GRAIL-Scores"

GRADE/X-GRADE generation
-------------

To generate the GRAIL and X-GRAIL scores, the following steps are required:

1. Install the package using the following command:
```bash
pip install .
```

2. Run the following command to generate the GRAIL and X-GRAIL scores:
```bash
python GRADE/grade.py -i <input_file> -o <output_file> -t <type> -m <model>
```

where:
`<input_file>`     is the input file that contains the protein-ligand complexes in PDB format.
`<output_file>`    is the output file that contains the GRAIL and X-GRAIL scores.
`<type>`           is the type of the score that should be calculated. It can be either `GRADE` or `X-GRADE`.

Repoducing the results
-------------

To reproduce the results of the paper, one has to install a package called PHANTOMDRAGON 
(PHarmacophore bAsed scoriNg funcTiOn iMplementations using DRug interAction Grail scOre calculatioNs). 
For your convinience, the package is available right here in this repository. To install the package, 
the following steps are required:

1. Install the package using the following command:
```bash
pip install .
```

2. Use one of the scripts in the `scripts` directory to reproduce the results.

The scripts are (in alphabetical order):

* 3DQSAR_GRADE.ipynb
* create_core_set_data.py
* create_general_set_data.py
* create_refined_set_data.py
* extract_data.ipynb
* pred_classes.py
* pred_PL-REX.py
* pred_test_set.py
* pred_validation_set.py
* test_time.py
* UMAP_generation.ipynb



Credits
-------

This package was created with the help of Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
