
# GRADE

This reposatory includes all scripts necessary to calculate the GRAil-based DEscriptor (GRADE) and the eXtended GRAil-based DEscriptor (X-GRADE) and the scripts that are necessary to reproduce the results of the paper "GRADE and X-GRADE: Unveiling novel Protein-Ligand Interaction Fingerprints based on GRAIL-Scores"

## Install

To reproduce the results of the paper, one has to install a package called PHANTOMDRAGON 
(PHarmacophore bAsed scoriNg funcTiOn iMplementations using DRug interAction Grail scOre calculatioNs). For your convinience, the package is available right here in this repository. To install the package, we recommend using a conda or mamba environment:

``conda create -n GRADE python==3.10``

This should give you a blank environment with only python and its environments installed. For this package python ``>3.7`` and ``<3.11`` is reqired.
Then activate the environment with:

``conda activate GRADE``

Then install the package with:

``pip install .``

Of course one has to be in the GRADE directory for this to work. This should also install all dependencies that you need.

## GRADE/X-GRADE generation

Depending on your data structure, run one of the following commands to generate GRADE or X-GRADE:

``calc_descr_pdb_bind.py [-h] -d COMPLEX_DATA_DIR -o OUT_CSV_FILE [-c] [-x]``

Calculates GRADE/X-GRADE for a set of input ligand-protein complexes. The Files have to be organized in PDBbind manner.

| Option               | Description                                                                                                    | Required | Default     |
|----------------------|----------------------------------------------------------------------------------------------------------------|----------|-------------|
| `-h` or `--help`       | Show this help message and exit                                                                               | No       | N/A         |
| `-d COMPLEX_DATA_DIR`| The directory containing the ligand-protein complexes to process, organized in PDBBind manner                   | Yes      | N/A         |
| `-o OUT_CSV_FILE`    | The path of the output CSV-file containing the descriptor values calculated for each input complex               | Yes      | N/A         |
| `-c`                 | Change protonation of acidic/basic groups to a state likely at pH7                                              | No       | false       |
| `-x`                 | Calculate extended GRAIL descriptor with subdivided HBA/HBD feature types                                       | No       | false       |

``calc_descr_PL_REX.py [-h] -d COMPLEX_DATA_DIR -o OUT_CSV_FILE [-c] [-x]``

Calculates GRADE/X-GRADE for a set of input ligand-protein complexes. The Files have to be organized in PL-REX manner.

| Option               | Description                                                                                                    | Required | Default     |
|----------------------|----------------------------------------------------------------------------------------------------------------|----------|-------------|
| `-h` or `--help`       | Show this help message and exit                                                                               | No       | N/A         |
| `-d COMPLEX_DATA_DIR`| The directory containing the ligand-protein complexes to process, organized in PL-REX manner                   | Yes      | N/A         |
| `-o OUT_CSV_FILE`    | The path of the output CSV-file containing the descriptor values calculated for each input complex               | Yes      | N/A         |
| `-c`                 | Change protonation of acidic/basic groups to a state likely at pH7                                              | No       | false       |
| `-x`                 | Calculate extended GRAIL descriptor with subdivided HBA/HBD feature types                                       | No       | false       |

``calc_descr_pdb_ligands.py [-h] -p PDB_FILE -l LIG_FILE -o OUT_CSV_FILE [-x] [-c]``

Calculates GRADE/X-GRADE for a PDB-file and set of input ligands.

| Option               | Description                                                                                                    | Required | Default     |
|----------------------|----------------------------------------------------------------------------------------------------------------|----------|-------------|
| `-h` or `--help`       | Show this help message and exit                                                                               | No       | N/A         |
| `-p PDB_FILE`        | The receptor PDB-file                                                                                           | Yes      | N/A         |
| `-l LIG_FILE`        | The file providing the ligands                                                                                  | Yes      | N/A         |
| `-o OUT_CSV_FILE`    | The path of the output CSV-file containing the descriptors calculated for each input ligand                     | Yes      | N/A         |
| `-x`                 | Calculate extended GRAIL descriptor with subdivided HBA/HBD feature types                                       | No       | false       |
| `-c`                 | Change protonation of acidic/basic groups to a state likely at pH7                                              | No       | false       |



## Repoducing the results

You can use one of the scripts in the `scripts` directory to reproduce the results.

The scripts are (in alphabetical order):

| Script Name                | Description                                                                                           | Type               |
|----------------------------|-------------------------------------------------------------------------------------------------------|--------------------|
| `3DQSAR_GRADE.ipynb`       | Performs 3D QSAR analysis using parts of GRADE and X-GRADE.    | Jupyter Notebook   |
| `calc_descr_pdb_bind.py`   | Calculates GRADE/X-GRADE for a set of input ligand-protein complexes. The Files have to be organized in PDBbind manner. (see above) | Python Script       |
| `calc_descr_pdb_ligands.py`| Calculates GRADE/X-GRADE for a PDB-file and set of input ligands. (see above) | Python Script       |
| `calc_descr_PL-REX.py`     | Calculates GRADE/X-GRADE for a set of input ligand-protein complexes. The Files have to be organized in PL-REX manner. (see above) | Python Script       |
| `create_core_set_data.py`  | Creates the PDBbind core dataset for model evaluation.    | Python Script       |
| `create_general_set_data.py`| Creates the PDBbind general set for model testing. | Python Script       |
| `create_refined_set_data.py`| Creates the PDBbind refined set for modle training. | Python Script       |
| `extract_data.ipynb`       | Extracts data from various sources and formats it for analysis.    | Jupyter Notebook   |
| `pred_classes.py`          | Estimates binding affinities for different EC numbers. | Python Script       |
| `pred_PL-REX.py`           | Estimates binding affinities for the PL-REX Dataset. | Python Script       |
| `pred_test_set.py`         | Evaluates estimations on a test set. | Python Script       |
| `pred_validation_set.py`   | Evaluates estimations on a validation set.     | Python Script       |
| `test_time.py`             | Measures and tests the runtime of the binding affinity estimations. | Python Script       |
| `UMAP_generation.ipynb`    | Generates UMAP visualizations for data analysis. Also includes t-SNE. | Jupyter Notebook   |


## Credits

CDPKit and the GRADE and X-GRADE implementation was done by Thomas Seidel (https://github.com/seidelt)

PHANTOMDRAGON, all testing and this Reposatory was created by Christian  Fellinger (https://github.com/Dragon3221)

This package was created with the help of Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
