# -*- mode: python; tab-width: 4 -*- 

## 
# calc_descr_pdb_bind.py 
#
# Copyright (C) 2023 Thomas A. Seidel <thomas.seidel@univie.ac.at>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; see the file COPYING. If not, write to
# the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
##

import argparse
import os
import sys

import CDPL.Chem as Chem
import CDPL.Biomol as Biomol
import CDPL.Math as Math
import CDPL.GRAIL as GRAIL


LIG_ENV_MAX_RADIUS = 21.0
REMOVE_NON_STD_RESIDUES = True
NON_STD_RESIDUE_MAX_ATOM_COUNT = 8


def parseArguments():
    parser = argparse.ArgumentParser(description='Calculates GRAIL affinity prediction descriptors for a set of input ligand-protein complexes.')
    
    parser.add_argument('-d',
                        dest='complex_data_dir',
                        required=True,
                        help='[Required] The directory containing the ligand-protein complexes to process, organized in PL-REX manner.',
                        nargs=1)
    parser.add_argument('-o',
                        dest='out_csv_file',
                        required=True,
                        help='[Required] The path of the output CSV-file containing the descriptor values calculated for each input complex.',
                        nargs=1)
    parser.add_argument('-c',
                        dest='norm_chgs',
                        help='[Optional] Change protonation of acidic/basic groups to a state likely at pH7 (default: false)·',
                        action='store_true',
                        default=False)
    parser.add_argument('-x',
                        dest='ext_descr',
                        help='[Optional] Calculate extended GRAIL descriptor with subdivided HBA/HBD feature types (default: false)·',
                        action='store_true',
                        default=False)
    return parser.parse_args()

def removeNonStdResidues(pdb_code, lig_env):
    residues = Biomol.ResidueList(lig_env)
    
    for res in residues:
        is_std_res = Biomol.ResidueDictionary.isStdResidue(Biomol.getResidueCode(res))
        
        if is_std_res and res.numAtoms < 5:
            print('!! While processing complex %s: isolated standard residue fragment of size %s found' % (pdb_code, str(res.numAtoms)), file=sys.stderr)
            lig_env -= res
            
        elif REMOVE_NON_STD_RESIDUES and not is_std_res and res.numAtoms <= NON_STD_RESIDUE_MAX_ATOM_COUNT:
            if res.numAtoms == 1 and Chem.AtomDictionary.isMetal(Chem.getType(res.atoms[0])):
                continue

            lig_env -= res

def checkProtein(pdb_code, protein):
    for atom in protein.atoms:
        if Chem.getType(atom) == Chem.AtomType.H and atom.numAtoms == 0:
            print('!! While processing complex %s: isolated hydrogen atom encountered' % pdb_code, file=sys.stderr)
             
        elif Chem.getType(atom) == Chem.AtomType.UNKNOWN:
            print('!! While processing complex %s: atom of unknown element encountered' % pdb_code, file=sys.stderr)

def processComplex(pdb_code, comp_data_dir, out_file, descr_calc, norm_chgs):
    print('Processing complex %s...' % pdb_code)

    sdf_reader = Chem.FileSDFMoleculeReader(comp_data_dir + '/' + 'ligand.sdf')
    ligand = Chem.BasicMolecule()

    if not sdf_reader.read(ligand):
        print('!! While processing complex %s: reading ligand SD-file failed' % pdb_code, file=sys.stderr)
        return
        
    pdb_reader = Biomol.FilePDBMoleculeReader(comp_data_dir + '/' + 'protein.pdb')
    protein = Chem.BasicMolecule()

    if not pdb_reader.read(protein):
        print('!! While processing complex %s: reading protein PDB-file failed' % pdb_code, file=sys.stderr)
        return
    
    checkProtein(pdb_code, protein)

    GRAIL.prepareForGRAILDescriptorCalculation(ligand, norm_chgs)
    GRAIL.prepareForGRAILDescriptorCalculation(protein, norm_chgs)

    lig_env = Chem.Fragment()

    Biomol.extractEnvironmentResidues(ligand, protein, lig_env, Chem.Atom3DCoordinatesFunctor(), LIG_ENV_MAX_RADIUS, False)
    removeNonStdResidues(pdb_code, lig_env)
    Chem.extractSSSRSubset(protein, lig_env, True)
        
    line = pdb_code
    descr = Math.DVector()
    lig_atom_coords = Math.Vector3DArray()

    Chem.get3DCoordinates(ligand, lig_atom_coords)

    descr_calc.initTargetData(lig_env, Chem.Atom3DCoordinatesFunctor())
    descr_calc.initLigandData(ligand)

    descr_calc.calculate(lig_atom_coords, descr)

    for i in range(0, descr_calc.TOTAL_DESCRIPTOR_SIZE):
        line += (', ' + str(descr(i)))

    out_file.write(line + '\n')
    out_file.flush()

def outputColNames(out_file, descr_calc):
    out_file.write('PDB code')

    for cn in descr_calc.ElementIndex.names.keys():
        out_file.write(', ' + cn)

    out_file.write('\n')
    out_file.flush()
    
def process(args):
    out_file = open(args.out_csv_file[0], 'w')

    if args.ext_descr:
        descr_calc = GRAIL.GRAILXDescriptorCalculator()
    else:
        descr_calc = GRAIL.GRAILDescriptorCalculator()

    outputColNames(out_file, descr_calc)
    
    for pdb_code in os.listdir(args.complex_data_dir[0]):
        comp_data_dir = os.path.join(args.complex_data_dir[0], pdb_code)

        if os.path.isfile(comp_data_dir): # sanity check
            continue

        try:
            df = processComplex(pdb_code, comp_data_dir, out_file, descr_calc, args.norm_chgs)

        except Exception as e:
            print('!! Processing complex %s failed: ' % pdb_code, e, file=sys.stderr)

    out_file.close()

    print('Done!')
    
if __name__ == '__main__':
    process(parseArguments())
