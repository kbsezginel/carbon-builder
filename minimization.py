# import os
import ase.io
# import subprocess
import numpy as np
# from rdkit import Chem
# from rdkit.Chem import AllChem
from openbabel import openbabel
from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.geometry.analysis import Analysis
# from rdkit2ase import rdkit2ase, ase2rdkit


FORCE_FIELDS = ['mmff94', 'mmff94s', 'ghemical', 'gaff', 'uff', 'EMT']


def minimize(pdb_file, force_field, steps=20, st=None):
    """
    Minimize pdb file with the selected force field using openbabel / ASE / rdKit.
    """
    if force_field in ['mmff94', 'mmff94s', 'ghemical', 'gaff', 'uff']:
        atoms = obff_minimize(pdb_file, steps=steps, ff=force_field, st=st)
        # atoms = obminimize(pdb_file, steps=steps, ff=force_field, st=st)
    elif force_field == 'EMT':
        atoms = emt_minimize(pdb_file)
    # elif force_field == 'MMFF (rdKit)':
    #     atoms = mmff_minimize(pdb_file)
    else:
        raise Exception(f'Force field {force_field} not available')
    return atoms


def emt_minimize(pdb_file):
    """ASE EMT minimization"""
    atoms_new = ase.io.read(pdb_file)
    atoms_new.calc = EMT()
    dyn = BFGS(atoms_new)
    dyn.run(fmax=0.5)
    return atoms_new


# def mmff_minimize(pdb_file):
#     """rdKit mmff minimization"""
#     mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)  # Set 
#     # Chem.AddHs(mol)
#     # mol = ase2rdkit(atoms)
#     mp = AllChem.MMFFGetMoleculeProperties(mol)
#     ff = AllChem.MMFFGetMoleculeForceField(mol, mp)
#     ff.Minimize()
#     return rdkit2ase(mol)


def obmol_to_ase_atoms(obmol):
      """Converts an Open Babel OBMol object to an ASE Atoms object."""
      # Get atom symbols and positions from OBMol
      symbols = ['C' for i in range(obmol.NumAtoms())]
      positions = []
      for i in range(obmol.NumAtoms()):
          vec = obmol.GetAtom(i+1).GetVector()
          positions.append([vec.GetX(), vec.GetY(), vec.GetZ()])
      # Create ASE Atoms object
      ase_atoms = Atoms(symbols=symbols, positions=positions)
      return ase_atoms


def atoms2obmol(atoms, symbol='C', bond_order=1):
    """
    Convert ASE atoms to openbabel mol.
    The bond order is set to 1 to get better minimization results.
    Bond length for the double bonds are too short.
    """
    ana = Analysis(atoms)
    bonds = ana.get_bonds(symbol, symbol, unique=True)[0]
    mol = openbabel.OBMol()
    for p in atoms.positions:
        obatom = mol.NewAtom()
        obatom.SetAtomicNum(6)
        obatom.SetVector(p[0], p[1], p[2])
    for b in bonds:
        mol.AddBond(int(b[0] + 1), int( b[1] + 1), bond_order)
    return mol


def obff_minimize(pdb_file, steps=20, ff='mmff94', algorithm='Conjugate Gradient', st=None):
    atoms = ase.io.read(pdb_file)
    mol = atoms2obmol(atoms)
    # na1 = mol.NumAtoms()
    # mol.AddHydrogens()
    # na2 = mol.NumAtoms()
    obff = openbabel.OBForceField.FindForceField(ff)
    ff_setup = obff.Setup(mol)
    e_init = obff.Energy()
    if algorithm == 'Conjugate Gradient':
        obff.ConjugateGradients(steps)
    elif algorithm == 'Steepest Descent':
        obff.SteepestDescent(steps)
    obff.GetCoordinates(mol)
    e_final = obff.Energy()
    if not ff_setup:
        st.text(f'{ff} setup : {ff_setup} | minimization completed | Energy: {e_init:.2f} > {e_final:.2f}')
    atoms = obmol_to_ase_atoms(mol)
    return atoms
