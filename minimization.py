import os
import ase.io
import subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel
from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.geometry.analysis import Analysis
from rdkit2ase import rdkit2ase, ase2rdkit


FORCE_FIELDS = ['mmff94', 'Ghemical', 'UFF', 'EMT', 'MMFF (rdKit)']


def minimize(pdb_file, force_field, steps=20, st=None):
    if force_field in ['mmff94', 'Ghemical', 'UFF']:
        atoms = obff_minimize(pdb_file, steps=steps, ff=force_field, st=st)
    elif force_field == 'EMT':
        atoms = emt_minimize(pdb_file)
    elif force_field == 'MMFF (rdKit)':
        atoms = mmff_minimize(pdb_file)
    else:
        raise Exception(f'Force field {force_field} not available')
    return atoms


def emt_minimize(pdb_file):
    """ASE EMT minimization"""
    atoms_new = ase.io.read(pdb_file)
    atoms_new.calc = EMT()
    # e_init = atoms_new.get_potential_energy()
    dyn = BFGS(atoms_new)
    dyn.run(fmax=0.5)
    # e_final = atoms_new.get_potential_energy()
    # st.text(f'Optimization completed: {e_init:.3f} -> {e_final:.3f}')
    return atoms_new


def mmff_minimize(pdb_file):
    """rdKit mmff minimization"""
    # atoms.write('rd.pdb')
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)  # Set 
    # Chem.AddHs(mol)

    # mol = ase2rdkit(atoms)
    mp = AllChem.MMFFGetMoleculeProperties(mol)
    ff = AllChem.MMFFGetMoleculeForceField(mol, mp)
    ff.Minimize()
    return rdkit2ase(mol)


def obminimize(pdb_file, steps=20, ff='MMFF94', st=None):
    """OpenBabel minimization"""
    OB_EXE = '/Users/kutay.sezginel/anaconda3/envs/flex/bin/obminimize'
    cmd = ['obminimize', '-ff', 'Ghemical', '-n', str(steps), 'test.mol2']

    exe = OB_EXE
    cmd = [exe, '-n', str(steps), '-ff', ff, pdb_file]
    # obabel rd.pdb -o pdb --minimize --steps 15 --ff MMFF94
    # obabel infile.xxx -O outfile.yyy --minimize --steps 1500 --sd
    cmd = [OB_EXE, pdb_file, '-Otmp_opt.pdb', '--minimize'] # , '--steps', str(steps), '--ff', ff]
    st.text(cmd)
    result = subprocess.run(cmd, capture_output=True, text=True)
    with open('tmp_opt.pdb', 'w') as f:
        f.write(result.stdout)
    # st.text(result.stdout)
    st.text(result.stdout)
    st.text(result.stderr)
    atoms = ase.io.read('tmp_opt.pdb')
    return atoms


def obmol_to_ase_atoms(obmol):
      """Converts an Open Babel OBMol object to an ASE Atoms object."""
      # Get atom symbols and positions from OBMol
      symbols = ['C' for i in range(obmol.NumAtoms())]
      positions = []
      for i in range(obmol.NumAtoms()):
          vec = obmol.GetAtom(i+1).GetVector()
          positions.append([vec.GetX(), vec.GetY(), vec.GetZ()])
      # positions = np.array([obmol.GetAtom(i+1).GetVector() for i in range(obmol.NumAtoms())])
      # Create ASE Atoms object
      ase_atoms = Atoms(symbols=symbols, positions=positions)
      return ase_atoms


def atoms2obmol(atoms):
    ana = Analysis(atoms)
    bonds = ana.get_bonds('C', 'C', unique=True)[0]
    mol = openbabel.OBMol()
    for p in atoms.positions:
        obatom = mol.NewAtom()
        obatom.SetAtomicNum(6)
        obatom.SetVector(p[0], p[1], p[2])
    for b in bonds:
        obbond = mol.NewBond()
        obbond.SetBegin( mol.GetAtom( int(b[0] + 1) ) )
        obbond.SetEnd(mol.GetAtom( int( b[1] + 1) ) )
    return mol


def obff_minimize(pdb_file, steps=20, ff='mmff94', algorithm='Conjugate Gradient', st=None):
    atoms = ase.io.read(pdb_file)
    mol = atoms2obmol(atoms)
    obff = openbabel.OBForceField.FindForceField(ff)
    obff.Setup(mol)
    e_init = obff.Energy()
    if algorithm == 'Conjugate Gradient':
        obff.ConjugateGradients(steps)
    elif algorithm == 'Steepest Descent':
        obff.SteepestDescent(steps)
    obff.GetCoordinates(mol)
    e_final = obff.Energy()
    st.text(f'{ff} minimization completed | Energy: {e_init:.2f} > {e_final:.2f}')
    atoms = obmol_to_ase_atoms(mol)
    return atoms


def obminimize(pdb_file, steps=20, ff='mmff94'):
    """OpenBabel minimization"""
    OB_EXE = '/Users/kutay.sezginel/anaconda3/envs/flex/bin/obminimize'
    cmd = [OB_EXE, '-n', str(steps), '-ff', ff, pdb_file]
    result = subprocess.run(cmd, capture_output=True, text=True)
    with open('tmp_min.pdb', 'w') as f:
        f.write(result.stdout)
    atoms = ase.io.read('tmp_min.pdb')
    return atoms