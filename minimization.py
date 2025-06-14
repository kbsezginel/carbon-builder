import subprocess
import ase.io
from rdkit import Chem
from rdkit.Chem import AllChem
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from rdkit2ase import rdkit2ase, ase2rdkit
import os

FORCE_FIELDS = ['MMFF94', 'Ghemical', 'UFF', 'EMT', 'MMFF (rdKit)']
# OB_EXE = 'obminimize'
# OB_EXE = '/Users/kutay.sezginel/anaconda3/envs/flex/bin/obminimize'

import os
from openbabel import openbabel
print('Openbabel Imported ', openbabel.__file__)
s = openbabel.__file__
lib_id = s.split('/').index('lib')
env_dir = '/'.join(s.split('/')[:lib_id])
print(env_dir)
print(os.listdir(env_dir))
print(f'{env_dir}/bin')
print(os.listdir(f'{env_dir}/bin'))
OB_EXE = f'{env_dir}/bin/obabel'
print(lib_id, env_dir, OB_EXE)


def minimize(pdb_file, force_field, steps=20, st=None):
    if force_field in ['MMFF94', 'Ghemical', 'UFF']:
        atoms = obminimize(pdb_file, steps=steps, ff=force_field, st=st)
    elif force_field == 'EMT':
        atoms = emt_minimize(pdb_file)
    elif force_field == 'MMFF (rdKit)':
        atoms = mmff_minimize(pdb_file)
    else:
        raise Exception(f'Force field {force_field} not available')
    return atoms


def emt_minimize(pdb_file):
    """ASE EMT minimization"""
    # atoms.write('opt.xyz')
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


# def obminimize(pdb_file, steps=20, ff='MMFF94', st=None):
#     """OpenBabel minimization"""
#     exe = OB_EXE
#     cmd = [exe, '-n', str(steps), '-ff', ff, pdb_file]
#     # obabel rd.pdb -o pdb --minimize --steps 15 --ff MMFF94
#     # obabel infile.xxx -O outfile.yyy --minimize --steps 1500 --sd
#     cmd = [OB_EXE, pdb_file, '-Otmp_opt.pdb', '--minimize'] # , '--steps', str(steps), '--ff', ff]
#     st.text(cmd)
#     result = subprocess.run(cmd, capture_output=True, text=True)
#     with open('tmp_opt.pdb', 'w') as f:
#         f.write(result.stdout)
#     # st.text(result.stdout)
#     st.text(result.stdout)
#     st.text(result.stderr)
#     atoms = ase.io.read('tmp_opt.pdb')
#     return atoms


def obminimize(pdb_file, steps=20, ff='MMFF94', st=None):
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat("pdb")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, pdb_file)
    st.text(mol)

    ff = openbabel.OBForceField.FindForceField(ff)
    ff.ConjugateGradients(steps)
    st.text(mol)
    st.text(ff)

    obConversion.SetOutFormat("pdb")
    s = obConversion.WriteString(mol)
    with open("ob_output.pdb", "w") as f:
        f.write(s)
    st.text(s)
    st.text(os.listdir())
    # obConversion.WriteFile(mol, "ob_output.pdb")
    
    atoms = ase.io.read("ob_output.pdb")
    return atoms