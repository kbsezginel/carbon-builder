import ase.io
from rdkit import Chem
from rdkit.Chem import AllChem
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from rdkit2ase import rdkit2ase, ase2rdkit


def emt_minimize(atoms):
    """ASE EMT minimization"""
    atoms.write('opt.xyz')
    atoms_new = ase.io.read('opt.xyz')
    atoms_new.calc = EMT()
    # e_init = atoms_new.get_potential_energy()
    dyn = BFGS(atoms_new)
    dyn.run(fmax=0.5)
    # e_final = atoms_new.get_potential_energy()
    # st.text(f'Optimization completed: {e_init:.3f} -> {e_final:.3f}')
    return atoms_new


def mmff_minimize(atoms):
    """rdKit mmff minimization"""
    atoms.write('rd.pdb')
    mol = Chem.MolFromPDBFile('rd.pdb', removeHs=False)  # Set 
    # Chem.AddHs(mol)

    # mol = ase2rdkit(atoms)
    mp = AllChem.MMFFGetMoleculeProperties(mol)
    ff = AllChem.MMFFGetMoleculeForceField(mol, mp)
    ff.Minimize()
    return rdkit2ase(mol)


def obminimize(atoms, steps=200, ff='mmff94'):
    """OpenBabel minimization"""
    atoms.write('tmp.pdb')
    exe = '/Users/kutay.sezginel/anaconda3/envs/flex/bin/obminimize'
    cmd = [exe, '-n', str(steps), '-ff', ff, 'tmp.pdb']
    result = subprocess.run(cmd, capture_output=True, text=True)
    with open('tmp_opt.pdb', 'w') as f:
        f.write(result.stdout)
    # st.text(result.stdout)
    atoms = ase.io.read('tmp_opt.pdb')
    return atoms
    # obminimize -n 200 -ff mmff94 seed2.pdb > seed2_opt.pdb
