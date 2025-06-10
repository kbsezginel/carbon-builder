import ase.io
from ase import Atom
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.build import graphene_nanoribbon
from ase.geometry.analysis import Analysis
from stmol import showmol, makeobj
import py3Dmol
import numpy as np
import streamlit as st
from io import StringIO
import subprocess

from coordinationNumbers import coordination_numbers


def detect_perimeter_by_coordination(atoms, cutoff=2):
    cnl, cnd = coordination_numbers(atoms)
    perimeter = []
    for i in cnl:
        if cnl[i] < cutoff + 1:
            perimeter.append(i)
    return perimeter


def detect_perimeter_by_bonds(atoms, cutoff=2):
    ana = Analysis(atoms)
    ccbonds = ana.get_bonds('C', 'C')
    bdict = {}
    for b in ccbonds[0]:
        b1, b2 = int(b[0]), int(b[1])
        if b1 not in bdict:
            bdict[b1] = []
        bdict[b1].append(b2)
        if b2 not in bdict:
            bdict[b2] = []
        bdict[b2].append(b1)
    bonds = []
    for b in bdict:
        if len(bdict[b]) < cutoff + 1:
            bonds.append(b)
    return bonds, bdict

def get_triangles(atoms):
    patoms, bdict = detect_perimeter_by_bonds(atoms, cutoff=2)
    triangles = {}
    for a in patoms:
        triangles[a] = bdict[a]
    return triangles

def get_lines(atoms):
    patoms, bdict = detect_perimeter_by_bonds(atoms, cutoff=1)
    lines = {}
    for a in patoms:
        lines[a] = bdict[a]
    return lines, bdict

def add_spokes(atoms):
    triangles = get_triangles(atoms)
    for t in triangles:
        p1 = atoms[t].position
        p2 = atoms[triangles[t][0]].position
        p3 = atoms[triangles[t][1]].position
        p23 = np.average([p2, p3], axis=0)
        v1 = p1 - p23
        pnew = p1 + v1 * 2
        a = Atom('C', pnew)
        atoms.append(a)
    return atoms

def add_dual_spokes(atoms):
    lines, bdict = get_lines(atoms)
    for l in lines:
        connected_atom = lines[l][0]
        tri = [i for i in bdict[connected_atom] if i != l]
        p1 = atoms[l].position
        p2 = atoms[connected_atom].position
        p3 = atoms[tri[0]].position
        p4 = atoms[tri[1]].position

        v1 = p1 - p2
        p3n = p3 + v1 * 2
        p4n = p4 + v1 * 2
        atoms.append(Atom('C', p3n))
        atoms.append(Atom('C', p4n))
    return atoms

def remove_duplicates(atoms, cutoff=0.1):
    delete = []
    for i1, a1 in enumerate(atoms):
        if i1 in delete:
            continue
        p1 = a1.position
        for i2, a2 in enumerate(atoms):
            if i1 >= i2:
                continue
            if i2 in delete:
                continue
            p2 = a2.position
            d = np.linalg.norm(p2 - p1)
            if d < cutoff:
                delete.append(i2)
    keep = [i for i in range(len(mol)) if i not in delete]
    # st.text(f'Total: {len(atoms)} Delete: {len(delete)} Keep: {len(keep)}')
    return atoms[keep]

def optimize(atoms):
    atoms.write('opt.xyz')
    atoms_new = ase.io.read('opt.xyz')
    atoms_new.calc = EMT()
    e_init = atoms_new.get_potential_energy()
    dyn = BFGS(atoms_new)
    dyn.run(fmax=0.5)
    e_final = atoms_new.get_potential_energy()
    st.text(f'Optimization completed: {e_init:.3f} -> {e_final:.3f}')
    return atoms_new

def obminimize(atoms, steps=200, ff='mmff94'):
    atoms.write('tmp.pdb')
    exe = '/Users/kutay.sezginel/anaconda3/envs/flex/bin/obminimize'
    cmd = [exe, '-n', str(steps), '-ff', ff, 'tmp.pdb', '>', 'tmp_opt.pdb']
    result = subprocess.run(cmd, capture_output=True, text=True)
    with open('tmp_opt.pdb', 'w') as f:
        f.write(result.stdout)
    # st.text(result.stdout)
    atoms = ase.io.read('tmp_opt.pdb')
    return atoms
    # obminimize -n 200 -ff mmff94 seed2.pdb > seed2_opt.pdb


def spoke(mol):
    mol_new = add_spokes(mol)
    # st.text(f'{len(mol)}')
    mol_new = add_dual_spokes(mol_new)
    # st.text(f'{len(mol)}')
    mol_new = remove_duplicates(mol_new, cutoff=0.5)
    # st.text(f'{len(mol)}')
    # mol_new = optimize(mol_new)
    return mol_new



st.set_page_config(page_title='Perimeter Builder')
st.title('Perimeter Builder')

if 'filename' not in st.session_state:
    st.session_state['filename'] = ''

uploaded_file = st.file_uploader('Upload molecule file')
if uploaded_file is not None:
    if str(uploaded_file.name) != st.session_state['filename']:
        st.session_state['count'] = 0
        st.session_state['filename'] = str(uploaded_file.name)
    # st.text(uploaded_file.name)
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
    xyz_str = stringio.read()
    fname = 'tmppp.xyz'
    with open(fname, 'w') as f:
        f.write(xyz_str)
    # st.write(stringio)
    mol = ase.io.read(fname)
    st.session_state['mol'] = mol
else:
    mol = ase.io.read('seed0.xyz')

if 'count' not in st.session_state:
    st.session_state['count'] = 0
if 'mol' not in st.session_state:
    st.session_state['mol'] = mol

cols = st.columns(3)
# tri_index = cols[0].number_input('Triangle', value=0, min_value=0, placeholder=str(3))
add_button = cols[0].button('Add Perimeter')
show_perimeters = cols[1].toggle('Show perimeter atoms', value=True)
# num = cols[1].number_input('Num', value=1, min_value=0, placeholder=str(3))
detection_type = cols[2].selectbox('Perimeter Detection', ('Coordination', 'Bonds'))


# get perimeter atoms
if detection_type == 'Bonds':
    patoms, bdict = detect_perimeter_by_bonds(mol)
elif detection_type == 'Coordination':
    patoms = detect_perimeter_by_coordination(mol)
    _, bdict = detect_perimeter_by_bonds(mol)

if add_button:
    st.session_state['count'] += 1
    for i in range(st.session_state['count']):
        mol = spoke(mol)
    # mol = obminimize(mol)
    st.session_state['mol'] = mol

    # st.session_state['mol'] = spoke(mol)
    # st.session_state['mol'] = spoke(st.session_state['mol'])
    # st.text('SPOKED: ' + str(len(st.session_state['mol'])) + ' ' + str(len(mol)))

st.text('Perimeter Count: ' + str(st.session_state['count']) + '  Num Atoms: ' + str(len(mol)))


# tri_str = ' '.join([str(i) for i in triangles])
# st.text(f'Triangles: {tri_str}')
# tri_str = ' '.join([str(i) for i in lines])
# st.text(f'Lines: {tri_str}')


patoms, bdict = detect_perimeter_by_bonds(st.session_state['mol'])

st.session_state['mol'].write('tmp.pdb')

pdbfile = open('tmp.pdb').read()
viewer = py3Dmol.view()
viewer.addModel(pdbfile,'pdb',{'doAssembly':True,'duplicateAssemblyAtoms':True})
viewer.setStyle({'sphere':{'scale':.3},'stick':{'colorscheme':'Jmol'}})

if show_perimeters:
    for a in patoms:
        viewer.addStyle({'model': -1, 'serial': a + 1}, {'sphere': {'color': 'red', 'radius': 2.}})

# triangles = get_triangles(mol)
# if tri_index < len(triangles):
#     t = list(triangles)[tri_index]
#     for b in triangles[t]:
#         viewer.addStyle({'model': -1, 'serial': b + 1}, {'sphere': {'color': 'blue', 'radius': 2.}})



# tri_str = ' '.join([str(i) for i in triangles])
# st.text(f'Triangles: {tri_str}')

# viewer.setStyle({'sym':2},{'sphere':{'scale':.5,'color':'blue'},'stick':{'color':'cyan'}})
# viewer.addStyle({'model': -1, 'serial': 1}, {'sphere': {'color': 'red', 'radius': 2.5}})
viewer.zoomTo()

showmol(viewer, height=500, width=800)


# @st.cache_data
def convert_for_download(mol):
    mol.write('tmpp.xyz')
    with open('tmpp.xyz', 'r') as f:
        s = f.read()
    return s


mol_download = convert_for_download(mol)

st.download_button(
    label="Download Molecule",
    data=mol_download,
    file_name="mol.xyz",
    icon=":material/download:",
)

