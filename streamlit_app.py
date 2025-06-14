import math
import ase.io
from ase import Atom
from ase.geometry.analysis import Analysis
from stmol import showmol, makeobj
import py3Dmol
import numpy as np
import streamlit as st
from io import StringIO
import subprocess
from pathlib import Path
from copy import deepcopy

from minimization import minimize as ffmin
from minimization import FORCE_FIELDS
from coordinationNumbers import coordination_numbers


def detect_perimeter_by_coordination(atoms, cutoff=2):
    cnl, cnd = coordination_numbers(atoms)
    perimeter = []
    for i in cnl:
        if cnl[i] < cutoff + 1:
            perimeter.append(i)
    return perimeter


def detect_perimeter_by_bonds(atoms, cutoff=2, skin=0.3):
    ana = Analysis(atoms, skin=skin)
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

def get_triangles(atoms, skin=0.3):
    patoms, bdict = detect_perimeter_by_bonds(atoms, cutoff=2, skin=skin)
    triangles = {}
    for a in patoms:
        if len(bdict[a]) == 2:
            triangles[a] = bdict[a]
    return triangles

def get_lines(atoms, skin=0.3):
    patoms, bdict = detect_perimeter_by_bonds(atoms, cutoff=1, skin=skin)
    lines = {}
    for a in patoms:
        lines[a] = bdict[a]
    return lines, bdict

def is_perimeter_complete(atoms, skin=0.3):
    patoms, bdict = detect_perimeter_by_bonds(atoms, cutoff=1, skin=skin)
    if len(patoms) == 0:
        return True
    else:
        return False

def add_single_spokes(atoms, angle=60):
    triangles = get_triangles(atoms)
    # st.text(f'Triangles: {len(triangles)} [num atoms: {len(atoms)}]')
    angle_radians = math.radians(angle)
    for t in triangles:
        p1 = atoms[t].position
        p2 = atoms[triangles[t][0]].position
        p3 = atoms[triangles[t][1]].position
        p23 = np.average([p2, p3], axis=0)
        v1 = p1 - p23
        # pnew = p1 + v1 * 2
        pnew = p1 + v1 / math.cos(angle_radians)
        a = Atom('C', pnew)
        atoms.append(a)
    return atoms


def add_dual_spokes(atoms, angle=60):
    lines, bdict = get_lines(atoms)
    # st.text(f'Lines: {len(lines)} [num atoms: {len(atoms)}]')
    angle_radians = math.radians(angle)
    # pnew = p1 + v1 
    for l in lines:
        connected_atom = lines[l][0]
        tri = [i for i in bdict[connected_atom] if i != l]
        if len(tri) > 1:
            p1 = atoms[l].position
            p2 = atoms[connected_atom].position
            p3 = atoms[tri[0]].position
            p4 = atoms[tri[1]].position

            v1 = p1 - p2
            p3n = p3 + v1 / math.cos(angle_radians)
            p4n = p4 + v1 / math.cos(angle_radians)
            atoms.append(Atom('C', p3n))
            atoms.append(Atom('C', p4n))
    return atoms

def remove_duplicates(atoms, cutoff=0.1, average=False):
    delete = []
    new_atoms = []
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
                if average:
                    delete.append(i1)
                    delete.append(i2)
                    new_atoms.append((p1 + p2) / 2)
                else:
                    delete.append(i2)
    keep = [i for i in range(len(atoms)) if i not in delete]
    atoms = atoms[keep]
    for a in new_atoms:
        atoms.append(Atom('C', a))
    # st.text(f'Total: {len(atoms)} Delete: {len(delete)} Keep: {len(keep)}')
    return atoms

def add_perimeter(
    atoms,
    dup_cutoff=1.0,
    dup_average=True,
    angle=60,
    run_minimization=False,
    force_field='MMFF94',
    min_steps=20
):
    atoms_new = add_single_spokes(atoms, angle=angle)
    if run_minimization:
        atoms_new = minimize(atoms_new, force_field)
    atoms_new = add_dual_spokes(atoms_new, angle=angle)
    if run_minimization:
        atoms_new = minimize(atoms_new, force_field)
    atoms_new = remove_duplicates(atoms_new, cutoff=dup_cutoff, average=dup_average)
    if run_minimization:
        atoms_new = minimize(atoms_new, force_field)
    return atoms_new

def complete_perimeter(
    atoms,
    dup_cutoff=1.0,
    dup_average=True,
    angle=60,
    max_attempts=3, 
    run_minimization=False,
    force_field='MMFF94',
    min_steps=20
):
    n = 0
    perimeter_completed = is_perimeter_complete(atoms)
    while not perimeter_completed:
        atoms = add_spoke(atoms, dup_cutoff=dup_cutoff, dup_average=dup_average, angle=angle, single_only=True)
        if run_minimization:
            mol_new = minimize(atoms, force_field)
        perimeter_completed = is_perimeter_complete(atoms)
        n += 1
        if n > max_attempts:
            break
    return atoms

def add_spoke(
    atoms,
    dup_cutoff=1.0,
    dup_average=True,
    angle=60,
    skin=0.3,
    single_only=False,
    run_minimization=False,
    force_field='MMFF94',
    min_steps=20
):
    angle_radians = math.radians(angle)
    patoms, bdict = detect_perimeter_by_bonds(atoms, cutoff=2, skin=skin)
    for a in patoms:
        if len(bdict[a]) == 1:
            connected_atom = bdict[a][0]
            tri = [i for i in bdict[connected_atom] if i != a]
            if len(tri) < 2:
                continue
            p1 = atoms[a].position
            p2 = atoms[connected_atom].position
            p3 = atoms[tri[0]].position
            p4 = atoms[tri[1]].position

            v1 = p1 - p2
            p3n = p3 + v1 / math.cos(angle_radians)
            p4n = p4 + v1 / math.cos(angle_radians)
            atoms.append(Atom('C', p3n))
            atoms.append(Atom('C', p4n))
        elif len(bdict[a]) == 2:
            if single_only:
                continue
            p1 = atoms[a].position
            p2 = atoms[bdict[a][0]].position
            p3 = atoms[bdict[a][1]].position
            p23 = np.average([p2, p3], axis=0)
            v1 = p1 - p23
            pnew = p1 + v1 / math.cos(angle_radians)
            a = Atom('C', pnew)
            atoms.append(a)
    atoms = remove_duplicates(atoms, cutoff=dup_cutoff, average=dup_average)
    if run_minimization:
        mol_new = minimize(atoms, force_field, steps=min_steps)
    return atoms

def remove_spokes(atoms, skin=0.3):
    patoms, bdict = detect_perimeter_by_bonds(atoms, cutoff=1, skin=skin)
    keep_atoms = [i for i in range(len(atoms)) if i not in patoms]
    return atoms[keep_atoms]

def add_bonds_pdb(pdb_file, atoms, skin=0.5):
    ana = Analysis(atoms, skin=skin)
    bonds = ana.get_bonds('C', 'C', unique=True)[0]

    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    new_lines = []
    for l in lines:
        if 'ENDMDL' in l:
            break
        new_lines.append(l)

    for atom in range(1, len(atoms) + 1):
        atom_bonds = [atom]
        for b in bonds:
            if atom == b[0] + 1:
                atom_bonds.append(b[1] + 1)
            elif atom == b[1] + 1:
                atom_bonds.append(b[0] + 1)
        new_lines.append('CONECT' + ' %4i' * len(atom_bonds) % tuple(atom_bonds) + '\n')

    with open(pdb_file, 'w') as f:
        for l in new_lines:
            f.write(l)
        f.write('ENDMDL\n')

def minimize(atoms, force_field='MMFF94', steps=20, skin=0.3, st=st):
    atoms.write('min.pdb')
    add_bonds_pdb('min.pdb', atoms, skin=skin)
    atoms_min = ffmin('min.pdb', force_field, steps=steps, st=st)
    return atoms_min


st.set_page_config(page_title='Perimeter Builder', layout='wide')
st.title('Perimeter Builder')

if 'filename' not in st.session_state:
    st.session_state['filename'] = ''
if 'operations' not in st.session_state:
    st.session_state['operations'] = {}
if 'count' not in st.session_state:
    st.session_state['count'] = []


moldir = Path('molecules')
init_xyz = moldir / 'seed0.xyz'

uploaded_file = st.file_uploader('Upload molecule file')
if uploaded_file is not None:
    if str(uploaded_file.name) != st.session_state['filename']:
        st.session_state['count'] = []
        st.session_state['filename'] = str(uploaded_file.name)
        stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
        mol_str = stringio.read()
        mol_ext = uploaded_file.name.split('.')[-1]
        fname = f'tmppp.{mol_ext}'
        with open(fname, 'w') as f:
            f.write(mol_str)
        st.session_state['init_mol'] = ase.io.read(fname)
        st.session_state['mol'] = ase.io.read(fname)
else:
    if st.session_state['filename'] != 'default':
        # init_mol = ase.io.read(init_xyz)
        st.session_state['mol'] = ase.io.read(init_xyz)
        st.session_state['init_mol'] = ase.io.read(init_xyz)
        st.session_state['filename'] = 'default'
        # st.session_state['reload_file'] = 'default'

bond_skin_help = """Skin distance used to calculate bonding.
The larger the distance the atoms that are further away will be considered bonded."""

spoke_angle_help = """The angle for the new spoke atoms. Default is 60 degrees based on hexogonal arrangement."""
merge_help = """When two atoms are closer to each other than merge cutoff, only one of those atoms will be kept.
If merge to the middle is selected, the average position between these atoms will be kept.
If merge to the middle is not selected, only the first atom will be kept. The positions will not be averaged."""

merge_cutoff_help = """If two atoms are closer to each other than this value the atoms will be merged.
"""

cols1 = st.columns(7, vertical_alignment='center')
# tri_index = cols[0].number_input('Triangle', value=0, min_value=0, placeholder=str(3))
add_perimeter_btn = cols1[0].button('Add Perimeter', use_container_width=True)
add_spoke_btn = cols1[1].button('Add Spokes', use_container_width=True)
remove_spokes_btn = cols1[2].button('Remove spokes', use_container_width=True)
complete_perimeter_btn = cols1[3].button('Complete Perimeter', use_container_width=True)
minimize_btn = cols1[4].button('Minimize', use_container_width=True)
undo_btn = cols1[5].button('Undo', use_container_width=True)
reload_btn = cols1[6].button('Reload', use_container_width=True)

# cols1b = st.columns(7)

with st.expander('Settings'):

    cols2 = st.columns(5, vertical_alignment='center')
    # show_perimeters = cols2[0].toggle('Show perimeter atoms', value=True)
    spoke_angle = 60
    # spoke_angle = cols2[0].number_input('Spoke angle', value=60, help=spoke_angle_help)
    force_field = cols2[0].selectbox('Force Field', FORCE_FIELDS, index=0)
    min_steps = cols2[1].number_input('Minimization Steps', value=20, min_value=0)
    bond_skin_dist = cols2[2].number_input('Bond skin distance', value=0.3, min_value=0.0, help=bond_skin_help)
    dup_cutoff = cols2[3].number_input('Merge cutoff (Å)', value=1.0, min_value=0.0, help=merge_help)
    minimize_every_step = cols2[4].toggle('Minimize every step', value=False)
    dup_average = True
    # dup_average = cols2[4].toggle('Merge nearby atoms to the middle', value=True, help=merge_help)

    cols3 = st.columns(6)
    perimeter_color = cols3[0].color_picker("Perimeter Color", "#E20000")
    atom_color = cols3[1].color_picker("Atom Color", "#373737")
    bond_color = cols3[2].color_picker("Bond Color", "#373737")
    bg_color = cols3[3].color_picker("Background Color", "#ffffff")
    mol_viewer_width = cols3[4].number_input('Width', value=1100, min_value=0)
    mol_viewer_height = cols3[5].number_input('Height', value=600, min_value=0)


# if add_perimeter_btn or add_spoke_btn or undo_btn:
if reload_btn:
    st.session_state['mol'] = st.session_state['init_mol']
    st.session_state['count'] = []
if add_perimeter_btn:
    st.session_state['count'].append('P')
    st.session_state['mol'] = add_perimeter(
        st.session_state['mol'],
        dup_cutoff=dup_cutoff,
        dup_average=dup_average,
        angle=spoke_angle,
        run_minimization=minimize_every_step,
        force_field=force_field,
        min_steps=min_steps
    )
if add_spoke_btn:
    st.session_state['count'].append('S')
    st.session_state['mol'] = add_spoke(
        st.session_state['mol'],
        dup_cutoff=dup_cutoff,
        dup_average=dup_average,
        angle=spoke_angle,
        run_minimization=minimize_every_step,
        force_field=force_field,
        min_steps=min_steps
    )
if remove_spokes_btn:
    st.session_state['count'].append('RS')
    st.session_state['mol'] = remove_spokes(st.session_state['mol'])
if minimize_btn:
    st.session_state['count'].append('M')
    # st.session_state['mol'].write('min.pdb')
    # add_bonds_pdb('min.pdb', st.session_state['mol'], skin=bond_skin_dist)
    # st.session_state['mol'] = mmff_minimize('min.pdb')
    # st.session_state['mol'] = minimize('min.pdb', force_field, steps=min_steps)
    st.session_state['mol'] = minimize(st.session_state['mol'], force_field, steps=min_steps)
if complete_perimeter_btn:
    st.session_state['count'].append('C')
    st.session_state['mol'] = complete_perimeter(
        st.session_state['mol'],
        dup_cutoff=dup_cutoff,
        dup_average=dup_average,
        angle=spoke_angle,
        run_minimization=minimize_every_step,
        force_field=force_field,
        min_steps=min_steps
    )
if undo_btn:
    st.session_state['count'].pop()
    st.session_state['mol'] = deepcopy(st.session_state['init_mol'])
    for c in st.session_state['count']:
        if c == 'P':
            st.session_state['mol'] = add_perimeter(
                st.session_state['mol'],
                dup_cutoff=dup_cutoff,
                dup_average=dup_average,
                angle=spoke_angle,
                run_minimization=minimize_every_step,
                force_field=force_field,
                min_steps=min_steps
            )
        if c == 'S':
            st.session_state['mol'] = add_spoke(
                st.session_state['mol'],
                dup_cutoff=dup_cutoff,
                dup_average=dup_average,
                angle=spoke_angle,
                run_minimization=minimize_every_step,
                force_field=force_field,
                min_steps=min_steps
            )
        if c == 'RS':
            st.session_state['mol'] = remove_spokes(st.session_state['mol'])
        if c == 'C':
            st.session_state['mol'] = complete_perimeter(
                st.session_state['mol'],
                dup_cutoff=dup_cutoff,
                dup_average=dup_average,
                angle=spoke_angle,
                run_minimization=minimize_every_step,
                force_field=force_field,
                min_steps=min_steps
            )
        if c == 'M':
            st.session_state['mol'].write('min.pdb')
            add_bonds_pdb('min.pdb', st.session_state['mol'], skin=bond_skin_dist)
            st.session_state['mol'] = minimize('min.pdb', force_field, steps=min_steps, st=st)


ana = Analysis(st.session_state['mol'])
ccbonds = ana.get_bonds('C', 'C', unique=True)
nb = len(ccbonds[0])
pc = st.session_state['count'].count('P')
sc = st.session_state['count'].count('S')
nc = st.session_state['count'].count('C')
na = len(st.session_state['mol'])
disp_txt = f'Perimeters: {pc} - Spokes: {sc} - Num. Atoms: {na} - Num Bonds: {nb}'
st.text(disp_txt)

# Save pdb file with bonds
st.session_state['mol'].write('tmp.pdb')
add_bonds_pdb('tmp.pdb', st.session_state['mol'], skin=bond_skin_dist)
pdbfile = open('tmp.pdb').read()

viewer = py3Dmol.view(height=mol_viewer_height, width=mol_viewer_width)
viewer.addModel(pdbfile,'pdb',{'doAssembly':True,'duplicateAssemblyAtoms':True})
viewer.setStyle({'sphere':{'scale':.3, 'color': atom_color},'stick':{'color': bond_color}})
viewer.setBackgroundColor(bg_color)

# Display molecular and perimeter atoms
patoms, bdict = detect_perimeter_by_bonds(st.session_state['mol'], skin=bond_skin_dist)
for a in patoms:
    viewer.addStyle({'model': -1, 'serial': a + 1}, {'sphere': {'color': perimeter_color, 'radius': 2.}})

# viewer.setStyle({'sym':2},{'sphere':{'scale':.5,'color':'blue'},'stick':{'color':'cyan'}})
# viewer.addStyle({'model': -1, 'serial': 1}, {'sphere': {'color': 'red', 'radius': 2.5}})
viewer.zoomTo()

showmol(viewer, height=mol_viewer_height, width=mol_viewer_width)


# @st.cache_data
def convert_for_download(mol):
    mol.write('tmpp.xyz')
    with open('tmpp.xyz', 'r') as f:
        s = f.read()
    return s


mol_download = convert_for_download(st.session_state['mol'])

st.download_button(
    label="Download Molecule",
    data=mol_download,
    file_name="mol.xyz",
    icon=":material/download:",
)
