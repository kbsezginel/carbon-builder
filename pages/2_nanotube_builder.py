import py3Dmol
import streamlit as st
from stmol import showmol
from ase.build import nanotube

st.set_page_config(page_title='Nanotube Builder')
st.title('Nanotube Builder')

nt_cols = st.columns(5)
nt_n = nt_cols[0].number_input('n', value=6, min_value=0, placeholder=str(6))
nt_m = nt_cols[1].number_input('m', value=6, min_value=0, placeholder=str(6))
nt_length = nt_cols[2].number_input('Length', value=6, min_value=0, placeholder=str(6))
nt_bond = nt_cols[3].number_input('Bond Length', value=1.42, min_value=0.0, placeholder=str(1.42))
nt_symbol = nt_cols[4].text_input('Atomic symbol', 'C')

nt = nanotube(nt_n, nt_m, nt_length, nt_bond, nt_symbol)
# gnr.rotate(90, 'x')
tmp_file = 'tmp_nt.pdb'
nt.write(tmp_file)

mol = open(tmp_file).read()
viewer = py3Dmol.view()
viewer.addModel(mol,'pdb',{'doAssembly':True,'duplicateAssemblyAtoms':True})
viewer.setStyle({'sphere':{'scale':.3},'stick':{'colorscheme':'Jmol'}})
viewer.zoomTo()


# rotate_cols = st.columns(3)

# if rotate_cols[0].button('x'):
#     viewer.rotate(90, 'x')
# if rotate_cols[1].button('y'):
#     viewer.rotate(90, 'y')
# if rotate_cols[2].button('z'):
#     viewer.rotate(90, 'z')


# xr = rotate_cols[0].number_input('x', value=90, min_value=0)
# viewer.rotate(xr, 'x')
# yr = rotate_cols[1].number_input('y', value=90, min_value=0)
# viewer.rotate(yr, 'y')
# zr = rotate_cols[2].number_input('z', value=90, min_value=0)
# viewer.rotate(zr, 'z')

showmol(viewer, height=500, width=800)





# @st.cache_data
def convert_for_download(mol):
    mol.write('tmpp.xyz')
    with open('tmpp.xyz', 'r') as f:
        s = f.read()
    return s


mol_download = convert_for_download(nt)

st.download_button(
    label="Download Molecule",
    data=mol_download,
    file_name="nanotube.xyz",
    icon=":material/download:",
)