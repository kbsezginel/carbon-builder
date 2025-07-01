from ase.build import graphene_nanoribbon
from stmol import showmol, makeobj
import py3Dmol

import streamlit as st


st.set_page_config(page_title='Graphene Builder', layout='wide')
st.title('Graphene Builder')

tmp_file = 'tmp_gn.pdb'

gnr_cols = st.columns(4)
gnr_width = gnr_cols[0].number_input('Width', value=3, min_value=0, placeholder=str(3))
gnr_length = gnr_cols[1].number_input('Length', value=3, min_value=0, placeholder=str(3))
gnr_type = gnr_cols[2].selectbox('Type', ('armchair', 'zigzag'))
gnr_sat = gnr_cols[3].toggle('Saturated', value=False)

gnr = graphene_nanoribbon(gnr_width, gnr_length, type=gnr_type, saturated=gnr_sat, vacuum=0)
gnr.rotate(90, 'x')
gnr.write(tmp_file)

mol = open(tmp_file).read()
viewer = py3Dmol.view()
viewer.addModel(mol,'pdb',{'doAssembly':True,'duplicateAssemblyAtoms':True})
viewer.setStyle({'sphere':{'scale':.3},'stick':{'colorscheme':'Jmol'}})
viewer.zoomTo()

showmol(viewer, height=500, width=800)


# @st.cache_data
def convert_for_download(mol):
    mol.write('tmpp.xyz')
    with open('tmpp.xyz', 'r') as f:
        s = f.read()
    return s


mol_download = convert_for_download(gnr)

st.download_button(
    label="Download Molecule",
    data=mol_download,
    file_name="graphene.xyz",
    icon=":material/download:",
)