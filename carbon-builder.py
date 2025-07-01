import numpy as np
import streamlit as st
from pathlib import Path
import streamlit.components.v1 as components

def ChangeWidgetFontSize(wgt_txt, wch_font_size = '12px'):
    htmlstr = """<script>var elements = window.parent.document.querySelectorAll('*'), i;
                    for (i = 0; i < elements.length; ++i) { if (elements[i].innerText == |wgt_txt|) 
                        { elements[i].style.fontSize='""" + wch_font_size + """'; elements[i].style.fontWeight ='bold'; } } </script>  """

    htmlstr = htmlstr.replace('|wgt_txt|', "'" + wgt_txt + "'")
    components.html(f"{htmlstr}", height=0, width=0)

st.set_page_config(page_title='Carbon Builder', layout='wide')
st.title('Carbon Builder')

"""
Welcome to Carbon Builder tools by Geoff Clark and Kutay Sezginel.

See [nanotubeintersections.com](https://nanotubeintersections.com/) for molecular structures and models designed by Geoff.

# Documentation
"""

with st.expander('Perimeter Builder'):
    """
    Perimeter builder is intended to expand the perimeter of the structure by adding new atoms.
    The perimeter atoms are defined as atoms two 2 bonds or less.
    Currently the perimeter builder is designed to work with hexangonal structures (graphene, CNT, etc.).
    By adjusting the settings and using minimization other geometries can be created in some situations.

    The new atom positions are calculated based on the perimeter atom's number of bonds.

    If the perimeter atom is connected to two atoms:
    - One new atom (B) is added

    If the perimeter atom is connected to one atom:
    - Two new atoms (B1, B2) are added

    While building a perimeter this would cause some atoms to overlap.
    In these cases overlapping atoms are merged by averaging their position.
    By default atoms that are closer than 1 Angstrom are merged.
    The merge cutoff can be adjusted from the settings.
    """

    st.image('pages/spoke_additions_graphic.png', width=700)

    """
    ### Buttons

    #### Add Perimeter
    First adds spokes for two bonded atoms.
    Runs minimization in between if selected.
    Adds spokes for single bonded atoms.
    Removes duplicate atoms.

    #### Add Spokes
    For every atom with a single bond:
    - Pick the atom it's bonded to, let's call this atom B
    - If atom B is also only connected to one atom, skip
    - If B is connected to at least two atoms, 
    create 2 new spokes using the triangle geometry

    For every atom with two bonds:
    - Using the positions of the two bonded atoms (B1, B2) calculate spoke vector
    - Add new atom to the calculated position

    #### Remove Spokes

    #### Complete Parameter
    Until each atom has at least two bonds it keeps trying to add single spokes.
    If the perimeter is not completed after 3 iterations it stops trying.

    #### Minimize
    Minimize molecule geometry. Minimization parameters can be adjusted from settings.

    #### Undo
    Undo the last operation. The perimeter builder keeps a history of each operation performed for a unique molecule.

    #### Reload
    Reloads the molecule and deletes operation history.

    ### Settings

    #### Force Field
    3 minimization engines and 7 force fields are available:
    - 'mmff94', 'mmff94s', 'ghemical', 'gaff', 'uff' (uses openbabel which is what Avogadro uses, seems to work quite well)
    - 'EMT' -> uses ASE (doesn't work well with carbon based structures)
    - 'MMFF (rdKit)' -> uses rdKit (doesn't work well with carbon based structures)

    #### Minimization Steps
    Number of steps to run during minimization.

    #### Bond Skin Distance
    Skin distance used to calculate bonding.
    The larger the distance the atoms that are further away will be considered bonded.

    #### Merge Cutoff
    When two atoms are closer to each other than the merge cutoff, only one of those atoms will be kept.

    #### Min Bond Length (Å)
    The minimum length of the new bonds. It applies to all operations.
    In some cases the spoke vectors can be short and this helps with keeping the new bons at a more consistent length.

    #### Minimize Every Step
    Minimizes geometry after every operation. If the minimization parameters are changed after performing operations
    """

with st.expander('Graphene Builder'):
    """
    Creates a graphene nanoribbon in the x-z plane, with the nanoribbon running along the z axis.

    #### Width
    The width of the nanoribbon. For armchair nanoribbons, this n may be half-integer to repeat by half a cell.

    #### Length
    The length of the nanoribbon.

    #### Type
    The orientation of the ribbon. Must be either ‘zigzag’ or ‘armchair’.

    #### Saturated
    If true, hydrogen atoms are placed along the edge.
    """


with st.expander('Nanotube Builder'):
    """
    Creates a single-walled nanotube whose structure is specified using the standardized (n, m) notation.
    """

    st.image('pages/cnt_nm.png', width=700)

    """
    #### n
    n in the (n, m) notation.

    #### m
    m in the (n, m) notation.

    #### Length
    Length (axial repetitions) of the nanotube.

    #### Bond Length
    Bond length between neighboring atoms.

    #### Atomic Symbol
    Chemical element to construct the nanotube from.

    """

ChangeWidgetFontSize('Nanotube Builder', '26px')
ChangeWidgetFontSize('Graphene Builder', '26px')
ChangeWidgetFontSize('Perimeter Builder', '26px')
