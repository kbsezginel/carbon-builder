import numpy as np
from ase.data import covalent_radii as radii, atomic_numbers

class idxListEmpty(Exception):
    pass

class rcDictEmpty(Exception):
    pass

class normValMissing(Exception):
    pass

class invalidProbe(Exception):
    pass

class typeError(Exception):
    pass

def generateIdxList(atoms, probe):
    """Generate list of atomic indices from 'probe' kwarg."""

    # If None / dafault, use all atoms.
    if probe == None:
        idxList = [atom.index for atom in atoms]
    else:
        idxList = []
        # Always convert to list for easier handling.
        if not isinstance(probe, list):
            probe = [probe]
        for item in probe:
            # int: user provided atomic index.
            if isinstance(item, int):
                idxList.append(item)
            # str: user provided element symbol.
            elif isinstance(item, str):
                idxList.extend([atom.index for atom in atoms
                                if atom.symbol == item])
            # Otherwise break (invalid type for probe.)
            else:
                raise invalidProbe(" ".join(
                        ["The supplied probe argument",
                         "'{:s}' is not of".format(str(probe)),
                        "either of the allowed types int or str."]))

    # Sanity check.
    if not idxList:
        raise idxListEmpty(" ".join(["The supplied 'probe' argument",
                                     "produced an empty selection of",
                                     "atomic indices; please refine.",
                                     "Exiting."]))

    # Remove duplicates.
    idxList = list(set(idxList))

    return idxList

def generateRcDict(atoms, rc, factor):
    """Generate dict of cutoff values for each element."""

    # Sanity check.
    if not isinstance(factor, int) and not isinstance(factor, float):
        raise typeError(
            "Argument 'factor' must be of type int or float.")

    rcDict = {}
    # If user provided {element: cutoff} dict:
    if isinstance(rc, dict):
        rcDict = rc
    # If None / default, calculate rc from covalent radii.
    elif rc == None:
        for sym in set([atom.symbol for atom in atoms]):
            rcDict.update({sym: radii[atomic_numbers[sym]] * factor})
    # If int / float, use this value for all elements.
    elif isinstance(rc, float) or isinstance(rc, int):
        for sym in set([atom.symbol for atom in atoms]):
            rcDict.update({sym: rc})
    # Sanity check.
    else:
        raise rcDictEmpty(" ".join(["Supplied 'rc' cutoff must match",
                                    "either of the following formats:",
                                    "{'element_symbol': rc} dict",
                                    "or float value or None (triggers",
                                    "automatic estimation based on",
                                    "covalent radii). Default: None."]))
    return rcDict

def generateNeighborList(atoms, idxList, rcDict):
    """Generate dict of neighbors for supplied idxList."""

    neighborList = {}
    for idx in idxList:
        idxRest = [atom.index for atom in atoms if atom.index != idx]
        dists = atoms.get_distances(idx, idxRest, mic=True)
        rc1 = rcDict[atoms[idx].symbol]
        for i, d in enumerate(dists):
            rc2 = rcDict[atoms[idxRest[i]].symbol]
            if d <= (rc1 + rc2):
                if idx not in neighborList:
                    neighborList.update({idx: [idxRest[i]]})
                else:
                    neighborList[idx].append(idxRest[i])

    return neighborList

def coordination_numbers(atoms, probe=None, rc=None, factor=1.2,
                         generalized=False, norm=None):

    """Calculate the classical coordination number of atoms.

    Parameters:

    atoms : obj
        ASE atoms object containing the to-be-analyzed structure.

    probe : int or str or list or None
        Which atom or atoms to probe for coordination numbers.
        If int, the atom with this index will be analyzed.
        If str, all atoms with this elemental symbol will be analyzed.
        If None (default), all atoms will be analyzed.
        It is also possible to supply a mixed list of indices and
        atomic symbols. Automatically excludes duplicates if selections
        overlap.

    rc : float or dict or None
        Cutoff radius for the estimation of nearest-neighbors.
        If float, use this value for all atoms.
        If dict, specify a cutoff for each element, for example
        rc = {'Pt': 2.0, 'O': 1.7}.
        If None, estimate `rc` for each pair of atoms based on their
        covalent radii. This can often fail for unconventional
        and strained structures with unusual bond lengths.

    factor : float
        A multiplication factor that will be applied to rc if
        rc == None. Default: 1.2, i.e. slightly larger than the
        combined atomic radii.

    generalized : bool
        Calculated generalized coordination number according to
            F. Calle-Vallejo, et al.,
            Angew. Chemie Int. Ed. 2014, 53, 8316-8319.
        instead of the classic coordination number (based on
        counting nearest-neighbors according to `rc`).

    norm : float (optional)
        If generalized == True, supply a norming value. The norm
        value is the maximum number of neighbors an atom in
        the crystal structure can have (fcc = 12, bcc = 8, etc.)

    Return values:

    cnList : dict
        Coordination number for every selected atoms.

    cnDist : dict
        Frequency distribution of coordination numbers.

    Usage example:

        >>> import matplotlib.pyplot as plt
        >>> from ase.cluster import Octahedron
        >>> from asetools.analysis.coordinationNumbers import coordination_numbers
        >>>
        >>> atoms = Octahedron("Pt", 10)
        >>> cnList, cnDist = coordination_numbers(atoms)
        >>>
        >>> keys, values = zip(*cnDist.items())
        >>> plt.bar(keys, values)
        >>> plt.xlabel("Coordination Number")
        >>> plt.ylabel("Frequence of Occurrence")
        >>> plt.show()
    """

    # Check if norm value is supplied when GCNs are calculated.
    if generalized:
        if norm == None:
            error = " ".join(["Norm value is required when calculating",
                             "generalized coordination numbers. Set",
                             "norm to the maximum number of nearest",
                             "neighbors possibe in the crystal",
                             "structure (fcc: 12, bcc: 8, etc.)."])
            raise normValMissing(error)

    # Generate idx list of to-be-analyzed atoms.
    idxList = generateIdxList(atoms, probe)

    # Generate dict of cutoff radii.
    rcDict = generateRcDict(atoms, rc, factor)

    # Make list of neighbors for each atom in idxList.
    neighborList = generateNeighborList(atoms, idxList, rcDict)

    # Calculate coordination number.
    if generalized:
        cnList = {}
        for idx, neighbors in neighborList.items():
            gcn = sum([len(neighborList[nn]) for nn in neighbors]) / norm
            cnList.update({idx: gcn})
    else:
        cnList = dict([(idx, len(neighbors)) for idx, neighbors
                                             in neighborList.items()])

    # Make distribution of coordination numbers.
    cnAll = np.array([*cnList.values()])
    cnUnique = set(cnAll)
    cnDist = dict([(cn, np.count_nonzero(cnAll == cn)) for cn in cnUnique])

    return cnList, cnDist
