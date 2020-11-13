"""Pre-SmileMerge_mcs_Filter

This script should take an input of a randomly selected file containing a list
of smiles.

A random number function will be used to pick 2 non-identical numbers for 0 to
the len(smile_list) Then those numbers are used to grab 2 non-identical
molecules from smile_list Those 2 smiles are tested using the MCS function to
find the most common structure (MCS). If MCS returns None (ie. no shared
structures) then mol2 is reassigned using the random function generator. This
iterates until a shared structure is returned.
"""

import copy
import random
import rdkit
from rdkit import Chem
from rdkit.Chem import rdFMCS

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

try:
    import autogrow4.autogrow.operators.crossover.smiles_merge.smiles_merge as smiles_merge
    import autogrow4.autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling as MOH
except Exception as e:
    import autogrow.operators.crossover.smiles_merge.smiles_merge as smiles_merge
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling as MOH




def test_for_mcs(vars, mol_1, mol_2):
    """
    Takes a ligand and a random selected molecule and looks for the Most
    common structure. rdFMCS.FindMCS(mols) flags an error if there are no
    common structures being compared. Try/except statement used to prevent
    program crash when 2 unlike molecules are compared. mols is a list of the
    molecules to be compare using rdFMCS.FindMCS

    This can function with mol_1 and mol_2 having H's but we recommend using
    this function with molecules with H's removed. If implicit H's are added
    they will be recoginized as part of MCS. This means 1 atom in common with
    3H's in common would pass a atom similarity count of 4 atoms shared, but
    really its only 1 non-H atom...

    Inputs:
    :param dict vars: User variables which will govern how the programs runs
    :param rdkit.Chem.rdchem.Mol mol_1: the 1st rdkit molecule
    :param rdkit.Chem.rdchem.Mol mol_2: the 2nd rdkit molecule

    Returns:
    :returns: <class 'rdkit.Chem.rdFMCS.MCSResult'> results: an MCSResults
        object returns None if it fails to find MCS sufficient with the User
        defined parameters.
    """

    mols = [mol_1, mol_2]
    time_timeout = vars["max_time_mcs_prescreen"]
    min_number_atoms_matched = vars["min_atom_match_mcs"]

    try:
        result = rdFMCS.FindMCS(
            mols,
            matchValences=False,
            ringMatchesRingOnly=True,
            completeRingsOnly=False,
            timeout=time_timeout,
        )
    except:
        return None

    # could be used for a theoretical timeout prefilter was to be implement
    # (but this isn't implemented) (ie. if it times out the prefilter dont use
    # in thorough MCS ligmerge) canceled: if True, the MCS calculation did not
    # finish

    # filter by UserDefined minimum number of atoms found. The higher the
    # number the more similar 2 ligands are but the more restrictive for
    # finding mergable ligands number of atoms in common found
    if result.numAtoms < min_number_atoms_matched:
        return None
    if result.canceled is True:
        return None

    return result


def find_random_lig2(vars, ligands_list, ligand1_pair):
    """
    Pick a random molecule from the list and check that it can be converted
    into a rdkit mol object and then test for a satistifactory Most common
    substructure (MCS) which satisifies the User specified minimum shared
    substructure

    NECESSARY INCASE THE SMILE CANNOT BE USED (ie. valence issue)

    Inputs:
    :param dict vars: User variables which will govern how the programs runs
    :param list ligands_list: list of all the lignads to chose from
    :param list ligand1_pair: information for the Ligand 1. This info includes
        the name and SMILES string

    Returns:
    :returns: list mol2_pair: a set of information for a 2nd ligand (Lig2)
        This includes the name and SMILES string this mol is from the ligands_list
    :returns: bool bool: returns False if no satistifactory matches were found
        it returns False
    """

    count = 0
    shuffled_num_list = list(range(0, len(ligands_list) - 1))
    random.shuffle(shuffled_num_list)

    # Convert lig_1 into an RDkit mol
    lig_1_string = ligand1_pair[0]
    lig1_mol = convert_mol_from_smiles(lig_1_string)

    while count < len(ligands_list) - 1:
        rand_num = shuffled_num_list[count]
        mol2_pair = ligands_list[rand_num]

        if mol2_pair[0] == lig_1_string:
            count = count + 1
            continue

        # Convert lig_1 into an RDkit mol
        lig_2_string = mol2_pair[0]
        lig2_mol = convert_mol_from_smiles(lig_2_string)

        if lig2_mol is None:
            count = count + 1
            continue

        # it converts and it is not Ligand1. now lets test for a common
        # substructure
        if test_for_mcs(vars, lig1_mol, lig2_mol) is None:
            count = count + 1
            continue

        # We found a good pair of Ligands
        return mol2_pair

    return False


def convert_mol_from_smiles(smiles_string):
    """
    Test a SMILES string can be converted into an rdkit molecule
    (rdkit.Chem.rdchem.Mol) and be sanitize. This also deprotanates them

    Inputs:
    :param str smiles_string: a single SMILES String

    Returns:
    :returns: rdkit.Chem.rdchem.Mol mol: an rdkit molecule object if it
        properly converts from the SMILE and None
    """

    try:
        mol = Chem.MolFromSmiles(smiles_string, sanitize=False)
    except:
        return None

    mol = MOH.check_sanitization(mol)
    if mol is None:
        return None

    mol = MOH.try_deprotanation(mol)
    if mol is None:
        return False
    return mol


def run_smiles_merge_prescreen(vars, ligands_list, ligand1_pair):
    """
    This function runs a series of functions to find two molecules with a
    sufficient amount of shared common structure (most common structure = MCS)

    These two parent molecules will be derived from a list of molecules from
    the parent generation.

    Inputs:
    :param dict vars: User variables which will govern how the programs runs
    :param list ligands_list: list of all the lignads to chose from
    :param list ligand1_pair: information for the Ligand 1. This info includes
        the name and SMILES string

    Returns:
    :returns: list lig_2_pair: a set of information for a 2nd ligand (Lig2)
        This includes the name and SMILES string this mol is from the
        ligands_list. returns None if a ligand with a sufficient MCS can not be
        found return None.
    """

    ligand_1_string = ligand1_pair[0]

    # check if ligand_1 can be converted to an rdkit mol
    lig_1 = convert_mol_from_smiles(ligand_1_string)
    if lig_1 is False:
        # Ligand1_string failed to be converted to rdkit mol format
        return None

    # GET TWO UNIQUE LIGANDS TO WITH A SHARED SUBSTRUCTURE
    lig_2_pair = find_random_lig2(vars, ligands_list, ligand1_pair)

    if lig_2_pair is False:
        # ligand_1 has no matches
        return None

    return lig_2_pair
