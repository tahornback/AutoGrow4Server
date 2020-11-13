"""
Top level for running AutoGrow.
Runs all population generation (operations) and docking.
Runs plotting at end.
"""

import glob
import os
import shutil
import sys

import autogrow4.autogrow.docking.concatenate_files as concatenate_files
import autogrow4.autogrow.docking.execute_docking as DockingClass
import autogrow4.autogrow.operators.operations as operations


def main_execute(vars):
    """
    This function takes the user variables and runs Autogrow

    Inputs:
    :param dict vars: dict of user variables which will govern how the
        programs runs
    """

    # Unpack necessary variables
    # output_directory is the root output folder for the run
    output_directory = vars["output_directory"]
    num_gens_to_make = vars["num_generations"]

    # Determine what was the last completed generation in the Run directory
    starting_generation_num = current_generation_number = 0

    if starting_generation_num > num_gens_to_make:
        print("This simulation has already been completed to the user defined number \
                of generations. Please check your user variables.")
        raise Exception("This simulation has already been completed to the user defined number \
                of generations. Please check your user variables.")

    # This is the main loop which will control and execute all commands This
    # is broken into 3 main sections:
    # 1)  operations which populating the new generation with ligands which
    #     both pass the userdefined filter and convert from 1D smiles to 3D
    #     PDB
    # 2)  Docking which handles converting from PDBs to Docking specific
    #     formats and running the actual Docking simulations
    # 3)  Ranking the generation based on the Docking scores
    sys.stdout.flush()

    # Get directory for smi to go
    current_generation_dir = vars["output_directory"] + "generation_{}{}".format(current_generation_number, os.sep)
    print(current_generation_dir)
    sys.stdout.flush()

    already_docked, smile_file_new_gen, new_gen_ligands_list = operations.populate_generation_zero(vars, generation_num=0)
    sys.stdout.flush()

    # Run file conversions of PDB to docking specific file type
    # and Begin Docking unweighted_ranked_smile_file is the file
    # name where the unweighted ranked but score .smi file resides
    unweighted_ranked_smile_file = DockingClass.run_docking_common(
        vars, current_generation_number,
        current_generation_dir, smile_file_new_gen)

    # Delete all temporary files; Skip if in Debugging Mode
    if vars["debug_mode"] is False:
        print("Deleting temporary files and directories")
        files_to_del = []
        folders_to_del = ["{}{}3D_SDFs{}".format(current_generation_dir, os.sep, os.sep), "{}{}3D_SDFs{}log{}".format(current_generation_dir, os.sep, os.sep, os.sep), "{}{}gypsum_submission_files{}".format(current_generation_dir, os.sep, os.sep)]
        for folder in folders_to_del:
            if os.path.exists(folder) is False:
                continue
            files_to_del.extend(glob.glob(folder+"*"))

        job_input = tuple([tuple([x]) for x in files_to_del if os.path.isfile(x) is True])
        vars["parallelizer"].run(job_input, delete_temporary_files_and_folders)
        # Delete Folders in an ordered manor incase folders are nested
        for i in range(0, len(folders_to_del)):
            delete_temporary_files_and_folders(folders_to_del[i])

    sys.stdout.flush()
    if vars["reduce_files_sizes"] is True:
        # Reduce the files in the PDBs folder to a single compiled file.
        # This reduces the data size And makes it easier to transfer the
        # data
        pdbs_folder = "{}{}PDBs{}".format(current_generation_dir, os.sep, os.sep)
        if os.path.exists(pdbs_folder) is True:
            concatenate_files.run_concatenation(vars["parallelizer"], pdbs_folder)
        else:
            print("\nNo PDB folder to concatenate and compress. This is likely generation 0 seeded with a Ranked .smi file.\n")
    print("")
    print("Finished generation ", current_generation_number)

    sys.stdout.flush()
#

def delete_temporary_files_and_folders(file_or_folder):
    """
    This deletes all temporary files.

    Inputs:
    :param str file_or_folder: the file or folder to delete

    """
    if os.path.exists(file_or_folder) is True:
        if os.path.isdir(file_or_folder) is True:
            try:
                shutil.rmtree(file_or_folder)
            except:
                pass
        else:
            try:
                os.remove(file_or_folder)
            except:
                pass

        # If it failed to delete try via bash command
        if os.path.exists(file_or_folder) is True:
            command = "rm -rf {}".format(file_or_folder)
            try:
                os.system(command)
            except:
                pass
    else:
        pass
#
