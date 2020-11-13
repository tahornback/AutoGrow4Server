"""user_vars
This should contain the functions for defining input variables.
Both the default variables and the user input variables.
This should also validate them.
"""

import copy
import json
import os
import platform
import sys
from shutil import copyfile

def multiprocess_handling(vars):
    """
    This function handles the multiprocessing functions. It establishes a Paralellizer object
    and adds it to the vars dictionary.
    Inputs:
    :param dict vars: dict of user variables which will govern how the programs runs
    Returns:
    :returns: dict vars: dict of user variables which will govern how the programs runs
    """

    # Handle Serial overriding number_of_processors
    # serial fixes it to 1 processor
    if vars["multithread_mode"].lower() == "serial":
        vars["multithread_mode"] = "serial"
        if vars["number_of_processors"] != 1:
            print(
                "Because --multithread_mode was set to serial, "
                + "this will be run on a single processor."
            )
        vars["number_of_processors"] = 1

    # Handle mpi errors if mpi4py isn't installed
    if vars["multithread_mode"].lower() == "mpi":
        vars["multithread_mode"] = "mpi"
        try:
            import mpi4py
        except:
            printout = "mpi4py not installed but --multithread_mode is set to"
            printout = printout + " mpi. \n Either install mpi4py or switch "
            printout = printout + "multithread_mode to multithreading or serial"
            raise ImportError(printout)

        try:
            import func_timeout
            from func_timeout import func_timeout, FunctionTimedOut
        except:
            printout = "func_timeout not installed but --multithread_mode is "
            printout = printout + "set to mpi. \n Either install func_timeout "
            printout = printout + "or switch multithread_mode to"
            printout = printout + " multithreading or serial"
            raise ImportError(printout)

    # # # launch mpi workers
    if vars["multithread_mode"] == "mpi":
        # Avoid EOF error
        try:
            from autogrow4.autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Parallelizer import (
                Parallelizer,
            )
        except Exception as e:
            from autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Parallelizer import (
                Parallelizer,
            )
        vars["parallelizer"] = Parallelizer(
            vars["multithread_mode"], vars["number_of_processors"]
        )

        if vars["parallelizer"] is None:
            printout = "EOF ERRORS FAILED TO CREATE A PARALLIZER OBJECT"
            print(printout)
            raise Exception(printout)

    else:
        # Lower level mpi (ie making a new Parallelizer within an mpi)
        #   has problems with importing the MPI environment and mpi4py
        #   So we will flag it to skip the MPI mode and just go to multithread/serial
        # This is a saftey precaution
        try:
            from autogrow4.autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Parallelizer import Parallelizer
        except Exception as e:
            from autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Parallelizer import Parallelizer
        vars["parallelizer"] = Parallelizer(
            vars["multithread_mode"], vars["number_of_processors"], True
        )

    return vars

def determine_bash_timeout_vs_gtimeout():
    """
    This function tests whether we should use the BASH command "timeout" (for linux)
     or the coreutils function "gtimeout" for MacOS which can be obtained
     through homebrew

    Returns:
    :returns: str timeout_option: A string either "timeout" or "gtimeout" describing
     whether the bash terminal is able to use the bash function timeout or gtimeout
    """

    if sys.platform.lower() in ["linux", "linux2"]:
        # Should be true and default installed in all Linux machines
        return "timeout"

    command = 'timeout 1 echo " "'
    # Running the os.system command for command will return 0,1, or 32512
    # 0 means that the timeout function works (most likely this is a linux os)
    # 32512 means that the timeout function DOES NOT Work (most likely this is MacOS)

    try:  # timeout or gtimeout
        timeout_result = os.system("g" + command)
    except:
        raise Exception(
            "Something is very wrong. This OS may not be supported \
            by Autogrow or you may need to execute through Bash."
        )
    if timeout_result == 0:
        timeout_option = "gtimeout"
        return timeout_option
    print("gtimeout failed to run, we will check timeout")

    try:  # timeout or gtimeout
        timeout_result = os.system(command)
    except:
        raise Exception(
            "Something is very wrong. This OS may not be supported by \
            Autogrow or you may need to execute through Bash."
        )

    if timeout_result == 0:
        timeout_option = "timeout"
        return timeout_option

    printout = "Need to install GNU tools for Bash to work. \n"
    printout = (
        printout
        + "This is essential to use Bash Timeout function in Autogrow. \n"
    )
    printout = printout + "\t This will require 1st installing homebrew. \n"
    printout = printout + "\t\t Instructions found at: https://brew.sh/ \n"
    printout = printout + "\t Once brew is installed, please run:"
    printout = printout + " sudo brew install coreutils \n\n"
    print(printout)
    raise Exception(printout)

def define_defaults():
    """
    Sets the command-line parameters to their default values.

    Returns:
    :returns: dict vars: a dictionary of all default variables
    """

    vars = {}

    # where we are currently (absolute filepath from route)
    # used for relative pathings
    script_dir = os.path.dirname(os.path.realpath(__file__))

    # Some variables which can be manually replaced but defaults
    # point to prepackaged locations.
    ## Neural Network executable for scoring binding
    vars["nn1_script"] = os.path.join(
        script_dir, "docking", "scoring", "nn_score_exe", "nnscore1", "NNScore.py"
    )
    # Example: vars['nn1_script'] =
    #    "/PATH/autogrow4/autogrow/docking/scoring/nn_score_exe/nnscore1/NNScore.py"

    vars["nn2_script"] = os.path.join(
        script_dir, "docking", "scoring", "nn_score_exe", "nnscore2", "NNScore2.py"
    )
    # Example: vars['nn2_script'] =
    #    "/PATH/autogrow4/autogrow/docking/scoring/nnscore2/NNScore2.py"

    #### OPTIONAL FILE-LOCATION VARIABLES ####
    # (RECOMMEND SETTING TO "" SO AUTOGROW CAN AUTOLOCATE THESE FILES)#

    # PARSER.add_argument('--conversion_choice', choices
    #    = ["MGLTools","obabel"], default="MGLTools",
    vars["conversion_choice"] = "MGLToolsConversion"
    vars["obabel_path"] = "obabel"
    vars["custom_conversion_script"] = ""
    # vars['prepare_ligand4.py'] =
    #   "/PATH/MGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
    vars["prepare_ligand4.py"] = ""
    # vars['prepare_receptor4.py'] =
    #   "/PATH/MGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
    vars["prepare_receptor4.py"] = ""
    # vars['mgl_python'] = "/PATH/MGLTools-1.5.4/bin/pythonsh"
    vars["mgl_python"] = ""

    # Crossover function
    vars["start_a_new_run"] = False
    vars["max_time_mcs_prescreen"] = 1
    vars["max_time_mcs_thorough"] = 1
    vars["min_atom_match_mcs"] = 4
    vars["protanate_step"] = False

    # Mutation Settings
    vars["rxn_library"] = "click_chem_rxns"
    vars["rxn_library_file"] = ""
    vars["function_group_library"] = ""
    vars["complementary_mol_directory"] = ""

    # processors
    vars["number_of_processors"] = 1
    vars["multithread_mode"] = "multithreading"

    # Genetic Algorithm Components
    vars["selector_choice"] = "Roulette_Selector"
    vars["tourn_size"] = 0.1

    # Seeding next gen and diversity
    vars["top_mols_to_seed_next_generation_first_generation"] = 10
    vars["top_mols_to_seed_next_generation"] = 10
    vars["diversity_mols_to_seed_first_generation"] = 10
    vars["diversity_seed_depreciation_per_gen"] = 2

    # Populations settings
    vars["filter_source_compounds"] = True
    vars["use_docked_source_compounds"] = True
    vars["num_generations"] = 10
    vars["number_of_crossovers_first_generation"] = 10
    vars["number_of_mutants_first_generation"] = 10
    vars["number_of_crossovers"] = 10
    vars["number_of_mutants"] = 10
    vars["number_elitism_advance_from_previous_gen"] = 10
    vars["number_elitism_advance_from_previous_gen_first_generation"] = 10
    vars["redock_elite_from_previous_gen"] = False

    # Filters
    vars["LipinskiStrictFilter"] = False
    vars["LipinskiLenientFilter"] = False
    vars["GhoseFilter"] = False
    vars["GhoseModifiedFilter"] = False
    vars["MozziconacciFilter"] = False
    vars["VandeWaterbeemdFilter"] = False
    vars["PAINSFilter"] = False
    vars["NIHFilter"] = False
    vars["BRENKFilter"] = False
    vars["No_Filters"] = False
    vars["alternative_filter"] = None

    # docking
    vars["dock_choice"] = "QuickVina2Docking"
    vars["docking_executable"] = None
    vars["docking_exhaustiveness"] = None
    vars["docking_num_modes"] = None
    vars["docking_timeout_limit"] = 120
    vars["custom_docking_script"] = ""

    # scoring
    vars["scoring_choice"] = "VINA"
    vars["rescore_lig_efficiency"] = False
    vars["custom_scoring_script"] = ""

    # gypsum # max variance is the number of conformers made per ligand
    vars["max_variants_per_compound"] = 3
    vars["gypsum_thoroughness"] = 3
    vars["min_ph"] = 6.4
    vars["max_ph"] = 8.4
    vars["pka_precision"] = 1.0
    vars["gypsum_timeout_limit"] = 10

    # Other vars
    vars["debug_mode"] = False
    vars["reduce_files_sizes"] = False
    vars["generate_plot"] = True
    # Check Bash Timeout function (There's a difference between MacOS and linux)
    # Linux uses timeout while MacOS uses gtimeout
    timeout_option = determine_bash_timeout_vs_gtimeout()
    if timeout_option in  ["timeout", "gtimeout"]:
        vars["timeout_vs_gtimeout"] = timeout_option
    else:
        raise Exception(
            "Something is very wrong. This OS may not be supported by \
             Autogrow or you may need to execute through Bash."
        )

    return vars
