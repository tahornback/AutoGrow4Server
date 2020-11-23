import base64
import os
import random
import shutil
import string

import autogrow4.autogrow.docking.execute_docking as execute_docking
import autogrow4.autogrow.operators.convert_files.conversion_to_3d as conversion_to_3d
import autogrow4.autogrow.operators.operations as operations
import autogrow4.autogrow.user_vars as user_vars
from flask import Flask, request, jsonify
from rdkit import Chem
import rdkit.Chem.Crippen as Crippen
import rdkit.Chem.GraphDescriptors as GraphDescriptors

application = Flask(__name__)


@application.route("/")
def hello():
    return "Hello World!"


@application.route("/env")
def env():
    os.system("which python > whichOutput.txt")
    os.system("which gunicorn >> whichOutput.txt")
    os.system("env | more > envOutput.txt")
    return "which: {}<br/>env|more: {}".format(
        open("whichOutput.txt").read(), open("envOutput.txt").read()
    )


@application.route("/updateProperties", methods=["POST"])
def updateProperties():
    # Input: 2Dmol sdf as string
    # Return: Smiles file, new 2dmol, 3dmol, descriptors
    # return obj {
    #     smiles: { STRING },
    #     mol2D: { STRING },
    #     mol3D: { STRING },
    #     mw: { INTEGER },
    #     formula: { STRING },
    #     donors: { INTEGER },
    #     acceptors: { INTEGER },
    #     rotatable: { INTEGER },
    #     rings: { INTEGER },
    #     tpsa: { DECIMAL },
    #     logP: { DECIMAL },
    #     heteroAtoms: { INTEGER },
    #     heavyAtoms: { INTEGER },
    #     complexity: { DECIMAL }
    # }

    # Set up some folders in case they don't exist yet
    CWD = os.getcwd()
    if not os.path.exists(CWD + "/smiles_dir"):
        os.makedirs(CWD + "/smiles_dir")
    if not os.path.exists(CWD + "/smiles_dir/generation_0"):
        os.makedirs(CWD + "/smiles_dir/generation_0")

    # Grab JSON from request
    json = request.get_json()

    # Get the 2dsdf from the request, save it, and turn it into an rdkit mol
    # referred_mol_name is just to make sure that multiple requests do not conflict with each other
    mol2d = json.get("mol")
    referred_mol_name = "".join(random.choice(string.ascii_letters) for i in range(8))
    temp_file_name = referred_mol_name + ".sdf"
    temp_file = open(temp_file_name, "w")
    temp_file.write(mol2d)
    temp_file.close()
    supplier = Chem.SDMolSupplier(temp_file_name)
    mol = supplier[0]
    smiles = Chem.MolToSmiles(mol)

    # This will return the smile if sanitized correctly or some other stuff we don't care about if it fails
    sanitized_smiles = operations.test_source_smiles_convert_update_properties(smiles)

    if smiles != sanitized_smiles:
        return "smiles did not sanitize"

    # vars (in this method) will contain a lot of stuff we don't care about.  That is okay because
    # the method calls will only use what they need

    vars = user_vars.define_defaults()
    vars = user_vars.multiprocess_handling(vars)
    vars["max_variants_per_compound"] = 1

    # Set up for conversion to 3d and convert
    smiles_list = [(sanitized_smiles, referred_mol_name)]
    smiles_to_convert_file, new_gen_folder_path = operations.save_generation_smi(
        CWD + "/smiles_dir/", 0, smiles_list, None,
    )
    conversion_to_3d.convert_to_3d(vars, smiles_to_convert_file, new_gen_folder_path)

    # the rdkit 2dsdf of the mol, for consistent drawing of molecule
    rdkit_mol_sdf = Chem.MolToMolBlock(mol)

    # return output
    threed_sdf = open(
        new_gen_folder_path + "/3D_SDFs/{}__input1.sdf".format(referred_mol_name)
    ).read()
    threed_sdf = threed_sdf.split(">")[0] + "$$$$"
    os.remove(new_gen_folder_path + "/3D_SDFs/{}__input1.sdf".format(referred_mol_name))

    dictionary = {
        "mw": Chem.rdMolDescriptors.CalcExactMolWt(mol),
        "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
        "donors": Chem.rdMolDescriptors.CalcNumHBD(mol),
        "acceptors": Chem.rdMolDescriptors.CalcNumHBA(mol),
        "rotatable": Chem.rdMolDescriptors.CalcNumRotatableBonds(mol),
        "rings": Chem.rdMolDescriptors.CalcNumRings(mol),
        "tpsa": Chem.rdMolDescriptors.CalcTPSA(mol),
        "logP": Crippen.MolLogP(mol),
        "heteroAtoms": Chem.rdMolDescriptors.CalcNumHeteroatoms(mol),
        "heavyAtoms": mol.GetNumHeavyAtoms(),
        "complexity": GraphDescriptors.BertzCT(mol),
        "smiles": sanitized_smiles,
        "mol2D": rdkit_mol_sdf,
        "mol3D": threed_sdf,
    }
    return jsonify(dictionary)


@application.route("/submitDocking", methods=["POST"])
def dock():
    # {
    #     "center_x": "-70.76",
    #     "center_y": "21.82",
    #     "center_z": "28.33",
    #     "size_x": "25.0",
    #     "size_y": "16.0",
    #     "size_z": "25.0",
    #     "docking_exhaustiveness": "1",
    #     "inputFile":[
    #         {
    #             "name":"`receptor.pdbqt`",
    #             "contents":"encodedReceptor"
    #         },
    #         {
    #             "name":"`ligand.pdbqt`",
    #             "contents":"encodedLigand"
    #         }
    #     ]
    # }

    # File/directory setup
    CWD = os.getcwd()
    suffix = "".join(random.choice(string.ascii_letters) for i in range(8))
    temp_folder = "/Output_{}".format(suffix)
    os.mkdir(CWD + "/" + temp_folder)

    # vars setup for autogrow methods
    vars = user_vars.define_defaults()
    vars = user_vars.multiprocess_handling(vars)
    vars["output_directory"] = CWD + temp_folder
    vars["root_output_folder"] = CWD + temp_folder
    current_generation_dir = (
            vars["output_directory"] + "/" + "generation_{}{}".format(0, os.sep)
    )
    vars["num_generations"] = 0
    vars["number_of_processors"] = -1
    vars["docking_num_modes"] = 1

    os.mkdir(current_generation_dir)
    os.mkdir(current_generation_dir + "PDBs/")

    # chunk out json stuff
    json = request.get_json()
    inputFile = json.get("inputFile")
    encodedReceptor = inputFile[0].get("contents")
    encodedLigand = inputFile[1].get("contents")

    # Write out files
    receptor_temp_file_name = "receptor_{}".format(suffix + ".pdbqt")
    receptor_temp_file = open(receptor_temp_file_name, "w")
    receptor_file_data = base64.b64decode(encodedReceptor).decode("utf-8")
    receptor_temp_file.write(receptor_file_data)
    receptor_temp_file.close()

    ligand_temp_file_name = "ligand_{}".format(suffix + ".pdbqt")
    ligand_temp_file_path = current_generation_dir + "PDBs/" + ligand_temp_file_name
    ligand_temp_file = open(
        ligand_temp_file_path, "w"
    )
    ligand_file_data = base64.b64decode(encodedLigand).decode("utf-8")
    ligand_temp_file.write(ligand_file_data)
    ligand_temp_file.close()

    # add opts from sample_submit_autogrow.json to vars, like mgltools_directory etc.
    # some of these will come from json
    vars["mgltools_directory"] = "/mgltools_x86_64Linux2_1.5.6/"
    vars["center_x"] = float(json.get("center_x"))
    vars["center_y"] = float(json.get("center_y"))
    vars["center_z"] = float(json.get("center_z"))
    vars["size_x"] = float(json.get("size_x"))
    vars["size_y"] = float(json.get("size_y"))
    vars["size_z"] = float(json.get("size_z"))
    vars["docking_exhaustiveness"] = int(json.get("docking_exhaustiveness"))
    vars["filename_of_receptor"] = (
            CWD + "/" + receptor_temp_file_name[:-2]
    )  # Cut off qt part, comparison later will check for `filename`+"qt" so we need to cut off the end
    vars["source_compound_file"] = CWD + "/" + ligand_temp_file_name

    execute_docking.run_docking_common(vars, 0, current_generation_dir, None)

    try:  # If everything was successful in docking, .vina file will be returned
        file_contents = open(
            current_generation_dir
            + "PDBs/"
            + ligand_temp_file_name
            + ".vina",
            "r",
            ).read()
    except:  # If something failed during docking, the output file will be returned, hopefully with good feedback
        file_contents = open(
            current_generation_dir
            + "PDBs/"
            + ligand_temp_file_name
            + "_docking_output.txt",
            "r",
                ).read()

    # Clean up files
    shutil.rmtree(vars["output_directory"])
    os.remove(CWD + "/" + receptor_temp_file_name)

    return file_contents


if __name__ == "__main__":
    application.debug = True
    application.run()
