import time

from flask import Flask, request, jsonify
import os
import autogrow4.autogrow.operators.operations as operations
import autogrow4.autogrow.operators.convert_files.conversion_to_3d as conversion_to_3d
import autogrow4.autogrow.user_vars as user_vars
import autogrow4.autogrow.autogrow_main_execute as autogrow_main_execute
import autogrow4.autogrow.docking.execute_docking as execute_docking

operations.test()
from rdkit import Chem
import random
import string
import tempfile
import base64

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
    if not os.path.exists(os.getcwd() + "/smiles_dir"):
        os.makedirs(os.getcwd() + "/smiles_dir")
    if not os.path.exists(os.getcwd() + "/smiles_dir/generation_0"):
        os.makedirs(os.getcwd() + "/smiles_dir/generation_0")
    json = request.get_json()
    mol2d = json.get("mol")
    referred_mol_name = "".join(random.choice(string.ascii_letters) for i in range(8))
    temp_file_name = referred_mol_name + ".sdf"
    temp_file = open(temp_file_name, "w")
    temp_file.write(mol2d)
    temp_file.close()
    supplier = Chem.SDMolSupplier(temp_file_name)
    mol = supplier[0]
    smiles = Chem.MolToSmiles(mol)

    sanitized_smiles = operations.test_source_smiles_convert_update_properties(smiles)

    # if smiles != sanitized_smiles:
    #     return "smiles did not sanitize"

    # vars (in this method) will contain a lot of stuff we don't care about.  That is okay because
    # the method calls will only use what they need

    vars = user_vars.define_defaults()
    vars = user_vars.multiprocess_handling(vars)
    vars["max_variants_per_compound"] = 1
    smiles_list = [(sanitized_smiles, referred_mol_name)]
    smiles_to_convert_file, new_gen_folder_path = operations.save_generation_smi(
        os.getcwd() + "/smiles_dir/", 0, smiles_list, None,
    )
    conversion_to_3d.convert_to_3d(vars, smiles_to_convert_file, new_gen_folder_path)
    file_in_3d_folder = os.system(
        "dir " + new_gen_folder_path + "/3D_SDFs/ > 3d_folder_contents.txt"
    )

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
        "logP": Chem.Crippen.MolLogP(mol),
        "heteroAtoms": Chem.rdMolDescriptors.CalcNumHeteroatoms(mol),
        "heavyAtoms": mol.GetNumHeavyAtoms(),
        "complexity": Chem.GraphDescriptors.BertzCT(mol),
        "smiles": sanitized_smiles,
        "mol2D": rdkit_mol_sdf,
        "mol3D": threed_sdf,
    }
    return jsonify(dictionary)


@application.route("/submitDocking", methods=["POST"])
def dock():
    # {
    #     "argList":"`--receptor receptor.pdbqt --ligand ligand.pdbqt --energy_range ${energy_range} --center_x ${receptorX} --center_y ${receptorY} --center_z ${receptorZ} --size_x ${sizeX} --size_y ${sizeY} --size_z ${sizeZ} --out output`",
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
    suffix = "".join(random.choice(string.ascii_letters) for i in range(8))
    temp_folder = "/Output_{}".format(suffix)
    os.mkdir(os.getcwd() + "/" + temp_folder)
    receptor_temp_file_name = "receptor_{}".format(suffix + ".pdbqt")
    receptor_temp_file = open(receptor_temp_file_name, "w")
    ligand_temp_file_name = "ligand_{}".format(suffix + ".pdbqt")
    ligand_temp_file = open(ligand_temp_file_name, "w")

    json = request.get_json()
    inputFile = json.get("inputFile")
    encodedReceptor = inputFile[0].get("contents")
    encodedLigand = inputFile[1].get("contents")

    receptor_file_data = base64.b64decode(encodedReceptor).decode("utf-8")
    ligand_file_data = base64.b64decode(encodedLigand).decode("utf-8")

    # print(receptor_file_data)
    # print(ligand_file_data)

    receptor_temp_file.write(receptor_file_data)
    ligand_temp_file.write(ligand_file_data)
    receptor_temp_file.close()
    ligand_temp_file.close()

    args = json.get("argList").split(" ")
    vars = user_vars.define_defaults()
    vars = user_vars.multiprocess_handling(vars)
    # add opts from sample_submit_autogrow.json to vars, like mgltools_directory etc.
    # some of these will come from json.get("argList")
    vars["mgltools_directory"] = "/mgltools_x86_64Linux2_1.5.6/"
    vars["center_x"] = float(args[args.index("--center_x") + 1])
    vars["center_y"] = float(args[args.index("--center_y") + 1])
    vars["center_z"] = float(args[args.index("--center_z") + 1])
    vars["size_x"] = float(args[args.index("--size_x") + 1])
    vars["size_y"] = float(args[args.index("--size_y") + 1])
    vars["size_z"] = float(args[args.index("--size_z") + 1])
    vars["num_generations"] = 0
    vars["output_directory"] = os.getcwd() + temp_folder
    vars["root_output_folder"] = os.getcwd() + temp_folder
    vars["filename_of_receptor"] = (
        os.getcwd() + "/" + receptor_temp_file_name[:-2]
    )  # Cut off qt part
    vars["source_compound_file"] = os.getcwd() + "/" + ligand_temp_file_name
    vars["docking_exhaustiveness"] = 1
    vars["number_of_processors"] = -1
    print(vars)
    # what is energy range??

    current_generation_dir = (
        vars["output_directory"] + "/" + "generation_{}{}".format(0, os.sep)
    )
    os.mkdir(current_generation_dir)
    os.mkdir(current_generation_dir + "PDBs/")
    ligand_temp_file_name = "ligand_{}".format(suffix + ".pdbqt")
    ligand_temp_file = open(
        current_generation_dir + "PDBs/" + ligand_temp_file_name, "w"
    )
    ligand_file_data = base64.b64decode(encodedLigand).decode("utf-8")
    ligand_temp_file.write(ligand_file_data)
    ligand_temp_file.close()

    # Should probably run user_vars.check_for_required_inputs(vars)
    # just to make sure we got everything
    # user_vars.check_for_required_inputs(vars)

    # autogrow_main_execute.main_execute(vars)
    execute_docking.run_docking_common(vars, 0, current_generation_dir, None)

    # os.system("ls -lR > filetree.txt")
    return "{}".format(
        open(
            current_generation_dir
            + "PDBs/"
            + ligand_temp_file_name
            + ".vina",
            "r",
        ).read()
    )


@application.route("/execute")
def execute():
    full_path = os.getcwd() + "/"
    sub_folder = "autogrow4/"
    py_file = "RunAutogrow.py"
    params = "sample_sub_scripts/sample_submit_autogrow.json"
    output_file = "current_run_output.txt"
    command = "/conda/bin/python {} -j {}".format(
        full_path + sub_folder + py_file, full_path + sub_folder + params
    )
    os.system("{} > {}".format(command, full_path + sub_folder + output_file))
    return "command: {}<br/>json file: {}<br/>output:<br/>{}".format(
        command,
        open(sub_folder + params).read(),
        open(sub_folder + output_file).read().replace(os.linesep, "<br/>"),
    )


@application.route("/execute-time-trial")
def executeTimeTrial():
    start_time = time.time()
    times_to_execute = request.args.get("times")
    printout = ""
    for i in range(int(times_to_execute)):
        printout += execute() + "<br/>"
    end_time = time.time()
    time_taken = end_time - start_time
    return "{} executions finished in {} seconds, averaging {} seconds per run<br/><br/><br/>{}".format(
        times_to_execute,
        str(time_taken),
        str(time_taken / int(times_to_execute)),
        printout,
    )


if __name__ == "__main__":
    application.debug = True
    application.run()
