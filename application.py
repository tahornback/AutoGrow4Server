import time

from flask import Flask, request
import os
import autogrow4.autogrow.operators.operations as operations
import autogrow4.autogrow.operators.convert_files.conversion_to_3d as conversion_to_3d
import autogrow4.autogrow.user_vars as user_vars
operations.test()
from rdkit import Chem
import random
import string

application = Flask(__name__)


@application.route("/")
def hello():
    return "Hello World!"


@application.route("/env")
def env():
    os.system("which python > whichOutput.txt")
    os.system("which gunicorn >> whichOutput.txt")
    os.system("env | more > envOutput.txt")
    return "which: {}<br/>env|more: {}".format(open("whichOutput.txt").read(), open("envOutput.txt").read())


@application.route("/updateProperties", methods=["POST"])
def updateProperties():
    # Input: 2Dmol sdf as string
    # Return: Smiles file, new 2dmol, 3dmol, descriptors
    if not os.path.exists(os.getcwd()+'/smiles_dir'):
        os.makedirs(os.getcwd()+'/smiles_dir')
    if not os.path.exists(os.getcwd()+'/smiles_dir/generation_0'):
        os.makedirs(os.getcwd()+'/smiles_dir/generation_0')
    json = request.get_json()
    #{game: 1, mol: "asdfas"}
    mol2d = json.get("mol")
    temp_file_name = "".join(random.choice(string.ascii_letters) for i in range(8)).join(".sdf")
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
    smiles_list = [sanitized_smiles]
    smiles_to_convert_file, new_gen_folder_path = operations.save_generation_smi(
        os.getcwd()+"/smiles_dir/",
        0,
        smiles_list,
        None,
    )
    conversion_to_3d.convert_to_3d(vars, smiles_to_convert_file, new_gen_folder_path)
    file_in_3d_folder = os.system("dir "+new_gen_folder_path+"/3D_SDFs/ > 3d_folder_contents.txt")
    return "{}\n{}\n{}\n{}\n{}\n{}\n{}".format(mol2d, mol.GetNumAtoms(), mol, smiles, sanitized_smiles, open("3d_folder_contents.txt").read(),
                                               open(new_gen_folder_path+"/3D_SDFs/C__input1.sdf").read()
                                               )
    # Gypsum
    # RDKit
    # RDKit
    # Gypsum
    # RDKit
    pass


@application.route("/execute")
def execute():
    full_path = os.getcwd() + "/"
    sub_folder = "autogrow4/"
    py_file = "RunAutogrow.py"
    params = "sample_sub_scripts/sample_submit_autogrow.json"
    output_file = "current_run_output.txt"
    command = "/conda/bin/python {} -j {}".format(full_path + sub_folder + py_file, full_path + sub_folder + params)
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
        times_to_execute, str(time_taken), str(time_taken / int(times_to_execute)), printout
    )


if __name__ == "__main__":
    application.debug = True
    application.run()
