import time

from flask import Flask, request
import os
from rdkit import Chem
import random
import string
import autogrow4.autogrow.operators.operations as operations

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
    json = request.get_json()
    #{game: 1, mol: "asdfas"}
    mol2d = json.get("mol")
    temp_file_name = "".join(random.choice(string.ascii_letters) for i in range(8)).join(".sdf")
    temp_file = open(temp_file_name, "w")
    temp_file.write(mol2d)
    temp_file.close()

    supplier = Chem.SDMolSupplier(temp_file_name)
    mol = supplier[0]
    smiles = Chem.molToSmiles(mol)

    sanitized_smiles = operations.test_source_smiles_convert_update_properties(smiles)

    # if smiles != sanitized_smiles:
    #     return "smiles did not sanitize"

    return "{} {} {} {} {}".format(mol2d, mol.GetNumAtoms(), mol, smiles, sanitized_smiles)
    # Gypsum
    # RDKit
    # RDKit
    # Gypsum
    # RDKit
    pass


@application.route("/execute")
def execute():
    full_path = os.getcwd() + "/"
    sub_folder = "autogrow4-4.0.2/"
    py_file = "RunAutogrow.py"
    params = "sample_sub_scripts/sample_submit_autogrow.json"
    output_file = "current_run_output.txt"
    command = "/conda/bin/python {} -j {}".format(full_path + sub_folder + py_file, full_path + sub_folder + params)
    os.system("{} > {}".format(command, full_path + sub_folder + output_file))
    return "command: {}<br/>json file: {}<br/>output:<br/>{}".format(
        command,
        open("autogrow4-4.0.2/" + params).read(),
        open("autogrow4-4.0.2/" + output_file).read().replace(os.linesep, "<br/>"),
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
