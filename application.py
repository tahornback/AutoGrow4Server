import time

from flask import Flask, request
import os
# from rdkit import Chem
import random
import string

application = Flask(__name__)


@application.route("/")
def hello():
    return "Hello World!"


@application.route("/env")
def env():
    os.system("which python > whichOutput.txt")
    os.system("env | more > envOutput.txt")
    return "which: {}<br/>env|more: {}".format(open("whichOutput.txt").read(), open("envOutput.txt"))


# @application.route("/updateProperties")
# def updateProperties():
#     # Input: 2Dmol sdf as string
#     # Return: Smiles file, new 2dmol, 3dmol, descriptors
#     json = request.get_json()
#     #{game: 1, mol: "asdfas"}
#     mol2d = json.get("mol")
#     temp_file_name = "".join(random.choice(string.ascii_letters) for i in range(8)).join(".sdf")
#     temp_file = open(temp_file_name, "rw")
#     temp_file.write(mol2d)
#     temp_file.close()
#
#     mol = Chem.SDMolSupplier(temp_file_name)
#     return mol.GetNumAtoms()
#     # Gypsum
#     # RDKit
#     # RDKit
#     # Gypsum
#     # RDKit
#     pass


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
