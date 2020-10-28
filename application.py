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
    total_start = time.time()
    output = ""
    start = time.time()
    if not os.path.exists(os.getcwd()+'/smiles_dir'):
        os.makedirs(os.getcwd()+'/smiles_dir')
    if not os.path.exists(os.getcwd()+'/smiles_dir/generation_0'):
        os.makedirs(os.getcwd()+'/smiles_dir/generation_0')
    end = time.time()
    output += "directory making: {} ".format(str(end-start))
    start = time.time()
    json = request.get_json()
    #{game: 1, mol: "asdfas"}
    mol2d = json.get("mol")
    referred_mol_name = "".join(random.choice(string.ascii_letters) for i in range(8))
    temp_file_name = referred_mol_name.join(".sdf")
    temp_file = open(temp_file_name, "w")
    temp_file.write(mol2d)
    temp_file.close()
    end = time.time()
    output += "parse input: {} ".format(str(end-start))
    start = time.time()
    supplier = Chem.SDMolSupplier(temp_file_name)
    mol = supplier[0]
    smiles = Chem.MolToSmiles(mol)

    sanitized_smiles = operations.test_source_smiles_convert_update_properties(smiles)
    end = time.time()
    output += "set up and sanitize: {} ".format(str(end-start))
    start = time.time()

    # if smiles != sanitized_smiles:
    #     return "smiles did not sanitize"

    # vars (in this method) will contain a lot of stuff we don't care about.  That is okay because
    # the method calls will only use what they need

    vars = user_vars.define_defaults()
    end = time.time()
    output += "define defaults: {} ".format(str(end-start))
    start = time.time()
    vars = user_vars.multiprocess_handling(vars)
    end = time.time()
    output += "multiprocess handling: {} ".format(str(end-start))
    vars["max_variants_per_compound"] = 1
    start = time.time()
    smiles_list = [(sanitized_smiles, referred_mol_name)]
    smiles_to_convert_file, new_gen_folder_path = operations.save_generation_smi(
        os.getcwd()+"/smiles_dir/",
        0,
        smiles_list,
        None,
        )
    end = time.time()
    output += "save smi: {} ".format(str(end-start))
    start = time.time()
    conversion_to_3d.convert_to_3d(vars, smiles_to_convert_file, new_gen_folder_path)
    end = time.time()
    output += "convert to 3d: {} ".format(str(end-start))
    start = time.time()
    file_in_3d_folder = os.system("dir "+new_gen_folder_path+"/3D_SDFs/ > 3d_folder_contents.txt")

    rdkit_mol_sdf = Chem.MolToMolBlock(mol)
    end = time.time()
    output += "mol to mol block: {} ".format(str(end-start))
    start = time.time()

    folder_contents = open("3d_folder_contents.txt").read()
    print(folder_contents)
    total_end = time.time()
    output += "total process time: {}".format(total_end-total_start)
    # return output
    threed_sdf = open(new_gen_folder_path+"/3D_SDFs/{}__input1.sdf".format(referred_mol_name)).read()
    os.remove(new_gen_folder_path+"/3D_SDFs/{}__input1.sdf".format(referred_mol_name))
    mw = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
    donors = Chem.rdMolDescriptors.CalcNumHBD(mol)
    acceptors = Chem.rdMolDescriptors.CalcNumHBA(mol)
    rotatable = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
    rings = Chem.rdMolDescriptors.CalcNumRings(mol)
    tpsa = Chem.rdMolDescriptors.CalcTPSA(mol)
    logP = Chem.Crippen.MolLogP(mol)
    heteroAtoms = Chem.rdMolDescriptors.CalcNumHeteroatoms(mol)
    heavyAtoms = mol.getNumHeavyAtoms()
    complexity = Chem.GraphDescriptors.BertzCT(mol)
    return "{}\n{}\n{}\n{}\n{}\n{}\n{}".format(mol2d, rdkit_mol_sdf, mol, smiles, sanitized_smiles,
                                               threed_sdf,
                                               output
                                               )
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
