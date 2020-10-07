from flask import Flask
import os

application = Flask(__name__)


@application.route("/")
def hello():
    return "Hello World!"


@application.route("/execute")
def execute():
    py_file = "RunAutogrow.py"
    params = "./sample_sub_scripts{}sample_submit_autogrow.json".format(os.sep)
    output_file = "current_run_output.txt"
    os.chdir("autogrow4-4.0.2")
    command = "/conda/bin/python {} -j {}".format(py_file, params)
    os.system("{} > {}".format(command, output_file))
    os.chdir("..")
    return "command: {}<br/>json file: {}<br/>output:<br/>{}".format(
        command,
        open("autogrow4-4.0.2/" + params).read(),
        open("autogrow4-4.0.2/" + output_file).read().replace(os.linesep, "<br/>"),
    )


if __name__ == "__main__":
    application.debug = True
    application.run()
