import time

from flask import Flask, request
import os

application = Flask(__name__)


@application.route("/")
def hello():
    return "Hello World!"


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
