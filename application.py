import time

from flask import Flask, request
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


@application.route("/execute-time-trial")
def executeTimeTrial():
    start_time = time.time()
    times_to_execute = request.args.get("times")
    printout = ""
    for i in range(times_to_execute):
        printout += execute() + "<br/>"
    end_time = time.time()
    time_taken = end_time - start_time
    return "Ten executions finished in {} ms, averaging {} per run<br/><br/><br/>{}".format(
        time_taken, time_taken // times_to_execute, printout
    )


if __name__ == "__main__":
    application.debug = True
    application.run()
