from flask import Flask
import os
application = Flask(__name__)


@application.route('/')
def hello():
    return "Hello World!"

@application.route('/execute')
def execute():
    # py_file = "{}{}autogrow4-4.0.2{}RunAutogrow.py".format(os.getcwd(), os.sep, os.sep)
    # params = "{}{}autogrow4-4.0.2{}sample_sub_scripts{}sample_submit_autogrow.json".format(os.getcwd(), os.sep, os.sep, os.sep)
    py_file = "RunAutogrow.py"
    params = "sample_sub_scripts{}sample_submit_autogrow.json".format(os.sep)
    os.chdir("autogrow4-4.0.2")
    command = "python {} -j {}".format(py_file, params)
    os.system(command)
    return "os.getcwd: {} \ncommand: {}".format(os.getcwd(), command)

if __name__ == '__main__':
    application.debug = True
    application.run()