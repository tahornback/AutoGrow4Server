from flask import Flask
import os
application = Flask(__name__)


@application.route('/')
def hello():
    return "Hello World!"

@application.route('/execute')
def execute():
    os.system("python autogrow4-4.0.2{}RunAutogrow.py -j autogrow4-4.0.2{}sample_sub_scripts{}sample_submit_autogrow.json".format(os.sep, os.sep, os.sep))
    return "Success?"

if __name__ == '__main__':
    application.debug = True
    application.run()