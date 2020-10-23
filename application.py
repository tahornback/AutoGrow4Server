import os
if __name__ == "__main__":
    os.system("/conda/bin/python ./flask_app.py")
else:
    raise Exception("application.py not started as __main__")