import os
import sys
if sys.executable.endswith('pythonw.exe'):
    sys.stdout = open(os.devnull, 'w')
    sys.stderr = open(os.path.join(os.getenv('TEMP'), 'stderr-{}'.format(os.path.basename(sys.argv[0]))), "w")

from pathlib import Path
    
from flask import Flask, redirect, request, render_template
from werkzeug.utils import secure_filename

from GUI import input_for_modules
from modules.main import run_modules

from GUI.flaskgui import FlaskUI

app = Flask(__name__)

# Define the path to the upload folder
UPLOAD_FOLDER = os.path.join("static", "uploads")
Path(UPLOAD_FOLDER).mkdir(parents=True, exist_ok=True)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

IDLE_INTERVAL_SEC = 3 * 60 * 60
ui = FlaskUI(app, width=500, height=500, idle_interval=IDLE_INTERVAL_SEC, host="0.0.0.0", port="5000")


@app.route('/')
def index():
    return render_template("tls_index.html")


@app.route('/success', methods=['POST', 'GET'])
def success():
    if request.method == "POST":
        global data
        data = request.form.to_dict()
        print(data)
        if request.files:
            uploaded_files = []
            uploaded_data = {}
            for key, f in request.files.items():
                filename = secure_filename(f.filename)
                if filename == '':
                    continue
                file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                f.save(file_path)
                uploaded_data[key] = file_path
                uploaded_files.append(filename)
            data['uploaded files'] = uploaded_files
            data["uploaded_data"] = uploaded_data
            print('files uploaded successfully: ', uploaded_files)
            print('Data uploaded organized: ', uploaded_data)

            global processed_user_input, model_preferences
            processed_user_input, model_preferences = input_for_modules.process_input_for_modules(data)
            return render_template("success.html", data={**processed_user_input, **model_preferences})
        else:
            print('No files uploaded')
    return redirect("/")


@app.route('/output', methods=['POST', 'GET'])
def output():
    if request.method == "POST":
        user_output, zip_file_path = run_modules(user_inp_raw=processed_user_input, model_preferences=model_preferences)
        return render_template("user_output.html", user_output=user_output, zip_file_path=zip_file_path)
    return redirect("/")


# @app.route('/communique', methods=['POST', 'GET'])
# def communique():
#     return render_template("optimize_form.html")

