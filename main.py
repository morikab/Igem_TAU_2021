import os
import sys
if sys.executable.endswith('pythonw.exe'):
    sys.stdout = open(os.devnull, 'w')
    sys.stderr = open(os.path.join(os.getenv('TEMP'), 'stderr-{}'.format(os.path.basename(sys.argv[0]))), "w")
    
import os
import sys
from pathlib import Path
    
from flask import Flask, redirect, request, render_template
from flaskwebgui import FlaskUI
from werkzeug.utils import secure_filename

from GUI import input_for_modules
from modules.main import run_modules


if sys.executable.endswith('pythonw.exe'):
    sys.stdout = open(os.devnull, 'w')
    sys.stderr = open(os.path.join(os.getenv('TEMP'), 'stderr-{}'.format(os.path.basename(sys.argv[0]))), "w")

app = Flask(__name__)
# Define the path to the upload folder
UPLOAD_FOLDER = os.path.join("static", "uploads")
Path(UPLOAD_FOLDER).mkdir(parents=True, exist_ok=True)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

ui = FlaskUI(app, width=500, height=500, idle_interval=100)


@app.route('/')
def index():
    return render_template("index.html")


@app.route('/success', methods=['POST', 'GET'])
def success():
    if request.method == "POST":
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
            data['uploaded files'] = uploaded_files     # TODO - can we get rid of uploaded_files?
            data["uploaded_data"] = uploaded_data
            print('files uploaded successfully: ', uploaded_files)
            print('Data uploaded organized: ', uploaded_data)
        else:
            print('No files uploaded')
        # TODO - need to create another screen with summarized info, and only then run the analysis
        processed_user_input = input_for_modules.process_input_for_modules(data)
        user_output = run_modules(processed_user_input)
        return render_template("success.html", data=data)
    return redirect("/")


@app.route('/communique', methods=['POST', 'GET'])
def communique():
    return render_template("optimize_form.html")


if __name__ == '__main__':
    # app.run(debug=True)
    # Run in order to make a standalone windows application and comment app.run
    ui.run()
