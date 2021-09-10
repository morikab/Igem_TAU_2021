import sys, os
if sys.executable.endswith('pythonw.exe'):
    sys.stdout = open(os.devnull, 'w')
    sys.stderr = open(os.path.join(os.getenv('TEMP'), 'stderr-{}'.format(os.path.basename(sys.argv[0]))), "w")
    
from flask import Flask, redirect, url_for, request, render_template, flash
from flaskwebgui import FlaskUI
from werkzeug.utils import secure_filename


app = Flask(__name__)
# Define the path to the upload folder
UPLOAD_FOLDER = 'static/uploads'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

ui = FlaskUI(app, width=500, height=500, idle_interval=100)


@app.route('/')
def index():
    return render_template("index.html")


@app.route('/success', methods=['POST', 'GET'])
def success():
    if request.method == "POST":
        data = request.form.to_dict(flat=False)
        display_data = request.form.to_dict()
        uploaded_files = []
        print(data)
        if request.files:
            for f in request.files.values():
                filename = secure_filename(f.filename)
                if filename == '':
                    continue
                f.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
                uploaded_files.append(filename)
            display_data['uploaded files'] = uploaded_files
            print('files uploaded successfully', uploaded_files)
        else:
            print('No files uploaded')
        return render_template("success.html", data=display_data)
    return redirect("/")


@app.route('/communique', methods=['POST', 'GET'])
def communique():
    return render_template("optimize_form.html")


if __name__ == '__main__':
    app.run(debug=True)
    # Run in order to make a standalone windows application and comment app.run
    # ui.run()
