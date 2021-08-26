import sys, os
if sys.executable.endswith('pythonw.exe'):
    sys.stdout = open(os.devnull, 'w')
    sys.stderr = open(os.path.join(os.getenv('TEMP'), 'stderr-{}'.format(os.path.basename(sys.argv[0]))), "w")
    
from flask import Flask, redirect, url_for, request, render_template, flash
from flaskwebgui import FlaskUI


app = Flask(__name__)
ui = FlaskUI(app, width=500, height=500, idle_interval=100)


@app.route('/')
def index():
    return render_template("index.html")


@app.route('/success', methods=['POST', 'GET'])
def success():
    return render_template("success.html")


@app.route('/communique', methods=['POST', 'GET'])
def communique():
    return render_template("communique.html")


if __name__ == '__main__':
    app.run(debug=True)
    # ui.run() // Run in order to make a standalone windows application and comment app.run
