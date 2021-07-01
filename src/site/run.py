import os
import sqlite3

from flask import abort, Flask, render_template, request, Response
from gizmos.tree import tree
from gizmos.search import search
from markdown import markdown

# Set up the app
app = Flask(__name__)


@app.route("/")
def index():
    with open("src/site/index.md", "r") as f:
        html = markdown(f.read())
    return render_template("base.html", content=html)


@app.route("/export/<export_file>.tsv")
def download_export(export_file):
    filepath = os.path.join("build", export_file + ".tsv")
    if not os.path.exists(filepath):
        abort(400, f"Unknown file: {export_file}.tsv")
    with open(filepath, "r") as f:
        return Response(f.read(), mimetype="text/tab-separated-values")


@app.route("/ontology")
def browse_mro():
    return render_tree(None)


@app.route("/ontology/<term_id>")
def browse_mro_at(term_id):
    return render_tree(term_id)


@app.route("/search")
def search_mro():
    text = request.args.get("text", "")
    db = "build/mro.db"
    if not os.path.exists(db):
        abort(500, "A database for MRO (build/mro.db) does not exist!")
    conn = sqlite3.connect(db)
    return search(conn, text)


def render_tree(term_id):
    db = "build/mro.db"
    if not os.path.exists(db):
        abort(500, "A database for MRO (build/mro.db) does not exist!")
    conn = sqlite3.connect(db)
    html = tree(conn, "MRO", term_id, href="/ontology/{curie}", include_search=True, standalone=False)
    return render_template("base.html", content=html)
