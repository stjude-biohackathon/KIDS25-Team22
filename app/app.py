from flask import Flask, render_template, request, jsonify
import os

app = Flask(__name__)
app.config["TEMPLATES_AUTO_RELOAD"] = True

@app.after_request
def add_no_cache(resp):
    # prevent aggressive browser/proxy caching during dev
    resp.headers["Cache-Control"] = "no-store"
    return resp

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/chat", methods=["POST"])
def chat():
    data = request.get_json(force=True) or {}
    msg = (data.get("message") or "").strip()
    return jsonify({"reply": f"Dummy response for: {msg}"})

if __name__ == "__main__":
    # ensure you actually run on 5001 if that’s what your browser opens
    app.run(debug=True, use_reloader=True, port=5001)