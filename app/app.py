from flask import Flask, render_template, request, jsonify, Response, stream_with_context
import contextlib
import io
import os
import queue
import sys
import threading
import uuid
from pathlib import Path
from werkzeug.utils import secure_filename

ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.append(str(ROOT_DIR))

import main as orchestrator_main

UPLOADS_DIR = ROOT_DIR / "uploads"
UPLOADS_DIR.mkdir(parents=True, exist_ok=True)


def _is_allowed_vcf(filename: str) -> bool:
    lowered = filename.lower()
    return lowered.endswith(".vcf") or lowered.endswith(".vcf.gz")


def _unique_destination(directory: Path, filename: str) -> Path:
    candidate = directory / filename
    if not candidate.exists():
        return candidate

    original = Path(filename)
    suffix = "".join(original.suffixes)
    base = original.name[: -len(suffix)] if suffix else original.name
    counter = 1
    while True:
        candidate = directory / f"{base}_{counter}{suffix}"
        if not candidate.exists():
            return candidate
        counter += 1


LOG_STREAMS = {}
LOG_LOCK = threading.Lock()


def _ensure_log_queue(log_id: str) -> queue.Queue:
    with LOG_LOCK:
        existing = LOG_STREAMS.get(log_id)
        if existing is not None:
            return existing
        q = queue.Queue()
        LOG_STREAMS[log_id] = q
        return q


def _publish_log(log_id: str, message: str) -> None:
    if message is None:
        return
    _ensure_log_queue(log_id).put(str(message))


def _finalize_log(log_id: str) -> None:
    with LOG_LOCK:
        q = LOG_STREAMS.pop(log_id, None)
    if q is not None:
        q.put(None)


class QueueStream(io.TextIOBase):
    def __init__(self, original: io.TextIOBase, log_id: str):
        self._original = original
        self._log_id = log_id
        self._buffer = ""

    def write(self, data: str) -> int:
        if not data:
            return 0
        self._original.write(data)
        self._original.flush()
        combined = self._buffer + data
        parts = combined.split("\n")
        self._buffer = parts.pop() if parts else ""
        for line in parts:
            _publish_log(self._log_id, line.rstrip("\r"))
        return len(data)

    def flush(self) -> None:
        self._original.flush()
        if self._buffer:
            _publish_log(self._log_id, self._buffer.rstrip("\r"))
            self._buffer = ""


@contextlib.contextmanager
def stream_logs_to_queue(log_id: str):
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    stdout_proxy = QueueStream(original_stdout, log_id)
    stderr_proxy = QueueStream(original_stderr, log_id)
    sys.stdout = stdout_proxy
    sys.stderr = stderr_proxy
    try:
        yield
    finally:
        stdout_proxy.flush()
        stderr_proxy.flush()
        sys.stdout = original_stdout
        sys.stderr = original_stderr

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

@app.route("/logs/<log_id>")
def stream_logs(log_id: str):
    log_queue = _ensure_log_queue(log_id)

    def generate():
        while True:
            message = log_queue.get()
            if message is None:
                yield "event: done\ndata:\n\n"
                break
            yield f"data: {message}\n\n"

    return Response(stream_with_context(generate()), mimetype="text/event-stream")

@app.route("/chat", methods=["POST"])
def chat():
    if request.content_type and request.content_type.startswith("multipart/form-data"):
        log_id = (request.form.get("log_id") or str(uuid.uuid4())).strip()
        _ensure_log_queue(log_id)

        phenotype = (request.form.get("phenotype") or "").strip()
        _ = (request.form.get("message") or "").strip()  # optional free-text notes (currently unused)
        uploaded_file = request.files.get("vcf_file")

        if not phenotype and not uploaded_file:
            error_msg = "Provide a phenotype description or attach a VCF file."
            _publish_log(log_id, f"Error: {error_msg}")
            _finalize_log(log_id)
            return jsonify({"error": error_msg, "log_id": log_id}), 400

        saved_path = None
        if uploaded_file and uploaded_file.filename:
            if not _is_allowed_vcf(uploaded_file.filename):
                error_msg = "Only .vcf or .vcf.gz files are supported."
                _publish_log(log_id, f"Error: {error_msg}")
                _finalize_log(log_id)
                return jsonify({"error": error_msg, "log_id": log_id}), 400

            filename = secure_filename(uploaded_file.filename)
            if not filename:
                error_msg = "Invalid filename."
                _publish_log(log_id, f"Error: {error_msg}")
                _finalize_log(log_id)
                return jsonify({"error": error_msg, "log_id": log_id}), 400

            destination = _unique_destination(UPLOADS_DIR, filename)
            uploaded_file.save(destination)
            saved_path = str(destination)

        status = 200
        payload = {"log_id": log_id}
        try:
            _publish_log(log_id, "Starting variant prioritization run…")
            with stream_logs_to_queue(log_id):
                result = orchestrator_main.main(
                    file_path=saved_path,
                    phenotype=phenotype or None,
                    use_cli_args=False,
                )
        except Exception as exc:
            status = 500
            error_text = str(exc)
            _publish_log(log_id, f"Error: {error_text}")
            payload["error"] = error_text
        else:
            _publish_log(log_id, "Run complete.")
            payload["reply"] = str(result)
        finally:
            _finalize_log(log_id)

        return jsonify(payload), status

    data = request.get_json(force=True) or {}
    msg = (data.get("message") or "").strip()
    return jsonify({"reply": f"Dummy response for: {msg}"})

if __name__ == "__main__":
    # ensure you actually run on 5001 if that’s what your browser opens
    app.run(debug=True, use_reloader=True, port=5001)
