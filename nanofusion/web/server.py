"""Local web server for the NanoFusion user interface (standard library only).

The simulation runs in this Python process; the browser page in ``static/``
edits parameters, starts runs, and draws the frames it fetches from the JSON
API below.

    GET  /                 the page
    GET  /static/<file>    page assets
    GET  /api/state        current job (with progress) and a short result summary
    GET  /api/defaults     default form values
    GET  /api/result       result metadata (properties, particle types, times, fusion events)
    GET  /api/frames       all frame positions: little-endian float32, frames x particles x 3, in ångström
    GET  /api/download     the result as a JSON file (marks it saved)
    POST /api/run          start a simulation; body = form values as JSON
    POST /api/open         load a saved simulation; body = the file's bytes
    POST /api/cancel       cancel the running job

POST requests must carry the header ``X-NanoFusion: 1``. Browsers will not
send that header cross-origin without a CORS preflight, which this server
never approves, so other web sites cannot drive it. When bound to a loopback
address the ``Host`` header must also be local, which blocks DNS rebinding.
"""

from __future__ import annotations

import argparse
import dataclasses
import io
import ipaddress
import json
import mimetypes
import re
import sys
import threading
import time
import webbrowser
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Callable, Optional
from urllib.parse import quote, unquote, urlsplit

import numpy as np

from .. import __version__
from ..engine import SimulationCancelled, detect_fusion, run_simulation
from ..fileio import load_simulation, load_simulation_text, write_simulation
from ..forms import (FormError, all_defaults, form_from_improved, form_from_properties, improved_from_form,
                     properties_from_form)
from ..improved import Cancelled as ImprovedCancelled
from ..improved import ImprovedProperties, run_improved
from ..model import SimulationResult

STATIC_DIR = Path(__file__).parent / "static"
ANGSTROM = 1e10
MAX_UPLOAD = 2 * 1024**3
MAX_EVENTS_LISTED = 200
DEFAULT_PORT = 8765


class Job:
    def __init__(self, kind: str, label: str):
        self.kind = kind
        self.label = label
        self.done = 0
        self.total = 0
        self.state = "running"  # running | done | failed | cancelled
        self.message = ""
        self.cancel_requested = False
        self.started = time.monotonic()
        self.elapsed = 0.0

    def as_json(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "label": self.label,
            "done": self.done,
            "total": self.total,
            "state": self.state,
            "message": self.message,
            "elapsed": round(self.elapsed if self.state != "running" else time.monotonic() - self.started, 2),
        }


class AppState:
    """Everything the server knows. One shared session: this is a local, single-user app."""

    def __init__(self) -> None:
        self.lock = threading.Lock()
        self.result: Optional[SimulationResult] = None
        self.result_id = 0
        self.origin = ""  # "run" or "file"
        self._frames_blob: Optional[bytes] = None
        self.job: Optional[Job] = None

    # --- jobs -------------------------------------------------------------------

    def start_job(self, kind: str, label: str, work: Callable[[Job], SimulationResult],
                  done_message: Callable[[SimulationResult, Job], str]) -> Job:
        with self.lock:
            if self.job is not None and self.job.state == "running":
                raise RuntimeError("Another job is still running.")
            job = Job(kind, label)
            self.job = job

        def runner() -> None:
            try:
                result = work(job)
            except (SimulationCancelled, ImprovedCancelled):
                job.state, job.message = "cancelled", f"{label} cancelled."
            except Exception as exc:  # reported to the page
                job.state, job.message = "failed", f"{type(exc).__name__}: {exc}"
            else:
                with self.lock:
                    self.result = result
                    self.result_id += 1
                    self.origin = kind
                    self._frames_blob = None
                job.message = done_message(result, job)
                job.state = "done"
            finally:
                job.elapsed = time.monotonic() - job.started

        threading.Thread(target=runner, name=f"nanofusion-{kind}", daemon=True).start()
        return job

    def start_run(self, form: dict[str, Any]) -> Job:
        if form.get("model") == "improved":
            return self.start_improved(improved_from_form(form))
        props = properties_from_form(form)

        def work(job: Job) -> SimulationResult:
            job.total = props.iterations

            def progress(done: int, total: int) -> None:
                job.done, job.total = done, total

            return run_simulation(props, progress=progress, cancelled=lambda: job.cancel_requested)

        def message(result: SimulationResult, job: Job) -> str:
            return (f"Simulated {props.particle_count:,} particles for {props.iterations:,} steps "
                    f"in {time.monotonic() - job.started:.1f} s.")

        return self.start_job("run", f"Simulating {props.particle_count:,} particles", work, message)

    def start_improved(self, props: ImprovedProperties) -> Job:
        def work(job: Job) -> SimulationResult:
            def progress(done: int, total: int) -> None:
                job.done, job.total = done, total

            result, report = run_improved(props, progress=progress, cancelled=lambda: job.cancel_requested)
            result.report = dataclasses.asdict(report)
            return result

        def message(result: SimulationResult, job: Job) -> str:
            rep = result.report
            return (f"Improved model: {rep['particles']:,} ions, {rep['simulated_fs']:.0f} fs in "
                    f"{rep['steps']:,} steps, {time.monotonic() - job.started:.1f} s.")

        n_label = "crystal" if props.shape == "sphere" else f"{int(props.size)}³ crystal"
        return self.start_job("run", f"Simulating a laser-stripped LiD {n_label}", work, message)

    def start_open(self, name: str, loader: Callable[[Callable[[int, int], None]], SimulationResult]) -> Job:
        def work(job: Job) -> SimulationResult:
            def progress(done: int, total: int) -> None:
                if job.cancel_requested:
                    raise SimulationCancelled()
                job.done, job.total = done, total

            result = loader(progress)
            result.fusion_events = detect_fusion(result)
            return result

        def message(result: SimulationResult, job: Job) -> str:
            return f"Opened {name} ({len(result.frames):,} frames of {len(result.types):,} particles)."

        return self.start_job("file", f"Opening {name}", work, message)

    # --- views of the result ----------------------------------------------------------

    def summary(self) -> Optional[dict[str, Any]]:
        r = self.result
        if r is None:
            return None
        return {"id": self.result_id, "name": r.properties.name, "canSave": r.can_save or r.report is not None,
                "saved": r.saved}

    def metadata(self) -> Optional[dict[str, Any]]:
        with self.lock:
            r, rid, origin = self.result, self.result_id, self.origin
        if r is None:
            return None
        lo, hi = r.bounds()
        events = r.fusion_events
        fusion = None
        if events is not None:
            fusion = {
                "count": len(events),
                "events": [
                    {"frame": e.frame, "p1": e.particle_1, "p2": e.particle_2, "distance": f"{float(e.distance):.4e}"}
                    for e in events[:MAX_EVENTS_LISTED]
                ],
            }
        if r.report is not None:
            form = form_from_improved(ImprovedProperties(**r.report["inputs"]))
        else:
            form = {"model": "original", **form_from_properties(r.properties)}
        return {
            "id": rid,
            "origin": origin,
            "model": "improved" if r.report is not None else "original",
            "report": r.report,
            "form": form,
            "name": r.properties.name,
            "particles": len(r.types),
            "frames": len(r.frames),
            "isLi": [1 if t == "Li" else 0 for t in r.types],
            "layer": r.lattice_index[:, 2].tolist(),
            "layers": int(r.lattice_index[:, 2].max()) + 1 if len(r.types) else 0,
            "times": [str(f.time) for f in r.frames],
            "executionMs": [f.execution_time_ms for f in r.frames],
            "bounds": {"lo": (lo * ANGSTROM).tolist(), "hi": (hi * ANGSTROM).tolist()},
            "fusion": fusion,
            "threshold": f"{float(r.properties.threshold):.3g}",
            "canSave": r.can_save or r.report is not None,
            "saved": r.saved,
            "source": Path(r.source_path).name if r.source_path else None,
            "warning": r.warning,
        }

    def frames_blob(self) -> Optional[bytes]:
        with self.lock:
            if self.result is None:
                return None
            if self._frames_blob is None:
                stacked = np.stack([f.positions for f in self.result.frames]) * ANGSTROM
                self._frames_blob = stacked.astype("<f4").tobytes()
            return self._frames_blob


def _is_loopback(host: str) -> bool:
    if host == "localhost":
        return True
    try:
        return ipaddress.ip_address(host).is_loopback
    except ValueError:
        return False


def _safe_filename(name: str) -> str:
    cleaned = re.sub(r'[\\/:*?"<>|\x00-\x1f]+', "_", name).strip(" .") or "simulation"
    return cleaned[:120]


class Handler(BaseHTTPRequestHandler):
    server_version = f"NanoFusion/{__version__}"
    app: AppState  # set on the subclass created in make_server()
    loopback_only: bool = True

    # --- plumbing -----------------------------------------------------------------

    def log_message(self, format: str, *args: Any) -> None:  # noqa: A002 (stdlib signature)
        pass  # keep the console quiet

    def _send(self, status: int, body: bytes, content_type: str, headers: Optional[dict[str, str]] = None) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        for key, value in (headers or {}).items():
            self.send_header(key, value)
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def _json(self, data: Any, status: int = HTTPStatus.OK) -> None:
        self._send(status, json.dumps(data).encode(), "application/json")

    def _error(self, status: int, message: str, errors: Optional[list[str]] = None) -> None:
        self._json({"error": message, "errors": errors or [message]}, status)

    def _host_ok(self) -> bool:
        if not self.loopback_only:
            return True
        host = (self.headers.get("Host") or "").strip()
        if host.startswith("["):
            host = host[1 : host.find("]")] if "]" in host else host
        else:
            host = host.split(":", 1)[0]
        return _is_loopback(host)

    def _guard(self, mutating: bool) -> bool:
        if not self._host_ok():
            self._error(HTTPStatus.FORBIDDEN, "Unexpected Host header.")
            return False
        if mutating and self.headers.get("X-NanoFusion") != "1":
            self._error(HTTPStatus.FORBIDDEN, "Missing X-NanoFusion header.")
            return False
        return True

    def _read_body(self) -> Optional[bytes]:
        try:
            length = int(self.headers.get("Content-Length", ""))
        except ValueError:
            self._error(HTTPStatus.LENGTH_REQUIRED, "Content-Length required.")
            return None
        if length < 0 or length > MAX_UPLOAD:
            self._error(HTTPStatus.REQUEST_ENTITY_TOO_LARGE, "File too large.")
            return None
        return self.rfile.read(length)

    # --- routes -------------------------------------------------------------------

    def do_HEAD(self) -> None:  # noqa: N802
        self.do_GET()

    def do_GET(self) -> None:  # noqa: N802
        if not self._guard(mutating=False):
            return
        path = urlsplit(self.path).path
        app = self.app
        if path in ("/", "/index.html"):
            self._static("index.html")
        elif path.startswith("/static/"):
            self._static(path[len("/static/"):])
        elif path == "/api/state":
            job = app.job
            self._json({"job": job.as_json() if job else None, "result": app.summary(), "version": __version__})
        elif path == "/api/defaults":
            self._json(all_defaults())
        elif path == "/api/result":
            meta = app.metadata()
            if meta is None:
                self._error(HTTPStatus.NOT_FOUND, "No result yet.")
            else:
                self._json(meta)
        elif path == "/api/frames":
            blob = app.frames_blob()
            if blob is None:
                self._error(HTTPStatus.NOT_FOUND, "No result yet.")
            else:
                self._send(HTTPStatus.OK, blob, "application/octet-stream",
                           {"X-Result-Id": str(app.result_id)})
        elif path == "/api/download":
            self._download()
        elif path == "/favicon.ico":
            self._send(HTTPStatus.NO_CONTENT, b"", "image/x-icon")
        else:
            self._error(HTTPStatus.NOT_FOUND, "Not found.")

    def do_POST(self) -> None:  # noqa: N802
        if not self._guard(mutating=True):
            return
        path = urlsplit(self.path).path
        app = self.app
        if path == "/api/run":
            body = self._read_body()
            if body is None:
                return
            try:
                form = json.loads(body or b"{}")
                if not isinstance(form, dict):
                    raise ValueError
            except ValueError:
                self._error(HTTPStatus.BAD_REQUEST, "Body must be a JSON object.")
                return
            try:
                job = app.start_run(form)
            except FormError as exc:
                self._error(HTTPStatus.BAD_REQUEST, "Check the parameters.", exc.errors)
                return
            except RuntimeError as exc:
                self._error(HTTPStatus.CONFLICT, str(exc))
                return
            self._json({"job": job.as_json()}, HTTPStatus.ACCEPTED)
        elif path == "/api/open":
            name = unquote(self.headers.get("X-Filename", "simulation.json"))[:200]
            body = self._read_body()
            if body is None:
                return
            try:
                text = body.decode("utf-8")
            except UnicodeDecodeError:
                self._error(HTTPStatus.BAD_REQUEST, f"{name} is not a UTF-8 text file.")
                return
            del body
            try:
                job = app.start_open(name, lambda progress: load_simulation_text(text, name, progress))
            except RuntimeError as exc:
                self._error(HTTPStatus.CONFLICT, str(exc))
                return
            self._json({"job": job.as_json()}, HTTPStatus.ACCEPTED)
        elif path == "/api/cancel":
            job = app.job
            if job is not None and job.state == "running":
                job.cancel_requested = True
            self._json({"job": job.as_json() if job else None})
        else:
            self._error(HTTPStatus.NOT_FOUND, "Not found.")

    def _static(self, name: str) -> None:
        allowed = {p.name: p for p in STATIC_DIR.iterdir() if p.is_file()}
        target = allowed.get(name)
        if target is None:
            self._error(HTTPStatus.NOT_FOUND, "Not found.")
            return
        content_type = {
            ".js": "text/javascript; charset=utf-8",
            ".css": "text/css; charset=utf-8",
            ".html": "text/html; charset=utf-8",
            ".svg": "image/svg+xml",
        }.get(target.suffix) or mimetypes.guess_type(target.name)[0] or "application/octet-stream"
        self._send(HTTPStatus.OK, target.read_bytes(), content_type)

    def _download(self) -> None:
        app = self.app
        with app.lock:
            result = app.result
        if result is None or not (result.can_save or result.report is not None):
            self._error(HTTPStatus.NOT_FOUND, "Nothing to save: run a simulation first.")
            return
        if result.report is not None:  # improved model: the report, not the frames
            body = json.dumps(result.report, indent=1).encode("utf-8")
            filename = _safe_filename(result.properties.name) + " report.json"
        else:
            buffer = io.StringIO()
            write_simulation(result, buffer)
            body = buffer.getvalue().encode("utf-8")
            del buffer
            filename = _safe_filename(result.properties.name) + ".json"
        disposition = f"attachment; filename=\"{filename.encode('ascii', 'replace').decode()}\"; " \
                      f"filename*=UTF-8''{quote(filename)}"
        self._send(HTTPStatus.OK, body, "application/json", {"Content-Disposition": disposition})
        result.saved = True


def make_server(host: str = "127.0.0.1", port: int = DEFAULT_PORT,
                state: Optional[AppState] = None) -> ThreadingHTTPServer:
    app = state or AppState()
    handler = type("NanoFusionHandler", (Handler,), {"app": app, "loopback_only": _is_loopback(host)})
    server = ThreadingHTTPServer((host, port), handler)
    server.daemon_threads = True
    server.app = app  # type: ignore[attr-defined]
    return server


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="nanofusion",
        description="Start the NanoFusion web interface (use 'nanofusion run --help' for the command line).",
    )
    parser.add_argument("file", nargs="?", help="saved simulation to open at start-up")
    parser.add_argument("--host", default="127.0.0.1", help="address to listen on (default: 127.0.0.1)")
    parser.add_argument("--port", type=int, default=None, help=f"port (default: {DEFAULT_PORT}, or a free one)")
    parser.add_argument("--no-browser", action="store_true", help="do not open a browser window")
    args = parser.parse_args(argv)

    if args.port is not None:
        server = make_server(args.host, args.port)
    else:
        try:
            server = make_server(args.host, DEFAULT_PORT)
        except OSError:
            server = make_server(args.host, 0)  # default port busy: let the OS pick one
    app: AppState = server.app  # type: ignore[attr-defined]

    if args.file:
        path = Path(args.file)
        if not path.is_file():
            print(f"error: {path} does not exist", file=sys.stderr)
            return 2
        app.start_open(path.name, lambda progress: load_simulation(path, progress))

    host, port = server.server_address[:2]
    shown_host = "localhost" if _is_loopback(str(host)) else host
    url = f"http://{shown_host}:{port}/"
    print(f"NanoFusion {__version__} is running at {url}  (Ctrl+C to stop)", flush=True)
    if not args.no_browser:
        threading.Timer(0.3, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopped.")
    finally:
        server.server_close()
    return 0
