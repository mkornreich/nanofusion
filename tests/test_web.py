"""The web API, exercised against a real server on a free port."""

from __future__ import annotations

import gzip
import http.client
import json
import threading
import time
from pathlib import Path

import numpy as np
import pytest

from nanofusion.fileio import load_simulation_text
from nanofusion.forms import FormError, properties_from_form
from nanofusion.web.server import make_server

FIXTURES = Path(__file__).parent / "fixtures"
SMALL_RUN = {
    "name": "Web test", "iterations": "8", "timeStep": "2.0E-17", "precision": "6",
    "sizeX": "3", "sizeY": "2", "sizeZ": "2", "fieldX": "1e12", "fieldY": "", "fieldZ": "",
    "threshold": "", "lorentz": True, "logTime": False, "transition": "0", "damping": "",
}


@pytest.fixture()
def server():
    srv = make_server("127.0.0.1", 0)
    thread = threading.Thread(target=srv.serve_forever, daemon=True)
    thread.start()
    yield srv
    srv.shutdown()
    srv.server_close()


def request(srv, method, path, body=None, headers=None):
    conn = http.client.HTTPConnection(*srv.server_address[:2], timeout=30)
    hdrs = {"X-NanoFusion": "1"} if method == "POST" else {}
    hdrs.update(headers or {})
    conn.request(method, path, body=body, headers=hdrs)
    response = conn.getresponse()
    data = response.read()
    conn.close()
    return response, data


def wait_for_job(srv, timeout=60):
    deadline = time.monotonic() + timeout
    while True:
        _, data = request(srv, "GET", "/api/state")
        state = json.loads(data)
        if state["job"] and state["job"]["state"] != "running":
            return state
        assert time.monotonic() < deadline, "job did not finish"
        time.sleep(0.05)


def test_page_and_assets(server):
    response, body = request(server, "GET", "/")
    assert response.status == 200 and b"NanoFusion" in body
    for name, kind in (("app.js", "javascript"), ("viewer.js", "javascript"), ("style.css", "css"),
                       ("about.html", "html")):
        response, _ = request(server, "GET", f"/static/{name}")
        assert response.status == 200 and kind in response.getheader("Content-Type")
    response, _ = request(server, "GET", "/static/../server.py")
    assert response.status == 404


def test_run_frames_and_download(server):
    response, data = request(server, "POST", "/api/run", json.dumps(SMALL_RUN), {"Content-Type": "application/json"})
    assert response.status == 202, data
    state = wait_for_job(server)
    assert state["job"]["state"] == "done"
    assert state["result"]["canSave"] and not state["result"]["saved"]

    _, data = request(server, "GET", "/api/result")
    meta = json.loads(data)
    assert (meta["particles"], meta["frames"], meta["layers"]) == (12, 9, 2)
    assert meta["times"][-1] == "1.60E-16"
    assert meta["fusion"]["count"] == 0
    assert meta["form"]["fieldX"] == "1E+12"

    response, blob = request(server, "GET", "/api/frames")
    frames = np.frombuffer(blob, dtype="<f4").reshape(9, 12, 3)
    assert response.getheader("X-Result-Id") == str(meta["id"])
    assert frames[0, 1, 0] == pytest.approx(2.04)  # ångström

    response, body = request(server, "GET", "/api/download")
    assert response.status == 200
    assert 'filename="Web test.json"' in response.getheader("Content-Disposition")
    saved = json.loads(body)  # standard JSON
    assert len(saved["results"]) == 9
    assert json.loads(request(server, "GET", "/api/state")[1])["result"]["saved"]


def test_validation_errors(server):
    bad = dict(SMALL_RUN, timeStep="abc", sizeX="0", logTime=True, damping="")
    response, data = request(server, "POST", "/api/run", json.dumps(bad), {"Content-Type": "application/json"})
    assert response.status == 400
    errors = json.loads(data)["errors"]
    assert any("Time step" in e for e in errors) and any("Lattice X" in e for e in errors)


def test_open_upload(server):
    with gzip.open(FIXTURES / "fusion_threshold.json.gz", "rb") as f:
        content = f.read()  # original Java-format file
    response, _ = request(server, "POST", "/api/open", content, {"X-Filename": "old%20run.json"})
    assert response.status == 202
    state = wait_for_job(server)
    assert state["job"]["state"] == "done", state["job"]["message"]
    meta = json.loads(request(server, "GET", "/api/result")[1])
    assert meta["origin"] == "file" and meta["source"] == "old run.json"
    assert meta["fusion"]["count"] > 0 and not meta["canSave"]
    assert request(server, "GET", "/api/download")[0].status == 404


def test_cancel(server):
    big = dict(SMALL_RUN, sizeX="10", sizeY="10", sizeZ="10", iterations="10000")
    assert request(server, "POST", "/api/run", json.dumps(big))[0].status == 202
    response, _ = request(server, "POST", "/api/run", json.dumps(SMALL_RUN))
    assert response.status == 409  # one job at a time
    request(server, "POST", "/api/cancel", "")
    assert wait_for_job(server)["job"]["state"] == "cancelled"


def test_rejects_cross_site_and_foreign_hosts(server):
    response, _ = request(server, "POST", "/api/run", json.dumps(SMALL_RUN), {"X-NanoFusion": ""})
    assert response.status == 403
    response, _ = request(server, "GET", "/api/state", headers={"Host": "evil.example:8765"})
    assert response.status == 403
    response, _ = request(server, "GET", "/api/state", headers={"Host": "localhost:8765"})
    assert response.status == 200


def test_form_round_trip():
    props = properties_from_form(SMALL_RUN)
    assert props.electric_field[0] == 1_000_000_000_000 and props.include_lorentz
    with pytest.raises(FormError):
        properties_from_form(dict(SMALL_RUN, sizeX="1000", sizeY="1000"))


def test_load_text_with_progress():
    with gzip.open(FIXTURES / "basic_4x4x4.json.gz", "rt") as f:
        text = f.read()
    seen = []
    result = load_simulation_text(text, "x.json", lambda done, total: seen.append((done, total)))
    assert len(result.frames) == 5 and seen[-1][0] == seen[-1][1]


def test_improved_model_run_and_report(server):
    form = {"model": "improved", "name": "Tiny", "impShape": "sphere", "impSize": "0.6", "impLithium": "Li6",
            "impIntensity": "1e19", "impWavelength": "3.6", "impFwhm": "6", "impPhase": "0", "impCentre": "0",
            "impDuration": "10", "impFrames": "5", "impTarget": "15"}
    response, data = request(server, "POST", "/api/run", json.dumps(form))
    assert response.status == 202, data
    state = wait_for_job(server)
    assert state["job"]["state"] == "done", state["job"]["message"]
    assert state["result"]["canSave"]
    meta = json.loads(request(server, "GET", "/api/result")[1])
    assert meta["model"] == "improved" and meta["form"]["model"] == "improved"
    assert meta["frames"] == 5 and meta["report"]["fusion"]["total"] >= 0
    blob = request(server, "GET", "/api/frames")[1]
    assert len(blob) == 5 * meta["particles"] * 3 * 4
    response, body = request(server, "GET", "/api/download")
    assert response.status == 200 and "report" in response.getheader("Content-Disposition")
    assert json.loads(body)["inputs"]["intensity_w_cm2"] == 1e19
    bad = dict(form, impSize="-1", impIntensity="lots")
    response, data = request(server, "POST", "/api/run", json.dumps(bad))
    assert response.status == 400 and len(json.loads(data)["errors"]) >= 1
