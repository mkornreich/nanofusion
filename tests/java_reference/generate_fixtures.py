"""Regenerate the reference outputs in tests/fixtures/ from the ORIGINAL Java engine.

Needs a JDK (javac/java on PATH) and the extracted Java project in ./java.

    python3 tests/java_reference/generate_fixtures.py
"""

from __future__ import annotations

import gzip
import json
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ENGINE = ROOT / "java" / "FusionSimulation" / "SimulationEngine"
GSON = ENGINE / "release" / "modules" / "ext" / "gson-2.2.4.jar"
FIXTURES = ROOT / "tests" / "fixtures"


def split_args(spec: str) -> list[str]:
    """'a=1 name=Long run' -> ['a=1', 'name=Long run'] (values may contain spaces)."""
    parts: list[str] = []
    for token in shlex.split(spec):
        if "=" in token or not parts:
            parts.append(token)
        else:
            parts[-1] += " " + token
    return parts


def main() -> int:
    if not ENGINE.is_dir():
        print(f"missing {ENGINE}; extract the original nanofusion.adding archive (not in this repository) into ./java first", file=sys.stderr)
        return 1
    cases = json.loads((HERE / "cases.json").read_text())
    FIXTURES.mkdir(parents=True, exist_ok=True)
    expected = {}
    with tempfile.TemporaryDirectory() as tmp:
        classes = Path(tmp) / "classes"
        sources = sorted(str(p) for p in (ENGINE / "src").rglob("*.java"))
        subprocess.run(
            ["javac", "-nowarn", "-cp", str(GSON), "-d", str(classes), *sources, str(HERE / "ReferenceRun.java")],
            check=True,
            stderr=subprocess.DEVNULL,
        )
        for name, spec in cases.items():
            out = Path(tmp) / f"{name}.json"
            proc = subprocess.run(
                ["java", "-cp", f"{classes}:{GSON}", "ReferenceRun", str(out), *split_args(spec)],
                check=True,
                capture_output=True,
                text=True,
            )
            fused = "fusion=true" in proc.stderr
            expected[name] = {"args": spec, "fusion": fused}
            with gzip.open(FIXTURES / f"{name}.json.gz", "wb", compresslevel=9) as gz:
                gz.write(out.read_bytes())
            print(f"{name}: fusion={fused}")
    (FIXTURES / "expected.json").write_text(json.dumps(expected, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
