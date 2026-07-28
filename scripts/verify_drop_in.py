#!/usr/bin/env python3
"""Run the complete local, offline verification for the lightweight Pages build."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run(command: list[str]) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def main() -> int:
    if sys.version_info < (3, 10):
        print("ERROR: Python 3.10 or newer is required.", file=sys.stderr)
        return 2

    node = shutil.which("node")
    if node is None:
        print("ERROR: Node.js 18 or newer is required and 'node' was not found.", file=sys.stderr)
        return 2

    required = [
        ROOT / "data" / "lookup_records_unified.json",
        ROOT / "DROPLET_RISK_unified_summary.json",
        ROOT / "site" / "index.html",
        ROOT / "site" / "assets" / "app.js",
        ROOT / "site" / "assets" / "score.js",
    ]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.is_file()]
    if missing:
        print("ERROR: required repository files are missing:", file=sys.stderr)
        for path in missing:
            print(f"  - {path}", file=sys.stderr)
        return 2

    try:
        run([sys.executable, "scripts/build_web_assets.py"])
        run([node, "--check", "site/assets/app.js"])
        run([node, "--check", "site/assets/score.js"])
        run([node, "tests/test_score.mjs"])
        run([sys.executable, "scripts/validate_web_build.py"])
    except subprocess.CalledProcessError as error:
        print(f"ERROR: verification failed with exit code {error.returncode}.", file=sys.stderr)
        return error.returncode or 1

    print("\nDrop-in verification passed. Serve the site with:")
    print(f"  {sys.executable} -m http.server 8000 --directory site")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
