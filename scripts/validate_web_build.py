#!/usr/bin/env python3
"""Validate generated Pages assets without third-party dependencies."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Any

ALLOWED_BANDS = {"lower", "leaky", "high", "very_high", "unscored"}
REQUIRED_COLUMNS = {
    "dataset", "name", "id", "riskClass", "score", "index", "mw",
    "logp", "psa", "hbd", "hba", "rtb", "aromaticRings",
}
REQUIRED_HTML_IDS = {
    "loading-banner", "instant-count", "full-row-count", "full-scored-count",
    "shown-count", "model-version", "search-input", "dataset-filter",
    "category-filter", "risk-filter", "sort-filter", "clear-button",
    "export-button", "chembl-fallback", "chembl-fallback-query",
    "chembl-search-button", "local-meta", "local-results-body",
    "chembl-status", "chembl-results", "record-detail",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site", type=Path, default=Path("site"))
    parser.add_argument("--source-lookup", type=Path, default=Path("data/lookup_records_unified.json"))
    return parser.parse_args()


def fail(message: str) -> None:
    raise AssertionError(message)


def load_json(path: Path) -> Any:
    if not path.is_file():
        fail(f"missing file: {path}")
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def decode_value(column: str, value: Any, dictionaries: dict[str, list[Any]]) -> Any:
    if value is None:
        return None
    dictionary = dictionaries.get(column)
    if dictionary is None:
        return value
    if not isinstance(value, int) or isinstance(value, bool):
        fail(f"dictionary column {column!r} contains non-integer code {value!r}")
    if not 0 <= value < len(dictionary):
        fail(f"dictionary code out of range for {column!r}: {value}")
    return dictionary[value]


def validate_lookup(
    payload: Any,
    summary: dict[str, Any],
    source_lookup: Path,
    compact_lookup: Path,
) -> dict[str, Any]:
    if not isinstance(payload, dict):
        fail("compact lookup must be a JSON object")
    if payload.get("version") != 2:
        fail(f"unexpected compact schema version: {payload.get('version')!r}")

    columns = payload.get("columns")
    rows = payload.get("rows")
    dictionaries = payload.get("dictionaries")
    if not isinstance(columns, list) or not all(isinstance(column, str) for column in columns):
        fail("columns must be a list of strings")
    if len(columns) != len(set(columns)):
        fail("columns contain duplicates")
    if not REQUIRED_COLUMNS.issubset(columns):
        fail(f"missing required columns: {sorted(REQUIRED_COLUMNS - set(columns))}")
    if "searchBlob" in columns or "extra" in columns:
        fail("searchBlob or extra leaked into the compact schema")
    if not isinstance(rows, list):
        fail("rows must be a list")
    if not isinstance(dictionaries, dict):
        fail("dictionaries must be an object")

    expected = summary.get("lookup_records")
    if expected is not None and len(rows) != expected:
        fail(f"row count {len(rows):,} does not match summary lookup_records {expected:,}")
    if payload.get("meta", {}).get("rowCount") != len(rows):
        fail("meta.rowCount does not match rows length")

    index_by_column = {column: index for index, column in enumerate(columns)}
    risk_index = index_by_column["riskClass"]
    score_index = index_by_column["score"]
    bounded_index = index_by_column["index"]
    name_index = index_by_column["name"]

    empty_names = 0
    band_counts: dict[str, int] = {band: 0 for band in ALLOWED_BANDS}
    for row_number, row in enumerate(rows, start=1):
        if not isinstance(row, list):
            fail(f"row {row_number} is not a list")
        if len(row) > len(columns):
            fail(f"row {row_number} is longer than the column list")

        padded = row + [None] * (len(columns) - len(row))
        name = decode_value("name", padded[name_index], dictionaries)
        if name is None or str(name).strip() == "":
            empty_names += 1

        band = decode_value("riskClass", padded[risk_index], dictionaries)
        if band not in ALLOWED_BANDS:
            fail(f"row {row_number} has invalid risk band {band!r}")
        band_counts[band] += 1

        score = padded[score_index]
        index = padded[bounded_index]
        if score is not None and (not isinstance(score, (int, float)) or not math.isfinite(float(score))):
            fail(f"row {row_number} has invalid score {score!r}")
        if index is not None:
            if not isinstance(index, (int, float)) or not math.isfinite(float(index)):
                fail(f"row {row_number} has invalid index {index!r}")
            if not 0 <= float(index) <= 1:
                fail(f"row {row_number} has out-of-range bounded index {index!r}")

    if empty_names:
        fail(f"{empty_names} records have empty names")

    size_result: dict[str, Any] = {}
    if source_lookup.is_file():
        compact_bytes = compact_lookup.stat().st_size
        source_bytes = source_lookup.stat().st_size
        if compact_bytes >= source_bytes:
            fail(f"compact lookup is not smaller than source ({compact_bytes:,} >= {source_bytes:,})")
        size_result = {
            "source_bytes": source_bytes,
            "compact_bytes": compact_bytes,
            "reduction_percent": round(100 * (1 - compact_bytes / source_bytes), 2),
        }

    return {"rows": len(rows), "band_counts": band_counts, **size_result}


def validate_site_files(site: Path) -> None:
    index_path = site / "index.html"
    html = index_path.read_text(encoding="utf-8")
    ids = set(re.findall(r'\bid="([^"]+)"', html))
    missing_ids = REQUIRED_HTML_IDS - ids
    if missing_ids:
        fail(f"index.html is missing required IDs: {sorted(missing_ids)}")

    local_refs = re.findall(r'(?:href|src)="(\./[^"]+)"', html)
    for reference in local_refs:
        clean = reference.split("#", 1)[0].split("?", 1)[0]
        path = site / clean.removeprefix("./")
        if not path.exists():
            fail(f"index.html references missing local asset: {reference}")

    app_path = site / "assets" / "app.js"
    score_path = site / "assets" / "score.js"
    if not app_path.is_file() or not score_path.is_file():
        fail("JavaScript assets are missing")
    app_text = app_path.read_text(encoding="utf-8")
    if "https://www.ebi.ac.uk/chembl/api/data" not in app_text:
        fail("app.js does not contain the ChEMBL API base URL")
    if 'from "./score.js"' not in app_text:
        fail("app.js does not import the score module")
    if not (site / ".nojekyll").is_file():
        fail("site/.nojekyll was not generated")


def main() -> int:
    args = parse_args()
    try:
        summary = load_json(args.site / "data" / "atlas-summary.json")
        lookup = load_json(args.site / "data" / "lookup-16k.compact.json")
        result = validate_lookup(
            lookup,
            summary,
            args.source_lookup,
            args.site / "data" / "lookup-16k.compact.json",
        )
        validate_site_files(args.site)
        chembl = summary.get("chembl", {})
        if chembl.get("rows", 0) <= result["rows"]:
            fail("full ChEMBL row count should exceed the local lookup row count")
        print(json.dumps({"status": "ok", **result}, indent=2))
        return 0
    except (AssertionError, OSError, ValueError, TypeError, json.JSONDecodeError) as error:
        print(f"VALIDATION ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
