#!/usr/bin/env python3
"""Build lightweight GitHub Pages assets from the existing 16k lookup JSON.

The source repository stores a rich object per lookup record. This builder:

* removes ``searchBlob`` and ``extra``;
* omits standard interpretation text by default (the browser regenerates it);
* stores records as rows under one shared column list;
* dictionary-encodes repeated strings such as dataset, source, pathway, and note;
* writes minified JSON for the static Pages application.

No third-party Python packages are required.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

COLUMNS = [
    "dataset",
    "name",
    "id",
    "riskClass",
    "score",
    "index",
    "rawClass",
    "category",
    "subclass",
    "source",
    "inchikey",
    "mw",
    "logp",
    "psa",
    "hbd",
    "hba",
    "rtb",
    "aromaticRings",
    "maxPhase",
    "therapeuticFlag",
    "confidence",
    "pathway",
    "note",
    "observed",
    "interpretation",
]

# These values repeat heavily across the ChEMBL subset and benefit from a
# shared dictionary. Name, identifier, and InChIKey remain direct values.
DICTIONARY_COLUMNS = {
    "dataset",
    "riskClass",
    "rawClass",
    "category",
    "subclass",
    "source",
    "confidence",
    "pathway",
    "note",
    "interpretation",
}

SCORE_MODEL = {
    "version": "droplet-risk-unified-v2",
    "coefficients": {
        "logp": 0.70,
        "psa": -0.012,
        "hbd": -0.18,
        "hba": -0.06,
        "mw": -0.0010,
        "rtb": 0.035,
        "aromaticRings": 0.10,
    },
    "cutoffs": {
        "leaky": -0.2563230876636397,
        "high": 0.1453605434241663,
        "veryHigh": 2.1,
    },
    "indexOffset": -2.1,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--lookup",
        type=Path,
        default=Path("data/lookup_records_unified.json"),
        help="Existing rich lookup JSON source.",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("DROPLET_RISK_unified_summary.json"),
        help="Existing atlas summary JSON.",
    )
    parser.add_argument(
        "--site",
        type=Path,
        default=Path("site"),
        help="GitHub Pages site directory.",
    )
    parser.add_argument(
        "--keep-interpretation",
        action="store_true",
        help="Keep record interpretation strings. The default regenerates standard interpretations by risk band in the client.",
    )
    parser.add_argument(
        "--copy-figures",
        type=Path,
        default=None,
        help="Optional figures directory to copy into site/figures.",
    )
    return parser.parse_args()


def json_safe(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def retained_record(record: dict[str, Any], columns: list[str]) -> dict[str, Any]:
    name = record.get("name") or record.get("id") or "Unlabeled record"
    output = {column: json_safe(record.get(column)) for column in columns}
    output["name"] = json_safe(name)
    return output


def build_dictionaries(records: list[dict[str, Any]], columns: list[str]) -> dict[str, list[str]]:
    dictionaries: dict[str, list[str]] = {}
    for column in columns:
        if column not in DICTIONARY_COLUMNS:
            continue
        counts: Counter[str] = Counter()
        for record in records:
            value = record.get(column)
            if value is not None and value != "":
                counts[str(value)] += 1
        # Most frequent first makes common values small integer indices and keeps
        # output deterministic by using the text as the secondary sort key.
        dictionaries[column] = [
            value for value, _ in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
        ]
    return dictionaries


def encode_rows(
    records: list[dict[str, Any]],
    columns: list[str],
    dictionaries: dict[str, list[str]],
) -> list[list[Any]]:
    dictionary_indexes = {
        column: {value: index for index, value in enumerate(values)}
        for column, values in dictionaries.items()
    }
    rows: list[list[Any]] = []
    for record in records:
        row: list[Any] = []
        for column in columns:
            value = record.get(column)
            if column in dictionary_indexes and value is not None and value != "":
                value = dictionary_indexes[column][str(value)]
            row.append(value)
        # Removing trailing nulls is valid because the decoder treats absent row
        # positions as null. It saves space for sparse biological/calibration rows.
        while row and row[-1] is None:
            row.pop()
        rows.append(row)
    return rows


def write_minified_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
        handle.write("\n")
    temporary.replace(path)


def main() -> int:
    args = parse_args()

    if not args.lookup.is_file():
        print(f"ERROR: lookup source not found: {args.lookup}", file=sys.stderr)
        return 2
    if not args.summary.is_file():
        print(f"ERROR: summary source not found: {args.summary}", file=sys.stderr)
        return 2
    if not (args.site / "index.html").is_file():
        print(f"ERROR: site template missing: {args.site / 'index.html'}", file=sys.stderr)
        return 2

    with args.lookup.open("r", encoding="utf-8") as handle:
        source_records = json.load(handle)
    if not isinstance(source_records, list):
        raise TypeError("lookup source must contain a top-level JSON array")
    if not all(isinstance(record, dict) for record in source_records):
        raise TypeError("every lookup record must be a JSON object")

    columns = [column for column in COLUMNS if args.keep_interpretation or column != "interpretation"]
    retained = [retained_record(record, columns) for record in source_records]
    dictionaries = build_dictionaries(retained, columns)
    rows = encode_rows(retained, columns, dictionaries)
    generated_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()

    output_lookup = args.site / "data" / "lookup-16k.compact.json"
    lookup_payload = {
        "version": 2,
        "generatedAt": generated_at,
        "columns": columns,
        "dictionaries": dictionaries,
        "rows": rows,
        "meta": {
            "rowCount": len(rows),
            "source": str(args.lookup).replace("\\", "/"),
            "omittedFields": ["searchBlob", "extra"] + ([] if args.keep_interpretation else ["interpretation"]),
            "encoding": "column rows with per-column dictionaries and trimmed trailing nulls",
            "scoreModel": SCORE_MODEL,
        },
    }
    write_minified_json(output_lookup, lookup_payload)

    with args.summary.open("r", encoding="utf-8") as handle:
        summary = json.load(handle)
    if not isinstance(summary, dict):
        raise TypeError("summary source must contain a top-level JSON object")
    summary["web_build"] = {
        "generated_at": generated_at,
        "lookup_records": len(rows),
        "lookup_asset": "data/lookup-16k.compact.json",
        "lookup_asset_bytes": output_lookup.stat().st_size,
        "architecture": "compact local atlas plus explicit live ChEMBL fallback",
        "live_scores_are_archived_chembl36": False,
    }
    output_summary = args.site / "data" / "atlas-summary.json"
    write_minified_json(output_summary, summary)

    (args.site / ".nojekyll").touch()

    if args.copy_figures:
        if not args.copy_figures.is_dir():
            print(f"ERROR: figures directory not found: {args.copy_figures}", file=sys.stderr)
            return 2
        destination = args.site / "figures"
        if destination.exists():
            shutil.rmtree(destination)
        shutil.copytree(args.copy_figures, destination)

    source_bytes = args.lookup.stat().st_size
    output_bytes = output_lookup.stat().st_size
    reduction = 100.0 * (1.0 - output_bytes / source_bytes) if source_bytes else 0.0
    print(
        json.dumps(
            {
                "records": len(rows),
                "source_bytes": source_bytes,
                "compact_bytes": output_bytes,
                "size_reduction_percent": round(reduction, 2),
                "dictionary_columns": {key: len(values) for key, values in dictionaries.items()},
                "lookup_output": str(output_lookup),
                "summary_output": str(output_summary),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
