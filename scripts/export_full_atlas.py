#!/usr/bin/env python3
"""Export the complete ChEMBL property table as compact scored JSONL.gz parts.

This script is optional for the website. The Pages interface scores fallback
records from the live ChEMBL API. Use this exporter to publish the archived
ChEMBL36 row-level atlas as GitHub Release assets without loading it into the
website.

The exporter uses only the Python standard library and processes the input one
row at a time. It accepts CSV, CSV.GZ, or a ZIP containing a CSV.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import math
import sys
import zipfile
from collections import Counter
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator, TextIO

LEAKY_CUTOFF = -0.2563230876636397
HIGH_CUTOFF = 0.1453605434241663
VERY_HIGH_CUTOFF = 2.1
INDEX_OFFSET = -2.1

ALIASES = {
    "chembl_id": ["chembl_id", "molecule_chembl_id"],
    "pref_name": ["pref_name"],
    "inchikey": [
        "compound_structures__standard_inchi_key",
        "molecule_structures__standard_inchi_key",
        "standard_inchi_key",
    ],
    "molecule_type": ["molecule_type"],
    "max_phase": ["max_phase"],
    "therapeutic_flag": ["therapeutic_flag"],
    "logp": ["compound_properties__alogp", "molecule_properties__alogp", "alogp", "logp"],
    "psa": ["compound_properties__psa", "molecule_properties__psa", "psa", "tpsa"],
    "hbd": ["compound_properties__hbd", "molecule_properties__hbd", "hbd"],
    "hba": ["compound_properties__hba", "molecule_properties__hba", "hba"],
    "mw": ["compound_properties__full_mwt", "molecule_properties__full_mwt", "full_mwt", "mw"],
    "rtb": ["compound_properties__rtb", "molecule_properties__rtb", "rtb"],
    "aromatic_rings": [
        "compound_properties__aromatic_rings",
        "molecule_properties__aromatic_rings",
        "aromatic_rings",
    ],
}

SCHEMA = {
    "format": "newline-delimited JSON compressed with gzip",
    "record_keys": {
        "i": "ChEMBL molecule identifier",
        "n": "preferred name; omitted when missing",
        "k": "standard InChIKey; omitted when missing",
        "t": "molecule type; omitted when missing",
        "p": "maximum clinical phase; omitted when missing",
        "tf": "therapeutic flag; omitted when missing",
        "s": "DROPLET-RISK score; omitted when unscored",
        "x": "bounded DROPLET-RISK index; omitted when unscored",
        "b": "risk band: lower, leaky, high, very_high, or unscored",
        "d": "descriptor array in descriptor_order; null values indicate missing descriptors",
    },
    "descriptor_order": ["alogp", "psa", "hbd", "hba", "full_mwt", "rtb", "aromatic_rings"],
    "score_equation": "0.70*alogp - 0.012*psa - 0.18*hbd - 0.06*hba - 0.0010*full_mwt + 0.035*rtb + 0.10*aromatic_rings",
    "index_equation": "sigmoid(-2.1 + score)",
    "cutoffs": {
        "lower_below": LEAKY_CUTOFF,
        "leaky_below": HIGH_CUTOFF,
        "high_below": VERY_HIGH_CUTOFF,
        "very_high_at_or_above": VERY_HIGH_CUTOFF,
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="ChEMBL property CSV, CSV.GZ, or ZIP containing a CSV.")
    parser.add_argument("--output-dir", type=Path, default=Path("release-assets"))
    parser.add_argument("--rows-per-part", type=int, default=500_000)
    parser.add_argument("--expected-rows", type=int, default=2_854_815)
    parser.add_argument("--allow-row-count-mismatch", action="store_true")
    return parser.parse_args()


def clean_text(value: Any) -> str | None:
    text = str(value or "").strip()
    if not text or text.lower() in {"nan", "none", "null", "<na>"}:
        return None
    return text


def number(value: Any) -> float | None:
    text = clean_text(value)
    if text is None:
        return None
    try:
        result = float(text)
    except ValueError:
        return None
    return result if math.isfinite(result) else None


def integer_or_float(value: Any) -> int | float | None:
    result = number(value)
    if result is None:
        return None
    return int(result) if result.is_integer() else result


def boolean_or_value(value: Any) -> bool | int | str | None:
    text = clean_text(value)
    if text is None:
        return None
    lowered = text.lower()
    if lowered in {"true", "t", "yes", "y", "1"}:
        return True
    if lowered in {"false", "f", "no", "n", "0"}:
        return False
    return text


def first_present(row: dict[str, str], names: list[str]) -> Any:
    for name in names:
        if name in row:
            value = row.get(name)
            if clean_text(value) is not None:
                return value
    return None


def sigmoid(value: float) -> float:
    if value >= 0:
        z = math.exp(-value)
        return 1.0 / (1.0 + z)
    z = math.exp(value)
    return z / (1.0 + z)


def score(descriptors: list[float | None]) -> tuple[float | None, float | None, str]:
    if any(value is None for value in descriptors):
        return None, None, "unscored"
    logp, psa, hbd, hba, mw, rtb, aromatic_rings = (float(value) for value in descriptors)
    raw = (
        0.70 * logp
        - 0.012 * psa
        - 0.18 * hbd
        - 0.06 * hba
        - 0.0010 * mw
        + 0.035 * rtb
        + 0.10 * aromatic_rings
    )
    comparison_score = round(raw, 10)
    if comparison_score < LEAKY_CUTOFF:
        band = "lower"
    elif comparison_score < HIGH_CUTOFF:
        band = "leaky"
    elif comparison_score < VERY_HIGH_CUTOFF:
        band = "high"
    else:
        band = "very_high"
    return raw, sigmoid(INDEX_OFFSET + raw), band


@contextmanager
def open_input(path: Path) -> Iterator[TextIO]:
    suffixes = [suffix.lower() for suffix in path.suffixes]
    if suffixes and suffixes[-1] == ".zip":
        archive = zipfile.ZipFile(path)
        candidates = [name for name in archive.namelist() if name.lower().endswith(".csv") and not name.endswith("/")]
        if not candidates:
            archive.close()
            raise FileNotFoundError("The ZIP archive contains no CSV file.")
        selected = max(candidates, key=lambda name: archive.getinfo(name).file_size)
        binary = archive.open(selected, "r")
        text = io.TextIOWrapper(binary, encoding="utf-8-sig", newline="")
        try:
            yield text
        finally:
            text.close()
            archive.close()
    elif suffixes[-2:] == [".csv", ".gz"] or (suffixes and suffixes[-1] == ".gz"):
        with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
            yield handle
    else:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            yield handle


class PartWriter:
    def __init__(self, output_dir: Path, rows_per_part: int):
        self.output_dir = output_dir
        self.rows_per_part = rows_per_part
        self.part_number = 0
        self.rows_in_part = 0
        self.total_rows = 0
        self.handle: TextIO | None = None
        self.parts: list[dict[str, Any]] = []

    def _open_next(self) -> None:
        self.close_current()
        self.part_number += 1
        filename = f"droplet-risk-chembl36.compact.part-{self.part_number:03d}.jsonl.gz"
        path = self.output_dir / filename
        self.handle = gzip.open(path, "wt", encoding="utf-8", newline="\n", compresslevel=9)
        self.rows_in_part = 0
        self.parts.append({"filename": filename, "rows": 0})

    def write(self, record: dict[str, Any]) -> None:
        if self.handle is None or self.rows_in_part >= self.rows_per_part:
            self._open_next()
        assert self.handle is not None
        self.handle.write(json.dumps(record, ensure_ascii=False, separators=(",", ":"), allow_nan=False))
        self.handle.write("\n")
        self.rows_in_part += 1
        self.total_rows += 1
        self.parts[-1]["rows"] = self.rows_in_part

    def close_current(self) -> None:
        if self.handle is not None:
            self.handle.close()
            self.handle = None
            path = self.output_dir / self.parts[-1]["filename"]
            self.parts[-1]["bytes"] = path.stat().st_size

    def close(self) -> None:
        self.close_current()


def build_record(row: dict[str, str]) -> tuple[dict[str, Any], list[str]]:
    chembl_id = clean_text(first_present(row, ALIASES["chembl_id"]))
    if not chembl_id:
        raise ValueError("A row is missing its ChEMBL identifier.")

    descriptors = [
        number(first_present(row, ALIASES["logp"])),
        number(first_present(row, ALIASES["psa"])),
        number(first_present(row, ALIASES["hbd"])),
        number(first_present(row, ALIASES["hba"])),
        number(first_present(row, ALIASES["mw"])),
        number(first_present(row, ALIASES["rtb"])),
        number(first_present(row, ALIASES["aromatic_rings"])),
    ]
    missing = [SCHEMA["descriptor_order"][index] for index, value in enumerate(descriptors) if value is None]
    raw_score, risk_index, band = score(descriptors)

    record: dict[str, Any] = {"i": chembl_id, "b": band, "d": descriptors}
    optional = {
        "n": clean_text(first_present(row, ALIASES["pref_name"])),
        "k": clean_text(first_present(row, ALIASES["inchikey"])),
        "t": clean_text(first_present(row, ALIASES["molecule_type"])),
        "p": integer_or_float(first_present(row, ALIASES["max_phase"])),
        "tf": boolean_or_value(first_present(row, ALIASES["therapeutic_flag"])),
    }
    record.update({key: value for key, value in optional.items() if value is not None})
    if raw_score is not None and risk_index is not None:
        record["s"] = round(raw_score, 7)
        record["x"] = round(risk_index, 8)
    return record, missing


def main() -> int:
    args = parse_args()
    if not args.input.is_file():
        print(f"ERROR: input file not found: {args.input}", file=sys.stderr)
        return 2
    if args.rows_per_part < 1:
        print("ERROR: --rows-per-part must be positive", file=sys.stderr)
        return 2

    args.output_dir.mkdir(parents=True, exist_ok=True)
    writer = PartWriter(args.output_dir, args.rows_per_part)
    risk_counts: Counter[str] = Counter()
    missing_counts: Counter[str] = Counter()
    named = 0
    max_phase_nonmissing = 0
    molecule_types: Counter[str] = Counter()
    started = datetime.now(timezone.utc)

    try:
        with open_input(args.input) as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                raise ValueError("Input CSV has no header row.")

            available = set(reader.fieldnames)
            unresolved = {
                key: names
                for key, names in ALIASES.items()
                if key in {"chembl_id", "logp", "psa", "hbd", "hba", "mw", "rtb", "aromatic_rings"}
                and not any(name in available for name in names)
            }
            if unresolved:
                formatted = "; ".join(f"{key}: {', '.join(names)}" for key, names in unresolved.items())
                raise KeyError(f"Required columns could not be resolved: {formatted}")

            for row_number, row in enumerate(reader, start=2):
                try:
                    record, missing = build_record(row)
                except Exception as exc:
                    raise RuntimeError(f"Failed at CSV row {row_number}: {exc}") from exc
                writer.write(record)
                risk_counts[record["b"]] += 1
                missing_counts.update(missing)
                if record.get("n"):
                    named += 1
                if record.get("p") is not None:
                    max_phase_nonmissing += 1
                molecule_types[record.get("t") or "missing"] += 1
                if writer.total_rows % 100_000 == 0:
                    print(f"Processed {writer.total_rows:,} rows…", flush=True)
    finally:
        writer.close()

    finished = datetime.now(timezone.utc)
    manifest = {
        "created_at": finished.replace(microsecond=0).isoformat(),
        "source_file": args.input.name,
        "total_rows": writer.total_rows,
        "named": named,
        "max_phase_nonmissing": max_phase_nonmissing,
        "risk_band_counts": dict(risk_counts),
        "molecule_type_counts": dict(molecule_types),
        "missing_descriptor_counts": dict(missing_counts),
        "parts": writer.parts,
        "elapsed_seconds": round((finished - started).total_seconds(), 2),
        "schema": SCHEMA,
        "usage_note": "Downloadable archived atlas. The GitHub Pages application does not load these files during normal use.",
    }
    manifest_path = args.output_dir / "droplet-risk-chembl36.compact.manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if writer.total_rows != args.expected_rows and not args.allow_row_count_mismatch:
        print(
            f"ERROR: exported {writer.total_rows:,} rows, expected {args.expected_rows:,}. "
            "Inspect the input or rerun with --allow-row-count-mismatch.",
            file=sys.stderr,
        )
        return 3

    print(json.dumps({"manifest": str(manifest_path), **manifest}, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
