#!/usr/bin/env python3
"""Optional live check that ChEMBL exposes the seven scoring descriptors."""

from __future__ import annotations

import argparse
import json
import math
import sys
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

FIELDS = [
    "molecule_chembl_id",
    "pref_name",
    "molecule_properties",
    "molecule_hierarchy",
]
DESCRIPTORS = {
    "alogp": "alogp",
    "psa": "psa",
    "hbd": "hbd",
    "hba": "hba",
    "full_mwt": "full_mwt",
    "rtb": "rtb",
    "aromatic_rings": "aromatic_rings",
}
LEAKY_CUTOFF = -0.2563230876636397
HIGH_CUTOFF = 0.1453605434241663
VERY_HIGH_CUTOFF = 2.1


def number(value):
    if value is None or value == "":
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def sigmoid(value: float) -> float:
    if value >= 0:
        z = math.exp(-value)
        return 1 / (1 + z)
    z = math.exp(value)
    return z / (1 + z)


def band(score: float) -> str:
    comparison = round(score, 10)
    if comparison < LEAKY_CUTOFF:
        return "lower"
    if comparison < HIGH_CUTOFF:
        return "leaky"
    if comparison < VERY_HIGH_CUTOFF:
        return "high"
    return "very_high"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("chembl_id", nargs="?", default="CHEMBL25")
    parser.add_argument("--timeout", type=float, default=20.0)
    args = parser.parse_args()

    chembl_id = args.chembl_id.strip().upper()
    query = urlencode({"only": ",".join(FIELDS)})
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json?{query}"
    request = Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "Droplet-Risk-Atlas-live-check/1.0",
        },
    )

    try:
        with urlopen(request, timeout=args.timeout) as response:
            molecule = json.load(response)
    except HTTPError as error:
        print(f"ERROR: ChEMBL returned HTTP {error.code} for {chembl_id}.", file=sys.stderr)
        return 3
    except URLError as error:
        print(f"ERROR: could not reach ChEMBL: {error.reason}", file=sys.stderr)
        return 3

    properties = molecule.get("molecule_properties") or {}
    values = {name: number(properties.get(api_name)) for name, api_name in DESCRIPTORS.items()}
    missing = [name for name, value in values.items() if value is None]

    output = {
        "molecule_chembl_id": molecule.get("molecule_chembl_id"),
        "pref_name": molecule.get("pref_name"),
        "descriptor_values": values,
        "missing_descriptors": missing,
        "parent_chembl_id": (molecule.get("molecule_hierarchy") or {}).get("parent_chembl_id"),
    }

    if not missing:
        score = (
            0.70 * values["alogp"]
            - 0.012 * values["psa"]
            - 0.18 * values["hbd"]
            - 0.06 * values["hba"]
            - 0.0010 * values["full_mwt"]
            + 0.035 * values["rtb"]
            + 0.10 * values["aromatic_rings"]
        )
        output.update({
            "score": score,
            "index": sigmoid(-2.1 + score),
            "risk_band": band(score),
        })

    print(json.dumps(output, indent=2, ensure_ascii=False))
    if missing:
        print("The API record is valid but cannot be scored because descriptors are missing.", file=sys.stderr)
        return 4
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
