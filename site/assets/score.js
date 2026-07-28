export const SCORE_MODEL = Object.freeze({
  version: "droplet-risk-unified-v2",
  coefficients: Object.freeze({
    logp: 0.70,
    psa: -0.012,
    hbd: -0.18,
    hba: -0.06,
    mw: -0.0010,
    rtb: 0.035,
    aromaticRings: 0.10,
  }),
  cutoffs: Object.freeze({
    leaky: -0.2563230876636397,
    high: 0.1453605434241663,
    veryHigh: 2.1,
  }),
  indexOffset: -2.1,
});

export const REQUIRED_DESCRIPTOR_LABELS = Object.freeze({
  logp: "ALogP",
  psa: "PSA",
  hbd: "HBD",
  hba: "HBA",
  mw: "full molecular weight",
  rtb: "rotatable bonds",
  aromaticRings: "aromatic rings",
});

function finiteNumber(value) {
  if (value === null || value === undefined || value === "") return null;
  const number = Number(value);
  return Number.isFinite(number) ? number : null;
}

export function sigmoid(value) {
  if (!Number.isFinite(value)) return null;
  if (value >= 0) {
    const z = Math.exp(-value);
    return 1 / (1 + z);
  }
  const z = Math.exp(value);
  return z / (1 + z);
}

export function riskBand(score) {
  if (!Number.isFinite(score)) return "unscored";
  // The archived Python pipeline rounds to 10 decimal places before assigning
  // a band. Reproduce that behavior so live calculations follow the same code.
  const comparisonScore = Number(score.toFixed(10));
  if (comparisonScore < SCORE_MODEL.cutoffs.leaky) return "lower";
  if (comparisonScore < SCORE_MODEL.cutoffs.high) return "leaky";
  if (comparisonScore < SCORE_MODEL.cutoffs.veryHigh) return "high";
  return "very_high";
}

export function riskInterpretation(band) {
  switch (band) {
    case "lower":
      return "Below the Payne-anchored leaky proxy cutoff; lower partition or material-sink concern under this descriptor-only screen.";
    case "leaky":
      return "Above the Payne-anchored leaky cutoff but below the high-transfer cutoff; prioritize controls when the molecule is assay-critical.";
    case "high":
      return "Above the Payne-anchored high-transfer proxy cutoff; screen for leakage, oil partition, interfacial transport, and material sinks.";
    case "very_high":
      return "Very-high proxy score with an index at or above 0.50; treat as a red-list candidate until experimentally cleared in the exact droplet system.";
    default:
      return "One or more descriptors required by the unified proxy are missing.";
  }
}

export function scoreDescriptors(input) {
  const descriptors = {
    logp: finiteNumber(input?.logp),
    psa: finiteNumber(input?.psa),
    hbd: finiteNumber(input?.hbd),
    hba: finiteNumber(input?.hba),
    mw: finiteNumber(input?.mw),
    rtb: finiteNumber(input?.rtb),
    aromaticRings: finiteNumber(input?.aromaticRings),
  };

  const missing = Object.entries(descriptors)
    .filter(([, value]) => value === null)
    .map(([key]) => key);

  if (missing.length) {
    return {
      scored: false,
      score: null,
      index: null,
      band: "unscored",
      descriptors,
      missing,
    };
  }

  const c = SCORE_MODEL.coefficients;
  const score =
    c.logp * descriptors.logp +
    c.psa * descriptors.psa +
    c.hbd * descriptors.hbd +
    c.hba * descriptors.hba +
    c.mw * descriptors.mw +
    c.rtb * descriptors.rtb +
    c.aromaticRings * descriptors.aromaticRings;

  return {
    scored: true,
    score,
    index: sigmoid(SCORE_MODEL.indexOffset + score),
    band: riskBand(score),
    descriptors,
    missing: [],
  };
}

export function descriptorsFromChembl(molecule) {
  const properties = molecule?.molecule_properties ?? {};
  return {
    logp: finiteNumber(properties.alogp),
    psa: finiteNumber(properties.psa),
    hbd: finiteNumber(properties.hbd),
    hba: finiteNumber(properties.hba),
    mw: finiteNumber(properties.full_mwt),
    rtb: finiteNumber(properties.rtb),
    aromaticRings: finiteNumber(properties.aromatic_rings),
  };
}

export function scoreChemblMolecule(molecule) {
  return scoreDescriptors(descriptorsFromChembl(molecule));
}

export function missingDescriptorLabels(missing) {
  return missing.map((key) => REQUIRED_DESCRIPTOR_LABELS[key] ?? key);
}
