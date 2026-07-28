import assert from "node:assert/strict";
import {
  riskBand,
  scoreChemblMolecule,
  scoreDescriptors,
} from "../site/assets/score.js";

const aspirin = {
  molecule_chembl_id: "CHEMBL25",
  molecule_properties: {
    alogp: "1.31",
    aromatic_rings: 1,
    full_mwt: "180.16",
    hba: 3,
    hbd: 1,
    psa: "63.60",
    rtb: 2,
  },
};

const result = scoreChemblMolecule(aspirin);
assert.equal(result.scored, true);
assert.ok(Math.abs(result.score - (-0.21636)) < 1e-10);
assert.ok(Math.abs(result.index - 0.08977706591742973) < 1e-12);
assert.equal(result.band, "leaky");

const missing = scoreDescriptors({
  logp: 1,
  psa: 30,
  hbd: 1,
  hba: 2,
  mw: 200,
  rtb: 2,
  aromaticRings: null,
});
assert.equal(missing.scored, false);
assert.deepEqual(missing.missing, ["aromaticRings"]);
assert.equal(missing.band, "unscored");

// The archived Python pipeline rounds the score to 10 decimal places before
// assigning a band. Test values safely on both sides of each boundary.
assert.equal(riskBand(-0.2563230878), "lower");
assert.equal(riskBand(-0.2563230876), "leaky");
assert.equal(riskBand(0.1453605433), "leaky");
assert.equal(riskBand(0.1453605435), "high");
assert.equal(riskBand(2.0999999999), "high");
assert.equal(riskBand(2.1), "very_high");

console.log("Score model tests passed.");
