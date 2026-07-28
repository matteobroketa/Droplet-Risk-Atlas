import {
  SCORE_MODEL,
  missingDescriptorLabels,
  riskInterpretation,
  scoreChemblMolecule,
} from "./score.js";

const CONFIG = Object.freeze({
  lookupUrl: "./data/lookup-16k.compact.json",
  summaryUrl: "./data/atlas-summary.json",
  chemblApiBase: "https://www.ebi.ac.uk/chembl/api/data",
  localRenderLimit: 350,
  chemblResultLimit: 12,
  requestTimeoutMs: 15000,
  chemblSearchFields: [
    "molecule_chembl_id",
    "pref_name",
    "molecule_type",
    "max_phase",
    "therapeutic_flag",
  ],
  chemblDetailFields: [
    "molecule_chembl_id",
    "pref_name",
    "molecule_type",
    "max_phase",
    "therapeutic_flag",
    "molecule_properties",
    "molecule_hierarchy",
  ],
});

const state = {
  records: [],
  filtered: [],
  localQueryMatchCount: 0,
  selectedLocalIndex: 0,
  query: "",
  dataset: "all",
  risk: "all",
  category: "all",
  sort: "score_desc",
  summary: null,
  chemblResults: [],
  chemblAbortController: null,
  chemblDetailCache: new Map(),
};

const $ = (id) => document.getElementById(id);

function escapeHtml(value) {
  return String(value ?? "").replace(/[&<>"']/g, (character) => ({
    "&": "&amp;",
    "<": "&lt;",
    ">": "&gt;",
    '"': "&quot;",
    "'": "&#039;",
  })[character]);
}

function normalizeText(value) {
  return String(value ?? "")
    .normalize("NFKD")
    .replace(/[\u0300-\u036f]/g, "")
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, " ")
    .trim();
}

function formatNumber(value, digits = 3) {
  if (value === null || value === undefined || !Number.isFinite(Number(value))) return "—";
  return Number(value).toFixed(digits).replace(/\.?0+$/, "");
}

function formatInteger(value) {
  if (value === null || value === undefined || !Number.isFinite(Number(value))) return "—";
  return Number(value).toLocaleString();
}

function debounce(callback, delayMs) {
  let timer = null;
  return (...args) => {
    window.clearTimeout(timer);
    timer = window.setTimeout(() => callback(...args), delayMs);
  };
}

async function fetchJson(url, options = {}) {
  const response = await fetch(url, {
    cache: "default",
    headers: { Accept: "application/json" },
    ...options,
  });
  if (!response.ok) {
    const error = new Error(`${response.status} ${response.statusText}`);
    error.status = response.status;
    throw error;
  }
  return response.json();
}

function decodeCompactLookup(payload) {
  if (!payload || !Array.isArray(payload.columns) || !Array.isArray(payload.rows)) {
    throw new Error("The compact lookup file has an invalid schema.");
  }
  const columns = payload.columns;
  const dictionaries = payload.dictionaries ?? {};
  return payload.rows.map((row) => {
    const record = {};
    for (let index = 0; index < columns.length; index += 1) {
      const column = columns[index];
      let value = row[index];
      if (value !== null && value !== undefined && Array.isArray(dictionaries[column])) {
        value = dictionaries[column][value];
      }
      if (value !== null && value !== undefined && value !== "") record[column] = value;
    }
    record.searchText = normalizeText([
      record.name,
      record.id,
      record.inchikey,
      record.dataset,
      record.category,
      record.subclass,
      record.rawClass,
    ].filter(Boolean).join(" "));
    return record;
  });
}

function uniqueSorted(values) {
  return [...new Set(values.filter(Boolean))].sort((a, b) => String(a).localeCompare(String(b)));
}

function fillSelect(select, values, allLabel) {
  select.innerHTML = ["all", ...values]
    .map((value) => `<option value="${escapeHtml(value)}">${escapeHtml(value === "all" ? allLabel : value)}</option>`)
    .join("");
}

function buildFilters() {
  fillSelect($("dataset-filter"), uniqueSorted(state.records.map((record) => record.dataset)), "All datasets");
  fillSelect($("category-filter"), uniqueSorted(state.records.map((record) => record.category)), "All categories");

  const risks = ["all", "lower", "leaky", "high", "very_high", "unscored"];
  $("risk-filter").innerHTML = risks
    .map((value) => `<option value="${value}">${value === "all" ? "All risk bands" : value.replace("_", " ")}</option>`)
    .join("");
}

function compareRecords(a, b) {
  switch (state.sort) {
    case "score_asc":
      return (a.score ?? Number.POSITIVE_INFINITY) - (b.score ?? Number.POSITIVE_INFINITY);
    case "name_asc":
      return String(a.name ?? "").localeCompare(String(b.name ?? ""));
    case "dataset_asc":
      return String(a.dataset ?? "").localeCompare(String(b.dataset ?? ""));
    case "score_desc":
    default:
      return (b.score ?? Number.NEGATIVE_INFINITY) - (a.score ?? Number.NEGATIVE_INFINITY);
  }
}

function filterLocalRecords() {
  const tokens = normalizeText(state.query).split(/\s+/).filter(Boolean);
  const queryMatches = state.records.filter((record) => (
    !tokens.length || tokens.every((token) => record.searchText.includes(token))
  ));
  state.localQueryMatchCount = queryMatches.length;
  state.filtered = queryMatches
    .filter((record) => {
      if (state.dataset !== "all" && record.dataset !== state.dataset) return false;
      if (state.risk !== "all" && record.riskClass !== state.risk) return false;
      if (state.category !== "all" && record.category !== state.category) return false;
      return true;
    })
    .sort(compareRecords);

  if (state.selectedLocalIndex >= state.filtered.length) state.selectedLocalIndex = 0;
}

function riskBadge(band) {
  const safeBand = ["lower", "leaky", "high", "very_high", "unscored"].includes(band) ? band : "unscored";
  return `<span class="badge badge-${safeBand}">${escapeHtml(safeBand.replace("_", " "))}</span>`;
}

function localRowHtml(record, index) {
  return `<tr data-local-index="${index}" class="${index === state.selectedLocalIndex ? "selected" : ""}">
    <td><strong>${escapeHtml(record.name ?? record.id ?? "Unlabeled record")}</strong><div class="muted small">${escapeHtml(record.dataset ?? "")}</div></td>
    <td>${riskBadge(record.riskClass)}</td>
    <td>${formatNumber(record.score)}</td>
    <td>${formatNumber(record.index)}</td>
    <td>${formatNumber(record.logp)}</td>
    <td>${formatNumber(record.psa, 1)}</td>
    <td>${formatNumber(record.mw, 1)}</td>
    <td>${escapeHtml(record.id ?? "—")}</td>
  </tr>`;
}

function renderLocalResults() {
  filterLocalRecords();
  const visibleRows = state.filtered.slice(0, CONFIG.localRenderLimit);
  const filterSuffix = state.localQueryMatchCount !== state.filtered.length
    ? ` (${state.localQueryMatchCount.toLocaleString()} before active filters)`
    : "";
  $("local-meta").textContent = `${state.filtered.length.toLocaleString()} local match${state.filtered.length === 1 ? "" : "es"}${filterSuffix}; showing ${visibleRows.length.toLocaleString()}.`;
  $("local-results-body").innerHTML = visibleRows.map(localRowHtml).join("");
  $("shown-count").textContent = state.filtered.length.toLocaleString();

  document.querySelectorAll("[data-local-index]").forEach((row) => {
    row.addEventListener("click", () => {
      state.selectedLocalIndex = Number(row.dataset.localIndex);
      renderLocalResults();
      renderLocalDetail();
    });
  });

  renderFallbackPrompt();
}

function detailCell(label, value) {
  return `<div class="detail-cell"><div class="detail-label">${escapeHtml(label)}</div><div class="detail-value">${escapeHtml(value ?? "—")}</div></div>`;
}

function renderLocalDetail() {
  const record = state.filtered[state.selectedLocalIndex];
  if (!record) {
    $("record-detail").innerHTML = `<div class="empty-state">No local record selected.</div>`;
    return;
  }

  $("record-detail").innerHTML = `
    <div class="detail-heading">
      <div><h2>${escapeHtml(record.name ?? record.id ?? "Unlabeled record")}</h2><p>${escapeHtml(record.dataset ?? "Local atlas")}</p></div>
      ${riskBadge(record.riskClass)}
    </div>
    <div class="detail-grid">
      ${detailCell("Identifier", record.id)}
      ${detailCell("InChIKey", record.inchikey)}
      ${detailCell("Score", formatNumber(record.score, 5))}
      ${detailCell("Index", formatNumber(record.index, 6))}
      ${detailCell("ALogP", formatNumber(record.logp))}
      ${detailCell("PSA", formatNumber(record.psa, 2))}
      ${detailCell("HBD", formatNumber(record.hbd, 0))}
      ${detailCell("HBA", formatNumber(record.hba, 0))}
      ${detailCell("Molecular weight", formatNumber(record.mw, 3))}
      ${detailCell("Rotatable bonds", formatNumber(record.rtb, 0))}
      ${detailCell("Aromatic rings", formatNumber(record.aromaticRings, 0))}
      ${detailCell("Maximum phase", record.maxPhase)}
      ${detailCell("Observed transfer", formatNumber(record.observed, 6))}
    </div>
    <div class="detail-section"><h3>Interpretation</h3><p>${escapeHtml(record.interpretation ?? riskInterpretation(record.riskClass))}</p></div>
    ${record.pathway ? `<div class="detail-section"><h3>Pathway</h3><p>${escapeHtml(record.pathway)}</p></div>` : ""}
    ${record.note ? `<div class="detail-section"><h3>Note</h3><p>${escapeHtml(record.note)}</p></div>` : ""}
    ${record.confidence ? `<div class="detail-section"><h3>Confidence</h3><p>${escapeHtml(record.confidence)}</p></div>` : ""}
    ${record.source ? `<div class="detail-section"><h3>Source</h3><p>${escapeHtml(record.source)}</p></div>` : ""}
  `;
}

function renderFallbackPrompt() {
  const prompt = $("chembl-fallback");
  const query = state.query.trim();
  const shouldShow = query.length >= 2 && state.localQueryMatchCount === 0;
  prompt.hidden = !shouldShow;
  if (!shouldShow) return;
  $("chembl-fallback-query").textContent = query;
  $("chembl-search-button").textContent = `Search ChEMBL for “${query}”`;
}

function setChemblStatus(message, kind = "neutral") {
  const status = $("chembl-status");
  status.textContent = message;
  status.dataset.kind = kind;
}

function chemblDisplayName(molecule) {
  return molecule?.pref_name || molecule?.molecule_chembl_id || "Unnamed ChEMBL molecule";
}

function chemblSearchRank(molecule, normalizedQuery) {
  const id = normalizeText(molecule?.molecule_chembl_id);
  const preferred = normalizeText(molecule?.pref_name);
  let score = 0;
  if (id === normalizedQuery) score += 1000;
  if (preferred === normalizedQuery) score += 900;
  if (preferred.startsWith(normalizedQuery)) score += 400;
  if (preferred.includes(normalizedQuery)) score += 150;
  score += Number(molecule?.max_phase ?? 0) * 5;
  if (molecule?.therapeutic_flag) score += 10;
  score -= Number(molecule?._apiOrder ?? 0) / 1000;
  return score;
}

async function chemblRequest(url) {
  if (state.chemblAbortController) state.chemblAbortController.abort();
  const controller = new AbortController();
  state.chemblAbortController = controller;
  const timeout = window.setTimeout(() => controller.abort(), CONFIG.requestTimeoutMs);
  try {
    return await fetchJson(url, { signal: controller.signal, cache: "no-store" });
  } finally {
    window.clearTimeout(timeout);
    if (state.chemblAbortController === controller) state.chemblAbortController = null;
  }
}

function withOnly(url, fields) {
  const separator = url.includes("?") ? "&" : "?";
  return `${url}${separator}only=${encodeURIComponent(fields.join(","))}`;
}

async function searchChemblCandidates(query) {
  const base = `${CONFIG.chemblApiBase}/molecule/search.json?q=${encodeURIComponent(query)}&limit=${CONFIG.chemblResultLimit}`;
  try {
    return await chemblRequest(withOnly(base, CONFIG.chemblSearchFields));
  } catch (error) {
    if (error?.name === "AbortError") throw error;
    return chemblRequest(base);
  }
}

async function fetchChemblMolecule(moleculeId) {
  if (state.chemblDetailCache.has(moleculeId)) return state.chemblDetailCache.get(moleculeId);
  const base = `${CONFIG.chemblApiBase}/molecule/${encodeURIComponent(moleculeId)}.json`;
  let molecule;
  try {
    molecule = await chemblRequest(withOnly(base, CONFIG.chemblDetailFields));
  } catch (error) {
    if (error?.name === "AbortError") throw error;
    molecule = await chemblRequest(base);
  }
  state.chemblDetailCache.set(moleculeId, molecule);
  return molecule;
}

async function searchChembl() {
  const query = state.query.trim();
  if (query.length < 2) return;

  const button = $("chembl-search-button");
  button.disabled = true;
  setChemblStatus("Searching the live ChEMBL service…");
  $("chembl-results").innerHTML = "";

  try {
    if (/^CHEMBL\d+$/i.test(query)) {
      const id = query.toUpperCase();
      state.chemblResults = [{ molecule_chembl_id: id, pref_name: id, _apiOrder: 0 }];
      renderChemblResults();
      await loadChemblDetailById(id);
      return;
    }

    const payload = await searchChemblCandidates(query);
    const molecules = Array.isArray(payload?.molecules) ? payload.molecules : [];
    const seen = new Set();
    const normalizedQuery = normalizeText(query);
    state.chemblResults = molecules
      .map((molecule, index) => ({ ...molecule, _apiOrder: index }))
      .filter((molecule) => {
        const id = molecule?.molecule_chembl_id;
        if (!id || seen.has(id)) return false;
        seen.add(id);
        return true;
      })
      .sort((a, b) => chemblSearchRank(b, normalizedQuery) - chemblSearchRank(a, normalizedQuery));

    if (!state.chemblResults.length) {
      setChemblStatus("ChEMBL returned no matching molecules.", "warning");
      $("chembl-results").innerHTML = `<div class="empty-state">Try a ChEMBL ID, a preferred name, or a common synonym.</div>`;
      return;
    }

    setChemblStatus(
      `${state.chemblResults.length} ChEMBL candidate${state.chemblResults.length === 1 ? "" : "s"}. Select the intended molecule or form before calculating.`,
      "success",
    );
    renderChemblResults();
  } catch (error) {
    const message = error?.name === "AbortError"
      ? "The ChEMBL request timed out or was replaced by a newer request."
      : `ChEMBL search failed: ${error.message}`;
    setChemblStatus(message, "error");
    $("chembl-results").innerHTML = `<div class="empty-state">The local 16k atlas remains available. Retry the external search later.</div>`;
  } finally {
    button.disabled = false;
  }
}

function chemblResultHtml(molecule, index) {
  const moleculeId = molecule?.molecule_chembl_id;
  return `<button class="chembl-result" data-chembl-index="${index}" type="button">
    <span><strong>${escapeHtml(chemblDisplayName(molecule))}</strong><small>${escapeHtml(moleculeId)} · ${escapeHtml(molecule?.molecule_type ?? "molecule")} · phase ${escapeHtml(molecule?.max_phase ?? "—")}</small></span>
    <span class="chembl-score"><strong>Calculate</strong><small>Load exact molecule descriptors</small></span>
  </button>`;
}

function renderChemblResults() {
  $("chembl-results").innerHTML = state.chemblResults.map(chemblResultHtml).join("");
  document.querySelectorAll("[data-chembl-index]").forEach((button) => {
    button.addEventListener("click", () => {
      const molecule = state.chemblResults[Number(button.dataset.chemblIndex)];
      if (molecule?.molecule_chembl_id) loadChemblDetailById(molecule.molecule_chembl_id);
    });
  });
}

function descriptorGrid(result) {
  const d = result.descriptors;
  return `
    ${detailCell("ALogP", formatNumber(d.logp))}
    ${detailCell("PSA", formatNumber(d.psa, 2))}
    ${detailCell("HBD", formatNumber(d.hbd, 0))}
    ${detailCell("HBA", formatNumber(d.hba, 0))}
    ${detailCell("Full molecular weight", formatNumber(d.mw, 3))}
    ${detailCell("Rotatable bonds", formatNumber(d.rtb, 0))}
    ${detailCell("Aromatic rings", formatNumber(d.aromaticRings, 0))}
  `;
}

function liveResultDetailHtml(molecule, result, contextNote = "") {
  const moleculeId = molecule?.molecule_chembl_id;
  const parentId = molecule?.molecule_hierarchy?.parent_chembl_id;
  const missingText = result.scored ? "" : missingDescriptorLabels(result.missing).join(", ");
  const parentButton = !result.scored && parentId && parentId !== moleculeId
    ? `<button type="button" class="button secondary" id="score-parent-button" data-parent-id="${escapeHtml(parentId)}">Score parent ${escapeHtml(parentId)}</button>`
    : "";

  return `
    <div class="detail-heading">
      <div><h2>${escapeHtml(chemblDisplayName(molecule))}</h2><p>Live ChEMBL API · ${escapeHtml(moleculeId)}</p></div>
      ${riskBadge(result.band)}
    </div>
    <div class="live-warning">
      Calculated on demand from the current ChEMBL API. This is not an archived ChEMBL36 row and may change when ChEMBL descriptors or compound forms change.
    </div>
    ${contextNote ? `<div class="notice">${escapeHtml(contextNote)}</div>` : ""}
    <div class="detail-grid">
      ${detailCell("Score", formatNumber(result.score, 5))}
      ${detailCell("Index", formatNumber(result.index, 6))}
      ${detailCell("Risk band", result.band.replace("_", " "))}
      ${detailCell("Maximum phase", molecule?.max_phase)}
      ${detailCell("Molecule type", molecule?.molecule_type)}
      ${detailCell("Therapeutic flag", molecule?.therapeutic_flag)}
      ${descriptorGrid(result)}
    </div>
    <div class="detail-section"><h3>Interpretation</h3><p>${escapeHtml(riskInterpretation(result.band))}</p></div>
    ${!result.scored ? `<div class="detail-section"><h3>Why it could not be scored</h3><p>Missing: ${escapeHtml(missingText)}.</p>${parentButton}</div>` : ""}
    <div class="detail-section"><h3>Source record</h3><p><a href="${CONFIG.chemblApiBase}/molecule/${encodeURIComponent(moleculeId)}.json" target="_blank" rel="noopener noreferrer">Open the raw ChEMBL JSON response</a></p></div>
  `;
}

async function loadChemblDetailById(moleculeId) {
  $("record-detail").innerHTML = `<div class="empty-state">Loading ${escapeHtml(moleculeId)} descriptors from ChEMBL…</div>`;
  try {
    const molecule = await fetchChemblMolecule(moleculeId);
    const result = scoreChemblMolecule(molecule);
    $("record-detail").innerHTML = liveResultDetailHtml(molecule, result);
    setChemblStatus(
      result.scored
        ? `${moleculeId} calculated successfully from all seven live descriptors.`
        : `${moleculeId} could not be calculated because required descriptors are missing.`,
      result.scored ? "success" : "warning",
    );
    const parentButton = $("score-parent-button");
    if (parentButton) {
      parentButton.addEventListener("click", () => scoreParentMolecule(parentButton.dataset.parentId, molecule));
    }
  } catch (error) {
    const message = error?.name === "AbortError"
      ? "The ChEMBL detail request timed out or was replaced by a newer request."
      : `Unable to load ${moleculeId}: ${error.message}`;
    setChemblStatus(message, "error");
    $("record-detail").innerHTML = `<div class="empty-state">${escapeHtml(message)}</div>`;
  }
}

async function scoreParentMolecule(parentId, originalMolecule) {
  const button = $("score-parent-button");
  if (button) button.disabled = true;
  try {
    const parent = await fetchChemblMolecule(parentId);
    const result = scoreChemblMolecule(parent);
    $("record-detail").innerHTML = liveResultDetailHtml(
      parent,
      result,
      `The selected record ${originalMolecule?.molecule_chembl_id} lacked complete descriptors. This score uses its ChEMBL parent ${parentId}; it is not a score for the original salt or form.`,
    );
  } catch (error) {
    setChemblStatus(`Unable to load parent ${parentId}: ${error.message}`, "error");
    if (button) button.disabled = false;
  }
}

function exportFilteredCsv() {
  const columns = [
    "dataset", "name", "id", "riskClass", "score", "index", "logp", "psa", "hbd", "hba", "mw", "rtb", "aromaticRings", "maxPhase", "inchikey",
  ];
  const quote = (value) => `"${String(value ?? "").replace(/"/g, '""')}"`;
  const lines = [columns.join(",")];
  for (const record of state.filtered) {
    lines.push(columns.map((column) => quote(record[column])).join(","));
  }
  const blob = new Blob([lines.join("\n")], { type: "text/csv;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = "droplet-risk-local-results.csv";
  anchor.click();
  window.setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function clearSearch() {
  state.query = "";
  state.dataset = "all";
  state.risk = "all";
  state.category = "all";
  state.sort = "score_desc";
  state.chemblResults = [];
  $("search-input").value = "";
  $("dataset-filter").value = "all";
  $("risk-filter").value = "all";
  $("category-filter").value = "all";
  $("sort-filter").value = "score_desc";
  $("chembl-results").innerHTML = "";
  setChemblStatus("ChEMBL is queried only after you request an external search.");
  renderLocalResults();
  renderLocalDetail();
}

function renderSummary() {
  const chembl = state.summary?.chembl ?? {};
  $("instant-count").textContent = state.records.length.toLocaleString();
  $("full-row-count").textContent = formatInteger(chembl.rows);
  $("full-scored-count").textContent = formatInteger(chembl.scored);
  $("model-version").textContent = SCORE_MODEL.version;
}

async function initialize() {
  try {
    const [lookupPayload, summary] = await Promise.all([
      fetchJson(CONFIG.lookupUrl),
      fetchJson(CONFIG.summaryUrl),
    ]);
    state.records = decodeCompactLookup(lookupPayload);
    state.summary = summary;
    buildFilters();
    renderSummary();
    renderLocalResults();
    renderLocalDetail();
    $("loading-banner").hidden = true;
  } catch (error) {
    $("loading-banner").textContent = `Unable to load the local atlas: ${error.message}`;
    $("loading-banner").dataset.kind = "error";
    console.error(error);
  }
}

const onSearchInput = debounce((value) => {
  state.query = value;
  state.selectedLocalIndex = 0;
  state.chemblResults = [];
  $("chembl-results").innerHTML = "";
  setChemblStatus("ChEMBL is queried only after you request an external search.");
  renderLocalResults();
  renderLocalDetail();
}, 120);

$("search-input").addEventListener("input", (event) => onSearchInput(event.target.value));
$("clear-button").addEventListener("click", clearSearch);
$("export-button").addEventListener("click", exportFilteredCsv);
$("chembl-search-button").addEventListener("click", searchChembl);
$("dataset-filter").addEventListener("change", (event) => {
  state.dataset = event.target.value;
  state.selectedLocalIndex = 0;
  renderLocalResults();
  renderLocalDetail();
});
$("risk-filter").addEventListener("change", (event) => {
  state.risk = event.target.value;
  state.selectedLocalIndex = 0;
  renderLocalResults();
  renderLocalDetail();
});
$("category-filter").addEventListener("change", (event) => {
  state.category = event.target.value;
  state.selectedLocalIndex = 0;
  renderLocalResults();
  renderLocalDetail();
});
$("sort-filter").addEventListener("change", (event) => {
  state.sort = event.target.value;
  state.selectedLocalIndex = 0;
  renderLocalResults();
  renderLocalDetail();
});

document.querySelectorAll("[data-example-query]").forEach((button) => {
  button.addEventListener("click", () => {
    const query = button.dataset.exampleQuery;
    $("search-input").value = query;
    state.query = query;
    state.selectedLocalIndex = 0;
    renderLocalResults();
    renderLocalDetail();
  });
});

initialize();
