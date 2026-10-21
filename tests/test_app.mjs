import test from "node:test";
import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import vm from "node:vm";
import * as score from "../site/assets/score.js";

const source = (await readFile(new URL("../site/assets/app.js", import.meta.url), "utf8"))
  .replace(/^import\s*\{[\s\S]*?\}\s*from\s*"\.\/score\.js";\s*/, "");

async function harness() {
  const elements = new Map();
  function element(id) {
    if (!elements.has(id)) elements.set(id, {
      innerHTML: "", textContent: "", value: "", hidden: false, disabled: false,
      dataset: {}, listeners: {},
      addEventListener(type, listener) { this.listeners[type] = listener; },
    });
    return elements.get(id);
  }
  const timers = new Map();
  let timerId = 0;
  const requests = [];
  const downloads = [];
  const example = element("example");
  example.dataset.exampleQuery = "Aspirin";
  const context = vm.createContext({
    ...score, console, AbortController, Blob,
    URL: { createObjectURL(blob) { downloads.push(blob); return "blob:test"; }, revokeObjectURL() {} },
    window: {
      setTimeout(callback, delay) { timers.set(++timerId, { callback, delay }); return timerId; },
      clearTimeout(id) { timers.delete(id); },
    },
    document: {
      getElementById: element,
      createElement: () => ({ click() {} }),
      querySelectorAll: selector => selector === "[data-example-query]" ? [example] : [],
    },
    fetch: async url => {
      if (url.includes("lookup-16k")) return { ok: true, json: async () => ({
        columns: ["name", "id", "dataset", "riskClass", "score", "hbd", "hba"],
        rows: [["Aspirin", "CHEMBL25", "Local", "leaky", -0.21636, 0, 10]],
      }) };
      if (url.includes("atlas-summary")) return { ok: true, json: async () => ({}) };
      // Deliberately allow a late response even after abort, as cached or already
      // decoded responses can finish after a newer user action.
      return new Promise(resolve => requests.push({ url, resolve }));
    },
  });
  vm.runInContext(source + "\nglobalThis.actions = {formatNumber, loadChemblDetailById, scoreParentMolecule};", context);
  await settle();
  return {
    element, timers, requests, downloads, actions: context.actions,
    click(id) { return element(id).listeners.click(); },
    input(value) {
      element("search-input").value = value;
      element("search-input").listeners.input({ target: element("search-input") });
    },
    flushInput() {
      for (const [id, timer] of [...timers]) {
        if (timer.delay === 120) { timers.delete(id); timer.callback(); }
      }
    },
    respond(index, molecule) {
      requests[index].resolve({ ok: true, json: async () => molecule });
    },
  };
}

async function settle() {
  for (let i = 0; i < 15; i += 1) await Promise.resolve();
}

const molecule = id => ({
  molecule_chembl_id: id, pref_name: "Live " + id,
  molecule_properties: { alogp: 1.31, aromatic_rings: 1, full_mwt: 180.16, hba: 3, hbd: 1, psa: 63.6, rtb: 2 },
});

test("integer descriptors retain zeros and zero remains visible", async () => {
  const app = await harness();
  for (const [value, digits, expected] of [[0, 0, "0"], [10, 0, "10"], [20, 0, "20"], [100, 0, "100"], [10, 3, "10"], [1.2, 3, "1.2"], [null, 0, "—"]]) {
    assert.equal(app.actions.formatNumber(value, digits), expected);
  }
  assert.match(app.element("record-detail").innerHTML, /HBA[\s\S]*?detail-value">10</);
});

test("Clear cancels a pending input instead of restoring the old query", async () => {
  const app = await harness();
  app.input("not-present");
  app.click("clear-button");
  app.flushInput();
  assert.equal(app.element("search-input").value, "");
  assert.equal(app.element("shown-count").textContent, "1");
  assert.equal(app.element("chembl-fallback").hidden, true);
});

test("example queries cancel older pending input", async () => {
  const app = await harness();
  app.input("not-present");
  app.click("example");
  app.flushInput();
  assert.equal(app.element("search-input").value, "Aspirin");
  assert.equal(app.element("shown-count").textContent, "1");
});

test("late ChEMBL details cannot replace the local record after Clear", async () => {
  const app = await harness();
  const pending = app.actions.loadChemblDetailById("CHEMBL999");
  app.click("clear-button");
  app.respond(0, molecule("CHEMBL999"));
  await pending;
  assert.match(app.element("record-detail").innerHTML, /Aspirin/);
  assert.doesNotMatch(app.element("record-detail").innerHTML, /Live CHEMBL999/);
  assert.equal(app.element("chembl-status").textContent, "ChEMBL is queried only after you request an external search.");
});

test("a late ChEMBL candidate search cannot repopulate a cleared search", async () => {
  const app = await harness();
  app.input("not-present");
  app.flushInput();
  const pending = app.click("chembl-search-button");
  app.click("clear-button");
  app.respond(0, { molecules: [molecule("CHEMBL999")] });
  await pending;
  assert.equal(app.element("chembl-results").innerHTML, "");
  assert.equal(app.element("chembl-search-button").disabled, false);
});

test("a newer live molecule selection wins when responses arrive out of order", async () => {
  const app = await harness();
  const first = app.actions.loadChemblDetailById("CHEMBL111");
  const second = app.actions.loadChemblDetailById("CHEMBL222");
  app.respond(1, molecule("CHEMBL222"));
  await second;
  app.respond(0, molecule("CHEMBL111"));
  await first;
  assert.match(app.element("record-detail").innerHTML, /Live CHEMBL222/);
  assert.doesNotMatch(app.element("record-detail").innerHTML, /Live CHEMBL111/);
});


test("export immediately after typing uses the current query", async () => {
  const app = await harness();
  app.input("not-present");
  app.click("export-button");
  const csv = await app.downloads[0].text();
  assert.equal(csv.split("\n").length, 1);
  assert.doesNotMatch(csv, /Aspirin/);
});

test("live search immediately after typing uses the current query", async () => {
  const app = await harness();
  app.input("older missing");
  app.flushInput();
  app.input("newer missing");
  const pending = app.click("chembl-search-button");
  assert.match(app.requests[0].url, /q=newer%20missing/);
  app.respond(0, { molecules: [] });
  await pending;
});


test("obsolete failed requests do not retry or replace the current local view", async () => {
  const app = await harness();
  const pending = app.actions.loadChemblDetailById("CHEMBL999");
  app.click("clear-button");
  app.requests[0].resolve({ ok: false, status: 503, statusText: "Unavailable" });
  await pending;
  assert.equal(app.requests.length, 1);
  assert.match(app.element("record-detail").innerHTML, /Aspirin/);
  assert.equal(app.element("chembl-status").dataset.kind, "neutral");
});

test("late parent scoring cannot replace a newer local search", async () => {
  const app = await harness();
  const pending = app.actions.scoreParentMolecule("CHEMBL999", molecule("CHEMBL111"));
  app.input("Aspirin");
  app.flushInput();
  app.respond(0, molecule("CHEMBL999"));
  await pending;
  assert.match(app.element("record-detail").innerHTML, /Aspirin/);
  assert.doesNotMatch(app.element("record-detail").innerHTML, /Live CHEMBL999/);
});
