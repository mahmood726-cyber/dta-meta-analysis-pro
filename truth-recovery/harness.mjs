// ============================================================
// harness.mjs -- Truth-recovery yardstick for dta-meta-analysis-pro.
//
// Wires the app's OWN improvedBivariatePool (copied verbatim into engine.mjs)
// to the known-truth bivariate DTA DGP and measures how often the pooled
// sensitivity / specificity CIs cover the TRUE Se / Sp.
//
// It also evaluates ONE candidate improvement. The shipped method labels its
// interval "HKSJ-t" but only swaps z -> t_{k-2}; it does NOT apply the HKSJ
// variance adjustment (the Q/(k-1) inflation with a max(1,.) floor). We add a
// genuine-HKSJ variant computed from the SAME RE weights and compare coverage.
//
// Truth-first: every number printed comes from seeded simulation here.
// Run:  node truth-recovery/harness.mjs --reps 400
// ============================================================

import { improvedBivariatePool, invLogit, tQuantile, normalQuantile } from './engine.mjs';
import { generate, makeRng, SCENARIOS } from './dgp-dta.mjs';

const BASE_SEED = 20260613;
const logit = (p) => Math.log(p / (1 - p));

// Recompute logit(Se)/logit(Sp) per study with the app's zero-cell CC, then
// build BOTH the shipped interval (ordinary RE SE x t_{k-2}) and a genuine-HKSJ
// interval (Q/(k-1) variance inflation, floored at 1). Point estimates identical.
function prep(studies) {
  const k = studies.length;
  const data = studies.map(s => {
    const hasZero = s.tp === 0 || s.fp === 0 || s.fn === 0 || s.tn === 0;
    const cc = hasZero ? 0.5 : 0;
    const tp = s.tp + cc, fp = s.fp + cc, fn = s.fn + cc, tn = s.tn + cc;
    const sens = tp / (tp + fn), spec = tn / (tn + fp);
    return { y1: Math.log(sens / (1 - sens)), y2: Math.log(spec / (1 - spec)),
             v1: 1 / tp + 1 / fn, v2: 1 / tn + 1 / fp };
  });
  return { k, data };
}

function reAxis(data, key, vkey, k) {
  const w = data.map(d => 1 / d[vkey]);
  const sumW = w.reduce((a, b) => a + b, 0);
  const muFE = data.reduce((s, d, i) => s + d[key] * w[i], 0) / sumW;
  const Q = data.reduce((s, d, i) => s + w[i] * (d[key] - muFE) ** 2, 0);
  const C = sumW - w.reduce((s, x) => s + x * x, 0) / sumW;
  const tau2 = Math.max(0, (Q - (k - 1)) / C);
  const wRE = data.map(d => 1 / (d[vkey] + tau2));
  const sumWRE = wRE.reduce((a, b) => a + b, 0);
  const mu = data.reduce((s, d, i) => s + d[key] * wRE[i], 0) / sumWRE;
  const seStd = 1 / Math.sqrt(sumWRE);
  const qRE = data.reduce((s, d, i) => s + wRE[i] * (d[key] - mu) ** 2, 0);
  const seHKSJ = Math.sqrt(Math.max(1, qRE / (k - 1)) / sumWRE);    // genuine HKSJ
  return { mu, seStd, seHKSJ };
}

// shipped: improvedBivariatePool (t_{k-2} crit on ordinary RE SE)
function methodShipped(studies) {
  const r = improvedBivariatePool(studies);
  return { sens: r.sens, spec: r.spec, sensCI: r.sensCI, specCI: r.specCI };
}
// genuine HKSJ variant (same point estimate; HKSJ-adjusted SE x t_{k-2})
function methodHKSJ(studies) {
  const { k, data } = prep(studies);
  const a = reAxis(data, 'y1', 'v1', k);
  const b = reAxis(data, 'y2', 'v2', k);
  const crit = k >= 30 ? normalQuantile(0.975) : tQuantile(0.975, Math.max(1, k - 2));
  return {
    sens: invLogit(a.mu), spec: invLogit(b.mu),
    sensCI: [invLogit(a.mu - crit * a.seHKSJ), invLogit(a.mu + crit * a.seHKSJ)],
    specCI: [invLogit(b.mu - crit * b.seHKSJ), invLogit(b.mu + crit * b.seHKSJ)],
  };
}

const METHODS = { shipped: methodShipped, 'genuine-HKSJ': methodHKSJ };
const mean = (a) => a.reduce((x, y) => x + y, 0) / a.length;

export function runCell(SeTrue, SpTrue, k, scenario, reps, rng) {
  const acc = {};
  for (const n of Object.keys(METHODS))
    acc[n] = { covSe: 0, covSp: 0, biasSe: 0, biasSp: 0, wSe: 0, wSp: 0, n: 0 };
  for (let r = 0; r < reps; r++) {
    const { studies } = generate(SeTrue, SpTrue, k, scenario, rng);
    for (const [name, fn] of Object.entries(METHODS)) {
      let o; try { o = fn(studies); } catch { continue; }
      if (!o || !isFinite(o.sens) || !isFinite(o.sensCI[0]) || !isFinite(o.specCI[0])) continue;
      const a = acc[name];
      a.n++;
      a.biasSe += o.sens - SeTrue; a.biasSp += o.spec - SpTrue;
      a.wSe += o.sensCI[1] - o.sensCI[0]; a.wSp += o.specCI[1] - o.specCI[0];
      if (o.sensCI[0] <= SeTrue && SeTrue <= o.sensCI[1]) a.covSe++;
      if (o.specCI[0] <= SpTrue && SpTrue <= o.specCI[1]) a.covSp++;
    }
  }
  const res = {};
  for (const [name, a] of Object.entries(acc)) {
    res[name] = {
      n: a.n,
      covSe: a.n ? +(a.covSe / a.n).toFixed(4) : null,
      covSp: a.n ? +(a.covSp / a.n).toFixed(4) : null,
      biasSe: a.n ? +(a.biasSe / a.n).toFixed(4) : null,
      biasSp: a.n ? +(a.biasSp / a.n).toFixed(4) : null,
      widthSe: a.n ? +(a.wSe / a.n).toFixed(4) : null,
    };
  }
  return res;
}

export function runGrid({ reps = 400, ks = [5, 10, 20], Se = 0.85, Sp = 0.80,
                          scenarios = SCENARIOS } = {}) {
  const rng = makeRng(BASE_SEED);
  const grid = [];
  for (const scen of scenarios)
    for (const k of ks) grid.push({ scen, k, results: runCell(Se, Sp, k, scen, reps, rng) });
  return grid;
}

export function summarize(grid, filter = () => true) {
  const names = Object.keys(grid[0].results);
  const out = {};
  for (const name of names) {
    const cSe = [], cSp = [];
    for (const c of grid) {
      if (!filter(c)) continue;
      if (c.results[name].covSe != null) cSe.push(c.results[name].covSe);
      if (c.results[name].covSp != null) cSp.push(c.results[name].covSp);
    }
    out[name] = { meanCovSe: +mean(cSe).toFixed(4), meanCovSp: +mean(cSp).toFixed(4) };
  }
  return out;
}

const isMain = process.argv[1]?.endsWith('harness.mjs');
if (isMain) {
  const i = process.argv.indexOf('--reps');
  const reps = i >= 0 ? Number(process.argv[i + 1]) : 400;
  const t0 = Date.now();
  const grid = runGrid({ reps });
  const s = summarize(grid);
  console.log(`\n# Truth-recovery yardstick -- dta-meta-analysis-pro`);
  console.log(`reps=${reps}/cell  Se=0.85 Sp=0.80  seed=${BASE_SEED}\n`);
  console.log('## Mean coverage of TRUE Se / Sp over all heterogeneity cells\n');
  console.log('method         meanCovSe  meanCovSp');
  for (const [n, v] of Object.entries(s))
    console.log(n.padEnd(13), String(v.meanCovSe).padStart(9), String(v.meanCovSp).padStart(10));
  console.log('\n## Per-cell coverage of true Se (shipped vs genuine-HKSJ)\n');
  console.log('scenario    k   ship.covSe  HKSJ.covSe   ship.covSp  HKSJ.covSp');
  for (const c of grid) {
    console.log(c.scen.padEnd(11), String(c.k).padStart(2),
      String(c.results.shipped.covSe).padStart(11),
      String(c.results['genuine-HKSJ'].covSe).padStart(11),
      String(c.results.shipped.covSp).padStart(12),
      String(c.results['genuine-HKSJ'].covSp).padStart(11));
  }
  console.log(`\n(${(Date.now() - t0) / 1000}s)`);
}
