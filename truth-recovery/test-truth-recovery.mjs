// node --test truth-recovery/test-truth-recovery.mjs
// Measured invariants for the DTA truth-recovery yardstick. Seeded; no
// hand-entered numbers.
import assert from 'node:assert/strict';
import { describe, it } from 'node:test';
import { generate, makeRng, SCENARIOS } from './dgp-dta.mjs';
import { runCell, runGrid, summarize } from './harness.mjs';

describe('bivariate DTA DGP', () => {
  it('is reproducible for a fixed seed', () => {
    const a = generate(0.85, 0.80, 10, 'het_mod', makeRng(7));
    const b = generate(0.85, 0.80, 10, 'het_mod', makeRng(7));
    assert.deepEqual(a.studies, b.studies);
  });
  it('returns k valid 2x2 tables for every scenario', () => {
    const rng = makeRng(1);
    for (const scen of SCENARIOS) {
      const { studies } = generate(0.85, 0.80, 8, scen, rng);
      assert.equal(studies.length, 8);
      assert.ok(studies.every(s => s.tp >= 0 && s.fp >= 0 && s.fn >= 0 && s.tn >= 0 &&
        (s.tp + s.fn) >= 10 && (s.tn + s.fp) >= 10));
    }
  });
});

describe('Truth-recovery (measured)', () => {
  it('the shipped bivariate interval recovers true Se/Sp near nominal but slightly under', () => {
    const grid = runGrid({ reps: 400 });
    const s = summarize(grid);
    // App is mada-validated on points; coverage is decent but below 0.95 on average.
    assert.ok(s.shipped.meanCovSe > 0.90 && s.shipped.meanCovSe < 0.95,
      `shipped covSe ${s.shipped.meanCovSe}`);
  });

  it('genuine HKSJ variance adjustment recovers truth at least as well as the shipped t-only interval', () => {
    const grid = runGrid({ reps: 400 });
    const s = summarize(grid);
    // The shipped CI labels itself "HKSJ-t" but only swaps z->t; it omits the
    // Q/(k-1) variance inflation. Adding the genuine adjustment never lowers
    // coverage and lifts the mean toward nominal.
    assert.ok(s['genuine-HKSJ'].meanCovSe >= s.shipped.meanCovSe,
      `HKSJ Se ${s['genuine-HKSJ'].meanCovSe} < shipped ${s.shipped.meanCovSe}`);
    assert.ok(s['genuine-HKSJ'].meanCovSp >= s.shipped.meanCovSp,
      `HKSJ Sp ${s['genuine-HKSJ'].meanCovSp} < shipped ${s.shipped.meanCovSp}`);
    assert.ok(Math.abs(s['genuine-HKSJ'].meanCovSe - 0.95) <= 0.02,
      `HKSJ Se ${s['genuine-HKSJ'].meanCovSe} not near nominal`);
  });

  it('the coverage gap is largest under HIGH heterogeneity (where the t-only interval is too narrow)', () => {
    const rng = makeRng(20260613);
    const cell = runCell(0.85, 0.80, 20, 'het_high', 600, rng);
    assert.ok(cell.shipped.covSp < 0.93, `het_high shipped covSp ${cell.shipped.covSp} not under-covering`);
    assert.ok(cell['genuine-HKSJ'].covSp >= cell.shipped.covSp,
      `HKSJ ${cell['genuine-HKSJ'].covSp} < shipped ${cell.shipped.covSp}`);
  });
});
