// engine.mjs -- pure DTA bivariate-pool core EXTRACTED VERBATIM from
// dta-pro-v3.7.html. invLogit, tQuantile, normalQuantile and
// improvedBivariatePool copied unchanged so the SAME math the app ships is
// what the harness measures. qnorm bound to self-contained normalQuantile.

const invLogit = x => 1 / (1 + Math.exp(-x));

function tQuantile(p, df) {
  if (df >= 120) return normalQuantile(p);
  if (df <= 0) return NaN;
  const z = normalQuantile(p);
  const g1 = (z * z * z + z) / 4;
  const g2 = (5 * Math.pow(z, 5) + 16 * Math.pow(z, 3) + 3 * z) / 96;
  const g3 = (3 * Math.pow(z, 7) + 19 * Math.pow(z, 5) + 17 * Math.pow(z, 3) - 15 * z) / 384;
  const g4 = (79 * Math.pow(z, 9) + 776 * Math.pow(z, 7) + 1482 * Math.pow(z, 5) - 1920 * Math.pow(z, 3) - 945 * z) / 92160;
  return z + g1/df + g2/(df*df) + g3/(df*df*df) + g4/(df*df*df*df);
}
function normalQuantile(p) {
  if (p <= 0) return -Infinity;
  if (p >= 1) return Infinity;
  if (p === 0.5) return 0;
  const a = [
    -3.969683028665376e+01, 2.209460984245205e+02,
    -2.759285104469687e+02, 1.383577518672690e+02,
    -3.066479806614716e+01, 2.506628277459239e+00
  ];
  const b = [
    -5.447609879822406e+01, 1.615858368580409e+02,
    -1.556989798598866e+02, 6.680131188771972e+01,
    -1.328068155288572e+01
  ];
  const c = [
    -7.784894002430293e-03, -3.223964580411365e-01,
    -2.400758277161838e+00, -2.549732539343734e+00,
    4.374664141464968e+00, 2.938163982698783e+00
  ];
  const d = [
    7.784695709041462e-03, 3.224671290700398e-01,
    2.445134137142996e+00, 3.754408661907416e+00
  ];
  const pLow = 0.02425;
  const pHigh = 1 - pLow;
  let q, r;
  if (p < pLow) {
    q = Math.sqrt(-2 * Math.log(p));
    return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
           ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1);
  } else if (p <= pHigh) {
    q = p - 0.5;
    r = q * q;
    return (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5])*q /
           (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1);
  } else {
    q = Math.sqrt(-2 * Math.log(1 - p));
    return -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
            ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1);
  }
}

const qnorm = normalQuantile;

function improvedBivariatePool(studies) {
  const n = studies.length;
  const data = studies.map(s => {
    const hasZero = s.tp === 0 || s.fp === 0 || s.fn === 0 || s.tn === 0;
    const cc = hasZero ? 0.5 : 0;
    const tp = s.tp + cc;
    const fp = s.fp + cc;
    const fn = s.fn + cc;
    const tn = s.tn + cc;
    const sens = tp / (tp + fn);
    const spec = tn / (tn + fp);
    const y1 = Math.log(sens / (1 - sens));
    const y2 = Math.log(spec / (1 - spec));
    const v1 = 1/tp + 1/fn;
    const v2 = 1/tn + 1/fp;
    return { sens, spec, y1, y2, v1, v2, hasZero };
  });
  const w1 = data.map(d => 1/d.v1);
  const w2 = data.map(d => 1/d.v2);
  const sumW1 = w1.reduce((a,b) => a+b, 0);
  const sumW2 = w2.reduce((a,b) => a+b, 0);
  const mu1_fe = data.reduce((s,d,i) => s + d.y1*w1[i], 0) / sumW1;
  const mu2_fe = data.reduce((s,d,i) => s + d.y2*w2[i], 0) / sumW2;
  const Q1 = data.reduce((s,d,i) => s + w1[i] * Math.pow(d.y1 - mu1_fe, 2), 0);
  const Q2 = data.reduce((s,d,i) => s + w2[i] * Math.pow(d.y2 - mu2_fe, 2), 0);
  const C1 = sumW1 - w1.reduce((s,w) => s + w*w, 0) / sumW1;
  const C2 = sumW2 - w2.reduce((s,w) => s + w*w, 0) / sumW2;
  const tau2_1 = Math.max(0, (Q1 - (n-1)) / C1);
  const tau2_2 = Math.max(0, (Q2 - (n-1)) / C2);
  const w1_re = data.map(d => 1/(d.v1 + tau2_1));
  const w2_re = data.map(d => 1/(d.v2 + tau2_2));
  const sumW1_re = w1_re.reduce((a,b) => a+b, 0);
  const sumW2_re = w2_re.reduce((a,b) => a+b, 0);
  const mu1 = data.reduce((s,d,i) => s + d.y1*w1_re[i], 0) / sumW1_re;
  const mu2 = data.reduce((s,d,i) => s + d.y2*w2_re[i], 0) / sumW2_re;
  const pooledSens = invLogit(mu1);
  const pooledSpec = invLogit(mu2);
  const seMu1 = 1 / Math.sqrt(sumW1_re);
  const seMu2 = 1 / Math.sqrt(sumW2_re);
  const k = n; // Number of studies for HKSJ adjustment
  const confLevelPool = (typeof document !== 'undefined' && document.getElementById('confLevel')) ? parseFloat(document.getElementById('confLevel').value) : 0.95;
  const alphaPool = 1 - confLevelPool;
  const critValue = k >= 30 ? qnorm(1 - alphaPool / 2) : (typeof tQuantile === 'function' ? tQuantile(1 - alphaPool / 2, Math.max(1, k - 2)) : qnorm(1 - alphaPool / 2));
  return {
    sens: pooledSens,
    spec: pooledSpec,
    sensCI: [invLogit(mu1 - critValue*seMu1), invLogit(mu1 + critValue*seMu1)],
    specCI: [invLogit(mu2 - critValue*seMu2), invLogit(mu2 + critValue*seMu2)],
    tau2_sens: tau2_1,
    tau2_spec: tau2_2,
    criticalValue: critValue,
    df: k - 2,
    ciMethod: k >= 30 ? 'z-normal' : 'HKSJ-t'
  };
}

export { improvedBivariatePool, invLogit, tQuantile, normalQuantile };
