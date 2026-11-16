const { Sim, defaultParams } = require('../src/sabo2d-solver.js');
const PRE = {
  A: { L: 600, B: 60, slope: 1 / 18, crossSlope: 0.03, sideGrad: 1.0,
    bend: { mode: 'table', table: [[0,0,52],[70,6,58],[150,-6,66],[230,-14,70],[310,-4,62],[380,10,74],[460,16,68],[530,6,60],[600,0,56]] },
    valley: { chD: 0.7, chW: 14, chWl: 0, relief: 0.3, wallH: 16, wallRough: 3, fans: [[120,-1,34,3.2],[340,1,40,3.5],[505,-1,28,2.5]] },
    dams: [75,135,195,255,315,375,435,495,555].map(x => ({ x, H: 3, b: 0 })) },
  B: { L: 660, B: 44, slope: 1 / 28, crossSlope: 0.02, sideGrad: 1.0,
    bend: { mode: 'table', table: [[0,0,40],[70,10,38],[140,-6,42],[200,-2,30],[235,0,24],[270,4,36],[340,18,84],[410,24,110],[480,18,112],[540,8,70],[575,2,34],[610,-2,40],[660,-6,42]] },
    valley: { chD: 0.8, chW: 12, chWl: 0, relief: 0.35, wallH: 18, wallRough: 3, fans: [[110,1,30,3],[300,-1,36,3.5],[440,1,46,4]] },
    dams: [{ x: 235, H: 5, b: 16 }, { x: 577, H: 9, b: 22 }] }
};
if (process.env.FILL) { PRE.A.prefill = 0.5; PRE.B.prefill = 0.5; }
const which = process.argv[2] || 'A', dx = +(process.argv[3] || 2), tEnd = +(process.argv[4] || 0);
const P = Object.assign(defaultParams(), PRE[which]); P.dx = dx; P.hydro.duration = 3600; P.hydro.tPeak = 1200; P.hydro.Qpeak = 80;
const s = new Sim(P);
console.log(which, 'dx', dx, 'nx', s.nx, 'ny', s.ny, 'cells', s.N, 'ns', s.ns, 'warn', s.warn);
s.dams.forEach(d => console.log(`  dam x=${d.x.toFixed(0)} th=${(d.th*180/Math.PI).toFixed(0)} cells=${d.cells.length} b=${d.b} spill=[${d.spMin.toFixed(1)},${d.spMax.toFixed(1)}] ext=[${d.sMin.toFixed(1)},${d.sMax.toFixed(1)}]`));
let leaks = 0;
for (let c = 0; c < s.N; c++) { if (s.kind[c] === 2 || !s.fl[c]) continue; const i = c % s.nx; for (const o of [1, s.nx]) { const c2 = c + o; if (c2 >= s.N || (o === 1 && i === s.nx - 1)) continue; if (s.kind[c2] === 2) continue; if (s.region[c] !== s.region[c2] && s.fl[c2]) { const d = s.dams[Math.min(s.region[c], s.region[c2])]; if (d && s.z[c] < d.crest && s.z[c2] < d.crest) leaks++; } } }
console.log('leaks', leaks, 'floor cells', s.fl.reduce((a, b) => a + b, 0));
const t0 = Date.now(); let n = 0;
while (s.phase !== 'done' && (tEnd ? s.t < tEnd : n < 3000)) { s.step(); n++; }
let bad = 0; for (let c = 0; c < s.N; c++) if (!isFinite(s.h[c]) || !isFinite(s.z[c])) bad++;
console.log(`steps ${n} in ${(Date.now()-t0)/1000}s  t=${s.t.toFixed(0)} th=${s.th.toFixed(0)} dt=${s.dt.toFixed(4)} bad ${bad} sedBal ${((s.sedIn - s.sedOut) / 0.6 - s.storage()).toExponential(1)} watBal ${(s.Vw0 + s.Vwin - s.Vwout - s.waterVolume()).toExponential(1)} Qout ${s.Qout.toFixed(1)} byDam ${s.storageByDam().map(v => v.toFixed(0)).join(' ')}`);
module.exports = PRE;
