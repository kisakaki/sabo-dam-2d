const { Sim, defaultParams } = require('../src/sabo2d-solver.js');
const P = defaultParams(); P.bend = { mode: 'sine', amp: 30, wl: 200, x0: 40, x1: 400 };
P.dams = [{ x: 90, H: 1.5, b: 0 }, { x: 140, H: 1.5, b: 12 }, { x: 190, H: 1.5, b: 0 }, { x: 240, H: 1.5, b: 0 }];
P.hydro.duration = 3600; P.hydro.tPeak = 1200;
const s = new Sim(P);
s.dams.forEach(d => console.log(`dam x=${d.x} angle=${(d.th * 180 / Math.PI).toFixed(1)}° cells=${d.cells.length} spill=[${d.spMin.toFixed(1)},${d.spMax.toFixed(1)}] ext=[${d.sMin.toFixed(1)},${d.sMax.toFixed(1)}]`));
console.log('warn', s.warn);
// 漏れの確認：各堰堤の上流側と下流側の floor セルが 4 近傍で直接つながっていないか
let leaks = 0;
for (let c = 0; c < s.N; c++) { if (s.kind[c] === 2 || !s.fl[c]) continue; const i = c % s.nx; for (const o of [1, s.nx]) { const c2 = c + o; if (c2 >= s.N || (o === 1 && i === s.nx - 1)) continue; if (s.kind[c2] === 2) continue; if (Math.abs(s.region[c] - s.region[c2]) >= 1 && s.fl[c2]) { const d = s.dams[Math.min(s.region[c], s.region[c2])]; if (d && s.z[c] < d.crest && s.z[c2] < d.crest) leaks++; } } }
console.log('leaky faces below crest:', leaks);
let n = 0; while (s.phase !== 'done' && n < 2e6) { s.step(); n++; }
let bad = 0; for (let c = 0; c < s.N; c++) if (!isFinite(s.h[c]) || !isFinite(s.z[c])) bad++;
console.log(`steps ${n} bad ${bad} sedBal ${((s.sedIn - s.sedOut) / 0.6 - s.storage()).toExponential(1)} watBal ${(s.Vw0 + s.Vwin - s.Vwout - s.waterVolume()).toExponential(1)} byDam ${s.storageByDam().map(v => v.toFixed(0)).join(' ')} sum=${s.storageByDam().reduce((a, b) => a + b).toFixed(0)} stor=${s.storage().toFixed(0)}`);
