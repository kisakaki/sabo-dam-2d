const { Sim, defaultParams } = require('../src/sabo2d-solver.js');
const P = defaultParams(), L = 400, B = 50;
P.bend = { mode: 'table', table: [[0, 0, B], [100, 15, B], [168, 0, 22], [200, -10, 22], [280, -20, 80], [L, 0, B]] };
P.dams = [{ x: 184, H: 4, b: 12 }, { x: 330, H: 2, b: 0 }];
P.hydro.duration = 2400; P.hydro.tPeak = 900;
const s = new Sim(P);
let w = []; for (let i = 0; i < s.nx; i += 20) w.push(`${((i + .5) * s.dx).toFixed(0)}:${s.Bx[i].toFixed(0)}`);
console.log('nx', s.nx, 'ny', s.ny, 'width', w.join(' '));
s.dams.forEach(d => console.log(`dam x=${d.x} th=${(d.th * 180 / Math.PI).toFixed(0)} cells=${d.cells.length} spill=[${d.spMin.toFixed(1)},${d.spMax.toFixed(1)}] ext=[${d.sMin.toFixed(1)},${d.sMax.toFixed(1)}]`));
let n = 0; while (s.phase !== 'done') { s.step(); n++; }
let bad = 0; for (let c = 0; c < s.N; c++) if (!isFinite(s.h[c]) || !isFinite(s.z[c])) bad++;
console.log(`steps ${n} bad ${bad} sedBal ${((s.sedIn - s.sedOut) / 0.6 - s.storage()).toExponential(1)} watBal ${(s.Vw0 + s.Vwin - s.Vwout - s.waterVolume()).toExponential(1)} byDam ${s.storageByDam().map(v => v.toFixed(0)).join(' ')}`);
