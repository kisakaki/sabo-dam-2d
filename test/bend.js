const { Sim, defaultParams } = require('../src/sabo2d-solver.js');
function run(P, tag) {
  const s = new Sim(P);
  let n = 0; while (s.phase !== 'done' && n < 2e6) { s.step(); n++; }
  let bad = 0; for (let c = 0; c < s.N; c++) if (!isFinite(s.h[c]) || !isFinite(s.z[c])) bad++;
  console.log(tag, `grid ${s.nx}x${s.ny} steps ${n} bad ${bad} sedBal ${((s.sedIn - s.sedOut) / (1 - P.porosity) - s.storage()).toExponential(1)} watBal ${(s.Vw0 + s.Vwin - s.Vwout - s.waterVolume()).toExponential(1)} Qout ${s.Qout.toFixed(1)} I ${s.I.toFixed(4)} warn [${s.warn.join('|')}] spill ${s.dams.map(d => d.j0 + '-' + d.j1).join(',')}`);
  return s;
}
let P = defaultParams(); P.bend = { mode: 'sine', amp: 30, wl: 200, x0: 40, x1: 400 }; P.dams = [{ x: 150, H: 1.2, b: 0 }, { x: 250, H: 1.2, b: 10 }]; P.hydro.duration = 3600; P.hydro.tPeak = 1200;
const sb = run(P, 'sine');
const rows = []; for (let j = sb.ny - 1; j >= 0; j--) { const r = []; for (let i = 0; i < sb.nx; i++) { const c = j * sb.nx + i; r.push(sb.kind[c] === 2 ? 0.5 * (sb.z0[c - 1] + sb.z0[c + 1]) : sb.z0[c]); } rows.push(r); }
P = defaultParams(); P.terrain = { rows, cs: 2 }; P.floorH = 2; P.perturb = 0; P.dams = [{ x: 150, H: 1.2, b: 0 }]; P.hydro.duration = 3600; P.hydro.tPeak = 1200;
run(P, 'csv ');
P = defaultParams(); P.bend = { mode: 'table', table: [[0, 0], [80, 0], [160, 25], [240, -25], [320, 0], [400, 0]] }; P.dams = []; P.hydro.duration = 1800; P.hydro.tPeak = 600;
run(P, 'tab ');
