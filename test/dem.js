// DEM の流域：例の DEM で、谷頭の供給点・格子の雨・斜面崩壊・天然ダムを計算し、水と土砂の収支を確かめる
// node test/dem.js
const M = require('../src/dem.js'), HY = require('../src/hydro.js');
const { Sim, defaultParams } = require('../src/sabo2d-solver.js');
const D = M.synthDEM({}), A = M.analyze(D, { thrArea: 20000 });
const b0 = A.basins.slice().sort((a, b) => b.area - a.area)[0], mask = M.upstream(D, A, b0.outlet);
const dom = M.prepareDomain(D, A, mask, b0.outlet, {}), heads = A.heads.filter(c => mask[c]);
console.log('DEM', D.nx, 'x', D.ny, '谷頭', A.heads.length, '流域', A.basins.length, '選んだ流域の谷頭', heads.length);
function base() {
  const P = defaultParams();
  P.terrain = { rows: dom.rows, cs: dom.cs, basin: dom.basin, chan: dom.chan, Iout: dom.Iout };
  P.dx = 10; P.csvSlope = 0; P.transport = 'debris'; P.dams = []; P.tribs = []; P.snapEvery = 0; P.dem.duration = 300;
  P.grains = [{ d: 0.01, bed: 0.3, feed: 0.3 }, { d: 0.05, bed: 0.4, feed: 0.4 }, { d: 0.2, bed: 0.3, feed: 0.3 }];
  return P;
}
function run(name, P) {
  const s = new Sim(P); while (s.phase !== 'done') s.step();
  const a = s.dx * s.dx, Cst = 1 - P.porosity; let hs = 0, sw = 0;
  for (let c = 0; c < s.N; c++) { hs += s.h[c]; for (let k = 0; k < s.K; k++) sw += s.Cc[c * s.K + k]; }
  hs *= a; sw *= a; const st = s.storage();
  const eb = st * Cst + sw + s.sedOut - s.sedIn;
  const ew = P.gw.on ? (hs - sw + s.groundWater()) - (s.Vw0 + s.W0 + (s.Vwin - s.sedIn) - (s.Vwout - s.sedOut) - s.gwOut) : (hs + st - s.voidAir) - (s.Vw0 + s.Vwin - s.Vwout);
  console.log(name.padEnd(14), 'Vwin', s.Vwin.toFixed(0).padStart(7), 'sedIn', s.sedIn.toFixed(0).padStart(5), '土砂の収支', eb.toExponential(1), '水の収支', ew.toExponential(1), s.slideLog.map(l => l.mode + ' ' + l.vol.toFixed(0)).join(' '));
}
const src = heads.slice(0, 3).map(c => { const p = dom.toXY(c), u = dom.dirXY(c); return { x: p[0], y: p[1], ux: u[0], uy: u[1], I: M.localSlope(D, A, c, 3), hydro: { mode: 'debris', Qbase: 0.5, Qpeak: 10, t0: 10, tRise: 30, tFall: 120, Cd: 0.3, Cs: 0 } }; });
let P = base(); P.sources = src; run('谷頭の供給点', P);
P = base(); P.rain = { on: true, method: 'grid', mode: 'uniform', table: [[0, 60], [5, 0]], loss: { type: 'f', f: 0.8 } }; run('格子の雨', P);
P = base(); const L = HY.lateralInflow(M, D, A, mask, { mode: 'uniform', table: [[0, 50], [10, 0]] }, { loss: { type: 'f1rsa', f1: 0.6, Rsa: 20 }, route: { type: 'storage', K: 8, P: 0.6 }, dt: 60, tEnd: 300 });
P.latIn = { dt: L.dt, pts: L.cells.map(c => { const xy = dom.toXY(c), u = dom.dirXY(c); return [xy[0], xy[1], u[0], u[1], 0.05]; }), Q: L.Q.map(a => Array.from(a)) }; run('雨の横流入', P);
const hp = dom.toXY(heads[1]);
for (const mv of ['fluid', 'block']) { P = base(); P.slides = [{ x: hp[0] + 15, y: hp[1] + 15, R: 35, D: 2, type: 'shallow', trigger: 'time', t: 20, move: mv }]; run('崩壊 ' + mv, P); }
let sc = heads[0]; for (let q = 0; q < 100; q++) { const d = M.downOf(D, A.dir, sc); if (d < 0) break; sc = d; }
const lp = dom.toXY(sc), lu = dom.dirXY(sc);
P = base(); P.ldams = [{ x: lp[0], y: lp[1], ux: lu[0], uy: lu[1], H: 12, Wc: 10, up: 30, dn: 25 }]; P.gw.on = true; P.ldam.k = 3e-2; P.debris.erosion = 'auto';
P.sources = [{ x: lp[0] - 30 * lu[0], y: lp[1] - 30 * lu[1], ux: lu[0], uy: lu[1], I: 0.05, hydro: { mode: 'const', Qpeak: 8, Qbase: 8, Cs: 0, duration: 1e9, tPeak: 1 } }];
run('天然ダム', P);
