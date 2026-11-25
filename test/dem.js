// DEM の流域：例の DEM で、谷頭の供給点・格子の雨・斜面崩壊・天然ダムを計算し、水と土砂の収支を確かめる
// node test/dem.js
const M = require('../src/dem.js'), HY = require('../src/hydro.js');
const { Sim, defaultParams } = require('../src/sabo2d-solver.js');
const D = M.synthDEM({}), A = M.analyze(D, { thrArea: 20000 });
const b0 = A.basins.slice().sort((a, b) => b.area - a.area)[0], mask = M.upstream(D, A, b0.outlet);
const dom = M.prepareDomain(D, A, mask, b0.outlet, { dx: 10 }), heads = A.heads.filter(c => mask[c]);
console.log('DEM', D.nx, 'x', D.ny, '谷頭', A.heads.length, '流域', A.basins.length, '選んだ流域の谷頭', heads.length);
function base() {
  const P = defaultParams();
  P.terrain = { rows: dom.rows, cs: dom.cs, basin: dom.basin, chan: dom.chan, corr: dom.corr, outXY: dom.outXY, Iout: dom.Iout };
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
  return s;
}
const src = heads.slice(0, 3).map(c => { const p = dom.toXY(c), u = dom.dirXY(c); return { x: p[0], y: p[1], ux: u[0], uy: u[1], I: M.localSlope(D, A, c, 3), hydro: { mode: 'debris', Qbase: 0.5, Qpeak: 10, t0: 10, tRise: 30, tFall: 120, Cd: 0.3, Cs: 0 } }; });
let P = base(); P.sources = src; run('谷頭の供給点', P);
P = base(); P.rain = { on: true, method: 'grid', mode: 'uniform', table: [[0, 60], [5, 0]], loss: { type: 'f', f: 0.8 } }; run('格子の雨', P);
P = base(); const L = HY.lateralInflow(M, D, A, mask, { mode: 'uniform', table: [[0, 50], [10, 0]] }, { loss: { type: 'f1rsa', f1: 0.6, Rsa: 20 }, route: { type: 'storage', K: 8, P: 0.6 }, dt: 60, tEnd: 300 });
P.latIn = { dt: L.dt, pts: L.cells.map(c => { const xy = dom.toXY(c), u = dom.dirXY(c); return [xy[0], xy[1], u[0], u[1], 0.05]; }), Q: L.Q.map(a => Array.from(a)) }; run('雨の横流入', P);
const hp = dom.toXY(heads[1]);
{ P = base(); const v = [[-10, -25], [40, -20], [55, 20], [20, 45], [-15, 25]].map(q => [hp[0] + 15 + q[0], hp[1] + 15 + q[1]]);
  P.slides = [{ poly: v, x: 0, y: 0, D: 6, type: 'deep', trigger: 'time', t: 20, move: 'block' }]; run('崩壊 多角形', P); }
for (const mv of ['fluid', 'block']) { P = base(); P.slides = [{ x: hp[0] + 15, y: hp[1] + 15, R: 35, D: 2, type: 'shallow', trigger: 'time', t: 20, move: mv }]; run('崩壊 ' + mv, P); }
let sc = heads[0]; for (let q = 0; q < 100; q++) { const d = M.downOf(D, A.dir, sc); if (d < 0) break; sc = d; }
const lp = dom.toXY(sc), lu = dom.dirXY(sc);
P = base(); P.ldams = [{ x: lp[0], y: lp[1], ux: lu[0], uy: lu[1], H: 12, Wc: 10, up: 30, dn: 25 }]; P.gw.on = true; P.ldam.k = 3e-2; P.debris.erosion = 'auto';
P.sources = [{ x: lp[0] - 30 * lu[0], y: lp[1] - 30 * lu[1], ux: lu[0], uy: lu[1], I: 0.05, hydro: { mode: 'const', Qpeak: 8, Qbase: 8, Cs: 0, duration: 1e9, tPeak: 1 } }];
run('天然ダム', P);
// 砂防堰堤（DEM）：本川の途中に不透過とスリットの 2 基。上流の範囲と収支
{
  let c1 = heads[0]; for (let q = 0; q < 60; q++) { const d = M.downOf(D, A.dir, c1); if (d < 0) break; c1 = d; }
  let c2 = c1; for (let q = 0; q < 40; q++) { const d = M.downOf(D, A.dir, c2); if (d < 0) break; c2 = d; }
  const mk = (c, H, ty) => { const p = dom.toXY(c), u = dom.dirXY(c); const d = { x: p[0], y: p[1], ang: Math.atan2(u[1], u[0]) * 180 / Math.PI, H, b: 15, so: 0, len: 200 }; if (ty) Object.assign(d, { type: ty, l: 2, n: 1 }); return d; };
  P = base(); P.sources = src; P.dem.duration = 400; P.dams = [mk(c1, 6), mk(c2, 8, 'slit')];
  const s0 = new Sim(P); const cnt = new Array(s0.dams.length + 1).fill(0); for (let c = 0; c < s0.N; c++) cnt[s0.region[c]]++;
  let pond = 0; for (let c = 0; c < s0.N; c++) pond += s0.h[c] * s0.dx * s0.dx;
  console.log('堰堤', s0.dams.length, '基', s0.dams.map(d => 'crest ' + d.crest.toFixed(1) + ' cells ' + d.cells.length + (d.op ? ' 開口 ' + d.op.cells.length : '')).join(' / '), '範囲のセル数', cnt.join(','), '初期湛水', pond.toFixed(0), 'm³', s0.warn.join(' '));
  run('砂防堰堤', P);
}
// 出口：自由流出（既定）と水路を延ばす（以前）。どちらも収支が合い、自由流出の方が出口の先に水を残さない
for (const mode of ['sink', 'channel']) {
  P = base(); P.sources = src.map(q => Object.assign({}, q, { hydro: { mode: 'const', Qbase: 4, Qpeak: 4, Cs: 0, duration: 1e12, tPeak: 1 } })); P.dem.duration = 900; P.dem.outlet = mode;
  const s = run('出口 ' + mode, P);
  if (s) { let vc = 0; (s.sinkCells || []).forEach(c => { vc += s.h[c]; }); console.log('   出口の先の水路のセル', (s.sinkCells || []).length, '流出 Vwout', s.Vwout.toFixed(0), 'm³', '穴埋め', s.nPit || 0, '窪地埋め', s.nFill || 0); }
}
