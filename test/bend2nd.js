// 2次流の効果の確認：蛇行水路の頂部で、流速最大の位置と外岸・内岸の河床高を比べる
const { Sim, defaultParams } = require('../src/sabo2d-solver.js');
const tag = process.argv[2] || 'none';
const P = defaultParams();
P.L = 360; P.B = 12; P.slope = 1 / 100; P.crossSlope = 0; P.dx = 1.0; P.sideCells = 6; P.sideGrad = 1.5;
P.valley = { chD: 0, chW: 0, chWl: 0, relief: 0.0, wallH: 3, wallRough: 0, fans: [] };
P.bend = { mode: 'sine', amp: 12, wl: 120, x0: 60, x1: 300 };
P.grains = [{ d: 0.01, bed: 0.5, feed: 0.5 }, { d: 0.03, bed: 0.5, feed: 0.5 }];
P.hydro = { mode: 'const', Qbase: 20, Qpeak: 20, tPeak: 0, duration: 1800, Cs: 0.002, table: [] };
P.dams = []; P.morfac = 4;
P.sec = { bed: tag.includes('bed'), mom: tag.includes('mom'), Nstar: 11.5, hrMax: 0.05 };
const s = new Sim(P);
while (s.phase !== 'done') s.step();
// 頂部（yo が極大・極小になる x = 60+30, 60+90, 60+150, 60+210）
let out = [];
for (const xa of [150, 210, 270]) {
  const i = Math.round(xa / s.dx - 0.5);
  let ys = [], vs = [], zs = [];
  for (let j = 0; j < s.ny; j++) { const c = j * s.nx + i; if (s.h[c] > 0.05) { const f = j * (s.nx + 1) + i; const u = 0.5 * (s.FxM[f] + s.FxM[f + 1]) / s.h[c], v = 0.5 * (s.GyM[c] + s.GyM[c + s.nx]) / s.h[c]; ys.push((j + 0.5) * s.dx); vs.push(Math.hypot(u, v)); zs.push(s.z[c] - s.z0[c]); } }
  const yc = (ys[0] + ys[ys.length - 1]) / 2, sgn = Math.sin(2 * Math.PI * (xa - 60) / 120) > 0 ? 1 : -1;   // 外岸 = sgn 側
  let im = 0; vs.forEach((v, k) => { if (v > vs[im]) im = k; });
  const n = zs.length, q = Math.max(1, Math.floor(n / 4));
  const zo = sgn > 0 ? zs.slice(n - q) : zs.slice(0, q), zi = sgn > 0 ? zs.slice(0, q) : zs.slice(n - q);
  const avg = a => a.reduce((x, y) => x + y, 0) / a.length;
  out.push(`x=${xa}: 流速最大の外岸側へのずれ ${(sgn * (ys[im] - yc)).toFixed(1)} m, 外岸−内岸の河床変化 ${(avg(zo) - avg(zi)).toFixed(2)} m`);
}
let mk = 0; for (let c = 0; c < s.N; c++) mk = Math.max(mk, Math.abs(s.kap[c] * s.h[c]));
console.log(tag, '| max h/r', mk.toFixed(3), '|', out.join(' | '), '| 収支', ((s.sedIn - s.sedOut) - 0.6 * s.storage() - s.sedimentInWater()).toExponential(1));
