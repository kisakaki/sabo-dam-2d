// 谷底いっぱいに流れる狭い水路で、側岸（段丘）の侵食による拡幅を比べる
const { Sim, defaultParams } = require('../src/sabo2d-solver.js');
const [mode, bankH, dry, wet, drya] = [process.argv[2] || 'eq', +(process.argv[3] || 0), +(process.argv[4] || 0), +(process.argv[5] || 30), +(process.argv[6] || 50)];
const P = defaultParams();
P.L = 240; P.B = 10; P.slope = 1 / 40; P.crossSlope = 0.02; P.dx = 1.0; P.sideCells = 8; P.sideGrad = 1.2;
P.valley = { chD: 0, chW: 0, chWl: 0, relief: 0.1, wallH: 3, wallRough: 0, fans: [] };
P.bank = { erodH: bankH, reposeWet: wet, reposeDry: drya, dryEro: dry };
P.transport = mode;
P.hydro = { mode: 'const', Qbase: 25, Qpeak: 25, tPeak: 0, duration: 1800, Cs: 0.0, table: [] };
P.dams = []; P.morfac = 4;
const s = new Sim(P);
while (s.phase !== 'done') s.step();
// 水面幅（x = 60〜200 m の平均）と、初期の谷底より外で削れた体積
let wsum = 0, nw = 0;
for (let i = Math.round(60 / P.dx); i < Math.round(200 / P.dx); i++) { let w = 0; for (let j = 0; j < s.ny; j++) if (s.h[j * s.nx + i] > 0.05) w++; wsum += w * s.dx; nw++; }
let vB = 0; for (let c = 0; c < s.N; c++) if (s.bankC[c] && s.z[c] < s.z0[c]) vB += (s.z0[c] - s.z[c]) * s.dx * s.dx;
const susp = s.sedimentInWater();
console.log(`${mode} erodH=${bankH} dryEro=${dry} wet=${wet} dry=${drya}: 水面幅 ${(wsum / nw).toFixed(1)} m  側岸の侵食 ${vB.toFixed(0)} m3  流出 ${s.sedOut.toFixed(0)} m3  収支 ${((s.sedIn - s.sedOut) - 0.6 * s.storage() - susp).toExponential(1)}`);
