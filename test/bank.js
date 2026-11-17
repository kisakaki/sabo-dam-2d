// 側岸侵食・非平衡流砂・巨礫の検証：体積の収支と、谷底幅の広がり
const { Sim, defaultParams } = require('../src/sabo2d-solver.js');
const mode = process.argv[2] || 'eq', bankH = +(process.argv[3] || 0), dry = +(process.argv[4] || 0), bL = +(process.argv[5] || 0);
const P = defaultParams();
P.L = 300; P.B = 30; P.slope = 1 / 30; P.crossSlope = 0; P.dx = 1.5; P.sideCells = 6; P.sideGrad = 1.0;
P.valley = { chD: 0.8, chW: 10, chWl: 0, relief: 0.15, wallH: 4, wallRough: 0, fans: [] };
P.bank = { erodH: bankH, reposeWet: 30, reposeDry: 50, dryEro: dry };
P.transport = mode; P.boulders = { L: bL, mode: 'parallel', size: 1 };
P.hydro = { mode: 'const', Qbase: 20, Qpeak: 20, tPeak: 0, duration: 2400, Cs: 0.002, table: [] };
P.dams = [];
const s = new Sim(P);
const wet0 = (() => { let n = 0; for (let c = 0; c < s.N; c++) if (s.kind[c] === 0 && !s.bankC[c]) n++; return n; })();
const t0 = Date.now(); let n = 0;
while (s.phase !== 'done') { s.step(); n++; }
let bad = 0; for (let c = 0; c < s.N; c++) if (!isFinite(s.h[c]) || !isFinite(s.z[c])) bad++;
// 谷底（初期の谷底より 0.3 m 以上低い所も含む）の幅：初期の岸より外で 0.3 m 以上削れたセル数
let eroBank = 0, vBank = 0; for (let c = 0; c < s.N; c++) if (s.bankC[c] && s.z[c] - s.z0[c] < -0.3) { eroBank++; vBank += (s.z0[c] - s.z[c]) * s.dx * s.dx; }
const susp = s.sedimentInWater ? s.sedimentInWater() : 0;
console.log(`${mode} bankH=${bankH} dry=${dry} boulders=${bL} (${s.nBoulder || 0} cells) grid ${s.nx}x${s.ny} steps ${n} ${((Date.now() - t0) / 1000).toFixed(1)}s bad ${bad}`);
console.log(`  sedBal ${((s.sedIn - s.sedOut) - (1 - P.porosity) * s.storage() - susp).toExponential(2)} m3  watBal ${(s.Vw0 + s.Vwin - s.Vwout - s.waterVolume()).toExponential(1)}  in ${s.sedIn.toFixed(0)} out ${s.sedOut.toFixed(0)} susp ${susp.toFixed(2)}  bank cells eroded ${eroBank} (${vBank.toFixed(0)} m3)  storage ${s.storage().toFixed(0)}`);
