const { Sim, defaultParams } = require('../src/sabo2d-solver.js');
for (const [m, hd] of [['atm', 'atm'], ['iwagaki', 'egiazaroff'], ['const', 'egiazaroff']]) {
  const P = defaultParams(); P.tcMethod = m; P.hiding = hd; P.hydro.duration = 3600; P.hydro.tPeak = 1200;
  P.dams = [60, 105, 150, 195, 240, 285, 330].map(x => ({ x, H: 1.2, b: 0 }));
  const s = new Sim(P); let n = 0; while (s.phase !== 'done' && n < 2e6) { s.step(); n++; }
  const tc = []; for (let c = 0; c < s.N; c++) if (s.tcs[c] === s.tcs[c] && s.h[c] > 0.02) tc.push(s.tcs[c]); tc.sort((a, b) => a - b);
  let bad = 0; for (let c = 0; c < s.N; c++) if (!isFinite(s.z[c])) bad++;
  console.log(`${m.padEnd(8)}/${hd.padEnd(10)} τ*cm 10/50/90% = ${[0.1, 0.5, 0.9].map(q => tc[Math.floor(q * (tc.length - 1))].toFixed(4)).join(' / ')}  流出土砂 ${(s.sedOut / 0.6).toFixed(0)} m3  区間内堆積 ${s.storage().toFixed(0)} m3  収支 ${((s.sedIn - s.sedOut) / 0.6 - s.storage()).toExponential(1)} bad ${bad}`);
}
