const { Sim, defaultParams } = require('../src/sabo2d-solver.js');
function base() { const P = defaultParams(); P.dams = []; return P; }

// A. 静水（lake at rest）
(function () {
  const P = base(); P.L = 120; P.B = 30; P.perturb = 0.3; P.hydro.mode = 'const'; P.hydro.Qpeak = 0; P.hydro.Cs = 0;
  P.dams = [{ x: 60, H: 2, b: 10 }];
  for (const order of [1, 2]) {
    P.order = order;
    const s = new Sim(P);
    const eta0 = 3.2;
    for (let c = 0; c < s.N; c++) { s.h[c] = Math.max(0, eta0 - s.z[c]); s.qx[c] = 0; s.qy[c] = 0; }
    s.phase = 'spinup'; s.Uref = 1e-9;  // 河床変動なし
    s.phase = 'x';
    for (let n = 0; n < 300; n++) { s.phase = 'done'; s.step(); }
    let um = 0, de = 0;
    for (let c = 0; c < s.N; c++) { um = Math.max(um, Math.abs(s.qx[c]), Math.abs(s.qy[c])); if (s.h[c] > 0.01) de = Math.max(de, Math.abs(s.z[c] + s.h[c] - eta0)); }
    console.log(`A 静水 order=${order}: max|q|=${um.toExponential(2)}  max|η-η0|=${de.toExponential(2)}`);
  }
})();

// B. 等流（Manning の等流水深との比較）
(function () {
  for (const I of [0.02, 0.05, 0.1]) for (const order of [1, 2]) {
    const P = base(); P.L = 200; P.B = 20; P.crossSlope = 0; P.perturb = 0; P.slope = I; P.order = order;
    P.hydro.mode = 'const'; P.hydro.Qpeak = 20; P.hydro.Cs = 0; P.sideGrad = 2;
    const s = new Sim(P);
    for (let n = 0; n < 6000; n++) { s.phase = 'done'; s.step(); }
    const q = 20 / s.Bf, hn = Math.pow(q * P.manning / Math.sqrt(I), 0.6);
    const j = (s.ny >> 1), i = (s.nx >> 1), c = j * s.nx + i;
    console.log(`B 等流 I=${I} order=${order}: h=${s.h[c].toFixed(4)} hn=${hn.toFixed(4)} 誤差=${((s.h[c] / hn - 1) * 100).toFixed(2)}%  u=${(s.qx[c] / s.h[c]).toFixed(3)} un=${(q / hn).toFixed(3)} Qout=${s.Qout.toFixed(3)}`);
  }
})();
