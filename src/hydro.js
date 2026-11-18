/* 降雨と流出（谷へ入る水）の計算。ブラウザと Node の両方で動く
 *  雨：一様（ハイエトグラフの表）か、格子（時刻ごとの雨量の格子、XRAIN などを変換したもの）
 *  損失（有効雨量）：
 *    'f'     流出係数 f：re = f·r（合理式の考え方）
 *    'f1rsa' 一次流出率 f1・飽和雨量 Rsa・初期損失 R0：累加雨量 ΣR < R0 で f = 0、R0 ≤ ΣR < R0+Rsa で f1、それ以上で 1（貯留関数法の標準の扱い）
 *    'ks'    浸透能（飽和透水係数 ks）と土層の空き容量 Smax = D·(n − θ0)：満ちるまでは r − min(r, ks) が表面流（ホートン型）、
 *            満ちたら r 全部が表面流（飽和余剰）
 *  流出（斜面から谷へ）：
 *    'rational' 合理式の考え方：到達時間 Tc 内の平均有効雨量強度 × 面積（時間面積法、面積は Tc の間に一様に到達）
 *               Tc は土研式（自然流域）T = 1.67×10⁻³ (L/√S)^0.7 [h]（L：斜面の最遠点からの流下距離 m、S：平均勾配）か指定値
 *    'storage'  貯留関数法 s = K q^P、ds/dt = re − q（q, re：mm/h、s：mm）、遅れ Tl
 *  谷筋のセルごとに、そこへ直接流れ込む斜面の面積から流量を出す（谷頭は上流の全面積）
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.SaboHydro = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // 時刻 t [s] の雨量強度 [mm/h]。表 [[時刻 分, mm/h], ...] は、その時刻から次の時刻までの強度（階段状）
  function rainUniform(tab, t) {
    if (!tab || !tab.length) return 0;
    var tm = t / 60, v = 0;
    for (var i = 0; i < tab.length; i++) { if (tm >= tab[i][0] - 1e-9) v = tab[i][1]; else break; }
    if (tab.length > 1 && tm >= tab[tab.length - 1][0] + (tab[tab.length - 1][0] - tab[tab.length - 2][0])) v = 0;   // 最後の区間の後は 0
    return Math.max(0, v);
  }
  // 格子の雨：frames = [{t: 分, r: 配列}]。時刻 t の直前の格子（階段状）
  function frameAt(frames, t) {
    if (!frames || !frames.length) return null;
    var tm = t / 60, f = null;
    for (var i = 0; i < frames.length; i++) { if (tm >= frames[i].t - 1e-9) f = frames[i]; else break; }
    if (f && frames.length > 1) { var dtf = frames[frames.length - 1].t - frames[frames.length - 2].t; if (tm >= frames[frames.length - 1].t + dtf) f = null; }
    return f;
  }

  // 損失モデル（1 セル分の状態を st に持つ）。dt [s]、r [mm/h] → 有効雨量 re [mm/h]
  function lossStep(L, st, r, dt) {
    var dth = dt / 3600;
    if (L.type === 'f1rsa') {
      st.R = (st.R || 0) + r * dth;
      var R0 = L.R0 || 0, f = st.R < R0 ? 0 : (st.R < R0 + (L.Rsa || 0) ? (L.f1 != null ? L.f1 : 0.5) : 1);
      return f * r;
    }
    if (L.type === 'ks') {
      var Smax = Math.max(0, (L.D || 1) * ((L.n || 0.45) - (L.th0 || 0.25)) * 1000);   // mm
      var ks = (L.ks || 1e-5) * 3.6e6;    // m/s → mm/h
      st.S = st.S || 0;
      if (st.S >= Smax) return r;
      var inf = Math.min(r, ks), add = inf * dth;
      if (st.S + add > Smax) { var over = (st.S + add - Smax) / dth; st.S = Smax; return r - inf + over; }
      st.S += add; return r - inf;
    }
    return (L.f != null ? L.f : 0.8) * r;
  }

  // 土研式（自然流域）の到達時間 [s]
  function tcPWRI(L, S) { return 3600 * 1.67e-3 * Math.pow(Math.max(1, L) / Math.sqrt(Math.max(1e-3, S)), 0.7); }

  /* 谷筋への流入の計算
   *  D, A：DEM と解析結果（SaboDEM.analyze）、mask：流域、rain：{mode:'uniform'|'grid', table, frames（DEM のセルに対応する配列）}
   *  opt：{loss:{type,...}, route:{type:'rational'|'storage', tc（分、0 で自動）, K, P, Tl（分）}, dt（s）, tEnd（s）}
   *  返り値：{cells: 谷筋のセルの番号, area: 直接の集水面積 m², Q: [Float64Array（時刻ごと m³/s）], dt, n, Qsum, Rmean}
   */
  function lateralInflow(SD, D, A, mask, rain, opt) {
    var dt = opt.dt || 60, tEnd = opt.tEnd || 3600 * 6, nT = Math.max(2, Math.ceil(tEnd / dt) + 1);
    var N = D.nx * D.ny, a1 = D.cs * D.cs, c, i, k;
    var HS = SD.hillslopes(D, A, mask);
    // 受け手（谷筋のセル）ごとに：面積、最遠の流下距離、平均勾配、雨の平均（格子の雨はセルごとに損失を計算して面積平均）
    var idx = new Int32Array(N).fill(-1), recv = [];
    for (c = 0; c < N; c++) if (mask[c] && A.stream[c]) { idx[c] = recv.length; recv.push({ c: c, n: 0, Lmax: 0, Ssum: 0 }); }
    var members = [];
    for (c = 0; c < N; c++) {
      if (!mask[c]) continue;
      var rc = HS.recv[c]; if (rc < 0 || idx[rc] < 0) continue;
      var R = recv[idx[rc]]; R.n++; if (HS.dist[c] > R.Lmax) R.Lmax = HS.dist[c];
      var dn = SD.downOf(D, A.dir, c);
      var sl = dn >= 0 ? Math.max(0, (A.zf[c] - A.zf[dn]) / (SD.DL[A.dir[c]] * D.cs)) : 0.01; R.Ssum += sl;
      members.push(c);
    }
    var nR = recv.length, Q = [], re = new Float64Array(nR), st = [], sto = new Float64Array(nR), Tc = new Float64Array(nR);
    var L = opt.loss || { type: 'f', f: 0.8 }, RT = opt.route || { type: 'rational' };
    for (i = 0; i < nR; i++) {
      Q.push(new Float64Array(nT));
      var Sm = recv[i].n ? recv[i].Ssum / recv[i].n : 0.1;
      Tc[i] = RT.tc > 0 ? RT.tc * 60 : Math.max(300, Math.min(7200, tcPWRI(recv[i].Lmax + D.cs / 2, Sm)));
    }
    var grid = rain.mode === 'grid' && rain.frames && rain.frames.length;
    var lossSt = grid ? members.map(function () { return {}; }) : [{}];
    var hist = []; for (i = 0; i < nR; i++) hist.push(new Float64Array(nT));   // 有効雨量の履歴（合理式の移動平均用）
    var K = RT.K > 0 ? RT.K : 10, Pp = RT.P > 0 ? RT.P : 0.6, lagN = Math.round((RT.Tl || 0) * 60 / dt);
    var Rsum = 0, Asum = 0;
    for (i = 0; i < nR; i++) Asum += recv[i].n * a1;
    for (k = 0; k < nT; k++) {
      var t = k * dt;
      // 有効雨量 re（受け手ごとの面積平均）
      for (i = 0; i < nR; i++) re[i] = 0;
      if (grid) {
        var fr = frameAt(rain.frames, t);
        for (var m = 0; m < members.length; m++) {
          c = members[m]; var r = fr ? Math.max(0, fr.r[c] || 0) : 0;
          var e = lossStep(L, lossSt[m], r, dt); re[idx[HS.recv[c]]] += e / recv[idx[HS.recv[c]]].n;
          Rsum += r * dt / 3600 / 1000 * a1;
        }
      } else {
        var ru = rainUniform(rain.table, t), eu = lossStep(L, lossSt[0], ru, dt);
        for (i = 0; i < nR; i++) re[i] = eu;
        Rsum += ru * dt / 3600 / 1000 * Asum;
      }
      for (i = 0; i < nR; i++) {
        hist[i][k] = (k > 0 ? hist[i][k - 1] : 0) + re[i];   // 有効雨量の累計（移動平均を速く求める）
        var qmm;
        if (RT.type === 'storage') {
          // ds/dt = re − q、q = (s/K)^(1/P)。2 次の Runge–Kutta
          var dth = dt / 3600, s0 = sto[i];
          var q0 = Math.pow(Math.max(0, s0) / K, 1 / Pp), s1 = s0 + dth * (re[i] - q0);
          var q1 = Math.pow(Math.max(0, s1) / K, 1 / Pp);
          sto[i] = Math.max(0, s0 + 0.5 * dth * (2 * re[i] - q0 - q1));
          qmm = Math.pow(sto[i] / K, 1 / Pp);
        } else {
          // 合理式（時間面積法）：到達時間 Tc の間の有効雨量強度の平均
          var nb = Math.max(1, Math.round(Tc[i] / dt));
          qmm = (hist[i][k] - (k - nb >= 0 ? hist[i][k - nb] : 0)) / nb;
        }
        var kk = k + lagN; if (kk < nT) Q[i][kk] = qmm * recv[i].n * a1 / 3.6e6;
      }
    }
    var Qsum = new Float64Array(nT);
    for (i = 0; i < nR; i++) for (k = 0; k < nT; k++) Qsum[k] += Q[i][k];
    return { cells: recv.map(function (R) { return R.c; }), area: recv.map(function (R) { return R.n * a1; }), Tc: Tc, Q: Q, dt: dt, n: nT, Qsum: Qsum, rainVol: Rsum, area0: Asum };
  }

  return { rainUniform: rainUniform, frameAt: frameAt, lossStep: lossStep, tcPWRI: tcPWRI, lateralInflow: lateralInflow };
});
