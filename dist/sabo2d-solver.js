/*
 * sabo2d-solver.js  —  連続砂防堰堤 平面2次元河床変動計算ソルバー
 * ---------------------------------------------------------------------------
 * 支配方程式
 *   流れ    : 平面2次元浅水流方程式（保存形）
 *               ∂h/∂t + ∂(hu)/∂x + ∂(hv)/∂y = 0
 *               ∂(hu)/∂t + ∂(hu²+gh²/2)/∂x + ∂(huv)/∂y = -gh ∂z/∂x - g n² u√(u²+v²)/h^(1/3)
 *               ∂(hv)/∂t + ∂(huv)/∂x + ∂(hv²+gh²/2)/∂y = -gh ∂z/∂y - g n² v√(u²+v²)/h^(1/3)
 *   河床    : 粒径階 k ごとの流砂の連続式（Exner）
 *               (1-λ) ∂z_k/∂t + ∂q_bx,k/∂x + ∂q_by,k/∂y = 0 ,  z = Σ z_k
 * 離散化
 *   流れ    : セル中心有限体積法。HLL近似リーマン解法 + 静水圧再構築
 *             (Audusse et al. 2004 / Chen & Noelle 2017)。水深 h と水位 η=z+h を
 *             minmod 制限付きで区分線形に再構築（等流と静水の両方を厳密に保つ）。
 *             時間積分は陽的Euler、摩擦項は半陰的。
 *   流砂    : 芦田・道上式（1972）。限界掃流力は芦田・高橋・水山（1978）式(24)
 *             τ*c = 0.04×10^(1.72 i) または岩垣式（1956）。粒径別の限界掃流力は
 *             芦田・高橋・水山の結果（全粒径で平均粒径の値）または Egiazaroff 式の芦田・道上修正。
 *             流砂の向きは流速方向＋河床勾配による重力効果（長谷川 1981 型の係数を
 *             等方的な勾配拡散として適用）。
 *   粒度    : 交換層（平野 1971）＋堆積層（多層）で粒度の鉛直履歴を保存。
 *   側岸    : 安息角を超えた斜面の崩落。
 *   堰堤    : 1セル幅の固定床（非侵食）として天端高・水通し・袖部を与える。
 * 適用範囲
 *   掃流砂が卓越する領域（概ね河床勾配 1/10 以下）。土石流・掃流状集合流動、
 *   浮遊砂・ウォッシュロード、二次流は扱わない。
 * ---------------------------------------------------------------------------
 */
(function (root) {
  'use strict';
  var G = 9.81;
  var LN19 = Math.log(19);

  function defaultParams() {
    return {
      // ---- 河道地形
      L: 400,            // 計算区間長 [m]
      B: 50,             // 谷底（河床）幅 [m]
      widthTable: [],    // 川幅の変化点 [[x(m), 幅(m)], ...]（線形補間。空なら B で一様）
      dx: 2,             // 格子幅 [m]（dx = dy）
      slope: 0.05,       // 縦断勾配 I [-]
      crossSlope: 0.02,  // 横断勾配（谷中央へ向かう勾配） [-]
      perturb: 0.05,     // 初期河床の擾乱振幅 [m]
      seed: 1,           // 乱数の種
      erodDepth: 3.0,    // 侵食可能深（移動床厚） [m]
      sideCells: 5,      // 谷壁（非侵食斜面）のセル数（片側）
      sideGrad: 0.8,     // 谷壁の勾配 [-]
      // ---- 渓流の作り込み（数値で指定する地形のみ。0 のとき使わない）
      valley: {
        chD: 0,          // 澪筋（低水路）の深さ [m]
        chW: 0,          // 澪筋の幅 [m]（0 のとき谷底幅の 3 割）
        chWl: 0,         // 澪筋の蛇行波長 [m]（0 のとき谷底幅の 4 倍）
        relief: 0,       // 谷底の起伏（砂州・淵）の振幅 [m]
        wallH: 0,        // 谷壁の高さ [m]（0 のとき堰堤の高さから決める）
        wallRough: 0,    // 岸の出入り（左右の岸線の不規則な凹凸）の振幅 [m]
        fans: []         // 支渓の扇状地 [[x(m), 岸(+1 左岸, -1 右岸), 半径(m), 高さ(m)], ...]
      },
      // ---- 流れ
      manning: 0.050,    // Manning の粗度係数 n [m^-1/3 s]
      cfl: 0.4,
      order: 2,          // 空間精度（1: 区分一定, 2: h と η を区分線形再構築）
      // ---- 土砂
      grains: [          // d: 粒径 [m], bed: 河床材料の割合, feed: 供給土砂の割合
        { d: 0.005, bed: 0.20, feed: 0.30 },
        { d: 0.020, bed: 0.30, feed: 0.35 },
        { d: 0.060, bed: 0.30, feed: 0.25 },
        { d: 0.150, bed: 0.20, feed: 0.10 }
      ],
      porosity: 0.4,     // 空隙率 λ
      sg: 2.65,          // 砂礫の比重 σ/ρ
      tcMethod: 'atm',   // 平均粒径の限界掃流力：'atm' 芦田・高橋・水山式 | 'iwagaki' 岩垣式 | 'const' 一定値 tauCm
      hiding: 'atm',     // 粒径別の限界掃流力：'atm' 全粒径で平均粒径の値（芦田・高橋・水山）| 'egiazaroff' Egiazaroff 式の芦田・道上修正
      tauCm: 0.05,       // tcMethod = 'const' のときの τ*cm
      reposeDeg: 35,     // 安息角 [deg]（流砂の勾配項の上限。崩落は下の bank の水中・陸上安息角で判定）
      // ---- 側岸侵食
      bank: {
        erodH: 0,        // 谷壁のうち、谷底から erodH [m] の高さまでを侵食できる側岸（段丘・崖錐）にする。0 で谷壁は固定床
        reposeWet: 30,   // 水中の安息角 [deg]（高橋・里深・奥村 1996 は 30°）
        reposeDry: 50,   // 陸上の安息角 [deg]（乾いた側岸が立っていられる角度。根や見かけの粘着力で急でも立つ）
        dryEro: 0,       // 水際の濡れたセルの侵食のうち、隣の乾いた側岸から削る割合（Delft3D の THETSD と同じ考え方）。0〜1
        collapse: 'repose' // 'repose' 安息角を超えた分を少しずつ崩す | 'hasegawa' 長谷川の崩落モデル（長田ら 1996 式17）
      },
      // ---- 2次流（湾曲部のらせん流）。長田・細田・村本・Rahman（1996）
      sec: {
        bed: false,      // 河床付近の流れを曲率の内側へ向け、掃流砂の向きを変える（式11・14：u_nb/u_sb = N* h/r）
        mom: false,      // 2次流による運動量の分散項（式19〜23、Kalkwijk・de Vriend）
        Nstar: 11.5,     // 式14の係数 N*
        hrMax: 0.05      // h/r の上限（弱い曲率を前提とする式のため。これを超える所は上限値で頭打ち）
      },
      // ---- 流砂の計算法
      transport: 'eq',   // 'eq' 平衡流砂量（芦田・道上式）| 'noneq' 侵食・堆積速度式（高橋・里深 1999、高橋ら 1998）
      noneq: {
        de: 0.6,         // 侵食係数 δe
        dd: 0.6,         // 堆積係数 δd
        hEdge: 0,        // 水際とみなす水深 [m]（0 のとき平均粒径の 2.5 倍）。これより浅い所で係数を割り増す
        fe: 2,           // 水際での侵食係数の倍率（高橋・里深 1999：0.6→1.2）
        fd: 10           // 水際での堆積係数の倍率（同：0.6→6.0）
      },
      // ---- 巨礫（動かない礫）を水みちの両岸に置く（高橋・里深・奥村 1996）
      boulders: { L: 0, mode: 'parallel', size: 1 },   // L：配置間隔 [m]（0 で置かない）、mode：'parallel' つがい｜'zigzag' 千鳥、size：セル数
      activeLayer: 0,    // 交換層厚 [m]（0 のとき自動：最大粒径）
      layerTh: 0.25,     // 堆積層 1 層の厚さ [m]
      morfac: 4,         // 河床変動の加速係数（1 = 加速なし）
      // ---- ハイドログラフ
      hydro: {
        mode: 'tri',     // 'const' | 'tri' | 'table'
        Qbase: 10,       // 基底流量 [m3/s]
        Qpeak: 60,       // ピーク流量 [m3/s]
        tPeak: 3600,     // ピーク時刻 [s]
        duration: 10800, // 継続時間 [s]
        Cs: 0.005,       // 供給土砂の流砂量濃度 Qs/Q（実質体積） [-]
        table: []        // [[t(s), Q(m3/s), Qs(m3/s, 省略可)], ...]
      },
      // ---- 堰堤
      dams: [],          // [{x: 位置[m], H: 堤高[m], b: 水通し幅[m]}]
      wingHeight: 1.5,   // 袖部の高さ（水通し天端から） [m]
      prefill: 0         // >0 のとき堰堤の上流をはじめから満砂にする。値は堆砂勾配と元の勾配の比
    };
  }

  // 再現性のある乱数（mulberry32）
  function rng(seed) {
    var a = (seed | 0) + 0x6D2B79F5;
    return function () {
      a |= 0; a = (a + 0x6D2B79F5) | 0;
      var t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  // 滑らかな乱数場（値ノイズ）。-1〜1
  function hashI(a, b, s) { var h = (Math.imul(a | 0, 374761393) + Math.imul(b | 0, 668265263) + Math.imul(s | 0, 2147483647)) | 0; h = Math.imul(h ^ (h >>> 13), 1274126177); return (((h ^ (h >>> 16)) >>> 0) / 4294967296) * 2 - 1; }
  function vn2(x, y, s) {
    var i = Math.floor(x), j = Math.floor(y), fx = x - i, fy = y - j;
    fx = fx * fx * (3 - 2 * fx); fy = fy * fy * (3 - 2 * fy);
    var a = hashI(i, j, s), b = hashI(i + 1, j, s), c = hashI(i, j + 1, s), d = hashI(i + 1, j + 1, s);
    return (a * (1 - fx) + b * fx) * (1 - fy) + (c * (1 - fx) + d * fx) * fy;
  }
  function vn1(x, s) { return 0.7 * vn2(x, 0.5, s) + 0.3 * vn2(x * 2.7 + 11.3, 3.5, s + 17); }

  // 岩垣式（1956）：限界摩擦速度の2乗 u*c² [m²/s²]。原式は cm・s 単位、σ/ρ = 2.65 の値なので (σ/ρ−1)/1.65 で補正する
  function iwagakiUc2(dm, sg1) {
    var d = dm * 100, u;
    if (d >= 0.303) u = 80.9 * d;
    else if (d >= 0.118) u = 134.6 * Math.pow(d, 31 / 22);
    else if (d >= 0.0565) u = 55.0 * d;
    else if (d >= 0.0065) u = 8.41 * Math.pow(d, 11 / 32);
    else u = 226 * d;
    return u * 1e-4 * sg1 / 1.65;
  }

  function Sim(P) {
    this.P = P;
    this.init();
  }

  // ------------------------------------------------------------------ 初期化
  // 中心線の横ずれ yo(x)（流下方向に向かって左が正） [m]
  function centerline(P, nx, dx) {
    var yo = new Float64Array(nx), Bd = P.bend || {}, i, x;
    if (Bd.mode === 'sine' && Bd.amp > 0 && Bd.wl > 0) {
      var x0 = Bd.x0 || 0, x1 = Bd.x1 > x0 ? Bd.x1 : 1e30, tp = Bd.wl / 2;
      for (i = 0; i < nx; i++) {
        x = (i + 0.5) * dx; if (x <= x0 || x >= x1) continue;
        var a = Math.min(1, (x - x0) / tp, (x1 - x) / tp); a = a * a * (3 - 2 * a);
        yo[i] = Bd.amp * a * Math.sin(2 * Math.PI * (x - x0) / Bd.wl);
      }
    } else if (Bd.mode === 'table' && Bd.table && Bd.table.length) {
      var T = Bd.table.filter(function (r) { return isFinite(r[0]) && isFinite(r[1]); }).sort(function (a, b) { return a[0] - b[0]; });
      var n = T.length;
      for (i = 0; i < nx; i++) {
        x = (i + 0.5) * dx;
        if (n === 1 || x <= T[0][0]) { yo[i] = T[0][1]; continue; }
        if (x >= T[n - 1][0]) { yo[i] = T[n - 1][1]; continue; }
        var k = 1; while (T[k][0] < x) k++;
        // Catmull-Rom（不等間隔）で滑らかにつなぐ
        var p0 = T[Math.max(0, k - 2)], p1 = T[k - 1], p2 = T[k], p3 = T[Math.min(n - 1, k + 1)];
        var hseg = p2[0] - p1[0], t = (x - p1[0]) / hseg;
        var m1 = k - 1 > 0 ? (p2[1] - p0[1]) / (p2[0] - p0[0]) * hseg : (p2[1] - p1[1]);
        var m2 = k < n - 1 ? (p3[1] - p1[1]) / (p3[0] - p1[0]) * hseg : (p2[1] - p1[1]);
        var t2 = t * t, t3 = t2 * t;
        yo[i] = (2 * t3 - 3 * t2 + 1) * p1[1] + (t3 - 2 * t2 + t) * m1 + (-2 * t3 + 3 * t2) * p2[1] + (t3 - t2) * m2;
      }
    }
    return yo;
  }

  // CSV 地形を (x, y) で双一次補間（rows[0] が左岸側 = y の大きい側）
  function sampleTerrain(T, x, y, Ly) {
    var nr = T.rows.length, nc = T.rows[0].length, cs = T.cs;
    var fx = x / cs - 0.5, fy = (Ly - y) / cs - 0.5;
    if (fx < 0) fx = 0; if (fx > nc - 1) fx = nc - 1; if (fy < 0) fy = 0; if (fy > nr - 1) fy = nr - 1;
    var i0 = Math.min(nc - 2, Math.floor(fx)), j0 = Math.min(nr - 2, Math.floor(fy));
    if (nc < 2) i0 = 0; if (nr < 2) j0 = 0;
    var tx = fx - i0, ty = fy - j0, r0 = T.rows[j0], r1 = T.rows[Math.min(nr - 1, j0 + 1)];
    var a = r0[i0] * (1 - tx) + r0[Math.min(nc - 1, i0 + 1)] * tx, b = r1[i0] * (1 - tx) + r1[Math.min(nc - 1, i0 + 1)] * tx;
    return a * (1 - ty) + b * ty;
  }

  Sim.prototype.init = function () {
    var P = this.P;
    var dx = this.dx = P.dx;
    var T = (P.terrain && P.terrain.rows && P.terrain.rows.length >= 3 && P.terrain.rows[0].length >= 3) ? P.terrain : null;
    this.fromCSV = !!T;
    var nx, ny, ns, i, j, c, k, q, Bx, yc, cosT;
    var maxH = 0; (P.dams || []).forEach(function (D) { if (D.H > maxH) maxH = D.H; });
    this.warn = [];
    var wallMax = Math.max((P.valley && P.valley.wallH > 0) ? P.valley.wallH : 6, maxH + P.wingHeight + 3);
    if (T) {
      nx = Math.max(20, Math.round(T.rows[0].length * T.cs / dx));
      ny = Math.max(6, Math.round(T.rows.length * T.cs / dx));
      ns = 0;
      Bx = new Float64Array(nx); yc = new Float64Array(nx); cosT = new Float64Array(nx).fill(1);
    } else {
      nx = Math.max(20, Math.round(P.L / dx));
      // 川幅（法線方向）の縦断変化：変化点の間は線形補間、セル数に丸める。
      // 平面形を表で与え、その表に川幅の列があるときは、そちらを使う（湾曲と狭窄を一緒に与える）
      var bdT = (P.bend && P.bend.mode === 'table' && P.bend.table) ? P.bend.table : [];
      var wt = bdT.filter(function (r) { return isFinite(r[0]) && r[2] > 0; }).map(function (r) { return [r[0], r[2]]; });
      if (!wt.length) wt = (P.widthTable || []).filter(function (r) { return isFinite(r[0]) && r[1] > 0; });
      wt.sort(function (a, b) { return a[0] - b[0]; });
      Bx = new Float64Array(nx);
      for (i = 0; i < nx; i++) {
        var bw = wt.length ? (wt.length > 1 ? interp(wt, (i + 0.5) * dx, 1) : wt[0][1]) : P.B;
        Bx[i] = Math.max(2, Math.round(bw / dx)) * dx;
      }
      // 岸の出入り：左右の岸線を別々に、滑らかな乱数でずらす
      var VL = P.valley || {}, rough = Math.max(0, VL.wallRough || 0);
      var hL = new Float64Array(nx), hR = new Float64Array(nx);
      for (i = 0; i < nx; i++) {
        var xr = (i + 0.5) * dx / 28;
        hL[i] = Math.max(dx, Bx[i] / 2 + (rough > 0 ? rough * vn1(xr, P.seed * 7 + 1) : 0));
        hR[i] = Math.max(dx, Bx[i] / 2 + (rough > 0 ? rough * vn1(xr, P.seed * 7 + 2) : 0));
        if (rough > 0) Bx[i] = hL[i] + hR[i];
      }
      this.hL = hL; this.hR = hR;
      // 平面形（湾曲）：中心線の横ずれと、その傾きによる y 方向幅の補正
      var yo = centerline(P, nx, dx), ymin = 1e30, ymax = -1e30, smax = 0;
      cosT = new Float64Array(nx);
      for (i = 0; i < nx; i++) {
        var sl = (yo[Math.min(nx - 1, i + 1)] - yo[Math.max(0, i - 1)]) / ((Math.min(nx - 1, i + 1) - Math.max(0, i - 1)) * dx);
        if (Math.abs(sl) > smax) smax = Math.abs(sl);
        cosT[i] = 1 / Math.sqrt(1 + sl * sl);
        if (yo[i] - hR[i] / cosT[i] < ymin) ymin = yo[i] - hR[i] / cosT[i]; if (yo[i] + hL[i] / cosT[i] > ymax) ymax = yo[i] + hL[i] / cosT[i];
      }
      if (smax > 1.2) this.warn.push('中心線の傾きが ' + (Math.atan(smax) * 180 / Math.PI).toFixed(0) + '° に達します。50° を超える曲がりは格子の向きと合わず、精度が落ちます。');
      ns = Math.max(1, P.sideCells | 0);
      if (VL.wallH > 0) ns = Math.max(ns, Math.ceil(wallMax / Math.max(0.2, P.sideGrad) / dx) + 1);
      ny = Math.round((ymax - ymin) / dx) + 2 * ns;
      yc = new Float64Array(nx);
      for (i = 0; i < nx; i++) yc[i] = ns * dx + (yo[i] - ymin);
    }
    this.nx = nx; this.ny = ny; this.ns = ns; this.nyb = ny - 2 * ns;
    var N = this.N = nx * ny;
    this.Lx = nx * dx; this.Ly = ny * dx;
    this.Bx = Bx; this.yc = yc; this.cosT = cosT;

    // 粒径階
    var K = this.K = P.grains.length;
    var d = this.d = new Float64Array(K), f0 = this.f0 = new Float64Array(K), ff = this.ff = new Float64Array(K);
    var sb = 0, sf = 0;
    for (k = 0; k < K; k++) { d[k] = P.grains[k].d; sb += Math.max(0, P.grains[k].bed); sf += Math.max(0, P.grains[k].feed); }
    if (sb <= 0) sb = 1; if (sf <= 0) sf = 1;
    var dm0 = 0, dmax = 0;
    for (k = 0; k < K; k++) {
      f0[k] = Math.max(0, P.grains[k].bed) / sb; ff[k] = Math.max(0, P.grains[k].feed) / sf;
      dm0 += f0[k] * d[k]; if (d[k] > dmax) dmax = d[k];
    }
    this.dm0 = dm0; this.dmax = dmax;
    this.La = P.activeLayer > 0 ? P.activeLayer : Math.max(dmax, 0.05);
    this.Ls = Math.max(0.05, P.layerTh);
    this.tanPhi = Math.tan(P.reposeDeg * Math.PI / 180);
    this.sg1 = P.sg - 1;
    this.n2 = P.manning * P.manning;

    // 流れの配列
    this.h = new Float64Array(N); this.qx = new Float64Array(N); this.qy = new Float64Array(N);
    this.u = new Float64Array(N); this.v = new Float64Array(N); this.eta = new Float64Array(N);
    this.z = new Float64Array(N); this.zfix = new Float64Array(N); this.z0 = new Float64Array(N);
    this.dhx = new Float64Array(N); this.dex = new Float64Array(N);
    this.dhy = new Float64Array(N); this.dey = new Float64Array(N);
    this.sc = new Float64Array(N);
    this.kind = new Uint8Array(N);   // 0: 河床, 1: 谷壁, 2: 堰堤
    this.fl = new Uint8Array(N);     // 1: 谷底の範囲（堰堤設置前の区分）
    this.noSl = new Uint8Array(N);
    var NFx = (nx + 1) * ny, NFy = nx * (ny + 1);
    this.FxM = new Float64Array(NFx); this.FxU = new Float64Array(NFx); this.FxV = new Float64Array(NFx);
    this.FxCL = new Float64Array(NFx); this.FxCR = new Float64Array(NFx);
    this.GyM = new Float64Array(NFy); this.GyU = new Float64Array(NFy); this.GyV = new Float64Array(NFy);
    this.GyCL = new Float64Array(NFy); this.GyCR = new Float64Array(NFy);
    this.qin = new Float64Array(ny); this.uin = new Float64Array(ny);

    // ---- 初期河床の擾乱
    var z = this.z, zfix = this.zfix, kind = this.kind, fl = this.fl;
    var rnd = rng(P.seed);
    var noise = new Float64Array(N), tmp = new Float64Array(N);
    for (c = 0; c < N; c++) tmp[c] = rnd() * 2 - 1;
    for (j = 0; j < ny; j++) for (i = 0; i < nx; i++) {   // 3x3 平滑化で格子スケールの凹凸を除く
      var s = 0, cnt = 0;
      for (var jj = -1; jj <= 1; jj++) for (var ii = -1; ii <= 1; ii++) {
        var i2 = i + ii, j2 = j + jj;
        if (i2 >= 0 && i2 < nx && j2 >= 0 && j2 < ny) { s += tmp[j2 * nx + i2]; cnt++; }
      }
      noise[j * nx + i] = s / cnt * 1.8;
    }

    // ---- 地形
    var colmin = this.colmin = new Float64Array(nx);
    if (T) {
      var floorH = P.floorH > 0 ? P.floorH : 3;
      for (j = 0; j < ny; j++) for (i = 0; i < nx; i++) z[j * nx + i] = sampleTerrain(T, (i + 0.5) * dx, (j + 0.5) * dx, this.Ly);
      if (P.csvSlope > 0) {   // 縦断勾配を指定した値に傾け直す（各横断の最低点の回帰直線で判定）
        var ax = 0, ay = 0, axx = 0, axy = 0;
        for (i = 0; i < nx; i++) { var mn0 = 1e30; for (j = 0; j < ny; j++) if (z[j * nx + i] < mn0) mn0 = z[j * nx + i]; var xx = (i + 0.5) * dx; ax += xx; ay += mn0; axx += xx * xx; axy += xx * mn0; }
        var b0 = (nx * axy - ax * ay) / (nx * axx - ax * ax), xm = ax / nx, add = -P.csvSlope - b0;
        for (j = 0; j < ny; j++) for (i = 0; i < nx; i++) z[j * nx + i] += add * ((i + 0.5) * dx - xm);
      }
      var sx = 0, sy = 0, sxx = 0, sxy = 0;
      for (i = 0; i < nx; i++) {
        var mn = 1e30, jm = 0;
        for (j = 0; j < ny; j++) if (z[j * nx + i] < mn) { mn = z[j * nx + i]; jm = j; }
        colmin[i] = mn; yc[i] = (jm + 0.5) * dx;
        var xi = (i + 0.5) * dx; sx += xi; sy += mn; sxx += xi * xi; sxy += xi * mn;
      }
      var bb = (nx * sxy - sx * sy) / (nx * sxx - sx * sx);
      this.trB = bb; this.trA = (sy - bb * sx) / nx;
      this.I = Math.max(1e-3, -bb);
      if (bb > 0) this.warn.push('読み込んだ地形は下流ほど高くなっています。「流下方向を逆にする」を確認してください。');
      // 谷底（侵食できる所）の印が付いていれば、それを使う（計算後の河床から続けるとき。初期の擾乱も足さない）
      var M = (T.mask && T.mask.length === T.rows.length && T.mask[0] && T.mask[0].length === T.rows[0].length) ? T.mask : null;
      for (i = 0; i < nx; i++) {
        var nf = 0;
        for (j = 0; j < ny; j++) {
          c = j * nx + i;
          var isF = M ? M[Math.min(M.length - 1, Math.max(0, Math.floor((this.Ly - (j + 0.5) * dx) / T.cs)))][Math.min(M[0].length - 1, Math.floor((i + 0.5) * dx / T.cs))] === 1
            : z[c] <= colmin[i] + floorH;
          if (isF) {
            fl[c] = 1; kind[c] = 0; nf++;
            var fixedEndC = (i < 2 || i >= nx - 2);
            if (!fixedEndC && !M) z[c] += P.perturb * noise[c];
            zfix[c] = fixedEndC ? z[c] : z[c] - P.erodDepth;
          } else { kind[c] = 1; zfix[c] = z[c]; }
        }
        Bx[i] = Math.max(1, nf) * dx;
      }
    } else {
      var I0 = this.I = P.slope;
      this.trA = this.Lx * I0; this.trB = -I0;
      var VV = P.valley || {}, chD = Math.max(0, VV.chD || 0), relief = Math.max(0, VV.relief || 0), wR = Math.max(0, VV.wallRough || 0);
      var fans = (VV.fans || []).filter(function (f) { return isFinite(f[0]) && f[2] > 0 && f[3] > 0; });
      this.shaped = chD > 0 || relief > 0 || fans.length > 0 || wR > 0;
      var bare = this.bare = new Uint8Array(N);
      // 澪筋の位置（中心線からの法線方向の距離、左が正）：谷底の中を滑らかに蛇行させる
      var ych = new Float64Array(nx), meanB = 0;
      for (i = 0; i < nx; i++) meanB += Bx[i] / nx;
      var chW = VV.chW > 0 ? VV.chW : 0.3 * meanB;
      if (chD > 0) {
        var r2 = rng(P.seed * 31 + 7), ph1 = r2() * 6.283, ph2 = r2() * 6.283, wlc = VV.chWl > 0 ? VV.chWl : Math.max(60, 4 * meanB);
        for (i = 0; i < nx; i++) {
          var xs_ = (i + 0.5) * dx, sv = 0.72 * Math.sin(2 * Math.PI * xs_ / wlc + ph1) + 0.28 * Math.sin(2 * Math.PI * xs_ * 2.37 / wlc + ph2);
          var lo_ = -this.hR[i] + chW / 2 + dx, hi_ = this.hL[i] - chW / 2 - dx;
          ych[i] = hi_ > lo_ ? (lo_ + hi_) / 2 + (hi_ - lo_) / 2 * 0.85 * sv : (this.hL[i] - this.hR[i]) / 2;
        }
        var wsm = Math.max(1, Math.round(chW / dx)), tmpY = new Float64Array(nx);   // 幅が急に変わる所でも澪筋がつながるように平滑化
        for (var pass = 0; pass < 2; pass++) {
          for (i = 0; i < nx; i++) { var sy_ = 0, ny_ = 0; for (k = Math.max(0, i - wsm); k <= Math.min(nx - 1, i + wsm); k++) { sy_ += ych[k]; ny_++; } tmpY[i] = sy_ / ny_; }
          ych.set(tmpY);
        }
      }
      this.ych = ych;
      var self = this, cs0 = P.crossSlope;
      // 谷底の面（澪筋へ向かう横断勾配）
      var plane = function (i_, o_) { return (self.Lx - (i_ + 0.5) * dx) * I0 + cs0 * Math.abs(o_ - ych[i_]); };
      for (i = 0; i < nx; i++) {
        var zc = (this.Lx - (i + 0.5) * dx) * I0, hLi = this.hL[i], hRi = this.hR[i];
        colmin[i] = zc;
        for (j = 0; j < ny; j++) {
          c = j * nx + i;
          var offs = ((j + 0.5) * dx - yc[i]) * cosT[i];   // 中心線からの法線方向距離（左が正）
          var xw = (i + 0.5) * dx, yw = (j + 0.5) * dx;
          if (offs < hLi && offs > -hRi) {
            fl[c] = 1; kind[c] = 0;
            var fixedEnd = (i < 2 || i >= nx - 2);
            var zz = plane(i, offs);
            if (chD > 0) {
              var uu_ = Math.abs(offs - ych[i]) / (chW / 2);
              zz -= chD * (uu_ < 0.55 ? 1 : uu_ < 1 ? 0.5 * (1 + Math.cos(Math.PI * (uu_ - 0.55) / 0.45)) : 0);
            }
            if (relief > 0) zz += relief * (0.7 * vn2(xw / 16, yw / 9, P.seed * 13 + 3) + 0.3 * vn2(xw / 6, yw / 6, P.seed * 13 + 4));
            z[c] = zz + (fixedEnd ? 0 : P.perturb * noise[c]);
          } else {
            var dw = offs >= hLi ? offs - hLi : -hRi - offs, eo = offs >= hLi ? hLi : -hRi;
            var zw_ = plane(i, eo) + Math.min(P.sideGrad * dw, wallMax);
            if (wR > 0) zw_ += wR * 0.6 * vn2(xw / 11, yw / 11, P.seed * 13 + 5) * Math.min(1, dw / (3 * dx));
            z[c] = zw_; kind[c] = 1;
          }
        }
      }
      // 支渓の扇状地（崖錐）：岸に頂点を置いた円錐を、谷底の面の上に重ねる
      for (q = 0; q < fans.length; q++) {
        var Fa = fans[q], fi = Math.max(0, Math.min(nx - 1, Math.round(Fa[0] / dx - 0.5))), sd = Fa[1] >= 0 ? 1 : -1;
        var eo2 = sd > 0 ? this.hL[fi] : -this.hR[fi], ax_ = (fi + 0.5) * dx, ay_ = yc[fi] + eo2 / cosT[fi], RR = Fa[2], HH = Fa[3];
        var i0f = Math.max(0, Math.floor((ax_ - RR) / dx)), i1f = Math.min(nx - 1, Math.ceil((ax_ + RR) / dx));
        for (i = i0f; i <= i1f; i++) for (j = 0; j < ny; j++) {
          c = j * nx + i;
          var rr_ = Math.hypot((i + 0.5) * dx - ax_, (j + 0.5) * dx - ay_); if (rr_ >= RR) continue;
          var off2 = ((j + 0.5) * dx - yc[i]) * cosT[i], oc = Math.max(-this.hR[i], Math.min(this.hL[i], off2));
          var zf = plane(i, oc) + HH * Math.pow(1 - rr_ / RR, 1.5);
          if (zf > z[c]) { z[c] = zf; if (kind[c] === 1) bare[c] = 1; }
        }
      }
      for (i = 0; i < nx; i++) {
        var fixedE = (i < 2 || i >= nx - 2), mnc = 1e30;
        for (j = 0; j < ny; j++) {
          c = j * nx + i;
          if (kind[c] === 0) { zfix[c] = fixedE ? z[c] : z[c] - P.erodDepth; if (z[c] < mnc) mnc = z[c]; }
          else zfix[c] = z[c];
        }
        if (this.shaped && mnc < 1e29) colmin[i] = mnc;
      }
    }
    var I = this.I;

    // ---- 堰堤（河道の中心線に直角＝両岸に直角に置く）
    this.dams = [];
    // 中心線の向き：CSV では最深部の位置を平滑化して使う
    var ycs = new Float64Array(nx);
    var win = T ? Math.max(3, Math.round(15 / dx * 2)) : 0;
    for (i = 0; i < nx; i++) {
      if (!win) { ycs[i] = yc[i]; continue; }
      var sw = 0, nw = 0; for (k = Math.max(0, i - win); k <= Math.min(nx - 1, i + win); k++) { sw += yc[k]; nw++; } ycs[i] = sw / nw;
    }
    this.ycs = ycs;
    var dl = (P.dams || []).slice().sort(function (a, b) { return a.x - b.x; });
    var sideG = T ? 0.3 : Math.max(0.2, P.sideGrad);
    for (q = 0; q < dl.length; q++) {
      var D = dl[q];
      if (!(D.H > 0)) continue;
      var id = Math.round(D.x / dx - 0.5);
      if (D.x < 0 || D.x > this.Lx) { this.warn.push('x = ' + D.x + ' m の堰堤は、地形の範囲（0〜' + Math.round(this.Lx) + ' m）の外なので置きませんでした。'); continue; }
      if (id < 4) id = 4; if (id > nx - 5) id = nx - 5;
      // 中心 y・回転角・天端標高・堤長・水通しの中心は、与えられていればその値を使う（計算後の河床から続けるときに堰堤を動かさない）
      var cx = (id + 0.5) * dx, cy = isFinite(D.y) ? D.y : ycs[id];
      var i1d = Math.min(nx - 1, id + 1), i0d = Math.max(0, id - 1);
      var th = isFinite(D.ang) ? D.ang * Math.PI / 180 : Math.atan2(ycs[i1d] - ycs[i0d], (i1d - i0d) * dx);
      var tx = Math.cos(th), ty = Math.sin(th), nxv = -ty, nyv = tx;      // t：流下方向、n：堰堤の軸方向（左岸が正）
      var zbase = T ? colmin[id] : (this.Lx - cx) * I;
      var crest = isFinite(D.zc) ? D.zc : zbase + D.H, wingTop = crest + P.wingHeight;
      zbase = crest - D.H;
      var spFix = isFinite(D.so) && D.b > 0;                                // 水通しの位置を固定
      var bn = D.b > 0 ? (spFix ? D.b : Math.min(D.b, Bx[id])) : Bx[id];   // 水通し幅（堰堤軸に沿う長さ）
      var hw0 = this.hL ? Math.max(this.hL[id], this.hR[id]) : Bx[id] / 2;
      var ext = D.len > 0 ? D.len / 2 : hw0 + (D.H + P.wingHeight) / sideG + 2 * dx;   // 袖を谷壁へ差し込む範囲
      var thick = (Math.abs(tx) + Math.abs(ty)) * dx / 2 + 1e-9;           // 4近傍で水が漏れない最小の厚さ
      var cells = [], spill = [], clash = false, sMin = 1e30, sMax = -1e30, spMin = 1e30, spMax = -1e30, zlow = 1e30;
      var iA = Math.max(1, Math.floor((cx - ext - dx) / dx)), iB = Math.min(nx - 2, Math.ceil((cx + ext + dx) / dx));
      for (j = 0; j < ny; j++) for (i = iA; i <= iB; i++) {
        c = j * nx + i;
        var px = (i + 0.5) * dx - cx, py = (j + 0.5) * dx - cy;
        var al = px * tx + py * ty, sl = px * nxv + py * nyv;
        if (Math.abs(al) > thick || Math.abs(sl) > ext) continue;
        var inSp = spFix ? Math.abs(sl - D.so) <= bn / 2 + 1e-9 : fl[c] === 1 && (!(D.b > 0) || Math.abs(sl) <= bn / 2 + 1e-9);
        var top = inSp ? crest : wingTop;
        if (top <= z[c] && !(spFix && inSp)) continue;    // 固定した水通しは、地盤が天端より高くても置く（穴をあけない）
        if (kind[c] === 2) clash = true;
        cells.push(c); spill.push(inSp ? 1 : 0);
        if (sl < sMin) sMin = sl; if (sl > sMax) sMax = sl;
        if (inSp) { if (sl < spMin) spMin = sl; if (sl > spMax) spMax = sl; }
        if (z[c] < zlow) zlow = z[c];
      }
      if (clash || !cells.length) { this.warn.push('x = ' + D.x + ' m の堰堤は、ほかの堰堤と重なるので置きませんでした。'); continue; }
      var j0 = ny, j1 = -1;
      for (k = 0; k < cells.length; k++) {
        c = cells[k]; var tp = spill[k] ? crest : wingTop;
        z[c] = tp; zfix[c] = tp; kind[c] = 2;
        this.noSl[c] = 1; this.noSl[c - 1] = 1; this.noSl[c + 1] = 1; if (c - nx >= 0) this.noSl[c - nx] = 1; if (c + nx < N) this.noSl[c + nx] = 1;
        if (spill[k] && c % nx === id) { var jj2 = (c - id) / nx; if (jj2 < j0) j0 = jj2; if (jj2 > j1) j1 = jj2; }
      }
      if (spMin > spMax) { spMin = 0; spMax = 0; }
      this.dams.push({ i: id, x: cx, cx: cx, cy: cy, th: th, tx: tx, ty: ty, nx: nxv, ny: nyv, ext: ext, thick: thick, H: D.H, b: bn,
        crest: crest, wingTop: wingTop, zbase: zbase, zlow: zlow, j0: j0, j1: j1, cells: cells, spill: spill,
        sMin: sMin - dx / 2, sMax: sMax + dx / 2, spMin: spMin - dx / 2, spMax: spMax + dx / 2 });
    }
    // 各セルが何番目の堰堤の上流側にあるか（堆砂量の集計・初期湛水に使う）
    var region = this.region = new Int16Array(N), nd = this.dams.length;
    for (j = 0; j < ny; j++) for (i = 0; i < nx; i++) {
      c = j * nx + i; var rq = nd;
      for (q = 0; q < nd; q++) {
        var Dq = this.dams[q], ppx = (i + 0.5) * dx - Dq.cx, ppy = (j + 0.5) * dx - Dq.cy;
        var slq = ppx * Dq.nx + ppy * Dq.ny, alq = ppx * Dq.tx + ppy * Dq.ty;
        if (Math.abs(slq) <= Dq.ext ? alq <= Dq.thick : ppx <= 0) { rq = q; break; }
      }
      region[c] = rq;
    }
    // 堰堤の上流をはじめから満砂にする（堆砂勾配 = 元の勾配 × prefill）
    if (P.prefill > 0) {
      var relF = (!T && P.valley && P.valley.relief > 0) ? P.valley.relief : 0;
      for (c = 0; c < N; c++) {
        if (kind[c] === 2) continue;
        var rqf = region[c]; if (rqf >= nd) continue;
        var Df = this.dams[rqf], pxf = ((c % nx) + 0.5) * dx - Df.cx, pyf = (((c / nx) | 0) + 0.5) * dx - Df.cy;
        var up = -(pxf * Df.tx + pyf * Df.ty); if (up < 0) continue;
        var zt = Df.crest + up * I * P.prefill;
        if (relF > 0) zt += relF * 0.6 * vn2(((c % nx) + 0.5) * dx / 14, (((c / nx) | 0) + 0.5) * dx / 8, P.seed * 13 + 9) * Math.min(1, up / 6);
        if (z[c] < zt) {
          if (kind[c] === 1) { kind[c] = 0; fl[c] = 1; }   // 谷壁の裾も埋める（その下の元の斜面は侵食しない）
          z[c] = zt;
        }
      }
      for (i = 0; i < nx; i++) { var mnp = 1e30; for (j = 0; j < ny; j++) { c = j * nx + i; if (kind[c] === 0 && z[c] < mnp) mnp = z[c]; } if (mnp < 1e29) colmin[i] = mnp; }
      this.shaped = true;
    }
    // ---- 侵食できる側岸：谷壁のうち、谷底の高さから erodH までの部分（段丘・崖錐）
    var BK = P.bank || {}, eH = BK.erodH > 0 ? BK.erodH : 0;
    if (typeof meanB === 'undefined' || !(meanB > 0)) { var meanB = 0; for (i = 0; i < nx; i++) meanB += Bx[i] / nx; }
    this.bankC = new Uint8Array(N);
    if (eH > 0) {
      for (i = 0; i < nx; i++) {
        var fmax = -1e30, fmin = 1e30;
        for (j = 0; j < ny; j++) { c = j * nx + i; if (kind[c] === 0) { if (z[c] > fmax) fmax = z[c]; if (z[c] < fmin) fmin = z[c]; } }
        if (fmax < -1e29) continue;
        for (j = 0; j < ny; j++) {
          c = j * nx + i;
          if (kind[c] !== 1 || z[c] > fmax + eH) continue;
          kind[c] = 0; this.bankC[c] = 1;
          zfix[c] = Math.min(z[c], (i < 2 || i >= nx - 2) ? z[c] : fmin - P.erodDepth);
        }
      }
    }
    // ---- 巨礫：水みち（初期の澪筋、なければ谷底）の両岸に、間隔 L で動かない礫を置く
    var BL = P.boulders || {}, bL = BL.L > 0 ? BL.L : 0;
    this.fixedC = new Uint8Array(N);
    if (bL > 0) {
      var bsz = Math.max(1, Math.min(4, BL.size | 0 || 1)), nb_ = 0, simB = this;
      var edgeAt = function (i_, side) {     // 横断 i_ の水みちの岸のセル（side +1：左岸、-1：右岸）
        var jj, best = -1;
        if (simB.ych && P.valley && P.valley.chD > 0 && !T) {
          var chW2 = (P.valley.chW > 0 ? P.valley.chW : 0.3 * meanB) / 2, target = simB.ych[i_] + side * chW2;
          var yt = yc[i_] + target / cosT[i_]; jj = Math.round(yt / dx - 0.5);
          return (jj >= 0 && jj < ny && kind[jj * nx + i_] === 0) ? jj : -1;
        }
        if (side > 0) { for (jj = ny - 1; jj >= 0; jj--) if (kind[jj * nx + i_] === 0 && fl[jj * nx + i_]) { best = jj; break; } }
        else { for (jj = 0; jj < ny; jj++) if (kind[jj * nx + i_] === 0 && fl[jj * nx + i_]) { best = jj; break; } }
        return best;
      };
      for (var side = -1; side <= 1; side += 2) {
        var x0b = bL / 2 + (BL.mode === 'zigzag' && side < 0 ? bL / 2 : 0);
        for (var xb = x0b; xb < this.Lx - 3 * dx; xb += bL) {
          var ib = Math.round(xb / dx - 0.5); if (ib < 3 || ib > nx - 4) continue;
          var jb = edgeAt(ib, side); if (jb < 0) continue;
          for (var di = 0; di < bsz; di++) for (var dj = 0; dj < bsz; dj++) {
            var ii = ib + di - (bsz >> 1), jj2 = jb + side * dj;
            if (ii < 1 || ii >= nx - 1 || jj2 < 0 || jj2 >= ny) continue;
            var cb = jj2 * nx + ii; if (kind[cb] !== 0) continue;
            this.fixedC[cb] = 1; zfix[cb] = z[cb]; nb_++;
          }
        }
      }
      this.nBoulder = nb_;
    }
    // 表示用の基準高さ
    var zs = 0, ns_ = 0, zlo = 1e30, zhi = -1e30;
    for (c = 0; c < N; c++) { if (fl[c]) { zs += z[c]; ns_++; } if (z[c] < zlo) zlo = z[c]; if (z[c] > zhi) zhi = z[c]; }
    this.zmid = ns_ ? zs / ns_ : 0; this.zlo = zlo; this.zhi = zhi;
    this.Bf = Bx[0];
    this.z0.set(z);

    // ---- 土砂の層構造
    var La = this.La, Ls = this.Ls;
    maxH = 0; for (q = 0; q < this.dams.length; q++) if (this.dams[q].H > maxH) maxH = this.dams[q].H;
    var thMax = 0; for (c = 0; c < N; c++) if (z[c] - zfix[c] > thMax) thMax = z[c] - zfix[c];
    var NL = this.NL = Math.min(100, Math.max(Math.ceil((P.erodDepth + maxH + 3) / Ls) + 2, Math.ceil(thMax / Ls) + 2));
    this.A = new Float64Array(N * K);            // 交換層内の粒径別体積（層厚換算） [m]
    this.SF = new Float32Array(N * NL * K);      // 堆積層の粒径別割合
    this.nl = new Int16Array(N);                 // 堆積層の層数
    this.topTh = new Float32Array(N);            // 最上層の厚さ [m]
    this.qbx = new Float64Array(N * K); this.qby = new Float64Array(N * K); this.gk = new Float64Array(N * K);
    this.act = new Uint8Array(N);
    this.tcs = new Float32Array(N).fill(NaN);     // 平均粒径の無次元限界掃流力 τ*cm（表示用）
    this.SFx = new Float64Array(NFx * K); this.SFy = new Float64Array(NFy * K);
    this.scS = new Float64Array(N * K);
    this.dA = new Float64Array(K);
    this.kap = new Float64Array(N); this.kapW = new Float64Array(N); this.secMask = new Uint8Array(N); this.uf = new Float64Array(N); this.vf = new Float64Array(N);   // 流線の曲率 1/r（左曲がりが正） [1/m]
    this.Txx = new Float64Array(N); this.Txy = new Float64Array(N);    // 2次流による運動量の分散項（セル中心）
    this.nHase = 0;                               // 長谷川モデルの崩落回数
    this.dzs = new Float64Array(N);               // 1 回の河床変動計算での河床変化（水際の側岸侵食に使う）
    this.Cc = new Float64Array(N * K);            // 侵食・堆積速度式：流れの中の粒径別土砂量 C_k h（実質体積の層厚） [m]
    this.Cn = new Float64Array(N * K);
    this.sedWaterInit = 0;
    for (c = 0; c < N; c++) {
      var del = z[c] - zfix[c];
      if (del <= 1e-9) continue;
      var a = Math.min(del, La);
      for (k = 0; k < K; k++) this.A[c * K + k] = a * f0[k];
      var rem = del - a;
      if (rem > 1e-9) {
        var nl = Math.min(NL, Math.ceil(rem / Ls - 1e-9));
        this.nl[c] = nl; this.topTh[c] = rem - (nl - 1) * Ls;
        for (var l = 0; l < nl; l++) for (k = 0; k < K; k++) this.SF[(c * NL + l) * K + k] = f0[k];
      }
    }

    // ---- 初期の水（等流の薄層流＋堰堤上流の湛水）
    var Q0 = this.Qat(0);
    var h = this.h, qx = this.qx;
    if (Q0 > 0.01) {
      for (i = 0; i < nx; i++) {
        var q0 = Q0 / Bx[i];
        var hn = Math.pow(q0 * P.manning / Math.sqrt(Math.max(I, 1e-4)), 0.6), un = q0 / hn;
        for (j = 0; j < ny; j++) {
          c = j * nx + i; if (kind[c] !== 0) continue;
          if (T || this.shaped) { var hw = colmin[i] + hn - z[c]; if (hw > 0) { h[c] = hw; qx[c] = hw * un; } }
          else { h[c] = hn; qx[c] = q0; }
        }
      }
    }
    for (c = 0; c < N; c++) {   // 堰堤の上流を天端まで湛水させておく
      var rq2 = this.region[c]; if (rq2 >= this.dams.length || kind[c] === 2) continue;
      var cr = this.dams[rq2].crest;
      if (z[c] < cr && cr - z[c] > h[c]) { h[c] = cr - z[c]; qx[c] = 0; this.qy[c] = 0; }
    }

    // ---- 時刻・収支
    this.t = 0;          // 出水の時刻（河床変動の時間軸） [s]
    this.th = 0;         // 流れの計算時間 [s]
    this.nstep = 0;
    this.dt = 0.01;
    this.phase = Q0 > 0.01 ? 'spinup' : 'run';
    this.spinT = 0;
    this.Qin = Q0; this.Qout = 0; this.QsIn = 0; this.QsOut = 0;
    this.Vwin = 0; this.Vwout = 0;
    this.sedIn = 0; this.sedOut = 0;     // 実質体積 [m3]
    this.sedOutK = new Float64Array(K);
    this.nLimW = 0; this.nLimS = 0;
    this.series = { t: [], Qin: [], Qout: [], QsIn: [], QsOut: [], stor: [] };
    this._nextSample = 0;
    this.Uref = Q0 > 0.01 ? Math.max(0.5, (Q0 / this.Bf) / Math.pow((Q0 / this.Bf) * P.manning / Math.sqrt(Math.max(I, 1e-4)), 0.6)) : 1;
  };

  // ------------------------------------------------------------ ハイドログラフ
  Sim.prototype.duration = function () {
    var H = this.P.hydro;
    if (H.mode === 'table' && H.table && H.table.length > 1) return H.table[H.table.length - 1][0];
    return H.duration;
  };
  function interp(tab, t, col) {
    var n = tab.length;
    if (t <= tab[0][0]) return tab[0][col];
    if (t >= tab[n - 1][0]) return tab[n - 1][col];
    for (var i = 1; i < n; i++) if (t <= tab[i][0]) {
      var w = (t - tab[i - 1][0]) / Math.max(1e-9, tab[i][0] - tab[i - 1][0]);
      return tab[i - 1][col] * (1 - w) + tab[i][col] * w;
    }
    return tab[n - 1][col];
  }
  Sim.prototype.Qat = function (t) {
    var H = this.P.hydro;
    if (H.mode === 'table' && H.table && H.table.length > 1) return Math.max(0, interp(H.table, t, 1));
    if (H.mode === 'const') return H.Qpeak;
    if (t <= 0) return H.Qbase;
    if (t < H.tPeak) return H.Qbase + (H.Qpeak - H.Qbase) * t / H.tPeak;
    if (t < H.duration) return H.Qpeak + (H.Qbase - H.Qpeak) * (t - H.tPeak) / Math.max(1e-9, H.duration - H.tPeak);
    return H.Qbase;
  };
  // 供給土砂量（実質体積） [m3/s]
  Sim.prototype.QsAt = function (t) {
    var H = this.P.hydro;
    if (H.mode === 'table' && H.table && H.table.length > 1 && H.table[0].length > 2 && isFinite(H.table[0][2]))
      return Math.max(0, interp(H.table, t, 2));
    return H.Cs * this.Qat(t);
  };

  // ------------------------------------------------------------------ 1 ステップ
  Sim.prototype.step = function () {
    var nx = this.nx, ny = this.ny, N = this.N, dx = this.dx, P = this.P;
    var h = this.h, qx = this.qx, qy = this.qy, u = this.u, v = this.v, eta = this.eta, z = this.z;
    var dhx = this.dhx, dex = this.dex, dhy = this.dhy, dey = this.dey, noSl = this.noSl;
    var FxM = this.FxM, FxU = this.FxU, FxV = this.FxV, FxCL = this.FxCL, FxCR = this.FxCR;
    var GyM = this.GyM, GyU = this.GyU, GyV = this.GyV, GyCL = this.GyCL, GyCR = this.GyCR;
    var HDRY = 1e-3, HS = 0.01, UMAX = 30, g = G, hg = 0.5 * G;
    var i, j, c, f, L, R, row, frow;
    var second = P.order !== 1;

    // ---- 流入境界（上流端の谷底セルに流量を配分）
    var Q = this.Qat(this.t);
    this.Qin = Q;
    var ns = this.ns, nyb = this.nyb, qin = this.qin, uin = this.uin;
    var wsum = 0, sqI = Math.sqrt(Math.max(this.I, 1e-4)), nman = P.manning;
    var fl = this.fl;
    for (j = 0; j < ny; j++) { var w = fl[j * nx] ? Math.pow(h[j * nx] + 0.05, 5 / 3) : 0; qin[j] = w; wsum += w; }
    for (j = 0; j < ny; j++) {
      var qj = Q * qin[j] / wsum / dx;                 // 単位幅流量 [m2/s]
      qin[j] = qj;
      uin[j] = qj > 0 ? qj / Math.pow(qj * nman / sqI, 0.6) : 0;
    }

    // ---- 2次流：流線の曲率（前のステップの面の流束から求める）
    var SEC = P.sec || {}, secOn = !!(SEC.bed || SEC.mom);
    if (secOn) this.curvature();

    // ---- 原始変数と勾配
    for (c = 0; c < N; c++) {
      var hh = h[c];
      if (hh > HDRY) {
        var uu = qx[c] / hh, vv = qy[c] / hh;
        if (uu > UMAX) { uu = UMAX; qx[c] = uu * hh; } else if (uu < -UMAX) { uu = -UMAX; qx[c] = uu * hh; }
        if (vv > UMAX) { vv = UMAX; qy[c] = vv * hh; } else if (vv < -UMAX) { vv = -UMAX; qy[c] = vv * hh; }
        u[c] = uu; v[c] = vv;
      } else { u[c] = 0; v[c] = 0; qx[c] = 0; qy[c] = 0; }
      eta[c] = z[c] + hh;
    }
    if (second) {
      var a, b, Idx = this.I * dx;
      for (j = 0; j < ny; j++) {
        row = j * nx;
        for (i = 0; i < nx; i++) {
          c = row + i;
          var sx = 0, ex = 0, sy = 0, ey = 0, hc = h[c];
          if (hc > HS && !noSl[c]) {
            if (i > 0 && i < nx - 1 && h[c - 1] > HS && h[c + 1] > HS) {
              a = hc - h[c - 1]; b = h[c + 1] - hc;
              sx = a * b <= 0 ? 0 : (Math.abs(a) < Math.abs(b) ? a : b);
              a = eta[c] - eta[c - 1]; b = eta[c + 1] - eta[c];
              ex = a * b <= 0 ? 0 : (Math.abs(a) < Math.abs(b) ? a : b);
            } else if (i === nx - 1 && h[c - 1] > HS) {
              // 下流端：初期勾配で続く仮想セル（同じ水深）を仮定
              a = eta[c] - eta[c - 1]; b = -Idx;
              ex = a * b <= 0 ? 0 : (Math.abs(a) < Math.abs(b) ? a : b);
            }
            if (j > 0 && j < ny - 1 && h[c - nx] > HS && h[c + nx] > HS) {
              a = hc - h[c - nx]; b = h[c + nx] - hc;
              sy = a * b <= 0 ? 0 : (Math.abs(a) < Math.abs(b) ? a : b);
              a = eta[c] - eta[c - nx]; b = eta[c + nx] - eta[c];
              ey = a * b <= 0 ? 0 : (Math.abs(a) < Math.abs(b) ? a : b);
            }
          }
          dhx[c] = sx; dex[c] = ex; dhy[c] = sy; dey[c] = ey;
        }
      }
    }

    // ---- x 方向の数値流束（HLL + 静水圧再構築）
    var smax = 1e-6;
    var hL, hR, eL, eR, zL, zR, zs, mw, hLs, hRs, uL, uR, cL, cR, SL, SR, mL, mR, fL1, fR1, F0, F1, inv, sa;
    for (j = 0; j < ny; j++) {
      row = j * nx; frow = j * (nx + 1);
      // 上流端：壁（流入は生成項として与える）
      FxM[frow] = 0; FxU[frow] = 0; FxV[frow] = 0; FxCL[frow] = 0; FxCR[frow] = hg * h[row] * h[row];
      for (i = 1; i < nx; i++) {
        L = row + i - 1; R = L + 1; f = frow + i;
        hL = h[L] + 0.5 * dhx[L]; hR = h[R] - 0.5 * dhx[R];
        if (hL < 1e-9 && hR < 1e-9) { FxM[f] = 0; FxU[f] = 0; FxV[f] = 0; FxCL[f] = 0; FxCR[f] = 0; continue; }
        eL = eta[L] + 0.5 * dex[L]; eR = eta[R] - 0.5 * dex[R];
        zL = eL - hL; zR = eR - hR;
        zs = zL > zR ? zL : zR; mw = eL < eR ? eL : eR; if (mw < zs) zs = mw;
        hLs = eL - zs; if (hLs > hL) hLs = hL;
        hRs = eR - zs; if (hRs > hR) hRs = hR;
        uL = u[L]; uR = u[R];
        if (hLs < 1e-9 && hRs < 1e-9) { F0 = 0; F1 = 0; }
        else {
          cL = Math.sqrt(g * hLs); cR = Math.sqrt(g * hRs);
          if (hLs < 1e-9) { SL = uR - 2 * cR; SR = uR + cR; }
          else if (hRs < 1e-9) { SL = uL - cL; SR = uL + 2 * cL; }
          else { SL = uL - cL < uR - cR ? uL - cL : uR - cR; SR = uL + cL > uR + cR ? uL + cL : uR + cR; }
          mL = hLs * uL; mR = hRs * uR;
          fL1 = mL * uL + hg * hLs * hLs; fR1 = mR * uR + hg * hRs * hRs;
          if (SL >= 0) { F0 = mL; F1 = fL1; }
          else if (SR <= 0) { F0 = mR; F1 = fR1; }
          else {
            inv = 1 / (SR - SL);
            F0 = (SR * mL - SL * mR + SL * SR * (hRs - hLs)) * inv;
            F1 = (SR * fL1 - SL * fR1 + SL * SR * (mR - mL)) * inv;
          }
          sa = -SL > SR ? -SL : SR; if (sa > smax) smax = sa;
        }
        FxM[f] = F0; FxU[f] = F1; FxV[f] = F0 * (F0 >= 0 ? v[L] : v[R]);
        FxCL[f] = hg * (hL + hLs) * (zs - zL);
        FxCR[f] = hg * (hR + hRs) * (zs - zR);
      }
      // 下流端：自由流出（流出のみ）
      L = row + nx - 1; f = frow + nx;
      hL = h[L]; uL = u[L];
      if (hL > HDRY && uL > 0) {
        FxM[f] = hL * uL; FxU[f] = hL * uL * uL + hg * hL * hL; FxV[f] = FxM[f] * v[L];
        sa = uL + Math.sqrt(g * hL); if (sa > smax) smax = sa;
      } else { FxM[f] = 0; FxU[f] = hg * hL * hL; FxV[f] = 0; }
      FxCL[f] = 0; FxCR[f] = 0;
    }

    // ---- y 方向の数値流束
    var vL, vR;
    for (i = 0; i < nx; i++) {              // 側方境界：壁
      GyM[i] = 0; GyU[i] = 0; GyV[i] = 0; GyCL[i] = 0; GyCR[i] = hg * h[i] * h[i];
      f = ny * nx + i; c = (ny - 1) * nx + i;
      GyM[f] = 0; GyU[f] = 0; GyV[f] = 0; GyCR[f] = 0; GyCL[f] = hg * h[c] * h[c];
    }
    for (j = 1; j < ny; j++) {
      row = j * nx;
      for (i = 0; i < nx; i++) {
        R = row + i; L = R - nx; f = R;
        hL = h[L] + 0.5 * dhy[L]; hR = h[R] - 0.5 * dhy[R];
        if (hL < 1e-9 && hR < 1e-9) { GyM[f] = 0; GyU[f] = 0; GyV[f] = 0; GyCL[f] = 0; GyCR[f] = 0; continue; }
        eL = eta[L] + 0.5 * dey[L]; eR = eta[R] - 0.5 * dey[R];
        zL = eL - hL; zR = eR - hR;
        zs = zL > zR ? zL : zR; mw = eL < eR ? eL : eR; if (mw < zs) zs = mw;
        hLs = eL - zs; if (hLs > hL) hLs = hL;
        hRs = eR - zs; if (hRs > hR) hRs = hR;
        vL = v[L]; vR = v[R];
        if (hLs < 1e-9 && hRs < 1e-9) { F0 = 0; F1 = 0; }
        else {
          cL = Math.sqrt(g * hLs); cR = Math.sqrt(g * hRs);
          if (hLs < 1e-9) { SL = vR - 2 * cR; SR = vR + cR; }
          else if (hRs < 1e-9) { SL = vL - cL; SR = vL + 2 * cL; }
          else { SL = vL - cL < vR - cR ? vL - cL : vR - cR; SR = vL + cL > vR + cR ? vL + cL : vR + cR; }
          mL = hLs * vL; mR = hRs * vR;
          fL1 = mL * vL + hg * hLs * hLs; fR1 = mR * vR + hg * hRs * hRs;
          if (SL >= 0) { F0 = mL; F1 = fL1; }
          else if (SR <= 0) { F0 = mR; F1 = fR1; }
          else {
            inv = 1 / (SR - SL);
            F0 = (SR * mL - SL * mR + SL * SR * (hRs - hLs)) * inv;
            F1 = (SR * fL1 - SL * fR1 + SL * SR * (mR - mL)) * inv;
          }
          sa = -SL > SR ? -SL : SR; if (sa > smax) smax = sa;
        }
        GyM[f] = F0; GyV[f] = F1; GyU[f] = F0 * (F0 >= 0 ? u[L] : u[R]);
        GyCL[f] = hg * (hL + hLs) * (zs - zL);
        GyCR[f] = hg * (hR + hRs) * (zs - zR);
      }
    }

    // ---- 時間刻み（CFL 条件）
    var dt = P.cfl * dx / smax;
    if (dt > 0.5) dt = 0.5;
    if (this.dt > 0 && dt > 1.5 * this.dt) dt = 1.5 * this.dt;
    this.dt = dt;
    var cf = dt / dx;

    // ---- 水深の非負性を保つ流出制限
    var sc = this.sc, any = false, out, fw, fe, gs, gn;
    for (j = 0; j < ny; j++) {
      row = j * nx; frow = j * (nx + 1);
      for (i = 0; i < nx; i++) {
        c = row + i; fw = frow + i; fe = fw + 1; gs = c; gn = c + nx;
        out = 0;
        if (FxM[fe] > 0) out += FxM[fe];
        if (FxM[fw] < 0) out -= FxM[fw];
        if (GyM[gn] > 0) out += GyM[gn];
        if (GyM[gs] < 0) out -= GyM[gs];
        out *= cf;
        if (out > h[c]) { sc[c] = out > 0 ? h[c] / out : 0; any = true; } else sc[c] = 1;
      }
    }
    if (any) {
      this.nLimW++;
      var r;
      for (j = 0; j < ny; j++) {
        row = j * nx; frow = j * (nx + 1);
        for (i = 1; i <= nx; i++) {
          f = frow + i;
          r = FxM[f] > 0 ? sc[row + i - 1] : (i < nx ? sc[row + i] : 1);
          if (r < 1) { FxM[f] *= r; FxU[f] *= r; FxV[f] *= r; }
        }
      }
      for (j = 1; j < ny; j++) {
        row = j * nx;
        for (i = 0; i < nx; i++) {
          f = row + i;
          r = GyM[f] > 0 ? sc[f - nx] : sc[f];
          if (r < 1) { GyM[f] *= r; GyU[f] *= r; GyV[f] *= r; }
        }
      }
    }

    // ---- 2次流による運動量の分散項（面の値は両側のセル中心の平均。どちらかが乾いていれば 0）
    var secM = !!SEC.mom, Txx = this.Txx, Txy = this.Txy, kmk = this.secMask, tE, tW, tN, tS, sE, sW, sN, sS;
    if (secM) this.secStress();

    // ---- 更新（勾配項・流入・摩擦）
    var n2 = this.n2, hn, qxn, qyn, Vv, fr, qoutSum = 0;
    for (j = 0; j < ny; j++) {
      row = j * nx; frow = j * (nx + 1);
      qoutSum += FxM[frow + nx];
      for (i = 0; i < nx; i++) {
        c = row + i; fw = frow + i; fe = fw + 1; gs = c; gn = c + nx;
        hh = h[c];
        hn = hh - cf * (FxM[fe] - FxM[fw] + GyM[gn] - GyM[gs]);
        qxn = qx[c] - cf * ((FxU[fe] + FxCL[fe]) - (FxU[fw] + FxCR[fw]) + GyU[gn] - GyU[gs]);
        qyn = qy[c] - cf * (FxV[fe] - FxV[fw] + (GyV[gn] + GyCL[gn]) - (GyV[gs] + GyCR[gs]));
        if (second && hh > 0) {
          qxn -= cf * g * hh * (dex[c] - dhx[c]);
          qyn -= cf * g * hh * (dey[c] - dhy[c]);
        }
        if (secM && kmk[c] > 0) {
          // x 運動量：∂Txx/∂x + ∂Txy/∂y、y 運動量：∂Txy/∂x + ∂Tyy/∂y（Tyy = −Txx）
          tE = (i < nx - 1 && kmk[c + 1] > 0) ? 0.5 * (Txx[c] + Txx[c + 1]) : 0; sE = (i < nx - 1 && kmk[c + 1] > 0) ? 0.5 * (Txy[c] + Txy[c + 1]) : 0;
          tW = (i > 0 && kmk[c - 1] > 0) ? 0.5 * (Txx[c] + Txx[c - 1]) : 0; sW = (i > 0 && kmk[c - 1] > 0) ? 0.5 * (Txy[c] + Txy[c - 1]) : 0;
          tN = (j < ny - 1 && kmk[c + nx] > 0) ? 0.5 * (Txx[c] + Txx[c + nx]) : 0; sN = (j < ny - 1 && kmk[c + nx] > 0) ? 0.5 * (Txy[c] + Txy[c + nx]) : 0;
          tS = (j > 0 && kmk[c - nx] > 0) ? 0.5 * (Txx[c] + Txx[c - nx]) : 0; sS = (j > 0 && kmk[c - nx] > 0) ? 0.5 * (Txy[c] + Txy[c - nx]) : 0;
          qxn -= cf * ((tE - tW) + (sN - sS));
          qyn -= cf * ((sE - sW) - (tN - tS));
        }
        if (i === 0 && qin[j] > 0) { hn += cf * qin[j]; qxn += cf * qin[j] * uin[j]; }
        if (hn < 0) hn = 0;
        if (hn > HDRY) {
          // 摩擦項は陰的（|q| の2次方程式を解く）
          Vv = Math.sqrt(qxn * qxn + qyn * qyn);
          fr = 0.5 * (1 + Math.sqrt(1 + 4 * dt * g * n2 * Vv / (hn * hn * Math.cbrt(hn))));
          qxn /= fr; qyn /= fr;
        } else { qxn = 0; qyn = 0; }
        h[c] = hn; qx[c] = qxn; qy[c] = qyn;
      }
    }
    this.Qout = qoutSum * dx;
    this.Vwin += Q * dt; this.Vwout += this.Qout * dt;
    this.th += dt; this.nstep++;

    // ---- 助走計算（水だけを流して流れを立ち上げる）
    if (this.phase === 'spinup') {
      this.spinT += dt;
      var tmin = 0.5 * this.Lx / this.Uref, tmax = 3 * this.Lx / this.Uref;
      if ((this.spinT > tmin && Math.abs(this.Qout - Q) < 0.08 * Q) || this.spinT > tmax) {
        this.phase = 'run'; this.Vwin = 0; this.Vwout = 0; this.Vw0 = this.waterVolume();
      }
      return;
    }
    if (this.phase === 'done') return;

    // ---- 河床変動
    var dtm = dt * Math.max(1, P.morfac);
    this.morph(dtm);
    this.t += dtm;
    if (this.t >= this._nextSample) {
      var S = this.series;
      S.t.push(this.t); S.Qin.push(Q); S.Qout.push(this.Qout); S.QsIn.push(this.QsIn); S.QsOut.push(this.QsOut);
      S.stor.push(this.storage());
      this._nextSample = this.t + this.duration() / 400;
    }
    if (this.t >= this.duration()) this.phase = 'done';
  };

  // --------------------------------------------------------- 河床変動 1 ステップ
  Sim.prototype.morph = function (dtm) {
    this.dzs.fill(0);
    if (this.P.transport === 'noneq') this.morphNoneq(dtm); else this.morphEq(dtm);
    this.dryBankErosion();
    this.hasegawaCollapse();
    this.collapse();
  };

  // 平衡流砂量（芦田・道上式）による河床変動
  Sim.prototype.morphEq = function (dtm) {
    var nx = this.nx, ny = this.ny, N = this.N, dx = this.dx, K = this.K, P = this.P, fxC = this.fixedC;
    var h = this.h, u = this.u, v = this.v, qxw = this.qx, qyw = this.qy, z = this.z, zfix = this.zfix;
    var A = this.A, qbx = this.qbx, qby = this.qby, gk = this.gk, act = this.act, d = this.d;
    var SFx = this.SFx, SFy = this.SFy, FxM = this.FxM, GyM = this.GyM;
    var g = G, n2 = this.n2, sg1 = this.sg1, La = this.La, tanPhi = this.tanPhi;
    var HSED = 0.02, invMu = (P.gravK == null ? 1 : P.gravK) / Math.sqrt(1.0 * 0.45);   // 1/sqrt(μs μk)
    var i, j, c, k, ck, row, frow, f, L, R;
    var cGrid = P.gridDamp == null ? 0.3 : P.gridDamp;

    var tcs = this.tcs, hideEg = P.hiding === 'egiazaroff';
    var SEC = P.sec || {}, secBed = !!SEC.bed, kapA = this.kap, Nst = SEC.Nstar || 11.5, hrMx = SEC.hrMax || 0.05;
    // ---- 各セルの粒径別掃流砂量（芦田・道上式）
    for (c = 0; c < N; c++) {
      act[c] = 0;
      var hh = h[c];
      if (hh < HSED || fxC[c]) continue;      // 巨礫のセルは動かない
      ck = c * K;
      var As = 0, dms = 0;
      for (k = 0; k < K; k++) { As += A[ck + k]; dms += A[ck + k] * d[k]; }
      if (As < 1e-6) continue;
      // 流速は面の質量流束の平均から求める（セル中心の運動量は凹凸のある河床上で
      // 2Δx の振動を持ちうるが、面の流束は連続式を満たすので滑らか）
      var ic = c % nx, jc = (c - ic) / nx, fwq = jc * (nx + 1) + ic;
      var uu = 0.5 * (FxM[fwq] + FxM[fwq + 1]) / hh, vv = 0.5 * (GyM[c] + GyM[c + nx]) / hh;
      var V2 = uu * uu + vv * vv;
      if (V2 < 1e-4) continue;
      var V = Math.sqrt(V2);
      var us2 = g * n2 * V2 / Math.cbrt(hh);          // u*²
      var dm = dms / As;
      var taum = us2 / (sg1 * g * dm);
      var arg = hh / (dm * (1 + 2 * taum));
      var den = 6.0 + 5.75 * Math.log10(arg > 0.2 ? arg : 0.2);
      var ue = V / (den > 2 ? den : 2);               // 有効摩擦速度 u*e
      var ue2 = ue * ue; if (ue2 > us2) { ue2 = us2; ue = Math.sqrt(us2); }
      var uss = Math.sqrt(us2), cap = 17 * ue2 * ue / (sg1 * g) / La;
      var anyq = false;
      var uscm2 = this.critUs2(c, dm, uu, vv);              // 平均粒径の限界摩擦速度の2乗
      tcs[c] = uscm2 / (sg1 * g * dm);
      // 格子スケールの凹凸（2Δx の振動）を抑える下限。Δx→0 で物理的な係数に戻る
      var epsMin = cGrid * dx / hh; if (epsMin > 4) epsMin = 4;
      // 2次流で河床付近の流れが曲率の内側へ向く（長田ら 1996 式14：u_nb/u_sb = N* h/r。n は流れの左向き）
      var secD = 0;
      if (secBed) { var hk0 = hh * kapA[c]; if (hk0 > hrMx) hk0 = hrMx; else if (hk0 < -hrMx) hk0 = -hrMx; secD = Nst * hk0; }
      for (k = 0; k < K; k++) {
        var rr = d[k] / dm, usc2 = uscm2;
        if (hideEg) { if (rr >= 0.4) { var lg = LN19 / Math.log(19 * rr); usc2 = uscm2 * lg * lg * rr; } else usc2 = 0.85 * uscm2; }
        if (us2 > usc2 && A[ck + k] > 0) {
          var usc = Math.sqrt(usc2);
          var qb = cap * A[ck + k] * (1 - usc2 / us2) * (1 - usc / uss);
          qbx[ck + k] = qb * (uu - secD * vv) / V; qby[ck + k] = qb * (vv + secD * uu) / V;
          var epsg = (usc / uss) * invMu; if (epsg < epsMin) epsg = epsMin;
          gk[ck + k] = qb * epsg;
          anyq = true;
        } else { qbx[ck + k] = 0; qby[ck + k] = 0; gk[ck + k] = 0; }
      }
      if (anyq) act[c] = 1;
    }

    // ---- 面の流砂量（上流化＋勾配による重力項）
    var aL, aR, S, fa, gg, fk;
    for (j = 0; j < ny; j++) {
      row = j * nx; frow = j * (nx + 1);
      f = frow * K; for (k = 0; k < K; k++) SFx[f + k] = 0;
      for (i = 1; i < nx; i++) {
        L = row + i - 1; R = L + 1; fk = (frow + i) * K;
        aL = act[L]; aR = act[R];
        if ((!aL && !aR) || fxC[L] || fxC[R]) { for (k = 0; k < K; k++) SFx[fk + k] = 0; continue; }
        S = (z[R] - z[L]) / dx; if (S > tanPhi) S = tanPhi; else if (S < -tanPhi) S = -tanPhi;
        // 固定床（堰堤など）への段上がりは掃流砂を通さない
        var blockLR = (z[R] - zfix[R] < 1e-6) && (z[R] - z[L] > 0.5 * h[L]);
        var blockRL = (z[L] - zfix[L] < 1e-6) && (z[L] - z[R] > 0.5 * h[R]);
        for (k = 0; k < K; k++) {
          fa = 0; gg = 0;
          if (aL) { if (qbx[L * K + k] > 0 && !blockLR) fa += qbx[L * K + k]; gg += gk[L * K + k]; }
          if (aR) { if (qbx[R * K + k] < 0 && !blockRL) fa += qbx[R * K + k]; gg += gk[R * K + k]; }
          SFx[fk + k] = fa - 0.5 * gg * S;
        }
      }
      // 下流端：流出
      L = row + nx - 1; fk = (frow + nx) * K;
      for (k = 0; k < K; k++) SFx[fk + k] = (act[L] && qbx[L * K + k] > 0) ? qbx[L * K + k] : 0;
    }
    for (i = 0; i < nx; i++) { for (k = 0; k < K; k++) { SFy[i * K + k] = 0; SFy[(ny * nx + i) * K + k] = 0; } }
    for (j = 1; j < ny; j++) {
      row = j * nx;
      for (i = 0; i < nx; i++) {
        R = row + i; L = R - nx; fk = R * K;
        aL = act[L]; aR = act[R];
        if ((!aL && !aR) || fxC[L] || fxC[R]) { for (k = 0; k < K; k++) SFy[fk + k] = 0; continue; }
        S = (z[R] - z[L]) / dx; if (S > tanPhi) S = tanPhi; else if (S < -tanPhi) S = -tanPhi;
        var bLR = (z[R] - zfix[R] < 1e-6) && (z[R] - z[L] > 0.5 * h[L]);
        var bRL = (z[L] - zfix[L] < 1e-6) && (z[L] - z[R] > 0.5 * h[R]);
        for (k = 0; k < K; k++) {
          fa = 0; gg = 0;
          if (aL) { if (qby[L * K + k] > 0 && !bLR) fa += qby[L * K + k]; gg += gk[L * K + k]; }
          if (aR) { if (qby[R * K + k] < 0 && !bRL) fa += qby[R * K + k]; gg += gk[R * K + k]; }
          SFy[fk + k] = fa - 0.5 * gg * S;
        }
      }
    }

    // ---- 粒径別の持ち出し制限（交換層内の存在量を超えない）
    var cf = dtm / ((1 - P.porosity) * dx);
    var scS = this.scS, any = false, out, fwk, fek, gsk, gnk;
    for (j = 0; j < ny; j++) {
      row = j * nx; frow = j * (nx + 1);
      for (i = 0; i < nx; i++) {
        c = row + i; ck = c * K;
        fwk = (frow + i) * K; fek = fwk + K; gsk = ck; gnk = (c + nx) * K;
        for (k = 0; k < K; k++) {
          out = 0;
          if (SFx[fek + k] > 0) out += SFx[fek + k];
          if (SFx[fwk + k] < 0) out -= SFx[fwk + k];
          if (SFy[gnk + k] > 0) out += SFy[gnk + k];
          if (SFy[gsk + k] < 0) out -= SFy[gsk + k];
          out *= cf;
          if (out > A[ck + k]) { scS[ck + k] = out > 0 ? A[ck + k] / out : 0; any = true; } else scS[ck + k] = 1;
        }
      }
    }
    if (any) {
      this.nLimS++;
      var r;
      for (j = 0; j < ny; j++) {
        row = j * nx; frow = j * (nx + 1);
        for (i = 1; i <= nx; i++) {
          fk = (frow + i) * K;
          for (k = 0; k < K; k++) {
            r = SFx[fk + k] > 0 ? scS[(row + i - 1) * K + k] : (i < nx ? scS[(row + i) * K + k] : 1);
            if (r < 1) SFx[fk + k] *= r;
          }
        }
      }
      for (j = 1; j < ny; j++) {
        row = j * nx;
        for (i = 0; i < nx; i++) {
          fk = (row + i) * K;
          for (k = 0; k < K; k++) {
            r = SFy[fk + k] > 0 ? scS[(row + i - nx) * K + k] : scS[(row + i) * K + k];
            if (r < 1) SFy[fk + k] *= r;
          }
        }
      }
    }

    // ---- 上流端の給砂
    var Qs = this.phase === 'run' ? this.QsAt(this.t) : 0;
    this.QsIn = Qs;
    var ns = this.ns, nyb = this.nyb, ff = this.ff, qin = this.qin, Qw = this.Qin;
    var feedCoef = (Qs > 0 && Qw > 0) ? Qs / Qw * dtm / ((1 - P.porosity) * dx) : 0;  // × qin[j] で層厚 [m]
    this.sedIn += Qs * dtm;

    // ---- Exner 式による更新と交換層・堆積層の入れ替え
    var qsOut = 0, dz, dAk;
    for (j = 0; j < ny; j++) {
      row = j * nx; frow = j * (nx + 1);
      fk = (frow + nx) * K;
      for (k = 0; k < K; k++) { qsOut += SFx[fk + k]; this.sedOutK[k] += SFx[fk + k] * dx * dtm; }
      for (i = 0; i < nx; i++) {
        c = row + i; ck = c * K;
        fwk = (frow + i) * K; fek = fwk + K; gsk = ck; gnk = (c + nx) * K;
        dz = 0;
        var feed = (i === 0 && feedCoef > 0 && qin[j] > 0) ? feedCoef * qin[j] : 0;
        for (k = 0; k < K; k++) {
          dAk = -cf * (SFx[fek + k] - SFx[fwk + k] + SFy[gnk + k] - SFy[gsk + k]);
          if (feed > 0) dAk += feed * ff[k];
          if (dAk !== 0) { var an = A[ck + k] + dAk; if (an < 0) { dAk = -A[ck + k]; an = 0; } A[ck + k] = an; dz += dAk; }
        }
        if (dz !== 0) { z[c] += dz; this.dzs[c] += dz; this.fixLayers(c); }
      }
    }
    this.QsOut = qsOut * dx;
    this.sedOut += this.QsOut * dtm;
  };

  // 安息角を超えた斜面の崩落。水に浸かった斜面は水中安息角、乾いた斜面は陸上安息角で判定する
  // （高橋・里深・奥村 1996 の 30° 規則を、XBeach・Delft3D の wet/dry slope の考え方で拡張）
  Sim.prototype.collapse = function () {
    var nx = this.nx, ny = this.ny, K = this.K, dx = this.dx, A = this.A, z = this.z, h = this.h, fxC = this.fixedC, BK = this.P.bank || {};
    var limW = Math.tan((BK.reposeWet || this.P.reposeDeg) * Math.PI / 180) * dx, limD = Math.tan((BK.reposeDry || this.P.reposeDeg) * Math.PI / 180) * dx;
    var HW = 0.01, i, j, c, k, row, lim, diff, hi, lo, av, mv, hk, lk, frc, hase = BK.collapse === 'hasegawa', kind = this.kind;
    for (var pass = 0; pass < 2; pass++) {
      var off = pass === 0 ? 1 : nx;
      var iEnd = pass === 0 ? nx - 1 : nx, jEnd = pass === 0 ? ny : ny - 1;
      for (j = 0; j < jEnd; j++) {
        row = j * nx;
        for (i = 0; i < iEnd; i++) {
          c = row + i;
          if (fxC[c] || fxC[c + off]) continue;
          lim = (h[c] > HW || h[c + off] > HW) ? limW : limD;
          diff = z[c] - z[c + off];
          if (diff > lim) { hi = c; lo = c + off; } else if (-diff > lim) { hi = c + off; lo = c; diff = -diff; } else continue;
          // 長谷川モデルでは、水際の乾いた側岸（高い側が乾き、低い側が濡れている）は崩落モデルに任せる
          if (hase && h[hi] <= HW && h[lo] > HW && kind[hi] === 0 && kind[lo] === 0) continue;
          hk = hi * K; lk = lo * K; av = 0;
          for (k = 0; k < K; k++) av += A[hk + k];
          if (av < 1e-6) continue;
          mv = 0.25 * (diff - lim); if (mv > 0.5 * av) mv = 0.5 * av;
          if (mv < 1e-6) continue;
          frc = mv / av;
          for (k = 0; k < K; k++) { var t = A[hk + k] * frc; A[hk + k] -= t; A[lk + k] += t; }
          z[hi] -= mv; z[lo] += mv;
          this.fixLayers(hi); this.fixLayers(lo);
        }
      }
    }
  };

  // ---- 2次流（長田・細田・村本・Rahman 1996）------------------------------------------
  // 流線の曲率 κ = 1/r = (u(u v_x + v v_y) − v(u u_x + v u_y)) / V³（左に曲がるとき正）。
  // 流速は前のステップの面の流束の平均。格子規模のばらつきを抑えるため、濡れたセルだけで 3×3 平均を 2 回とる。
  // 弱い曲率を前提とする式なので、|h κ| を hrMax で頭打ちにする。
  Sim.prototype.curvature = function () {
    var nx = this.nx, ny = this.ny, N = this.N, dx = this.dx, h = this.h, FxM = this.FxM, GyM = this.GyM;
    var kap = this.kap, tmp = this.kapW, uf = this.uf, vf = this.vf, HS = 0.05, i, j, c, f;
    var hrMax = (this.P.sec && this.P.sec.hrMax) || 0.05;
    for (c = 0; c < N; c++) {
      var hh = h[c];
      if (hh > HS) { i = c % nx; j = (c - i) / nx; f = j * (nx + 1) + i; uf[c] = 0.5 * (FxM[f] + FxM[f + 1]) / hh; vf[c] = 0.5 * (GyM[c] + GyM[c + nx]) / hh; }
      else { uf[c] = 0; vf[c] = 0; }
      tmp[c] = 0;
    }
    for (j = 1; j < ny - 1; j++) for (i = 1; i < nx - 1; i++) {
      c = j * nx + i;
      if (!(h[c] > HS && h[c - 1] > HS && h[c + 1] > HS && h[c - nx] > HS && h[c + nx] > HS)) continue;
      var uu = uf[c], vv = vf[c], V2 = uu * uu + vv * vv;
      if (V2 < 0.01) continue;
      var ux = (uf[c + 1] - uf[c - 1]) / (2 * dx), uy = (uf[c + nx] - uf[c - nx]) / (2 * dx);
      var vx = (vf[c + 1] - vf[c - 1]) / (2 * dx), vy = (vf[c + nx] - vf[c - nx]) / (2 * dx);
      tmp[c] = (uu * (uu * vx + vv * vy) - vv * (uu * ux + vv * uy)) / (V2 * Math.sqrt(V2));
    }
    for (var pass = 0; pass < 2; pass++) {
      var src = pass === 0 ? tmp : kap, dst = pass === 0 ? kap : tmp;
      for (j = 0; j < ny; j++) for (i = 0; i < nx; i++) {
        c = j * nx + i; dst[c] = 0;
        if (!(h[c] > HS)) continue;
        var sm = 0, nn = 0;
        for (var dj = -1; dj <= 1; dj++) { var jj = j + dj; if (jj < 0 || jj >= ny) continue;
          for (var di = -1; di <= 1; di++) { var ii = i + di; if (ii < 0 || ii >= nx) continue; var cc = jj * nx + ii; if (h[cc] > HS) { sm += src[cc]; nn++; } } }
        dst[c] = sm / nn;
      }
    }
    for (c = 0; c < N; c++) {
      var lim = hrMax / Math.max(h[c], 1e-3), kk = tmp[c];
      kap[c] = kk > lim ? lim : (kk < -lim ? -lim : kk);
    }
  };
  // 2次流による運動量の分散項（Kalkwijk・de Vriend 1980、長田ら 1996 式19〜23）。
  // s–n 座標の T_sn = h ū_s ū_n の水深積分 = −k_sn h² V² κ（外岸向きに主流の運動量を運ぶ）を x–y に回す：
  //   Txx = −Tyy = 2 k_sn h² κ u v、Txy = −k_sn h² κ (u² − v²)
  //   k_sn = 5a − 15.6a² + 37.5a³、a = √g/(κ_K C)、C = h^(1/6)/n。外挿を避けるため a ≤ 0.4
  Sim.prototype.secStress = function () {
    var N = this.N, h = this.h, uf = this.uf, vf = this.vf, kap = this.kap, Txx = this.Txx, Txy = this.Txy, msk = this.secMask;
    var nman = this.P.manning, HS = 0.05, sg = Math.sqrt(G);
    for (var c = 0; c < N; c++) {
      Txx[c] = 0; Txy[c] = 0; msk[c] = 0;
      var hh = h[c]; if (!(hh > HS)) continue;
      msk[c] = 1;
      var uu = uf[c], vv = vf[c], kk = kap[c];
      if (kk === 0) continue;
      var a = sg * nman / (0.4 * Math.pow(hh, 1 / 6)); if (a > 0.4) a = 0.4;
      var ksn = 5 * a - 15.6 * a * a + 37.5 * a * a * a;
      var w = ksn * hh * hh * kk;
      Txx[c] = 2 * w * uu * vv; Txy[c] = -w * (uu * uu - vv * vv);
    }
  };

  // ---- 長谷川の側岸崩落モデル（長谷川 1983、長田ら 1996 式17）-----------------------------
  // 乾いた側岸セル b の隣が水に浸かっているとき、b の位置の水面から水中安息角 βk で下る線を引く。
  // 線より上の土砂（崩落する土塊 A）と、線より下の洗掘の不足分 B を、横断方向（4 方向）に数える。
  // B ≥ A になったら土塊 A が一度に崩れ、B の部分を埋め戻す（崩落限界幅 ΔBc は 1 セル）。
  Sim.prototype.hasegawaCollapse = function () {
    var BK = this.P.bank || {}; if (BK.collapse !== 'hasegawa') return;
    var nx = this.nx, ny = this.ny, N = this.N, dx = this.dx, h = this.h, z = this.z, kind = this.kind, fxC = this.fixedC;
    var tb = Math.tan((BK.reposeWet || this.P.reposeDeg) * Math.PI / 180) * dx;
    var HW = 0.01, MAXM = 40, exC = this._exC || (this._exC = new Int32Array(MAXM + 1)), exV = this._exV || (this._exV = new Float64Array(MAXM + 1));
    var dfC = this._dfC || (this._dfC = new Int32Array(MAXM + 1)), dfV = this._dfV || (this._dfV = new Float64Array(MAXM + 1));
    var DI = [1, -1, 0, 0], DJ = [0, 0, 1, -1];
    for (var c = 0; c < N; c++) {
      if (h[c] >= HW || kind[c] !== 0 || fxC[c]) continue;
      var i0 = c % nx, j0 = (c - i0) / nx;
      for (var q = 0; q < 4; q++) {
        var i1 = i0 + DI[q], j1 = j0 + DJ[q];
        if (i1 < 0 || i1 >= nx || j1 < 0 || j1 >= ny) continue;
        var w1 = j1 * nx + i1;
        if (h[w1] < HW || kind[w1] !== 0 || fxC[w1]) continue;
        var H = z[w1] + h[w1]; if (z[c] <= H) continue;
        var nE = 1, nD = 0, Aex = z[c] - H, Bdf = 0, stage = 0;
        exC[0] = c; exV[0] = Aex;
        for (var m = 1; m <= MAXM; m++) {
          var ii = i0 + m * DI[q], jj = j0 + m * DJ[q];
          if (ii < 0 || ii >= nx || jj < 0 || jj >= ny) break;
          var cc = jj * nx + ii;
          if (h[cc] < HW || kind[cc] !== 0 || fxC[cc]) break;
          var dzl = z[cc] - (H - m * tb);
          if (dzl > 0) { if (stage === 0) { exC[nE] = cc; exV[nE] = dzl; nE++; Aex += dzl; } else break; }
          else if (dzl < 0) { stage = 1; dfC[nD] = cc; dfV[nD] = -dzl; nD++; Bdf -= dzl; }
          else if (stage === 1) break;
        }
        if (nD === 0 || Bdf < Aex) continue;
        for (var e = 0; e < nE; e++) for (var fq = 0; fq < nD; fq++) this.moveSed(exC[e], dfC[fq], exV[e] * dfV[fq] / Bdf);
        this.nHase++;
        break;
      }
    }
  };
  // セル a の上から厚さ dz の土砂を（堆積層の粒度で）セル b へ移す。動かせた厚さを返す
  Sim.prototype.moveSed = function (a, b, dz) {
    var K = this.K, A = this.A, ak = a * K, bk = b * K, left = dz, k, it = 0;
    while (left > 1e-9 && it++ < 400) {
      var av = 0; for (k = 0; k < K; k++) av += A[ak + k];
      if (av < 1e-9) break;
      var take = left < av ? left : av, frc = take / av;
      for (k = 0; k < K; k++) { var t = A[ak + k] * frc; A[ak + k] -= t; A[bk + k] += t; }
      this.z[a] -= take; this.z[b] += take; left -= take;
      this.fixLayers(a); this.fixLayers(b);
    }
    return dz - left;
  };

  // 水際の側岸侵食：濡れたセルが侵食されたとき、その一部（割合 dryEro）を隣の乾いた側岸セルから削る
  // （Delft3D の dry-cell erosion と同じ考え方。流砂が乾いた側岸を直接削れない格子モデルの弱点を補う）
  Sim.prototype.dryBankErosion = function () {
    var th = (this.P.bank && this.P.bank.dryEro) || 0; if (!(th > 0)) return;
    var nx = this.nx, ny = this.ny, N = this.N, K = this.K, A = this.A, z = this.z, h = this.h, kind = this.kind, fxC = this.fixedC, dzs = this.dzs;
    var HW = 0.01, nb = [0, 0, 0, 0], c, k, q, n;
    for (c = 0; c < N; c++) {
      if (dzs[c] >= 0 || h[c] < HW || kind[c] !== 0 || fxC[c]) continue;
      var i = c % nx, cnt = 0;
      var cand = [i > 0 ? c - 1 : -1, i < nx - 1 ? c + 1 : -1, c - nx, c + nx];
      for (q = 0; q < 4; q++) {
        n = cand[q]; if (n < 0 || n >= N) continue;
        if (h[n] >= HW || kind[n] !== 0 || fxC[n] || z[n] <= z[c]) continue;
        nb[cnt++] = n;
      }
      if (!cnt) continue;
      var want = th * (-dzs[c]) / cnt;
      for (q = 0; q < cnt; q++) {
        n = nb[q]; var nk = n * K, ck = c * K, av = 0;
        for (k = 0; k < K; k++) av += A[nk + k];
        if (av < 1e-6) continue;
        var mv = Math.min(want, 0.5 * av, 0.5 * (z[n] - z[c]));
        if (mv < 1e-9) continue;
        var frc = mv / av;
        for (k = 0; k < K; k++) { var t = A[nk + k] * frc; A[nk + k] -= t; A[ck + k] += t; }
        z[n] -= mv; z[c] += mv;
        this.fixLayers(n); this.fixLayers(c);
      }
    }
  };

  // 侵食・堆積速度式による非平衡の河床変動（高橋・里深 1999、高橋・中川・里深・奥村・安本 1998）
  //   i_k = f_bk δe (C∞−C)/C* |u*−u*ck|  （C∞ > C かつ u* > u*ck、侵食）
  //   i_k = f_0k δd (C∞−C)/C* |u*−u*ck|  （C∞ ≤ C かつ u* < u*ck、堆積）
  //   C∞ = Σ f_bk q_bu,k / q、q_bu,k は粒径 k だけの河床に対する芦田・道上式
  //   ∂(C_k h)/∂t + ∇·(C_k h u') = C* i_k、∂z/∂t = −Σ i_k
  //   u' は長谷川（1983）の横断方向流砂で向きを β2 = arctan(−K ∂z/∂n) だけ回した流速
  Sim.prototype.morphNoneq = function (dtm) {
    var nsub = Math.max(1, Math.ceil(dtm / Math.max(this.dt, 1e-6) - 1e-9)), dts = dtm / nsub;
    for (var s_ = 0; s_ < nsub; s_++) this.noneqStep(dts);
  };
  Sim.prototype.noneqStep = function (dts) {
    var nx = this.nx, ny = this.ny, N = this.N, dx = this.dx, K = this.K, P = this.P, NE = P.noneq || {};
    var h = this.h, z = this.z, zfix = this.zfix, A = this.A, d = this.d, Cc = this.Cc, kind = this.kind, fxC = this.fixedC;
    var FxM = this.FxM, GyM = this.GyM, SFx = this.SFx, SFy = this.SFy, tu = this.qbx, tv = this.qby;   // 作業配列を流用
    var g = G, sg1 = this.sg1, Cst = 1 - P.porosity, HSED = 0.02, tcs = this.tcs, hideEg = P.hiding === 'egiazaroff';
    var de0 = NE.de > 0 ? NE.de : 0.6, dd0 = NE.dd > 0 ? NE.dd : 0.6, hE = NE.hEdge > 0 ? NE.hEdge : 2.5 * this.dm0;
    var fe = NE.fe > 0 ? NE.fe : 1, fd = NE.fd > 0 ? NE.fd : 1, gK = (P.gravK == null ? 1 : P.gravK) / Math.sqrt(1.0 * 0.45);
    var c, k, ck, i, j, row, frow;
    if (!this.dord) { this.dord = []; for (k = 0; k < K; k++) this.dord.push(k); var dd_ = d; this.dord.sort(function (a, b) { return dd_[a] - dd_[b]; }); }
    var ord = this.dord, sqd3 = this.sqd3 || (this.sqd3 = Array.prototype.map.call(d, function (dk) { return Math.sqrt(sg1 * g * dk * dk * dk); }));
    var iK = this.dA;
    var SEC = P.sec || {}, secBed = !!SEC.bed, kapA = this.kap, Nst = SEC.Nstar || 11.5, hrMx = SEC.hrMax || 0.05;
    // ---- 1a. 輸送速度と限界摩擦速度（河床を変える前の地形で全セル分を先に求める）
    var ucm = this.neUcm || (this.neUcm = new Float64Array(N));
    for (c = 0; c < N; c++) {
      ck = c * K; tu[ck] = 0; tv[ck] = 0; ucm[c] = 0;
      var h0 = h[c]; if (h0 < HSED) continue;
      var ic0 = c % nx, jc0 = (c - ic0) / nx, fw0 = jc0 * (nx + 1) + ic0;
      var u0 = 0.5 * (FxM[fw0] + FxM[fw0 + 1]) / h0, v0 = 0.5 * (GyM[c] + GyM[c + nx]) / h0, V0 = Math.sqrt(u0 * u0 + v0 * v0);
      var As0 = 0, dms0 = 0; for (k = 0; k < K; k++) { As0 += A[ck + k]; dms0 += A[ck + k] * d[k]; }
      var dm0 = As0 > 1e-9 ? dms0 / As0 : this.dm0;
      var ks0 = d[ord[K - 1]];
      if (As0 > 1e-9) { var cum0 = 0; for (var o0 = 0; o0 < K; o0++) { cum0 += A[ck + ord[o0]] / As0; if (cum0 >= 0.9 - 1e-9) { ks0 = d[ord[o0]]; break; } } }
      var den0 = 6.0 + 5.75 * Math.log10(Math.max(h0 / ks0, 0.2)); if (den0 < 2) den0 = 2;
      var us0 = V0 / den0, us20 = us0 * us0, uc0 = this.critUs2(c, dm0, u0, v0);
      ucm[c] = uc0; tcs[c] = uc0 / (sg1 * g * dm0);
      if (V0 > 1e-6) {
        var b1 = Math.atan2(v0, u0), nxn = -Math.sin(b1), nyn = Math.cos(b1);
        var zl = ic0 > 0 ? z[c - 1] : z[c], zr = ic0 < nx - 1 ? z[c + 1] : z[c], zd = jc0 > 0 ? z[c - nx] : z[c], zu = jc0 < ny - 1 ? z[c + nx] : z[c];
        var dzdn = ((zr - zl) / ((ic0 > 0 && ic0 < nx - 1) ? 2 * dx : dx)) * nxn + ((zu - zd) / ((jc0 > 0 && jc0 < ny - 1) ? 2 * dx : dx)) * nyn;
        var Kc = us20 > 1e-12 ? Math.sqrt(Math.min(1, uc0 / us20)) * gK : 0;
        var sd0 = 0;   // 2次流による河床付近の流向の偏り（長田ら 1996 式11・14）
        if (secBed) { var hk1 = h0 * kapA[c]; if (hk1 > hrMx) hk1 = hrMx; else if (hk1 < -hrMx) hk1 = -hrMx; sd0 = Nst * hk1; }
        var b2 = Math.atan(sd0 - Kc * dzdn); if (b2 > 1) b2 = 1; else if (b2 < -1) b2 = -1;
        tu[ck] = V0 * Math.cos(b1 + b2); tv[ck] = V0 * Math.sin(b1 + b2);
      }
    }
    // ---- 1b. 各セルの侵食・堆積
    for (c = 0; c < N; c++) {
      ck = c * K; var hh = h[c], dzc = 0;
      if (hh < HSED) {             // 水が引いたら、流れの中の土砂はその場に堆積させる
        var any = 0; for (k = 0; k < K; k++) if (Cc[ck + k] > 0) { var dep = Cc[ck + k] / Cst; A[ck + k] += dep; dzc += dep; Cc[ck + k] = 0; any = 1; }
        if (any) { z[c] += dzc; this.dzs[c] += dzc; this.fixLayers(c); }
        continue;
      }
      if (fxC[c]) continue;        // 巨礫のセルでは侵食も堆積もしない（土砂は流れて通る）
      var ic = c % nx, jc = (c - ic) / nx, fwq = jc * (nx + 1) + ic;
      var uu = 0.5 * (FxM[fwq] + FxM[fwq + 1]) / hh, vv = 0.5 * (GyM[c] + GyM[c + nx]) / hh;
      var V = Math.sqrt(uu * uu + vv * vv);
      var As = 0, dms = 0; for (k = 0; k < K; k++) { As += A[ck + k]; dms += A[ck + k] * d[k]; }
      var dm = As > 1e-9 ? dms / As : this.dm0;
      var ks = d[ord[K - 1]];
      if (As > 1e-9) { var cum = 0; for (var o = 0; o < K; o++) { cum += A[ck + ord[o]] / As; if (cum >= 0.9 - 1e-9) { ks = d[ord[o]]; break; } } }
      var den = 6.0 + 5.75 * Math.log10(Math.max(hh / ks, 0.2)); if (den < 2) den = 2;
      var us = V / den, us2 = us * us, uscm2 = ucm[c];
      var q = hh * V, Ct = 0; for (k = 0; k < K; k++) Ct += Cc[ck + k]; Ct /= hh;
      var Cinf = 0;
      for (k = 0; k < K; k++) {
        var rr = d[k] / dm, usc2 = uscm2;
        if (hideEg) { if (rr >= 0.4) { var lg = LN19 / Math.log(19 * rr); usc2 = uscm2 * lg * lg * rr; } else usc2 = 0.85 * uscm2; }
        iK[k] = usc2;           // 一時的に u*ck² を入れておく
        if (As > 1e-9 && q > 1e-9 && us2 > usc2) {
          var tk = us2 / (sg1 * g * d[k]), tck = usc2 / (sg1 * g * d[k]);
          var qbu = 17 * sqd3[k] * Math.pow(tk, 1.5) * (1 - tck / tk) * (1 - Math.sqrt(usc2) / us);
          Cinf += (A[ck + k] / As) * qbu / q;
        }
      }
      var de = hh < hE ? de0 * fe : de0, ddd = hh < hE ? dd0 * fd : dd0, diffC = Cinf - Ct;
      for (k = 0; k < K; k++) {
        var usck = Math.sqrt(iK[k]), rate = 0;
        if (diffC > 0) { if (us > usck && As > 1e-9) rate = (A[ck + k] / As) * de * diffC / Cst * (us - usck); }
        else if (Ct > 0 && us < usck) rate = (Cc[ck + k] / (Ct * hh)) * ddd * diffC / Cst * (usck - us);
        var dzk = -rate * dts;                       // 河床の変化（空隙込みの層厚）
        if (dzk < 0 && -dzk > A[ck + k]) dzk = -A[ck + k];               // 交換層にある量まで
        if (dzk > 0 && dzk * Cst > Cc[ck + k]) dzk = Cc[ck + k] / Cst;   // 流れの中にある量まで
        if (dzk !== 0) { A[ck + k] += dzk; Cc[ck + k] -= dzk * Cst; dzc += dzk; }
      }
      if (dzc !== 0) { z[c] += dzc; this.dzs[c] += dzc; this.fixLayers(c); }
    }
    // ---- 2. 流れの中の土砂の移流（風上差分）
    var Qs = this.phase === 'run' ? this.QsAt(this.t) : 0, Qw = this.Qin, qin = this.qin, ff = this.ff;
    this.QsIn = Qs; this.sedIn += Qs * dts;
    var Cin = (Qs > 0 && Qw > 0) ? Qs / Qw : 0, cf = dts / dx, L, R, fk, uf, wl, wr, dn;
    for (j = 0; j < ny; j++) {
      row = j * nx; frow = j * (nx + 1);
      for (k = 0; k < K; k++) SFx[frow * K + k] = 0;
      for (i = 1; i < nx; i++) {
        L = row + i - 1; R = L + 1; fk = (frow + i) * K;
        wl = h[L] >= HSED; wr = h[R] >= HSED;
        uf = wl && wr ? 0.5 * (tu[L * K] + tu[R * K]) : wl ? Math.max(0, tu[L * K]) : wr ? Math.min(0, tu[R * K]) : 0;
        dn = uf > 0 ? L : R;
        for (k = 0; k < K; k++) SFx[fk + k] = uf * Cc[dn * K + k];
      }
      L = row + nx - 1; fk = (frow + nx) * K; uf = h[L] >= HSED ? Math.max(0, tu[L * K]) : 0;
      for (k = 0; k < K; k++) SFx[fk + k] = uf * Cc[L * K + k];
    }
    for (i = 0; i < nx; i++) for (k = 0; k < K; k++) { SFy[i * K + k] = 0; SFy[(ny * nx + i) * K + k] = 0; }
    for (j = 1; j < ny; j++) {
      row = j * nx;
      for (i = 0; i < nx; i++) {
        R = row + i; L = R - nx; fk = R * K;
        wl = h[L] >= HSED; wr = h[R] >= HSED;
        uf = wl && wr ? 0.5 * (tv[L * K] + tv[R * K]) : wl ? Math.max(0, tv[L * K]) : wr ? Math.min(0, tv[R * K]) : 0;
        dn = uf > 0 ? L : R;
        for (k = 0; k < K; k++) SFy[fk + k] = uf * Cc[dn * K + k];
      }
    }
    // 持ち出し制限（セルにある量を超えて流さない）
    var scS = this.scS, lim = false;
    for (j = 0; j < ny; j++) {
      row = j * nx; frow = j * (nx + 1);
      for (i = 0; i < nx; i++) {
        c = row + i; ck = c * K; var fw = (frow + i) * K, fe_ = fw + K, gs = ck, gn = (c + nx) * K;
        for (k = 0; k < K; k++) {
          var out = (SFx[fe_ + k] > 0 ? SFx[fe_ + k] : 0) - (SFx[fw + k] < 0 ? SFx[fw + k] : 0) + (SFy[gn + k] > 0 ? SFy[gn + k] : 0) - (SFy[gs + k] < 0 ? SFy[gs + k] : 0);
          out *= cf;
          if (out > Cc[ck + k] + 1e-15) { scS[ck + k] = out > 0 ? Cc[ck + k] / out : 0; lim = true; } else scS[ck + k] = 1;
        }
      }
    }
    if (lim) {
      for (j = 0; j < ny; j++) { row = j * nx; frow = j * (nx + 1);
        for (i = 1; i <= nx; i++) { fk = (frow + i) * K; for (k = 0; k < K; k++) { var r = SFx[fk + k] > 0 ? scS[(row + i - 1) * K + k] : (i < nx ? scS[(row + i) * K + k] : 1); if (r < 1) SFx[fk + k] *= r; } } }
      for (j = 1; j < ny; j++) { row = j * nx;
        for (i = 0; i < nx; i++) { fk = (row + i) * K; for (k = 0; k < K; k++) { var r2 = SFy[fk + k] > 0 ? scS[(row + i - nx) * K + k] : scS[(row + i) * K + k]; if (r2 < 1) SFy[fk + k] *= r2; } } }
    }
    var qsOut = 0;
    for (j = 0; j < ny; j++) {
      row = j * nx; frow = j * (nx + 1);
      fk = (frow + nx) * K;
      for (k = 0; k < K; k++) { qsOut += SFx[fk + k]; this.sedOutK[k] += SFx[fk + k] * dx * dts; }
      for (i = 0; i < nx; i++) {
        c = row + i; ck = c * K; var fw2 = (frow + i) * K, fe2 = fw2 + K, gn2 = (c + nx) * K;
        var src = (i === 0 && Cin > 0 && qin[j] > 0) ? cf * Cin * qin[j] : 0;
        for (k = 0; k < K; k++) {
          var nv = Cc[ck + k] - cf * (SFx[fe2 + k] - SFx[fw2 + k] + SFy[gn2 + k] - SFy[ck + k]) + src * ff[k];
          Cc[ck + k] = nv > 0 ? nv : 0;
        }
      }
    }
    this.QsOut = qsOut * dx;
    this.sedOut += this.QsOut * dts;
  };
  // 流れの中にある土砂の量（実質体積） [m3]
  Sim.prototype.sedimentInWater = function () {
    var s = 0, Cc = this.Cc; for (var c = 0; c < Cc.length; c++) s += Cc[c];
    return s * this.dx * this.dx;
  };

  // 交換層厚を La に保つよう、堆積層との間で土砂を出し入れする
  Sim.prototype.fixLayers = function (c) {
    var K = this.K, A = this.A, La = this.La, Ls = this.Ls, NL = this.NL, SF = this.SF;
    var ck = c * K, k, As = 0;
    for (k = 0; k < K; k++) As += A[ck + k];
    var nl = this.nl[c], th = this.topTh[c], base;
    if (As > La + 1e-12) {
      // 堆積：余剰分を交換層の粒度のまま堆積層へ
      var e = As - La, fr = e / As;
      if (nl === 0) { nl = 1; th = 0; }
      base = (c * NL + nl - 1) * K;
      var tn = th + e;
      for (k = 0; k < K; k++) {
        var tk = A[ck + k] * fr;
        A[ck + k] -= tk;
        SF[base + k] = (SF[base + k] * th + tk) / tn;
      }
      th = tn;
      if (th > Ls && nl < NL) {
        var nb = base + K;
        for (k = 0; k < K; k++) SF[nb + k] = SF[base + k];
        nl++; th -= Ls;
      }
      this.nl[c] = nl; this.topTh[c] = th;
    } else if (As < La - 1e-12 && nl > 0) {
      // 侵食：不足分を堆積層最上層の粒度で補う
      var def = La - As;
      while (def > 1e-12 && nl > 0) {
        base = (c * NL + nl - 1) * K;
        var take = def < th ? def : th;
        for (k = 0; k < K; k++) A[ck + k] += take * SF[base + k];
        th -= take; def -= take;
        if (th <= 1e-9) { nl--; th = nl > 0 ? Ls : 0; }
      }
      this.nl[c] = nl; this.topTh[c] = th;
    }
  };

  // 流れの向きの河床勾配 i（下り勾配を正、0〜0.2 に制限）。堰堤のセルは勾配の計算から除く
  Sim.prototype.flowSlope = function (c, ux, uy) {
    var nx = this.nx, z = this.z, kind = this.kind, dx = this.dx, i = c % nx, a, b, gx, gy;
    a = (i > 0 && kind[c - 1] !== 2) ? c - 1 : -1; b = (i < nx - 1 && kind[c + 1] !== 2) ? c + 1 : -1;
    gx = (a >= 0 && b >= 0) ? (z[b] - z[a]) / (2 * dx) : a >= 0 ? (z[c] - z[a]) / dx : b >= 0 ? (z[b] - z[c]) / dx : -this.I;
    a = (c - nx >= 0 && kind[c - nx] !== 2) ? c - nx : -1; b = (c + nx < this.N && kind[c + nx] !== 2) ? c + nx : -1;
    gy = (a >= 0 && b >= 0) ? (z[b] - z[a]) / (2 * dx) : a >= 0 ? (z[c] - z[a]) / dx : b >= 0 ? (z[b] - z[c]) / dx : 0;
    var V = Math.sqrt(ux * ux + uy * uy); if (V < 1e-9) return 0;
    var sl = -(gx * ux + gy * uy) / V;
    return sl < 0 ? 0 : sl > 0.2 ? 0.2 : sl;
  };
  // 平均粒径 dm に対する限界摩擦速度の2乗 u*cm² [m²/s²]
  Sim.prototype.critUs2 = function (c, dm, ux, uy) {
    var P = this.P, sg1 = this.sg1;
    if (P.tcMethod === 'iwagaki') return iwagakiUc2(dm, sg1);
    if (P.tcMethod === 'const') return P.tauCm * sg1 * G * dm;
    // 芦田・高橋・水山（1978）式(24)：τ*c = 0.04 × 10^(1.72 i)
    return 0.04 * Math.pow(10, 1.72 * this.flowSlope(c, ux, uy)) * sg1 * G * dm;
  };

  // ------------------------------------------------------------------ 集計
  Sim.prototype.waterVolume = function () {
    var s = 0, h = this.h; for (var c = 0; c < this.N; c++) s += h[c];
    return s * this.dx * this.dx;
  };
  // 計算区間内の堆積土砂量の増分（空隙込み） [m3]
  Sim.prototype.storage = function () {
    var s = 0, z = this.z, z0 = this.z0; for (var c = 0; c < this.N; c++) s += z[c] - z0[c];
    return s * this.dx * this.dx;
  };
  // 堰堤ごとの堆砂量（各堰堤の上流側区間の河床上昇・低下の合計、空隙込み） [m3]
  Sim.prototype.storageByDam = function () {
    var res = [], nx = this.nx, ny = this.ny, z = this.z, z0 = this.z0, a2 = this.dx * this.dx;
    for (var q = 0; q <= this.dams.length; q++) res.push(0);
    for (var c = 0; c < this.N; c++) res[this.region[c]] += z[c] - z0[c];
    for (q = 0; q < res.length; q++) res[q] *= a2;
    return res;   // 最後の要素は最下流堰堤より下流の区間
  };
  // 表層（交換層）の平均粒径 [m]
  Sim.prototype.surfaceDm = function (c) {
    var K = this.K, A = this.A, d = this.d, s = 0, m = 0, ck = c * K;
    for (var k = 0; k < K; k++) { s += A[ck + k]; m += A[ck + k] * d[k]; }
    return s > 1e-9 ? m / s : NaN;
  };
  // 縦断分布（横断方向に集約）
  Sim.prototype.profile = function () {
    var nx = this.nx, ns = this.ns, nyb = this.nyb, dx = this.dx, K = this.K;
    var z = this.z, z0 = this.z0, h = this.h, kind = this.kind;
    var x = new Float64Array(nx), zmin = new Float64Array(nx), zmean = new Float64Array(nx);
    var z0min = new Float64Array(nx), z0mean = new Float64Array(nx), ws = new Float64Array(nx);
    var dm = new Float64Array(nx), wet = new Float64Array(nx), qs = new Float64Array(nx), hmax = new Float64Array(nx);
    for (var i = 0; i < nx; i++) {
      x[i] = (i + 0.5) * dx;
      var mn = 1e30, m0 = 1e30, sm = 0, s0 = 0, cnt = 0, wsum = 0, wcnt = 0, dsum = 0, dcnt = 0, hm = 0, isDam = false;
      for (var j = 0; j < this.ny; j++) {
        var c = j * nx + i;
        if (!this.fl[c]) continue;
        if (kind[c] === 2) { isDam = true; continue; }
        if (z[c] < mn) mn = z[c]; if (z0[c] < m0) m0 = z0[c];
        sm += z[c]; s0 += z0[c]; cnt++;
        if (h[c] > 0.03) {
          wsum += (z[c] + h[c]) * h[c]; wcnt += h[c]; if (h[c] > hm) hm = h[c];
          var dd = this.surfaceDm(c); if (dd === dd) { dsum += dd; dcnt++; }
        }
      }
      if (cnt) { zmin[i] = mn; z0min[i] = m0; zmean[i] = sm / cnt; z0mean[i] = s0 / cnt; } else { zmin[i] = z0min[i] = zmean[i] = z0mean[i] = NaN; }
      ws[i] = wcnt > 0 ? wsum / wcnt : NaN;
      wet[i] = 0; for (j = 0; j < this.ny; j++) if (this.fl[j * nx + i] && h[j * nx + i] > 0.03) wet[i] += dx;
      dm[i] = dcnt > 0 ? dsum / dcnt : NaN; hmax[i] = hm;
      var q = 0;
      for (j = 0; j < this.ny; j++) { var fk = (j * (nx + 1) + i + 1) * K; for (var k = 0; k < K; k++) q += this.SFx[fk + k]; }
      qs[i] = q * dx;
    }
    return { x: x, zmin: zmin, zmean: zmean, z0min: z0min, z0mean: z0mean, ws: ws, dm: dm, wet: wet, qs: qs, hmax: hmax };
  };

  var api = { Sim: Sim, defaultParams: defaultParams, G: G };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.Sabo2D = api;
})(typeof window !== 'undefined' ? window : globalThis);
