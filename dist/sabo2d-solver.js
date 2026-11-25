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
      // 支川：{ x: 合流点 [m], side: 'L'|'R'（左岸・右岸）, B: 幅 [m], I: 勾配 1/n, len: 計算区間に入れる長さ [m],
      //        r: 本川流量に対する流量の比, rs: 本川の供給土砂量に対する比, lag: 本川より遅れる時間 [s] }
      tribs: [],
      // 伏流（河床への浸透と地下水の流れ。Darcy 則・Dupuit 近似）
      gw: { on: false, k: 0, kMax: 0.05, s0: 0 },   // k：透水係数 [m/s]（0 なら粒度から Kozeny–Carman 式）、kMax：上限、s0：乾いた堆積層の初期の飽和度
      // ---- DEM の流域（地形を DEM から切り出したとき）
      dem: { soilD: 1.0, chW: 6, ks: 1e-5, s0: 0.3, duration: 3600, outlet: 'sink', fill: true, pit: 0.5, edgeOut: true },   // outlet：出口の扱い（sink 自由流出 | channel 水路を延ばす）、edgeOut：流域界の外へ流れる水を抜く、fill：計算格子で窪地を埋める、pit：これより深い 1 セルの穴を埋める [m]（0 でしない）   // 斜面の土層の厚さ [m]、侵食できる谷筋の幅 [m]、土層の飽和透水係数 [m/s]、土層の初期の飽和度、計算時間 [s]
      sources: [],       // 点の流入 [{x, y, ux, uy, I, hydro: {上流端と同じ形}}]
      latIn: null,       // 雨の横流入 {dt, pts: [[x, y, ux, uy, I]], Q: [[流量の時系列]]}
      rain: { on: false, method: 'lumped', mode: 'uniform', table: [], frames: [], loss: { type: 'f', f: 0.8 } },
      slides: [],        // 斜面崩壊 [{x, y, R, D, type: 'shallow'|'deep', trigger: 'fs'|'time', t, move: 'auto'|'block'|'fluid', c, phi, gt, gs, ne, m0, Tk, Sr0, CL, mu, thDep, mc}]
      slideChk: 1,       // 崩壊の判定の間隔 [s]
      ldams: [],         // 天然ダム [{x, y, ux, uy, H, Wc, up, dn, lat}]
      ldam: { k: 0, fail: 'auto', slipLevel: 0.9, slipDeg: 0, thDep: 10, phi: 35, c: 0, gs: 20, Ts: 0.5, s0: 0 },
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
      hiding: 'egiazaroff', // 粒径別の限界掃流力：'egiazaroff' Egiazaroff 式の芦田・道上修正（選択的な輸送が起きる）| 'atm' 全粒径で平均粒径の値（芦田・高橋・水山。全粒径が同じ割合で動く）
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
      transport: 'eq',   // 'eq' 平衡流砂量 | 'noneq' 侵食・堆積速度式（掃流） | 'debris' 土石流（高橋の式）
      // ---- 土石流（高橋 2007、Kanako 2D と同じ形の式）
      debris: {
        de: 0.0007,      // 侵食速度係数 δe
        dd: 0.05,        // 堆積速度係数 δd
        tanPhi: 0.7,     // 内部摩擦角の正接 tanφ
        rho: 1000,       // 間隙流体の密度 ρ [kg/m3]
        Vmin: 0.1,       // 堆積速度に使う流速の下限 [m/s]（止まりかけの流れでも堆積させる。本プログラムの仮定）
        manningFloor: true   // 抵抗は Manning 式の値を下限にする（本プログラムの仮定）
      },   // 'eq' 平衡流砂量（芦田・道上式）| 'noneq' 侵食・堆積速度式（高橋・里深 1999、高橋ら 1998）
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
      snapEvery: 600,    // 時間ごとの結果を記録する間隔 [s]（出水の時刻）。0 で記録しない
      // ---- 透過型堰堤（スリット・鋼製）の共通設定
      openDam: {
        w: 0.5,          // 鋼製の部材の幅・複数スリットの間の柱の幅 [m]
        beta: 1.5,       // 閉塞の判定：開口幅 l ≤ beta × d95 で、流れの土砂濃度が Cblk 以上なら閉塞（土石流対策技術指針 2000 の 1.5、2016 の指針は 1.0）
        Cblk: 0.2,       // 閉塞が起きる流れの土砂濃度（容積）の下限（本プログラムの仮定）
        alpha: 1.0,      // 粒径 d が l/alpha 以上の礫は開口を通れない
        phiBlock: 0.05   // 閉塞後の開口率（礫の隙間を水だけが通る）
      },
      // ---- ハイドログラフ
      hydro: {
        mode: 'tri',     // 'const' | 'tri' | 'table'
        Qbase: 10,       // 基底流量 [m3/s]
        Qpeak: 60,       // ピーク流量 [m3/s]
        tPeak: 3600,     // ピーク時刻 [s]
        duration: 10800, // 継続時間 [s]
        Cs: 0.005,       // 供給土砂の流砂量濃度 Qs/Q（実質体積） [-]
        // mode 'debris'：土石流の一波。基底流量（清水）に、t0 から三角形の波（ピーク Qpeak、立ち上がり tRise、減衰 tFall）を重ねる
        t0: 60, tRise: 30, tFall: 120, Cd: 0.4,   // 波の土砂濃度（容積）
        Iup: 0,          // 流入点より上流の勾配 1/n（流入する土石流の流速・水深を決める。0 なら上流端の勾配）
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
    this.tribs = [];
    if (T && P.tribs && P.tribs.length) (this.warn || (this.warn = [])).push('支川は、数値で指定した地形でだけ使えます（CSV の地形では使いません）。');
    var nx, ny, ns, i, j, c, k, q, Bx, yc, cosT;
    var maxH = 0; (P.dams || []).forEach(function (D) { if (D.H > maxH) maxH = D.H; });
    this.warn = this.warn || [];
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
      // 支川があれば、その側の谷壁の外側を支川の長さだけ広げる
      var nsL = ns, nsR = ns;
      (P.tribs || []).forEach(function (tb) { var nn = Math.ceil(Math.max(10, tb.len || 40) / dx) + 1; if (tb.side === 'R') nsR = Math.max(nsR, nn); else nsL = Math.max(nsL, nn); });
      ny = Math.round((ymax - ymin) / dx) + nsL + nsR;
      yc = new Float64Array(nx);
      for (i = 0; i < nx; i++) yc[i] = nsR * dx + (yo[i] - ymin);
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
    this.hOld = new Float64Array(N); this.dfC = new Float64Array(N); this.dfD = new Float64Array(N);   // 土石流：前の水深、濃度、代表粒径
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
    this.Ix = new Float64Array(nx);   // 各列の基準の縦断勾配（勾配の表があれば区間ごと）
    this.zbAt = null;
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
      for (i = 0; i < nx; i++) this.Ix[i] = this.I;
      var trA_ = this.trA, trB_ = this.trB; this.zbAt = function (x) { return trA_ + trB_ * x; };
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
      if (T.basin && T.basin.length === T.rows.length) this.setupBasin(T);
    } else {
      // 縦断：勾配の表（位置 x から先の勾配）があれば、区間ごとの勾配をつないだ折れ線にする
      var st = (P.slopeTable || []).filter(function (r) { return isFinite(r[0]) && r[1] > 0; }).slice().sort(function (a, b) { return a[0] - b[0]; });
      var Lx_ = this.Lx, segI = function (x) { var v = P.slope; for (var q2 = 0; q2 < st.length; q2++) if (x >= st[q2][0] - 1e-9) v = 1 / st[q2][1]; return v; };
      if (st.length && st[0][0] > 0) st.unshift([0, 1 / P.slope]);
      var brk = st.map(function (r) { return Math.max(0, Math.min(Lx_, r[0])); }).concat([Lx_]);
      var zbAt = this.zbAt = function (x) {      // 基準の河床高（下流端を 0 として上流へ積み上げる）
        if (!st.length) return (Lx_ - x) * P.slope;
        var zz = 0, xa = Math.max(0, Math.min(Lx_, x));
        for (var q2 = 0; q2 < st.length; q2++) { var a0 = brk[q2], b0 = brk[q2 + 1]; if (b0 <= xa) continue; zz += (b0 - Math.max(a0, xa)) / st[q2][1]; }
        return zz;
      };
      var I0 = this.I = st.length ? Math.max(1e-4, zbAt(0) / Lx_) : P.slope;   // 平均勾配
      this.trA = this.Lx * I0; this.trB = -I0;
      for (i = 0; i < nx; i++) this.Ix[i] = st.length ? segI((i + 0.5) * dx) : P.slope;
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
      var zbc = new Float64Array(nx); for (i = 0; i < nx; i++) zbc[i] = zbAt((i + 0.5) * dx);
      var plane = function (i_, o_) { return zbc[i_] + cs0 * Math.abs(o_ - ych[i_]); };
      for (i = 0; i < nx; i++) {
        var zc = zbc[i], hLi = this.hL[i], hRi = this.hR[i];
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
      // 支川：谷壁を切り込んだ水路（勾配 I）を、計算区間の側方の端から本川の谷底の縁までつなぐ。端のセルから流入させる
      this.tribs = [];
      var TR = (P.tribs || []).filter(function (tb) { return isFinite(tb.x) && tb.B > 0 && tb.I > 0; });
      for (q = 0; q < TR.length; q++) {
        // 合流角 ang：本川の流下方向と支川の流下方向のなす角（90° で直角、小さいほど本川の下流向きに斜めに入る）
        var tb = TR[q], sdT = tb.side === 'R' ? -1 : 1, It = 1 / tb.I, angD = isFinite(tb.ang) ? Math.max(30, Math.min(150, tb.ang)) : 90,
          aR = angD * Math.PI / 180, sA = angD === 90 ? 1 : Math.sin(aR), cA = angD === 90 ? 0 : Math.cos(aR), cot = cA / sA, hw = tb.B / 2, xj = tb.x;
        var inC = [], jbT = sdT > 0 ? ny - 1 : 0, zE = new Float64Array(nx);
        for (i = 0; i < nx; i++) zE[i] = plane(i, sdT > 0 ? this.hL[i] : -this.hR[i]);
        for (j = 0; j < ny; j++) for (i = 2; i < nx - 2; i++) {
          c = j * nx + i;
          var offT = ((j + 0.5) * dx - yc[i]) * cosT[i], dist = sdT > 0 ? offT - this.hL[i] : -this.hR[i] - offT;
          if (dist <= 0) continue;
          var xr = (i + 0.5) * dx - xj, pp = xr * sA + dist * cA;      // pp：支川の中心線からの距離（中心線に直角）
          var iw = Math.max(0, Math.min(nx - 1, Math.floor(((i + 0.5) * dx + dist * cot) / dx))), zb = zE[iw] + dist / sA * It;
          if (pp >= -hw && pp < hw) { z[c] = zb; kind[c] = 0; fl[c] = 1; bare[c] = 0; if (j === jbT) inC.push(c); }
          else if (kind[c] === 1) {
            // 支川の谷：支川の河床より低い谷壁を、水路から離れるほど高くなる V 字の谷（側壁の勾配 1:2）に盛る
            var lat = Math.abs(pp) - hw;
            if (lat < 30) { var ztb = zb + Math.min(1.5, dist * 0.5) + 0.5 * lat; if (z[c] < ztb) z[c] = ztb; }
          }
        }
        if (!inC.length) (this.warn || (this.warn = [])).push('支川（x = ' + tb.x + ' m）が計算区間の側方の端に届きません。合流角か位置を変えてください。');
        if (inC.length) this.tribs.push({ x: tb.x, side: sdT, I: It, B: tb.B, ang: angD, sA: sA, cA: cA, r: tb.r > 0 ? tb.r : 0, rs: tb.rs >= 0 ? tb.rs : (tb.r || 0), lag: tb.lag || 0, cells: inC, q: 0, v: 0, Q: 0, Qs: 0 });
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
    this.openC = new Int16Array(N);    // 透過型堰堤の開口のセル（堰堤の番号 + 1）
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
      var zbase = T ? colmin[id] : this.zbAt(cx);
      if (T && T.basin && !isFinite(D.zc)) {   // DEM：堰堤の中心のまわり（3×3）の最も低い地盤を基礎にする
        var ic = Math.floor(cx / dx), jc = Math.floor(cy / dx); zbase = 1e30;
        for (var dj0 = -1; dj0 <= 1; dj0++) for (var di0 = -1; di0 <= 1; di0++) { var ii = ic + di0, jj0 = jc + dj0; if (ii >= 0 && jj0 >= 0 && ii < nx && jj0 < ny && kind[jj0 * nx + ii] !== 2) zbase = Math.min(zbase, z[jj0 * nx + ii]); }
        if (!(zbase < 1e29)) zbase = colmin[id];
      }
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
      // ---- 透過型：水通しの中に開口を開ける（開口のセルは開口の底の高さの移動床。水と土砂は開口率に応じて通る）
      var OD = P.openDam || {}, typ = D.type || 'closed', op = null;
      if (typ !== 'closed' && OD) {
        var lw = D.l > 0 ? D.l : 1.0, wm = OD.w > 0 ? OD.w : 0.5, nO = Math.max(1, Math.round(D.n || 1));
        var zo = zbase + Math.max(0, Math.min(D.H - 0.1, D.hb || 0)), sc0 = 0.5 * (spMin + spMax), wz, phiW, phiS;
        if (typ === 'slit') { nO = 1; wz = lw; }
        else if (typ === 'mslit') wz = nO * lw + (nO - 1) * wm;
        else wz = spMax - spMin + dx;                                     // 鋼製：水通しの全幅
        var oc = [], kept = [], keptS = [];
        for (k = 0; k < cells.length; k++) {
          c = cells[k];
          if (spill[k]) {
            var pxo = ((c % nx) + 0.5) * dx - cx, pyo = (((c / nx) | 0) + 0.5) * dx - cy, slo = pxo * nxv + pyo * nyv;
            if (Math.abs(slo - sc0) <= wz / 2 + 1e-9) { oc.push(c); continue; }
          }
          kept.push(c); keptS.push(spill[k]);
        }
        if (!oc.length) {   // 開口が格子より狭いときは、水通しの中心に最も近いセルを 1 つ開ける
          var best = -1, bd = 1e30;
          for (k = 0; k < cells.length; k++) if (spill[k]) { c = cells[k]; var pxq = ((c % nx) + 0.5) * dx - cx, pyq = (((c / nx) | 0) + 0.5) * dx - cy, dd = Math.abs(pxq * nxv + pyq * nyv - sc0); if (dd < bd) { bd = dd; best = k; } }
          if (best >= 0) { oc.push(cells[best]); kept = []; keptS = []; for (k = 0; k < cells.length; k++) if (k !== best) { kept.push(cells[k]); keptS.push(spill[k]); } }
        }
        if (oc.length) {
          // 開口のセルが占める幅（堰堤軸方向）に対する、開口の幅の割合
          var cover = oc.length * dx / Math.max(1, Math.round(thick * 2 / dx));
          if (typ === 'slit' || typ === 'mslit') { phiW = Math.min(1, nO * lw / Math.max(cover, 1e-9)); phiS = phiW; }
          else { var fr1 = lw / (lw + wm); phiW = typ === 'grid' ? fr1 * fr1 : fr1; phiS = fr1; }
          for (k = 0; k < oc.length; k++) { c = oc[k]; z[c] = zo; zfix[c] = zo; kind[c] = 0; fl[c] = 1; this.openC[c] = this.dams.length + 1; }
          cells = kept; spill = keptS;
          var zs0 = (typ === 'slit' || typ === 'mslit') ? sc0 - wz / 2 : spMin - dx / 2, zs1 = (typ === 'slit' || typ === 'mslit') ? sc0 + wz / 2 : spMax + dx / 2;
          op = { type: typ, l: lw, n: nO, w: wm, zo: zo, cells: oc, phiW: phiW, phiS: phiS, blocked: false, tBlock: null, beta: D.beta > 0 ? D.beta : (OD.beta || 1.5), s0: zs0, s1: zs1, hb: zo - zbase };
        }
      }
      this.dams.push({ i: id, x: cx, cx: cx, cy: cy, th: th, tx: tx, ty: ty, nx: nxv, ny: nyv, ext: ext, thick: thick, H: D.H, b: bn,
        crest: crest, wingTop: wingTop, zbase: zbase, zlow: zlow, j0: j0, j1: j1, cells: cells, spill: spill,
        sMin: sMin - dx / 2, sMax: sMax + dx / 2, spMin: spMin - dx / 2, spMax: spMax + dx / 2, type: op ? typ : 'closed', op: op,
        hb: D.hb || 0 });
    }
    this.hasOpen = this.dams.some(function (Dd) { return !!Dd.op; });
    // 各セルが何番目の堰堤の上流側にあるか（堆砂量の集計・初期湛水に使う）
    var region = this.region = new Int16Array(N), nd = this.dams.length;
    if (T && T.basin && nd) this.damRegionDEM();
    else for (j = 0; j < ny; j++) for (i = 0; i < nx; i++) {
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
        var zt = Df.crest + up * this.Ix[Math.max(0, Df.i - 1)] * P.prefill;
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
      if (del <= 1e-9 || kind[c] === 3) continue;   // 斜面の土層（kind 3）は流れで侵食しない（層を持たない）
      var a = Math.min(del, La);
      for (k = 0; k < K; k++) this.A[c * K + k] = a * f0[k];
      var rem = del - a;
      if (rem > 1e-9) {
        var nl = Math.min(NL, Math.ceil(rem / Ls - 1e-9));
        this.nl[c] = nl; this.topTh[c] = rem - (nl - 1) * Ls;
        for (var l = 0; l < nl; l++) for (k = 0; k < K; k++) this.SF[(c * NL + l) * K + k] = f0[k];
      }
    }

    // 支川の流入セル
    this.tribIn = new Int16Array(N);
    for (q = 0; q < this.tribs.length; q++) for (k = 0; k < this.tribs[q].cells.length; k++) this.tribIn[this.tribs[q].cells[k]] = q + 1;
    this.setupPoints();
    // ---- 初期の水（等流の薄層流＋堰堤上流の湛水）
    var Q0 = this.demMode ? 0 : this.Qat(0);
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
      var rq2 = this.region[c]; if (rq2 >= this.dams.length || kind[c] === 2 || this.openC[c]) continue;
      if (this.dams[rq2].op) continue;   // 透過型の上流は湛水させない
      var cr = this.dams[rq2].crest;
      if (z[c] < cr && cr - z[c] > h[c]) { h[c] = cr - z[c]; qx[c] = 0; this.qy[c] = 0; }
    }

    // ---- 伏流：河床（堆積層）の中の水 W [m]。水に浸かったセルは飽和、乾いたセルは飽和度 s0
    this.W = new Float64Array(N); this.kgw = new Float64Array(N); this.gwOut = 0;
    var GW0 = P.gw || {};
    this.por = new Uint8Array(N); for (c = 0; c < N; c++) this.por[c] = (kind[c] === 0 || kind[c] === 3) ? 1 : 0;   // 伏流の層があるセル（河床と斜面の土層）
    if (GW0.on) for (c = 0; c < N; c++) if (this.por[c]) { var Tc = Math.max(0, z[c] - zfix[c]); this.W[c] = P.porosity * Tc * (h[c] > 1e-3 ? 1 : Math.max(0, Math.min(1, kind[c] === 3 && P.dem && P.dem.s0 >= 0 ? P.dem.s0 : (GW0.s0 || 0)))); }
    this.W0 = 0; for (c = 0; c < N; c++) this.W0 += this.W[c]; this.W0 *= dx * dx;
    this.voidAir = 0;   // 伏流なしの土石流で、堆積物の間隙を満たす水が足りなかった体積（空気）
    this.setupSlides();
    // ---- 時刻・収支
    this.t = 0;          // 出水の時刻（河床変動の時間軸） [s]
    this.th = 0;         // 流れの計算時間 [s]
    this.nstep = 0;
    this.dt = 0.01;
    this.phase = Q0 > 0.01 ? 'spinup' : 'run';
    this.spinT = 0;
    this.Qin = Q0; this.Qout = 0; this.QsIn = 0; this.QsOut = 0;
    this.Vwin = 0; this.Vwout = 0;
    this.rainVol = 0; this.rainLoss = 0; this.slideLog = [];
    this.sedIn = 0; this.sedOut = 0;     // 実質体積 [m3]
    this.sedOutK = new Float64Array(K);
    this.nLimW = 0; this.nLimS = 0;
    this.series = { t: [], Qin: [], Qout: [], QsIn: [], QsOut: [], stor: [], QinT: [], QsInT: [] };   // QinT・QsInT：全体の流入（上流端＋支川＋供給点＋雨）
    this.snaps = []; this.snapFull = false; this._nextSnap = 0;   // 時間ごとの結果（記録は出水の開始から）
    this._nextSample = 0;
    this.Uref = Q0 > 0.01 ? Math.max(0.5, (Q0 / this.Bf) / Math.pow((Q0 / this.Bf) * P.manning / Math.sqrt(Math.max(I, 1e-4)), 0.6)) : 1;
    // 天然ダムを置いたあとの地形を初期の地形とする
    this.z0.set(z); this.Vw0 = this.waterVolume(); this.W0 = this.groundWater();
  };

  // ------------------------------------------------------------ ハイドログラフ
  Sim.prototype.duration = function () {
    if (this.demMode && this.P.dem && this.P.dem.duration > 0) return this.P.dem.duration;
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
  // 土石流の一波（三角形）の流量 [m3/s]（混合物）
  Sim.prototype.waveQ = function (t) {
    var H = this.P.hydro, a = t - (H.t0 || 0), r = Math.max(1, H.tRise || 30), f = Math.max(1, H.tFall || 120);
    if (a <= 0 || a >= r + f) return 0;
    return (H.Qpeak - H.Qbase) * (a < r ? a / r : 1 - (a - r) / f);
  };
  Sim.prototype.Qat = function (t) { return this.demMode ? 0 : QatH(this.P.hydro, t); };
  Sim.prototype.QsAt = function (t) { return this.demMode ? 0 : QsAtH(this.P.hydro, t); };
  function waveQH(H, t) {
    var a = t - (H.t0 || 0), r = Math.max(1, H.tRise || 30), f = Math.max(1, H.tFall || 120);
    if (a <= 0 || a >= r + f) return 0;
    return (H.Qpeak - H.Qbase) * (a < r ? a / r : 1 - (a - r) / f);
  }
  // ハイドロ H（上流端と同じ形の設定）の時刻 t の流量 [m3/s]
  function QatH(H, t) {
    if (H.mode === 'debris') return H.Qbase + waveQH(H, t);
    if (H.mode === 'table' && H.table && H.table.length > 1) return Math.max(0, interp(H.table, t, 1));
    if (H.mode === 'const') return H.Qpeak;
    if (t <= 0) return H.Qbase;
    if (t < H.tPeak) return H.Qbase + (H.Qpeak - H.Qbase) * t / H.tPeak;
    if (t < H.duration) return H.Qpeak + (H.Qbase - H.Qpeak) * (t - H.tPeak) / Math.max(1e-9, H.duration - H.tPeak);
    return H.Qbase;
  };
  // 供給土砂量（実質体積） [m3/s]
  function QsAtH(H, t) {
    if (H.mode === 'table' && H.table && H.table.length > 1 && H.table[0].length > 2 && isFinite(H.table[0][2]))
      return Math.max(0, interp(H.table, t, 2));
    if (H.mode === 'debris') return (H.Cd || 0) * waveQH(H, t) + H.Cs * H.Qbase;
    return H.Cs * QatH(H, t);
  }

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
    var wsum = 0, sqI = Math.sqrt(Math.max(this.Ix[0], 1e-4)), nman = P.manning;
    var fl = this.fl;
    for (j = 0; j < ny; j++) { var w = fl[j * nx] ? Math.pow(h[j * nx] + 0.05, 5 / 3) : 0; qin[j] = w; wsum += w; }
    var dbIn = P.transport === 'debris', Cin0 = dbIn && Q > 0 ? Math.min(0.9 * (1 - P.porosity), this.QsAt(this.t) / Q) : 0;
    var Iu = (P.hydro.Iup > 0 ? 1 / P.hydro.Iup : this.Ix[0]);
    for (j = 0; j < ny; j++) {
      var qj = wsum > 0 ? Q * qin[j] / wsum / dx : 0;   // 単位幅流量 [m2/s]
      qin[j] = qj;
      uin[j] = qj > 0 ? (dbIn ? qj / this.debrisDepth(qj, Cin0, Iu) : qj / Math.pow(qj * nman / Math.sqrt(Math.max(Iu, 1e-4)), 0.6)) : 0;
    }
    // 支川の流入（本川の流量の r 倍、lag 秒遅れ）
    var TRB = this.tribs, QtSum = 0;
    for (var tq_ = 0; tq_ < TRB.length; tq_++) {
      var tb = TRB[tq_], tl = Math.max(0, this.t - tb.lag), Qt = tb.r * this.Qat(tl);
      tb.Q = Qt; tb.q = Qt / (tb.cells.length * dx); QtSum += Qt;
      var Ct = dbIn && Qt > 0 && this.phase === 'run' ? Math.min(0.9 * (1 - P.porosity), tb.rs * this.QsAt(tl) / Qt) : 0;
      var qw = tb.q * tb.sA;   // 支川の単位幅流量（斜めに入ると、端のセルの列は幅 B/sin になる）
      tb.v = qw > 0 ? (dbIn ? qw / this.debrisDepth(qw, Ct, tb.I) : qw / Math.pow(qw * nman / Math.sqrt(tb.I), 0.6)) : 0;
    }

    // 点の流入（谷頭の供給点・雨の横流入）と、格子に降る雨
    var PIN = this.pin, QpSum = PIN && PIN.length ? this.pinStep() : 0, rainH = null;

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
      var a, b, Idx = this.Ix[nx - 1] * dx;
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
    if (this.rainOn) rainH = this.rainStep(dt);   // 格子に降る雨（このステップの有効雨量の水深 [m]）
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

    // ---- 土石流：濃度と代表粒径（抵抗則に使う）
    var dbr = P.transport === 'debris', dfC = this.dfC, dfD = this.dfD, DB = P.debris || {}, Cst_ = 1 - P.porosity;
    if (dbr) { this.hOld.set(h); this.debrisParams(); }
    var rhoS = 1000 * P.sg, rhoW = DB.rho || 1000, mfl = DB.manningFloor !== false;
    // ---- 更新（勾配項・流入・摩擦）
    var pinIn = PIN && PIN.length ? this.pinIn : null;
    var n2 = this.n2, hn, qxn, qyn, Vv, fr, qoutSum = 0, tribIn = this.tribIn, openC = this.openC, dams = this.dams, phiB = (P.openDam && P.openDam.phiBlock > 0) ? P.openDam.phiBlock : 0.05, frX;
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
        if (tribIn[c]) { var tbc = TRB[tribIn[c] - 1]; hn += cf * tbc.q; qyn -= cf * tbc.q * tbc.v * tbc.sA * tbc.side; qxn += cf * tbc.q * tbc.v * tbc.cA; }
        if (pinIn && pinIn[c]) { hn += dt * this.pH[c]; qxn += dt * this.pMx[c]; qyn += dt * this.pMy[c]; }
        if (rainH && rainH[c] > 0) hn += rainH[c];
        if (hn < 0) hn = 0;
        if (hn > HDRY) {
          // 摩擦項は陰的（|q| の2次方程式を解く）
          Vv = Math.sqrt(qxn * qxn + qyn * qyn);
          frX = 4 * dt * g * n2 * Vv / (hn * hn * Math.cbrt(hn));
          if (dbr && dfC[c] >= 0.01) {
            // 高橋の抵抗則：τ/ρ_T = Ψ u|V|。石礫型（C ≥ 0.4C*）Ψ = d²/{8 h² (C + (1−C)ρ/σ)((C*/C)^(1/3) − 1)²}、乱流型（0.01 ≤ C < 0.4C*）Ψ = d²/(0.49 h²)
            var Cd_ = Math.min(dfC[c], 0.98 * Cst_), dd_ = dfD[c], psi;
            if (Cd_ >= 0.4 * Cst_) { var X_ = Math.cbrt(Cst_ / Cd_) - 1; psi = dd_ * dd_ / (8 * hn * hn * (Cd_ + (1 - Cd_) * rhoW / rhoS) * X_ * X_); }
            else psi = dd_ * dd_ / (0.49 * hn * hn);
            var frD = 4 * dt * psi * Vv / (hn * hn);
            frX = mfl ? Math.max(frX, frD) : frD;
          }
          if (openC[c]) {
            // 透過型堰堤の開口：多孔板の損失係数（Idelchik）ζ = {(1 + 0.707√(1−φ) − φ)/φ}²、損失水頭 ζV²/2g を 1 セルで与える。
            // 天端より上を流れる分は損失なし（φ_eff は水深で重み付け）
            var Dq_ = dams[openC[c] - 1], op_ = Dq_.op, ph = op_.blocked ? phiB : op_.phiW, Hs = z[c] + hn;
            var hsub = Math.max(0, Math.min(Hs, Dq_.crest) - z[c]), hov = Math.max(0, Hs - Dq_.crest), pe = (ph * hsub + hov) / hn;
            if (pe < 1) { if (pe < 1e-3) pe = 1e-3; var zt_ = (1 + 0.707 * Math.sqrt(1 - pe) - pe) / pe; frX += 4 * dt * (zt_ * zt_ / (2 * dx * hn)) * Vv; }
          }
          fr = 0.5 * (1 + Math.sqrt(1 + frX));
          qxn /= fr; qyn /= fr;
        } else { qxn = 0; qyn = 0; }
        h[c] = hn; qx[c] = qxn; qy[c] = qyn;
      }
    }
    this.Qout = qoutSum * dx;
    this.Vwin += (Q + QtSum + QpSum) * dt + (rainH ? this._rainV : 0); this.Vwout += this.Qout * dt;
    this.QinT = Q + QtSum + QpSum + (rainH ? this._rainV / dt : 0);
    var gwOn_ = !!(P.gw && P.gw.on);
    this.th += dt; this.nstep++;

    // ---- 助走計算（水だけを流して流れを立ち上げる）
    if (this.phase === 'spinup') {
      this.spinT += dt;
      var tmin = 0.5 * this.Lx / this.Uref, tmax = 3 * this.Lx / this.Uref;
      if ((this.spinT > tmin && Math.abs(this.Qout - Q) < 0.08 * Q) || this.spinT > tmax) {
        this.phase = 'run'; this.Vwin = 0; this.Vwout = 0; this.Vw0 = this.waterVolume();
        this.W0 = this.groundWater(); this.gwOut = 0; this.voidAir = 0;
        this.snapshot();
      }
      if (gwOn_) this.gwStep(dt);
      return;
    }
    if (this.phase === 'done') return;

    // ---- 河床変動
    if (this.t === 0 && !this.snaps.length) this.snapshot();   // 出水の開始時の状態
    var dtm = dt * (P.transport === 'debris' ? 1 : Math.max(1, P.morfac));   // 土石流は加速しない
    var sIn0 = this.sedIn;
    this.morph(dtm);
    if (gwOn_) this.gwStep(dt);      // 伏流は流れの時間刻みで（加速しない）
    if (this.slides.length || this.depOn) this.slideStep(dt);
    if (this.sinkCells && this.sinkCells.length) this.sinkStep(dt);
    this.QsInT = (this.sedIn - sIn0) / dtm;
    this.t += dtm;
    if (this.t >= this._nextSample) {
      var S = this.series;
      S.t.push(this.t); S.Qin.push(Q); S.Qout.push(this.Qout); S.QsIn.push(this.QsIn); S.QsOut.push(this.QsOut); S.QinT.push(this.QinT); S.QsInT.push(this.QsInT);
      S.stor.push(this.storage());
      this._nextSample = this.t + this.duration() / 400;
    }
    if (this.t >= this.duration()) this.phase = 'done';
    var ev = +P.snapEvery || 0;
    if (ev > 0 && (this.t >= this._nextSnap - 1e-9 || (this.phase === 'done' && this.snaps.length && this.snaps[this.snaps.length - 1].t < this.t))) this.snapshot();
  };

  // ---- 時間ごとの結果の記録 ---------------------------------------------------------
  // 河床変動量・水深・流速・表層の平均粒径を 4 バイトの実数で残す（河床高は z0 + 河床変動量で戻す）
  var SNAP_MAX = 25e6;   // 記録する値の総数の上限（約 100 MB）
  Sim.prototype.snapshot = function () {
    var ev = +this.P.snapEvery || 0;
    if (!(ev > 0)) return;
    this._nextSnap = (Math.floor(this.t / ev + 1e-9) + 1) * ev;
    if (this.snaps.length && Math.abs(this.snaps[this.snaps.length - 1].t - this.t) < 1e-9) return;
    var N = this.N, nx = this.nx;
    if ((this.snaps.length + 1) * N * 4 > SNAP_MAX) { this.snapFull = true; return; }
    var dz = new Float32Array(N), hh = new Float32Array(N), vv = new Float32Array(N), dm = new Float32Array(N);
    var z = this.z, z0 = this.z0, h = this.h, FxM = this.FxM, GyM = this.GyM;
    for (var c = 0; c < N; c++) {
      dz[c] = z[c] - z0[c]; hh[c] = h[c];
      if (h[c] > 1e-3) { var i = c % nx, j = (c - i) / nx, f = j * (nx + 1) + i, u = 0.5 * (FxM[f] + FxM[f + 1]) / h[c], v = 0.5 * (GyM[c] + GyM[c + nx]) / h[c]; vv[c] = Math.sqrt(u * u + v * v); }
      dm[c] = this.kind[c] === 0 ? this.surfaceDm(c) : NaN;
    }
    this.snaps.push({ t: this.t, dz: dz, h: hh, V: vv, dm: dm, Qin: this.Qin, Qout: this.Qout, QsIn: this.QsIn, QsOut: this.QsOut,
      sedIn: this.sedIn, sedOut: this.sedOut, stor: this.storage(), byDam: this.storageByDam() });
  };
  // 記録の間隔を途中で変える（次の記録は新しい間隔の倍数の時刻）
  Sim.prototype.setSnapEvery = function (sec) {
    this.P.snapEvery = sec > 0 ? sec : 0;
    this._nextSnap = sec > 0 ? (Math.floor(this.t / sec + 1e-9) + 1) * sec : Infinity;
    if (sec > 0 && this.phase !== 'spinup' && !this.snaps.length) this.snapshot();
  };

  // --------------------------------------------------------- 河床変動 1 ステップ
  Sim.prototype.morph = function (dtm) {
    this.dzs.fill(0);
    if (this.P.transport === 'debris') this.morphDebris(dtm);
    else if (this.P.transport === 'noneq') this.morphNoneq(dtm); else this.morphEq(dtm);
    this.dryBankErosion();
    this.hasegawaCollapse();
    this.collapse();
    if (this.demMode && this.P.dem && this.P.dem.pit > 0) this.pitFill(this.P.dem.pit);
  };
  // DEM：1 セルの穴（4 近傍の河床のどれよりも pit 以上低いセル）を、近傍の最も高い侵食できるセルの土砂で埋める（質量保存の数値処理）。
  // 1 セル幅の谷筋で侵食が 1 セルだけ進んだ穴に水がたまり、水深が突出するのを防ぐ。堰堤の隣は埋めない（水叩きの洗掘を残す）
  Sim.prototype.pitFill = function (tol) {
    var nx = this.nx, ny = this.ny, K = this.K, A = this.A, z = this.z, kind = this.kind, bas = this.basin, i, j, c, k, m, q;
    var self = this, cand = this._pitC || (this._pitC = []); cand.length = 0;
    // 候補：4 近傍（壁・流域の外を除く）のどれよりも tol 以上低い谷筋のセル。堰堤の隣は除く
    function zmin4(c) {
      var nb = [c - 1, c + 1, c - nx, c + nx], zm = 1e30;
      for (var m = 0; m < 4; m++) { var o = nb[m]; if (kind[o] === 2) return -1e30; if (kind[o] === 1 || (bas && !bas[o])) continue; if (z[o] < zm) zm = z[o]; }
      return zm;
    }
    for (j = 1; j < ny - 1; j++) for (i = 1; i < nx - 1; i++) {
      c = j * nx + i; if (kind[c] !== 0 || (bas && !bas[c])) continue;
      if (z[c] < zmin4(c) - tol) cand.push(c);
    }
    for (q = 0; q < cand.length; q++) {
      c = cand[q]; var zm = zmin4(c); if (!(z[c] < zm - tol)) continue;
      var nb = [c - 1, c + 1, c - nx, c + nx], don = -1, zd = -1e30, avd = 0;
      for (m = 0; m < 4; m++) {
        var o = nb[m]; if (kind[o] === 1 || (bas && !bas[o])) continue;
        var av = 0; for (k = 0; k < K; k++) av += A[o * K + k];
        if (av > 1e-6 && z[o] > zd) { zd = z[o]; don = o; avd = av; }
      }
      if (don < 0) continue;
      var mv = Math.min(0.5 * (zm - tol - z[c]), 0.5 * avd, 0.5 * (z[don] - z[c])); if (!(mv > 1e-6)) continue;
      var fr = mv / avd;
      for (k = 0; k < K; k++) { var t = A[don * K + k] * fr; A[don * K + k] -= t; A[c * K + k] += t; }
      z[don] -= mv; z[c] += mv; self.dzs[don] -= mv; self.dzs[c] += mv; self.fixLayers(don); self.fixLayers(c);
      self.nPit = (self.nPit || 0) + 1;
    }
  };

  // ---- 土石流（高橋 2007。Kanako 2D と同じ形の式）--------------------------------------
  // 平衡濃度 C∞（tanθ は水面勾配）。tanθ ≥ tanφ：0.9C*、> 0.138：石礫型、> 0.03：乱流型（6.7 倍の 2 乗）、≤ 0.03：掃流
  Sim.prototype.cinfT = function (t, h, d) {
    var P = this.P, DB = P.debris || {}, Cst = 1 - P.porosity, tp = DB.tanPhi > 0 ? DB.tanPhi : 0.7, rho = DB.rho || 1000, sig = 1000 * P.sg;
    if (!(t > 0)) return 0;
    if (t >= tp) return 0.9 * Cst;
    var c1 = rho * t / ((sig - rho) * (tp - t)), C;
    if (t > 0.138) C = c1;
    else if (t > 0.03) C = 6.7 * c1 * c1;
    else {
      var sgr = sig / rho, ts = h * t / ((sgr - 1) * d), tc = 0.04 * Math.pow(10, 1.72 * t);
      if (!(ts > tc)) return 0;
      var st = sig * t / (sig - rho), a2 = 2 * (0.425 - st) / (1 - st);
      if (!(a2 > 0)) return 0;
      C = (1 + 5 * t) * t / (sgr - 1) * (1 - a2 * tc / ts) * (1 - Math.sqrt(a2 * tc / ts));
    }
    return C < 0 ? 0 : C > 0.9 * Cst ? 0.9 * Cst : C;
  };
  // 土石流：開口で止めた土砂の体積を、混合物の水深として上流側に残す（下流側へは水だけが通る）
  Sim.prototype.holdVolume = function (dX, dY, cf) {
    var nx = this.nx, ny = this.ny, h = this.h, j, i, t, a, b;
    for (j = 0; j < ny; j++) for (i = 1; i < nx; i++) {
      if (dX[j * (nx + 1) + i] === 0) continue;
      a = j * nx + i - 1; b = a + 1; t = cf * dX[j * (nx + 1) + i];
      t = t > 0 ? Math.min(t, h[b]) : -Math.min(-t, h[a]);
      h[a] += t; h[b] -= t;
    }
    for (j = 1; j < ny; j++) for (i = 0; i < nx; i++) {
      if (dY[j * nx + i] === 0) continue;
      b = j * nx + i; a = b - nx; t = cf * dY[b];
      t = t > 0 ? Math.min(t, h[b]) : -Math.min(-t, h[a]);
      h[a] += t; h[b] -= t;
    }
  };
  // 各セルの流れの土砂濃度と代表粒径（流れている土砂の平均粒径。土砂がなければ河床表層の値）
  Sim.prototype.debrisParams = function () {
    var N = this.N, K = this.K, h = this.h, Cc = this.Cc, d = this.d, dfC = this.dfC, dfD = this.dfD, Cst = 1 - this.P.porosity;
    for (var c = 0; c < N; c++) {
      var S = 0, Sd = 0, ck = c * K;
      for (var k = 0; k < K; k++) { S += Cc[ck + k]; Sd += Cc[ck + k] * d[k]; }
      var hh = h[c];
      dfC[c] = hh > 1e-3 ? Math.min(0.98 * Cst, S / hh) : 0;
      if (S > 1e-9) dfD[c] = Sd / S; else { var dm = this.surfaceDm(c); dfD[c] = dm === dm ? dm : this.dm0; }
    }
  };
  // 流入する土石流の等流水深（単位幅流量 q、濃度 C、勾配 I）
  Sim.prototype.debrisDepth = function (q, C, I) {
    var P = this.P, DB = P.debris || {}, Cst = 1 - P.porosity, rhoS = 1000 * P.sg, rhoW = DB.rho || 1000, nman = P.manning;
    var sn = I / Math.sqrt(1 + I * I), ff = this.ff, d = this.d, dd = 0, k, mfl = DB.manningFloor !== false;
    for (k = 0; k < this.K; k++) dd += ff[k] * d[k];
    var Cd = Math.min(C, 0.98 * Cst);
    var vel = function (h) {
      var Vm = Math.pow(h, 2 / 3) * Math.sqrt(Math.max(I, 1e-4)) / nman;
      if (Cd < 0.01) return Vm;
      var psi;
      if (Cd >= 0.4 * Cst) { var X = Math.cbrt(Cst / Cd) - 1; psi = dd * dd / (8 * h * h * (Cd + (1 - Cd) * rhoW / rhoS) * X * X); }
      else psi = dd * dd / (0.49 * h * h);
      var Vd = Math.sqrt(G * h * sn / psi);
      return mfl ? Math.min(Vd, Vm) : Vd;
    };
    var lo = 1e-4, hi = 50;
    for (var it = 0; it < 60; it++) { var mid = Math.sqrt(lo * hi); if (mid * vel(mid) > q) hi = mid; else lo = mid; }
    return Math.sqrt(lo * hi);
  };
  // 土石流の河床変動：混合物の流束で土砂を運び（風上）、高橋の侵食・堆積速度式で河床と交換する
  //   侵食（C < C∞）：i = δe (C∞ − C)/(C* − C∞)·q/d、堆積（C ≥ C∞）：i = δd (C∞ − C)/C*·V
  //   ∂z/∂t = −i、∂(C h)/∂t + ∇·(C q) = C* i、∂h/∂t + ∇·q = i（河床は飽和）
  Sim.prototype.morphDebris = function (dt) {
    var nx = this.nx, ny = this.ny, N = this.N, dx = this.dx, K = this.K, P = this.P, DB = P.debris || {};
    var h = this.h, hO = this.hOld, z = this.z, zfix = this.zfix, A = this.A, Cc = this.Cc, d = this.d, FxM = this.FxM, GyM = this.GyM;
    var SFx = this.SFx, SFy = this.SFy, fxC = this.fixedC, Cst = 1 - P.porosity, HSED = 0.02, ff = this.ff, qin = this.qin;
    var de = DB.de > 0 ? DB.de : 0.0007, dd = DB.dd > 0 ? DB.dd : 0.05, Vmin = DB.Vmin >= 0 ? DB.Vmin : 0.1, gwOn = !!(P.gw && P.gw.on);
    var eroMode = DB.erosion || 'sat', Kt = DB.Kt > 0 ? DB.Kt : 0.0019, tanPd = DB.tanPhi > 0 ? DB.tanPhi : 0.7, sgm = P.sg * 1000 / (DB.rho || 1000);
    var i, j, c, k, ck, row, frow, fk, L, R, f, cu, hu;
    // ---- 1. 移流（面の混合物流束 × 風上セルの濃度）
    for (j = 0; j < ny; j++) {
      row = j * nx; frow = j * (nx + 1);
      for (k = 0; k < K; k++) SFx[frow * K + k] = 0;
      for (i = 1; i <= nx; i++) {
        f = frow + i; fk = f * K; L = row + i - 1; R = L + 1;
        var F = FxM[f];
        if (i === nx) cu = F > 0 ? L : -1; else cu = F > 0 ? L : R;
        hu = cu >= 0 ? hO[cu] : 0;
        if (cu < 0 || hu < 1e-6) { for (k = 0; k < K; k++) SFx[fk + k] = 0; continue; }
        for (k = 0; k < K; k++) SFx[fk + k] = F * Cc[cu * K + k] / hu;
      }
    }
    for (i = 0; i < nx; i++) for (k = 0; k < K; k++) { SFy[i * K + k] = 0; SFy[(ny * nx + i) * K + k] = 0; }
    for (j = 1; j < ny; j++) {
      row = j * nx;
      for (i = 0; i < nx; i++) {
        R = row + i; L = R - nx; fk = R * K;
        var Gf = GyM[R]; cu = Gf > 0 ? L : R; hu = hO[cu];
        if (hu < 1e-6) { for (k = 0; k < K; k++) SFy[fk + k] = 0; continue; }
        for (k = 0; k < K; k++) SFy[fk + k] = Gf * Cc[cu * K + k] / hu;
      }
    }
    var cf = dt / dx, scS = this.scS, lim = false;
    if (this.hasOpen) {
      var nfx = ny * (nx + 1), nfy = (ny + 1) * nx, sx0 = new Float64Array(nfx), sy0 = new Float64Array(nfy), q_;
      for (q_ = 0; q_ < nfx; q_++) for (k = 0; k < K; k++) sx0[q_] += SFx[q_ * K + k];
      for (q_ = 0; q_ < nfy; q_++) for (k = 0; k < K; k++) sy0[q_] += SFy[q_ * K + k];
      this.checkBlock(true); this.openFilter(SFx, SFy);
      for (q_ = 0; q_ < nfx; q_++) for (k = 0; k < K; k++) sx0[q_] -= SFx[q_ * K + k];
      for (q_ = 0; q_ < nfy; q_++) for (k = 0; k < K; k++) sy0[q_] -= SFy[q_ * K + k];
      this.holdVolume(sx0, sy0, cf);
    }
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
    var Qs = this.phase === 'run' ? this.QsAt(this.t) : 0, Qw = this.Qin;
    this.QsIn = Qs; this.sedIn += Qs * dt;
    var Cin = (Qs > 0 && Qw > 0) ? Qs / Qw : 0, qsOut = 0;
    for (j = 0; j < ny; j++) {
      row = j * nx; frow = j * (nx + 1);
      fk = (frow + nx) * K;
      for (k = 0; k < K; k++) { qsOut += SFx[fk + k]; this.sedOutK[k] += SFx[fk + k] * dx * dt; }
      for (i = 0; i < nx; i++) {
        c = row + i; ck = c * K; var fw2 = (frow + i) * K, fe2 = fw2 + K, gs2 = ck, gn2 = (c + nx) * K;
        for (k = 0; k < K; k++) {
          var nv = Cc[ck + k] - cf * (SFx[fe2 + k] - SFx[fw2 + k] + SFy[gn2 + k] - SFy[gs2 + k]);
          if (i === 0 && Cin > 0 && qin[j] > 0) nv += cf * Cin * qin[j] * ff[k];
          Cc[ck + k] = nv > 0 ? nv : 0;
        }
      }
    }
    this.QsOut = qsOut * dx; this.sedOut += qsOut * dx * dt;
    this.tribFeed(dt, false);

    // ---- 2. 河床との交換
    var eta = this.eta, wf = this._wf || (this._wf = new Uint8Array(N));
    for (c = 0; c < N; c++) { eta[c] = z[c] + h[c]; wf[c] = h[c] >= HSED ? 1 : 0; }   // 勾配は交換の前の水面・濡れ方で求める
    for (c = 0; c < N; c++) {
      ck = c * K; var hh = h[c], S = 0;
      for (k = 0; k < K; k++) S += Cc[ck + k];
      if (hh < HSED) {      // 乾いたら、流れの中の土砂はその場に堆積させる
        if (S > 0) {
          var dz0 = S / Cst;
          for (k = 0; k < K; k++) { A[ck + k] += Cc[ck + k] / Cst; Cc[ck + k] = 0; }
          z[c] += dz0; this.dzs[c] += dz0; this.fixLayers(c);
          if (gwOn) { var pw0 = Math.min((1 - Cst) * dz0, Math.max(0, hh - S)); this.W[c] += pw0; h[c] = Math.max(0, hh - S - pw0); }   // 間隙の水は手元の水の分だけ
          else { h[c] = hh - dz0; if (h[c] < 0) { this.voidAir -= h[c] * dx * dx; h[c] = 0; } }
        }
        continue;
      }
      if (fxC[c]) continue;
      var ic = c % nx, jc = (c - ic) / nx, fwq = jc * (nx + 1) + ic;
      var uu = 0.5 * (FxM[fwq] + FxM[fwq + 1]) / hh, vv = 0.5 * (GyM[c] + GyM[c + nx]) / hh, V = Math.sqrt(uu * uu + vv * vv);
      // 流れの向きの水面勾配
      var tq = 0;
      if (V > 1e-6) {
        var wl = ic > 0 && wf[c - 1], wr = ic < nx - 1 && wf[c + 1], ws = jc > 0 && wf[c - nx], wn = jc < ny - 1 && wf[c + nx];
        // 流れの下流側の片側差分（中心差分では 1 セルおきの凹凸が見えず、堆積面が市松模様になる）
        var gx = uu >= 0 ? (wr ? (eta[c + 1] - eta[c]) / dx : wl ? (eta[c] - eta[c - 1]) / dx : -this.Ix[ic])
                         : (wl ? (eta[c] - eta[c - 1]) / dx : wr ? (eta[c + 1] - eta[c]) / dx : -this.Ix[ic]);
        var gy = vv >= 0 ? (wn ? (eta[c + nx] - eta[c]) / dx : ws ? (eta[c] - eta[c - nx]) / dx : 0)
                         : (ws ? (eta[c] - eta[c - nx]) / dx : wn ? (eta[c + nx] - eta[c]) / dx : 0);
        tq = -(gx * uu + gy * vv) / V;
      }
      var C = S / hh, As = 0, dmB = 0;
      for (k = 0; k < K; k++) { As += A[ck + k]; dmB += A[ck + k] * d[k]; }
      var dRep = C > 1e-6 ? this.dfD[c] : (As > 1e-9 ? dmB / As : this.dm0);
      if (!(dRep > 0)) dRep = this.dm0;
      var Ci = this.cinfT(tq, hh, dRep), dz = 0;
      if (C < Ci) {
        if (As < 1e-9) continue;
        var dB = dmB / As;
        var unsat = eroMode === 'unsat' || (eroMode === 'auto' && gwOn && this.por[c] && this.W[c] < 0.9 * (1 - Cst) * Math.max(0, z[c] - zfix[c]));
        if (unsat && tq > 1e-4 && tq < tanPd) {
          // 不飽和河床の侵食速度式（高橋・里深）：i/√(gh) = Kt sin^{3/2}θ {1 − (σ−ρ)/ρ C (tanφ/tanθ − 1)}^{1/2} (tanφ/tanθ − 1)(C∞ − C) h/d
          var tr_ = tanPd / tq - 1, br_ = 1 - (sgm - 1) * C * tr_, sn_ = tq / Math.sqrt(1 + tq * tq);
          dz = br_ > 0 ? Kt * Math.sqrt(G * hh) * Math.pow(sn_, 1.5) * Math.sqrt(br_) * tr_ * (Ci - C) * hh / dB * dt : 0;
        } else dz = de * (Ci - C) / (Cst - Ci) * hh * V / dB * dt;
        var lim1 = (Ci * hh - S) / (Cst - Ci); if (dz > lim1) dz = lim1;
        if (dz > As) dz = As;
        if (dz > z[c] - zfix[c]) dz = Math.max(0, z[c] - zfix[c]);
        if (!(dz > 0)) continue;
        for (k = 0; k < K; k++) { var e = dz * A[ck + k] / As; A[ck + k] -= e; Cc[ck + k] += Cst * e; }
        if (gwOn) {   // 削った層に含まれる地下水だけが流れに入る（乾いた層なら土砂だけ）
          var Tq = z[c] - zfix[c], bq = this.W[c] / (1 - Cst), win = (1 - Cst) * Math.max(0, Math.min(dz, bq - (Tq - dz)));
          this.W[c] -= win; h[c] += Cst * dz + win;
        } else h[c] += dz;
        z[c] -= dz; this.dzs[c] -= dz;
      } else if (C > Ci && S > 0) {
        dz = dd * (C - Ci) / Cst * Math.max(V, Vmin) * dt;
        var lim2 = (S - Ci * hh) / (Cst - Ci); if (dz > lim2) dz = lim2;
        if (dz > S / Cst) dz = S / Cst;
        if (!(dz > 0)) continue;
        for (k = 0; k < K; k++) { var gq = dz * Cc[ck + k] / S; A[ck + k] += gq; Cc[ck + k] -= Cst * gq; if (Cc[ck + k] < 0) Cc[ck + k] = 0; }
        z[c] += dz; this.dzs[c] += dz;
        if (gwOn) { var pw = Math.min((1 - Cst) * dz, Math.max(0, hh - S)); this.W[c] += pw; h[c] = hh - Cst * dz - pw; }   // 堆積物の間隙の水は河床の中へ
        else { h[c] = hh - dz; if (h[c] < 0) { this.voidAir -= h[c] * dx * dx; h[c] = 0; } }
      } else continue;
      this.fixLayers(c);
    }
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

    if (this.hasOpen) { this.checkBlock(false); this.openFilter(SFx, SFy); }
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
    this.tribFeed(dtm, true);
  };

  // 安息角を超えた斜面の崩落。水に浸かった斜面は水中安息角、乾いた斜面は陸上安息角で判定する
  // （高橋・里深・奥村 1996 の 30° 規則を、XBeach・Delft3D の wet/dry slope の考え方で拡張）
  Sim.prototype.collapse = function () {
    var nx = this.nx, ny = this.ny, K = this.K, dx = this.dx, A = this.A, z = this.z, h = this.h, fxC = this.fixedC, BK = this.P.bank || {};
    var limW = Math.tan((BK.reposeWet || this.P.reposeDeg) * Math.PI / 180) * dx, limD = Math.tan((BK.reposeDry || this.P.reposeDeg) * Math.PI / 180) * dx;
    var HW = 0.01, i, j, c, k, row, lim, diff, hi, lo, av, mv, hk, lk, frc, hase = BK.collapse === 'hasegawa', kind = this.kind, noC = this.sinkMask;
    for (var pass = 0; pass < 2; pass++) {
      var off = pass === 0 ? 1 : nx;
      var iEnd = pass === 0 ? nx - 1 : nx, jEnd = pass === 0 ? ny : ny - 1;
      for (j = 0; j < jEnd; j++) {
        row = j * nx;
        for (i = 0; i < iEnd; i++) {
          c = row + i;
          if (fxC[c] || fxC[c + off] || (noC && (noC[c] || noC[c + off]))) continue;   // 抜く所（出口の先・流域界の外）へは崩さない
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

  // ---- 支川からの給砂（本川の供給土砂量の rs 倍、lag 秒遅れ）。bulk：河床に置く（平衡流砂量）、そうでなければ流れの中へ
  Sim.prototype.tribFeed = function (dtm, bulk) {
    if (this.phase === 'run' && this.pin && this.pin.length) this.pinFeed(dtm, bulk);
    var TRB = this.tribs; if (!TRB.length || this.phase !== 'run') return;
    var K = this.K, ff = this.ff, A = this.A, Cc = this.Cc, dx = this.dx, lam = 1 - this.P.porosity;
    for (var q = 0; q < TRB.length; q++) {
      var tb = TRB[q], Qs = tb.rs * this.QsAt(Math.max(0, this.t - tb.lag));
      if (this.P.transport === 'debris' && Qs > 0.9 * (1 - this.P.porosity) * tb.Q) Qs = 0.9 * (1 - this.P.porosity) * tb.Q;   // 混合物の濃度は 0.9C* まで
      tb.Qs = Qs; if (!(Qs > 0)) continue;
      this.sedIn += Qs * dtm;
      var per = Qs * dtm / (tb.cells.length * dx * dx);     // 実質体積の層厚 [m]
      for (var m = 0; m < tb.cells.length; m++) {
        var c = tb.cells[m], ck = c * K;
        if (bulk) { for (var k = 0; k < K; k++) A[ck + k] += per / lam * ff[k]; this.z[c] += per / lam; this.dzs[c] += per / lam; this.fixLayers(c); }
        else for (var k2 = 0; k2 < K; k2++) Cc[ck + k2] += per * ff[k2];
      }
    }
  };
  // ======================================================================================
  // DEM の流域（SaboDEM.prepareDomain で切り出した地形）
  //  流域の外：高い壁（水が出入りしない）。谷筋（幅 chW に広げる）：侵食できる河床（kind 0）。
  //  斜面：kind 3（流れでは侵食しない。土層の厚さ soilD は伏流と斜面崩壊の土層として使う）
  // ======================================================================================
  // 4 近傍の Priority-Flood（Barnes ら 2014）＋ε の勾配。流域のセルだけ。種は右端の列と出口の先のセル。同じ高さは番号の小さい順
  function fill4(z, bas, seed, nx, ny, eps) {
    var N = nx * ny, done = new Uint8Array(N), hz = [], hc = [], n = 0, cnt = 0;
    function less(a, b) { return hz[a] < hz[b] || (hz[a] === hz[b] && hc[a] < hc[b]); }
    function push(zz, c) {
      var i = n++; hz[i] = zz; hc[i] = c;
      while (i > 0) { var p = (i - 1) >> 1; if (!less(i, p)) break; var tz = hz[p], tc = hc[p]; hz[p] = hz[i]; hc[p] = hc[i]; hz[i] = tz; hc[i] = tc; i = p; }
    }
    function pop() {
      var rz = hz[0], rc = hc[0]; n--;
      if (n > 0) {
        hz[0] = hz[n]; hc[0] = hc[n]; var i = 0;
        for (;;) { var l = 2 * i + 1, r = l + 1, m = i; if (l < n && less(l, m)) m = l; if (r < n && less(r, m)) m = r; if (m === i) break; var tz = hz[m], tc = hc[m]; hz[m] = hz[i]; hc[m] = hc[i]; hz[i] = tz; hc[i] = tc; i = m; }
      }
      return [rz, rc];
    }
    for (var c = 0; c < N; c++) if (bas[c] && (c % nx === nx - 1 || seed[c])) { done[c] = 1; push(z[c], c); }
    while (n > 0) {
      var t = pop(), zc = t[0], cc = t[1], i = cc % nx, j = (cc - i) / nx;
      var nb = [i > 0 ? cc - 1 : -1, i < nx - 1 ? cc + 1 : -1, j > 0 ? cc - nx : -1, j < ny - 1 ? cc + nx : -1];
      for (var m = 0; m < 4; m++) {
        var o = nb[m]; if (o < 0 || done[o] || !bas[o]) continue;
        done[o] = 1; if (z[o] < zc + eps) { z[o] = zc + eps; cnt++; }
        push(z[o], o);
      }
    }
    return cnt;
  }
  Sim.prototype.setupBasin = function (T) {
    var P = this.P, nx = this.nx, ny = this.ny, dx = this.dx, N = this.N, z = this.z, zfix = this.zfix, kind = this.kind, fl = this.fl;
    var DM = P.dem || {}, soilD = DM.soilD > 0 ? DM.soilD : 1, chR = Math.max(0, Math.round(((DM.chW > 0 ? DM.chW : dx) / dx - 1) / 2));
    var Ly = this.Ly, cs = T.cs, nr = T.rows.length, nc = T.rows[0].length, i, j, c, r, q;
    var bas = this.basin = new Uint8Array(N), ch = new Uint8Array(N), zTop = -1e30;
    this.tIdx = new Int32Array(N);
    // DEM の格子 → 計算格子：計算格子のセルに中心が入る DEM のセルを集める。谷筋のセルがあればその最低値（狭い谷を
    // ふさがない）、なければ流域内のセルの平均。流域かどうかは、谷筋を含むか、半分以上が流域か。
    // 入る DEM のセルがない（計算格子が DEM より細かい）ときは、最も近い DEM のセル（谷筋なら値そのもの、ほかは双一次補間）
    var xA = new Float64Array(N).fill(-1e30), nK = new Float64Array(N), nT = new Float64Array(N), nB = new Float64Array(N), sB = new Float64Array(N), sA = new Float64Array(N), nC = new Float64Array(N), mC = new Float64Array(N).fill(1e30);
    for (r = 0; r < nr; r++) {
      var jy = Math.floor((Ly - (r + 0.5) * cs) / dx); if (jy < 0 || jy >= ny) continue;
      var rowZ = T.rows[r], rowB = T.basin[r], rowC = T.chan ? T.chan[r] : null, rowK = T.corr ? T.corr[r] : null;
      for (q = 0; q < nc; q++) {
        var ix = Math.floor((q + 0.5) * cs / dx); if (ix < 0 || ix >= nx) continue;
        c = jy * nx + ix; var v = rowZ[q]; nT[c]++; sA[c] += v; if (v > xA[c]) xA[c] = v; if (rowK && rowK[q]) nK[c]++;
        if (rowB[q]) { nB[c]++; sB[c] += v; if (rowC && rowC[q]) { nC[c]++; if (v < mC[c]) mC[c] = v; } }
      }
    }
    for (j = 0; j < ny; j++) for (i = 0; i < nx; i++) {
      c = j * nx + i;
      r = Math.min(nr - 1, Math.max(0, Math.floor((Ly - (j + 0.5) * dx) / cs))); q = Math.min(nc - 1, Math.max(0, Math.floor((i + 0.5) * dx / cs)));
      this.tIdx[c] = r * nc + q;
      if (nT[c] > 0) {
        ch[c] = nC[c] > 0 ? 1 : 0; bas[c] = (nC[c] > 0 || 2 * nB[c] >= nT[c]) ? 1 : 0;
        // 流域のセル：谷筋の最低値か流域内の平均。流域の外のセル：流域の外の DEM のセルの平均（谷底の値を混ぜない）
        z[c] = bas[c] ? (nC[c] > 0 ? mC[c] : sB[c] / nB[c]) : (nT[c] > nB[c] ? (sA[c] - sB[c]) / (nT[c] - nB[c]) : sA[c] / nT[c]);
      } else {
        bas[c] = T.basin[r][q] ? 1 : 0; ch[c] = T.chan && T.chan[r][q] ? 1 : 0;
        z[c] = ch[c] ? T.rows[r][q] : sampleTerrain(T, (i + 0.5) * dx, (j + 0.5) * dx, Ly);
      }
      if (bas[c] && z[c] > zTop) zTop = z[c];
    }
    // 出口の先の水路のセル：自由流出（P.dem.outlet = 'sink'、既定）ならここで抜く
    var sinkL = [];
    var isSink = new Uint8Array(N);
    if (T.corr && DM.outlet !== 'channel') for (c = 0; c < N; c++) {
      var isK = nT[c] > 0 ? nK[c] > 0 : !!T.corr[Math.min(nr - 1, Math.max(0, Math.floor((Ly - (((c - c % nx) / nx) + 0.5) * dx) / cs)))][Math.min(nc - 1, Math.max(0, Math.floor(((c % nx) + 0.5) * dx / cs)))];
      if (isK && bas[c]) { sinkL.push(c); isSink[c] = 1; }
    }
    this.sinkCells = sinkL;
    // 計算格子での窪地埋め（4 近傍の Priority-Flood。右端と出口の先から）。DEM を計算格子に移すときにできる窪地をなくす
    if (DM.fill !== false) this.nFill = fill4(z, bas, isSink, nx, ny, 1e-3);
    var chD = new Uint8Array(N);
    for (c = 0; c < N; c++) if (ch[c]) {
      var ci = c % nx, cj = (c - ci) / nx;
      for (var dj = -chR; dj <= chR; dj++) for (var di = -chR; di <= chR; di++) {
        var ii = ci + di, jj = cj + dj; if (ii < 0 || jj < 0 || ii >= nx || jj >= ny) continue;
        var o = jj * nx + ii; if (bas[o]) chD[o] = 1;
      }
    }
    var wallZ = zTop + 30, nf;
    this.zOut = new Float32Array(N);   // 流域の外の元の標高（表示用。計算では高い壁）
    for (i = 0; i < nx; i++) {
      nf = 0;
      for (j = 0; j < ny; j++) {
        c = j * nx + i;
        if (!bas[c]) { this.zOut[c] = z[c]; kind[c] = 1; fl[c] = 0; z[c] = wallZ; zfix[c] = wallZ; }
        else if (chD[c]) { kind[c] = 0; fl[c] = 1; nf++; zfix[c] = (i < 2 || i >= nx - 2) ? z[c] : z[c] - P.erodDepth; }
        else { kind[c] = 3; fl[c] = 0; zfix[c] = z[c] - soilD; }
      }
      this.Bx[i] = Math.max(1, nf) * dx;
    }
    for (var q2 = 0; q2 < sinkL.length; q2++) { c = sinkL[q2]; zfix[c] = z[c]; }   // 抜く所は侵食しない
    // 流域界の外側の 1 列（4 近傍に流域のセルがある流域の外のセル）：元の地形の高さに戻し、入った水と土砂を抜く。
    // 流域界は分水界なので、外向きに流れる水（流域界を越える水）だけが出ていく（P.dem.edgeOut、既定 true）
    this.nEdge = 0;
    if (DM.edgeOut !== false) for (j = 0; j < ny; j++) for (i = 0; i < nx; i++) {
      c = j * nx + i; if (bas[c] || isSink[c]) continue;
      if (!((i > 0 && bas[c - 1]) || (i < nx - 1 && bas[c + 1]) || (j > 0 && bas[c - nx]) || (j < ny - 1 && bas[c + nx]))) continue;
      // 高さは、そのセルに入る DEM のセルの最高値（細い尾根＝分水界を平均で消さない。外へ出るのは分水界を越える水だけ）
      var zo_ = nT[c] > 0 ? Math.max(this.zOut[c], xA[c]) : this.zOut[c]; if (!(zo_ === zo_) || zo_ >= wallZ) continue;
      z[c] = zo_; zfix[c] = zo_; sinkL.push(c); this.nEdge++;
    }
    if (sinkL.length) { this.sinkMask = new Uint8Array(N); for (q2 = 0; q2 < sinkL.length; q2++) this.sinkMask[sinkL[q2]] = 1; }
    this.demMode = true; this.zTop = zTop; this.wallZ = wallZ;
    var Io = T.Iout > 0 ? T.Iout : 0.01; this.I = Io;
    for (i = 0; i < nx; i++) { this.Ix[i] = Io; this.colmin[i] = 1e30; for (j = 0; j < ny; j++) { c = j * nx + i; if (bas[c] && z[c] < this.colmin[i]) this.colmin[i] = z[c]; } if (this.colmin[i] > 1e29) this.colmin[i] = wallZ; }
  };

  function cellAtXY(sim, x, y) {
    var i = Math.floor(x / sim.dx), j = Math.floor(y / sim.dx);
    if (i < 0 || j < 0 || i >= sim.nx || j >= sim.ny) return -1;
    return j * sim.nx + i;
  }
  // 壁のセルなら、近くの壁でないセルへずらす
  function snapCell(sim, c) {
    if (c < 0 || sim.kind[c] !== 1 || (sim.basin && sim.basin[c] && sim.kind[c] !== 1)) return c;
    var nx = sim.nx, ny = sim.ny, ci = c % nx, cj = (c - ci) / nx, best = -1, bd = 1e9;
    for (var dj = -3; dj <= 3; dj++) for (var di = -3; di <= 3; di++) {
      var ii = ci + di, jj = cj + dj; if (ii < 0 || jj < 0 || ii >= nx || jj >= ny) continue;
      var o = jj * nx + ii; if (sim.kind[o] === 1 || sim.kind[o] === 2) continue;
      var dd = di * di + dj * dj; if (dd < bd) { bd = dd; best = o; }
    }
    return best >= 0 ? best : c;
  }

  // 点の流入：谷頭などの供給点（それぞれのハイドロ）と、雨の横流入（谷筋のセルごとの流量の表）
  Sim.prototype.setupPoints = function () {
    var P = this.P, N = this.N, self = this, k;
    this.pin = [];
    (P.sources || []).forEach(function (S) {
      var c = snapCell(self, cellAtXY(self, S.x, S.y)); if (c < 0) return;
      var ux = isFinite(S.ux) ? S.ux : 1, uy = isFinite(S.uy) ? S.uy : 0, L = Math.sqrt(ux * ux + uy * uy) || 1;
      self.pin.push({ c: c, ux: ux / L, uy: uy / L, I: Math.max(1e-3, S.I || 0.05), H: S.hydro || null, tab: null, Q: 0, Qs: 0, v: 0, name: S.name || '' });
    });
    var LI = P.latIn;
    if (LI && LI.pts && LI.Q) for (k = 0; k < LI.pts.length; k++) {
      var pt = LI.pts[k], c2 = snapCell(self, cellAtXY(self, pt[0], pt[1])); if (c2 < 0) continue;
      self.pin.push({ c: c2, ux: pt[2], uy: pt[3], I: Math.max(1e-3, pt[4] || 0.05), H: null, tab: LI.Q[k], dtab: LI.dt > 0 ? LI.dt : 60, Q: 0, Qs: 0, v: 0, lat: true });
    }
    this.pinIn = new Uint8Array(N); this.pH = new Float64Array(N); this.pMx = new Float64Array(N); this.pMy = new Float64Array(N);
    for (k = 0; k < this.pin.length; k++) this.pinIn[this.pin[k].c] = 1;
    // 格子に降る雨
    var R = P.rain || {};
    this.rainOn = !!(R.on && R.method === 'grid' && this.demMode);
    if (this.rainOn) { this.rS = new Float64Array(N); this.rR = new Float64Array(N); }
  };
  Sim.prototype.pinStep = function () {
    var P = this.P, t = this.t, dx = this.dx, PIN = this.pin, sum = 0, dbIn = P.transport === 'debris', Cmax = 0.9 * (1 - P.porosity), nman = P.manning, m, p;
    for (m = 0; m < PIN.length; m++) { var c0 = PIN[m].c; this.pH[c0] = 0; this.pMx[c0] = 0; this.pMy[c0] = 0; }
    for (m = 0; m < PIN.length; m++) {
      p = PIN[m]; var Q, Qs;
      if (p.tab) {
        var f = t / p.dtab, k0 = Math.floor(f), w = f - k0, nT = p.tab.length;
        Q = k0 >= nT - 1 ? p.tab[nT - 1] : p.tab[k0] * (1 - w) + p.tab[k0 + 1] * w; Qs = 0;
      } else { Q = QatH(p.H, t); Qs = QsAtH(p.H, t); }
      if (!(Q > 0)) { p.Q = 0; p.Qs = 0; p.v = 0; continue; }
      if (dbIn && Qs > Cmax * Q) Qs = Cmax * Q;
      p.Q = Q; p.Qs = Qs > 0 ? Qs : 0; sum += Q;
      var q = Q / dx;
      p.v = dbIn ? q / this.debrisDepth(q, p.Qs / Q, p.I) : q / Math.pow(q * nman / Math.sqrt(p.I), 0.6);
      var a = Q / (dx * dx);
      this.pH[p.c] += a; this.pMx[p.c] += a * p.v * p.ux; this.pMy[p.c] += a * p.v * p.uy;
    }
    return sum;
  };
  // 供給点の土砂（粒度は「供給」の割合）
  Sim.prototype.pinFeed = function (dtm, bulk) {
    var PIN = this.pin; if (!PIN || !PIN.length) return;
    var K = this.K, ff = this.ff, A = this.A, Cc = this.Cc, dx = this.dx, lam = 1 - this.P.porosity, add = 0;
    for (var m = 0; m < PIN.length; m++) {
      var p = PIN[m]; if (!(p.Qs > 0)) continue;
      var per = p.Qs * dtm / (dx * dx), c = p.c, ck = c * K, k;
      add += p.Qs * dtm;
      if (bulk) { for (k = 0; k < K; k++) A[ck + k] += per / lam * ff[k]; this.z[c] += per / lam; this.dzs[c] += per / lam; this.fixLayers(c); }
      else for (k = 0; k < K; k++) Cc[ck + k] += per * ff[k];
    }
    this.sedIn += add; this.QsIn += add / Math.max(1e-9, dtm);
  };

  // ---- 雨（格子）
  function rainUni(tab, t) {
    if (!tab || !tab.length) return 0;
    var tm = t / 60, v = 0, i;
    for (i = 0; i < tab.length; i++) { if (tm >= tab[i][0] - 1e-9) v = tab[i][1]; else break; }
    if (tab.length > 1 && tm >= tab[tab.length - 1][0] + (tab[tab.length - 1][0] - tab[tab.length - 2][0])) v = 0;
    return Math.max(0, v);
  }
  function rainFrame(frames, t) {
    if (!frames || !frames.length) return null;
    var tm = t / 60, f = null, i;
    for (i = 0; i < frames.length; i++) { if (tm >= frames[i].t - 1e-9) f = frames[i]; else break; }
    if (f && frames.length > 1) { var dtf = frames[frames.length - 1].t - frames[frames.length - 2].t; if (tm >= frames[frames.length - 1].t + dtf) f = null; }
    return f;
  }
  // 損失（有効雨量）：流出係数 f、一次流出率 f1・飽和雨量 Rsa・初期損失 R0、浸透能 ks と土層の空き容量
  function lossCell(L, Sa, Ra, c, r, dt) {
    var dth = dt / 3600;
    if (L.type === 'f1rsa') {
      Ra[c] += r * dth;
      var R0 = L.R0 || 0, f = Ra[c] < R0 ? 0 : (Ra[c] < R0 + (L.Rsa || 0) ? (L.f1 != null ? L.f1 : 0.5) : 1);
      return f * r;
    }
    if (L.type === 'ks') {
      var Smax = Math.max(0, (L.D || 1) * ((L.n || 0.45) - (L.th0 || 0.25)) * 1000), ks = (L.ks || 1e-5) * 3.6e6;
      if (Sa[c] >= Smax) return r;
      var inf = Math.min(r, ks), add = inf * dth;
      if (Sa[c] + add > Smax) { var over = (Sa[c] + add - Smax) / dth; Sa[c] = Smax; return r - inf + over; }
      Sa[c] += add; return r - inf;
    }
    return (L.f != null ? L.f : 0.8) * r;
  }
  Sim.prototype.rainAtCell = function (c, t) {
    var R = this.P.rain || {};
    if (R.mode === 'grid' && R.frames && R.frames.length) { var fr = rainFrame(R.frames, t); return fr && this.tIdx ? Math.max(0, fr.r[this.tIdx[c]] || 0) : 0; }
    return rainUni(R.table, t);
  };
  Sim.prototype.rainStep = function (dt) {
    var P = this.P, R = P.rain, N = this.N, bas = this.basin, t = this.t, dx = this.dx, a = dx * dx;
    var rr = this._rainR || (this._rainR = new Float64Array(N)), L = R.loss || { type: 'f', f: 1 }, gwOn = !!(P.gw && P.gw.on);
    var fr = (R.mode === 'grid' && R.frames && R.frames.length) ? rainFrame(R.frames, t) : null, ru = fr ? 0 : rainUni(R.table, t);
    var vol = 0, tot = 0;
    for (var c = 0; c < N; c++) {
      if (!bas[c]) { rr[c] = 0; continue; }
      var r = fr ? Math.max(0, fr.r[this.tIdx[c]] || 0) : ru;
      if (!(r > 0)) { rr[c] = 0; continue; }
      var e = gwOn ? r : lossCell(L, this.rS, this.rR, c, r, dt);   // 伏流ありなら全部を地表へ（浸透は伏流の計算で）
      rr[c] = e / 3.6e6 * dt; vol += rr[c]; tot += r / 3.6e6 * dt;
    }
    this._rainV = vol * a; this.rainVol += tot * a; this.rainLoss += (tot - vol) * a;
    return rr;
  };

  // ======================================================================================
  // 斜面崩壊・天然ダム
  // ======================================================================================
  // 崩壊の範囲：中心 (x, y)、半径 R。表層：一様な深さ D（無限長斜面）。深層：深さ D のさじ形のすべり面（中心で D、縁で 0）
  // 点が多角形の中か（偶奇則）
  function inPoly(P, x, y) {
    var n = P.length, a = false;
    for (var i = 0, j = n - 1; i < n; j = i++) {
      var xi = P[i][0], yi = P[i][1], xj = P[j][0], yj = P[j][1];
      if ((yi > y) !== (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi) a = !a;
    }
    return a;
  }
  // 点から多角形の辺までの最短距離
  function edgeDist(P, x, y) {
    var n = P.length, m = 1e30;
    for (var i = 0, j = n - 1; i < n; j = i++) {
      var ax = P[j][0], ay = P[j][1], bx = P[i][0] - ax, by = P[i][1] - ay, L2 = bx * bx + by * by;
      var t = L2 > 0 ? Math.max(0, Math.min(1, ((x - ax) * bx + (y - ay) * by) / L2)) : 0;
      var ex = ax + t * bx - x, ey = ay + t * by - y; m = Math.min(m, Math.sqrt(ex * ex + ey * ey));
    }
    return m;
  }
  Sim.prototype.setupSlides = function () {
    var P = this.P, N = this.N, dx = this.dx, nx = this.nx, ny = this.ny, z = this.z, kind = this.kind, self = this;
    this.slides = []; this.dep = new Uint8Array(N); this.depOn = false; this.ldamList = []; this.depK = 0; this._chkT = 0;
    (P.slides || []).forEach(function (S, q) {
      var R = Math.max(dx, S.R || 20), D = Math.max(0.1, S.D || 2), deep = S.type === 'deep', cells = [], zs = [];
      var poly = S.poly && S.poly.length >= 3 ? S.poly : null, i0, i1, j0, j1, i, j, c;
      var zt = -1e30, sx = 0, sy = 0;
      if (poly) {
        // 多角形：中の格子を選ぶ。深さは表層なら一様、深層なら縁からの距離でスプーン形
        var xa = 1e30, xb = -1e30, ya = 1e30, yb = -1e30, pc = [], de = [], dmx = 0;
        poly.forEach(function (v) { xa = Math.min(xa, v[0]); xb = Math.max(xb, v[0]); ya = Math.min(ya, v[1]); yb = Math.max(yb, v[1]); });
        i0 = Math.max(0, Math.floor(xa / dx)); i1 = Math.min(nx - 1, Math.floor(xb / dx));
        j0 = Math.max(0, Math.floor(ya / dx)); j1 = Math.min(ny - 1, Math.floor(yb / dx));
        for (j = j0; j <= j1; j++) for (i = i0; i <= i1; i++) {
          c = j * nx + i; if (kind[c] === 2 || (self.basin && !self.basin[c])) continue;
          var px = (i + 0.5) * dx, py = (j + 0.5) * dx;
          if (!inPoly(poly, px, py)) continue;
          var d = edgeDist(poly, px, py); pc.push(c); de.push(d); if (d > dmx) dmx = d;
        }
        dmx = Math.max(dmx, 0.5 * dx);
        for (var k = 0; k < pc.length; k++) {
          var u = Math.min(1, de[k] / dmx), th = deep ? D * Math.sqrt(Math.max(0, 1 - (1 - u) * (1 - u))) : D;
          if (th < 0.05) th = 0.05;
          c = pc[k]; cells.push(c); zs.push(z[c] - th); if (z[c] > zt) zt = z[c]; sx += c % nx; sy += (c - c % nx) / nx;
        }
        R = Math.max(dx, Math.sqrt(cells.length * dx * dx / Math.PI));   // 等価半径（到達距離に使う）
      } else {
        i0 = Math.max(0, Math.floor((S.x - R) / dx)); i1 = Math.min(nx - 1, Math.floor((S.x + R) / dx));
        j0 = Math.max(0, Math.floor((S.y - R) / dx)); j1 = Math.min(ny - 1, Math.floor((S.y + R) / dx));
        for (j = j0; j <= j1; j++) for (i = i0; i <= i1; i++) {
          c = j * nx + i; if (kind[c] === 2 || (self.basin && !self.basin[c])) continue;
          var rx = (i + 0.5) * dx - S.x, ry = (j + 0.5) * dx - S.y, rr = Math.sqrt(rx * rx + ry * ry) / R;
          if (rr > 1) continue;
          var th2 = deep ? D * Math.sqrt(Math.max(0, 1 - rr * rr)) : D;
          if (th2 < 0.05) continue;
          cells.push(c); zs.push(z[c] - th2); if (z[c] > zt) zt = z[c]; sx += i; sy += j;
        }
      }
      if (!cells.length) return;
      var ne = S.ne > 0 ? S.ne : (deep ? 0.05 : 0.15), m0 = S.m0 >= 0 ? S.m0 : 0.2;
      self.slides.push({ id: q, S: S, cells: cells, zslip: Float64Array.from(zs), deep: deep, D: D, R: R, ne: ne, fired: false, done: new Uint8Array(cells.length),
        Sw: m0 * D * ne * 1000, Fs: NaN, zTop: zt, ci: Math.round(sx / cells.length), cj: Math.round(sy / cells.length), m: m0 });
    });
    (P.ldams || []).forEach(function (L) { self.placeLdam(L); });
  };

  // 表面の勾配 tanβ（中心差分、端は片側）
  Sim.prototype.surfSlope = function (c) {
    var nx = this.nx, ny = this.ny, z = this.z, dx = this.dx, i = c % nx, j = (c - i) / nx;
    var a = i > 0 ? c - 1 : c, b = i < nx - 1 ? c + 1 : c, s = j > 0 ? c - nx : c, n = j < ny - 1 ? c + nx : c;
    var gx = (z[b] - z[a]) / ((b - a) * dx || dx), gy = (z[n] - z[s]) / (((n - s) / nx) * dx || dx);
    if (b === a) gx = 0; if (n === s) gy = 0;
    return Math.sqrt(gx * gx + gy * gy);
  };
  // すべり面より上の飽和厚 hw [m]（伏流ありで表層なら土層の地下水、それ以外はタンク）
  Sim.prototype.slideHw = function (sl, k) {
    var P = this.P, c = sl.cells[k], T = this.z[c] - sl.zslip[k];
    if (!sl.deep && P.gw && P.gw.on && this.por[c]) {
      var n = P.porosity, b = this.W[c] / n, base = this.zfix[c];   // 地下水面の高さ = zfix + W/n
      return Math.max(0, Math.min(T, base + b - sl.zslip[k]));
    }
    return Math.max(0, Math.min(T, sl.Sw / 1000 / sl.ne));
  };
  // 安全率。表層：無限長斜面 Fs = [c' + (γt(D−hw) + γ'hw) cos²β tanφ] / [(γt(D−hw) + γsat hw) sinβ cosβ]（セルごと）
  // 深層：Hovland（1977）の 3 次元簡便法（柱ごとの重さ・間隙水圧、すべり面の傾きは移動方向に測る）
  Sim.prototype.slideFs = function (sl) {
    var S = sl.S, gt = S.gt > 0 ? S.gt : 18, gs = S.gs > 0 ? S.gs : 20, gw = 9.81, coh = S.c >= 0 ? S.c : (sl.deep ? 30 : 5), tp = Math.tan((S.phi > 0 ? S.phi : (sl.deep ? 30 : 35)) * Math.PI / 180);
    var z = this.z, dx = this.dx, nx = this.nx, k, c, n = sl.cells.length;
    if (!sl.deep) {
      var fsc = sl.fsc || (sl.fsc = new Float64Array(n)), mn = 1e9, msum = 0, cnt = 0;
      for (k = 0; k < n; k++) {
        if (sl.done[k]) { fsc[k] = 9; continue; }
        c = sl.cells[k]; var D = z[c] - sl.zslip[k]; if (D <= 0.01) { fsc[k] = 9; continue; }
        var hw = this.slideHw(sl, k), tb = this.surfSlope(c), b = Math.atan(tb), cb = Math.cos(b), sb = Math.sin(b);
        var res = coh + (gt * (D - hw) + (gs - gw) * hw) * cb * cb * tp, drv = (gt * (D - hw) + gs * hw) * sb * cb;
        fsc[k] = drv > 1e-9 ? res / drv : 9; if (fsc[k] < mn) mn = fsc[k]; msum += hw / D; cnt++;
      }
      sl.m = cnt ? msum / cnt : 0; return mn;
    }
    // 移動方向：崩壊域の表面勾配の平均の下り向き
    var gxs = 0, gys = 0, ny = this.ny;
    for (k = 0; k < n; k++) {
      c = sl.cells[k]; var i = c % nx, j = (c - i) / nx;
      if (i > 0 && i < nx - 1) gxs += (z[c + 1] - z[c - 1]) / (2 * dx);
      if (j > 0 && j < ny - 1) gys += (z[c + nx] - z[c - nx]) / (2 * dx);
    }
    var L = Math.sqrt(gxs * gxs + gys * gys); if (L < 1e-9) return 9;
    var ux = -gxs / L, uy = -gys / L, num = 0, den = 0, a1 = dx * dx, ms = 0;
    var zsl = this._zsl || (this._zsl = new Float64Array(this.N)), msk = this._zm || (this._zm = new Uint8Array(this.N));
    for (k = 0; k < n; k++) { zsl[sl.cells[k]] = sl.zslip[k]; msk[sl.cells[k]] = 1; }
    for (k = 0; k < n; k++) {
      c = sl.cells[k]; var T = z[c] - sl.zslip[k]; if (T <= 0.01) continue;
      var i2 = c % nx, j2 = (c - i2) / nx;
      var sxp = (i2 < nx - 1 && msk[c + 1]) ? zsl[c + 1] : z[Math.min(c + 1, this.N - 1)], sxm = (i2 > 0 && msk[c - 1]) ? zsl[c - 1] : z[Math.max(0, c - 1)];
      var syp = (j2 < ny - 1 && msk[c + nx]) ? zsl[c + nx] : (j2 < ny - 1 ? z[c + nx] : zsl[c]), sym = (j2 > 0 && msk[c - nx]) ? zsl[c - nx] : (j2 > 0 ? z[c - nx] : zsl[c]);
      var al = Math.atan(-((sxp - sxm) / (2 * dx) * ux + (syp - sym) / (2 * dx) * uy));   // すべり面の移動方向の傾き
      var hw2 = this.slideHw(sl, k), Wt = (gt * (T - hw2) + gs * hw2) * a1, Ab = a1 / Math.max(0.2, Math.cos(al)), u = gw * hw2 * Math.cos(al) * Math.cos(al);
      num += coh * Ab + Math.max(0, Wt * Math.cos(al) - u * Ab) * tp; den += Wt * Math.sin(al); ms += hw2 / T;
    }
    for (k = 0; k < n; k++) msk[sl.cells[k]] = 0;
    sl.m = n ? ms / n : 0; sl.ux = ux; sl.uy = uy;
    return den > 1e-9 ? num / den : 9;
  };

  Sim.prototype.slideStep = function (dt) {
    var P = this.P, t = this.t, m;
    for (m = 0; m < this.slides.length; m++) {
      var sl = this.slides[m], S = sl.S; if (sl.fired) continue;
      // 地下水のタンク：dS/dt = r − S/Tk（S：mm）
      var Tk = (S.Tk > 0 ? S.Tk : (sl.deep ? 48 : 6)) * 3600, r = P.rain && P.rain.on ? this.rainAtCell(sl.cells[0], t) : 0;
      sl.Sw += (r / 3600 - sl.Sw / Tk) * dt;
      if (sl.Sw > sl.D * sl.ne * 1000) sl.Sw = sl.D * sl.ne * 1000;
    }
    this._chkT += dt;
    var chk = P.slideChk > 0 ? P.slideChk : 1;
    if (this._chkT < chk) return;
    this._chkT = 0;
    for (m = 0; m < this.slides.length; m++) {
      var s2 = this.slides[m], S2 = s2.S; if (s2.fired) continue;
      var Fs = this.slideFs(s2); s2.Fs = Fs;
      if (S2.trigger === 'time') { if (t >= (S2.t || 0)) this.releaseSlide(s2, null); }
      else if (s2.deep) { if (Fs < 1) this.releaseSlide(s2, null); }
      else if (Fs < 1) {   // 表層：安全率が 1 を下回ったセルだけ崩れる
        var list = [];
        for (var k = 0; k < s2.cells.length; k++) if (!s2.done[k] && s2.fsc[k] < 1) list.push(k);
        if (list.length) this.releaseSlide(s2, list);
      }
    }
    if (this.depOn && P.gw && P.gw.on) this.progressive();
    this.ldamSlip();
  };

  // 層から厚さ dz（かさ）を取り去る。取った粒径別の実質体積を out に足す
  Sim.prototype.takeSed = function (c, dz, out) {
    var K = this.K, A = this.A, ck = c * K, left = dz, it = 0, k, lam = 1 - this.P.porosity;
    while (left > 1e-9 && it++ < 400) {
      var av = 0; for (k = 0; k < K; k++) av += A[ck + k];
      if (av < 1e-9) break;
      var take = left < av ? left : av, frc = take / av;
      for (k = 0; k < K; k++) { var tt = A[ck + k] * frc; A[ck + k] -= tt; if (out) out[k] += tt * lam; }
      left -= take; this.fixLayers(c);
    }
    return dz - left;
  };
  // 堆積物（かさ dz、粒度 fr）を足す
  Sim.prototype.addDeposit = function (c, dz, fr) {
    var K = this.K, A = this.A, ck = c * K;
    if (this.kind[c] !== 0) { this.kind[c] = 0; this.zfix[c] = this.z[c]; this.fl[c] = 1; if (this.fixedC) this.fixedC[c] = 0; }
    for (var k = 0; k < K; k++) A[ck + k] += dz * fr[k];
    this.z[c] += dz; this.dzs[c] += dz; this.fixLayers(c);
    this.por[c] = 1; this.dep[c] = 1; this.depOn = true;
  };

  // 崩壊：list（セルの番号の配列、null で全部）を崩す
  Sim.prototype.releaseSlide = function (sl, list) {
    var P = this.P, S = sl.S, n = P.porosity, Cst = 1 - n, K = this.K, dx = this.dx, a1 = dx * dx, z = this.z, h = this.h, Cc = this.Cc, f0 = this.f0;
    var gwOn = !!(P.gw && P.gw.on), idx = list || sl.cells.map(function (_, k) { return k; }), k, q;
    this.slideFs(sl);
    var mode = S.move || 'auto';
    if (mode === 'auto') mode = (sl.m >= (S.mc > 0 ? S.mc : 0.8) && P.transport === 'debris') ? 'fluid' : 'block';
    if (mode === 'fluid' && P.transport !== 'debris') mode = 'block';
    var CL = S.CL > 0 ? S.CL : 0.5, Sr0 = S.Sr0 >= 0 ? S.Sr0 : 0.5, vol = 0, wat = 0, sed = new Float64Array(K), added = 0, sedTot = 0;
    var frc = new Float64Array(K);
    for (q = 0; q < idx.length; q++) {
      k = idx[q]; if (sl.done[k]) continue; sl.done[k] = 1;
      var c = sl.cells[k], T = z[c] - sl.zslip[k]; if (T <= 1e-6) continue;
      var hw = this.slideHw(sl, k), wsl;
      if (gwOn && this.por[c] && !sl.deep) { wsl = Math.max(0, Math.min(this.W[c], this.W[c] - n * Math.max(0, sl.zslip[k] - this.zfix[c]))); this.W[c] -= wsl; }
      else if (gwOn) { wsl = n * (hw + Sr0 * Math.max(0, T - hw)); this.Vwin += wsl * a1; }   // 追っていなかった間隙の水を足す
      else wsl = mode === 'fluid' ? n * T : 0;   // 伏流なし：河床の間隙はいつも水で満ちている扱い（かさの中の水）
      // すべり面より上の土：河床の層があれば層から、なければ「河床材料」の粒度
      for (var kk = 0; kk < K; kk++) frc[kk] = 0;
      var got = 0;
      if (this.kind[c] === 0) got = this.takeSed(c, Math.min(T, Math.max(0, z[c] - this.zfix[c])), frc);
      for (kk = 0; kk < K; kk++) frc[kk] += (T - got) * Cst * f0[kk];
      z[c] = sl.zslip[k]; this.dzs[c] -= T;
      if (this.zfix[c] > z[c]) this.zfix[c] = z[c];
      if (this.kind[c] === 3) { this.zfix[c] = z[c]; }    // 崩れた跡は土層なし
      if (gwOn && this.por[c]) { var cap = n * Math.max(0, z[c] - this.zfix[c]); if (this.W[c] > cap) { h[c] += this.W[c] - cap; this.W[c] = cap; } }
      if (mode === 'fluid') {
        var ss = 0; for (kk = 0; kk < K; kk++) { Cc[c * K + kk] += frc[kk]; ss += frc[kk]; }
        h[c] += ss + wsl;
        var Sx = 0; for (kk = 0; kk < K; kk++) Sx += Cc[c * K + kk];
        if (Sx > CL * h[c]) { var ad = Sx / CL - h[c]; h[c] += ad; added += ad * a1; }   // 流動化する濃度まで水を足す（法面の流れ・湧水を取り込むと仮定）
        sedTot += ss * a1;
      } else { for (kk = 0; kk < K; kk++) sed[kk] += frc[kk] * a1; vol += T * a1; wat += wsl * a1; }
    }
    this.Vwin += added;
    if (mode === 'block' && vol > 0) {
      var ss2 = 0; for (kk = 0; kk < K; kk++) ss2 += sed[kk];
      for (kk = 0; kk < K; kk++) frc[kk] = ss2 > 0 ? sed[kk] / ss2 : f0[kk];
      this.runout(sl, ss2 / Cst, wat, frc);
    }
    if (!list || sl.done.every(function (v) { return v; })) sl.fired = true;
    this.slideLog.push({ t: this.t, id: sl.id, mode: mode, vol: mode === 'block' ? vol : sedTot / Cst, Fs: sl.Fs, m: sl.m, water: added });
  };

  // 土塊の瞬時移動：最急降下の経路を、エネルギー線（崩壊域の最高点から勾配 μ）が地表に届くまでたどり、そこに堆積させる。
  // μ：指定（等価摩擦角）か Scheidegger（1973）log10(H/L) = −0.15666 log10 V + 0.62419（V：m³）
  Sim.prototype.runout = function (sl, V, Wv, frc) {
    var S = sl.S, nx = this.nx, ny = this.ny, dx = this.dx, z = this.z, kind = this.kind;
    var mu = S.mu > 0 ? Math.tan(S.mu * Math.PI / 180) : Math.max(0.05, Math.min(0.9, Math.pow(10, -0.15666 * Math.log(Math.max(1, V)) / Math.LN10 + 0.62419)));
    var c = sl.cj * nx + sl.ci, L = 0, it = 0, top = sl.zTop;
    var DCX = [1, 1, 0, -1, -1, -1, 0, 1], DCY = [0, 1, 1, 1, 0, -1, -1, -1];
    while (it++ < 20000) {
      var ci = c % nx, cj = (c - ci) / nx, best = -1, bs = 0, bl = 0;
      for (var q = 0; q < 8; q++) {
        var ii = ci + DCX[q], jj = cj + DCY[q]; if (ii < 0 || jj < 0 || ii >= nx || jj >= ny) continue;
        var o = jj * nx + ii; if (kind[o] === 2 || (this.basin && !this.basin[o])) continue;
        var dl = (q % 2 ? Math.SQRT2 : 1) * dx, s = (z[c] - z[o]) / dl;
        if (s > bs) { bs = s; best = o; bl = dl; }
      }
      if (best < 0) break;
      L += bl; c = best;
      if (z[c] >= top - mu * (L + sl.R)) break;
    }
    sl.stop = c; sl.mu = mu; sl.L = L;
    this.depositCone(c, V, Math.tan((S.thDep > 0 ? S.thDep : 10) * Math.PI / 180), Wv, frc);
  };
  // 円錐形の堆積（頂点はセル c0、表面の勾配 tanT）。体積 V（かさ）に合う高さを二分法で求める
  Sim.prototype.depositCone = function (c0, V, tanT, Wv, frc) {
    var nx = this.nx, ny = this.ny, dx = this.dx, z = this.z, kind = this.kind, a1 = dx * dx, self = this;
    var ci = c0 % nx, cj = (c0 - ci) / nx, z0 = z[c0];
    function vol(H, apply) {
      var R = H / tanT, rr = Math.ceil(R / dx) + 1, s = 0, cells = apply ? [] : null;
      for (var j = Math.max(0, cj - rr); j <= Math.min(ny - 1, cj + rr); j++) for (var i = Math.max(0, ci - rr); i <= Math.min(nx - 1, ci + rr); i++) {
        var c = j * nx + i; if (kind[c] === 2 || (self.basin && !self.basin[c])) continue;
        var r = Math.sqrt((i - ci) * (i - ci) + (j - cj) * (j - cj)) * dx, zt = z0 + H - tanT * r;
        if (zt > z[c]) { s += (zt - z[c]) * a1; if (apply) cells.push([c, zt - z[c]]); }
      }
      return apply ? cells : s;
    }
    var lo = 0, hi = Math.max(1, Math.cbrt(V)), it;
    for (it = 0; it < 60 && vol(hi) < V; it++) hi *= 2;
    for (it = 0; it < 60; it++) { var mid = 0.5 * (lo + hi); if (vol(mid) < V) lo = mid; else hi = mid; }
    var cells = vol(hi, true), tot = 0, q;
    for (q = 0; q < cells.length; q++) tot += cells[q][1];
    var sc = tot > 0 ? V / (tot * a1) : 0;   // 体積をきっちり合わせる
    this.depositCells(cells, sc, Wv, frc);
  };
  Sim.prototype.depositCells = function (cells, sc, Wv, frc) {
    var P = this.P, n = P.porosity, a1 = this.dx * this.dx, gwOn = !!(P.gw && P.gw.on), tot = 0, q, lam = 1 - n;
    var fb = new Float64Array(this.K); for (var k = 0; k < this.K; k++) fb[k] = frc[k];   // 粒度（かさの割合）
    for (q = 0; q < cells.length; q++) tot += cells[q][1] * sc;
    for (q = 0; q < cells.length; q++) {
      var c = cells[q][0], t = cells[q][1] * sc; if (!(t > 0)) continue;
      this.addDeposit(c, t, fb);
      var w = tot > 0 ? Wv / a1 * t / tot : 0;         // 運ばれた水を厚さに比例して配る
      if (gwOn) { var cap = n * Math.max(0, this.z[c] - this.zfix[c]), room = Math.max(0, cap - this.W[c]), ww = Math.min(w, room); this.W[c] += ww; this.h[c] += w - ww; }
      else this.h[c] += w;
      this.Vdep = (this.Vdep || 0) + t * a1;
    }
    void lam;
  };

  // 天然ダムを直接置く：中心 (x, y)、流れの向き (ux, uy)、高さ H、天端幅 Wc、上流・下流の法面勾配（度）、横の範囲 lat
  // DEM：堰堤の上流の範囲。堰堤の上流面のすぐ上のセルから、上り（または平ら）の 4 近傍をたどって広げる（谷の集水域）。
  // 上流の堰堤から先に決め、決まったセルは取らない（いちばん近い下流の堰堤の範囲になる）
  Sim.prototype.damRegionDEM = function () {
    var nx = this.nx, ny = this.ny, N = this.N, dx = this.dx, z = this.z, kind = this.kind, nd = this.dams.length, region = this.region, bas = this.basin;
    region.fill(nd);
    var ord = this.dams.map(function (D, q) { return q; }).sort(function (a, b) { return this.dams[b].crest - this.dams[a].crest; }.bind(this));
    var st = new Int32Array(N);
    for (var r = 0; r < nd; r++) {
      var q = ord[r], D = this.dams[q], n = 0, i, j, c;
      for (var k = 0; k < D.cells.length; k++) {
        c = D.cells[k]; var nb = [c - 1, c + 1, c - nx, c + nx];
        for (var m = 0; m < 4; m++) {
          var o = nb[m]; if (o < 0 || o >= N || kind[o] === 2 || region[o] !== nd || (bas && !bas[o])) continue;
          var px = ((o % nx) + 0.5) * dx - D.cx, py = (((o / nx) | 0) + 0.5) * dx - D.cy;
          if (px * D.tx + py * D.ty >= 0) continue;   // 上流側だけ
          region[o] = q; st[n++] = o;
        }
      }
      while (n > 0) {
        c = st[--n]; i = c % nx; j = (c - i) / nx;
        var nb2 = [i > 0 ? c - 1 : -1, i < nx - 1 ? c + 1 : -1, j > 0 ? c - nx : -1, j < ny - 1 ? c + nx : -1];
        for (var m2 = 0; m2 < 4; m2++) {
          var o2 = nb2[m2]; if (o2 < 0 || kind[o2] === 2 || region[o2] !== nd || (bas && !bas[o2])) continue;
          if (z[o2] < z[c] - 1e-6) continue;   // 下りには進まない
          region[o2] = q; st[n++] = o2;
        }
      }
    }
  };
  // DEM の出口：出口の先の水路に入った水と土砂をその場で抜き、流出として数える（自由流出）。
  // 水の量は混合物の体積（境界の流出と同じ数え方）、土砂は実質体積
  Sim.prototype.sinkStep = function (dt) {
    var L = this.sinkCells, h = this.h, Cc = this.Cc, K = this.K, a = this.dx * this.dx, wv = 0, sv = 0, z = this.z, z0 = this.z0, A = this.A, Cst = 1 - this.P.porosity;
    for (var q = 0; q < L.length; q++) {
      var c = L[q];
      // 崩落や堆積でたまった土砂も流出として抜き、もとの高さに戻す（河床の変化は空隙込みの体積で数えるので、水も同じ量を流出に）
      var B = z[c] - z0[c];
      if (B > 1e-9) {
        var As = 0, k0; for (k0 = 0; k0 < K; k0++) As += A[c * K + k0];
        for (k0 = 0; k0 < K; k0++) { this.sedOutK[k0] += B * a * Cst * (As > 0 ? A[c * K + k0] / As : 1 / K); A[c * K + k0] = 0; }
        this.nl[c] = 0; this.topTh[c] = 0; z[c] = z0[c];
        wv += B * a; sv += B * a * Cst;
      }
      if (!(h[c] > 0)) continue;
      wv += h[c] * a; h[c] = 0; this.qx[c] = 0; this.qy[c] = 0;
      for (var k = 0; k < K; k++) { var v = Cc[c * K + k]; if (v > 0) { sv += v * a; this.sedOutK[k] += v * a; Cc[c * K + k] = 0; } }
    }
    this.Vwout += wv; this.sedOut += sv;
    this.Qout += wv / dt; this.QsOut += sv / dt;
  };
  Sim.prototype.placeLdam = function (L) {
    var nx = this.nx, ny = this.ny, dx = this.dx, z = this.z, kind = this.kind, c0 = snapCell(this, cellAtXY(this, L.x, L.y));
    if (c0 < 0) return;
    var ux = isFinite(L.ux) ? L.ux : 1, uy = isFinite(L.uy) ? L.uy : 0, nn = Math.sqrt(ux * ux + uy * uy) || 1; ux /= nn; uy /= nn;
    var H = Math.max(0.5, L.H || 10), Wc = Math.max(0, L.Wc || 0), tu = Math.tan((L.up > 0 ? L.up : 30) * Math.PI / 180), td = Math.tan((L.dn > 0 ? L.dn : 25) * Math.PI / 180);
    var lat = L.lat > 0 ? L.lat : Math.max(50, 8 * H), x0 = (c0 % nx + 0.5) * dx, y0 = (Math.floor(c0 / nx) + 0.5) * dx, zc = z[c0];
    var ext = Math.max(H / tu, H / td) + Wc + lat, cells = [], f0 = this.f0;
    for (var j = Math.max(0, Math.floor((y0 - ext) / dx)); j <= Math.min(ny - 1, Math.floor((y0 + ext) / dx)); j++)
      for (var i = Math.max(0, Math.floor((x0 - ext) / dx)); i <= Math.min(nx - 1, Math.floor((x0 + ext) / dx)); i++) {
        var c = j * nx + i; if (kind[c] === 2 || (this.basin && !this.basin[c])) continue;
        var px = (i + 0.5) * dx - x0, py = (j + 0.5) * dx - y0, s = px * ux + py * uy, w = -px * uy + py * ux;
        if (Math.abs(w) > lat) continue;
        var zt = zc + H - (s > Wc / 2 ? (s - Wc / 2) * td : s < -Wc / 2 ? (-Wc / 2 - s) * tu : 0);
        if (zt > z[c] + 1e-3) cells.push([c, zt - z[c]]);
      }
    var LD = this.P.ldam || {}, fr = new Float64Array(this.K), K = this.K;
    for (var k = 0; k < K; k++) fr[k] = f0[k];
    if (LD.grains && LD.grains.length === K) { var sm = 0; for (k = 0; k < K; k++) sm += Math.max(0, LD.grains[k]); for (k = 0; k < K; k++) fr[k] = sm > 0 ? Math.max(0, LD.grains[k]) / sm : f0[k]; }
    var n = this.P.porosity, gwOn = !!(this.P.gw && this.P.gw.on), s0 = LD.s0 >= 0 ? LD.s0 : 0;
    var V = 0;
    for (var q = 0; q < cells.length; q++) {
      var cc = cells[q][0], t = cells[q][1];
      this.addDeposit(cc, t, fr); V += t * dx * dx;
      if (gwOn) this.W[cc] += n * t * s0;
    }
    if (gwOn) { var wv = 0; for (q = 0; q < cells.length; q++) wv += n * cells[q][1] * s0; this.W0add = (this.W0add || 0) + wv * dx * dx; }
    this.ldamList.push({ c0: c0, x0: x0, y0: y0, zc: zc, H: H, ux: ux, uy: uy, crest: zc + H, lat: lat, V: V, cells: cells.map(function (a) { return a[0]; }), slipped: false, L: L });
  };

  // 進行性崩壊：天然ダム（堆積物）のうち、地表まで飽和して湧水が出ている（浸潤面が地表に達した）セルで、
  // 浸透流が斜面に平行なときの無限長斜面 Fs = c'/(γsat T sinβ cosβ) + (γ'/γsat) tanφ/tanβ が 1 を下回れば、
  // 下の隣のセルへ、勾配が限界 tanβc = (γ'/γsat) tanφ になるまで土を移す（崩れた土は飽和したまま。流れが運ぶ）
  Sim.prototype.progressive = function () {
    var P = this.P, LD = P.ldam || {}, N = this.N, nx = this.nx, ny = this.ny, dx = this.dx, z = this.z, zfix = this.zfix, h = this.h, W = this.W, n = P.porosity;
    var gs = LD.gs > 0 ? LD.gs : 20, gw = 9.81, tp = Math.tan((LD.phi > 0 ? LD.phi : 35) * Math.PI / 180), coh = LD.c > 0 ? LD.c : 0, Ts = LD.Ts > 0 ? LD.Ts : 0.5;
    var tanC = (gs - gw) / gs * tp, dep = this.dep, cnt = 0, vol = 0;
    for (var c = 0; c < N; c++) {
      if (!dep[c]) continue;
      var Tq = z[c] - zfix[c]; if (Tq < Ts) continue;        // すべりの厚さ Ts より薄い所（ダムの縁）は対象にしない
      if (W[c] < 0.98 * n * Tq || h[c] > 0.1) continue;   // 地表まで飽和、湛水の中ではない
      var i = c % nx, j = (c - i) / nx, best = -1, bs = 0;
      var nb = [i < nx - 1 ? c + 1 : -1, i > 0 ? c - 1 : -1, j < ny - 1 ? c + nx : -1, j > 0 ? c - nx : -1];
      for (var q = 0; q < 4; q++) { var o = nb[q]; if (o < 0 || this.kind[o] === 2 || (this.basin && !this.basin[o])) continue; var s = (z[c] - z[o]) / dx; if (s > bs) { bs = s; best = o; } }
      if (best < 0) continue;
      var b = Math.atan(bs), Fs = coh / (gs * Ts * Math.sin(b) * Math.cos(b)) + (gs - gw) / gs * tp / bs;
      if (Fs >= 1) continue;
      var dl = Math.min((z[c] - z[best] - tanC * dx) / 2, Tq);
      if (!(dl > 1e-4)) continue;
      var got = this.moveSed(c, best, dl);
      if (!(got > 0)) continue;
      this.dzs[c] -= got; this.dzs[best] += got;
      if (this.kind[best] !== 0) { this.kind[best] = 0; this.zfix[best] = z[best] - got; this.fl[best] = 1; }
      this.por[best] = 1; dep[best] = 1;
      var wm = Math.min(W[c], n * got); W[c] -= wm;
      var cap = n * Math.max(0, z[best] - zfix[best]), room = Math.max(0, cap - W[best]), ww = Math.min(wm, room); W[best] += ww; h[best] += wm - ww;
      cnt++; vol += got * dx * dx;
    }
    if (cnt) { this.progCount = (this.progCount || 0) + cnt; this.progVol = (this.progVol || 0) + vol; }
  };

  // すべり破壊（指定）：湛水位が指定の高さに達したら、ダムのうち、下流の法尻から角度 αs で上流へ上がる面より上の土を、
  // 法尻から下流へ堆積勾配 θdep の斜面になるように瞬時に動かす（Awal ら（2008）、高山ら（2020）の簡易モデルの 2 次元版）
  Sim.prototype.ldamSlip = function () {
    var LD = this.P.ldam || {}; if (LD.fail !== 'slip' || !this.ldamList.length) return;
    var nx = this.nx, dx = this.dx, z = this.z, h = this.h, zfix = this.zfix, a1 = dx * dx, K = this.K;
    for (var m = 0; m < this.ldamList.length; m++) {
      var D = this.ldamList[m]; if (D.slipped) continue;
      // 湛水位：ダムの上流側（s < 0）で、ダムのセルに接する濡れたセルのうち、河床がいちばん低いセルの水位
      var lev = LD.slipLevel > 0 ? D.zc + LD.slipLevel * D.H : D.crest - 0.05 * D.H;
      var wl = -1e30, zlo = 1e30, sTo = -1e30, cToe = -1, k, c, isD = this._isD || (this._isD = new Uint8Array(this.N));
      for (k = 0; k < D.cells.length; k++) {
        c = D.cells[k]; var px = (c % nx + 0.5) * dx - D.x0, py = (Math.floor(c / nx) + 0.5) * dx - D.y0, s = px * D.ux + py * D.uy;
        if (z[c] - zfix[c] > 0.05) { isD[c] = 1; if (s > sTo) { sTo = s; cToe = c; } }
      }
      for (k = 0; k < D.cells.length; k++) {
        c = D.cells[k]; if (!isD[c]) continue;
        var ci2 = c % nx, cj2 = Math.floor(c / nx), nb = [ci2 > 0 ? c - 1 : -1, ci2 < nx - 1 ? c + 1 : -1, cj2 > 0 ? c - nx : -1, cj2 < this.ny - 1 ? c + nx : -1];
        for (var q2 = 0; q2 < 4; q2++) {
          var o = nb[q2]; if (o < 0 || isD[o] || h[o] <= 0.05) continue;
          var so = ((o % nx + 0.5) * dx - D.x0) * D.ux + ((Math.floor(o / nx) + 0.5) * dx - D.y0) * D.uy;
          if (so < 0 && z[o] < zlo) { zlo = z[o]; wl = z[o] + h[o]; }   // いちばん河床の低い所の水位（湛水の水面）。斜面を流れ下る水は除く
        }
      }
      for (k = 0; k < D.cells.length; k++) isD[D.cells[k]] = 0;
      if (wl < lev || cToe < 0) continue;
      D.slipped = true;
      // すべり面：下流の法尻（元の地盤）から、指定の角度で上流へ上がる平面。角度 0（既定）は法尻と上流側の天端の端を結ぶ面
      var ztoe = zfix[cToe], sUp = -(D.L && D.L.Wc > 0 ? D.L.Wc : 0) / 2;
      var al = LD.slipDeg > 0 ? Math.tan(LD.slipDeg * Math.PI / 180) : Math.max(0.05, (D.crest - ztoe) / Math.max(dx, sTo - sUp));
      var tanD = Math.tan((LD.thDep > 0 ? LD.thDep : 10) * Math.PI / 180);
      var sed = new Float64Array(K), Vb = 0, Wm = 0, n = this.P.porosity;
      for (k = 0; k < D.cells.length; k++) {
        c = D.cells[k]; var px3 = (c % nx + 0.5) * dx - D.x0, py3 = (Math.floor(c / nx) + 0.5) * dx - D.y0, s3 = px3 * D.ux + py3 * D.uy;
        var zp = ztoe + (sTo - s3) * al, base = Math.max(zfix[c], zp), ex = z[c] - base;
        if (ex <= 1e-3) continue;
        var got = this.takeSed(c, ex, sed); z[c] -= got; this.dzs[c] -= got; Vb += got * a1;
        if (this.W) { var cap = n * Math.max(0, z[c] - zfix[c]); if (this.W[c] > cap) { Wm += (this.W[c] - cap) * a1; this.W[c] = cap; } }
      }
      if (Vb <= 0) continue;
      var ss = 0; for (k = 0; k < K; k++) ss += sed[k];
      var fr = new Float64Array(K); for (k = 0; k < K; k++) fr[k] = ss > 0 ? sed[k] / ss : this.f0[k];
      this.depositRamp(D, sTo, Vb, tanD, Wm, fr);
      this.slideLog.push({ t: this.t, id: 'ldam' + m, mode: 'slip', vol: Vb, water: 0 });
    }
  };
  // 法尻（s = s0）から下流へ、勾配 tanT で下がる斜面の形に体積 V（かさ）を置く
  Sim.prototype.depositRamp = function (D, s0, V, tanT, Wv, frc) {
    var nx = this.nx, ny = this.ny, dx = this.dx, z = this.z, kind = this.kind, a1 = dx * dx, self = this;
    var r0 = Math.ceil((D.H * 12 + D.lat) / dx), ci = D.c0 % nx, cj = Math.floor(D.c0 / nx);
    function vol(Z, apply) {
      var s = 0, cells = apply ? [] : null;
      for (var j = Math.max(0, cj - r0); j <= Math.min(ny - 1, cj + r0); j++) for (var i = Math.max(0, ci - r0); i <= Math.min(nx - 1, ci + r0); i++) {
        var c = j * nx + i; if (kind[c] === 2 || (self.basin && !self.basin[c])) continue;
        var px = (i + 0.5) * dx - D.x0, py = (j + 0.5) * dx - D.y0, sv = px * D.ux + py * D.uy, w = -px * D.uy + py * D.ux;
        if (sv < s0 - dx || Math.abs(w) > D.lat) continue;
        var zt = Z - tanT * Math.max(0, sv - s0);
        if (zt > z[c]) { s += (zt - z[c]) * a1; if (apply) cells.push([c, zt - z[c]]); }
      }
      return apply ? cells : s;
    }
    var lo = D.zc - 50, hi = D.crest + 50, it;
    for (it = 0; it < 60 && vol(hi) < V; it++) hi += 50;
    for (it = 0; it < 60; it++) { var mid = 0.5 * (lo + hi); if (vol(mid) < V) lo = mid; else hi = mid; }
    var cells = vol(hi, true), tot = 0, q;
    for (q = 0; q < cells.length; q++) tot += cells[q][1];
    this.depositCells(cells, tot > 0 ? V / (tot * a1) : 0, Wv, frc);
  };

  // ---- 伏流（河床への浸透・地下水の流れ・湧出）--------------------------------------------
  // 透水係数 k：指定がなければ交換層の粒度から Kozeny–Carman 式 k = (g/ν)·n³/(1−n)²·d_h²/180（d_h は調和平均径）、上限 kMax。
  // 浸透：水に浸かったセルで、堆積層の空隙が残っていれば単位動水勾配で k の速さ。
  // 地下水：堆積層の底（侵食の下限）の上の飽和厚 b = W/n、水頭 zfix + b の差で Darcy 流（Dupuit 近似）。
  // 湧出：W が堆積層の空隙量 n(z − zfix) を超えた分は地表へ戻す。下流端からは基準勾配で流出させる。
  Sim.prototype.gwStep = function (dt) {
    var P = this.P, GW = P.gw || {}, N = this.N, nx = this.nx, ny = this.ny, dx = this.dx, K = this.K, n = P.porosity;
    var h = this.h, z = this.z, zfix = this.zfix, kind = this.kind, W = this.W, kg = this.kgw, A = this.A, d = this.d, qx = this.qx, qy = this.qy;
    var Cst = 1 - n, dbr = P.transport === 'debris', Cc = this.Cc, kMax = GW.kMax > 0 ? GW.kMax : 0.05, HD = 5e-4, c, k, i, j;   // 浸透で残す水深は流れの乾湿の境（1e-3）より小さくする（境目に止まると丸め誤差で乾湿が入れ替わる）
    var kc0 = G / 1e-6 * n * n * n / ((1 - n) * (1 - n)) / 180, bmax = 0, kmx = 0, por = this.por, DMk = P.dem ? P.dem.ks : 0, LDk = P.ldam ? P.ldam.k : 0;
    for (c = 0; c < N; c++) {
      if (!por[c]) { kg[c] = 0; continue; }
      var kk = kind[c] === 3 ? (DMk > 0 ? DMk : 1e-5) : (this.dep[c] && LDk > 0 ? LDk : (GW.k > 0 ? GW.k : 0));
      if (!kk) { var ck = c * K, sA = 0, sh = 0; for (k = 0; k < K; k++) { sA += A[ck + k]; sh += A[ck + k] / d[k]; } var dh = sh > 0 ? sA / sh : this.dm0; kk = kc0 * dh * dh; }
      kg[c] = kk > kMax ? kMax : kk;
      var Wm = n * Math.max(0, z[c] - zfix[c]), hh = h[c];
      if (W[c] > Wm) { var ex = W[c] - Wm; W[c] = Wm; h[c] = hh + ex; hh = h[c]; }
      if (hh > HD && W[c] < Wm) {
        var av = hh - HD;
        if (dbr) { var S = 0, ck2 = c * K; for (k = 0; k < K; k++) S += Cc[ck2 + k]; av = Math.min(av, hh - S / (0.98 * Cst)); }
        var f = kg[c] * dt; if (f > Wm - W[c]) f = Wm - W[c]; if (f > av) f = av;
        if (f > 0) { var r = (hh - f) / hh; h[c] = hh - f; qx[c] *= r; qy[c] *= r; W[c] += f; }
      }
      var b = W[c] / n; if (b > bmax) bmax = b; if (kg[c] > kmx) kmx = kg[c];
    }
    // 地下水の流れ（陽解法。安定条件で小分けにする）
    var dtg = kmx * bmax > 0 ? 0.2 * n * dx * dx / (kmx * bmax) : dt, ns = Math.min(50, Math.max(1, Math.ceil(dt / dtg))), dts = dt / ns;
    var Hg = this._gwH || (this._gwH = new Float64Array(N)), dW = this._gwD || (this._gwD = new Float64Array(N)), Ib = this.Ix[nx - 1];
    for (var it = 0; it < ns; it++) {
      for (c = 0; c < N; c++) { Hg[c] = zfix[c] + W[c] / n; dW[c] = 0; }
      for (j = 0; j < ny; j++) for (i = 0; i < nx; i++) {
        c = j * nx + i; if (!por[c]) continue;
        for (var a = 0; a < 2; a++) {
          var o = a === 0 ? (i < nx - 1 ? c + 1 : -1) : (j < ny - 1 ? c + nx : -1);
          if (o < 0 || !por[o]) continue;
          var dH = Hg[c] - Hg[o]; if (dH === 0) continue;
          var up = dH > 0 ? c : o, kf = kg[c] > 0 && kg[o] > 0 ? 2 * kg[c] * kg[o] / (kg[c] + kg[o]) : 0;
          var vol = kf * (W[up] / n) * dH / dx * dts / dx;            // 移る水の量（層厚） [m]
          var cap = 0.25 * W[up]; if (vol > cap) vol = cap; else if (vol < -cap) vol = -cap;
          dW[c] -= vol; dW[o] += vol;
        }
        if (i === nx - 1) { var vo = kg[c] * (W[c] / n) * Ib * dts / dx; if (vo > 0.25 * W[c]) vo = 0.25 * W[c]; dW[c] -= vo; this.gwOut += vo * dx * dx; }
      }
      for (c = 0; c < N; c++) W[c] += dW[c];
    }
    for (c = 0; c < N; c++) {     // 湧出
      if (!por[c]) continue;
      var Wm2 = n * Math.max(0, z[c] - zfix[c]);
      if (W[c] > Wm2) { h[c] += W[c] - Wm2; W[c] = Wm2; }
    }
  };
  Sim.prototype.groundWater = function () { var s = 0, W = this.W; for (var c = 0; c < this.N; c++) s += W[c]; return s * this.dx * this.dx; };

  // ---- 透過型堰堤：開口を通る土砂（粒径別）と閉塞 ---------------------------------------
  // 開口のセルの面を通る流砂量に、通過率 f_k を掛ける。f_k = 0（閉塞後、または d_k ≥ l/α の礫）、それ以外は開口率 φs
  Sim.prototype.openFilter = function (SFx, SFy) {
    var D = this.dams, nx = this.nx, ny = this.ny, K = this.K, d = this.d, OD = this.P.openDam || {}, al = OD.alpha > 0 ? OD.alpha : 1;
    var stX = this._ofX || (this._ofX = new Int32Array(ny * (nx + 1))), stY = this._ofY || (this._ofY = new Int32Array((ny + 1) * nx));
    var fK = this._ofK || (this._ofK = new Float64Array(K)), tag = this._ofT = (this._ofT | 0) + 1, q, k, m, c, i, j, f;
    for (q = 0; q < D.length; q++) {
      var op = D[q].op; if (!op) continue;
      for (k = 0; k < K; k++) fK[k] = op.blocked ? 0 : (d[k] * al < op.l ? op.phiS : 0);
      for (m = 0; m < op.cells.length; m++) {
        c = op.cells[m]; i = c % nx; j = (c - i) / nx;
        f = j * (nx + 1) + i;
        if (stX[f] !== tag) { stX[f] = tag; for (k = 0; k < K; k++) SFx[f * K + k] *= fK[k]; }
        if (stX[f + 1] !== tag) { stX[f + 1] = tag; for (k = 0; k < K; k++) SFx[(f + 1) * K + k] *= fK[k]; }
        if (stY[c] !== tag) { stY[c] = tag; for (k = 0; k < K; k++) SFy[c * K + k] *= fK[k]; }
        if (stY[c + nx] !== tag) { stY[c + nx] = tag; for (k = 0; k < K; k++) SFy[(c + nx) * K + k] *= fK[k]; }
      }
    }
  };
  // 閉塞の判定：開口のすぐ上流の流れの土砂濃度 C（流砂量 ÷ 流量）が Cblk 以上で、流れている土砂の d95 に対し l ≤ β d95 なら閉塞。
  // 閉塞後は不透過型と同じに働く（開口の損失係数は φ = phiBlock、土砂は通さない）
  Sim.prototype.checkBlock = function (useCc) {
    var D = this.dams, nx = this.nx, ny = this.ny, K = this.K, h = this.h, kind = this.kind, openC = this.openC, FxM = this.FxM, GyM = this.GyM;
    if (!this.dord) { this.dord = []; for (var k0 = 0; k0 < this.K; k0++) this.dord.push(k0); var dd0 = this.d; this.dord.sort(function (a, b) { return dd0[a] - dd0[b]; }); }
    var OD = this.P.openDam || {}, Cb = OD.Cblk > 0 ? OD.Cblk : 0.2, ord = this.dord, d = this.d, dx = this.dx;
    var tk = this._cbK || (this._cbK = new Float64Array(K)), q, k, m;
    for (q = 0; q < D.length; q++) {
      var Dq = D[q], op = Dq.op; if (!op || op.blocked) continue;
      var sS = 0, sQ = 0; for (k = 0; k < K; k++) tk[k] = 0;
      for (m = 0; m < op.cells.length; m++) {
        var c0 = op.cells[m], i0 = c0 % nx, j0 = (c0 - i0) / nx, nbs = [i0 > 0 ? c0 - 1 : -1, i0 < nx - 1 ? c0 + 1 : -1, j0 > 0 ? c0 - nx : -1, j0 < ny - 1 ? c0 + nx : -1];
        for (var a = 0; a < 4; a++) {
          var nb = nbs[a]; if (nb < 0 || openC[nb] || kind[nb] !== 0 || h[nb] < 0.02) continue;
          var ib = nb % nx, jb = (nb - ib) / nx;
          if (((ib + 0.5) * dx - Dq.cx) * Dq.tx + ((jb + 0.5) * dx - Dq.cy) * Dq.ty >= 0) continue;   // 上流側だけ
          var fw = jb * (nx + 1) + ib, uu = 0.5 * (FxM[fw] + FxM[fw + 1]) / h[nb], vv = 0.5 * (GyM[nb] + GyM[nb + nx]) / h[nb], V = Math.sqrt(uu * uu + vv * vv);
          sQ += h[nb] * V;
          for (k = 0; k < K; k++) {
            var qk = useCc ? this.Cc[nb * K + k] * V : Math.sqrt(this.qbx[nb * K + k] * this.qbx[nb * K + k] + this.qby[nb * K + k] * this.qby[nb * K + k]);
            tk[k] += qk; sS += qk;
          }
        }
      }
      if (!(sQ > 1e-9) || sS / sQ < Cb) continue;
      var cum = 0, d95 = d[ord[K - 1]];
      for (k = 0; k < K; k++) { cum += tk[ord[k]]; if (cum >= 0.95 * sS) { d95 = d[ord[k]]; break; } }
      if (op.l <= op.beta * d95) { op.blocked = true; op.tBlock = this.t; op.Cblk = sS / sQ; op.d95 = d95; }
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
    if (this.hasOpen) { this.checkBlock(true); this.openFilter(SFx, SFy); }
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
    this.tribFeed(dts, false);
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
    gx = (a >= 0 && b >= 0) ? (z[b] - z[a]) / (2 * dx) : a >= 0 ? (z[c] - z[a]) / dx : b >= 0 ? (z[b] - z[c]) / dx : -this.Ix[i];
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
