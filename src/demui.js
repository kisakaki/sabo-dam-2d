/* DEM と流域の画面（地図・流域の選択・供給点・雨・斜面崩壊・天然ダム）
 * ui.js の readForm から SaboDEMUI.params(p) を呼んで、計算の条件 p を作る */
(function () {
  'use strict';
  var SD = window.SaboDEM, SH = window.SaboHydro;
  var $ = function (id) { return document.getElementById(id); };
  var ST = { D: null, A: null, mask: null, outlet: -1, bid: -1, heads: [], src: [], slides: [], ldams: [], dams: [], matz: [], hist: [], img: null, view: null, rainFrames: null, lat: null, dom: null };

  function fmt(v, n) { return isFinite(v) ? (+v).toFixed(n) : ''; }
  function info(msg, err) { var e = $('dem_info'); e.textContent = msg; e.style.color = err ? 'var(--err,#b3261e)' : ''; }

  // ------------------------------------------------------------ 読み込み
  function inflateBrowser(u8) {
    var ds = new DecompressionStream('deflate');
    var w = ds.writable.getWriter(); w.write(u8); w.close();
    return new Response(ds.readable).arrayBuffer().then(function (b) { return new Uint8Array(b); });
  }
  // ZIP の展開（基盤地図情報のダウンロードは ZIP。ZIP の中の ZIP も開く）
  function inflateRaw(u8) {
    var ds = new DecompressionStream('deflate-raw');
    var w = ds.writable.getWriter(); w.write(u8); w.close();
    return new Response(ds.readable).arrayBuffer().then(function (b) { return new Uint8Array(b); });
  }
  function unzip(buf, depth) {
    var u8 = new Uint8Array(buf), dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength), e = -1;
    for (var i = u8.length - 22; i >= Math.max(0, u8.length - 65557); i--) if (dv.getUint32(i, true) === 0x06054b50) { e = i; break; }
    if (e < 0) return Promise.reject(new Error('ZIP として読めません'));
    var n = dv.getUint16(e + 10, true), p = dv.getUint32(e + 16, true), jobs = [];
    for (var k = 0; k < n; k++) {
      if (dv.getUint32(p, true) !== 0x02014b50) break;
      var meth = dv.getUint16(p + 10, true), csz = dv.getUint32(p + 20, true), nl = dv.getUint16(p + 28, true), xl = dv.getUint16(p + 30, true), cl = dv.getUint16(p + 32, true), lo = dv.getUint32(p + 42, true);
      var name = new TextDecoder().decode(u8.subarray(p + 46, p + 46 + nl));
      p += 46 + nl + xl + cl;
      if (/\/$/.test(name)) continue;
      var d0 = lo + 30 + dv.getUint16(lo + 26, true) + dv.getUint16(lo + 28, true), raw = u8.subarray(d0, d0 + csz);
      (function (name, meth, raw) {
        var pr = meth === 0 ? Promise.resolve(raw.slice()) : meth === 8 ? inflateRaw(raw) : Promise.reject(new Error(name + '：未対応の圧縮方式 ' + meth));
        jobs.push(pr.then(function (b) {
          if (/\.zip$/i.test(name) && (depth || 0) < 3) return unzip(b.buffer, (depth || 0) + 1);
          return [{ name: name.split('/').pop(), data: b }];
        }));
      })(name, meth, raw);
    }
    return Promise.all(jobs).then(function (a) { return [].concat.apply([], a); });
  }
  function loadFiles(files) {
    if (!files || !files.length) return;
    info('読み込み中…');
    var arr = Array.prototype.slice.call(files);
    // ZIP は展開して中身のファイルに置き換える
    Promise.all(arr.map(function (f) {
      if (/\.zip$/i.test(f.name)) return f.arrayBuffer().then(function (b) { return unzip(b, 0); });
      return f.arrayBuffer().then(function (b) { return [{ name: f.name, data: new Uint8Array(b) }]; });
    })).then(function (a) {
      var list = [].concat.apply([], a), dec = new TextDecoder();
      var by = function (re) { return list.filter(function (f) { return re.test(f.name.toLowerCase()); }); };
      var xml = by(/\.xml$/), tif = by(/\.tiff?$/), asc = by(/\.(asc|txt|xyz|csv|dat)$/);
      var label = (arr.length === 1 ? arr[0].name : arr[0].name + ' ほか ' + (arr.length - 1) + ' 個');
      // 基盤地図情報の XML を優先。標高の XML でないもの（メタデータなど）は除く
      xml = xml.filter(function (f) { return /ElevationModel|DEM|tupleList/.test(dec.decode(f.data.subarray(0, 4000))) || /dem/i.test(f.name); });
      if (xml.length) { label += '（XML ' + xml.length + ' 個）'; return { D: SD.parseGSIXML(xml.map(function (f) { return dec.decode(f.data); })), label: label }; }
      if (tif.length) return SD.parseTIFF(tif[0].data.buffer.slice(tif[0].data.byteOffset, tif[0].data.byteOffset + tif[0].data.byteLength), inflateBrowser).then(function (D) { return { D: D, label: label }; });
      if (asc.length) { var t = dec.decode(asc[0].data); return { D: /^\s*ncols/i.test(t) ? SD.parseASC(t) : SD.parseXYZ(t), label: label }; }
      throw new Error('DEM のファイルが見つかりません（読めるのは 基盤地図情報の XML・GeoTIFF・ESRI ASCII・XYZ と、それらの ZIP）。中身：' + list.slice(0, 5).map(function (f) { return f.name; }).join('、'));
    }).then(function (r) {
      if (!r) return;
      setDEM(r.D, r.label);
    }).catch(function (e) { info('読み込めませんでした：' + (e && e.message || e), true); console.error(e); });
  }
  function setDEM(D, label) {
    var N = D.nx * D.ny, k = 1;
    while (N / (k * k) > 9e6) k++;
    if (k > 1) { D = SD.coarsen(D, k); label += '、大きいので ' + k + '×' + k + ' セルを平均して粗くしました'; }
    if (D.skipped) label += '、読めないタイル ' + D.skipped + ' 個を飛ばしました';
    var nh = SD.fillHoles(D); if (nh) label += '、内側の欠測 ' + nh.toLocaleString() + ' セルを周りの値で埋めました';
    var nv = 0; for (var c = 0; c < D.nx * D.ny; c++) if (isFinite(D.z[c])) nv++;
    if (nv < 100) { info('標高の値がほとんどありません（' + nv + ' セル）。ファイルを確かめてください。', true); return; }
    ST.D0 = D; ST.label = label; ST.mask = null; ST.bid = -1; ST.outlet = -1; ST.src = []; ST.slides = []; ST.ldams = []; ST.dams = []; ST.matz = []; ST.hist = []; ST.draft = []; ST.mode = 'basin'; applySmooth();
  }
  // 平滑化（元の DEM は ST.D0 に残す）。谷筋・流域を取り直す。出口を選んでいたら、その近くの谷筋に選び直す
  function applySmooth() {
    if (!ST.D0) return;
    var sg = Math.max(0, parseFloat($('dem_smooth').value) || 0), t0 = performance.now();
    ST.D = sg > 0 ? SD.smooth(ST.D0, sg) : ST.D0; ST.smoothMs = Math.round(performance.now() - t0);
    var keepOut = ST.mode === 'outlet' ? ST.outlet : -1;
    analyze();
    if (keepOut >= 0) { var s2 = nearStream(keepOut, 3); if (s2 >= 0) selectOutlet(s2); draw(); }
  }
  function analyze() {
    var D = ST.D; if (!D) return;
    var t0 = performance.now();
    ST.A = SD.analyze(D, { thrArea: Math.max(0.0005, parseFloat($('dem_thr').value) || 0.02) * 1e6 });
    var A = ST.A, big = A.basins.slice().sort(function (a, b) { return b.area - a.area; });
    ST.order = big;
    if (ST.bid >= 0 && ST.bid < A.basins.length && ST.mask && ST.mode !== 'outlet') selectBasin(ST.bid, true);
    else selectBasin(big.length ? big[0].id : -1, true);
    var zmin = Infinity, zmax = -Infinity; for (var c = 0; c < D.nx * D.ny; c++) { var v = D.z[c]; if (isFinite(v)) { if (v < zmin) zmin = v; if (v > zmax) zmax = v; } }
    info(ST.label + '：' + D.nx + ' × ' + D.ny + '（格子 ' + fmt(D.cs, 2) + ' m' + (D.smooth ? '、平滑化 ' + fmt(D.smooth, 0) + ' m' : '') + '、標高 ' + fmt(zmin, 0) + '〜' + fmt(zmax, 0) + ' m）。谷頭 ' + A.heads.length + ' 点、出口ごとの流域 ' + A.basins.length + '（' + Math.round(performance.now() - t0) + ' ms）');
    ST.img = null; ST.det = null; ST.zst = null; ST.cam = null; draw();
  }
  function selectBasin(id, quiet) {
    var A = ST.A; if (!A || id < 0) { ST.mask = null; return; }
    var N = ST.D.nx * ST.D.ny, m = new Uint8Array(N);
    for (var c = 0; c < N; c++) if (A.bid[c] === id) m[c] = 1;
    ST.mask = m; ST.outlet = A.basins[id].outlet; ST.bid = id; ST.mode = 'basin';
    afterBasin(quiet);
  }
  function selectOutlet(c) {
    ST.mask = SD.upstream(ST.D, ST.A, c); ST.outlet = c; ST.bid = ST.A.bid[c]; ST.mode = 'outlet';
    afterBasin(false);
  }
  function afterBasin(quiet) {
    var m = ST.mask, A = ST.A;
    ST.heads = A.heads.filter(function (c) { return m[c]; });
    ST.src = ST.src.filter(function (r) { return m[r.c]; });
    ST.slides = ST.slides.filter(function (r) { return m[r.c]; });
    ST.ldams = ST.ldams.filter(function (r) { return m[r.c]; });
    ST.dams = ST.dams.filter(function (r) { return m[r.c]; });
    if (!ST.src.length && ST.heads.length) topHeads(3);
    ST.lat = null; ST.mcv = null;
    if ($('dm_fit').checked && ST.D) fitBasin();
    renderTables(); draw();
    if (!quiet) selInfo();
  }
  function selInfo() {
    if (!ST.mask) { $('dem_sel').textContent = ''; return; }
    var n = 0; for (var c = 0; c < ST.mask.length; c++) n += ST.mask[c];
    var a = n * ST.D.cs * ST.D.cs / 1e6;
    $('dem_sel').textContent = '選んだ流域：面積 ' + fmt(a, a < 1 ? 3 : 2) + ' km²、谷頭 ' + ST.heads.length + ' 点、供給点 ' + ST.src.length + ' 点、斜面崩壊 ' + ST.slides.length + ' か所、天然ダム ' + ST.ldams.length + ' か所。「▶ 計算開始」で計算領域に切り出します。';
  }
  function topHeads(n) {
    var A = ST.A, hs = ST.heads.slice().sort(function (a, b) { return A.acc[b] - A.acc[a]; }).slice(0, n);
    ST.src = hs.map(function (c, i) { return srcDefault(c, i); });
  }
  function srcDefault(c) {
    var db = $('gen_type').value !== 'bed' ? false : ($('s_tr').value === 'debris');
    return { c: c, mode: db ? 'debris' : 'tri', Qb: 0.2, Qp: 5, t0: 30, tr: 60, tf: 300, C: db ? 0.3 : 0.01 };
  }

  // ------------------------------------------------------------ 表
  function opt(v, list) { return list.map(function (o) { return '<option value="' + o[0] + '"' + (o[0] === v ? ' selected' : '') + '>' + o[1] + '</option>'; }).join(''); }
  function renderTables() {
    var tb = $('srctab').tBodies[0]; tb.innerHTML = '';
    ST.src.forEach(function (r, i) {
      var tr = document.createElement('tr');
      tr.innerHTML = '<td>' + (i + 1) + '</td><td><select>' + opt(r.mode, [['tri', '三角形'], ['const', '一定'], ['debris', '土石流の一波']]) + '</select></td>' +
        ['Qb', 'Qp', 't0', 'tr', 'tf', 'C'].map(function (k) { return '<td><input type="number" step="any" value="' + r[k] + '" data-k="' + k + '"></td>'; }).join('') +
        '<td><button class="x" aria-label="外す">×</button></td>';
      tr.querySelector('select').onchange = function () { r.mode = this.value; };
      Array.prototype.forEach.call(tr.querySelectorAll('input'), function (inp) { inp.onchange = function () { r[inp.dataset.k] = parseFloat(inp.value) || 0; }; });
      tr.querySelector('.x').onclick = function () { snap(); ST.src.splice(i, 1); renderTables(); draw(); selInfo(); };
      tb.appendChild(tr);
    });
    tb = $('sltab').tBodies[0]; tb.innerHTML = '';
    ST.slides.forEach(function (r, i) {
      var tr = document.createElement('tr');
      tr.innerHTML = '<td>' + (i + 1) + '</td><td><select data-k="type">' + opt(r.type, [['shallow', '表層'], ['deep', '深層']]) + '</select></td>' +
        (r.poly ? '<td title="多角形（頂点 ' + r.poly.length + '）">' + fmt(polyArea(r.poly), 0) + ' m²</td>' : '<td><input type="number" step="any" value="' + r.R + '" data-k="R"></td>') + '<td><input type="number" step="any" value="' + r.D + '" data-k="D"></td>' +
        '<td><select data-k="trigger">' + opt(r.trigger, [['fs', '安全率'], ['time', '時刻']]) + '</select></td><td><input type="number" step="any" value="' + r.t + '" data-k="t"></td>' +
        '<td><select data-k="move">' + opt(r.move, [['auto', '自動'], ['fluid', '土石流化'], ['block', '土塊']]) + '</select></td>' +
        '<td><input type="number" step="any" value="' + r.phi + '" data-k="phi"></td><td><input type="number" step="any" value="' + r.coh + '" data-k="coh"></td><td><input type="number" step="any" value="' + r.m0 + '" data-k="m0"></td>' +
        '<td><button class="x" aria-label="消す">×</button></td>';
      Array.prototype.forEach.call(tr.querySelectorAll('input,select'), function (inp) {
        inp.onchange = function () {
          snap(); var k = inp.dataset.k; r[k] = inp.tagName === 'SELECT' ? inp.value : (parseFloat(inp.value) || 0);
          if (k === 'type') { var d = slideDefaults(r.type); r.R = d.R; r.D = d.D; r.phi = d.phi; r.coh = d.coh; renderTables(); }
          draw();
        };
      });
      tr.querySelector('.x').onclick = function () { snap(); ST.slides.splice(i, 1); renderTables(); draw(); selInfo(); };
      tb.appendChild(tr);
    });
    tb = $('ldtab').tBodies[0]; tb.innerHTML = '';
    ST.ldams.forEach(function (r, i) {
      var sh = ldShape(r), tr = document.createElement('tr');
      tr.innerHTML = '<td>' + (i + 1) + '</td>' + ['H', 'Wc', 'up', 'dn', 'W'].map(function (k) { return '<td><input type="number" step="any" value="' + (r[k] || 0) + '" data-k="' + k + '"></td>'; }).join('') +
        '<td><input type="number" step="any" value="' + Math.round(sh.V) + '" data-k="V" title="体積を入れると、それに合う高さにします"></td><td>' + fmt(sh.L, 0) + ' × ' + fmt(sh.B, 0) + '</td><td><button class="x" aria-label="消す">×</button></td>';
      Array.prototype.forEach.call(tr.querySelectorAll('input'), function (inp) {
        inp.onchange = function () {
          snap(); var k = inp.dataset.k, v = parseFloat(inp.value) || 0;
          if (k === 'V') fitLdamV(r, v); else r[k] = v;
          r.sh = null; renderTables(); draw();
        };
      });
      tr.querySelector('.x').onclick = function () { snap(); ST.ldams.splice(i, 1); renderTables(); draw(); selInfo(); };
      tb.appendChild(tr);
    });
    tb = $('dmtab').tBodies[0]; tb.innerHTML = '';
    ST.dams.forEach(function (r, i) {
      var tr = document.createElement('tr'), op = r.type !== 'closed';
      tr.innerHTML = '<td>' + (i + 1) + '</td><td><select data-k="type">' + opt(r.type, [['closed', '不透過'], ['slit', 'スリット'], ['mslit', '複数スリット'], ['grid', '鋼製格子型'], ['sslit', '鋼製スリット型']]) + '</select></td>' +
        ['H', 'b', 'len'].map(function (k) { return '<td><input type="number" step="any" value="' + (r[k] || 0) + '" data-k="' + k + '"></td>'; }).join('') +
        ['l', 'n', 'hb'].map(function (k) { return '<td><input type="number" step="any" value="' + (r[k] || 0) + '" data-k="' + k + '"' + (op ? '' : ' disabled') + '></td>'; }).join('') +
        '<td>' + fmt(damLen(r), 0) + '</td><td><button class="x" aria-label="消す">×</button></td>';
      Array.prototype.forEach.call(tr.querySelectorAll('input,select'), function (inp) {
        inp.onchange = function () {
          snap(); var k = inp.dataset.k; r[k] = inp.tagName === 'SELECT' ? inp.value : (parseFloat(inp.value) || 0);
          if (k === 'type') { r.l = r.type === 'grid' || r.type === 'sslit' ? 1 : 2; r.n = r.type === 'mslit' ? 3 : 1; }
          renderTables(); draw();
        };
      });
      tr.querySelector('.x').onclick = function () { snap(); ST.dams.splice(i, 1); renderTables(); draw(); selInfo(); };
      tb.appendChild(tr);
    });
    $('dm_back').disabled = !ST.hist.length;
  }
  // ------------------------------------------------------------ 取り消し（元に戻す）
  function snap() {
    ST.hist.push(JSON.stringify({ src: ST.src, slides: ST.slides, ldams: ST.ldams.map(function (r) { var o = Object.assign({}, r); delete o.sh; return o; }), dams: ST.dams, matz: ST.matz }));
    if (ST.hist.length > 50) ST.hist.shift();
  }
  function undo() {
    if (!ST.hist.length) return;
    var o = JSON.parse(ST.hist.pop()); ST.src = o.src; ST.slides = o.slides; ST.ldams = o.ldams; ST.dams = o.dams; ST.matz = o.matz || [];
    renderTables(); draw(); selInfo(); $('dem_sel').textContent = '1 つ前に戻しました。' + (ST.hist.length ? '' : '（これ以上は戻せません）');
  }
  // クリックした点の近く（画面で 14 px 以内）にある物を外す。返り値：外した物の名前
  function removeAt(c) {
    var D = ST.D, x = c % D.nx, y = (c / D.nx) | 0, sc = ST.view.sc, best = null, bd = 14 / sc;
    function cand(list, i, cc, name) { var d = Math.hypot(cc % D.nx - x, ((cc / D.nx) | 0) - y); if (d < bd) { bd = d; best = { list: list, i: i, name: name }; } }
    ST.src.forEach(function (r, i) { cand(ST.src, i, r.c, '供給点 ' + (i + 1)); });
    ST.slides.forEach(function (r, i) {
      if (r.poly && inPolyCells(r.poly, x, y)) { bd = 0; best = { list: ST.slides, i: i, name: '崩壊 ' + (i + 1) }; return; }
      var dd = Math.hypot(r.c % D.nx - x, ((r.c / D.nx) | 0) - y); if (!r.poly && dd <= Math.max(r.R / D.cs, 14 / sc) && dd < bd + r.R / D.cs) { bd = Math.min(bd, dd); best = { list: ST.slides, i: i, name: '崩壊 ' + (i + 1) }; }
    });
    ST.ldams.forEach(function (r, i) { cand(ST.ldams, i, r.c, '天然ダム ' + (i + 1)); });
    ST.dams.forEach(function (r, i) { cand(ST.dams, i, r.c, '堰堤 ' + (i + 1)); });
    if (!best) ST.matz.forEach(function (r, i) { if (inPolyCells(r.poly, x, y)) best = { list: ST.matz, i: i, name: '河床材料 ' + ('ABCDEFGHIJKL'[r.mat] || '?') + ' の範囲' }; });
    if (!best) return null;
    snap(); best.list.splice(best.i, 1); return best.name;
  }
  function inPolyCells(poly, x, y) {
    var D = ST.D, a = false;
    for (var i = 0, j = poly.length - 1; i < poly.length; j = i++) {
      var xi = poly[i] % D.nx, yi = (poly[i] / D.nx) | 0, xj = poly[j] % D.nx, yj = (poly[j] / D.nx) | 0;
      if ((yi > y) !== (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi) a = !a;
    }
    return a;
  }
  // ------------------------------------------------------------ 谷筋の向き・天然ダムの形・堰堤の長さ
  // 谷筋の向き：3 セル上流（集水の大きい方）から 3 セル下流へ。[上流のセル, 下流のセル]
  function streamEnds(c) {
    var D = ST.D, A = ST.A, dn = c, up = c, k;
    for (k = 0; k < 3; k++) { var d = SD.downOf(D, A.dir, dn); if (d < 0 || (ST.mask && !ST.mask[d])) break; dn = d; }
    for (k = 0; k < 3; k++) {
      var x0 = up % D.nx, y0 = (up / D.nx) | 0, bu = -1, ba = 0;
      for (var dy = -1; dy <= 1; dy++) for (var dx = -1; dx <= 1; dx++) {
        if (!dx && !dy) continue; var xx = x0 + dx, yy = y0 + dy; if (xx < 0 || yy < 0 || xx >= D.nx || yy >= D.ny) continue;
        var o = yy * D.nx + xx; if (SD.downOf(D, A.dir, o) === up && A.acc[o] > ba) { ba = A.acc[o]; bu = o; }
      }
      if (bu < 0) break; up = bu;
    }
    if (up === dn) { var q = A.dir[c]; return q >= 0 ? [c, -1, SD.DC[q], SD.DR[q]] : [c, -1, 1, 0]; }
    var ux = dn % D.nx - up % D.nx, uy = ((dn / D.nx) | 0) - ((up / D.nx) | 0), L = Math.hypot(ux, uy) || 1;
    return [up, dn, ux / L, uy / L];
  }
  // 解析座標の流下方向（DEM の向きを計算領域の回転に合わせる）
  function solverDir(dom, c) {
    var e = streamEnds(c);
    if (e[1] >= 0) { var a = dom.toXY(e[0]), b = dom.toXY(e[1]), vx = b[0] - a[0], vy = b[1] - a[1], L = Math.hypot(vx, vy) || 1; return [vx / L, vy / L]; }
    return dom.dirXY(c);
  }
  // 天然ダムの形（ソルバーと同じ式を DEM の上で）。台形の縦断：天端の標高 = 谷底 + H、天端幅 Wc（流れの向き）、
  // 上流・下流の法の角度。谷を横切る向きは幅 W（0 なら max(100, 16H)）まで、地盤より高い所だけ
  function ldShape(r) {
    var D = ST.D, key = [r.c, r.H, r.Wc, r.up, r.dn, r.W].join(',');
    if (r.sh && r.sh.key === key) return r.sh;
    var e = streamEnds(r.c), ux = e[2], uy = e[3], cs = D.cs, x0 = r.c % D.nx, y0 = (r.c / D.nx) | 0, zc = D.z[r.c];
    var H = Math.max(0.5, r.H || 10), Wc = Math.max(0, r.Wc || 0), tu = Math.tan((r.up > 0 ? r.up : 30) * Math.PI / 180), td = Math.tan((r.dn > 0 ? r.dn : 25) * Math.PI / 180);
    var lat = r.W > 0 ? r.W / 2 : Math.max(50, 8 * H), ext = (Math.max(H / tu, H / td) + Wc + lat) / cs, cells = [], th = [], V = 0, sA = 1e9, sB = -1e9, wA = 1e9, wB = -1e9;
    for (var j = Math.max(0, Math.floor(y0 - ext)); j <= Math.min(D.ny - 1, Math.ceil(y0 + ext)); j++)
      for (var i = Math.max(0, Math.floor(x0 - ext)); i <= Math.min(D.nx - 1, Math.ceil(x0 + ext)); i++) {
        var c = j * D.nx + i; if (ST.mask && !ST.mask[c]) continue;
        var px = (i - x0) * cs, py = (j - y0) * cs, sa = px * ux + py * uy, w = -px * uy + py * ux;
        if (Math.abs(w) > lat) continue;
        var zt = zc + H - (sa > Wc / 2 ? (sa - Wc / 2) * td : sa < -Wc / 2 ? (-Wc / 2 - sa) * tu : 0), t = zt - D.z[c];
        if (t > 1e-3) { cells.push(c); th.push(t); V += t * cs * cs; if (sa < sA) sA = sa; if (sa > sB) sB = sa; if (w < wA) wA = w; if (w > wB) wB = w; }
      }
    r.sh = { key: key, cells: cells, th: th, V: V, L: cells.length ? sB - sA + cs : 0, B: cells.length ? wB - wA + cs : 0, ux: ux, uy: uy };
    return r.sh;
  }
  // 体積 V に合う高さを二分法で
  function fitLdamV(r, V) {
    if (!(V > 0)) return; var a = 0.5, b = 300;
    for (var k = 0; k < 40; k++) { var m = 0.5 * (a + b); r.H = m; r.sh = null; if (ldShape(r).V < V) a = m; else b = m; }
    r.H = Math.round(0.5 * (a + b) * 10) / 10; r.sh = null;
  }
  // 堰堤の長さ（自動）：谷を横切る向きに、地盤が天端＋袖の高さを超えるまで。左右の和
  function damLen(r) {
    if (r.len > 0) return r.len;
    var D = ST.D, e = streamEnds(r.c), ux = e[2], uy = e[3], x0 = r.c % D.nx, y0 = (r.c / D.nx) | 0, zb = 1e30, wing = +($('d_wing') ? $('d_wing').value : 1) || 1;
    for (var dy = -1; dy <= 1; dy++) for (var dx = -1; dx <= 1; dx++) { var v = D.z[(y0 + dy) * D.nx + x0 + dx]; if (v < zb) zb = v; }
    var top = zb + r.H + wing, tot = 0;
    [-1, 1].forEach(function (sg) {
      for (var k = 1; k < 600; k++) {
        var s = 0.5 * k, xi = Math.round(x0 - sg * uy * s), yi = Math.round(y0 + sg * ux * s);
        if (xi < 0 || yi < 0 || xi >= D.nx || yi >= D.ny || !(D.z[yi * D.nx + xi] < top)) { tot += s * D.cs; return; }
      }
      tot += 300 * D.cs;
    });
    return Math.min(800, tot + 2 * D.cs);
  }
  function damDefault(c) { return { c: c, type: 'closed', H: 8, b: 15, len: 0, l: 2, n: 1, hb: 0 }; }
  // 例：本川の中流で、谷が深く勾配が 1.5 % 以上の点に天然ダム、その 20 セル上流の本川に一定流量の供給点
  function ldamExample() {
    var D = ST.D, A = ST.A, m = ST.mask; if (!D || !m) return;
    var accO = A.acc[ST.outlet], best = -1, bd = -1, c, k;
    for (c = 0; c < D.nx * D.ny; c++) {
      if (!m[c] || !A.stream[c] || A.acc[c] < 0.2 * accO || A.acc[c] > 0.8 * accO) continue;
      var q = A.dir[c]; if (q < 0 || SD.localSlope(D, A, c, 3) < 0.015) continue;
      var nxv = -SD.DR[q], nyv = SD.DC[q], dmin = 1e9;
      [-1, 1].forEach(function (sg) {
        var mx = -1e9;
        for (k = 1; k <= 6; k++) { var x = c % D.nx + Math.round(sg * k * nyv), y = ((c / D.nx) | 0) + Math.round(sg * k * nxv); if (x < 0 || y < 0 || x >= D.nx || y >= D.ny) continue; mx = Math.max(mx, D.z[y * D.nx + x]); }
        dmin = Math.min(dmin, mx - D.z[c]);
      });
      if (dmin > bd) { bd = dmin; best = c; }
    }
    if (best < 0) { $('dem_sel').textContent = '天然ダムに向いた谷が見つかりませんでした。地図でクリックして置いてください。'; return; }
    var up = best;
    for (k = 0; k < 20; k++) {
      var x0 = up % D.nx, y0 = (up / D.nx) | 0, bu = -1, ba = 0;
      for (var dy = -1; dy <= 1; dy++) for (var dx = -1; dx <= 1; dx++) {
        if (!dx && !dy) continue; var xx = x0 + dx, yy = y0 + dy; if (xx < 0 || yy < 0 || xx >= D.nx || yy >= D.ny) continue;
        var o = yy * D.nx + xx; if (SD.downOf(D, A.dir, o) === up && A.acc[o] > ba) { ba = A.acc[o]; bu = o; }
      }
      if (bu < 0) break; up = bu;
    }
    snap(); ST.ldams = [{ c: best, H: Math.max(5, Math.min(12, Math.round(bd * 0.3))), Wc: 10, up: 30, dn: 25, W: 0 }];
    ST.src = [{ c: up, mode: 'const', Qb: 40, Qp: 40, t0: 0, tr: 60, tf: 300, C: 0 }];
    $('gen_type').value = 'ldam'; genNote();
    renderTables(); draw(); selInfo();
    $('dem_sel').textContent += '　例：天然ダム（高さ ' + ST.ldams[0].H + ' m）と、その上流の本川に一定流量 40 m³/s の供給点を置きました（谷頭でなくても供給点にできます）。';
  }
  function slideDefaults(type) { return type === 'deep' ? { R: 80, D: 20, phi: 30, coh: 30 } : { R: 25, D: 1.5, phi: 35, coh: 3 }; }

  // ------------------------------------------------------------ 地図
  // 地形の画像（流域の選択によらない）。表示は 陰影＋段彩・赤色立体風・傾斜・段彩のみ。等高線も描ける
  // 全体は f×f セルごとに 1 画素（画素数 250 万以下）。拡大したときは見えている範囲を 1 セル 1 画素で描き直す（detail）
  function shadeStats(D) {
    if (ST.zst && ST.zst.D === D) return ST.zst;
    var zmin = Infinity, zmax = -Infinity; for (var c = 0; c < D.nx * D.ny; c++) { var v = D.z[c]; if (v === v) { if (v < zmin) zmin = v; if (v > zmax) zmax = v; } }
    var zr = Math.max(1e-6, zmax - zmin), ci = 1, nice = [1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000];
    for (var k = 0; k < nice.length; k++) { ci = nice[k]; if (zr / ci <= 30) break; }
    return (ST.zst = { D: D, zmin: zmin, zmax: zmax, zr: zr, ci: ci });
  }
  var TINT = [[0, 110, 165, 100], [0.25, 175, 200, 125], [0.5, 220, 200, 140], [0.75, 190, 150, 110], [1, 240, 236, 230]];
  function shade(D, A, x0, y0, x1, y1, f, ciIn) {
    var nx = D.nx, ny = D.ny, W = Math.ceil((x1 - x0) / f), H = Math.ceil((y1 - y0) / f);
    var style = $('dm_shade') ? $('dm_shade').value : 'hs', cont = $('dm_cont') ? $('dm_cont').checked : true;
    var st = shadeStats(D), zmin = st.zmin, zr = st.zr, ci = ciIn || st.ci, cs = D.cs * f, z = D.z, i, j, k;
    var R = Math.max(2, Math.min(30, Math.round(75 / cs))), M = R + 1, WW = W + 2 * M, HH = H + 2 * M;
    // 余白つきで間引いた標高（DEM の外は端の値）。欠測は NaN のまま
    var Z = new Float64Array(WW * HH);
    for (j = 0; j < HH; j++) {
      var y = y0 + (j - M) * f + (f >> 1); y = y < 0 ? 0 : y >= ny ? ny - 1 : y; var row = y * nx;
      for (i = 0; i < WW; i++) { var x = x0 + (i - M) * f + (f >> 1); x = x < 0 ? 0 : x >= nx ? nx - 1 : x; Z[j * WW + i] = z[row + x]; }
    }
    // 谷筋（ブロック内の最大次数）
    var S = new Uint8Array(W * H);
    if (A) for (j = y0; j < Math.min(ny, y0 + H * f); j++) for (i = x0; i < Math.min(nx, x0 + W * f); i++) { var c = j * nx + i; if (A.stream[c]) { k = ((j - y0) / f | 0) * W + ((i - x0) / f | 0); if (A.strahler[c] > S[k]) S[k] = A.strahler[c]; } }
    // TPI 用の積分画像
    var I = new Float64Array((WW + 1) * (HH + 1)), N = new Float64Array((WW + 1) * (HH + 1));
    for (j = 0; j < HH; j++) { var rs = 0, rn = 0, o1 = (j + 1) * (WW + 1), o0 = j * (WW + 1); for (i = 0; i < WW; i++) { var vv = Z[j * WW + i]; if (vv === vv) { rs += vv; rn++; } I[o1 + i + 1] = I[o0 + i + 1] + rs; N[o1 + i + 1] = N[o0 + i + 1] + rn; } }
    var img = new ImageData(W, H), px = img.data, L = Math.SQRT1_2;
    // 光源 4 方向（北西を主に）。重み × (sin, cos)
    var AZ = [[315, 0.55], [270, 0.2], [0, 0.15], [225, 0.1]].map(function (q) { var t = q[0] * Math.PI / 180; return [q[1] * L * Math.sin(t), q[1] * L * Math.cos(t), q[1] * L]; });
    var sr = 1 / (0.12 * R * cs), inv8 = 1 / (8 * cs);
    for (j = 0; j < H; j++) for (i = 0; i < W; i++) {
      var p = (j * W + i) * 4, q0 = (j + M) * WW + i + M, zz = Z[q0];
      if (zz !== zz) { px[p] = px[p + 1] = px[p + 2] = 232; px[p + 3] = 255; continue; }
      var a1 = Z[q0 - WW - 1], a2 = Z[q0 - WW], a3 = Z[q0 - WW + 1], a4 = Z[q0 - 1], a6 = Z[q0 + 1], a7 = Z[q0 + WW - 1], a8 = Z[q0 + WW], a9 = Z[q0 + WW + 1];
      if (a1 !== a1) a1 = zz; if (a2 !== a2) a2 = zz; if (a3 !== a3) a3 = zz; if (a4 !== a4) a4 = zz; if (a6 !== a6) a6 = zz; if (a7 !== a7) a7 = zz; if (a8 !== a8) a8 = zz; if (a9 !== a9) a9 = zz;
      // Horn の 3×3 勾配（gx：東向き、gy：北向き）
      var gx = ((a3 + 2 * a6 + a9) - (a1 + 2 * a4 + a7)) * inv8, gy = ((a1 + 2 * a2 + a3) - (a7 + 2 * a8 + a9)) * inv8;
      var g2 = gx * gx + gy * gy, nrm = 1 / Math.sqrt(1 + g2), hs = 0;
      for (k = 0; k < 4; k++) { var h1 = (-gx * AZ[k][0] - gy * AZ[k][1] + AZ[k][2]) * nrm; if (h1 > 0) hs += h1; }
      hs = hs > 0.72 ? 1 : hs / 0.72;
      var ii0 = i + M - R, ii1 = i + M + R + 1, jj0 = j + M - R, jj1 = j + M + R + 1;   // 半径 R の窓（余白込みの座標）
      var sm = I[jj1 * (WW + 1) + ii1] - I[jj0 * (WW + 1) + ii1] - I[jj1 * (WW + 1) + ii0] + I[jj0 * (WW + 1) + ii0];
      var sn = N[jj1 * (WW + 1) + ii1] - N[jj0 * (WW + 1) + ii1] - N[jj1 * (WW + 1) + ii0] + N[jj0 * (WW + 1) + ii0];
      var tpi = Math.tanh((zz - sm / (sn || 1)) * sr), r, g, b2;
      if (style === 'red') {   // 赤色立体風：明るさ＝尾根・谷（TPI）、赤み＝傾斜
        var slp = Math.atan(Math.sqrt(g2)) * 57.3, br = 0.5 + 0.5 * tpi, sat = slp > 50 ? 1 : slp / 50;
        r = 255 * br; g = 255 * br * (1 - 0.8 * sat); b2 = 255 * br * (1 - 0.9 * sat);
      } else if (style === 'slope') {   // 傾斜：急なほど濃い
        var sl2 = Math.atan(Math.sqrt(g2)) * 57.3, gr = 255 * (1 - (sl2 > 55 ? 1 : sl2 / 55)) * (0.8 + 0.2 * hs);
        r = g = b2 = gr;
      } else {
        var t = (zz - zmin) / zr, q = 1; while (q < 4 && t > TINT[q][0]) q++;
        var A0 = TINT[q - 1], B0 = TINT[q], u = Math.max(0, Math.min(1, (t - A0[0]) / (B0[0] - A0[0])));
        var cr = A0[1] + u * (B0[1] - A0[1]), cg = A0[2] + u * (B0[2] - A0[2]), cb = A0[3] + u * (B0[3] - A0[3]);
        if (style === 'plain') { r = cr; g = cg; b2 = cb; }
        else {   // 陰影は強め、尾根を少し明るく谷を少し暗く
          var sh = (0.12 + 1.0 * hs) * (1 + 0.22 * tpi), w = 0.55;   // 色を少し薄めて陰影を見やすく
          r = (w * cr + (1 - w) * 225) * sh; g = (w * cg + (1 - w) * 225) * sh; b2 = (w * cb + (1 - w) * 225) * sh;
        }
      }
      if (cont) {   // 等高線（5 本ごとに濃い計曲線）
        var lv = Math.floor(zz / ci), ln = Math.floor((a6 === a6 ? a6 : zz) / ci), ld = Math.floor((a8 === a8 ? a8 : zz) / ci);
        if (lv !== ln || lv !== ld) { var mx = lv > ln ? (lv > ld ? lv : ld) : (ln > ld ? ln : ld), dk = mx % 5 === 0 ? 0.45 : 0.7; r *= dk; g *= dk; b2 *= dk; }
      }
      var so = S[j * W + i]; if (so) { r = 25; g = 95 + 12 * (so > 8 ? 8 : so); b2 = 215; }
      px[p] = r; px[p + 1] = g; px[p + 2] = b2; px[p + 3] = 255;
    }
    var cv = document.createElement('canvas'); cv.width = W; cv.height = H; cv.getContext('2d').putImageData(img, 0, 0);
    return { cv: cv, f: f, W: W, H: H, ci: ci, x0: x0, y0: y0 };
  }
  function shadeAll(D, A) { var f = Math.max(1, Math.ceil(Math.sqrt(D.nx * D.ny / 2.5e6))); return shade(D, A, 0, 0, D.nx, D.ny, f); }
  // 拡大したとき、見えている範囲を細かく描き直す（動かし終わって 0.25 秒後）
  var detT = 0;
  function wantDetail() {
    clearTimeout(detT);
    var I = ST.img, v = ST.view, D = ST.D; if (!I || !v || I.f === 1 || v.sc * I.f < 1.6) return;
    detT = setTimeout(function () {
      var v2 = ST.view; if (!v2) return;
      var x0 = Math.max(0, Math.floor(v2.x0) - 2), y0 = Math.max(0, Math.floor(v2.y0) - 2), x1 = Math.min(D.nx, Math.ceil(v2.x0 + v2.W / v2.sc) + 2), y1 = Math.min(D.ny, Math.ceil(v2.y0 + v2.H / v2.sc) + 2);
      if (x1 <= x0 || y1 <= y0 || (x1 - x0) * (y1 - y0) > 3e6) return;
      var dt = ST.det; if (dt && dt.x0 <= x0 && dt.y0 <= y0 && dt.x0 + dt.W >= x1 && dt.y0 + dt.H >= y1) return;
      var mg = Math.round(0.25 * Math.max(x1 - x0, y1 - y0));   // 少し広めに描いて、小さな移動で描き直さない
      x0 = Math.max(0, x0 - mg); y0 = Math.max(0, y0 - mg); x1 = Math.min(D.nx, x1 + mg); y1 = Math.min(D.ny, y1 + mg);
      ST.det = shade(D, ST.A, x0, y0, x1, y1, 1, I.ci); ST.detM = null; draw();
    }, 250);
  }
  // 選んだ流域の外を白くかぶせる画像（画像 I と同じ範囲・間引き）
  function maskImg(D, I) {
    var W = I.W, H = I.H, f = I.f, im = new ImageData(W, H), px = im.data, m = ST.mask;
    if (m) for (var j = 0; j < H; j++) for (var i = 0; i < W; i++) {
      var c = Math.min(D.ny - 1, I.y0 + j * f + (f >> 1)) * D.nx + Math.min(D.nx - 1, I.x0 + i * f + (f >> 1));
      if (!m[c]) { var p = (j * W + i) * 4; px[p] = px[p + 1] = px[p + 2] = 255; px[p + 3] = 120; }
    }
    var cv = document.createElement('canvas'); cv.width = W; cv.height = H; cv.getContext('2d').putImageData(im, 0, 0); return cv;
  }
  // ------------------------------------------------------------ 表示範囲（ズーム・移動）
  // ST.cam = {cx, cy, sc}：中心のセル座標と 1 セルの画素数
  function mapSize() {
    var D = ST.D, Wc = $('demmap').parentNode.clientWidth || 600;
    return [Wc, Math.round(Math.min(620, Math.max(320, Wc * D.ny / D.nx)))];
  }
  function fitBox(x0, y0, x1, y1) {
    var s = mapSize(), sc = 0.98 * Math.min(s[0] / Math.max(1, x1 - x0), s[1] / Math.max(1, y1 - y0));
    ST.cam = { cx: (x0 + x1) / 2, cy: (y0 + y1) / 2, sc: sc };
  }
  function fitAll() { fitBox(0, 0, ST.D.nx, ST.D.ny); }
  function fitBasin() {
    var D = ST.D; if (!ST.mask) { fitAll(); return; }
    var x0 = D.nx, y0 = D.ny, x1 = 0, y1 = 0;
    for (var c = 0; c < ST.mask.length; c++) if (ST.mask[c]) { var x = c % D.nx, y = (c / D.nx) | 0; if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y; }
    var mg = Math.max(4, 0.06 * Math.max(x1 - x0, y1 - y0));
    fitBox(x0 - mg, y0 - mg, x1 + mg + 1, y1 + mg + 1);
  }
  function scLimits() { var s = mapSize(), D = ST.D; return [0.5 * Math.min(s[0] / D.nx, s[1] / D.ny), 40]; }
  function zoomAt(px, py, k) {
    var v = ST.view, cam = ST.cam; if (!v || !cam) return;
    var lim = scLimits(), sc = Math.max(lim[0], Math.min(lim[1], cam.sc * k));
    var gx = v.x0 + px / v.sc, gy = v.y0 + py / v.sc;   // カーソルの下のセルを動かさない
    cam.cx = gx - (px - v.W / 2) / sc; cam.cy = gy - (py - v.H / 2) / sc; cam.sc = sc; redraw();
  }
  var rq = 0; function redraw() { if (!rq) rq = requestAnimationFrame(function () { rq = 0; draw(); }); }
  function draw() {
    var cv = $('demmap'); if (!cv || $('demcard').hidden) return;
    var D = ST.D, dpr = window.devicePixelRatio || 1, Wc = cv.parentNode.clientWidth || 600;
    if (!D) { cv.width = Wc * dpr; cv.height = 120 * dpr; cv.style.height = '120px'; var c0 = cv.getContext('2d'); c0.setTransform(dpr, 0, 0, dpr, 0, 0); c0.clearRect(0, 0, Wc, 120); c0.fillStyle = '#888'; c0.font = '13px sans-serif'; c0.fillText('DEM を読み込むか「例の DEM」を押してください', 12, 60); return; }
    if (!ST.img) { ST.img = shadeAll(D, ST.A); ST.mcv = null; ST.det = null; }
    if (!ST.mcv) { ST.mcv = maskImg(D, ST.img); ST.detM = null; }
    if (!ST.cam) { if ($('dm_fit').checked) fitBasin(); else fitAll(); }
    var ms = mapSize(), Hc = ms[1]; Wc = ms[0];
    var sc = ST.cam.sc, vb = { x0: ST.cam.cx - Wc / 2 / sc, y0: ST.cam.cy - Hc / 2 / sc }; vb.x1 = vb.x0 + Wc / sc; vb.y1 = vb.y0 + Hc / sc;
    cv.width = Wc * dpr; cv.height = Hc * dpr; cv.style.height = Hc + 'px';
    var ctx = cv.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.fillStyle = '#eceff1'; ctx.fillRect(0, 0, Wc, Hc);
    // 画像 I の見えている部分を描く（DEM の外は描かない）
    function blit(I, mcv) {
      var ix1 = I.x0 + I.W * I.f, iy1 = I.y0 + I.H * I.f;
      var cx0 = Math.max(I.x0, vb.x0), cy0 = Math.max(I.y0, vb.y0), cx1 = Math.min(Math.min(D.nx, ix1), vb.x1), cy1 = Math.min(Math.min(D.ny, iy1), vb.y1);
      if (cx1 <= cx0 || cy1 <= cy0) return;
      ctx.imageSmoothingEnabled = sc * I.f < 12;
      var sx = (cx0 - I.x0) / I.f, sy = (cy0 - I.y0) / I.f, sw = (cx1 - cx0) / I.f, sh = (cy1 - cy0) / I.f;
      var dx0 = (cx0 - vb.x0) * sc, dy0 = (cy0 - vb.y0) * sc, dw = (cx1 - cx0) * sc, dh = (cy1 - cy0) * sc;
      ctx.drawImage(I.cv, sx, sy, sw, sh, dx0, dy0, dw, dh);
      ctx.drawImage(mcv, sx, sy, sw, sh, dx0, dy0, dw, dh);
    }
    blit(ST.img, ST.mcv);
    if (ST.det && sc * ST.img.f >= 1.6) { if (!ST.detM) ST.detM = maskImg(D, ST.det); blit(ST.det, ST.detM); }
    ST.view = { x0: vb.x0, y0: vb.y0, sc: sc, W: Wc, H: Hc };
    wantDetail();
    function P(c) { return [((c % D.nx) + 0.5 - vb.x0) * sc, (((c / D.nx) | 0) + 0.5 - vb.y0) * sc]; }
    // 流域界（選んだ流域の縁）
    if (ST.mask) {
      ctx.fillStyle = 'rgba(200,40,40,0.9)';
      var m = ST.mask, step = Math.max(1, Math.floor(1 / sc));
      for (var y = Math.max(0, Math.floor(vb.y0)); y < Math.min(D.ny, vb.y1); y += step) for (var x = Math.max(0, Math.floor(vb.x0)); x < Math.min(D.nx, vb.x1); x += step) {
        var c = y * D.nx + x; if (!m[c]) continue;
        if (x === 0 || y === 0 || x === D.nx - 1 || y === D.ny - 1 || !m[c - 1] || !m[c + 1] || !m[c - D.nx] || !m[c + D.nx]) ctx.fillRect((x - vb.x0) * sc, (y - vb.y0) * sc, Math.max(1.2, sc * step), Math.max(1.2, sc * step));
      }
      var po = P(ST.outlet); ctx.strokeStyle = '#c62828'; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(po[0], po[1], 7, 0, 7); ctx.stroke();
    }
    // 谷頭
    var selSet = {}; ST.src.forEach(function (r, i) { selSet[r.c] = i + 1; });
    ST.heads.forEach(function (c) {
      var p = P(c); if (p[0] < -5 || p[1] < -5 || p[0] > Wc + 5 || p[1] > Hc + 5) return;
      var s = selSet[c], hr = s ? 6 : Math.max(1.3, Math.min(3.2, 4 * Math.sqrt(sc)));   // 縮小したときは小さく（地形が隠れないように）
      ctx.beginPath(); ctx.arc(p[0], p[1], hr, 0, 7);
      ctx.fillStyle = s ? '#ef6c00' : '#fff'; ctx.fill(); ctx.lineWidth = hr < 2.5 && !s ? 0.6 : 1.2; ctx.strokeStyle = '#333'; ctx.stroke();
      if (s) { ctx.fillStyle = '#000'; ctx.font = '600 11px sans-serif'; ctx.fillText(String(s), p[0] + 7, p[1] - 6); }
    });
    // 谷頭でない供給点（本川の途中など）
    var hset = {}; ST.heads.forEach(function (c) { hset[c] = 1; });
    ST.src.forEach(function (r, i) {
      if (hset[r.c]) return; var p = P(r.c);
      ctx.beginPath(); ctx.rect(p[0] - 5, p[1] - 5, 10, 10); ctx.fillStyle = '#ef6c00'; ctx.fill(); ctx.lineWidth = 1.2; ctx.strokeStyle = '#333'; ctx.stroke();
      ctx.fillStyle = '#000'; ctx.font = '600 11px sans-serif'; ctx.fillText(String(i + 1), p[0] + 7, p[1] - 6);
    });
    // 斜面崩壊
    function polyPath(list) { ctx.beginPath(); list.forEach(function (c, k) { var q = P(c); if (k) ctx.lineTo(q[0], q[1]); else ctx.moveTo(q[0], q[1]); }); }
    // 河床材料の範囲（材料ごとの色）
    var MM = window.__saboMat || { sym: 'ABCDEFGHIJKL', col: ['#b89a6a', '#1baf7a', '#2a78d6'] };
    ST.matz.forEach(function (r) {
      var col = MM.col[r.mat] || '#1baf7a';
      polyPath(r.poly); ctx.closePath(); ctx.globalAlpha = 0.28; ctx.fillStyle = col; ctx.fill(); ctx.globalAlpha = 1;
      ctx.setLineDash([6, 3]); ctx.strokeStyle = col; ctx.lineWidth = 1.6; ctx.stroke(); ctx.setLineDash([]);
      var q = P(r.poly[0]); ctx.fillStyle = col; ctx.font = '600 11px sans-serif'; ctx.fillText('河床 ' + MM.sym[r.mat], q[0] + 4, q[1] - 4);
    });
    ST.slides.forEach(function (r, i) {
      var fillC = r.type === 'deep' ? 'rgba(121,85,72,0.45)' : 'rgba(255,152,0,0.4)', p, R;
      if (r.poly) {
        polyPath(r.poly); ctx.closePath(); ctx.fillStyle = fillC; ctx.fill(); ctx.strokeStyle = '#5d4037'; ctx.lineWidth = 1.5; ctx.stroke();
        var xm = -1e9; r.poly.forEach(function (c) { var q = P(c); if (q[0] > xm) xm = q[0]; }); p = P(r.c); R = xm - p[0];
      } else {
        p = P(r.c); R = Math.max(4, r.R / D.cs * sc);
        ctx.beginPath(); ctx.arc(p[0], p[1], R, 0, 7); ctx.fillStyle = fillC; ctx.fill(); ctx.strokeStyle = '#5d4037'; ctx.lineWidth = 1.5; ctx.stroke();
      }
      ctx.fillStyle = '#3e2723'; ctx.font = '600 11px sans-serif'; ctx.fillText('崩' + (i + 1), p[0] + R + 2, p[1]);
    });
    // 描きかけの多角形
    if (ST.draft && ST.draft.length) {
      polyPath(ST.draft); ctx.strokeStyle = '#d84315'; ctx.lineWidth = 2; ctx.setLineDash([5, 4]); if (ST.draft.length > 2) ctx.closePath(); ctx.stroke(); ctx.setLineDash([]);
      ST.draft.forEach(function (c) { var q = P(c); ctx.beginPath(); ctx.arc(q[0], q[1], 3.5, 0, 7); ctx.fillStyle = '#d84315'; ctx.fill(); });
    }
    // 天然ダム：堤体の範囲（厚いほど濃い）と番号
    ST.ldams.forEach(function (r, i) {
      var sh = ldShape(r), cw = Math.max(1, sc), Hm = Math.max(0.5, r.H);
      for (var k = 0; k < sh.cells.length; k++) {
        var q = P(sh.cells[k]); if (q[0] < -cw || q[1] < -cw || q[0] > Wc + cw || q[1] > Hc + cw) continue;
        ctx.fillStyle = 'rgba(93,64,55,' + (0.25 + 0.6 * Math.min(1, sh.th[k] / Hm)) + ')'; ctx.fillRect(q[0] - cw / 2, q[1] - cw / 2, cw, cw);
      }
      var p = P(r.c), wv = Math.max(8, sh.B / 2 / D.cs * sc);
      ctx.strokeStyle = '#3e2723'; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(p[0] + sh.uy * wv, p[1] - sh.ux * wv); ctx.lineTo(p[0] - sh.uy * wv, p[1] + sh.ux * wv); ctx.stroke();
      ctx.fillStyle = '#3e2723'; ctx.font = '600 11px sans-serif'; ctx.fillText('天然ダム' + (i + 1) + '（' + fmt(r.H, 1) + ' m, ' + fmt(sh.V / 1e4, 1) + ' 万 m³）', p[0] + wv * Math.abs(sh.uy) + 6, p[1] + 4);
    });
    // 砂防堰堤：谷を横切る線（濃い：袖、白：水通し）
    ST.dams.forEach(function (r, i) {
      var e = streamEnds(r.c), ux = e[2], uy = e[3], p = P(r.c), hl = damLen(r) / 2 / D.cs * sc, hb = Math.min(hl, (r.b || 10) / 2 / D.cs * sc), lw = Math.max(4, 1.6 * sc);
      ctx.lineCap = 'butt';
      ctx.strokeStyle = '#37474f'; ctx.lineWidth = lw; ctx.beginPath(); ctx.moveTo(p[0] + uy * hl, p[1] - ux * hl); ctx.lineTo(p[0] - uy * hl, p[1] + ux * hl); ctx.stroke();
      ctx.strokeStyle = r.type === 'closed' ? '#eceff1' : '#80deea'; ctx.lineWidth = Math.max(2, lw - 2.5); ctx.beginPath(); ctx.moveTo(p[0] + uy * hb, p[1] - ux * hb); ctx.lineTo(p[0] - uy * hb, p[1] + ux * hb); ctx.stroke();
      ctx.fillStyle = '#263238'; ctx.font = '600 11px sans-serif'; ctx.fillText('堰' + (i + 1) + '（' + fmt(r.H, 1) + ' m）', p[0] + Math.abs(uy) * hl + 6, p[1] - 6);
    });
    $('dem_key').innerHTML = '<span><i style="border-color:#1e6ec8"></i>谷筋</span><span><i class="b" style="background:#fff;border:1px solid #333"></i>谷頭</span><span><i class="b" style="background:#ef6c00"></i>供給点</span><span><i style="border-color:#c62828"></i>流域界・出口</span><span><i class="b" style="background:rgba(255,152,0,0.5)"></i>斜面崩壊</span><span><i class="b" style="background:#5d4037"></i>天然ダム（濃いほど厚い）</span><span><i style="border-color:#37474f;border-width:3px"></i>砂防堰堤（白：水通し）</span><span><i class="b" style="background:rgba(27,175,122,0.4)"></i>河床材料の範囲（材料ごとの色）</span>' + ($('dm_cont').checked ? '<span>等高線 ' + ST.img.ci + ' m 間隔（太線は ' + 5 * ST.img.ci + ' m ごと）</span>' : '') + '<span style="color:#666">ホイールで拡大・縮小、ドラッグで移動、右クリックで近くの物を外す</span>';
  }
  function cellAt(e) {
    var cv = $('demmap'), r = cv.getBoundingClientRect(), v = ST.view; if (!v || !ST.D) return -1;
    var x = Math.floor(v.x0 + (e.clientX - r.left) / v.sc), y = Math.floor(v.y0 + (e.clientY - r.top) / v.sc);
    if (x < 0 || y < 0 || x >= ST.D.nx || y >= ST.D.ny) return -1;
    return y * ST.D.nx + x;
  }
  function nearest(c, list, pxMax) {
    var D = ST.D, x = c % D.nx, y = (c / D.nx) | 0, best = -1, bd = Infinity, lim = pxMax / ST.view.sc;
    list.forEach(function (o) { var dx = o % D.nx - x, dy = ((o / D.nx) | 0) - y, d = Math.sqrt(dx * dx + dy * dy); if (d < bd) { bd = d; best = o; } });
    return bd <= lim ? best : -1;
  }
  // クリックした点の近くの谷筋のセル。画面上の距離（px）から、次数 1 ごとに 1.5 px、DEM の端へ出る末端の出口なら 4 px を引いた値が最小のもの
  // （本川の近くをクリックしたとき細い支流に吸い寄せられず、末端の出口の近くなら出口そのものを選ぶ）
  function nearStream(c, r) {
    var D = ST.D, A = ST.A, x = c % D.nx, y = (c / D.nx) | 0, best = -1, bs = Infinity, sc = ST.view ? ST.view.sc : 1;
    for (var dy = -r; dy <= r; dy++) for (var dx = -r; dx <= r; dx++) {
      var xx = x + dx, yy = y + dy; if (xx < 0 || yy < 0 || xx >= D.nx || yy >= D.ny) continue;
      var o = yy * D.nx + xx; if (!A.stream[o]) continue;
      var sco = Math.sqrt(dx * dx + dy * dy) * sc - 1.5 * A.strahler[o] - (SD.downOf(D, A.dir, o) < 0 ? 4 : 0);
      if (sco < bs || (sco === bs && A.acc[o] > A.acc[best])) { bs = sco; best = o; }
    }
    return best;
  }
  function newSlide(c, poly) {
    var tp = $('gen_type').value === 'ldam' ? 'deep' : 'shallow', d = slideDefaults(tp);
    var r = { c: c, type: tp, R: d.R, D: d.D, trigger: $('rn_on').checked ? 'fs' : 'time', t: 30, move: tp === 'deep' ? 'block' : 'auto', phi: d.phi, coh: d.coh, m0: $('rn_on').checked ? 0.2 : 0.9 };
    if (poly) { r.poly = poly; r.R = Math.round(Math.sqrt(polyArea(poly) / Math.PI)); }
    return r;
  }
  // 多角形の面積 [m²]（頂点は DEM のセル番号）
  function polyArea(poly) {
    var D = ST.D, a = 0;
    for (var i = 0, j = poly.length - 1; i < poly.length; j = i++) {
      var xi = poly[i] % D.nx, yi = (poly[i] / D.nx) | 0, xj = poly[j] % D.nx, yj = (poly[j] / D.nx) | 0;
      a += xj * yi - xi * yj;
    }
    return Math.abs(a) / 2 * D.cs * D.cs;
  }
  function polyMsg() {
    var n = (ST.draft || []).length;
    $('dem_sel').textContent = n < 3 ? '多角形の頂点をクリックしてください（いま ' + n + ' 点。3 点以上で閉じられます）。' : '頂点 ' + n + ' 点、面積 ' + fmt(polyArea(ST.draft), 0) + ' m²。ダブルクリックか「多角形を閉じる」で確定します。';
  }
  function closePoly() {
    var dr = ST.draft || []; if (dr.length < 3) { polyMsg(); return; }
    var D = ST.D, sx = 0, sy = 0; dr.forEach(function (c) { sx += c % D.nx; sy += (c / D.nx) | 0; });
    var cc = Math.round(sy / dr.length) * D.nx + Math.round(sx / dr.length);
    if (!ST.mask[cc]) cc = dr[0];
    snap();
    if ($('dm_mode').value === 'mat') { if (!window.__saboMat || window.__saboMat.n() < 2) { ST.hist.pop(); $('dem_sel').textContent = '先に「粒径と河床材料」で「＋ 河床材料を追加」をして、材料 B などの割合を入れてください。'; ST.draft = []; draw(); return; } ST.matz.push({ poly: dr.slice(), mat: +$('dm_matsel').value || 1 }); }
    else ST.slides.push(newSlide(cc, dr.slice()));
    ST.draft = [];
    renderTables(); draw(); selInfo();
  }
  function shapeUI() {
    var on = $('dm_mode').value === 'slide', mt = $('dm_mode').value === 'mat', pg = (on && $('dm_shape').value === 'poly') || mt;
    $('dm_shapew').hidden = !on; $('dm_matw').hidden = !mt; $('dm_close').hidden = !pg; $('dm_undo').hidden = !pg;
    if (!pg && ST.draft && ST.draft.length) { ST.draft = []; draw(); }
    if (pg) polyMsg();
  }
  function onClick(e) {
    if (ST.dragged) { ST.dragged = false; return; }
    var c = cellAt(e); if (c < 0) return;
    var mode = $('dm_mode').value, A = ST.A;
    if (mode === 'basin') { if (A.bid[c] >= 0) selectBasin(A.bid[c], false); return; }
    if (mode === 'outlet') { var s = nearStream(c, Math.max(2, Math.min(40, Math.round(10 / ST.view.sc)))); if (s >= 0) selectOutlet(s); else $('dem_sel').textContent = '谷筋（青い線）の近くをクリックしてください。'; return; }
    if (mode === 'del') { var nm = removeAt(c); $('dem_sel').textContent = nm ? nm + ' を外しました（「元に戻す」で戻せます）。' : '近くに外せる物がありません。'; renderTables(); draw(); return; }
    if (!ST.mask || !ST.mask[c]) { $('dem_sel').textContent = '選んだ流域の中をクリックしてください。'; return; }
    if (!(mode === 'slide' && $('dm_shape').value === 'poly') && mode !== 'mat') snap();
    if (mode === 'src') {
      var h = nearest(c, ST.heads, 14); if (h < 0) { ST.hist.pop(); return; }
      var k = -1; ST.src.forEach(function (r, i) { if (r.c === h) k = i; });
      if (k >= 0) ST.src.splice(k, 1); else ST.src.push(srcDefault(h));
    } else if (mode === 'mat') {
      ST.hist.pop(); ST.draft = ST.draft || [];
      if (ST.draft[ST.draft.length - 1] !== c) ST.draft.push(c);
      polyMsg(); draw(); return;
    } else if (mode === 'slide') {
      if ($('dm_shape').value === 'poly') {
        ST.draft = ST.draft || [];
        if (ST.draft[ST.draft.length - 1] !== c) ST.draft.push(c);   // ダブルクリックの 2 回目は足さない
        polyMsg(); draw(); return;
      }
      ST.slides.push(newSlide(c, null));
    } else if (mode === 'ldam') {
      var s2 = nearStream(c, Math.max(2, Math.min(40, Math.round(10 / ST.view.sc)))); if (s2 < 0) { ST.hist.pop(); $('dem_sel').textContent = '谷筋（青い線）の近くをクリックしてください。'; return; }
      ST.ldams.push({ c: s2, H: 15, Wc: 10, up: 30, dn: 25, W: 0 });
    } else if (mode === 'dam') {
      var s3 = nearStream(c, Math.max(2, Math.min(40, Math.round(10 / ST.view.sc)))); if (s3 < 0) { ST.hist.pop(); $('dem_sel').textContent = '谷筋（青い線）の近くをクリックしてください。'; return; }
      ST.dams.push(damDefault(s3));
    }
    renderTables(); draw(); selInfo();
  }
  function onMove(e) {
    var c = cellAt(e), tip = $('dem_tip'); if (c < 0 || !ST.A) { tip.hidden = true; return; }
    var D = ST.D, A = ST.A, a = A.acc[c] * D.cs * D.cs, sl = SD.localSlope(D, A, c, 2);
    tip.hidden = false; tip.style.left = (e.offsetX + 14) + 'px'; tip.style.top = (e.offsetY + 10) + 'px';
    tip.textContent = '標高 ' + fmt(D.z[c], 1) + ' m　集水面積 ' + (a >= 1e5 ? fmt(a / 1e6, 3) + ' km²' : fmt(a, 0) + ' m²') + '　勾配 ' + fmt(sl, 3) + (A.stream[c] ? '　谷筋（次数 ' + A.strahler[c] + '）' : '');
  }

  // ------------------------------------------------------------ 雨
  function lossOf() {
    var t = $('rn_loss').value;
    if (t === 'f1rsa') return { type: 'f1rsa', f1: +$('rn_f1').value, Rsa: +$('rn_rsa').value, R0: +$('rn_r0').value };
    if (t === 'ks') return { type: 'ks', ks: +$('rn_ks').value, D: +$('rn_D').value, n: +$('rn_n').value, th0: +$('rn_th').value };
    return { type: 'f', f: +$('rn_f').value };
  }
  function rainTable() {
    var tab = [];
    $('rn_tab').value.split(/\n/).forEach(function (ln) { var a = ln.split(/[,\s\t]+/).filter(function (s) { return s !== ''; }).map(Number); if (a.length >= 2 && isFinite(a[0]) && isFinite(a[1])) tab.push([a[0], Math.max(0, a[1])]); });
    tab.sort(function (a, b) { return a[0] - b[0]; });
    return tab;
  }
  // 格子の雨の CSV（t[分], X, Y, mm/h）→ 時刻ごとに DEM のセルへ（規則格子の点とみなして最も近い点）
  function parseRainCSV(text) {
    var D = ST.D, rows = [];
    text.split(/\r?\n/).forEach(function (ln) { var a = ln.split(/[,\s\t]+/).filter(function (s) { return s !== ''; }).map(Number); if (a.length >= 4 && a.every(isFinite)) rows.push(a); });
    if (!rows.length) throw new Error('雨の CSV に「時刻, X, Y, 雨量」の行がありません');
    var geo = D.src === 'gsi' && D.geo, mx = 1, my = 1;
    if (geo) { var latm = (geo.lat0 + geo.lat1) / 2 * Math.PI / 180; my = 111132.92 - 559.82 * Math.cos(2 * latm) + 1.175 * Math.cos(4 * latm); mx = 111412.84 * Math.cos(latm) - 93.5 * Math.cos(3 * latm); }
    var xs = rows.map(function (r) { return r[1] * mx; }), ys = rows.map(function (r) { return r[2] * my; });
    function stepOf(v) { var u = v.slice().sort(function (a, b) { return a - b; }), d = []; for (var i = 1; i < u.length; i++) if (u[i] - u[i - 1] > 1e-6) d.push(u[i] - u[i - 1]); d.sort(function (a, b) { return a - b; }); return { min: u[0], d: d.length ? d[d.length >> 1] : 1 }; }
    var SX = stepOf(xs), SY = stepOf(ys), byT = {};
    rows.forEach(function (r, i) { var k = String(r[0]); (byT[k] || (byT[k] = {}))[Math.round((xs[i] - SX.min) / SX.d) + ',' + Math.round((ys[i] - SY.min) / SY.d)] = r[3]; });
    var ts = Object.keys(byT).map(Number).sort(function (a, b) { return a - b; }), N = D.nx * D.ny;
    var gi = new Int32Array(N), gj = new Int32Array(N);
    for (var c = 0; c < N; c++) { var x = D.x0 + ((c % D.nx) + 0.5) * D.cs, y = D.yTop - (((c / D.nx) | 0) + 0.5) * D.cs; gi[c] = Math.round((x - SX.min) / SX.d); gj[c] = Math.round((y - SY.min) / SY.d); }
    var frames = ts.map(function (t) {
      var M = byT[String(t)], r = new Float32Array(N), hit = 0;
      for (var c2 = 0; c2 < N; c2++) { var v = M[gi[c2] + ',' + gj[c2]]; if (v != null) { r[c2] = v; hit++; } }
      return { t: t, r: r, hit: hit };
    });
    var cover = frames.length ? frames[0].hit / N : 0;
    return { frames: frames, info: '時刻 ' + ts.length + ' 個（' + ts[0] + '〜' + ts[ts.length - 1] + ' 分）、点の間隔 ' + fmt(SX.d, 0) + ' × ' + fmt(SY.d, 0) + ' m、DEM のセルの ' + fmt(cover * 100, 0) + ' % に値がある' };
  }
  function rainSpec(forSolverDom) {
    var mode = $('rn_mode').value, R = { on: $('rn_on').checked, method: $('rn_meth').value, mode: mode, table: rainTable(), frames: [], loss: lossOf() };
    if (mode === 'grid') {
      if (!ST.rainFrames) throw new Error('格子の雨の CSV を読み込んでください');
      R.frames = ST.rainFrames;
      if (forSolverDom) {   // 計算領域の格子（切り出して回した DEM）へ
        var dom = forSolverDom, nr = dom.ny, nc = dom.nx;
        R.frames = ST.rainFrames.map(function (f) {
          var r = new Array(nr * nc);
          for (var j = 0; j < nr; j++) for (var i = 0; i < nc; i++) { var oc = dom.cellOf[j][i]; r[j * nc + i] = oc >= 0 ? +f.r[oc].toFixed(2) : 0; }
          return { t: f.t, r: r };
        });
      }
    }
    return R;
  }
  function lateral(durS) {
    var R = rainSpec(null), rt = $('rn_route').value;
    var route = rt === 'storage' ? { type: 'storage', K: +$('rn_K').value, P: +$('rn_P').value, Tl: +$('rn_tl').value } : { type: 'rational', tc: +$('rn_tc').value };
    return SH.lateralInflow(SD, ST.D, ST.A, ST.mask, R, { loss: R.loss, route: route, dt: 60, tEnd: durS });
  }
  function rainChart(L) {
    var w = $('rn_chartw'); w.hidden = false;
    var cv = $('rn_chart'), dpr = window.devicePixelRatio || 1, W = w.clientWidth || 600, H = 200;
    cv.width = W * dpr; cv.height = H * dpr; cv.style.height = H + 'px';
    var ctx = cv.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, W, H);
    var n = L.n, Qmax = 0, k; for (k = 0; k < n; k++) if (L.Qsum[k] > Qmax) Qmax = L.Qsum[k];
    var tab = rainTable(), rmax = 0; tab.forEach(function (r) { if (r[1] > rmax) rmax = r[1]; });
    var x0 = 48, x1 = W - 44, y0 = 14, y1 = H - 28, T = (n - 1) * L.dt;
    ctx.strokeStyle = '#999'; ctx.lineWidth = 1; ctx.strokeRect(x0, y0, x1 - x0, y1 - y0);
    if ($('rn_mode').value === 'uniform' && rmax > 0) {   // 雨（上から下向きの棒）
      ctx.fillStyle = 'rgba(30,110,200,0.25)';
      for (k = 0; k < tab.length; k++) { var ta = tab[k][0] * 60, tb = k + 1 < tab.length ? tab[k + 1][0] * 60 : ta + 600; if (ta > T) break; var xa = x0 + (x1 - x0) * ta / T, xb = x0 + (x1 - x0) * Math.min(T, tb) / T; ctx.fillRect(xa, y0, xb - xa, (y1 - y0) * 0.45 * tab[k][1] / rmax); }
    }
    ctx.strokeStyle = '#c62828'; ctx.lineWidth = 2; ctx.beginPath();
    for (k = 0; k < n; k++) { var x = x0 + (x1 - x0) * k / (n - 1), y = y1 - (y1 - y0) * L.Qsum[k] / Math.max(1e-9, Qmax * 1.1); if (k) ctx.lineTo(x, y); else ctx.moveTo(x, y); }
    ctx.stroke();
    ctx.fillStyle = '#444'; ctx.font = '11px sans-serif';
    ctx.fillText(fmt(Qmax * 1.1, 1) + ' m³/s', 2, y0 + 10); ctx.fillText('0', x0 - 12, y1); ctx.fillText(fmt(T / 60, 0) + ' 分', x1 - 20, H - 10); ctx.fillText('0', x0, H - 10);
    if (rmax > 0) { ctx.fillStyle = 'rgba(30,110,200,0.9)'; ctx.fillText(fmt(rmax, 0) + ' mm/h', x1 + 2, y0 + 10); }
    ctx.fillStyle = '#c62828'; ctx.fillText('谷筋へ入る流量の合計（河道の流れの遅れは含まない）', x0 + 6, y1 - 6);
  }
  function calcRain() {
    try {
      if (!ST.mask) throw new Error('流域を選んでください');
      var dur = Math.max(1, +$('dem_dur').value || 30) * 60, L = lateral(dur);
      ST.lat = L; rainChart(L);
      var Vq = 0; for (var k = 0; k < L.n; k++) Vq += L.Qsum[k] * L.dt;
      var pk = 0, tk = 0; for (k = 0; k < L.n; k++) if (L.Qsum[k] > pk) { pk = L.Qsum[k]; tk = k * L.dt; }
      $('rn_info').textContent = '谷筋のセル ' + L.cells.length + '、流域面積 ' + fmt(L.area0 / 1e6, 3) + ' km²。計算時間内の雨 ' + fmt(L.rainVol, 0) + ' m³、谷へ入る水 ' + fmt(Vq, 0) + ' m³（' + fmt(100 * Vq / Math.max(1, L.rainVol), 0) + ' %）、ピーク ' + fmt(pk, 2) + ' m³/s（' + fmt(tk / 60, 0) + ' 分）。到達時間（土研式）は ' + fmt(Math.min.apply(null, Array.from(L.Tc)) / 60, 0) + '〜' + fmt(Math.max.apply(null, Array.from(L.Tc)) / 60, 0) + ' 分。';
    } catch (e) { $('rn_info').textContent = e.message; }
  }

  // ------------------------------------------------------------ 計算条件へ
  // ------------------------------------------------------------ 縦断の線（本川・支川）
  // 本川：出口から集水の大きい方へ谷頭まで。支川：本川（と 1 次の支川）に合流する谷筋で、集水面積が流域の 3 % 以上のもの。
  // 上流 → 下流の順に、計算領域の座標で返す。本川は出口から右端まで（出口の水路）を足す
  function upMain(c) {
    var D = ST.D, A = ST.A, x0 = c % D.nx, y0 = (c / D.nx) | 0, bu = -1, ba = 0;
    for (var dy = -1; dy <= 1; dy++) for (var dx = -1; dx <= 1; dx++) {
      if (!dx && !dy) continue; var xx = x0 + dx, yy = y0 + dy; if (xx < 0 || yy < 0 || xx >= D.nx || yy >= D.ny) continue;
      var o = yy * D.nx + xx; if (ST.mask[o] && SD.downOf(D, A.dir, o) === c && A.acc[o] > ba) { ba = A.acc[o]; bu = o; }
    }
    return bu;
  }
  function traceUp(c) { var p = [c], k = 0; while (k++ < 1e6) { var u = upMain(c); if (u < 0) break; p.push(u); c = u; } return p; }   // 下流 → 上流
  function buildPaths(dom) {
    var D = ST.D, A = ST.A, a1 = D.cs * D.cs, aB = A.acc[ST.outlet] * a1, minA = Math.max(0.03 * aB, 4 * (parseFloat($('dem_thr').value) || 0.02) * 1e6);
    var main = traceUp(ST.outlet), out = [], s0 = 0;
    function addTribs(path, level, parent) {   // path：下流 → 上流（支川は先頭が合流点。合流点は見ない）
      var dist = 0;
      for (var k = 0; k < path.length; k++) {
        if (level > 1 && k === 0) continue;
        var c = path[k], nxt = path[k + 1], x0 = c % D.nx, y0 = (c / D.nx) | 0;
        if (k) dist += Math.hypot(c % D.nx - path[k - 1] % D.nx, ((c / D.nx) | 0) - ((path[k - 1] / D.nx) | 0)) * D.cs;
        for (var dy = -1; dy <= 1; dy++) for (var dx = -1; dx <= 1; dx++) {
          if (!dx && !dy) continue; var xx = x0 + dx, yy = y0 + dy; if (xx < 0 || yy < 0 || xx >= D.nx || yy >= D.ny) continue;
          var o = yy * D.nx + xx; if (o === nxt || !ST.mask[o] || SD.downOf(D, A.dir, o) !== c || A.acc[o] * a1 < minA) continue;
          var tp = traceUp(o); tp.unshift(c);
          var t = { cells: tp, area: A.acc[o] * a1, level: level, at: dist, parent: parent }; out.push(t);
          if (level < 2) addTribs(tp, level + 1, t);
        }
      }
    }
    addTribs(main, 1, null);
    out.sort(function (a, b) { return b.area - a.area; }); out = out.slice(0, 20);
    var paths = [{ name: '本川（流域 ' + fmt(aB / 1e6, aB < 1e6 ? 3 : 2) + ' km²）', cells: main }];
    out.forEach(function (t, i) { t.no = i + 1; });
    out.forEach(function (t) {
      var pn = t.parent ? (t.parent.no ? '支川 ' + t.parent.no : '支川') : '本川';
      t.name = '支川 ' + t.no + '（流域 ' + fmt(t.area / 1e6, t.area < 1e6 ? 3 : 2) + ' km²、' + pn + 'の下流端から ' + fmt(t.at, 0) + ' m で合流）'; paths.push(t);
    });
    ST.paths = paths.map(function (p, q) {
      var pts = p.cells.slice().reverse().map(function (c) { return dom.toXY(c); });   // 上流 → 下流
      if (q === 0) { var e = pts[pts.length - 1]; pts.push([dom.nx * dom.cs - 0.5 * dom.cs, e[1]]); }   // 出口の水路
      return { name: p.name, pts: pts };
    });
  }
  function params(p) {
    if (!ST.D || !ST.mask) throw new Error('DEM を読み込み、流域を選んでください');
    var D = ST.D, A = ST.A;
    if (!(p.dx >= D.cs)) p.dx = D.cs;
    var dom = SD.prepareDomain(D, A, ST.mask, ST.outlet, { fill: $('dem_fill').checked, dx: p.dx });
    ST.dom = dom;
    try { buildPaths(dom); } catch (e) { ST.paths = null; }
    p.terrain = { rows: dom.rows, cs: dom.cs, basin: dom.basin, chan: dom.chan, corr: dom.corr, outXY: dom.outXY, Iout: dom.Iout };
    p.csvSlope = 0; p.L = dom.nx * dom.cs; p.sideCells = 0; p.dams = []; p.tribs = []; p.prefill = 0;
    var dur = Math.max(1, +$('dem_dur').value || 30) * 60;
    p.dem = { soilD: Math.max(0.1, +$('dem_soil').value || 1), chW: Math.max(1, +$('dem_chw').value || 10), ks: Math.max(0, +$('dem_ks').value || 0), s0: 0.3, duration: dur, outlet: $('dem_out').value, smooth: ST.D && ST.D.smooth || 0, fill: $('dem_fill').checked, edgeOut: $('dem_edge').checked, pit: Math.max(0, +$('dem_pit').value || 0) };
    var gen = $('gen_type').value;
    if (gen !== 'bed') p.transport = 'debris';
    p.sources = ST.src.map(function (r) {
      var xy = dom.toXY(r.c), u = dom.dirXY(r.c), H;
      if (r.mode === 'debris') H = { mode: 'debris', Qbase: r.Qb, Qpeak: Math.max(r.Qb, r.Qp), t0: r.t0, tRise: Math.max(1, r.tr), tFall: Math.max(1, r.tf), Cd: Math.min(0.55, Math.max(0, r.C)), Cs: 0 };
      else if (r.mode === 'const') H = { mode: 'const', Qbase: r.Qp, Qpeak: r.Qp, Cs: Math.max(0, r.C), duration: 1e12, tPeak: 1 };
      else {   // 三角形：開始までは基底、立上り・減衰で三角形、そのあと基底
        H = { mode: 'table', table: [[0, r.Qb, r.Qb * r.C], [Math.max(0, r.t0), r.Qb, r.Qb * r.C], [r.t0 + Math.max(1, r.tr), r.Qp, r.Qp * r.C], [r.t0 + r.tr + Math.max(1, r.tf), r.Qb, r.Qb * r.C], [1e9, r.Qb, r.Qb * r.C]] };
      }
      return { x: xy[0], y: xy[1], ux: u[0], uy: u[1], I: SD.localSlope(D, A, r.c, 3), hydro: H };
    });
    p.rain = { on: false }; p.latIn = null;
    if ($('rn_on').checked) {
      var R = rainSpec($('rn_meth').value === 'grid' ? dom : null);
      if (R.method === 'grid') p.rain = R;
      else {
        var L = lateral(dur); ST.lat = L;
        p.latIn = { dt: L.dt, pts: L.cells.map(function (c) { var xy = dom.toXY(c), u = dom.dirXY(c); return [xy[0], xy[1], u[0], u[1], SD.localSlope(D, A, c, 3)]; }), Q: L.Q.map(function (a) { return Array.prototype.map.call(a, function (v) { return +v.toFixed(5); }); }) };
        p.rain = { on: true, method: 'lumped', mode: R.mode, table: R.table, frames: [], loss: R.loss };   // 斜面崩壊のタンク用に雨を残す（谷へは latIn で入れる）
        if (R.mode === 'grid') p.rain.frames = rainSpec(dom).frames;
      }
    }
    var thD = +$('sl_dep').value || 10, mu = +$('sl_mu').value || 0, CL = +$('sl_cl').value || 0.5;
    p.slides = ST.slides.map(function (r) {
      var xy = dom.toXY(r.c);
      return { x: xy[0], y: xy[1], poly: r.poly ? r.poly.map(function (c) { return dom.toXY(c); }) : undefined, R: r.R, D: r.D, type: r.type, trigger: r.trigger, t: r.t, move: r.move, phi: r.phi, c: r.coh, m0: r.m0, thDep: thD, mu: mu, CL: CL, Tk: r.type === 'deep' ? +$('sl_tk2').value : +$('sl_tk1').value };
    });
    p.ldams = ST.ldams.map(function (r) { var xy = dom.toXY(r.c), u = solverDir(dom, r.c); return { x: xy[0], y: xy[1], ux: u[0], uy: u[1], H: r.H, Wc: r.Wc, up: r.up, dn: r.dn, lat: r.W > 0 ? r.W / 2 : 0 }; });
    // 場所ごとの河床材料（多角形）。DEM では x の区間の表は使わない
    p.bedZones = ST.matz.map(function (r) { return { mat: r.mat, poly: r.poly.map(function (c) { return dom.toXY(c); }) }; });
    // 砂防堰堤：中心・流下方向・長さを DEM から。水通しは中心に固定（so = 0）
    p.dams = ST.dams.map(function (r) {
      var xy = dom.toXY(r.c), u = solverDir(dom, r.c), d = { x: xy[0], y: xy[1], ang: Math.atan2(u[1], u[0]) * 180 / Math.PI, H: r.H, b: Math.max(dom.cs, r.b || 10), so: 0, len: damLen(r) };
      if (r.type !== 'closed') { d.type = r.type; d.l = r.l > 0 ? r.l : 1; d.n = Math.max(1, Math.round(r.n || 1)); if (r.hb > 0) d.hb = r.hb; }
      return d;
    });
    p.ldam = { k: Math.max(0, +$('ld_k').value || 0), fail: $('ld_fail').value, slipLevel: +$('ld_lev').value || 0.9, slipDeg: 0, thDep: thD, phi: +$('ld_phi').value || 35, c: +$('ld_c').value || 0, gs: 20, Ts: 0.5, s0: 0 };
    if (p.ldams.length) p.gw.on = true;
    p.debris.erosion = (p.ldams.length || ST.slides.some(function (r) { return r.move !== 'fluid'; })) ? $('ld_ero').value : 'sat';
    p.debris.Kt = +$('ld_kt').value || 0.0019;
    return { Bmax: dom.ny * dom.cs, dom: dom };
  }
  function genNote() {
    var g = $('gen_type').value;
    $('gen_note').textContent = g === 'bed' ? '谷頭の供給点から水（または土石流）を流し、渓床の侵食で土石流が発達する場合です。流砂モデルを「土石流」にすると、高橋の侵食・堆積速度式で発達・減衰します。'
      : g === 'slide' ? '斜面崩壊の土が流動化して土石流になる場合です。地図の「崩壊を置く」で崩壊の場所を選び、雨を降らせて安全率で崩すか、時刻で崩します。流砂モデルは自動で「土石流」になります。'
      : '斜面崩壊（土塊）か、直接置いた天然ダムが谷をふさぎ、湛水したあと決壊する場合です。決壊の形（越流侵食・進行性崩壊）は堤体の透水係数で決まり、すべり破壊は湛水位で指定できます。伏流は自動で有効になります。';
    if (g === 'slide') $('dm_mode').value = 'slide';
    if (g === 'ldam') $('dm_mode').value = 'ldam';
    shapeUI();
  }
  function rainSwitch() {
    $('rn_box').hidden = !$('rn_on').checked;
    var lt = $('rn_loss').value;
    Array.prototype.forEach.call(document.querySelectorAll('#rn_box .lf'), function (e) { e.hidden = lt !== 'f'; });
    Array.prototype.forEach.call(document.querySelectorAll('#rn_box .lr'), function (e) { e.hidden = lt !== 'f1rsa'; });
    Array.prototype.forEach.call(document.querySelectorAll('#rn_box .lk'), function (e) { e.hidden = lt !== 'ks'; });
    var rt = $('rn_route').value;
    Array.prototype.forEach.call(document.querySelectorAll('#rn_box .rr'), function (e) { e.hidden = rt !== 'rational'; });
    Array.prototype.forEach.call(document.querySelectorAll('#rn_box .rs'), function (e) { e.hidden = rt !== 'storage'; });
    $('rn_lump').hidden = $('rn_meth').value !== 'lumped';
    $('rn_uni').hidden = $('rn_mode').value !== 'uniform'; $('rn_grd').hidden = $('rn_mode').value !== 'grid';
  }
  function show(on) {
    $('g_dem').hidden = !on; $('demcard').hidden = !on; $('dem_more').hidden = !on; ['dem_dams'].forEach(function (k) { if ($(k)) $(k).hidden = !on; }); ['hyd_sec', 'dam_sec', 'g_intro'].forEach(function (k) { if ($(k)) $(k).hidden = on; });
    if (on) { draw(); genNote(); }
  }

  // ------------------------------------------------------------ つなぐ
  function wire() {
    $('dem_file').onchange = function () { loadFiles(this.files); };
    $('dem_sample').onclick = function () { setDEM(SD.synthDEM({}), '例の DEM（架空）'); };
    $('dem_thr').onchange = analyze;
    $('dem_smooth').onchange = applySmooth;
    $('dm_fit').onchange = function () { if (!ST.D) return; if ($('dm_fit').checked) fitBasin(); else fitAll(); draw(); };
    $('dm_shade').onchange = $('dm_cont').onchange = function () { ST.img = null; ST.det = null; draw(); };
    $('dm_zin').onclick = function () { var v = ST.view; if (v) zoomAt(v.W / 2, v.H / 2, 1.6); };
    $('dm_zout').onclick = function () { var v = ST.view; if (v) zoomAt(v.W / 2, v.H / 2, 1 / 1.6); };
    $('dm_zall').onclick = function () { if (ST.D) { fitAll(); draw(); } };
    $('dm_zbas').onclick = function () { if (ST.D) { fitBasin(); draw(); } };
    // ホイールで拡大縮小（カーソルの位置を中心に）、ドラッグで移動。ドラッグのあとのクリックは無視
    var mv = $('demmap'), drag = null;
    mv.addEventListener('wheel', function (e) {
      if (!ST.D) return; e.preventDefault();
      var r = mv.getBoundingClientRect(), d = e.deltaMode === 1 ? e.deltaY * 33 : e.deltaY;
      zoomAt(e.clientX - r.left, e.clientY - r.top, Math.exp(-d * 0.0015));
    }, { passive: false });
    mv.addEventListener('pointerdown', function (e) { if (!ST.cam || e.button !== 0) return; drag = { x: e.clientX, y: e.clientY, cx: ST.cam.cx, cy: ST.cam.cy, moved: false, id: e.pointerId }; ST.dragged = false; });
    window.addEventListener('pointermove', function (e) {
      if (!drag || e.pointerId !== drag.id) return;
      var dx = e.clientX - drag.x, dy = e.clientY - drag.y;
      if (!drag.moved && dx * dx + dy * dy < 25) return;
      if (!drag.moved) { drag.moved = true; try { mv.setPointerCapture(drag.id); } catch (er) {} mv.style.cursor = 'grabbing'; $('dem_tip').hidden = true; }
      ST.cam.cx = drag.cx - dx / ST.cam.sc; ST.cam.cy = drag.cy - dy / ST.cam.sc; redraw();
    });
    window.addEventListener('pointerup', function (e) { if (drag && e.pointerId === drag.id) { ST.dragged = drag.moved; drag = null; mv.style.cursor = 'crosshair'; } });
    $('demmap').addEventListener('click', onClick);
    $('demmap').addEventListener('dblclick', function (e) { if (($('dm_mode').value === 'slide' && $('dm_shape').value === 'poly') || $('dm_mode').value === 'mat') { e.preventDefault(); closePoly(); } });
    $('dm_mode').addEventListener('change', shapeUI); $('dm_shape').addEventListener('change', shapeUI);
    $('dm_close').onclick = closePoly;
    $('dm_undo').onclick = function () { if (ST.draft && ST.draft.length) { ST.draft.pop(); polyMsg(); draw(); } };
    $('demmap').addEventListener('mousemove', onMove);
    $('demmap').addEventListener('mouseleave', function () { $('dem_tip').hidden = true; });
    $('src_all').onclick = function () { if (!ST.src.length) return; var f = ST.src[0]; ST.src.forEach(function (r) { ['mode', 'Qb', 'Qp', 't0', 'tr', 'tf', 'C'].forEach(function (k) { r[k] = f[k]; }); }); renderTables(); };
    $('src_top').onclick = function () { topHeads(3); renderTables(); draw(); selInfo(); };
    $('src_clr').onclick = function () { ST.src = []; renderTables(); draw(); selInfo(); };
    $('sl_clr').onclick = function () { ST.slides = []; renderTables(); draw(); selInfo(); };
    $('ld_clr').onclick = function () { snap(); ST.ldams = []; renderTables(); draw(); selInfo(); };
    $('dmt_clr').onclick = function () { snap(); ST.dams = []; renderTables(); draw(); selInfo(); };
    $('dm_back').onclick = undo;
    document.addEventListener('keydown', function (e) {
      if ((e.ctrlKey || e.metaKey) && !e.shiftKey && (e.key === 'z' || e.key === 'Z') && !$('demcard').hidden && !/INPUT|TEXTAREA|SELECT/.test((document.activeElement || {}).tagName || '')) { e.preventDefault(); undo(); }
      if (e.key === 'Escape' && ST.draft && ST.draft.length) { ST.draft = []; polyMsg(); draw(); }
    });
    $('demmap').addEventListener('contextmenu', function (e) {
      if (!ST.D) return; e.preventDefault(); var c = cellAt(e); if (c < 0) return;
      var nm = removeAt(c); $('dem_sel').textContent = nm ? nm + ' を外しました（「元に戻す」か Ctrl+Z で戻せます）。' : '近くに外せる物がありません。'; renderTables(); draw();
    });
    $('ld_ex').onclick = ldamExample;
    ['rn_on', 'rn_loss', 'rn_route', 'rn_meth', 'rn_mode'].forEach(function (id) { $(id).addEventListener('change', rainSwitch); });
    $('rn_calc').onclick = calcRain;
    $('rn_file').onchange = function () {
      var f = this.files[0]; if (!f) return; var r = new FileReader();
      r.onload = function () { try { if (!ST.D) throw new Error('先に DEM を読み込んでください'); var o = parseRainCSV(r.result); ST.rainFrames = o.frames; $('rn_ginfo').textContent = o.info; } catch (e) { $('rn_ginfo').textContent = e.message; } };
      r.readAsText(f);
    };
    $('gen_type').onchange = genNote;
    $('g_src').addEventListener('change', function () { show(this.value === 'dem'); });
    window.addEventListener('resize', function () { if (!$('demcard').hidden && ST.D) { if ($('dm_fit').checked) fitBasin(); else fitAll(); draw(); } });
    rainSwitch(); show($('g_src').value === 'dem');
  }
  window.SaboDEMUI = { params: params, state: ST, draw: draw, setDEM: setDEM, paths: function () { return ST.paths; }, active: function () { return $('g_src').value === 'dem'; } };
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', wire); else wire();
})();
