/* DEM と流域の画面（地図・流域の選択・供給点・雨・斜面崩壊・天然ダム）
 * ui.js の readForm から SaboDEMUI.params(p) を呼んで、計算の条件 p を作る */
(function () {
  'use strict';
  var SD = window.SaboDEM, SH = window.SaboHydro;
  var $ = function (id) { return document.getElementById(id); };
  var ST = { D: null, A: null, mask: null, outlet: -1, bid: -1, heads: [], src: [], slides: [], ldams: [], img: null, view: null, rainFrames: null, lat: null, dom: null };

  function fmt(v, n) { return isFinite(v) ? (+v).toFixed(n) : ''; }
  function info(msg, err) { var e = $('dem_info'); e.textContent = msg; e.style.color = err ? 'var(--err,#b3261e)' : ''; }

  // ------------------------------------------------------------ 読み込み
  function inflateBrowser(u8) {
    var ds = new DecompressionStream('deflate');
    var w = ds.writable.getWriter(); w.write(u8); w.close();
    return new Response(ds.readable).arrayBuffer().then(function (b) { return new Uint8Array(b); });
  }
  function loadFiles(files) {
    if (!files || !files.length) return;
    info('読み込み中…');
    var arr = Array.prototype.slice.call(files), first = arr[0], name = first.name.toLowerCase();
    var asText = function (f) { return new Promise(function (ok, ng) { var r = new FileReader(); r.onload = function () { ok(r.result); }; r.onerror = ng; r.readAsText(f); }); };
    var p;
    if (/\.tiff?$/.test(name)) p = first.arrayBuffer().then(function (b) { return SD.parseTIFF(b, inflateBrowser); });
    else if (/\.xml$/.test(name)) p = Promise.all(arr.map(asText)).then(function (ts) { return SD.parseGSIXML(ts); });
    else p = asText(first).then(function (t) { return /^\s*ncols/i.test(t) ? SD.parseASC(t) : SD.parseXYZ(t); });
    p.then(function (D) { setDEM(D, first.name + (arr.length > 1 ? ' ほか ' + (arr.length - 1) + ' 個' : '')); })
      .catch(function (e) { info('読み込めませんでした：' + e.message, true); });
  }
  function setDEM(D, label) {
    var N = D.nx * D.ny;
    if (N > 9e6) { info('DEM が大きすぎます（' + D.nx + '×' + D.ny + '）。900 万セル以下に切り出してください。', true); return; }
    ST.D = D; ST.label = label; analyze();
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
    info(ST.label + '：' + D.nx + ' × ' + D.ny + '（格子 ' + fmt(D.cs, 2) + ' m、標高 ' + fmt(zmin, 0) + '〜' + fmt(zmax, 0) + ' m）。谷頭 ' + A.heads.length + ' 点、出口ごとの流域 ' + A.basins.length + '（' + Math.round(performance.now() - t0) + ' ms）');
    ST.img = null; draw();
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
    if (!ST.src.length && ST.heads.length) topHeads(3);
    ST.lat = null; ST.img = null;
    renderTables(); draw();
    if (!quiet) selInfo();
  }
  function selInfo() {
    if (!ST.mask) { $('dem_sel').textContent = ''; return; }
    var n = 0; for (var c = 0; c < ST.mask.length; c++) n += ST.mask[c];
    var a = n * ST.D.cs * ST.D.cs / 1e6;
    $('dem_sel').textContent = '選んだ流域：面積 ' + fmt(a, a < 1 ? 3 : 2) + ' km²、谷頭 ' + ST.heads.length + ' 点、供給点 ' + ST.src.length + ' 点、斜面崩壊 ' + ST.slides.length + ' か所、天然ダム ' + ST.ldams.length + ' か所。「適用して最初から」で計算領域に切り出します。';
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
      tr.querySelector('.x').onclick = function () { ST.src.splice(i, 1); renderTables(); draw(); selInfo(); };
      tb.appendChild(tr);
    });
    tb = $('sltab').tBodies[0]; tb.innerHTML = '';
    ST.slides.forEach(function (r, i) {
      var tr = document.createElement('tr');
      tr.innerHTML = '<td>' + (i + 1) + '</td><td><select data-k="type">' + opt(r.type, [['shallow', '表層'], ['deep', '深層']]) + '</select></td>' +
        '<td><input type="number" step="any" value="' + r.R + '" data-k="R"></td><td><input type="number" step="any" value="' + r.D + '" data-k="D"></td>' +
        '<td><select data-k="trigger">' + opt(r.trigger, [['fs', '安全率'], ['time', '時刻']]) + '</select></td><td><input type="number" step="any" value="' + r.t + '" data-k="t"></td>' +
        '<td><select data-k="move">' + opt(r.move, [['auto', '自動'], ['fluid', '土石流化'], ['block', '土塊']]) + '</select></td>' +
        '<td><input type="number" step="any" value="' + r.phi + '" data-k="phi"></td><td><input type="number" step="any" value="' + r.coh + '" data-k="coh"></td><td><input type="number" step="any" value="' + r.m0 + '" data-k="m0"></td>' +
        '<td><button class="x" aria-label="消す">×</button></td>';
      Array.prototype.forEach.call(tr.querySelectorAll('input,select'), function (inp) {
        inp.onchange = function () {
          var k = inp.dataset.k; r[k] = inp.tagName === 'SELECT' ? inp.value : (parseFloat(inp.value) || 0);
          if (k === 'type') { var d = slideDefaults(r.type); r.R = d.R; r.D = d.D; r.phi = d.phi; r.coh = d.coh; renderTables(); }
          draw();
        };
      });
      tr.querySelector('.x').onclick = function () { ST.slides.splice(i, 1); renderTables(); draw(); selInfo(); };
      tb.appendChild(tr);
    });
    tb = $('ldtab').tBodies[0]; tb.innerHTML = '';
    ST.ldams.forEach(function (r, i) {
      var tr = document.createElement('tr');
      tr.innerHTML = '<td>' + (i + 1) + '</td>' + ['H', 'Wc', 'up', 'dn'].map(function (k) { return '<td><input type="number" step="any" value="' + r[k] + '" data-k="' + k + '"></td>'; }).join('') + '<td><button class="x" aria-label="消す">×</button></td>';
      Array.prototype.forEach.call(tr.querySelectorAll('input'), function (inp) { inp.onchange = function () { r[inp.dataset.k] = parseFloat(inp.value) || 0; draw(); }; });
      tr.querySelector('.x').onclick = function () { ST.ldams.splice(i, 1); renderTables(); draw(); selInfo(); };
      tb.appendChild(tr);
    });
  }
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
    ST.ldams = [{ c: best, H: Math.max(5, Math.min(12, Math.round(bd * 0.3))), Wc: 10, up: 30, dn: 25 }];
    ST.src = [{ c: up, mode: 'const', Qb: 40, Qp: 40, t0: 0, tr: 60, tf: 300, C: 0 }];
    $('gen_type').value = 'ldam'; genNote();
    renderTables(); draw(); selInfo();
    $('dem_sel').textContent += '　例：天然ダム（高さ ' + ST.ldams[0].H + ' m）と、その上流の本川に一定流量 40 m³/s の供給点を置きました（谷頭でなくても供給点にできます）。';
  }
  function slideDefaults(type) { return type === 'deep' ? { R: 80, D: 20, phi: 30, coh: 30 } : { R: 25, D: 1.5, phi: 35, coh: 3 }; }

  // ------------------------------------------------------------ 地図
  function shade(D, A) {
    var nx = D.nx, ny = D.ny, f = Math.max(1, Math.ceil(Math.max(nx, ny) / 1100)), W = Math.ceil(nx / f), H = Math.ceil(ny / f);
    var img = new ImageData(W, H), px = img.data, z = D.z, m = ST.mask, cs = D.cs * f;
    var zmin = Infinity, zmax = -Infinity, c;
    for (c = 0; c < nx * ny; c++) { var v = z[c]; if (isFinite(v)) { if (v < zmin) zmin = v; if (v > zmax) zmax = v; } }
    for (var j = 0; j < H; j++) for (var i = 0; i < W; i++) {
      var x = Math.min(nx - 1, i * f), y = Math.min(ny - 1, j * f); c = y * nx + x;
      var p = (j * W + i) * 4, zz = z[c];
      if (!isFinite(zz)) { px[p] = px[p + 1] = px[p + 2] = 235; px[p + 3] = 255; continue; }
      var xe = Math.min(nx - 1, x + f), yn = Math.max(0, y - f);
      var dzdx = (z[y * nx + xe] - zz) / cs, dzdy = (z[yn * nx + x] - zz) / cs;
      if (!isFinite(dzdx)) dzdx = 0; if (!isFinite(dzdy)) dzdy = 0;
      var hs = Math.max(0, Math.min(1, 0.55 + 0.45 * (-dzdx * 0.6 + dzdy * 0.6) / Math.sqrt(1 + dzdx * dzdx + dzdy * dzdy)));
      var t = (zz - zmin) / Math.max(1e-6, zmax - zmin);
      var r = 120 + 110 * t, g = 150 + 70 * t, b = 110 + 90 * t;
      var inB = m ? m[c] : 1;
      if (!inB) { var gr = (r + g + b) / 3; r = g = b = 0.6 * gr + 70; }
      if (A && A.stream[c]) { var o = A.strahler[c]; r = inB ? 30 : 140; g = inB ? 110 + 10 * o : 160; b = inB ? 200 : 190; hs = 1; }
      px[p] = r * hs; px[p + 1] = g * hs; px[p + 2] = b * hs; px[p + 3] = 255;
    }
    return { img: img, f: f, W: W, H: H };
  }
  function viewBox() {
    var D = ST.D, x0 = 0, y0 = 0, x1 = D.nx, y1 = D.ny;
    if ($('dm_fit').checked && ST.mask) {
      x0 = D.nx; y0 = D.ny; x1 = 0; y1 = 0;
      for (var c = 0; c < ST.mask.length; c++) if (ST.mask[c]) { var x = c % D.nx, y = (c / D.nx) | 0; if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y; }
      var mg = Math.max(4, 0.06 * Math.max(x1 - x0, y1 - y0)); x0 = Math.max(0, x0 - mg); y0 = Math.max(0, y0 - mg); x1 = Math.min(D.nx, x1 + mg + 1); y1 = Math.min(D.ny, y1 + mg + 1);
    }
    return { x0: x0, y0: y0, x1: x1, y1: y1 };
  }
  function draw() {
    var cv = $('demmap'); if (!cv || $('demcard').hidden) return;
    var D = ST.D, dpr = window.devicePixelRatio || 1, Wc = cv.parentNode.clientWidth || 600;
    if (!D) { cv.width = Wc * dpr; cv.height = 120 * dpr; cv.style.height = '120px'; var c0 = cv.getContext('2d'); c0.setTransform(dpr, 0, 0, dpr, 0, 0); c0.clearRect(0, 0, Wc, 120); c0.fillStyle = '#888'; c0.font = '13px sans-serif'; c0.fillText('DEM を読み込むか「例の DEM」を押してください', 12, 60); return; }
    if (!ST.img) ST.img = shade(D, ST.A);
    var vb = viewBox(), vw = vb.x1 - vb.x0, vh = vb.y1 - vb.y0, Hc = Math.min(560, Math.max(220, Wc * vh / vw));
    var sc = Math.min(Wc / vw, Hc / vh); Hc = vh * sc;
    cv.width = Wc * dpr; cv.height = Hc * dpr; cv.style.height = Hc + 'px';
    var ctx = cv.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, Wc, Hc);
    var I = ST.img, tmp = ST.tmpCv || (ST.tmpCv = document.createElement('canvas'));
    if (tmp.width !== I.W || tmp.height !== I.H || ST.tmpImg !== I) { tmp.width = I.W; tmp.height = I.H; tmp.getContext('2d').putImageData(I.img, 0, 0); ST.tmpImg = I; }
    ctx.imageSmoothingEnabled = I.f > 1;
    ctx.drawImage(tmp, vb.x0 / I.f, vb.y0 / I.f, vw / I.f, vh / I.f, 0, 0, vw * sc, vh * sc);
    ST.view = { x0: vb.x0, y0: vb.y0, sc: sc, W: Wc, H: Hc };
    function P(c) { return [((c % D.nx) + 0.5 - vb.x0) * sc, (((c / D.nx) | 0) + 0.5 - vb.y0) * sc]; }
    // 流域界（選んだ流域の縁）
    if (ST.mask) {
      ctx.fillStyle = 'rgba(200,40,40,0.9)';
      var m = ST.mask, step = Math.max(1, Math.floor(1 / sc));
      for (var y = Math.floor(vb.y0); y < vb.y1; y += step) for (var x = Math.floor(vb.x0); x < vb.x1; x += step) {
        var c = y * D.nx + x; if (!m[c]) continue;
        if (x === 0 || y === 0 || x === D.nx - 1 || y === D.ny - 1 || !m[c - 1] || !m[c + 1] || !m[c - D.nx] || !m[c + D.nx]) ctx.fillRect((x - vb.x0) * sc, (y - vb.y0) * sc, Math.max(1.2, sc * step), Math.max(1.2, sc * step));
      }
      var po = P(ST.outlet); ctx.strokeStyle = '#c62828'; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(po[0], po[1], 7, 0, 7); ctx.stroke();
    }
    // 谷頭
    var selSet = {}; ST.src.forEach(function (r, i) { selSet[r.c] = i + 1; });
    ST.heads.forEach(function (c) {
      var p = P(c); if (p[0] < -5 || p[1] < -5 || p[0] > Wc + 5 || p[1] > Hc + 5) return;
      var s = selSet[c];
      ctx.beginPath(); ctx.arc(p[0], p[1], s ? 6 : 3.2, 0, 7);
      ctx.fillStyle = s ? '#ef6c00' : '#fff'; ctx.fill(); ctx.lineWidth = 1.2; ctx.strokeStyle = '#333'; ctx.stroke();
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
    ST.slides.forEach(function (r, i) {
      var p = P(r.c), R = r.R / D.cs * sc;
      ctx.beginPath(); ctx.arc(p[0], p[1], Math.max(4, R), 0, 7); ctx.fillStyle = r.type === 'deep' ? 'rgba(121,85,72,0.45)' : 'rgba(255,152,0,0.4)'; ctx.fill();
      ctx.strokeStyle = '#5d4037'; ctx.lineWidth = 1.5; ctx.stroke(); ctx.fillStyle = '#3e2723'; ctx.font = '600 11px sans-serif'; ctx.fillText('崩' + (i + 1), p[0] + Math.max(4, R) + 2, p[1]);
    });
    // 天然ダム
    ST.ldams.forEach(function (r, i) {
      var p = P(r.c), q = ST.A.dir[r.c], dxv = q >= 0 ? SD.DC[q] : 1, dyv = q >= 0 ? SD.DR[q] : 0, L = Math.sqrt(dxv * dxv + dyv * dyv), wv = Math.max(8, r.H * 3 / D.cs * sc);
      var nxv = -dyv / L, nyv = dxv / L;
      ctx.strokeStyle = '#4e342e'; ctx.lineWidth = Math.max(4, r.Wc / D.cs * sc); ctx.beginPath(); ctx.moveTo(p[0] - nxv * wv, p[1] - nyv * wv); ctx.lineTo(p[0] + nxv * wv, p[1] + nyv * wv); ctx.stroke();
      ctx.fillStyle = '#3e2723'; ctx.font = '600 11px sans-serif'; ctx.fillText('ダム' + (i + 1), p[0] + wv + 3, p[1] + 4);
    });
    $('dem_key').innerHTML = '<span><i style="border-color:#1e6ec8"></i>谷筋</span><span><i class="b" style="background:#fff;border:1px solid #333"></i>谷頭</span><span><i class="b" style="background:#ef6c00"></i>供給点</span><span><i style="border-color:#c62828"></i>流域界・出口</span><span><i class="b" style="background:rgba(255,152,0,0.5)"></i>斜面崩壊</span><span><i class="b" style="background:#4e342e"></i>天然ダム</span>';
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
  function nearStream(c, r) {
    var D = ST.D, A = ST.A, x = c % D.nx, y = (c / D.nx) | 0, best = -1, bd = Infinity;
    for (var dy = -r; dy <= r; dy++) for (var dx = -r; dx <= r; dx++) {
      var xx = x + dx, yy = y + dy; if (xx < 0 || yy < 0 || xx >= D.nx || yy >= D.ny) continue;
      var o = yy * D.nx + xx; if (!A.stream[o]) continue;
      var d = dx * dx + dy * dy; if (d < bd || (d === bd && A.acc[o] > A.acc[best])) { bd = d; best = o; }
    }
    return best;
  }
  function onClick(e) {
    var c = cellAt(e); if (c < 0) return;
    var mode = $('dm_mode').value, A = ST.A;
    if (mode === 'basin') { if (A.bid[c] >= 0) selectBasin(A.bid[c], false); return; }
    if (mode === 'outlet') { var s = nearStream(c, Math.max(3, Math.round(14 / ST.view.sc))); if (s >= 0) selectOutlet(s); else $('dem_sel').textContent = '谷筋（青い線）の近くをクリックしてください。'; return; }
    if (!ST.mask || !ST.mask[c]) { $('dem_sel').textContent = '選んだ流域の中をクリックしてください。'; return; }
    if (mode === 'src') {
      var h = nearest(c, ST.heads, 14); if (h < 0) return;
      var k = -1; ST.src.forEach(function (r, i) { if (r.c === h) k = i; });
      if (k >= 0) ST.src.splice(k, 1); else ST.src.push(srcDefault(h));
    } else if (mode === 'slide') {
      var d = slideDefaults($('gen_type').value === 'ldam' ? 'deep' : 'shallow');
      ST.slides.push({ c: c, type: $('gen_type').value === 'ldam' ? 'deep' : 'shallow', R: d.R, D: d.D, trigger: $('rn_on').checked ? 'fs' : 'time', t: 30, move: $('gen_type').value === 'ldam' ? 'block' : 'auto', phi: d.phi, coh: d.coh, m0: $('rn_on').checked ? 0.2 : 0.9 });
    } else if (mode === 'ldam') {
      var s2 = nearStream(c, Math.max(3, Math.round(14 / ST.view.sc))); if (s2 < 0) { $('dem_sel').textContent = '谷筋（青い線）の近くをクリックしてください。'; return; }
      ST.ldams.push({ c: s2, H: 15, Wc: 10, up: 30, dn: 25 });
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
  function params(p) {
    if (!ST.D || !ST.mask) throw new Error('DEM を読み込み、流域を選んでください');
    var D = ST.D, A = ST.A;
    if (!(p.dx >= D.cs)) p.dx = D.cs;
    var dom = SD.prepareDomain(D, A, ST.mask, ST.outlet, { fill: $('dem_fill').checked, dx: p.dx });
    ST.dom = dom;
    p.terrain = { rows: dom.rows, cs: dom.cs, basin: dom.basin, chan: dom.chan, Iout: dom.Iout };
    p.csvSlope = 0; p.L = dom.nx * dom.cs; p.sideCells = 0; p.dams = []; p.tribs = []; p.prefill = 0;
    var dur = Math.max(1, +$('dem_dur').value || 30) * 60;
    p.dem = { soilD: Math.max(0.1, +$('dem_soil').value || 1), chW: Math.max(1, +$('dem_chw').value || 10), ks: Math.max(0, +$('dem_ks').value || 0), s0: 0.3, duration: dur };
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
      return { x: xy[0], y: xy[1], R: r.R, D: r.D, type: r.type, trigger: r.trigger, t: r.t, move: r.move, phi: r.phi, c: r.coh, m0: r.m0, thDep: thD, mu: mu, CL: CL, Tk: r.type === 'deep' ? +$('sl_tk2').value : +$('sl_tk1').value };
    });
    p.ldams = ST.ldams.map(function (r) { var xy = dom.toXY(r.c), u = dom.dirXY(r.c); return { x: xy[0], y: xy[1], ux: u[0], uy: u[1], H: r.H, Wc: r.Wc, up: r.up, dn: r.dn }; });
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
    $('g_dem').hidden = !on; $('demcard').hidden = !on; $('dem_more').hidden = !on;
    if (on) { draw(); genNote(); }
  }

  // ------------------------------------------------------------ つなぐ
  function wire() {
    $('dem_file').onchange = function () { loadFiles(this.files); };
    $('dem_sample').onclick = function () { setDEM(SD.synthDEM({}), '例の DEM（架空）'); };
    $('dem_thr').onchange = analyze;
    $('dm_fit').onchange = draw;
    $('demmap').addEventListener('click', onClick);
    $('demmap').addEventListener('mousemove', onMove);
    $('demmap').addEventListener('mouseleave', function () { $('dem_tip').hidden = true; });
    $('src_all').onclick = function () { if (!ST.src.length) return; var f = ST.src[0]; ST.src.forEach(function (r) { ['mode', 'Qb', 'Qp', 't0', 'tr', 'tf', 'C'].forEach(function (k) { r[k] = f[k]; }); }); renderTables(); };
    $('src_top').onclick = function () { topHeads(3); renderTables(); draw(); selInfo(); };
    $('src_clr').onclick = function () { ST.src = []; renderTables(); draw(); selInfo(); };
    $('sl_clr').onclick = function () { ST.slides = []; renderTables(); draw(); selInfo(); };
    $('ld_clr').onclick = function () { ST.ldams = []; renderTables(); draw(); selInfo(); };
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
    window.addEventListener('resize', function () { if (!$('demcard').hidden) draw(); });
    rainSwitch(); show($('g_src').value === 'dem');
  }
  window.SaboDEMUI = { params: params, state: ST, draw: draw, setDEM: setDEM, active: function () { return $('g_src').value === 'dem'; } };
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', wire); else wire();
})();
