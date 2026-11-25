(function () {
  'use strict';
  var S2 = window.Sabo2D, G = 9.81;
  /*TESTCSV*/
  var $ = function (id) { return document.getElementById(id); };
  var sim = null, running = false, ref = null, P = S2.defaultParams();
  P.dams = [60, 105, 150, 195, 240, 285, 330].map(function (x) { return { x: x, H: 1.2, b: 0 }; });

  // ------------------------------------------------------------ 入力欄
  // 堰堤の行：位置・堤高・水通し幅と、任意で回転角・中心 y・天端標高・堤長（空欄は自動）
  function fv(v, n) { return isFinite(v) && v !== null && v !== '' ? +(+v).toFixed(n) : ''; }
  var DTYPES = [['closed', '不透過'], ['slit', 'スリット'], ['mslit', '複数スリット'], ['grid', '鋼製格子型'], ['sslit', '鋼製スリット型']];
  var DCODE = { closed: 0, slit: 1, mslit: 2, grid: 3, sslit: 4 }, DNAME = ['closed', 'slit', 'mslit', 'grid', 'sslit'];
  function damRow(d) {
    var tr = document.createElement('tr');
    tr.innerHTML = '<td><input type="number" step="5" min="0" value="' + d.x + '" aria-label="位置"></td>' +
      '<td><input type="number" step="0.1" min="0.1" value="' + d.H + '" aria-label="堤高"></td>' +
      '<td><input type="number" step="1" min="0" value="' + d.b + '" aria-label="水通し幅"></td>' +
      '<td><select aria-label="形式">' + DTYPES.map(function (t) { return '<option value="' + t[0] + '"' + ((d.type || 'closed') === t[0] ? ' selected' : '') + '>' + t[1] + '</option>'; }).join('') + '</select></td>' +
      '<td><input type="number" step="0.1" min="0.05" value="' + (d.l > 0 ? d.l : '') + '" placeholder="–" aria-label="開口幅"></td>' +
      '<td><input type="number" step="1" min="1" value="' + (d.n > 0 ? d.n : '') + '" placeholder="–" aria-label="開口数"></td>' +
      '<td class="adv"><input type="number" step="1" value="' + fv(d.ang, 4) + '" placeholder="自動" aria-label="回転角"></td>' +
      '<td class="adv"><input type="number" step="1" value="' + fv(d.y, 4) + '" placeholder="自動" aria-label="中心 y"></td>' +
      '<td class="adv"><input type="number" step="0.1" value="' + fv(d.zc, 3) + '" placeholder="自動" aria-label="天端標高"></td>' +
      '<td class="adv"><input type="number" step="2" min="2" value="' + fv(d.len, 2) + '" placeholder="自動" aria-label="堤長"></td>' +
      '<td><button class="x" aria-label="削除">×</button></td>';
    if (isFinite(d.so) && d.so !== null && d.b > 0) { tr.dataset.so = d.so; tr.dataset.b0 = d.b; }   // 水通しの位置（固定したときだけ）
    if (d.hb > 0) tr.dataset.hb = d.hb;
    tr.querySelector('button').onclick = function () { tr.remove(); };
    return tr;
  }
  function advSwitch() { $('damtab').classList.toggle('simple', !$('d_adv').checked); $('dam_adv_note').hidden = !$('d_adv').checked; }
  // いまの計算の堰堤を、位置・向き・天端標高・堤長まで表に書き込む（地形が変わっても同じ場所に置くため）
  function damsFixed() {
    return sim.dams.map(function (D) {
      var sp = D.spMax - D.spMin;
      var o = { x: +D.x.toFixed(3), H: D.H, b: +Math.max(0, sp).toFixed(3), ang: +(D.th * 180 / Math.PI).toFixed(4), y: +D.cy.toFixed(4),
        zc: +D.crest.toFixed(3), len: +(2 * Math.max(Math.abs(D.sMin), Math.abs(D.sMax))).toFixed(2), so: +((D.spMin + D.spMax) / 2).toFixed(4) };
      if (D.op) { o.type = D.op.type; o.l = D.op.l; o.n = D.op.n; if (D.op.hb > 0) o.hb = +D.op.hb.toFixed(3); }
      return o;
    });
  }
  function fixDams() { if (!sim || !sim.dams.length) return 0; setDams(damsFixed()); $('d_adv').checked = true; advSwitch(); presetDams = false; return sim.dams.length; }
  // CSV のコメント行に書いた堰堤（# dam x=… H=… …）を読む
  function damsFromComments(text) {
    var out = [], fh = null, mask = [];
    String(text).split(/\r?\n/).forEach(function (ln) {
      var m = ln.match(/^\s*#\s*dam\s+(.*)$/i);
      if (m) { var d = {}; m[1].replace(/(\w+)=(-?[\d.eE+-]+)/g, function (_, k, v) { d[k] = parseFloat(v); }); if (isFinite(d.ty) && DNAME[d.ty]) { d.type = DNAME[d.ty]; delete d.ty; } if (isFinite(d.x) && d.H > 0) out.push(d); }
      var f = ln.match(/^\s*#\s*floorH=(-?[\d.]+)/i); if (f) fh = parseFloat(f[1]);
      var fm = ln.match(/^\s*#\s*floor\s+([01]+)\s*$/); if (fm) mask.push(fm[1].split('').map(Number));
    });
    return { dams: out, floorH: fh, mask: mask.length ? mask : null };
  }
  // 計算後の河床を格子 CSV にする（堰堤のセルは周りの地盤の高さ。堰堤の諸元はコメント行に書く）
  function bedCSV(zA, tt) {
    zA = zA || sim.z; if (tt == null) tt = sim.t;
    var nx = sim.nx, ny = sim.ny, lines = [], fh = 0;
    for (var i = 0; i < nx; i++) {
      var mn = 1e30, mx = -1e30;
      for (var j = 0; j < ny; j++) { var c = j * nx + i; if (sim.kind[c] === 0) { if (zA[c] < mn) mn = zA[c]; if (zA[c] > mx) mx = zA[c]; } }
      if (mx > mn && mx - mn > fh) fh = mx - mn;
    }
    fh = Math.min(20, Math.max(0.5, Math.ceil((fh + 0.05) * 10) / 10));
    lines.push('# 計算後の河床（出水の時刻 ' + (tt / 60).toFixed(1) + ' 分）　1行目が左岸側、左の列が上流　格子間隔 ' + sim.dx + ' m');
    lines.push('# floorH=' + fh);
    damsFixed().forEach(function (d) { lines.push('# dam x=' + d.x + ' H=' + d.H + ' b=' + d.b + ' ang=' + d.ang + ' y=' + d.y + ' zc=' + d.zc + ' len=' + d.len + ' so=' + d.so + (d.type && d.type !== 'closed' ? ' ty=' + DCODE[d.type] + ' l=' + d.l + ' n=' + d.n + (d.hb ? ' hb=' + d.hb : '') : '')); });
    // 谷底（侵食できるセル）の印。1 が谷底、0 が谷壁・堰堤。行の並びは数値の行と同じ
    for (var jm = ny - 1; jm >= 0; jm--) { var sm = ''; for (var im = 0; im < nx; im++) { var cm = jm * nx + im; sm += sim.kind[cm] === 0 || (sim.kind[cm] === 2 && sim.fl[cm]) ? '1' : '0'; } lines.push('# floor ' + sm); }
    for (var j2 = ny - 1; j2 >= 0; j2--) { var r = []; for (var i2 = 0; i2 < nx; i2++) { var c2 = j2 * nx + i2, v = sim.kind[c2] === 2 ? Math.min(damGround(c2, zA), zA[c2] - 0.01) : zA[c2]; r.push(v.toFixed(3)); } lines.push(r.join(',')); }   // 堰堤の下は天端より少し低くして、読み込み直しても同じセルに堰堤が立つように
    return { text: lines.join('\r\n'), floorH: fh };
  }
  // ファイルとして保存（公開ページでは保存の確認が出る。オフライン版はブラウザのダウンロード）
  function saveText(name, text, msgId, raw) {
    var mid = msgId || 'csv_msg', data = raw ? text : '\uFEFF' + text, done = function () { $(mid).hidden = false; $(mid).textContent = name + ' を保存しました。'; };
    var mime = /\.zip$/.test(name) ? 'application/zip' : /\.json$/.test(name) ? 'application/json' : 'text/csv';
    var fallback = function () { try { var a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([data], { type: mime })); a.download = name; document.body.appendChild(a); a.click(); setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 1000); done(); } catch (e) { $(mid).hidden = false; $(mid).textContent = '保存できませんでした。' + (mid === 'csv_msg' ? '「CSV をコピー」を使ってください。' : ''); } };
    if (window.claude && window.claude.use) {
      window.claude.use('downloads').then(function (dl) {
        if (!dl) { fallback(); return; }
        dl.save({ filename: name, data: data }).then(done, function (e) { if (e && e.code === 'declined') return; fallback(); });
      }, fallback);
    } else fallback();
  }
  // ---- 時間ごとの結果の書き出し -----------------------------------------------------
  function snapLabel(t) {   // ファイル名の時刻（分）。間隔が分の整数倍なら整数、そうでなければ小数 1 桁
    var ev = sim.P.snapEvery || 600, w = String(Math.round(sim.duration() / 60)).length;
    if (Math.abs(ev / 60 - Math.round(ev / 60)) < 1e-9) { var m = String(Math.round(t / 60)); while (m.length < w) m = '0' + m; return m; }
    return (t / 60).toFixed(1);
  }
  function snapZ(sn) { var z0 = sim.z0, N = sim.N, z = new Float64Array(N); for (var c = 0; c < N; c++) z[c] = z0[c] + sn.dz[c]; return z; }
  function snapInfo() {
    if (!sim) return;
    var S = sim.snaps, ev = sim.P.snapEvery || 0, el = $('o_info');
    if (!(ev > 0)) { el.textContent = '記録していません（記録の間隔が 0）。'; }
    else if (!S.length) { el.textContent = '計算を始めると、出水の時刻 ' + +(ev / 60).toFixed(2) + ' 分ごとに記録します。'; }
    else el.textContent = '記録 ' + S.length + ' 回（' + (S[0].t / 60).toFixed(1) + '〜' + (S[S.length - 1].t / 60).toFixed(1) + ' 分）。' + (sim.snapFull ? '　記録の量が上限（約 100 MB）に達したため、これより後は記録していません。間隔を広げてください。' : '');
    $('o_long').disabled = $('o_zip').disabled = !S.length;
  }
  function longCSV() {
    var nx = sim.nx, ny = sim.ny, dx = sim.dx, z0 = sim.z0, S = sim.snaps, out = [];
    out.push('# 時間ごとの結果（記録の間隔 ' + +((sim.P.snapEvery || 0) / 60).toFixed(2) + ' 分、格子間隔 ' + dx + ' m）。移動床のセルだけ');
    out.push('# x：上流端からの流下方向の距離。y：右岸側の端からの距離（左岸側が正）。河床変動量は計算開始時の河床からの変化（＋が堆積）。水深 0 は乾いたセル');
    out.push('時刻_分,x_m,y_m,河床高_m,河床変動量_m,水深_m,水位_m,流速_m/s,表層平均粒径_mm');
    for (var q = 0; q < S.length; q++) {
      var sn = S[q], tm = (sn.t / 60).toFixed(2);
      for (var i = 0; i < nx; i++) for (var j = 0; j < ny; j++) {
        var c = j * nx + i; if (sim.kind[c] !== 0) continue;
        var zz = z0[c] + sn.dz[c], hh = sn.h[c] > 1e-3 ? sn.h[c] : 0, dm = sn.dm[c];
        out.push(tm + ',' + ((i + 0.5) * dx).toFixed(2) + ',' + ((j + 0.5) * dx).toFixed(2) + ',' + zz.toFixed(3) + ',' + sn.dz[c].toFixed(4) + ',' +
          hh.toFixed(3) + ',' + (hh > 0 ? (zz + hh).toFixed(3) : '') + ',' + (hh > 0 ? sn.V[c].toFixed(3) : '0.000') + ',' + (isFinite(dm) ? (dm * 1000).toFixed(2) : ''));
      }
    }
    return out.join('\r\n');
  }
  function matCSV(title, arr, fmt) {
    var nx = sim.nx, ny = sim.ny, lines = ['# ' + title + '　1行目が左岸側、左の列が上流　格子間隔 ' + sim.dx + ' m'];
    for (var j = ny - 1; j >= 0; j--) { var r = []; for (var i = 0; i < nx; i++) r.push(fmt(arr[j * nx + i], j * nx + i)); lines.push(r.join(',')); }
    return lines.join('\r\n');
  }
  function profileCSV() {
    if (pathMode()) {   // DEM：選んだ線（本川・支川）に沿って、谷底の帯（±谷筋の幅の半分）の最深値と平均値
      var PSc = pathSamples(), Sn = sim.snaps, hw = Math.max(sim.dx, ((sim.P.dem && sim.P.dem.chW) || 6) / 2 + 0.5 * sim.dx), hd2 = ['s_m', 'x_m', 'y_m'], rw = [];
      Sn.forEach(function (sn) { hd2.push('最深河床高_' + snapLabel(sn.t) + '分'); }); Sn.forEach(function (sn) { hd2.push('平均河床高_' + snapLabel(sn.t) + '分'); });
      for (var k = 0; k < PSc.x.length; k++) {
        var cl = [], seen = {}, nxv = -PSc.ty[k], nyv = PSc.tx[k];
        for (var u = -hw; u <= hw + 1e-9; u += sim.dx / 2) { var c0 = cellXY(PSc.px[k] + nxv * u, PSc.py[k] + nyv * u); if (inBasin(c0) && !seen[c0] && sim.kind[c0] !== 2) { seen[c0] = 1; cl.push(c0); } }
        var mi = [], me = [];
        Sn.forEach(function (sn) { var mn = Infinity, sm = 0; cl.forEach(function (c1) { var zz = sim.z0[c1] + sn.dz[c1]; if (zz < mn) mn = zz; sm += zz; }); mi.push(cl.length ? mn.toFixed(3) : ''); me.push(cl.length ? (sm / cl.length).toFixed(3) : ''); });
        rw.push([PSc.x[k].toFixed(2), PSc.px[k].toFixed(2), PSc.py[k].toFixed(2)].concat(mi, me).join(','));
      }
      return ['# 時刻別の縦断（' + DP[PI].name + '）。s は線の上流端からの距離、x・y は計算領域の座標。谷底の帯のセルの最深値と平均値。堰堤のセルは除く', hd2.join(',')].concat(rw).join('\r\n');
    }
    var nx = sim.nx, ny = sim.ny, S = sim.snaps, z0 = sim.z0, hd = ['x_m'], rows = [];
    S.forEach(function (sn) { hd.push('最深河床高_' + snapLabel(sn.t) + '分'); });
    S.forEach(function (sn) { hd.push('平均河床高_' + snapLabel(sn.t) + '分'); });
    for (var i = 0; i < nx; i++) {
      var r = [((i + 0.5) * sim.dx).toFixed(2)], mins = [], means = [];
      for (var q = 0; q < S.length; q++) {
        var mn = Infinity, sm = 0, n = 0;
        for (var j = 0; j < ny; j++) { var c = j * nx + i; if (!sim.fl[c] || sim.kind[c] !== 0) continue; var zz = z0[c] + S[q].dz[c]; if (zz < mn) mn = zz; sm += zz; n++; }
        mins.push(n ? mn.toFixed(3) : ''); means.push(n ? (sm / n).toFixed(3) : '');
      }
      rows.push(r.concat(mins, means).join(','));
    }
    return ['# 時刻別の縦断。谷底（移動床）のセルの最深値と平均値。堰堤の列は空欄', hd.join(',')].concat(rows).join('\r\n');
  }
  function seriesCSV() {
    var S = sim.snaps, nd = sim.dams.length, hd = ['時刻_分', '流入流量_m3/s', '流出流量_m3/s', '流入土砂量_m3/s', '流出土砂量_m3/s', '累計流入土砂_m3', '累計流出土砂_m3', '堆砂量合計_m3'];
    for (var q = 0; q < nd; q++) hd.push('堰堤' + (q + 1) + 'の上流_m3');
    hd.push('最下流堰堤より下流_m3');
    var rows = S.map(function (sn) {
      var r = [(sn.t / 60).toFixed(2), sn.Qin.toFixed(3), sn.Qout.toFixed(3), sn.QsIn.toExponential(4), sn.QsOut.toExponential(4), sn.sedIn.toFixed(2), sn.sedOut.toFixed(2), sn.stor.toFixed(2)];
      for (var k = 0; k <= nd; k++) r.push(sn.byDam[k].toFixed(2));
      return r.join(',');
    });
    return ['# 時系列の集計。土砂量は実質体積、堆砂量は空隙込み（河床の上昇・低下の体積）', hd.join(',')].concat(rows).join('\r\n');
  }
  // 圧縮しない zip（ファイル名は UTF-8）
  var CRC = null;
  function crc32(b) {
    if (!CRC) { CRC = new Uint32Array(256); for (var n = 0; n < 256; n++) { var c = n; for (var k = 0; k < 8; k++) c = c & 1 ? 0xEDB88320 ^ (c >>> 1) : c >>> 1; CRC[n] = c >>> 0; } }
    var x = 0xFFFFFFFF; for (var i = 0; i < b.length; i++) x = CRC[(x ^ b[i]) & 0xFF] ^ (x >>> 8);
    return (x ^ 0xFFFFFFFF) >>> 0;
  }
  function makeZip(files) {
    var enc = new TextEncoder(), parts = [], central = [], off = 0;
    var d = new Date(), dt = ((d.getHours() << 11) | (d.getMinutes() << 5) | (d.getSeconds() >> 1)) & 0xFFFF, dd = (((d.getFullYear() - 1980) << 9) | ((d.getMonth() + 1) << 5) | d.getDate()) & 0xFFFF;
    files.forEach(function (f) {
      var nm = enc.encode(f.name), data = typeof f.data === 'string' ? enc.encode(f.data) : f.data, crc = crc32(data);
      var h = new DataView(new ArrayBuffer(30));
      h.setUint32(0, 0x04034b50, true); h.setUint16(4, 20, true); h.setUint16(6, 0x0800, true); h.setUint16(8, 0, true);
      h.setUint16(10, dt, true); h.setUint16(12, dd, true); h.setUint32(14, crc, true); h.setUint32(18, data.length, true); h.setUint32(22, data.length, true);
      h.setUint16(26, nm.length, true); h.setUint16(28, 0, true);
      parts.push(new Uint8Array(h.buffer), nm, data);
      var c = new DataView(new ArrayBuffer(46));
      c.setUint32(0, 0x02014b50, true); c.setUint16(4, 20, true); c.setUint16(6, 20, true); c.setUint16(8, 0x0800, true); c.setUint16(10, 0, true);
      c.setUint16(12, dt, true); c.setUint16(14, dd, true); c.setUint32(16, crc, true); c.setUint32(20, data.length, true); c.setUint32(24, data.length, true);
      c.setUint16(28, nm.length, true); c.setUint32(42, off, true);
      central.push(new Uint8Array(c.buffer), nm);
      off += 30 + nm.length + data.length;
    });
    var csz = central.reduce(function (s, a) { return s + a.length; }, 0), e = new DataView(new ArrayBuffer(22));
    e.setUint32(0, 0x06054b50, true); e.setUint16(8, files.length, true); e.setUint16(10, files.length, true); e.setUint32(12, csz, true); e.setUint32(16, off, true);
    var all = parts.concat(central, [new Uint8Array(e.buffer)]), tot = all.reduce(function (s, a) { return s + a.length; }, 0), out = new Uint8Array(tot), p = 0;
    all.forEach(function (a) { out.set(a, p); p += a.length; });
    return out;
  }
  function snapZip() {
    var S = sim.snaps, files = [], B = '\uFEFF';
    files.push({ name: '読み方.txt', data: '時間ごとの結果（連続砂防堰堤 河床変動シミュレータ）\r\n\r\n' +
      '集計_時系列.csv：記録した時刻ごとの流量・流砂量・累計土砂量・堰堤ごとの堆砂量\r\n' +
      '縦断_時刻別.csv：時刻別の縦断（谷底の最深河床高と平均河床高）\r\n' +
      '河床高/河床高_XXX分.csv：その時刻の河床高の行列（1行目が左岸側、左の列が上流）。堰堤の諸元をコメント行に含むので、「格子の CSV」で読み込むとその時刻の河床から計算できます\r\n' +
      '河床変動量/河床変動量_XXX分.csv：計算開始時の河床からの変化 [m]（＋が堆積）\r\n' +
      '水深/水深_XXX分.csv：水深 [m]\r\n\r\n時刻は出水の時刻（河床変動の時間軸）です。記録の間隔 ' + +((sim.P.snapEvery || 0) / 60).toFixed(2) + ' 分、格子間隔 ' + sim.dx + ' m。\r\n' });
    files.push({ name: '集計_時系列.csv', data: B + seriesCSV() });
    files.push({ name: '縦断_時刻別.csv', data: B + profileCSV() });
    S.forEach(function (sn) {
      var lb = snapLabel(sn.t), tm = (sn.t / 60).toFixed(1);
      files.push({ name: '河床高/河床高_' + lb + '分.csv', data: B + bedCSV(snapZ(sn), sn.t).text });
      files.push({ name: '河床変動量/河床変動量_' + lb + '分.csv', data: B + matCSV('河床変動量 [m]（出水の時刻 ' + tm + ' 分、＋が堆積）', sn.dz, function (v) { return v.toFixed(4); }) });
      files.push({ name: '水深/水深_' + lb + '分.csv', data: B + matCSV('水深 [m]（出水の時刻 ' + tm + ' 分）', sn.h, function (v) { return (v > 1e-3 ? v : 0).toFixed(3); }) });
    });
    return makeZip(files);
  }
  $('o_ev').onchange = function () {
    var v = Math.max(0, +this.value || 0); this.value = v;
    if (sim) { sim.setSnapEvery(v * 60); if (P) P.snapEvery = v * 60; snapInfo(); }
  };
  $('o_long').onclick = function () {
    if (!sim || !sim.snaps.length) return;
    $('o_msg').hidden = false; $('o_msg').textContent = '書き出しています…';
    setTimeout(function () { saveText('時間ごとの結果_' + (sim.P.snapEvery / 60) + '分間隔.csv', longCSV(), 'o_msg'); }, 30);
  };
  $('o_zip').onclick = function () {
    if (!sim || !sim.snaps.length) return;
    $('o_msg').hidden = false; $('o_msg').textContent = 'まとめています…';
    setTimeout(function () { saveText('時間ごとの結果_' + (sim.P.snapEvery / 60) + '分間隔.zip', snapZip(), 'o_msg', true); }, 30);
  };

  // ------------------------------------------------------------ 粒度分布：代表的な河床材料、平均粒径と幅から作る、要約と粒径加積曲線
  var GR_PRE = { sand: [0.5, 1.5, 3, 1], sg: [5, 4, 5, 1], gravel: [30, 3.5, 5, 1], debris: [50, 5, 6, 1], uni: [5, 1, 1, 1] };
  function erf(x) { var t = 1 / (1 + 0.3275911 * Math.abs(x)), y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x); return x >= 0 ? y : -y; }
  function lnCdf(x, mu, sg) { return 0.5 * (1 + erf((x - mu) / (Math.SQRT2 * sg))); }
  // 対数正規分布を K 階に分ける（ln d の ±2.5σ を等分）。粒径は各階の幾何中心、割合は区間の確率。算術平均が dm になるよう粒径をそろえて縮める
  function makeGrains(dm, sgG, K, fr) {
    K = Math.max(1, Math.min(6, K | 0)); dm = Math.max(0.05, dm) / 1000;
    if (K === 1 || !(sgG > 1.0001)) return [{ d: dm, bed: 1, feed: 1 }];
    var sg = Math.log(sgG), lo = -2.5 * sg, w = 5 * sg / K, d = [], pb = [], pf = [], k, sb = 0, sf = 0, mf = Math.log(fr > 0 ? fr : 1);
    for (k = 0; k < K; k++) {
      var a = lo + k * w, b = a + w; d.push(Math.exp(a + w / 2));
      var e0 = k === 0 ? -1e9 : a, e1 = k === K - 1 ? 1e9 : b;   // 両端の階は裾も含める
      pb.push(lnCdf(e1, 0, sg) - lnCdf(e0, 0, sg)); pf.push(lnCdf(e1, mf, sg) - lnCdf(e0, mf, sg));
      sb += pb[k]; sf += pf[k];
    }
    var mean = 0; for (k = 0; k < K; k++) { pb[k] /= sb; pf[k] /= sf; mean += pb[k] * d[k]; }
    return d.map(function (dk, q) { return { d: dk * dm / mean, bed: pb[q], feed: pf[q] }; });
  }
  function readGrTable() {
    var g = [];
    Array.prototype.forEach.call($('grtab').tBodies[0].rows, function (tr) { var v = tr.querySelectorAll('input'), d = +v[0].value, b = +v[1].value, f = +v[2].value; if (d > 0) g.push({ d: d / 1000, bed: b > 0 ? b : 0, feed: f > 0 ? f : 0 }); });
    return g.sort(function (a, b) { return a.d - b.d; });
  }
  // 割合 p（合計 1）の d_x：累積割合を各階の粒径で ln d に直線補間（階の粒径を代表値とする）
  function dPct(g, key, x) {
    var tot = 0; g.forEach(function (q) { tot += q[key]; }); if (!(tot > 0)) return NaN;
    var cum = 0, pc = [], ld = g.map(function (q) { return Math.log(q.d); });
    g.forEach(function (q) { cum += q[key] / tot; pc.push(cum - 0.5 * q[key] / tot); });
    if (x <= pc[0]) return g[0].d; if (x >= pc[pc.length - 1]) return g[g.length - 1].d;
    for (var k = 1; k < pc.length; k++) if (x <= pc[k]) return Math.exp(ld[k - 1] + (ld[k] - ld[k - 1]) * (x - pc[k - 1]) / (pc[k] - pc[k - 1] || 1));
    return g[g.length - 1].d;
  }
  function grSummary() {
    var g = readGrTable(); if (!g.length) { $('gr_sum').textContent = ''; return; }
    function line(key, name) {
      var tot = 0, m = 0; g.forEach(function (q) { tot += q[key]; m += q[key] * q.d; }); if (!(tot > 0)) return name + '：割合が 0 です';
      var mm = function (v) { v *= 1000; return v >= 10 ? v.toFixed(0) : v >= 1 ? v.toFixed(1) : v.toFixed(2); };
      return name + '：平均粒径 dm ' + mm(m / tot) + ' mm、d16 ' + mm(dPct(g, key, 0.16)) + '、d50 ' + mm(dPct(g, key, 0.5)) + '、d84 ' + mm(dPct(g, key, 0.84)) + '、d90 ' + mm(dPct(g, key, 0.9)) + ' mm（合計 ' + tot.toFixed(1) + ' %）';
    }
    $('gr_sum').textContent = line('bed', '河床') + '　／　' + line('feed', '供給') + '。割合の合計は 100 % でなくても比で使います。';
    // 粒径加積曲線（横軸は対数）
    var cv = $('gr_cv'), dpr = window.devicePixelRatio || 1, W = cv.clientWidth || 300, H = 150; cv.width = W * dpr; cv.height = H * dpr;
    var ctx = cv.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, W, H);
    var lmin = Math.log10(g[0].d * 1000) - 0.3, lmax = Math.log10(g[g.length - 1].d * 1000) + 0.3, L = 34, R = 8, T = 8, B = 22;
    var X = function (dmm) { return L + (Math.log10(dmm) - lmin) / (lmax - lmin) * (W - L - R); }, Y = function (p) { return T + (1 - p) * (H - T - B); };
    ctx.strokeStyle = css('--line'); ctx.fillStyle = css('--ink3'); ctx.font = '10px sans-serif'; ctx.lineWidth = 1;
    for (var e = Math.ceil(lmin); e <= Math.floor(lmax); e++) { var xx = X(Math.pow(10, e)); ctx.beginPath(); ctx.moveTo(xx, T); ctx.lineTo(xx, H - B); ctx.stroke(); ctx.fillText((e < 0 ? Math.pow(10, e).toFixed(-e) : Math.pow(10, e)) + ' mm', xx - 10, H - 6); }
    [0, 0.5, 1].forEach(function (p) { ctx.beginPath(); ctx.moveTo(L, Y(p)); ctx.lineTo(W - R, Y(p)); ctx.stroke(); ctx.fillText((p * 100) + '%', 2, Y(p) + 3); });
    [['bed', css('--sed')], ['feed', css('--s2')]].forEach(function (kc) {
      var tot = 0; g.forEach(function (q) { tot += q[kc[0]]; }); if (!(tot > 0)) return;
      ctx.strokeStyle = kc[1]; ctx.lineWidth = 2; ctx.beginPath(); var cum = 0;
      g.forEach(function (q, k) { cum += q[kc[0]] / tot; var xx = X(q.d * 1000), yy = Y(cum); if (k) ctx.lineTo(xx, yy); else ctx.moveTo(xx, yy); });
      ctx.stroke();
    });
    ctx.fillStyle = css('--sed'); ctx.fillText('— 河床', W - 110, H - B - 6); ctx.fillStyle = css('--s2'); ctx.fillText('— 供給', W - 56, H - B - 6);
  }
  function grRow(g) {
    var tr = document.createElement('tr');
    tr.innerHTML = '<td><input type="number" step="1" min="0.1" value="' + (+(g.d * 1000).toFixed(2)) + '" aria-label="粒径"></td>' +
      '<td><input type="number" step="5" min="0" value="' + (+(g.bed * 100).toFixed(1)) + '" aria-label="河床の割合"></td>' +
      '<td><input type="number" step="5" min="0" value="' + (+(g.feed * 100).toFixed(1)) + '" aria-label="供給の割合"></td>' +
      '<td><button class="x" aria-label="削除">×</button></td>';
    tr.querySelector('button').onclick = function () { if ($('grtab').tBodies[0].rows.length > 1) { tr.remove(); grSummary(); } };
    Array.prototype.forEach.call(tr.querySelectorAll('input'), function (i) { i.addEventListener('input', grSummary); });
    return tr;
  }
  function wRow(w) {
    var tr = document.createElement('tr');
    tr.innerHTML = '<td><input type="number" step="10" min="0" value="' + w[0] + '" aria-label="位置"></td><td><input type="number" step="2" min="4" value="' + w[1] + '" aria-label="川幅"></td><td><button class="x" aria-label="削除">×</button></td>';
    tr.querySelector('button').onclick = function () { tr.remove(); };
    return tr;
  }
  function iRow(w) {
    var tr = document.createElement('tr');
    tr.innerHTML = '<td><input type="number" step="10" min="0" value="' + w[0] + '" aria-label="位置"></td><td><input type="number" step="1" min="2" value="' + w[1] + '" aria-label="勾配 1/n"></td><td><button class="x" aria-label="削除">×</button></td>';
    tr.querySelector('button').onclick = function () { tr.remove(); };
    return tr;
  }
  function tRow(t) {
    var tr = document.createElement('tr');
    tr.innerHTML = '<td><input type="number" step="10" min="0" value="' + t.x + '" aria-label="合流点"></td>' +
      '<td><select aria-label="岸"><option value="L"' + (t.side !== 'R' ? ' selected' : '') + '>左岸</option><option value="R"' + (t.side === 'R' ? ' selected' : '') + '>右岸</option></select></td>' +
      '<td><input type="number" step="1" min="2" value="' + t.B + '" aria-label="幅"></td>' +
      '<td><input type="number" step="1" min="1.5" value="' + t.I + '" aria-label="勾配 1/n"></td>' +
      '<td><input type="number" step="0.1" min="0" value="' + t.r + '" aria-label="流量比"></td>' +
      '<td><input type="number" step="0.1" min="0" value="' + t.rs + '" aria-label="土砂量比"></td>' +
      '<td><input type="number" step="10" min="0" value="' + (t.lag || 0) + '" aria-label="遅れ"></td>' +
      '<td><input type="number" step="5" min="30" max="150" value="' + (isFinite(t.ang) ? t.ang : 90) + '" aria-label="合流角"></td>' +
      '<td><button class="x" aria-label="削除">×</button></td>';
    tr.querySelector('button').onclick = function () { tr.remove(); };
    return tr;
  }
  function setT(list) { var tb = $('ttab').tBodies[0]; tb.innerHTML = ''; list.forEach(function (t) { tb.appendChild(tRow(t)); }); }
  function setI(list) { var tb = $('itab').tBodies[0]; tb.innerHTML = ''; list.forEach(function (w) { tb.appendChild(iRow(w)); }); }
  function bRow(w) {
    var tr = document.createElement('tr');
    tr.innerHTML = '<td><input type="number" step="10" value="' + w[0] + '" aria-label="位置"></td><td><input type="number" step="5" value="' + w[1] + '" aria-label="中心線の横ずれ"></td>' +
      '<td><input type="number" step="2" min="4" value="' + (w[2] > 0 ? w[2] : '') + '" placeholder="–" aria-label="川幅"></td><td><button class="x" aria-label="削除">×</button></td>';
    tr.querySelector('button').onclick = function () { tr.remove(); };
    return tr;
  }
  function fRow(f) {
    var tr = document.createElement('tr');
    tr.innerHTML = '<td><input type="number" step="10" min="0" value="' + f[0] + '" aria-label="位置"></td>' +
      '<td style="min-width:3.6em"><select aria-label="岸" style="width:100%;min-width:3.6em;padding-left:4px;padding-right:2px"><option value="1"' + (f[1] >= 0 ? ' selected' : '') + '>左</option><option value="-1"' + (f[1] < 0 ? ' selected' : '') + '>右</option></select></td>' +
      '<td><input type="number" step="2" min="2" value="' + f[2] + '" aria-label="半径"></td><td><input type="number" step="0.5" min="0.2" value="' + f[3] + '" aria-label="高さ"></td>' +
      '<td><button class="x" aria-label="削除">×</button></td>';
    tr.querySelector('button').onclick = function () { tr.remove(); };
    return tr;
  }
  function setF(list) { var tb = $('ftab').tBodies[0]; tb.innerHTML = ''; list.forEach(function (f) { tb.appendChild(fRow(f)); }); }
  function setValley(V, wallGrad) {
    $('v_chD').value = V.chD || 0; $('v_chW').value = V.chW || 0; $('v_chL').value = V.chWl || 0; $('v_rel').value = V.relief || 0;
    $('v_wH').value = V.wallH || 0; $('v_wR').value = V.wallRough || 0; $('v_wG').value = wallGrad; setF(V.fans || []);
  }
  // モデル渓流（写真の河道を参考にした架空の地形）
  var PRESETS = {
    // 床固工群が続く、谷底の広い渓流（早月川の床固工群を参考）
    A: { L: 600, B: 60, I: 18, S: 3, grad: 1.0, wing: 1.5, dx: 2.5,
      bend: [[0, 0, 52], [70, 6, 58], [150, -6, 66], [230, -14, 70], [310, -4, 62], [380, 10, 74], [460, 16, 68], [530, 6, 60], [600, 0, 56]],
      V: { chD: 0.7, chW: 14, chWl: 0, relief: 0.4, wallH: 16, wallRough: 3, fans: [[120, -1, 34, 3.2], [340, 1, 40, 3.5], [505, -1, 28, 2.5]] },
      dams: [75, 135, 195, 255, 315, 375, 435, 495, 555].map(function (x) { return { x: x, H: 3, b: 36 }; }), fill: true },
    // 狭窄部に堰堤を置き、その上流に広い堆砂域がある渓流（南股入の堰堤群を参考）
    B: { L: 660, B: 44, I: 28, S: 2, grad: 1.0, wing: 1.5, dx: 2.5,
      bend: [[0, 0, 40], [70, 10, 38], [140, -6, 42], [200, -2, 30], [235, 0, 24], [270, 4, 36], [340, 18, 84], [410, 24, 110], [480, 18, 112], [540, 8, 70], [575, 2, 34], [610, -2, 40], [660, -6, 42]],
      V: { chD: 0.8, chW: 12, chWl: 0, relief: 0.35, wallH: 18, wallRough: 3, fans: [[110, 1, 30, 3], [300, -1, 36, 3.5], [440, 1, 46, 4]] },
      dams: [{ x: 235, H: 5, b: 16 }, { x: 577, H: 9, b: 22 }], fill: true },
    // 以前の単純な地形（直線・一様幅）
    S: { L: 400, B: 50, I: 20, S: 2, grad: 0.8, wing: 1.5, dx: 2, bend: null, V: {},
      dams: [60, 105, 150, 195, 240, 285, 330].map(function (x) { return { x: x, H: 1.2, b: 0 }; }) }
  };
  function trSwitch() { $('ne_box').hidden = $('s_tr').value !== 'noneq'; $('db_box').hidden = $('s_tr').value !== 'debris'; $('c_mf').disabled = $('s_tr').value === 'debris'; }
  function gwSwitch() { $('gw_box').hidden = !$('gw_on').checked; }
  function fillSwitch() { $('d_fill_row').hidden = $('d_fill_note').hidden = !$('d_fill').checked; }
  function applyPreset(k) {
    var R = PRESETS[k]; if (!R) return;
    $('g_src').value = 'param';
    $('g_L').value = R.L; $('g_B').value = R.B; $('g_I').value = R.I; $('g_S').value = R.S; setW([]); setI([]);
    if (R.bend) { $('b_mode').value = 'table'; setB(R.bend); } else { $('b_mode').value = 'none'; setB([[0, 0], [R.L, 0]]); }
    $('b_x1').value = R.L; $('c_dx').value = String(R.dx);
    setValley(R.V, R.grad); setDams(R.dams); $('d_wing').value = R.wing;
    $('d_fill').checked = !!R.fill; fillSwitch(); presetDams = true;
    geoSwitch();
  }
  function setB(list) { var tb = $('btab').tBodies[0]; tb.innerHTML = ''; list.forEach(function (w) { tb.appendChild(bRow(w)); }); }
  function setW(list) { var tb = $('wtab').tBodies[0]; tb.innerHTML = ''; list.forEach(function (w) { tb.appendChild(wRow(w)); }); }
  function setDams(list) { var tb = $('damtab').tBodies[0]; tb.innerHTML = ''; list.forEach(function (d) { tb.appendChild(damRow(d)); }); }
  function fillForm() {
    var H = P.hydro;
    $('h_mode').value = H.mode; $('h_qb').value = H.Qbase; $('h_qp').value = H.Qpeak;
    $('h_tp').value = H.tPeak / 60; $('h_td').value = H.duration / 60; $('h_cs').value = +(H.Cs * 100).toFixed(3); $('h_cs2').value = +(H.Cs * 100).toFixed(3);
    $('h_table').value = '0, 10, 0.05\n60, 60, 0.30\n120, 30, 0.15\n180, 10, 0.05';
    $('g_L').value = P.L; $('g_B').value = P.B; $('g_I').value = +(1 / P.slope).toFixed(2); $('g_S').value = P.crossSlope * 100;
    $('g_P').value = P.perturb * 100; $('g_D').value = P.erodDepth; $('g_seed').value = P.seed;
    setW(P.widthTable || []); setI(P.slopeTable || []);
    $('b_mode').value = 'none'; $('b_amp').value = 25; $('b_wl').value = 160; $('b_x0').value = 40; $('b_x1').value = P.L;
    setB([[0, 0], [P.L, 0]]);
    setDams(P.dams); $('d_wing').value = P.wingHeight;
    setValley(P.valley || {}, P.sideGrad);
    var tb = $('grtab').tBodies[0]; tb.innerHTML = ''; P.grains.forEach(function (g) { tb.appendChild(grRow(g)); }); setTimeout(grSummary, 0);
    $('s_n').value = P.manning; $('s_tc').value = P.tauCm; $('s_tcm').value = P.tcMethod || 'atm'; $('s_hid').value = P.hiding || 'egiazaroff'; tcSwitch(); $('s_por').value = P.porosity; $('s_rep').value = P.reposeDeg;
    var BK = P.bank || {}, NE = P.noneq || {}, BO = P.boulders || {};
    $('s_tr').value = P.transport || 'eq'; $('ne_de').value = NE.de; $('ne_dd').value = NE.dd; $('ne_he').value = NE.hEdge; $('ne_fe').value = NE.fe; $('ne_fd').value = NE.fd;
    $('bk_h').value = BK.erodH; $('bk_w').value = BK.reposeWet; $('bk_d').value = BK.reposeDry; $('bk_e').value = BK.dryEro; $('bk_cm').value = BK.collapse || 'repose';
    var SC = P.sec || {}; $('sc_bed').checked = !!SC.bed; $('sc_mom').checked = !!SC.mom; $('sc_N').value = SC.Nstar || 11.5; $('sc_hr').value = SC.hrMax || 0.05;
    $('bo_L').value = BO.L; $('bo_m').value = BO.mode || 'parallel'; $('bo_s').value = BO.size || 1; trSwitch();
    $('c_dx').value = String(P.dx); $('c_mf').value = P.morfac; $('c_ord').value = String(P.order);
    $('o_ev').value = P.snapEvery == null ? 10 : +(P.snapEvery / 60).toFixed(2);
    $('hd_t0').value = H.t0 == null ? 60 : H.t0; $('hd_tr').value = H.tRise || 30; $('hd_tf').value = H.tFall || 120; $('hd_c').value = H.Cd == null ? 0.4 : H.Cd; $('hd_iu').value = H.Iup || 0;
    $('hd_td').value = H.mode === 'debris' ? H.duration / 60 : 10;
    var OD = P.openDam || {}; $('od_w').value = OD.w || 0.5; $('od_b').value = OD.beta || 1.5; $('od_c').value = OD.Cblk || 0.2; $('od_a').value = OD.alpha || 1; $('od_p').value = OD.phiBlock || 0.05;
    var DB = P.debris || {}; $('db_de').value = DB.de || 0.0007; $('db_dd').value = DB.dd || 0.05; $('db_tp').value = DB.tanPhi || 0.7; $('db_rho').value = DB.rho || 1000;
    var GW = P.gw || {}; $('gw_on').checked = !!GW.on; $('gw_k').value = GW.k || 0; $('gw_kx').value = GW.kMax || 0.05; $('gw_s0').value = GW.s0 || 0; gwSwitch();
    setT(P.tribs || []); trSwitch();
    modeSwitch();
  }
  function num(id, def, lo, hi) { var v = parseFloat($(id).value); if (!isFinite(v)) v = def; if (lo != null && v < lo) v = lo; if (hi != null && v > hi) v = hi; return v; }
  function readForm() {
    var p = S2.defaultParams(), H = p.hydro;
    H.mode = $('h_mode').value;
    H.Qbase = num('h_qb', 10, 0); H.Qpeak = num('h_qp', 60, 0);
    H.duration = num('h_td', 180, 5) * 60; H.tPeak = Math.min(num('h_tp', 60, 1) * 60, H.duration * 0.98);
    H.Cs = (H.mode === 'table' ? num('h_cs2', 0.5, 0) : num('h_cs', 0.5, 0)) / 100;
    if (H.mode === 'debris') {
      H.t0 = num('hd_t0', 60, 0); H.tRise = num('hd_tr', 30, 1); H.tFall = num('hd_tf', 120, 1); H.Cd = num('hd_c', 0.4, 0, 0.55); H.Iup = num('hd_iu', 0, 0);
      H.duration = num('hd_td', 10, 1) * 60; H.Cs = 0; H.tPeak = H.t0 + H.tRise;
    }
    H.table = [];
    if (H.mode === 'table') {
      $('h_table').value.split(/\n/).forEach(function (ln) {
        var a = ln.split(/[,\s\t]+/).filter(function (s) { return s !== ''; }).map(Number);
        if (a.length >= 2 && isFinite(a[0]) && isFinite(a[1])) H.table.push(a.length >= 3 && isFinite(a[2]) ? [a[0] * 60, a[1], a[2]] : [a[0] * 60, a[1]]);
      });
      H.table.sort(function (a, b) { return a[0] - b[0]; });
      var all3 = H.table.every(function (r) { return r.length > 2; });
      if (!all3) H.table = H.table.map(function (r) { return [r[0], r[1]]; });
      if (H.table.length < 2) { H.mode = 'tri'; }
    }
    p.dx = parseFloat($('c_dx').value) || 2;   // 測量値の格子化に使うので先に決める
    p.L = num('g_L', 400, 60, 2000); p.B = num('g_B', 50, 8, 400); p.slope = 1 / num('g_I', 20, 3, 1000);
    p.crossSlope = num('g_S', 2, 0, 30) / 100; p.perturb = num('g_P', 5, 0, 100) / 100; p.erodDepth = num('g_D', 3, 0.2, 30); p.seed = num('g_seed', 1, 1) | 0;
    p.widthTable = [];
    Array.prototype.forEach.call($('wtab').tBodies[0].rows, function (tr) {
      var a = tr.querySelectorAll('input'), x = parseFloat(a[0].value), w = parseFloat(a[1].value);
      if (isFinite(x) && w > 0) p.widthTable.push([x, Math.min(400, Math.max(4, w))]);
    });
    p.slopeTable = [];
    Array.prototype.forEach.call($('itab').tBodies[0].rows, function (tr) {
      var a = tr.querySelectorAll('input'), x = parseFloat(a[0].value), n = parseFloat(a[1].value);
      if (isFinite(x) && n >= 2) p.slopeTable.push([x, Math.min(5000, n)]);
    });
    var Bmax = p.B; if (p.widthTable.length) { Bmax = 0; p.widthTable.forEach(function (r) { if (r[1] > Bmax) Bmax = r[1]; }); }
    p.bend = { mode: $('b_mode').value, amp: num('b_amp', 25, 0, 500), wl: num('b_wl', 160, 20, 5000), x0: num('b_x0', 0, 0), x1: num('b_x1', p.L, 0), table: [] };
    Array.prototype.forEach.call($('btab').tBodies[0].rows, function (tr) {
      var a = tr.querySelectorAll('input'), x = parseFloat(a[0].value), y = parseFloat(a[1].value), w = parseFloat(a[2].value);
      if (isFinite(x) && isFinite(y)) p.bend.table.push(w > 0 ? [x, y, Math.min(400, Math.max(4, w))] : [x, y]);
    });
    if (p.bend.mode === 'table') { var bws = p.bend.table.filter(function (r) { return r[2] > 0; }); if (bws.length) { Bmax = 0; bws.forEach(function (r) { if (r[2] > Bmax) Bmax = r[2]; }); } }
    p.sideGrad = num('v_wG', 0.8, 0.2, 3);
    p.valley = { chD: num('v_chD', 0, 0, 10), chW: num('v_chW', 0, 0, 400), chWl: num('v_chL', 0, 0, 10000), relief: num('v_rel', 0, 0, 5),
      wallH: num('v_wH', 0, 0, 200), wallRough: num('v_wR', 0, 0, 50), fans: [] };
    Array.prototype.forEach.call($('ftab').tBodies[0].rows, function (tr) {
      var a = tr.querySelectorAll('input'), x = parseFloat(a[0].value), sd = +tr.querySelector('select').value, R = parseFloat(a[1].value), Hf = parseFloat(a[2].value);
      if (isFinite(x) && R > 0 && Hf > 0) p.valley.fans.push([x, sd, R, Hf]);
    });
    Bmax += 2 * p.valley.wallRough;
    var spread = 0;
    if (p.bend.mode === 'sine') spread = 2 * p.bend.amp;
    else if (p.bend.mode === 'table' && p.bend.table.length) { var ys = p.bend.table.map(function (r) { return r[1]; }); spread = Math.max.apply(null, ys) - Math.min.apply(null, ys); }
    Bmax += spread * 1.2;
    p.terrain = null;
    if ($('g_src').value === 'survey') {
      var sr = surveyRows(p.dx); if (sr.err) throw new Error(sr.err);
      p.terrain = { rows: sr.rows, cs: sr.cs }; p.floorH = num('sv_fh', 3, 0.1, 100);
      p.L = sr.nx * sr.cs; Bmax = sr.ny * sr.cs; p.sideCells = 0;
    }
    if ($('g_src').value === 'csv') {
      var cr = csvParsed(); if (cr.err) throw new Error(cr.err);
      p.terrain = { rows: cr.rows, cs: cr.cs }; p.floorH = num('csv_fh', 3, 0.1, 100);
      if (!$('csv_tr').checked && !$('csv_flip').checked) { var mk = damsFromComments($('csv_text').value).mask; if (mk && mk.length === cr.nr && mk[0].length === cr.nc) p.terrain.mask = mk; }
      var cI = parseFloat($('csv_I').value); p.csvSlope = cI >= 3 ? 1 / cI : 0;
      p.L = cr.nc * cr.cs; Bmax = cr.nr * cr.cs; p.sideCells = 0;
    }
    p.dams = [];
    Array.prototype.forEach.call($('damtab').tBodies[0].rows, function (tr) {
      var a = tr.querySelectorAll('input'); var x = parseFloat(a[0].value), h = parseFloat(a[1].value), b = parseFloat(a[2].value);
      if (!(isFinite(x) && h > 0)) return;
      var D = { x: x, H: h, b: isFinite(b) && b > 0 ? b : 0 };
      var ty = tr.querySelector('select').value;
      if (ty && ty !== 'closed') { D.type = ty; var lv = parseFloat(a[3].value), nv = parseFloat(a[4].value); D.l = lv > 0 ? lv : (ty === 'grid' || ty === 'sslit' ? 1.0 : 2.0); D.n = nv >= 1 ? Math.round(nv) : (ty === 'mslit' ? 3 : 1); if (tr.dataset.hb) D.hb = +tr.dataset.hb; }
      if ($('d_adv').checked) {
        var ang = parseFloat(a[5].value), y = parseFloat(a[6].value), zc = parseFloat(a[7].value), len = parseFloat(a[8].value);
        if (isFinite(ang)) D.ang = Math.max(-80, Math.min(80, ang)); if (isFinite(y)) D.y = y; if (isFinite(zc)) D.zc = zc; if (len > 0) D.len = len;
        var so = parseFloat(tr.dataset.so); if (isFinite(so) && Math.abs(D.b - parseFloat(tr.dataset.b0)) < 1e-6) D.so = so;   // 水通し幅を変えたら位置の固定は外す
      }
      p.dams.push(D);
    });
    p.wingHeight = num('d_wing', 1.5, 0, 20);
    p.prefill = $('d_fill').checked ? num('d_fill_r', 0.5, 0.05, 1) : 0;
    var mxH = 0; p.dams.forEach(function (d) { if (d.H > mxH) mxH = d.H; });
    p.grains = [];
    Array.prototype.forEach.call($('grtab').tBodies[0].rows, function (tr) {
      var a = tr.querySelectorAll('input'); var d = parseFloat(a[0].value), b = parseFloat(a[1].value), f = parseFloat(a[2].value);
      if (d > 0) p.grains.push({ d: d / 1000, bed: b > 0 ? b / 100 : 0, feed: f > 0 ? f / 100 : 0 });
    });
    if (!p.grains.length) p.grains = S2.defaultParams().grains;
    p.grains.sort(function (a, b) { return a.d - b.d; });
    p.manning = num('s_n', 0.05, 0.01, 0.3); p.tauCm = num('s_tc', 0.05, 0.01, 0.2); p.tcMethod = $('s_tcm').value; p.hiding = $('s_hid').value; p.porosity = num('s_por', 0.4, 0.1, 0.7); p.reposeDeg = num('s_rep', 35, 15, 60);
    p.transport = $('s_tr').value === 'noneq' ? 'noneq' : $('s_tr').value === 'debris' ? 'debris' : 'eq';
    if (H.mode === 'debris') p.transport = 'debris';
    p.debris = { de: num('db_de', 0.0007, 0, 1), dd: num('db_dd', 0.05, 0, 10), tanPhi: num('db_tp', 0.7, 0.3, 1.2), rho: num('db_rho', 1000, 1000, 1600), Vmin: 0.1, manningFloor: true };
    p.openDam = { w: num('od_w', 0.5, 0.05, 5), beta: num('od_b', 1.5, 0.5, 5), Cblk: num('od_c', 0.2, 0.01, 0.6), alpha: num('od_a', 1, 0.5, 3), phiBlock: num('od_p', 0.05, 0.001, 0.5) };
    p.gw = { on: $('gw_on').checked, k: num('gw_k', 0, 0, 1), kMax: num('gw_kx', 0.05, 0.0001, 1), s0: num('gw_s0', 0, 0, 1) };
    p.tribs = [];
    Array.prototype.forEach.call($('ttab').tBodies[0].rows, function (tr) {
      var a = tr.querySelectorAll('input'), x = parseFloat(a[0].value), B = parseFloat(a[1].value), I = parseFloat(a[2].value);
      if (!(isFinite(x) && B > 0 && I > 0)) return;
      p.tribs.push({ x: x, side: tr.querySelector('select').value, B: B, I: I, len: 40, r: Math.max(0, parseFloat(a[3].value) || 0), rs: Math.max(0, parseFloat(a[4].value) || 0), lag: Math.max(0, parseFloat(a[5].value) || 0), ang: Math.max(30, Math.min(150, parseFloat(a[6].value) || 90)) });
    });
    p.noneq = { de: num('ne_de', 0.6, 0.001, 100), dd: num('ne_dd', 0.6, 0.001, 100), hEdge: num('ne_he', 0, 0, 10), fe: num('ne_fe', 2, 1, 100), fd: num('ne_fd', 10, 1, 1000) };
    p.bank = { erodH: num('bk_h', 0, 0, 100), reposeWet: num('bk_w', 30, 15, 60), reposeDry: num('bk_d', 50, 15, 80), dryEro: num('bk_e', 0, 0, 1), collapse: $('bk_cm').value };
    p.sec = { bed: $('sc_bed').checked, mom: $('sc_mom').checked, Nstar: num('sc_N', 11.5, 0, 50), hrMax: num('sc_hr', 0.05, 0.005, 0.2) };
    p.boulders = { L: num('bo_L', 0, 0, 10000), mode: $('bo_m').value, size: num('bo_s', 1, 1, 4) };
    p.morfac = num('c_mf', 4, 1, 50); p.order = +$('c_ord').value || 2;
    p.snapEvery = num('o_ev', 10, 0, 100000) * 60;
    if ($('g_src').value === 'dem') { var dr = window.SaboDEMUI.params(p); Bmax = dr.Bmax; }
    // 格子数の上限（ブラウザで動く範囲）
    var sideN = function () { return p.terrain ? 0 : Math.max(p.sideCells, p.valley.wallH > 0 ? Math.ceil(Math.max(p.valley.wallH, mxH + p.wingHeight + 3) / p.sideGrad / p.dx) + 1 : 0); };
    var nc = (p.L / p.dx) * (Bmax / p.dx + 2 * sideN());
    while (nc > ($('g_src').value === 'dem' ? 60000 : 40000)) { p.dx *= 1.25; nc = (p.L / p.dx) * (Bmax / p.dx + 2 * sideN()); }
    return p;
  }
  function geoSwitch() {
    var src = $('g_src').value, csv = src === 'csv', sv = src === 'survey', dm = src === 'dem', bm = $('b_mode').value;
    $('g_param').hidden = csv || sv || dm; $('g_csv').hidden = !csv; $('g_survey').hidden = !sv;
    $('b_sine').hidden = bm !== 'sine'; $('b_tab').hidden = bm !== 'table'; $('b_note').hidden = bm === 'none';
  }
  // ---- CSV 地形の読み込み
  function uniq(a) {
    a = a.slice().sort(function (p, q) { return p - q; });
    var tol = 1e-6 * Math.max(1, Math.abs(a[a.length - 1] - a[0])), out = [a[0]];
    for (var i = 1; i < a.length; i++) if (a[i] - out[out.length - 1] > tol) out.push(a[i]);
    return out;
  }
  function medDiff(u) { var d = []; for (var i = 1; i < u.length; i++) d.push(u[i] - u[i - 1]); d.sort(function (p, q) { return p - q; }); return d[d.length >> 1]; }
  function fillHoles(g, nr, nc) {
    for (var it = 0; it < 200; it++) {
      var left = 0, nxt = g.map(function (r) { return r.slice(); });
      for (var j = 0; j < nr; j++) for (var i = 0; i < nc; i++) {
        if (g[j][i] === g[j][i]) continue;
        var s = 0, n = 0;
        [[0, 1], [0, -1], [1, 0], [-1, 0]].forEach(function (o) { var jj = j + o[0], ii = i + o[1]; if (jj >= 0 && jj < nr && ii >= 0 && ii < nc && g[jj][ii] === g[jj][ii]) { s += g[jj][ii]; n++; } });
        if (n) nxt[j][i] = s / n; else left++;
      }
      g = nxt; if (!left) break;
    }
    return g;
  }
  function parseCSV(text, cs, transpose, flip) {
    var data = [];
    text.split(/\r?\n/).forEach(function (ln) {
      ln = ln.trim(); if (!ln) return;
      var v = ln.split(/[,\t; ]+/).filter(function (t) { return t !== ''; }).map(Number);
      if (!v.length || v.some(function (x) { return !isFinite(x); })) return;   // 見出し行などは飛ばす
      data.push(v);
    });
    if (data.length < 3) return { err: '数値の行が 3 行以上必要です。' };
    var rows = null, kind = '', note = '';
    if (data.length >= 9 && data.every(function (r) { return r.length === 3; })) {
      var ux = uniq(data.map(function (r) { return r[0]; })), uy = uniq(data.map(function (r) { return r[1]; }));
      if (ux.length >= 3 && uy.length >= 3 && ux.length * uy.length <= data.length * 1.6 + 4) {
        var csx = medDiff(ux), csy = medDiff(uy), x0 = ux[0], y0 = uy[0];
        var nc = Math.round((ux[ux.length - 1] - x0) / csx) + 1, nr = Math.round((uy[uy.length - 1] - y0) / csy) + 1;
        if (nc * nr > 4e6) return { err: '点の数が多すぎます（' + nc + ' × ' + nr + '）。' };
        var g = []; for (var j = 0; j < nr; j++) { g.push(new Array(nc).fill(NaN)); }
        data.forEach(function (r) { var i = Math.round((r[0] - x0) / csx), jj = nr - 1 - Math.round((r[1] - y0) / csy); if (i >= 0 && i < nc && jj >= 0 && jj < nr) g[jj][i] = r[2]; });
        g = fillHoles(g, nr, nc);
        if (Math.abs(csy - csx) > 0.05 * csx) {   // y 方向の間隔を x 方向にそろえる
          var nr2 = Math.max(3, Math.round((nr - 1) * csy / csx) + 1), g2 = [];
          for (var j2 = 0; j2 < nr2; j2++) { var f = j2 * (nr - 1) / (nr2 - 1), a = Math.floor(f), t = f - a, b = Math.min(nr - 1, a + 1), row = []; for (var i2 = 0; i2 < nc; i2++) row.push(g[a][i2] * (1 - t) + g[b][i2] * t); g2.push(row); }
          g = g2; nr = nr2; note = '（y 方向の間隔 ' + csy + ' m を ' + csx + ' m にそろえました）';
        }
        rows = g; cs = csx; kind = 'x, y, z 形式';
      }
    }
    if (!rows) {
      var lens = {}; data.forEach(function (r) { lens[r.length] = (lens[r.length] || 0) + 1; });
      var L = +Object.keys(lens).sort(function (p, q) { return lens[q] - lens[p]; })[0];
      rows = data.filter(function (r) { return r.length >= L; }).map(function (r) { return r.slice(0, L); });
      kind = '行列形式';
    }
    if (transpose) { var t2 = []; for (var i3 = 0; i3 < rows[0].length; i3++) { var r3 = []; for (var j3 = 0; j3 < rows.length; j3++) r3.push(rows[j3][i3]); t2.push(r3); } rows = t2; }
    if (flip) rows = rows.map(function (r) { return r.slice().reverse(); });
    var nrF = rows.length, ncF = rows[0].length;
    if (ncF < 10 || nrF < 3) return { err: '流下方向に 10 点以上、横断方向に 3 点以上が必要です（いまは ' + ncF + ' × ' + nrF + '）。「行と列を入れ替える」も試してください。' };
    if (!(cs > 0)) return { err: '格子間隔を入力してください。' };
    return { rows: rows, cs: cs, kind: kind, nr: nrF, nc: ncF, note: note };
  }
  var csvCache = { key: '', res: null };
  function textHash(t) { var h = 2166136261 >>> 0; for (var i = 0; i < t.length; i++) { h ^= t.charCodeAt(i); h = Math.imul(h, 16777619) >>> 0; } return t.length + ':' + h; }
  function csvParsed() {
    var key = textHash($('csv_text').value) + '|' + $('csv_cs').value + $('csv_tr').checked + $('csv_flip').checked;   // 全文のハッシュ（先頭と末尾が同じ別の地形を取り違えない）
    if (key !== csvCache.key) { csvCache.key = key; csvCache.res = $('csv_text').value.trim() ? parseCSV($('csv_text').value, parseFloat($('csv_cs').value), $('csv_tr').checked, $('csv_flip').checked) : { err: 'CSV のファイルを選ぶか、貼り付けてください。' }; }
    return csvCache.res;
  }
  // ---- 地形ファイルの種類を中身から見分けて、正しい欄へ入れる
  // 横断（x, s, z。同じ x が何行も続く）／縦断（x, z）／格子（行列、または x, y, z の格子点）
  var FILES = { xs: '', long: '', grid: '' }, presetDams = false;
  function classify(text) {
    var rows = [];
    String(text).split(/\r?\n/).forEach(function (ln) {
      ln = ln.trim(); if (!ln || ln.charAt(0) === '#') return;
      var v = ln.split(/[,\t; ]+/).filter(function (t) { return t !== ''; }).map(Number);
      if (v.length && v.every(isFinite)) rows.push(v);
    });
    if (rows.length < 2) return 'none';
    var cnt = {}; rows.forEach(function (r) { cnt[r.length] = (cnt[r.length] || 0) + 1; });
    var m = +Object.keys(cnt).sort(function (a, b) { return cnt[b] - cnt[a]; })[0];
    if (m === 2) return 'long';
    if (m >= 4) return 'grid';
    if (m === 3) {
      var r3 = parseCSV(text, 1, false, false);
      if (!r3.err && r3.kind === 'x, y, z 形式') return 'grid';
      return 'xs';
    }
    return 'none';
  }
  var KNAME = { xs: '横断 CSV（x, s, z）', long: '縦断 CSV（x, z）', grid: '格子 CSV' };
  function routeFile(text, name, asked) {
    var k = classify(text), msg = '';
    if (k === 'none') { var box = asked === 'grid' ? 'csv_info' : 'sv_info'; $(box).textContent = name + '：数値の行が見つかりませんでした。'; return; }
    if (k !== asked) msg = name + ' は中身から' + KNAME[k] + 'と判断したので、' + (k === 'grid' ? '「格子の CSV を読み込む」' : '「横断・縦断の測量値」の' + (k === 'xs' ? '横断' : '縦断')) + 'の欄に入れました。';
    if (k === 'grid') { $('csv_text').value = text; FILES.grid = name; $('g_src').value = 'csv'; }
    else { $(k === 'xs' ? 'sv_xs' : 'sv_long').value = text; FILES[k] = name; $('g_src').value = 'survey'; }
    if ($('d_fill').checked) { $('d_fill').checked = false; fillSwitch(); }
    var damMsg = '', cm = damsFromComments(text);
    if (k === 'grid' && cm.dams.length) {
      setDams(cm.dams); $('d_adv').checked = true; advSwitch(); presetDams = false;
      if (cm.floorH) $('csv_fh').value = cm.floorH;
      damMsg = 'ファイルに書かれた堰堤 ' + cm.dams.length + ' 基を、前と同じ位置・向き・天端標高で表に入れました。';
    } else if (presetDams) { setDams([]); presetDams = false; damMsg = 'モデル渓流の堰堤は、読み込んだ地形に合わないので外しました。堰堤は改めて入れてください。'; }
    geoSwitch();
    if (k === 'grid') { $('csv_msg').hidden = true; csvInfo(); matchDx(); if (msg) $('csv_info').textContent = msg + '　' + $('csv_info').textContent; }
    else { surveyInfo(); if (msg) $('sv_info').textContent = msg + '　' + $('sv_info').textContent; }
    damRangeNote();
    $('file_note').hidden = false; $('file_note').textContent = '読み込んだファイル：' + name + '。' + damMsg + '「適用して最初から」を押すと、この地形で計算をやり直します。';
  }
  // 格子 CSV の間隔が計算の格子幅の選択肢にあれば、そろえる（補間で地形がなまらないように）
  function matchDx() {
    var r = csvParsed(); if (r.err) return;
    var opts = Array.prototype.map.call($('c_dx').options, function (o) { return +o.value; });
    for (var i = 0; i < opts.length; i++) if (Math.abs(opts[i] - r.cs) < 1e-6) { $('c_dx').value = $('c_dx').options[i].value; return; }
  }
  function csvInfo() {
    var r = csvParsed();
    $('csv_info').textContent = r.err ? r.err : '読み込み：' + r.kind + '、' + r.nr + ' 行 × ' + r.nc + ' 列、格子間隔 ' + r.cs + ' m（区間長 ' + (r.nc * r.cs).toFixed(0) + ' m、幅 ' + (r.nr * r.cs).toFixed(0) + ' m）' + r.note;
  }
  // ---- 横断・縦断の測量値
  function surveyRows(dx) {
    return window.SaboSurvey.toRows({ xsText: $('sv_xs').value, longText: $('sv_long').value, dx: dx, align: $('sv_align').value, wall: num('sv_wall', 1, 0.1, 10),
      B: num('g_B', 50, 4, 400), cross: num('g_S', 2, 0, 30) / 100, side: 0.8, Hw: 6 });
  }
  function surveyInfo() {
    var r = surveyRows(parseFloat($('c_dx').value) || 2);
    var fn = (FILES.xs || FILES.long) ? '使う値：横断 ' + ($('sv_xs').value.trim() ? (FILES.xs || '貼り付け・テスト形状') : 'なし') + '、縦断 ' + ($('sv_long').value.trim() ? (FILES.long || '貼り付け・テスト形状') : 'なし') + '。' : '';
    $('sv_info').textContent = r.err ? fn + r.err : fn + '地形：' + r.ny + ' 行 × ' + r.nx + ' 列（区間長 ' + r.L.toFixed(0) + ' m、幅 ' + r.W.toFixed(0) + ' m）。断面 ' + (r.usedTemplate ? 'なし（数値指定の形を使用）' : r.nsec) + '、縦断 ' + (r.hasLong ? 'あり' : 'なし（断面の最深部を直線でつなぐ）') + '。x は最上流を 0 とします。';
    return r;
  }
  function curL() {
    if ($('g_src').value === 'csv') { var r = csvParsed(); if (!r.err) return r.nc * r.cs; }
    if ($('g_src').value === 'survey') { var s = surveyRows(parseFloat($('c_dx').value) || 2); if (!s.err) return s.L; }
    return num('g_L', 400);
  }
  function tcSwitch() {
    var m = $('s_tcm').value; $('s_tc_row').hidden = m !== 'const';
    $('s_tc_note').textContent = m === 'atm' ? 'τ*c = 0.04×10^(1.72 i)。i は各セルの流れの向きの河床勾配です。' : m === 'iwagaki' ? '岩垣式は勾配の効果を含みません。粒径別は「修正 Egiazaroff」と組み合わせるのが一般的です。' : '';
  }
  function modeSwitch() {
    var m = $('h_mode').value; $('h_par').hidden = (m === 'table'); $('h_tab').hidden = (m !== 'table'); $('h_deb').hidden = (m !== 'debris');
    $('h_qb').closest('.row').hidden = (m === 'const'); $('h_tp').closest('.row').hidden = (m === 'const' || m === 'debris');
    $('h_td').closest('.row').hidden = $('h_cs').closest('.row').hidden = (m === 'debris');
    if (m === 'debris' && $('s_tr').value !== 'debris') { $('s_tr').value = 'debris'; trSwitch(); }
  }

  // ------------------------------------------------------------ 色
  function css(n) { return getComputedStyle(document.documentElement).getPropertyValue(n).trim(); }
  function hex(h) { return [parseInt(h.substr(1, 2), 16), parseInt(h.substr(3, 2), 16), parseInt(h.substr(5, 2), 16)]; }
  function ramp(list) {
    var st = list.map(hex), n = st.length - 1;
    return function (t, o) { t = t < 0 ? 0 : t > 1 ? 1 : t; var x = t * n, i = Math.min(n - 1, x | 0), f = x - i, a = st[i], b = st[i + 1]; o[0] = a[0] + (b[0] - a[0]) * f; o[1] = a[1] + (b[1] - a[1]) * f; o[2] = a[2] + (b[2] - a[2]) * f; };
  }
  var seqBlue = ramp(['#e3effc', '#b7d3f6', '#86b6ef', '#5598e7', '#2a78d6', '#1c5cab', '#104281', '#0d366b']);
  var seqOr = ramp(['#fbe9df', '#f6c8b2', '#f0a583', '#eb8457', '#d95926', '#a8431b', '#7a3013']);
  var divBR = ramp(['#1c5cab', '#5598e7', '#b7d3f6', '#e9e8e4', '#f3b6ac', '#e0705c', '#b3302a']);
  var seqEarth = ramp(['#f4eedd', '#e5d5ab', '#cfb67e', '#ae9058', '#8b6d3f', '#674e2c', '#45331d']);
  var divTB = ramp(['#0d4f63', '#2f8a9a', '#9fcfd3', '#efeee8', '#e4c896', '#c18e4c', '#7d521f']);
  var VARS = {
    real: { name: '写実' },
    hs: { name: '陰影＋段彩', terr: true }, red: { name: '赤色立体風', terr: true }, slp: { name: '傾斜', terr: true },
    bed: { name: '河床位（標高）', unit: 'm', ramp: seqEarth },
    bedrel: { name: '河床位（初期勾配を除く）', unit: 'm', ramp: divTB, div: true },
    depth: { name: '水深', unit: 'm', water: true, ramp: seqBlue },
    vel: { name: '流速', unit: 'm/s', water: true, ramp: seqBlue },
    tau: { name: '無次元掃流力 τ*', unit: '', water: true, ramp: seqBlue },
    tcr: { name: '限界掃流力 τ*c（平均粒径）', unit: '', water: true, ramp: seqBlue },
    fr: { name: 'フルード数', unit: '', water: true, ramp: seqBlue },
    dz: { name: '河床変動量', unit: 'm', ramp: divBR, div: true },
    dm: { name: '表層の平均粒径', unit: 'mm', ramp: seqOr },
    conc: { name: '流れの土砂濃度', unit: '', water: true, ramp: seqOr },
    gw: { name: '伏流（堆積層の飽和度）', unit: '', ramp: seqBlue }
  };
  function hash(c) { var x = Math.sin(c * 12.9898 + 78.233) * 43758.5453; return x - Math.floor(x); }
  function trend(x) { return sim.zbAt ? sim.zbAt(x) : sim.trA + sim.trB * x; }
  var trendX = trend;
  function relZ(c) { return sim.z[c] - trend((c % sim.nx + 0.5) * sim.dx); }
  // 表層の粒度：out[0] 平均粒径, out[1] d84 [m]
  function gstat(c, out) {
    var K = sim.K, A = sim.A, d = sim.d, ck = c * K, s = 0, m = 0, k;
    for (k = 0; k < K; k++) { s += A[ck + k]; m += A[ck + k] * d[k]; }
    if (s < 1e-6) return false;
    out[0] = m / s; var acc = 0; out[1] = d[K - 1];
    for (k = 0; k < K; k++) { acc += A[ck + k] / s; if (acc >= 0.84) { out[1] = d[k]; break; } }
    return true;
  }
  function coarse(dm) { var k = Math.log(dm / 0.002) / Math.log(0.04 / 0.002); return k < 0 ? 0 : k > 1 ? 1 : k; }
  function contourInfo(mode) {
    if (!$('v_cont').checked || !sim) return null;
    var rel = mode === 'bedrel', lo = 1e30, hi = -1e30;
    for (var c = 0; c < sim.N; c++) if (sim.fl[c]) { var v = rel ? relZ(c) : sim.z[c]; if (v < lo) lo = v; if (v > hi) hi = v; }
    return { rel: rel, ci: nice((hi - lo) / (rel ? 10 : 22)) };
  }
  function nice(v) { if (!(v > 0)) return 1; var e = Math.pow(10, Math.floor(Math.log10(v))), m = v / e; return (m <= 1 ? 1 : m <= 2 ? 2 : m <= 5 ? 5 : 10) * e; }

  // 変数値（NaN = 対象外）
  function value(mode, c) {
    var h = sim.h[c];
    if (mode === 'depth') return h > 0.01 ? h : NaN;
    if (mode === 'vel') return h > 0.01 ? Math.hypot(sim.qx[c], sim.qy[c]) / h : NaN;
    if (mode === 'fr') return h > 0.01 ? Math.hypot(sim.qx[c], sim.qy[c]) / h / Math.sqrt(G * h) : NaN;
    if (mode === 'tau') {
      if (h < 0.02) return NaN; var dm = sim.surfaceDm(c); if (dm !== dm) dm = sim.dm0;
      var V = Math.hypot(sim.qx[c], sim.qy[c]) / h; return G * sim.n2 * V * V / Math.cbrt(h) / (sim.sg1 * G * dm);
    }
    if (mode === 'dz') return sim.kind[c] === 0 || Math.abs(sim.z[c] - sim.z0[c]) > 1e-3 ? sim.z[c] - sim.z0[c] : NaN;
    if (mode === 'dm') { var d = sim.surfaceDm(c); return d === d ? d * 1000 : NaN; }
    if (mode === 'tcr') { var tc = sim.tcs[c]; return h > 0.02 && tc === tc ? tc : NaN; }
    if (mode === 'bed') return sim.z[c];
    if (mode === 'bedrel') return relZ(c);
    if (mode === 'conc') {   // 土石流：混合物中の容積濃度。侵食・堆積速度式（掃流）：水の中の土砂の容積比
      if (h < 0.02 || !sim.Cc) return NaN; var K = sim.K, sc = 0; for (var k = 0; k < K; k++) sc += sim.Cc[c * K + k];
      return sc / h;
    }
    if (mode === 'gw') { if (!sim.W || sim.kind[c] !== 0) return NaN; var Tq = sim.z[c] - sim.zfix[c]; return Tq > 0.01 ? Math.min(1, sim.W[c] / (sim.P.porosity * Tq)) : NaN; }
    return NaN;
  }
  var vrange = { lo: 0, hi: 1 };
  function computeRange(mode) {
    var V = VARS[mode]; if (!V || !V.ramp) return;
    if (mode === 'dm') {   // 実際の分布に合わせる（谷底のセルの 2〜98 % 点、対数目盛）
      var vals = [];
      for (var c3 = 0; c3 < sim.N; c3++) if (sim.fl[c3] && sim.kind[c3] === 0) { var dv = value('dm', c3); if (dv === dv) vals.push(dv); }
      vals.sort(function (a, b) { return a - b; });
      var dlo = sim.d[0] * 1000, dhi = sim.dmax * 1000;
      if (vals.length) { dlo = Math.max(dlo, vals[Math.floor(vals.length * 0.02)]); dhi = Math.min(dhi, vals[Math.min(vals.length - 1, Math.floor(vals.length * 0.98))]); }
      if (dhi / dlo < 1.5) { var gm = Math.sqrt(dhi * dlo); dlo = gm / Math.sqrt(1.5); dhi = gm * Math.sqrt(1.5); }
      vrange.lo = dlo; vrange.hi = dhi; return;
    }
    if (mode === 'fr') { vrange.lo = 0; vrange.hi = 2; return; }
    if (mode === 'gw') { vrange.lo = 0; vrange.hi = 1; return; }
    if (mode === 'conc') { var mxc = 0; for (var c4 = 0; c4 < sim.N; c4++) { var cv = value('conc', c4); if (cv > mxc) mxc = cv; } vrange.lo = 0; vrange.hi = Math.max(0.01, Math.min(0.6, nice(mxc))); return; }
    if (mode === 'bed' || mode === 'bedrel') {
      var lo = 1e30, hi = -1e30, am = 0;
      for (var c2 = 0; c2 < sim.N; c2++) if (sim.fl[c2]) { var w = value(mode, c2); if (w < lo) lo = w; if (w > hi) hi = w; if (Math.abs(w) > am) am = Math.abs(w); }
      if (mode === 'bed') { vrange.lo = lo; vrange.hi = hi > lo ? hi : lo + 1; } else { var a2 = nice(am * 0.8); vrange.lo = -a2; vrange.hi = a2; }
      return;
    }
    var mx = 0;
    for (var c = 0; c < sim.N; c++) { var v = value(mode, c); if (v === v && Math.abs(v) > mx) mx = Math.abs(v); }
    if (V.div) { var a = Math.max(0.5, nice(mx * 0.8)); vrange.lo = -a; vrange.hi = a; }
    else { vrange.lo = 0; vrange.hi = nice(Math.max(mx * 0.85, mode === 'tau' ? 0.1 : 0.5)); }
  }
  function norm(mode, v) {
    if (mode === 'dm') return Math.log(v / vrange.lo) / Math.log(vrange.hi / vrange.lo || 2);
    return (v - vrange.lo) / (vrange.hi - vrange.lo);
  }
  // 地表の色（写実）
  var gcol = [0, 0, 0];
  // 堰堤のセルの下の地盤：堰堤でない4近傍の平均（なければそのまま）
  function damGround(c, zA) {
    zA = zA || sim.z;
    var nx = sim.nx, i = c % nx, s = 0, n = 0, nb = [i > 0 ? c - 1 : -1, i < nx - 1 ? c + 1 : -1, c - nx, c + nx];
    for (var k = 0; k < 4; k++) { var q = nb[k]; if (q >= 0 && q < sim.N && sim.kind[q] !== 2) { s += zA[q]; n++; } }
    return n ? s / n : zA[c];
  }
  function nonDam(c) {
    var nx = sim.nx, i = c % nx, nb = [i > 0 ? c - 1 : -1, i < nx - 1 ? c + 1 : -1, c - nx, c + nx];
    for (var k = 0; k < 4; k++) { var q = nb[k]; if (q >= 0 && q < sim.N && sim.kind[q] !== 2) return q; }
    return c;
  }
  // 堰堤を軸に沿って3つの部分（左袖・水通し・右袖）に分ける
  function damTypeLabel(D) {
    if (!D.op) return '不透過';
    var nm = { slit: 'スリット', mslit: '複数スリット', grid: '鋼製格子型', sslit: '鋼製スリット型' }[D.op.type] || D.op.type;
    return nm + ' l=' + D.op.l + ' m' + (D.op.type === 'mslit' ? '×' + D.op.n : '') + (D.op.blocked ? '　<b style="color:var(--err,#b3261e)">閉塞（' + (D.op.tBlock / 60).toFixed(1) + ' 分、C=' + D.op.Cblk.toFixed(2) + '、d95=' + (D.op.d95 * 1000).toFixed(0) + ' mm）</b>' : '　開口');
  }
  function damSegs(D) {
    var out = [], op = D.op;
    if (D.spMin - D.sMin > 1e-6) out.push({ s0: D.sMin, s1: D.spMin, top: D.wingTop, spill: false });
    if (op) {   // 透過型：開口の部分を抜き、柱（複数スリット）・鋼製の部材を置く
      var a0 = Math.max(D.spMin, op.s0), a1 = Math.min(D.spMax, op.s1);
      if (a0 - D.spMin > 1e-6) out.push({ s0: D.spMin, s1: a0, top: D.crest, spill: true });
      if (op.type === 'mslit') for (var m = 0; m < op.n - 1; m++) { var p0 = op.s0 + op.l + m * (op.l + op.w); out.push({ s0: p0, s1: p0 + op.w, top: D.crest, spill: true }); }
      if (op.type === 'grid' || op.type === 'sslit') {
        for (var p1 = a0 + op.l; p1 + op.w <= a1 + 1e-6; p1 += op.l + op.w) out.push({ s0: p1, s1: p1 + op.w, top: D.crest, bot: op.zo, steel: true });
        out.push({ s0: a0, s1: a1, top: D.crest, bot: D.crest - op.w, steel: true, horiz: true });                      // 天端の横部材
        if (op.type === 'grid') for (var hz = op.zo + op.l; hz + op.w < D.crest - op.w; hz += op.l + op.w) out.push({ s0: a0, s1: a1, top: hz + op.w, bot: hz, steel: true, horiz: true });
      }
      if (D.spMax - a1 > 1e-6) out.push({ s0: a1, s1: D.spMax, top: D.crest, spill: true });
    } else if (D.spMax - D.spMin > 1e-6) out.push({ s0: D.spMin, s1: D.spMax, top: D.crest, spill: true });
    if (D.sMax - D.spMax > 1e-6) out.push({ s0: D.spMax, s1: D.sMax, top: D.wingTop, spill: false });
    return out;
  }
  function ground(c, o) {
    if (sim.basin && !sim.basin[c]) { o[0] = o[1] = o[2] = 214; return; }   // DEM の流域の外
    var k = sim.kind[c], n = hash(c);
    if (k === 2) { var s = 178 + n * 14; o[0] = s; o[1] = s - 2; o[2] = s - 9; return; }
    if (k === 1 && sim.bare && sim.bare[c] && sim.z[c] - sim.z0[c] < 0.02) { var tb = 0.9 + 0.2 * n; o[0] = 184 * tb; o[1] = 176 * tb; o[2] = 160 * tb; return; }   // 扇状地（裸地）
    if (k === 1 && sim.z[c] - sim.z0[c] < 0.02) { var t = hash(c * 3 + 1); o[0] = 46 + 30 * t + 8 * n; o[1] = 78 + 22 * t + 10 * n; o[2] = 42 + 12 * t; return; }   // 森林
    if (sim.bankC && sim.bankC[c] && Math.abs(sim.z[c] - sim.z0[c]) < 0.05) { var t2 = hash(c * 3 + 1); o[0] = 52 + 30 * t2 + 8 * n; o[1] = 84 + 22 * t2 + 10 * n; o[2] = 44 + 12 * t2; return; }   // まだ削れていない側岸（段丘の森林）
    if (sim.fixedC && sim.fixedC[c]) { var tq = 0.75 + 0.2 * n; o[0] = 92 * tq; o[1] = 88 * tq; o[2] = 82 * tq; return; }   // 巨礫
    var dm = sim.surfaceDm(c);
    if (dm !== dm) { o[0] = 96 + 16 * n; o[1] = 90 + 14 * n; o[2] = 84 + 12 * n; return; }
    var f = Math.log(dm / 0.004) / Math.log(0.15 / 0.004); f = f < 0 ? 0 : f > 1 ? 1 : f;
    var b = 0.88 + 0.24 * n;
    o[0] = (208 - 72 * f) * b; o[1] = (194 - 64 * f) * b; o[2] = (166 - 40 * f) * b;
    if (sim.h[c] > 0.01) { o[0] *= 0.7; o[1] *= 0.72; o[2] *= 0.74; }
  }
  var WCLR = { sh: [96, 172, 170], dp: [22, 78, 104] }, WTURB = { sh: [168, 142, 100], dp: [112, 86, 54] };
  function waterCol() { return $('v_turb').checked ? WTURB : WCLR; }

  // ------------------------------------------------------------ 平面図
  var planOff = null, planImg = null, PB = null, xsX = null;
  function ih(x, y) { var h = (Math.imul(x, 374761393) + Math.imul(y, 668265263)) | 0; h = Math.imul(h ^ (h >>> 13), 1274126177); return ((h ^ (h >>> 16)) >>> 0) / 4294967296; }
  function planSize() {
    var cv = $('plan'), Wp = cv.parentNode.clientWidth, Wd = Math.min(Wp, 640 * sim.Lx / sim.Ly), Hh = Math.max(70, Math.round(Wd * sim.Ly / sim.Lx));
    cv.style.width = Wd + 'px'; return { W: Wd, H: Hh };
  }
  // 地形の色（陰影＋段彩・赤色立体風・傾斜）。いまの地形（DEM の流域の外は元の標高）から、セルごとに RGB を作る
  // 陰影：Horn の 3×3 勾配と 4 方向の光源（北西が主）。尾根・谷：まわり約 75 m の平均標高との差（TPI）
  var TC = null;
  function terrColors(mode) {
    var nx = sim.nx, ny = sim.ny, N = sim.N, dx = sim.dx, kind = sim.kind, i, j, c;
    if (!TC || TC.length !== N * 3) TC = new Float32Array(N * 3);
    var Z = new Float64Array(N), zmin = 1e30, zmax = -1e30;
    for (c = 0; c < N; c++) { var zz = kind[c] === 2 ? damGround(c) : (sim.zOut && !sim.basin[c] ? sim.zOut[c] : sim.z[c]); Z[c] = zz; if (zz < zmin) zmin = zz; if (zz > zmax) zmax = zz; }
    var zr = Math.max(1e-6, zmax - zmin), R = Math.max(2, Math.min(30, Math.round(75 / dx)));
    var W1 = nx + 1, I = new Float64Array(W1 * (ny + 1));
    for (j = 0; j < ny; j++) { var rs = 0; for (i = 0; i < nx; i++) { rs += Z[j * nx + i]; I[(j + 1) * W1 + i + 1] = I[j * W1 + i + 1] + rs; } }
    var L = Math.SQRT1_2, AZ = [[315, 0.55], [270, 0.2], [0, 0.15], [225, 0.1]].map(function (q) { var t = q[0] * Math.PI / 180; return [q[1] * L * Math.sin(t), q[1] * L * Math.cos(t), q[1] * L]; });
    var TT = [[0, 110, 165, 100], [0.25, 175, 200, 125], [0.5, 220, 200, 140], [0.75, 190, 150, 110], [1, 240, 236, 230]];
    function zc(ii, jj) { ii = ii < 0 ? 0 : ii >= nx ? nx - 1 : ii; jj = jj < 0 ? 0 : jj >= ny ? ny - 1 : jj; return Z[jj * nx + ii]; }
    for (j = 0; j < ny; j++) for (i = 0; i < nx; i++) {
      c = j * nx + i; var z0 = Z[c];
      // j が増える向きが北（平面図の上）
      var a1 = zc(i - 1, j + 1), a2 = zc(i, j + 1), a3 = zc(i + 1, j + 1), a4 = zc(i - 1, j), a6 = zc(i + 1, j), a7 = zc(i - 1, j - 1), a8 = zc(i, j - 1), a9 = zc(i + 1, j - 1);
      var gx = ((a3 + 2 * a6 + a9) - (a1 + 2 * a4 + a7)) / (8 * dx), gy = ((a1 + 2 * a2 + a3) - (a7 + 2 * a8 + a9)) / (8 * dx);
      var g2 = gx * gx + gy * gy, nr = 1 / Math.sqrt(1 + g2), hs = 0, k;
      for (k = 0; k < 4; k++) { var h1 = (-gx * AZ[k][0] - gy * AZ[k][1] + AZ[k][2]) * nr; if (h1 > 0) hs += h1; }
      hs = hs > 0.72 ? 1 : hs / 0.72;
      var i0 = Math.max(0, i - R), i1 = Math.min(nx, i + R + 1), j0 = Math.max(0, j - R), j1 = Math.min(ny, j + R + 1);
      var sm = I[j1 * W1 + i1] - I[j0 * W1 + i1] - I[j1 * W1 + i0] + I[j0 * W1 + i0], sn = (i1 - i0) * (j1 - j0);
      var tpi = Math.tanh((z0 - sm / sn) / (0.12 * R * dx)), slp = Math.atan(Math.sqrt(g2)) * 57.3, r, g, b;
      if (mode === 'red') { var br = 0.5 + 0.5 * tpi, sat = slp > 50 ? 1 : slp / 50; r = 255 * br; g = 255 * br * (1 - 0.8 * sat); b = 255 * br * (1 - 0.9 * sat); }
      else if (mode === 'slp') { r = g = b = 255 * (1 - (slp > 55 ? 1 : slp / 55)) * (0.8 + 0.2 * hs); }
      else {
        var t = (z0 - zmin) / zr, q = 1; while (q < 4 && t > TT[q][0]) q++;
        var A0 = TT[q - 1], B0 = TT[q], u = Math.max(0, Math.min(1, (t - A0[0]) / (B0[0] - A0[0]))), sh = (0.12 + hs) * (1 + 0.22 * tpi), w = 0.55;
        r = (w * (A0[1] + u * (B0[1] - A0[1])) + (1 - w) * 225) * sh; g = (w * (A0[2] + u * (B0[2] - A0[2])) + (1 - w) * 225) * sh; b = (w * (A0[3] + u * (B0[3] - A0[3])) + (1 - w) * 225) * sh;
      }
      if (kind[c] === 2) { r = g = b = 150; }
      // DEM の流域の外は薄く
      if (sim.zOut && !sim.basin[c]) { r = 0.55 * r + 110; g = 0.55 * g + 110; b = 0.55 * b + 110; }
      TC[c * 3] = r > 255 ? 255 : r; TC[c * 3 + 1] = g > 255 ? 255 : g; TC[c * 3 + 2] = b > 255 ? 255 : b;
    }
    return TC;
  }
  function drawPlan() {
    var cv = $('plan'), mode = $('v_mode').value, nx = sim.nx, ny = sim.ny, N = sim.N;
    var sz = planSize(), W = sz.W, Hh = sz.H, dpr = Math.min(2, window.devicePixelRatio || 1);
    if (cv.width !== Math.round(W * dpr) || cv.height !== Math.round(Hh * dpr)) { cv.width = Math.round(W * dpr); cv.height = Math.round(Hh * dpr); cv.style.height = Hh + 'px'; }
    var sub = Math.max(1, Math.min(6, Math.round(cv.width / nx))), pw = nx * sub, ph = ny * sub;
    if (!planOff || planOff.width !== pw || planOff.height !== ph) {
      planOff = document.createElement('canvas'); planOff.width = pw; planOff.height = ph; planImg = planOff.getContext('2d').createImageData(pw, ph);
      PB = { n1: new Float32Array(pw * ph), n3: new Float32Array(pw * ph), cr: new Float32Array(pw * ph), v: new Float32Array(pw * ph) };
      for (var py0 = 0; py0 < ph; py0++) for (var px0 = 0; px0 < pw; px0++) { var q0 = py0 * pw + px0; PB.n1[q0] = ih(px0, py0) - 0.5; PB.n3[q0] = ih((px0 / 3) | 0, ((py0 / 3) | 0) + 7919) - 0.5; PB.cr[q0] = ih(px0 * 7 + 3, py0 * 13 + 1); }
    }
    if (!PB.N || PB.N !== N) PB = Object.assign(PB, { N: N, gr: new Float32Array(N * 3), wr: new Float32Array(N * 3), wa: new Float32Array(N), k: new Float32Array(N), cv: new Float32Array(N), vg: new Float32Array(N) });
    // 樹冠のまだら（真上から見た森林）：樹冠の大きさ約 5 m のボロノイ模様
    var crown = 2.6 * sub / sim.dx;
    if (PB.crown !== crown) {
      PB.crown = crown; PB.cn = new Float32Array(pw * ph);
      var cg = crown * 1.5;
      for (var qy = 0; qy < ph; qy++) for (var qx = 0; qx < pw; qx++) {
        var gx = qx / cg, gy = qy / cg, ix = Math.floor(gx), iy = Math.floor(gy), f1 = 9, tone = 0;
        for (var oy = -1; oy <= 1; oy++) for (var ox = -1; ox <= 1; ox++) {
          var cx_ = ix + ox, cy_ = iy + oy, px_ = cx_ + 0.15 + 0.7 * ih(cx_, cy_ + 911), py_ = cy_ + 0.15 + 0.7 * ih(cx_ + 577, cy_);
          var dd = (gx - px_) * (gx - px_) + (gy - py_) * (gy - py_); if (dd < f1) { f1 = dd; tone = ih(cx_ + 33, cy_ + 71); }
        }
        var rr_ = Math.sqrt(f1) * cg / crown;
        PB.cn[qy * pw + qx] = (rr_ < 1 ? 1.12 - 0.42 * rr_ * rr_ : 0.55) * (0.82 + 0.36 * tone);
      }
    }
    var z = sim.z, h = sim.h, dx = sim.dx, wc = waterCol(), col = [0, 0, 0], V = VARS[mode], gs = [0, 0];
    var ci = contourInfo(mode), real = mode === 'real', terr = !!(V && V.terr), i, j, c, tcol = terr ? terrColors(mode) : null;
    // セルごとの色・水・粗さ・等高線の値
    for (j = 0; j < ny; j++) for (i = 0; i < nx; i++) {
      c = j * nx + i;
      var zw = z[i > 0 ? c - 1 : c], ze = z[i < nx - 1 ? c + 1 : c], zs = z[j > 0 ? c - nx : c], zn = z[j < ny - 1 ? c + nx : c];
      var sh = 1 + 0.9 * ((zw - ze) / (2 * dx) - (sim.Ix ? sim.Ix[i] : -sim.trB)) + 0.9 * (zn - zs) / (2 * dx); sh = sh < 0.55 ? 0.55 : sh > 1.35 ? 1.35 : sh;
      var r, g, b, wa = 0;
      if (real || terr) {
        if (terr) { r = tcol[c * 3]; g = tcol[c * 3 + 1]; b = tcol[c * 3 + 2]; }
        else { ground(c, gcol); r = gcol[0] * sh; g = gcol[1] * sh; b = gcol[2] * sh; }
        var hh = h[c];
        if (hh > 0.01) {
          wa = Math.min(0.88, 1 - Math.exp(-hh * 3.2)); var t = 1 - Math.exp(-hh * 1.4);
          var wr = wc.sh[0] + (wc.dp[0] - wc.sh[0]) * t, wg = wc.sh[1] + (wc.dp[1] - wc.sh[1]) * t, wb = wc.sh[2] + (wc.dp[2] - wc.sh[2]) * t;
          var fr = Math.hypot(sim.qx[c], sim.qy[c]) / hh / Math.sqrt(G * hh);
          var fo = Math.max(0, Math.min(0.75, (fr - 1.15) * 0.8)) * (0.5 + hash(c + sim.nstep % 7));
          PB.wr[c * 3] = wr + (245 - wr) * fo; PB.wr[c * 3 + 1] = wg + (247 - wg) * fo; PB.wr[c * 3 + 2] = wb + (246 - wb) * fo;
        }
        PB.k[c] = !terr && sim.kind[c] !== 2 && gstat(c, gs) ? coarse(gs[0]) : 0;
        PB.vg[c] = !terr && ((sim.kind[c] === 1 && !(sim.bare && sim.bare[c]) && sim.z[c] - sim.z0[c] < 0.02) || (sim.bankC && sim.bankC[c] && Math.abs(sim.z[c] - sim.z0[c]) < 0.05)) ? 1 : 0;
      } else {
        PB.vg[c] = 0;
        var base = sim.kind[c] === 2 ? 150 : sim.kind[c] === 1 ? 108 : 128, v = value(mode, c);
        if (v === v) { V.ramp(norm(mode, v), col); var s2 = V.water ? 1 : (0.75 + 0.25 * sh); r = col[0] * s2; g = col[1] * s2; b = col[2] * s2; }
        else { r = base * sh; g = base * sh; b = (base + 2) * sh; }
        PB.k[c] = 0;
      }
      PB.gr[c * 3] = r; PB.gr[c * 3 + 1] = g; PB.gr[c * 3 + 2] = b; PB.wa[c] = wa;
      if (ci) PB.cv[c] = ci.rel ? relZ(c) : z[c];
    }
    // 画素ごとに双一次補間し、粒径に応じたきめを付ける
    var d = planImg.data, gr = PB.gr, wr3 = PB.wr, wa3 = PB.wa, kk = PB.k, cvv = PB.cv, n1 = PB.n1, n3 = PB.n3, crv = PB.cr, pv = PB.v, vg = PB.vg, cn = PB.cn;
    for (var py = 0; py < ph; py++) {
      var fy = (ph - py - 0.5) / sub - 0.5; if (fy < 0) fy = 0; if (fy > ny - 1) fy = ny - 1;
      var j0 = Math.min(ny - 2, fy | 0), ty = fy - j0;
      for (var px = 0; px < pw; px++) {
        var fx = (px + 0.5) / sub - 0.5; if (fx < 0) fx = 0; if (fx > nx - 1) fx = nx - 1;
        var i0 = Math.min(nx - 2, fx | 0), tx = fx - i0;
        var c00 = j0 * nx + i0, c10 = c00 + 1, c01 = c00 + nx, c11 = c01 + 1;
        var w00 = (1 - tx) * (1 - ty), w10 = tx * (1 - ty), w01 = (1 - tx) * ty, w11 = tx * ty;
        var R = gr[c00 * 3] * w00 + gr[c10 * 3] * w10 + gr[c01 * 3] * w01 + gr[c11 * 3] * w11;
        var Gc = gr[c00 * 3 + 1] * w00 + gr[c10 * 3 + 1] * w10 + gr[c01 * 3 + 1] * w01 + gr[c11 * 3 + 1] * w11;
        var B = gr[c00 * 3 + 2] * w00 + gr[c10 * 3 + 2] * w10 + gr[c01 * 3 + 2] * w01 + gr[c11 * 3 + 2] * w11;
        var q = py * pw + px, o = q * 4;
        if (real || terr) {
          var K = kk[c00] * w00 + kk[c10] * w10 + kk[c01] * w01 + kk[c11] * w11;
          // 細粒：細かく弱いざらつき　粗粒：大きなまだらと礫間の陰
          var nv = (1 - K) * n1[q] * 0.10 + K * (n3[q] * 0.55 + n1[q] * 0.22);
          var f = terr ? 1 : 1 + nv; if (!terr && K > 0.25 && crv[q] < 0.14 * K) f *= 0.62;
          var Vg = vg[c00] * w00 + vg[c10] * w10 + vg[c01] * w01 + vg[c11] * w11;
          if (Vg > 0.01) f = f * (1 - Vg) + cn[q] * Vg;
          R *= f; Gc *= f; B *= f;
          var A = wa3[c00] * w00 + wa3[c10] * w10 + wa3[c01] * w01 + wa3[c11] * w11;
          if (A > 0.004) {
            var sw = w00 * wa3[c00] + w10 * wa3[c10] + w01 * wa3[c01] + w11 * wa3[c11];
            var WR = (wr3[c00 * 3] * w00 * wa3[c00] + wr3[c10 * 3] * w10 * wa3[c10] + wr3[c01 * 3] * w01 * wa3[c01] + wr3[c11 * 3] * w11 * wa3[c11]) / sw;
            var WG = (wr3[c00 * 3 + 1] * w00 * wa3[c00] + wr3[c10 * 3 + 1] * w10 * wa3[c10] + wr3[c01 * 3 + 1] * w01 * wa3[c01] + wr3[c11 * 3 + 1] * w11 * wa3[c11]) / sw;
            var WB = (wr3[c00 * 3 + 2] * w00 * wa3[c00] + wr3[c10 * 3 + 2] * w10 * wa3[c10] + wr3[c01 * 3 + 2] * w01 * wa3[c01] + wr3[c11 * 3 + 2] * w11 * wa3[c11]) / sw;
            var Ae = A < 0.04 ? A * 0.5 : A;
            R += (WR - R) * Ae; Gc += (WG - Gc) * Ae; B += (WB - B) * Ae;
          }
        }
        d[o] = R; d[o + 1] = Gc; d[o + 2] = B; d[o + 3] = 255;
        if (ci) pv[q] = cvv[c00] * w00 + cvv[c10] * w10 + cvv[c01] * w01 + cvv[c11] * w11;
      }
    }
    if (ci) {   // 等高線：隣の画素と区間が変わる所を暗くする
      var inv = 1 / ci.ci;
      for (py = 0; py < ph - 1; py++) for (px = 0; px < pw - 1; px++) {
        q = py * pw + px; var a0 = Math.floor(pv[q] * inv);
        if (a0 !== Math.floor(pv[q + 1] * inv) || a0 !== Math.floor(pv[q + pw] * inv)) { o = q * 4; var dk = real ? 0.55 : 0.45; d[o] *= dk; d[o + 1] *= dk; d[o + 2] *= dk; }
      }
    }
    planOff.getContext('2d').putImageData(planImg, 0, 0);
    var ctx = cv.getContext('2d'); ctx.imageSmoothingEnabled = true; ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(planOff, 0, 0, cv.width, cv.height);
    var sx = cv.width / sim.Lx, sy = cv.height / sim.Ly;
    // 堰堤（河道に直角。濃い色が袖部、白が水通し）
    ctx.save(); ctx.lineCap = 'butt';
    sim.dams.forEach(function (D) {
      damSegs(D).forEach(function (g) {
        if (g.horiz) return;   // 横部材は平面図では描かない
        ctx.strokeStyle = g.steel ? (D.op && D.op.blocked ? 'rgba(200,60,40,.95)' : 'rgba(70,80,96,.95)') : g.spill ? 'rgba(236,236,228,.96)' : 'rgba(40,40,40,.9)'; ctx.lineWidth = Math.max(2 * dpr, dx * sx * (g.steel ? 0.5 : 0.85));
        ctx.beginPath();
        ctx.moveTo((D.cx + D.nx * g.s0) * sx, cv.height - (D.cy + D.ny * g.s0) * sy);
        ctx.lineTo((D.cx + D.nx * g.s1) * sx, cv.height - (D.cy + D.ny * g.s1) * sy);
        ctx.stroke();
      });
    });
    ctx.restore();
    if ($('p_core').checked) drawCores(ctx, sx, sy, cv.height, dpr, real);
    // 縦断の線と横断の位置（DEM）
    if (pathMode()) {
      var PSm = pathSamples(), Hc = cv.height;
      ctx.save(); ctx.lineJoin = 'round';
      DP.forEach(function (pth, q) {
        if (q === PI) return; ctx.strokeStyle = 'rgba(255,255,255,.55)'; ctx.lineWidth = 1.2 * dpr; ctx.beginPath();
        pth.pts.forEach(function (pt, k) { var X1 = pt[0] * sx, Y1 = Hc - pt[1] * sy; if (k) ctx.lineTo(X1, Y1); else ctx.moveTo(X1, Y1); }); ctx.stroke();
      });
      ctx.strokeStyle = 'rgba(255,214,0,.95)'; ctx.lineWidth = 2.4 * dpr; ctx.beginPath();
      DP[PI].pts.forEach(function (pt, k) { var X1 = pt[0] * sx, Y1 = Hc - pt[1] * sy; if (k) ctx.lineTo(X1, Y1); else ctx.moveTo(X1, Y1); }); ctx.stroke();
      if (xsS != null) {
        var kx = Math.max(0, Math.min(PSm.x.length - 1, Math.round(xsS / (PSm.L / PSm.x.length) - 0.5))), ux = -PSm.ty[kx], uy = PSm.tx[kx], Lh = 150;
        ctx.setLineDash([6 * dpr, 4 * dpr]); ctx.lineWidth = 1.5 * dpr; ctx.strokeStyle = 'rgba(255,255,255,.95)'; ctx.beginPath();
        ctx.moveTo((PSm.px[kx] + ux * Lh) * sx, Hc - (PSm.py[kx] + uy * Lh) * sy); ctx.lineTo((PSm.px[kx] - ux * Lh) * sx, Hc - (PSm.py[kx] - uy * Lh) * sy); ctx.stroke();
        ctx.strokeStyle = 'rgba(0,0,0,.6)'; ctx.lineDashOffset = 5 * dpr; ctx.stroke();
      }
      ctx.restore();
    }
    // 横断の位置
    else if (xsX != null) {
      ctx.save(); ctx.setLineDash([6 * dpr, 4 * dpr]); ctx.lineWidth = 1.5 * dpr; ctx.strokeStyle = 'rgba(255,255,255,.95)';
      ctx.beginPath(); ctx.moveTo(xsX * sx, 0); ctx.lineTo(xsX * sx, cv.height); ctx.stroke();
      ctx.strokeStyle = 'rgba(0,0,0,.6)'; ctx.lineDashOffset = 5 * dpr; ctx.stroke(); ctx.restore();
    }
    // 凡例
    var lg = $('legend'), ctxt = ci ? '　等高線 ' + fmt(ci.ci) + ' m 間隔' + (ci.rel ? '（初期勾配を除いた高さ）' : '') : '';
    var coreTxt = $('p_core').checked ? '　流心線：各横断の単位幅流量の極大を結んだ線（太さは流量）' : '';
    if (terr) lg.innerHTML = '<span>' + V.name + (mode === 'red' ? '（明るい：尾根　暗い：谷　赤い：急な斜面）' : mode === 'slp' ? '（濃いほど急：55° 以上で黒）' : '（北西から光を当てた陰影と標高の色）') + '　水色：水' + coreTxt + ctxt + '</span>';
    else if (real) lg.innerHTML = '<span>白い線：堰堤の水通し　濃い線：袖部　白波：射流（Fr &gt; 1.2）　まだらの粗さ：表層の粒径' + coreTxt + ctxt + '</span>';
    else {
      lg.innerHTML = '<span>' + V.name + '</span><span id="lg0"></span><canvas width="180" height="10"></canvas><span id="lg1"></span><span>' + coreTxt + ctxt + '</span>';
      var lc = lg.querySelector('canvas'), lx = lc.getContext('2d');
      for (i = 0; i < 180; i++) { V.ramp(i / 179, col); lx.fillStyle = 'rgb(' + (col[0] | 0) + ',' + (col[1] | 0) + ',' + (col[2] | 0) + ')'; lx.fillRect(i, 0, 1, 10); }
      $('lg0').textContent = fmt(vrange.lo) + (V.unit ? ' ' + V.unit : ''); $('lg1').textContent = fmt(vrange.hi) + (V.unit ? ' ' + V.unit : '') + (mode === 'dz' ? '（赤：堆積　青：侵食）' : mode === 'dm' ? '（対数）' : mode === 'bedrel' ? '（青：低い　茶：高い）' : '');
    }
  }
  // ------------------------------------------------------------ 流心線
  // 各横断で単位幅流量 q = h|u| の極大（流心）を求め、隣の横断の流心どうしを結ぶ。
  // 砂州で流れが分かれると、流心が2本以上に分かれて見える。
  function flowCores() {
    var nx = sim.nx, ny = sim.ny, h = sim.h, kind = sim.kind, qx = sim.qx, qy = sim.qy, dx = sim.dx, N = sim.N;
    var q = new Float64Array(N), qa = new Float64Array(N), cols = [], gmax = 0, i, j, c;
    for (c = 0; c < N; c++) q[c] = (kind[c] !== 2 && h[c] > 0.03) ? Math.hypot(qx[c], qy[c]) : 0;
    // 縦横に平滑化（流下方向は ±2 列、横断方向は ±2 行の重み付き平均）して格子規模の凹凸を除く
    var wx = [1, 2, 3, 2, 1], wy = [1, 2, 3, 2, 1];
    for (j = 0; j < ny; j++) for (i = 0; i < nx; i++) {
      c = j * nx + i; if (q[c] <= 0) { qa[c] = 0; continue; }
      var s = 0, w = 0;
      for (var a = -2; a <= 2; a++) for (var bq = -2; bq <= 2; bq++) {
        var ii = i + a, jj = j + bq; if (ii < 0 || ii >= nx || jj < 0 || jj >= ny) continue;
        var cc = jj * nx + ii; if (kind[cc] === 2) continue;
        var ww = wx[a + 2] * wy[bq + 2]; s += q[cc] * ww; w += ww;
      }
      qa[c] = w ? s / w : 0;
    }
    for (i = 0; i < nx; i++) {
      var cmax = 0, pk = [];
      for (j = 0; j < ny; j++) if (qa[j * nx + i] > cmax) cmax = qa[j * nx + i];
      if (cmax > 0) {
        for (j = 0; j < ny; j++) {
          var v = qa[j * nx + i], l = j > 0 ? qa[(j - 1) * nx + i] : 0, r = j < ny - 1 ? qa[(j + 1) * nx + i] : 0;
          if (v >= 0.3 * cmax && v >= l && v > r) pk.push(j);
        }
        for (var m = pk.length - 1; m > 0; m--) {   // 谷の浅い隣どうしは1つにまとめる
          var ja = pk[m - 1], jb = pk[m], vmin = 1e30; for (j = ja; j <= jb; j++) vmin = Math.min(vmin, qa[j * nx + i]);
          if (vmin > 0.75 * Math.min(qa[ja * nx + i], qa[jb * nx + i])) pk.splice(qa[ja * nx + i] >= qa[jb * nx + i] ? m : m - 1, 1);
        }
        if (cmax > gmax) gmax = cmax;
      }
      cols.push(pk.map(function (jj) {
        var v0 = qa[jj * nx + i], l = jj > 0 ? qa[(jj - 1) * nx + i] : v0, r = jj < ny - 1 ? qa[(jj + 1) * nx + i] : v0, den = l - 2 * v0 + r, off = den < 0 ? 0.5 * (l - r) / den : 0;
        return { i: i, y: (jj + 0.5 + Math.max(-0.5, Math.min(0.5, off))) * dx, q: v0 };
      }));
    }
    // 流心をつないで線にする。1本の線に2つの流心が続くときは分岐、2本が1つに入るときは合流とする
    var lim = Math.max(5 * dx, 10), gap = 4, tracks = [], open = [];
    for (i = 0; i < nx; i++) {
      var P = cols[i], claimed = {};
      open = open.filter(function (t) { return i - t.pts[t.pts.length - 1].i <= gap; });
      var order = P.map(function (p, k) { return k; }).sort(function (a2, b2) { return P[b2].q - P[a2].q; });
      order.forEach(function (k) {
        var p = P[k], best = null, bd = lim;
        open.forEach(function (t) { var e = t.pts[t.pts.length - 1]; if (e.i >= i) return; var d = Math.abs(e.y - p.y); if (d < bd) { bd = d; best = t; } });
        if (best && !claimed[best.id]) { best.pts.push(p); claimed[best.id] = 1; }
        else if (best) { var nt = { id: tracks.length, pts: [best.pts[best.pts.length - 2] || best.pts[best.pts.length - 1], p], fixed: 1 }; tracks.push(nt); open.push(nt); claimed[nt.id] = 1; }
        else { var t0 = { id: tracks.length, pts: [p], fixed: 0 }; tracks.push(t0); open.push(t0); claimed[t0.id] = 1; }
      });
      // 続かなかった線は、近くの流心へつないで合流として終える
      open.forEach(function (t) {
        if (claimed[t.id]) return; var e = t.pts[t.pts.length - 1]; if (e.i !== i - 1 || t.joined) return;
        var best = null, bd = lim; P.forEach(function (p) { var d = Math.abs(e.y - p.y); if (d < bd) { bd = d; best = p; } });
        if (best) { t.pts.push(best); t.joined = 1; }
      });
    }
    // 各線の位置を流下方向に平滑化する（端は固定）
    tracks.forEach(function (t) {
      var n = t.pts.length; if (n < 5) return;
      var ys = t.pts.map(function (p) { return p.y; }), out = ys.slice();
      for (var k = 2; k < n - 2; k++) out[k] = (ys[k - 2] + 2 * ys[k - 1] + 3 * ys[k] + 2 * ys[k + 1] + ys[k + 2]) / 9;
      t.ys = out;
    });
    return { tracks: tracks.filter(function (t) { var qm = 0; t.pts.forEach(function (p) { qm += p.q; }); return t.pts.length >= 6 && qm / t.pts.length >= 0.15 * gmax; }), gmax: gmax };
  }
  function drawCores(ctx, sx, sy, H, dpr, real) {
    var F = flowCores(), dx = sim.dx;
    ctx.save(); ctx.lineCap = 'round'; ctx.lineJoin = 'round';
    for (var pass = 0; pass < 2; pass++) {
      F.tracks.forEach(function (t) {
        var ys = t.ys || t.pts.map(function (p) { return p.y; });
        for (var k = 1; k < t.pts.length; k++) {
          var a = t.pts[k - 1], b = t.pts[k], w = (1.1 + 3.4 * (a.q + b.q) / (2 * F.gmax || 1)) * dpr;
          ctx.strokeStyle = pass === 0 ? (real ? 'rgba(10,30,40,.45)' : 'rgba(255,255,255,.75)') : (real ? 'rgba(255,250,235,.96)' : 'rgba(20,28,26,.9)');
          ctx.lineWidth = pass === 0 ? w + 2 * dpr : w;
          ctx.beginPath(); ctx.moveTo((a.i + 0.5) * dx * sx, H - ys[k - 1] * sy); ctx.lineTo((b.i + 0.5) * dx * sx, H - ys[k] * sy); ctx.stroke();
        }
      });
    }
    ctx.restore();
  }
  // ------------------------------------------------------------ 横断図
  // 谷筋に直角な横断（下流を向いて左が左岸）。幅は流域の外に出るか ±400 m まで
  function drawXsecPath() {
    var S = pathSamples(), n = S.x.length, dx = sim.dx; if (xsS == null || xsS > S.L) xsS = S.L * 0.55;
    var k = Math.max(0, Math.min(n - 1, Math.round(xsS / (S.L / n) - 0.5))), nxv = -S.ty[k], nyv = S.tx[k], st = dx / 2, ext = [0, 0];
    [1, -1].forEach(function (sg, q) { var u = 0; while (u < 400) { var c = cellXY(S.px[k] + nxv * sg * (u + st), S.py[k] + nyv * sg * (u + st)); if (!inBasin(c)) break; u += st; } ext[q] = u; });
    var yy = [], zb = [], z0 = [], ws = [], dmA = [], gs = [0, 0], isDam = false;
    for (var u = ext[0]; u >= -ext[1] - 1e-9; u -= st) {
      var c = cellXY(S.px[k] + nxv * u, S.py[k] + nyv * u); if (c < 0) continue;
      yy.push(ext[0] - u); zb.push(sim.z[c]); z0.push(sim.z0[c]); ws.push(sim.h[c] > 0.01 ? sim.z[c] + sim.h[c] : NaN);
      dmA.push(sim.kind[c] !== 2 && gstat(c, gs) ? gs[0] : NaN); if (sim.kind[c] === 2) isDam = true;
    }
    var nm = DP[PI].name.replace(/（.*$/, '');
    $('xs_title').textContent = '横断（' + nm + '、上流端から ' + S.x[k].toFixed(0) + ' m' + (isDam ? '、堰堤の位置' : '') + '）';
    $('xs_val').textContent = S.x[k].toFixed(0) + ' m';
    var dlo = sim.d[0], dhi = sim.dmax, colr = [0, 0, 0], W = ext[0] + ext[1];
    chart($('xsec'), {
      series: [
        { x: yy, y: zb, color: css('--sed'), fill: css('--sed') + '55', label: '河床', w: 2 },
        { x: yy, y: z0, color: css('--ink3'), label: '初期河床', w: 1.2 },
        { x: yy, y: ws, color: css('--s1'), label: '水面', w: 1.6 }
      ], x0: 0, x1: Math.max(W, dx), zoom: 'xs', xlabel: '左岸側の端からの距離 [m]（谷筋は ' + ext[0].toFixed(0) + ' m）', xname: '左から', xunit: ' m', height: 220,
      over: function (ctx, X, Y) {
        for (var q = 0; q < yy.length; q++) {
          if (dmA[q] !== dmA[q]) continue;
          seqOr(Math.log(dmA[q] / dlo) / Math.log(dhi / dlo || 2), colr);
          ctx.fillStyle = 'rgb(' + (colr[0] | 0) + ',' + (colr[1] | 0) + ',' + (colr[2] | 0) + ')';
          var xa = X(yy[q] - st / 2), xb = X(yy[q] + st / 2); ctx.fillRect(xa, Y(zb[q]) + 1, xb - xa + 0.5, 5);
        }
        ctx.save(); ctx.setLineDash([3, 3]); ctx.strokeStyle = css('--ink3'); ctx.beginPath(); ctx.moveTo(X(ext[0]), 0); ctx.lineTo(X(ext[0]), 600); ctx.stroke(); ctx.restore();
      }
    });
    $('xs_key').innerHTML = keyHTML([{ color: css('--sed'), label: '河床' }, { color: css('--ink3'), label: '初期河床' }, { color: css('--s1'), label: '水面' }]) +
      '<span><i class="b" style="background:linear-gradient(90deg,#fbe9df,#eb8457,#7a3013)"></i>表層の平均粒径 ' + fmt(dlo * 1000) + '〜' + fmt(dhi * 1000) + ' mm（対数）</span><span>点線：谷筋（縦断の線）の位置</span>';
  }
  function drawXsec() {
    if (!sim) return;
    if (pathMode()) { drawXsecPath(); return; }
    if (xsX == null || xsX > sim.Lx) xsX = sim.Lx / 2;
    var nx = sim.nx, ny = sim.ny, dx = sim.dx, i = Math.min(nx - 1, Math.max(0, Math.floor(xsX / dx))), Ly = sim.Ly;
    var yy = new Float64Array(ny), zb = new Float64Array(ny), z0 = new Float64Array(ny), ws = new Float64Array(ny), dmA = new Float64Array(ny), gs = [0, 0], isDam = false;
    for (var k = 0; k < ny; k++) {
      var j = ny - 1 - k, c = j * nx + i;
      yy[k] = Ly - (j + 0.5) * dx; zb[k] = sim.z[c]; z0[k] = sim.z0[c];
      ws[k] = sim.h[c] > 0.01 ? sim.z[c] + sim.h[c] : NaN;
      dmA[k] = sim.kind[c] !== 2 && gstat(c, gs) ? gs[0] : NaN;
      if (sim.kind[c] === 2) isDam = true;
    }
    $('xs_title').textContent = '横断（x = ' + ((i + 0.5) * dx).toFixed(0) + ' m' + (isDam ? '、堰堤の位置' : '') + '）';
    $('xs_val').textContent = ((i + 0.5) * dx).toFixed(0) + ' m';
    var dlo = sim.d[0], dhi = sim.dmax, colr = [0, 0, 0];
    chart($('xsec'), {
      series: [
        { x: yy, y: zb, color: css('--sed'), fill: css('--sed') + '55', label: '河床', w: 2 },
        { x: yy, y: z0, color: css('--ink3'), label: '初期河床', w: 1.2 },
        { x: yy, y: ws, color: css('--s1'), label: '水面', w: 1.6 }
      ], x0: 0, x1: Ly, zoom: 'xs', xlabel: '左岸からの距離 [m]', xname: '左岸から', xunit: ' m', height: 220,
      over: function (ctx, X, Y) {   // 河床表層の平均粒径を河床線の下に色で示す
        for (var k = 0; k < ny; k++) {
          if (dmA[k] !== dmA[k]) continue;
          seqOr(Math.log(dmA[k] / dlo) / Math.log(dhi / dlo || 2), colr);
          ctx.fillStyle = 'rgb(' + (colr[0] | 0) + ',' + (colr[1] | 0) + ',' + (colr[2] | 0) + ')';
          var xa = X(yy[k] - dx / 2), xb = X(yy[k] + dx / 2); ctx.fillRect(xa, Y(zb[k]) + 1, xb - xa + 0.5, 5);
        }
      }
    });
    $('xs_key').innerHTML = keyHTML([{ color: css('--sed'), label: '河床' }, { color: css('--ink3'), label: '初期河床' }, { color: css('--s1'), label: '水面' }]) +
      '<span><i class="b" style="background:linear-gradient(90deg,#fbe9df,#eb8457,#7a3013)"></i>表層の平均粒径 ' + fmt(dlo * 1000) + '〜' + fmt(dhi * 1000) + ' mm（対数）</span>';
  }
  function fmt(v) { var a = Math.abs(v); return a >= 100 ? v.toFixed(0) : a >= 10 ? v.toFixed(1) : a >= 1 ? v.toFixed(2) : v.toFixed(3); }
  (function () {
    var cv = $('plan'), tip = $('plan_tip');
    cv.addEventListener('pointermove', function (e) {
      if (!sim) return; var r = cv.getBoundingClientRect(), x = (e.clientX - r.left) / r.width * sim.Lx, y = (1 - (e.clientY - r.top) / r.height) * sim.Ly;
      var i = Math.min(sim.nx - 1, Math.max(0, Math.floor(x / sim.dx))), j = Math.min(sim.ny - 1, Math.max(0, Math.floor(y / sim.dx))), c = j * sim.nx + i, h = sim.h[c];
      var dm = sim.surfaceDm(c), gs = [0, 0], d84 = gstat(c, gs) ? gs[1] : NaN;
      tip.innerHTML = 'x=' + x.toFixed(0) + ' m　河床位 ' + sim.z[c].toFixed(2) + ' m　Δz ' + (sim.z[c] - sim.z0[c]).toFixed(2) + ' m<br>水深 ' + h.toFixed(2) + ' m　流速 ' + (h > 0.01 ? (Math.hypot(sim.qx[c], sim.qy[c]) / h).toFixed(2) : '–') + ' m/s<br>表層 dm ' + (dm === dm ? (dm * 1000).toFixed(0) + ' mm' : '–') + '　d84 ' + (d84 === d84 ? (d84 * 1000).toFixed(0) + ' mm' : '–');
      tip.hidden = false; var tx = e.clientX - r.left + 12; if (tx > r.width - 230) tx = r.width - 230; tip.style.left = Math.max(0, tx) + 'px'; tip.style.top = (e.clientY - r.top + 14) + 'px';
    });
    cv.addEventListener('pointerleave', function () { tip.hidden = true; });
    cv.addEventListener('click', function (e) { if (!sim) return; var r = cv.getBoundingClientRect(); if (pathMode()) { setXsFromPlan((e.clientX - r.left) / r.width * sim.Lx, (1 - (e.clientY - r.top) / r.height) * sim.Ly); drawXsec(); drawPlan(); drawProfile(); return; } xsX = Math.max(0, Math.min(sim.Lx - 1e-6, (e.clientX - r.left) / r.width * sim.Lx)); $('xs_x').value = xsX; drawXsec(); drawPlan(); drawProfile(); });
  })();

  // ------------------------------------------------------------ 折れ線グラフ
  // 横軸の拡大範囲（キーごとに共有）
  var ZOOM = {};
  function zoomRange(key, f0, f1) { var z = ZOOM[key]; if (!z || Math.abs(z.f0 - f0) > 1e-9 || Math.abs(z.f1 - f1) > 1e-9) z = ZOOM[key] = { f0: f0, f1: f1, x0: f0, x1: f1 }; return z; }
  function zoomBy(key, f, xc) {
    var z = ZOOM[key]; if (!z) return; var w = z.x1 - z.x0, full = z.f1 - z.f0, nw = Math.max(full / 300, Math.min(full, w * f));
    if (xc == null) xc = (z.x0 + z.x1) / 2;
    var a = xc - (xc - z.x0) * nw / w; a = Math.max(z.f0, Math.min(z.f1 - nw, a)); z.x0 = a; z.x1 = a + nw;
  }
  function zoomPan(key, d) { var z = ZOOM[key]; if (!z) return; var w = z.x1 - z.x0, a = Math.max(z.f0, Math.min(z.f1 - w, z.x0 + d)); z.x0 = a; z.x1 = a + w; }
  function chart(cv, o) {
    var W = cv.clientWidth, Hh = o.height || 220, dpr = Math.min(2, window.devicePixelRatio || 1);
    if (cv.width !== Math.round(W * dpr) || cv.height !== Math.round(Hh * dpr)) { cv.width = Math.round(W * dpr); cv.height = Math.round(Hh * dpr); cv.style.height = Hh + 'px'; }
    var ctx = cv.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, W, Hh);
    var ml = 46, mr = 10, mt = 8, mb = 26, pw = W - ml - mr, ph = Hh - mt - mb;
    var x0 = o.x0, x1 = o.x1, y0 = o.y0, y1 = o.y1;
    if (o.zoom) { var zr = zoomRange(o.zoom, o.x0, o.x1); x0 = zr.x0; x1 = zr.x1; }
    if (y0 == null) {
      y0 = 1e30; y1 = -1e30;
      o.series.forEach(function (s) { for (var i = 0; i < s.x.length; i++) { var v = s.y[i]; if (v === v && s.x[i] >= x0 && s.x[i] <= x1) { if (v < y0) y0 = v; if (v > y1) y1 = v; } } });
      if (o.zero && y0 > 0) y0 = 0; if (y0 > y1) { y0 = 0; y1 = 1; } if (y1 - y0 < 1e-9) y1 = y0 + 1;
      var pad = (y1 - y0) * 0.06; y1 += pad; if (!(o.zero && y0 === 0)) y0 -= pad;
    }
    var X = function (x) { return ml + (x - x0) / (x1 - x0) * pw; }, Y = function (y) { return mt + (1 - (y - y0) / (y1 - y0)) * ph; };
    var ink3 = css('--ink3'), line = css('--line');
    ctx.font = '11px ' + css('--mono'); ctx.fillStyle = ink3; ctx.strokeStyle = line; ctx.lineWidth = 1;
    var ty = nice((y1 - y0) / 4), tx = nice((x1 - x0) / 7), v;
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    for (v = Math.ceil(y0 / ty) * ty; v <= y1 + 1e-9; v += ty) { ctx.beginPath(); ctx.moveTo(ml, Y(v) + .5); ctx.lineTo(ml + pw, Y(v) + .5); ctx.stroke(); ctx.fillText(+v.toFixed(6) + '', ml - 6, Y(v)); }
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    for (v = Math.ceil(x0 / tx) * tx; v <= x1 + 1e-9; v += tx) ctx.fillText(+v.toFixed(6) + '', X(v), mt + ph + 5);
    ctx.textAlign = 'right'; ctx.fillText(o.xlabel || '', ml + pw, mt + ph + 14);
    ctx.save(); ctx.beginPath(); ctx.rect(ml, mt, pw, ph); ctx.clip();
    if (o.under) o.under(ctx, X, Y, y0);
    o.series.forEach(function (s) {
      if (s.fill) {
        ctx.fillStyle = s.fill; ctx.beginPath(); var st = false;
        for (var i = 0; i < s.x.length; i++) { if (s.y[i] !== s.y[i]) continue; if (!st) { ctx.moveTo(X(s.x[i]), Y(y0)); st = true; } ctx.lineTo(X(s.x[i]), Y(s.y[i])); }
        if (st) { ctx.lineTo(X(s.x[s.x.length - 1]), Y(y0)); ctx.fill(); }
      }
      ctx.strokeStyle = s.color; ctx.lineWidth = s.w || 2; ctx.setLineDash(s.dash ? [5, 4] : []); ctx.lineJoin = 'round'; ctx.beginPath(); var pen = false;
      for (i = 0; i < s.x.length; i++) { if (s.y[i] !== s.y[i]) { pen = false; continue; } if (pen) ctx.lineTo(X(s.x[i]), Y(s.y[i])); else { ctx.moveTo(X(s.x[i]), Y(s.y[i])); pen = true; } }
      ctx.stroke(); ctx.setLineDash([]);
    });
    if (o.over) o.over(ctx, X, Y, y0);
    if (cv._hx != null) { ctx.strokeStyle = ink3; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(X(cv._hx) + .5, mt); ctx.lineTo(X(cv._hx) + .5, mt + ph); ctx.stroke(); }
    ctx.restore();
    cv._o = o; cv._g = { ml: ml, pw: pw, x0: x0, x1: x1 };
    if (!cv._bound) {
      cv._bound = true; var tip = cv.parentNode.querySelector('.tip'), drag = null;
      cv.addEventListener('wheel', function (e) {
        var g = cv._g, oo = cv._o; if (!g || !oo.zoom) return; e.preventDefault();
        var r = cv.getBoundingClientRect(), xc = g.x0 + (e.clientX - r.left - g.ml) / g.pw * (g.x1 - g.x0);
        var dy = e.deltaY * (e.deltaMode === 1 ? 33 : e.deltaMode === 2 ? 400 : 1);
        zoomBy(oo.zoom, Math.exp(Math.max(-0.5, Math.min(0.5, dy * 0.0015))), xc); redrawCharts();
      }, { passive: false });
      cv.addEventListener('pointerdown', function (e) { if (cv._o && cv._o.zoom) { drag = e.clientX; cv.setPointerCapture(e.pointerId); cv.style.cursor = 'grabbing'; } });
      cv.addEventListener('pointerup', function () { drag = null; cv.style.cursor = ''; });
      cv.addEventListener('dblclick', function () { var oo = cv._o; if (oo && oo.zoom && ZOOM[oo.zoom]) { var z = ZOOM[oo.zoom]; z.x0 = z.f0; z.x1 = z.f1; redrawCharts(); } });
      cv.addEventListener('pointermove', function (e) {
        var g = cv._g, oo = cv._o; if (!g) return;
        if (drag != null && oo.zoom) { zoomPan(oo.zoom, -(e.clientX - drag) / g.pw * (g.x1 - g.x0)); drag = e.clientX; } var r = cv.getBoundingClientRect(), x = g.x0 + (e.clientX - r.left - g.ml) / g.pw * (g.x1 - g.x0);
        if (x < g.x0 || x > g.x1) { tip.hidden = true; cv._hx = null; return; }
        cv._hx = x; var html = (oo.xname || 'x') + ' ' + fmt(x) + (oo.xunit || '');
        oo.series.forEach(function (s) { if (!s.label || !s.x.length) return; var bi = 0, bd = 1e30; for (var i = 0; i < s.x.length; i++) { var dd = Math.abs(s.x[i] - x); if (dd < bd) { bd = dd; bi = i; } } var v = s.y[bi]; if (v === v) html += '<br><span style="color:' + s.color + '">■</span> ' + s.label + ' ' + fmt(v); });
        tip.innerHTML = html; tip.hidden = false; var tx = e.clientX - r.left + 12; if (tx > r.width - 190) tx = e.clientX - r.left - 190; tip.style.left = Math.max(0, tx) + 'px'; tip.style.top = '8px'; redrawCharts();
      });
      cv.addEventListener('pointerleave', function () { tip.hidden = true; cv._hx = null; redrawCharts(); });
    }
  }
  function keyHTML(items) { return items.map(function (k) { return '<span><i class="' + (k.dash ? 'd' : k.box ? 'b' : '') + '" style="' + (k.box ? 'background:' : 'border-color:') + k.color + '"></i>' + k.label + '</span>'; }).join(''); }

  var prof = null;
  // ------------------------------------------------------------ DEM：谷に沿った縦断と、それに直角な横断
  var DP = null, PI = 0, xsS = null, PS = null;   // 線の一覧、選んだ線、横断の位置（上流端からの距離）、標本の控え
  function pathMode() { return !!(DP && DP.length && sim && sim.demMode); }
  // 選んだ線を dx おきの点に分ける。点ごとに位置と下流向きの単位ベクトル
  function pathSamples() {
    if (PS && PS.pi === PI && PS.sim === sim) return PS;
    var pts = DP[PI].pts, cum = [0], k;
    for (k = 1; k < pts.length; k++) cum.push(cum[k - 1] + Math.hypot(pts[k][0] - pts[k - 1][0], pts[k][1] - pts[k - 1][1]));
    var L = cum[cum.length - 1], dx = sim.dx, n = Math.max(2, Math.floor(L / dx)), S = { pi: PI, sim: sim, L: L, x: new Float64Array(n), px: [], py: [], tx: [], ty: [] }, q = 0;
    function at(s) { while (q < cum.length - 2 && cum[q + 1] < s) q++; while (q > 0 && cum[q] > s) q--; var a = pts[q], b = pts[q + 1], u = (s - cum[q]) / Math.max(1e-9, cum[q + 1] - cum[q]); return [a[0] + u * (b[0] - a[0]), a[1] + u * (b[1] - a[1])]; }
    for (k = 0; k < n; k++) {
      var sk = (k + 0.5) * L / n, p = at(sk), pa = at(Math.max(0, sk - 1.5 * dx)), pb = at(Math.min(L, sk + 1.5 * dx));
      var tx = pb[0] - pa[0], ty = pb[1] - pa[1], tl = Math.hypot(tx, ty) || 1;
      S.x[k] = sk; S.px.push(p[0]); S.py.push(p[1]); S.tx.push(tx / tl); S.ty.push(ty / tl);
    }
    return (PS = S);
  }
  function cellXY(x, y) { var i = Math.floor(x / sim.dx), j = Math.floor(y / sim.dx); return i < 0 || j < 0 || i >= sim.nx || j >= sim.ny ? -1 : j * sim.nx + i; }
  function inBasin(c) { return c >= 0 && (!sim.basin || sim.basin[c]); }
  // 縦断：各点で、谷筋に直角な線の上の谷底の帯（±谷筋の幅の半分）から最深・平均の河床、水面、水面幅、表層の粒径、流砂量
  function pathProfile() {
    var S = pathSamples(), n = S.x.length, dx = sim.dx, K = sim.K, nx = sim.nx, z = sim.z, z0 = sim.z0, h = sim.h, kind = sim.kind;
    var hw = Math.max(dx, ((sim.P.dem && sim.P.dem.chW) || 6) / 2 + 0.5 * dx), ww = Math.max(3 * hw, 30);
    var R = { x: S.x, zmin: new Float64Array(n), zmean: new Float64Array(n), z0min: new Float64Array(n), z0mean: new Float64Array(n), ws: new Float64Array(n), dm: new Float64Array(n), wet: new Float64Array(n), qs: new Float64Array(n), hmax: new Float64Array(n), isDam: new Uint8Array(n) };
    for (var k = 0; k < n; k++) {
      var nxv = -S.ty[k], nyv = S.tx[k], mn = 1e30, m0 = 1e30, sm = 0, s0 = 0, cnt = 0, wsum = 0, wcnt = 0, ds = 0, dc = 0, hm = 0, wet = 0, q = 0, seen = {}, dam = false, gs = [0, 0];
      for (var u = -ww; u <= ww + 1e-9; u += dx / 2) {
        var c = cellXY(S.px[k] + nxv * u, S.py[k] + nyv * u); if (!inBasin(c) || seen[c]) continue; seen[c] = 1;
        if (kind[c] === 2) { if (Math.abs(u) <= hw) dam = true; continue; }
        if (h[c] > 0.03) wet += dx;
        if (Math.abs(u) > hw) continue;
        if (z[c] < mn) mn = z[c]; if (z0[c] < m0) m0 = z0[c]; sm += z[c]; s0 += z0[c]; cnt++;
        if (h[c] > 0.03) { wsum += (z[c] + h[c]) * h[c]; wcnt += h[c]; if (h[c] > hm) hm = h[c]; if (gstat(c, gs)) { ds += gs[0]; dc++; } }
        var i = c % nx, j = (c - i) / nx, fx = 0, fy = 0;
        for (var kk = 0; kk < K; kk++) { fx += 0.5 * (sim.SFx[(j * (nx + 1) + i) * K + kk] + sim.SFx[(j * (nx + 1) + i + 1) * K + kk]); fy += 0.5 * (sim.SFy[(j * nx + i) * K + kk] + sim.SFy[((j + 1) * nx + i) * K + kk]); }
        q += (fx * S.tx[k] + fy * S.ty[k]) * dx;
      }
      if (cnt) { R.zmin[k] = mn; R.z0min[k] = m0; R.zmean[k] = sm / cnt; R.z0mean[k] = s0 / cnt; } else { R.zmin[k] = R.z0min[k] = R.zmean[k] = R.z0mean[k] = NaN; }
      R.ws[k] = wcnt > 0 ? wsum / wcnt : NaN; R.dm[k] = dc ? ds / dc : NaN; R.wet[k] = wet; R.qs[k] = q; R.hmax[k] = hm; R.isDam[k] = dam ? 1 : 0;
    }
    // 初期河床（最深）の直線の当てはめ：「初期勾配を差し引く」に使う
    var sx = 0, sy = 0, sxx = 0, sxy = 0, m = 0;
    for (k = 0; k < n; k++) if (R.z0min[k] === R.z0min[k]) { sx += S.x[k]; sy += R.z0min[k]; sxx += S.x[k] * S.x[k]; sxy += S.x[k] * R.z0min[k]; m++; }
    var bb = m > 1 ? (m * sxy - sx * sy) / Math.max(1e-9, m * sxx - sx * sx) : 0; R.trB = bb; R.trA = m ? (sy - bb * sx) / m : 0;
    // 線の近くの堰堤：線上の最も近い点
    R.dams = [];
    sim.dams.forEach(function (D) {
      var best = -1, bd = 1e30; for (var k2 = 0; k2 < n; k2++) { var d = Math.hypot(S.px[k2] - D.cx, S.py[k2] - D.cy); if (d < bd) { bd = d; best = k2; } }
      if (best >= 0 && bd <= Math.max(2 * dx, D.ext)) R.dams.push({ s: S.x[best], crest: D.crest });
    });
    return R;
  }
  function pathOptions() {
    var sel = $('l_path'), on = pathMode();
    $('l_pathw').hidden = !on; if (!on) return;
    sel.innerHTML = DP.map(function (p, q) { return '<option value="' + q + '"' + (q === PI ? ' selected' : '') + '>' + p.name + '</option>'; }).join('');
  }
  function setXsFromPlan(x, y) {   // 平面図のクリック：いちばん近い線の上の点へ（いまの線を優先）
    var best = -1, bs = 0, bd = 1e30, keep = PI;
    for (var q = 0; q < DP.length; q++) {
      PI = q; var S = pathSamples();
      for (var k = 0; k < S.x.length; k++) { var d = Math.hypot(S.px[k] - x, S.py[k] - y) - (q === keep ? 2 * sim.dx : 0); if (d < bd) { bd = d; best = q; bs = S.x[k]; } }
    }
    PI = best >= 0 ? best : keep; xsS = bs; PS = null; pathOptions(); var S2 = pathSamples(); $('xs_x').max = S2.L; $('xs_x').value = xsS;
  }
  function drawProfile() {
    var PM = pathMode();
    prof = PM ? pathProfile() : sim.profile();
    var det = $('l_det').checked, useMin = $('l_bed').value === 'min', nx = PM ? prof.x.length : sim.nx, I = sim.I, L = PM ? pathSamples().L : sim.Lx;
    var trend = PM ? function (x) { return prof.trA + prof.trB * x; } : trendX;
    var tr = function (zarr, i) { return det ? zarr[i] - trend(prof.x[i]) : zarr[i]; };
    var zb = new Float64Array(nx), z0 = new Float64Array(nx), ws = new Float64Array(nx), i;
    var isDam = PM ? prof.isDam : new Uint8Array(nx); if (!PM) sim.dams.forEach(function (D) { isDam[D.i] = 1; });
    var damMarks = PM ? prof.dams.map(function (d) { return { x: d.s, crest: d.crest }; }) : sim.dams.map(function (D) { return { x: D.x, crest: D.crest }; });
    var xsPos = PM ? xsS : xsX, xlab = PM ? '谷に沿った距離（' + DP[PI].name.replace(/（.*$/, '') + 'の上流端から）[m]' : '上流端からの距離 [m]';
    for (i = 0; i < nx; i++) {
      zb[i] = isDam[i] ? NaN : tr(useMin ? prof.zmin : prof.zmean, i); z0[i] = isDam[i] ? NaN : tr(useMin ? prof.z0min : prof.z0mean, i);
      ws[i] = prof.ws[i] === prof.ws[i] ? (det ? prof.ws[i] - trend(prof.x[i]) : prof.ws[i]) : NaN;
      if (useMin && ws[i] === ws[i]) ws[i] = (det ? prof.zmin[i] - trend(prof.x[i]) : prof.zmin[i]) + prof.hmax[i];
    }
    var series = [
      { x: prof.x, y: zb, color: css('--sed'), fill: css('--sed') + '55', label: '河床', w: 2 },
      { x: prof.x, y: z0, color: css('--ink3'), label: '初期河床', w: 1.2 },
      { x: prof.x, y: ws, color: css('--s1'), label: '水面', w: 1.6 }
    ];
    if (!PM && ref && ref.nx === nx && ref.dx === sim.dx) {
      var zr = new Float64Array(nx); for (i = 0; i < nx; i++) zr[i] = ref.isDam[i] ? NaN : (det ? (useMin ? ref.zmin[i] : ref.zmean[i]) - trend(prof.x[i]) : (useMin ? ref.zmin[i] : ref.zmean[i]));
      series.push({ x: prof.x, y: zr, color: css('--sref'), dash: true, label: '比較', w: 1.6 });
    }
    chart($('prof'), {
      series: series, x0: 0, x1: L, xlabel: xlab, xname: PM ? 's' : 'x', xunit: ' m', height: 260, zoom: 'prof',
      over: function (ctx, X, Y, y0) {
        ctx.fillStyle = css('--ink2');
        if (xsPos != null) { ctx.save(); ctx.setLineDash([5, 4]); ctx.strokeStyle = css('--ink2'); ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(X(xsPos), 0); ctx.lineTo(X(xsPos), 600); ctx.stroke(); ctx.restore(); }
        damMarks.forEach(function (D) { var top = det ? D.crest - trend(D.x) : D.crest, w = Math.max(3, X(sim.dx) - X(0)); ctx.fillRect(X(D.x) - w / 2, Y(top), w, Y(y0) - Y(top)); });
      }
    });
    var key = [{ color: css('--sed'), label: '河床（' + (useMin ? '最深' : PM ? '谷底の帯の平均' : '横断平均') + '）' }, { color: css('--ink3'), label: '初期河床' }, { color: css('--s1'), label: '水面' }, { color: css('--ink2'), box: true, label: '堰堤' }];
    if (series.length > 3) key.push({ color: css('--sref'), dash: true, label: '比較（' + ref.label + '）' });
    $('prof_key').innerHTML = keyHTML(key);
    // 縦断分布
    var lv = $('l_var').value, y = new Float64Array(nx), lab = { dm: '平均粒径 mm', qs: '流砂量 m³/s', wet: '水面幅 m', dz: 'Δz m' }[lv];
    for (i = 0; i < nx; i++) y[i] = lv === 'dm' ? prof.dm[i] * 1000 : lv === 'qs' ? prof.qs[i] : lv === 'wet' ? prof.wet[i] : (isDam[i] ? NaN : prof.zmean[i] - prof.z0mean[i]);
    var ss = [{ x: prof.x, y: y, color: lv === 'dz' ? css('--sed') : css('--s3'), label: lab }];
    if (!PM && ref && ref.nx === nx && ref.dx === sim.dx && ref[lv]) ss.push({ x: prof.x, y: ref[lv], color: css('--sref'), dash: true, label: '比較', w: 1.5 });
    chart($('strip'), { series: ss, x0: 0, x1: L, xlabel: xlab, xname: PM ? 's' : 'x', xunit: ' m', height: 150, zoom: 'prof', zero: lv !== 'dz' && lv !== 'dm',
      under: function (ctx, X, Y) { ctx.strokeStyle = css('--line'); ctx.lineWidth = 1; ctx.setLineDash([2, 3]); damMarks.forEach(function (D) { ctx.beginPath(); ctx.moveTo(X(D.x), 0); ctx.lineTo(X(D.x), 400); ctx.stroke(); }); ctx.setLineDash([]); } });
  }
  // ------------------------------------------------------------ 横断位置の流量・流砂量（いまの横断の線を通る量。計算の進みに合わせて記録）
  var XH = null;
  function xsKey() { return sim ? (pathMode() ? 'p' + PI + ':' + xsS : 'x' + xsX) : ''; }
  function xsFlux() {
    var dx = sim.dx, nx = sim.nx, K = sim.K, Q = 0, Qs = 0, st = dx / 2, px, py, tx, ty, ux, uy, u0, u1;
    if (pathMode()) {
      var S = pathSamples(), n = S.x.length, k = Math.max(0, Math.min(n - 1, Math.round(xsS / (S.L / n) - 0.5)));
      px = S.px[k]; py = S.py[k]; tx = S.tx[k]; ty = S.ty[k]; ux = -ty; uy = tx; u0 = -400; u1 = 400;
    } else { px = xsX; py = 0; tx = 1; ty = 0; ux = 0; uy = 1; u0 = 0; u1 = sim.Ly; }
    var inside = false;
    for (var u = u0 + st / 2; u < u1; u += st) {
      var c = cellXY(px + ux * u, py + uy * u); if (c < 0) continue;
      if (pathMode() && !inBasin(c)) { if (inside && u > 0) break; continue; } inside = true;
      Q += (sim.qx[c] * tx + sim.qy[c] * ty) * st;
      var i = c % nx, j = (c - i) / nx, fx = 0, fy = 0;
      for (var kk = 0; kk < K; kk++) { fx += 0.5 * (sim.SFx[(j * (nx + 1) + i) * K + kk] + sim.SFx[(j * (nx + 1) + i + 1) * K + kk]); fy += 0.5 * (sim.SFy[(j * nx + i) * K + kk] + sim.SFy[((j + 1) * nx + i) * K + kk]); }
      Qs += (fx * tx + fy * ty) * st;
    }
    return [Q, Qs];
  }
  function xsRecord() {
    if (!sim || sim.phase === 'spinup') return;
    var key = xsKey();
    if (!XH || XH.key !== key || XH.sim !== sim) XH = { key: key, sim: sim, t: [], Q: [], Qs: [] };
    if (XH.t.length && sim.t <= XH.t[XH.t.length - 1] + sim.duration() / 600) return;
    var f = xsFlux(); XH.t.push(sim.t); XH.Q.push(f[0]); XH.Qs.push(f[1]);
  }
  function xsName() { return pathMode() ? DP[PI].name.replace(/（.*$/, '') + ' ' + xsS.toFixed(0) + ' m' : 'x = ' + (+xsX).toFixed(0) + ' m'; }
  function drawSeries() {
    var S = sim.series, T = sim.duration() / 60, tm = S.t.map(function (t) { return t / 60; }), dem = !!sim.demMode && S.QinT;
    var px = [], pq = [], ps = [];
    if (dem) { px = tm; pq = S.QinT; ps = S.QsInT; }   // DEM：供給点・雨などの流入の合計（記録した値）
    else for (var i = 0; i <= 120; i++) { var t = sim.duration() * i / 120; px.push(t / 60); pq.push(sim.Qat(t)); ps.push(sim.QsAt(t)); }   // 入力ハイドログラフ（全期間）
    var xh = XH && XH.sim === sim && XH.t.length > 1 ? XH : null, xt = xh ? xh.t.map(function (t) { return t / 60; }) : null;
    var deb = sim.P.transport === 'debris', qlab = deb ? '（土石流は水と土砂の混合物）' : '';
    var sq = [{ x: px, y: pq, color: css('--s1'), label: '流入', w: 2 }, { x: tm, y: S.Qout, color: css('--s2'), label: '流出', w: 2 }];
    var kq = [{ color: css('--s1'), label: dem ? '流入（供給点・雨の合計）' : '上流端の流入' }, { color: css('--s2'), label: (dem ? '出口からの流出' : '下流端の流出') + qlab }];
    if (xh) { sq.push({ x: xt, y: xh.Q, color: css('--s3'), label: '横断', w: 1.6 }); kq.push({ color: css('--s3'), label: '横断位置（' + xsName() + '）を通る流量' }); }
    chart($('tsq'), { series: sq, x0: 0, x1: T, zero: true, xlabel: '時刻 [分]', xname: 't', xunit: ' 分', height: 200 });
    $('tsq_key').innerHTML = keyHTML(kq);
    var ss = [{ x: px, y: ps, color: css('--s1'), label: '供給', w: 2 }, { x: tm, y: S.QsOut, color: css('--s2'), label: '流出', w: 2 }], key = [{ color: css('--s1'), label: dem ? '供給（供給点の土砂）' : '上流端の供給' }, { color: css('--s2'), label: dem ? '出口からの流出' : '下流端の流出' }];
    if (xh) { ss.push({ x: xt, y: xh.Qs, color: css('--s3'), label: '横断', w: 1.6 }); key.push({ color: css('--s3'), label: '横断位置（' + xsName() + '）を通る流砂量' }); }
    if (ref) { ss.push({ x: ref.t, y: ref.QsOut, color: css('--sref'), dash: true, label: '比較', w: 1.6 }); key.push({ color: css('--sref'), dash: true, label: '比較の流出（' + ref.label + '）' }); }
    chart($('tss'), { series: ss, x0: 0, x1: T, zero: true, xlabel: '時刻 [分]', xname: 't', xunit: ' 分', height: 200 });
    $('tss_key').innerHTML = keyHTML(key);
  }
  function redrawCharts() { if (!sim) return; drawProfile(); drawSeries(); drawXsec(); }

  // ------------------------------------------------------------ 数値のまとめ
  function tile(k, v, u) { return '<div class="tile"><div class="k">' + k + '</div><div class="v">' + v + (u ? '<small>' + u + '</small>' : '') + '</div></div>'; }
  function drawStats() {
    var lam = sim.P.porosity, S = sim.series, pin = 0, pout = 0, i;
    if (sim.demMode && S.QsInT) { for (i = 0; i < S.QsInT.length; i++) if (S.QsInT[i] > pin) pin = S.QsInT[i]; }
    else for (i = 0; i < 200; i++) { var q = sim.QsAt(sim.t * i / 199); if (q > pin) pin = q; }
    for (i = 0; i < S.QsOut.length; i++) if (S.QsOut[i] > pout) pout = S.QsOut[i];
    var stor = sim.storage(), sIn = sim.sedIn / (1 - lam), sOut = sim.sedOut / (1 - lam);
    var h = tile('出水の時刻', (sim.t / 60).toFixed(0), '/ ' + (sim.duration() / 60).toFixed(0) + ' 分') +
      tile('流量（流入 → 流出）', (sim.demMode ? (sim.QinT || 0) : sim.Qin).toFixed(1) + ' → ' + sim.Qout.toFixed(1), 'm³/s') +
      tile('供給土砂量', Math.round(sIn).toLocaleString(), 'm³') +
      tile('流出土砂量', Math.round(sOut).toLocaleString(), 'm³') +
      tile('区間内の堆積量', (stor >= 0 ? '+' : '') + Math.round(stor).toLocaleString(), 'm³') +
      tile('流砂量ピーク（供給 → 流出）', pin.toFixed(2) + ' → ' + pout.toFixed(2), 'm³/s');
    if (ref) h += tile('比較との流出土砂量の差', ((sOut - ref.sOutAt(sim.t)) >= 0 ? '+' : '') + Math.round(sOut - ref.sOutAt(sim.t)).toLocaleString(), 'm³');
    $('tiles').innerHTML = h;
    var by = sim.storageByDam(), rows = '<tr><th>区間</th><th>位置 [m]</th><th>堤高 [m]</th><th>水通し幅 [m]</th><th>形式</th><th>堆砂量 [m³]</th><th>満砂時の目安 [m³]</th><th>充填率</th></tr>';
    sim.dams.forEach(function (D, q) {
      // 目安：元河床勾配 I、堆砂勾配 I/2 のくさび（幅は堰堤位置の川幅で一様と仮定）
      var Iloc = sim.Ix ? sim.Ix[Math.max(0, D.i - 1)] : sim.I, len = D.H / (Iloc / 2), up = q > 0 ? D.x - sim.dams[q - 1].x : D.x; var cap;
      if (len <= up) cap = 0.5 * D.H * len * sim.Bx[D.i]; else { var hu = D.H * (1 - up / len); cap = 0.5 * (D.H + hu) * up * sim.Bx[D.i]; }
      rows += '<tr><td>堰堤 ' + (q + 1) + ' の上流</td><td>' + D.x.toFixed(0) + '</td><td>' + D.H.toFixed(1) + '</td><td>' + (D.b >= sim.Bx[D.i] - 1e-6 ? '全幅' : D.b.toFixed(0)) + '</td><td>' + damTypeLabel(D) + '</td><td>' + Math.round(by[q]).toLocaleString() + '</td><td>' + Math.round(cap).toLocaleString() + '</td><td>' + (by[q] / cap * 100).toFixed(0) + '%</td></tr>';
    });
    rows += '<tr><td>' + (sim.dams.length ? '最下流の堰堤より下流' : '全区間') + '</td><td>–</td><td>–</td><td>–</td><td>–</td><td>' + Math.round(by[by.length - 1]).toLocaleString() + '</td><td>–</td><td>–</td></tr>';
    $('damres').innerHTML = rows;
    var st = sim.phase === 'spinup' ? '助走計算中（水だけを流して流れを立ち上げています）' : sim.phase === 'done' ? '計算終了' : (running ? '計算中' : '停止中');
    $('status').textContent = st + '　ステップ ' + sim.nstep.toLocaleString() + '　Δt ' + sim.dt.toFixed(3) + ' s' + (spd ? '　' + spd.toFixed(0) + ' ステップ/秒' : '');
    $('prog').style.width = (sim.phase === 'spinup' ? 0 : Math.min(100, sim.t / sim.duration() * 100)) + '%';
  }
  function saveRef() {
    if (!sim) return; var p = sim.profile(), lam = sim.P.porosity, S = sim.series, nx = sim.nx, i;
    var isDam = new Uint8Array(nx); sim.dams.forEach(function (D) { isDam[D.i] = 1; });
    var dz = new Float64Array(nx), dm = new Float64Array(nx); for (i = 0; i < nx; i++) { dz[i] = isDam[i] ? NaN : p.zmean[i] - p.z0mean[i]; dm[i] = p.dm[i] * 1000; }
    var tt = S.t.slice(), cum = [], acc = 0; for (i = 0; i < tt.length; i++) { acc += S.QsOut[i] * (i ? tt[i] - tt[i - 1] : tt[0]); cum.push(acc / (1 - lam)); }
    ref = { label: sim.dams.length ? '堰堤 ' + sim.dams.length + ' 基' : '堰堤なし', nx: nx, dx: sim.dx, zmean: p.zmean, zmin: p.zmin, isDam: isDam, dz: dz, dm: dm, qs: p.qs, wet: p.wet,
      t: tt.map(function (t) { return t / 60; }), QsOut: S.QsOut.slice(),
      sOutAt: function (t) { if (!tt.length) return 0; for (var i = tt.length - 1; i >= 0; i--) if (tt[i] <= t) return cum[i]; return 0; } };
    $('saveref').textContent = '比較用に保存済み（' + ref.label + '）'; redrawCharts(); drawStats();
  }

  // ------------------------------------------------------------ 鳥瞰（WebGL）
  var T3 = null;
  function init3D() {
    if (!window.THREE) { $('hint3d').textContent = '3D 表示のライブラリを読み込めませんでした。平面図と縦断図は使えます。'; $('c3d').hidden = true; return; }
    var THREE = window.THREE, cv = $('c3d'), r;
    try { r = new THREE.WebGLRenderer({ canvas: cv, antialias: true }); } catch (e) { $('hint3d').textContent = 'この環境では WebGL を使えません。平面図と縦断図は使えます。'; cv.hidden = true; return; }
    r.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
    var scene = new THREE.Scene(); scene.background = new THREE.Color(0xb9cfdd);
    var cam = new THREE.PerspectiveCamera(36, 2, 1, 6000);
    scene.add(new THREE.HemisphereLight(0xdfeaf2, 0x6b6250, 0.78));
    var sun = new THREE.DirectionalLight(0xfff3dd, 0.95); sun.position.set(-0.45, 0.8, 0.4); scene.add(sun);
    // 礫のきめ
    var tc = document.createElement('canvas'); tc.width = tc.height = 256; var tx = tc.getContext('2d'); tx.fillStyle = '#e4e4e4'; tx.fillRect(0, 0, 256, 256);
    var rr = 7; function rn() { rr = (rr * 16807) % 2147483647; return rr / 2147483647; }
    for (var q = 0; q < 2600; q++) { var g = 150 + rn() * 105 | 0, s = 1 + rn() * 4; tx.fillStyle = 'rgb(' + g + ',' + g + ',' + (g - 4) + ')'; tx.beginPath(); tx.ellipse(rn() * 256, rn() * 256, s, s * (0.6 + rn() * 0.5), rn() * 3, 0, 6.3); tx.fill(); }
    var tex = new THREE.CanvasTexture(tc); tex.wrapS = tex.wrapT = THREE.RepeatWrapping; tex.anisotropy = 4;
    var wmat = new THREE.ShaderMaterial({
      transparent: true, depthWrite: false, side: THREE.DoubleSide,
      uniforms: { uTime: { value: 0 }, uSun: { value: new THREE.Vector3(-0.45, 0.8, 0.4).normalize() }, uSh: { value: new THREE.Color() }, uDp: { value: new THREE.Color() }, uSky: { value: new THREE.Color(0xcfe0ea) }, uSci: { value: 0 } },
      vertexShader: 'attribute vec4 aD; attribute vec3 aC; varying vec3 vP; varying vec3 vN; varying vec4 vD; varying vec3 vC;' +
        'void main(){ vP=position; vN=normal; vD=aD; vC=aC; gl_Position=projectionMatrix*modelViewMatrix*vec4(position,1.0); }',
      fragmentShader: 'precision highp float; uniform float uTime,uSci; uniform vec3 uSun,uSh,uDp,uSky; varying vec3 vP; varying vec3 vN; varying vec4 vD; varying vec3 vC;' +
        'float hs(vec2 p){ return fract(sin(dot(p,vec2(127.1,311.7)))*43758.5453); }' +
        'float ns(vec2 p){ vec2 i=floor(p), f=fract(p); f=f*f*(3.0-2.0*f); return mix(mix(hs(i),hs(i+vec2(1,0)),f.x),mix(hs(i+vec2(0,1)),hs(i+vec2(1,1)),f.x),f.y); }' +
        'float fb(vec2 p){ return ns(p)*0.6+ns(p*2.3+7.1)*0.4; }' +
        'float flow(vec2 p, vec2 v){ float a=fract(uTime*0.35), b=fract(uTime*0.35+0.5); float w=abs(1.0-2.0*a); return mix(fb(p-v*a*2.2), fb(p-v*b*2.2+13.7), w); }' +
        'void main(){ float d=vD.x; if(d<0.006) discard; vec2 p=vP.xz*0.55; vec2 v=vD.zw;' +
        ' float e=0.35; float n0=flow(p,v), nx=flow(p+vec2(e,0.0),v), nz=flow(p+vec2(0.0,e),v);' +
        ' float amp=0.10+0.16*clamp(length(v)*0.25,0.0,1.0); vec3 N=normalize(normalize(vN)+vec3((n0-nx)*amp*6.0,0.0,(n0-nz)*amp*6.0));' +
        ' vec3 V=normalize(cameraPosition-vP); float fr=pow(1.0-max(dot(N,V),0.0),3.0)*0.55+0.06;' +
        ' vec3 base=mix(uSh,uDp,1.0-exp(-d*1.5)); float dif=0.72+0.28*max(dot(N,uSun),0.0);' +
        ' vec3 col=mix(base*dif,uSky,fr); float sp=pow(max(dot(N,normalize(uSun+V)),0.0),90.0)*0.7;' +
        ' float fo=clamp(vD.y*1.5-0.55+ (n0-0.5)*1.3,0.0,1.0)*step(0.01,vD.y); col=mix(col,vec3(0.96,0.97,0.96),fo); col+=sp;' +
        ' float al=clamp(0.42+d*1.1,0.0,0.93); al=max(al,fo*0.95); al*=smoothstep(0.006,0.05,d);' +
        ' if(uSci>0.5){ col=vC*(0.8+0.2*max(dot(normalize(vN),uSun),0.0)); al=(uSci>1.5?0.10:0.96)*smoothstep(0.006,0.04,d); }' +
        ' gl_FragColor=vec4(col,al); }'
    });
    // 河床の質感：表層の粒径で「さらさら（砂）」と「ごつごつ（礫）」を描き分ける
    var tmat = new THREE.ShaderMaterial({
      uniforms: { uSun: { value: new THREE.Vector3(-0.45, 0.8, 0.4).normalize() }, uReal: { value: 1 }, uCI: { value: 0 }, uFlat: { value: 0 } },
      vertexShader: 'attribute vec3 aCol; attribute vec4 aG; varying vec3 vCol; varying vec4 vG; varying vec3 vN; varying vec3 vW;' +
        'void main(){ vCol=aCol; vG=aG; vN=normal; vW=position; gl_Position=projectionMatrix*modelViewMatrix*vec4(position,1.0); }',
      fragmentShader: 'precision highp float; uniform vec3 uSun; uniform float uReal,uCI,uFlat; varying vec3 vCol; varying vec4 vG; varying vec3 vN; varying vec3 vW;' +
        'float h21(vec2 p){ p=fract(p*vec2(123.34,456.21)); p+=dot(p,p+45.32); return fract(p.x*p.y); }' +
        'vec2 h22(vec2 p){ float n=h21(p); return vec2(n,h21(p+n+17.0)); }' +
        'void main(){ vec3 N=normalize(vN); vec3 base=vCol; float k=vG.x;' +
        ' if(uReal>0.5){ float s=max(vG.y,0.004); vec2 p=vW.xz/s; vec2 fw=fwidth(p); float vis=1.0-smoothstep(0.35,1.0,max(fw.x,fw.y));' +
        '  vec2 ip=floor(p), fp=fract(p); float f1=8.0,f2=8.0,id=0.0; vec2 dv=vec2(0.0);' +
        '  for(int j=-1;j<=1;j++){ for(int i=-1;i<=1;i++){ vec2 g=vec2(float(i),float(j)); vec2 o=h22(ip+g)*0.75+0.125; vec2 r=g+o-fp; float d=dot(r,r);' +
        '   if(d<f1){ f2=f1; f1=d; dv=r; id=h21(ip+g+3.7);} else if(d<f2){ f2=d; } } }' +
        '  f1=sqrt(f1); f2=sqrt(f2); float crev=1.0-smoothstep(0.03,0.16,f2-f1);' +
        '  vec3 stone=base*(0.76+0.52*id)*vec3(1.0+0.08*(id-0.5),1.0,1.0-0.08*(id-0.5)); stone*=1.0-0.62*crev;' +
        '  vec2 q=vW.xz/(s*5.0); vec2 fq=fwidth(q); float vq=1.0-smoothstep(0.4,1.2,max(fq.x,fq.y)); vec3 far=base*(0.9+0.32*(h21(floor(q))-0.5)*vq);' +
        '  vec3 cob=mix(far,stone,vis);' +
        '  vec2 ps=vW.xz/0.04; vec2 fs=fwidth(ps); float vs=1.0-smoothstep(0.4,1.0,max(fs.x,fs.y)); vec3 sand=base*(1.0+0.09*(h21(floor(ps))-0.5)*vs);' +
        '  float kk=smoothstep(0.1,0.9,k); base=mix(sand,cob,kk);' +
        '  vec2 b=-dv*(0.6+1.8*f1)*vis*kk; N=normalize(N+vec3(b.x,0.0,b.y)*0.9); }' +
        ' float dif=max(dot(N,uSun),0.0); float hemi=0.5+0.5*N.y; vec3 col=base*mix(0.50*hemi+0.70*dif+0.06,0.55+0.45*(0.5*hemi+0.7*dif),uFlat);' +
        ' if(uCI>0.0){ float v=vG.z/uCI; float w=max(fwidth(v),1e-4); float f=abs(fract(v-0.5)-0.5)/w; col=mix(col,col*0.35,(1.0-min(f,1.0))*0.9); }' +
        ' gl_FragColor=vec4(col,1.0); }'
    });
    tmat.extensions.derivatives = true;
    // 礫（粗い粒径階を実体の石として置く）
    var rg = new THREE.IcosahedronGeometry(1, 0), rpos = rg.attributes.position;
    for (var vi = 0; vi < rpos.count; vi++) {
      var vx = rpos.getX(vi), vy = rpos.getY(vi), vz = rpos.getZ(vi), key = Math.round(vx * 97) * 7 + Math.round(vy * 89) * 13 + Math.round(vz * 83) * 31;
      var jf = 0.78 + 0.4 * hash(key); rpos.setXYZ(vi, vx * jf, vy * jf, vz * jf);
    }
    rg.computeVertexNormals();
    var RMAX = 40000, rocks = new THREE.InstancedMesh(rg, new THREE.MeshLambertMaterial({ color: 0xffffff }), RMAX);
    rocks.instanceMatrix.setUsage(THREE.DynamicDrawUsage); rocks.frustumCulled = false;
    rocks.setColorAt(0, new THREE.Color(1, 1, 1));   // 色の配列は count が RMAX の間に作る
    rocks.instanceColor.setUsage(THREE.DynamicDrawUsage); rocks.count = 0;
    scene.add(rocks);
    var RV = [];   // 形のばらつき（回転＋扁平）
    for (var q = 0; q < 64; q++) {
      var yaw = hash(q * 3.1) * 6.283, tl = (hash(q * 5.7) - 0.5) * 0.5, sx = 0.8 + 0.4 * hash(q * 7.3), sy = 0.45 + 0.25 * hash(q * 9.1), sz = 0.8 + 0.4 * hash(q * 11.9);
      var m = new THREE.Matrix4().makeRotationY(yaw).multiply(new THREE.Matrix4().makeRotationX(tl)).multiply(new THREE.Matrix4().makeScale(sx, sy, sz));
      RV.push({ e: m.elements.slice(), sy: sy });
    }
    // 樹木（谷壁の森林。表示だけ）：針葉樹は円錐、広葉樹は丸い樹冠
    var conG = new THREE.ConeGeometry(1, 1, 7); conG.translate(0, 0.5, 0);
    var brG = new THREE.IcosahedronGeometry(1, 1), bp = brG.attributes.position;
    for (var vi2 = 0; vi2 < bp.count; vi2++) { var bx = bp.getX(vi2), by = bp.getY(vi2), bz = bp.getZ(vi2), jf2 = 0.85 + 0.3 * hash(Math.round(bx * 50) * 7 + Math.round(by * 50) * 13 + Math.round(bz * 50) * 31); bp.setXYZ(vi2, bx * jf2, by * jf2 * 0.8 + 0.8, bz * jf2); }
    brG.computeVertexNormals();
    var TMAX = 12000, mkTrees = function (geo) { var m = new THREE.InstancedMesh(geo, new THREE.MeshLambertMaterial({ color: 0xffffff }), TMAX); m.setColorAt(0, new THREE.Color(1, 1, 1)); m.count = 0; m.frustumCulled = false; scene.add(m); return m; };
    var treesC = mkTrees(conG), treesB = mkTrees(brG);
    T3 = { THREE: THREE, r: r, scene: scene, treesC: treesC, treesB: treesB, TMAX: TMAX, treeKey: '', cam: cam, tex: tex, wmat: wmat, tmat: tmat, rocks: rocks, RMAX: RMAX, RV: RV, th: -0.62, ph: 0.62, dist: 1, tgt: new THREE.Vector3(), group: null, t0: performance.now() };
    // 視点操作：目標の視点（g〜）へ毎フレーム滑らかに近づける
    T3.gth = T3.th; T3.gph = T3.ph; T3.gdist = T3.dist; T3.gtgt = new THREE.Vector3();
    var ray = new THREE.Raycaster(), ndc = new THREE.Vector2(), ptr = {}, pinch = null;
    function hitGround(cx, cy) {   // 画面上の点に当たる地表の位置
      var r = cv.getBoundingClientRect(); ndc.set((cx - r.left) / r.width * 2 - 1, -(cy - r.top) / r.height * 2 + 1);
      ray.setFromCamera(ndc, cam);
      var hit = T3.terr ? ray.intersectObject(T3.terr, false) : [];
      if (hit.length) return hit[0].point.clone();
      var o = ray.ray.origin, d = ray.ray.direction; if (Math.abs(d.y) < 1e-6) return null;
      var t = (T3.tgt.y - o.y) / d.y; return t > 0 ? o.clone().add(d.clone().multiplyScalar(t)) : null;
    }
    function clampGoal() {
      if (!sim) return;
      var big = Math.max(sim.Lx, sim.Ly);
      T3.gdist = Math.max(4, Math.min(big * 3, T3.gdist));
      T3.gph = Math.max(0.04, Math.min(1.55, T3.gph));
      T3.gtgt.x = Math.max(-sim.Lx * 0.6, Math.min(sim.Lx * 0.6, T3.gtgt.x));
      T3.gtgt.z = Math.max(-sim.Ly * 0.6 - 20, Math.min(sim.Ly * 0.6 + 20, T3.gtgt.z));
    }
    function pan(dxp, dyp) {      // 地図のように地面をつかんで動かす
      var k = T3.gdist / cv.clientHeight * 1.15, kf = k / Math.max(0.35, Math.sin(T3.gph));
      var rx = Math.sin(T3.gth), rz = -Math.cos(T3.gth), fx = -Math.cos(T3.gth), fz = -Math.sin(T3.gth);
      T3.gtgt.x += -rx * dxp * k + fx * dyp * kf; T3.gtgt.z += -rz * dxp * k + fz * dyp * kf; clampGoal();
    }
    function rotate(dxp, dyp) { T3.gth -= dxp * 0.006; T3.gph += dyp * 0.005; clampGoal(); }
    function zoomAt(f, cx, cy) {  // カーソルの位置へ向かって拡大・縮小
      var P = cx != null ? hitGround(cx, cy) : null;
      if (P) { T3.gtgt.x = P.x + (T3.gtgt.x - P.x) * f; T3.gtgt.z = P.z + (T3.gtgt.z - P.z) * f; T3.gtgt.y = P.y + (T3.gtgt.y - P.y) * f; }
      T3.gdist *= f; clampGoal();
    }
    T3.zoomAt = zoomAt;
    cv.addEventListener('contextmenu', function (e) { e.preventDefault(); });
    cv.addEventListener('pointerdown', function (e) {
      cv.setPointerCapture(e.pointerId); cv.focus({ preventScroll: true });
      ptr[e.pointerId] = { x: e.clientX, y: e.clientY, pan: e.button === 2 || e.button === 1 || e.shiftKey || e.ctrlKey };
      cv.style.cursor = ptr[e.pointerId].pan ? 'move' : 'grabbing';
      var ids = Object.keys(ptr); if (ids.length === 2) { var a = ptr[ids[0]], b = ptr[ids[1]]; pinch = { d: Math.hypot(a.x - b.x, a.y - b.y), mx: (a.x + b.x) / 2, my: (a.y + b.y) / 2 }; }
    });
    cv.addEventListener('pointermove', function (e) {
      var p = ptr[e.pointerId]; if (!p) return;
      var ids = Object.keys(ptr), dxp = e.clientX - p.x, dyp = e.clientY - p.y;
      p.x = e.clientX; p.y = e.clientY;
      if (ids.length >= 2) {      // 2本指：拡大縮小と移動
        var a = ptr[ids[0]], b = ptr[ids[1]], d = Math.hypot(a.x - b.x, a.y - b.y), mx = (a.x + b.x) / 2, my = (a.y + b.y) / 2;
        if (pinch && d > 10) { zoomAt(pinch.d / d, mx, my); pan(mx - pinch.mx, my - pinch.my); }
        pinch = { d: d, mx: mx, my: my }; return;
      }
      if (p.pan || e.shiftKey) pan(dxp, dyp); else rotate(dxp, dyp);
    });
    var up = function (e) { delete ptr[e.pointerId]; if (Object.keys(ptr).length < 2) pinch = null; cv.style.cursor = 'grab'; };
    cv.addEventListener('pointerup', up); cv.addEventListener('pointercancel', up);
    cv.addEventListener('wheel', function (e) {
      e.preventDefault();
      var dy = e.deltaY * (e.deltaMode === 1 ? 33 : e.deltaMode === 2 ? 400 : 1);
      zoomAt(Math.exp(Math.max(-0.5, Math.min(0.5, dy * 0.0015))), e.clientX, e.clientY);
    }, { passive: false });
    cv.addEventListener('dblclick', function (e) { var P = hitGround(e.clientX, e.clientY); if (P) { T3.gtgt.copy(P); T3.gdist *= 0.6; clampGoal(); } });
    cv.addEventListener('keydown', function (e) {
      var st = 40, used = true;
      if (e.key === 'ArrowLeft') e.shiftKey ? rotate(-st, 0) : pan(st, 0);
      else if (e.key === 'ArrowRight') e.shiftKey ? rotate(st, 0) : pan(-st, 0);
      else if (e.key === 'ArrowUp') e.shiftKey ? rotate(0, -st) : pan(0, st);
      else if (e.key === 'ArrowDown') e.shiftKey ? rotate(0, st) : pan(0, -st);
      else if (e.key === '+' || e.key === '=' || e.key === ';') zoomAt(0.8);
      else if (e.key === '-' || e.key === '_') zoomAt(1.25);
      else if (e.key === 'Home') home3D();
      else used = false;
      if (used) e.preventDefault();
    });
    $('v_zi').onclick = function () { zoomAt(0.8); }; $('v_zo').onclick = function () { zoomAt(1.25); };
    $('v_home').onclick = function () { $('v_view').value = 'home'; home3D(); };
    $('v_view').onchange = function () { setView(this.value); };
  }
  // 決まった向きの視点（instant = true ですぐに切り替える）
  function setView(v, instant) {
    if (!T3 || !sim) return;
    var L = sim.Lx, W = sim.Ly, f = 2 * Math.tan(T3.cam.fov * Math.PI / 360), asp = T3.cam.aspect || 2;
    var fit = Math.max(L / asp, W) / f * 1.08;   // 区間全体が収まる距離（真上から）
    var V = { home: [0.42, 0.78, Math.max(L * 0.78, W * 2.2)], top: [Math.PI / 2, 1.55, fit], up: [Math.PI, 0.32, L * 0.55],
      down: [0, 0.32, L * 0.55], left: [-Math.PI / 2, 0.42, Math.max(L * 0.75, W * 1.6)] }[v] || [0.42, 0.78, L * 0.78];
    T3.gth = V[0]; T3.gph = V[1]; T3.gdist = V[2]; T3.gtgt.set(0, 0, 0);
    if (instant) { T3.th = T3.gth; T3.ph = T3.gph; T3.dist = T3.gdist; T3.tgt.copy(T3.gtgt); }
  }
  function home3D(instant) { setView('home', instant); }
  function build3D() {
    if (!T3) return; var THREE = T3.THREE, nx = sim.nx, ny = sim.ny, N = sim.N, dx = sim.dx;
    if (T3.group) { T3.scene.remove(T3.group); T3.group.traverse(function (o) { if (o.geometry) o.geometry.dispose(); }); }
    var grp = T3.group = new THREE.Group(); T3.scene.add(grp);
    var idx = new Uint32Array((nx - 1) * (ny - 1) * 6), k = 0, i, j;
    for (j = 0; j < ny - 1; j++) for (i = 0; i < nx - 1; i++) { var a = j * nx + i, b = a + 1, c = a + nx, d = c + 1; idx[k++] = a; idx[k++] = b; idx[k++] = c; idx[k++] = b; idx[k++] = d; idx[k++] = c; }
    var tg = new THREE.BufferGeometry();
    tg.setAttribute('position', new THREE.BufferAttribute(new Float32Array(N * 3), 3)); tg.setAttribute('aCol', new THREE.BufferAttribute(new Float32Array(N * 3), 3));
    tg.setAttribute('aG', new THREE.BufferAttribute(new Float32Array(N * 4), 4)); tg.setIndex(new THREE.BufferAttribute(idx, 1));
    T3.terr = new THREE.Mesh(tg, T3.tmat); grp.add(T3.terr);
    var wg = new THREE.BufferGeometry();
    wg.setAttribute('position', new THREE.BufferAttribute(new Float32Array(N * 3), 3)); wg.setAttribute('normal', new THREE.BufferAttribute(new Float32Array(N * 3), 3));
    wg.setAttribute('aD', new THREE.BufferAttribute(new Float32Array(N * 4), 4)); wg.setAttribute('aC', new THREE.BufferAttribute(new Float32Array(N * 3), 3));
    wg.setIndex(new THREE.BufferAttribute(idx, 1));
    T3.water = new THREE.Mesh(wg, T3.wmat); T3.water.frustumCulled = false; T3.water.renderOrder = 2; grp.add(T3.water);
    // 側面（地盤の断面）
    var per = 2 * (nx + ny) - 4, sg = new THREE.BufferGeometry(), si = new Uint32Array(per * 6);
    for (i = 0; i < per; i++) { var n = (i + 1) % per; si.set([i * 2, i * 2 + 1, n * 2, n * 2, i * 2 + 1, n * 2 + 1], i * 6); }
    sg.setAttribute('position', new THREE.BufferAttribute(new Float32Array(per * 6), 3)); sg.setIndex(new THREE.BufferAttribute(si, 1));
    T3.skirt = new THREE.Mesh(sg, new THREE.MeshLambertMaterial({ color: 0x6f6252, side: THREE.DoubleSide })); grp.add(T3.skirt);
    T3.ring = []; for (i = 0; i < nx; i++) T3.ring.push(i); for (j = 1; j < ny; j++) T3.ring.push(j * nx + nx - 1); for (i = nx - 2; i >= 0; i--) T3.ring.push((ny - 1) * nx + i); for (j = ny - 2; j >= 1; j--) T3.ring.push(j * nx);
    // 堰堤（コンクリート）
    var dg = new THREE.BufferGeometry(); T3.damSegs = []; T3.steelSegs = [];
    sim.dams.forEach(function (D) { damSegs(D).forEach(function (g) { (g.steel ? T3.steelSegs : T3.damSegs).push({ D: D, g: g }); }); });
    var sgm = new THREE.BufferGeometry();
    sgm.setAttribute('position', new THREE.BufferAttribute(new Float32Array(Math.max(1, T3.steelSegs.length) * 24 * 3), 3));
    var si2 = new Uint32Array(Math.max(1, T3.steelSegs.length) * 36), F2 = [[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11], [12, 13, 14, 15], [16, 17, 18, 19], [20, 21, 22, 23]];
    T3.steelSegs.forEach(function (c, q) { F2.forEach(function (f, m) { si2.set([q * 24 + f[0], q * 24 + f[1], q * 24 + f[2], q * 24 + f[0], q * 24 + f[2], q * 24 + f[3]], q * 36 + m * 6); }); });
    sgm.setIndex(new THREE.BufferAttribute(si2, 1));
    T3.steel = new THREE.Mesh(sgm, new THREE.MeshLambertMaterial({ color: 0x55606c, side: THREE.DoubleSide })); grp.add(T3.steel);
    dg.setAttribute('position', new THREE.BufferAttribute(new Float32Array(T3.damSegs.length * 24 * 3), 3));
    var di = new Uint32Array(T3.damSegs.length * 36), F = [[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11], [12, 13, 14, 15], [16, 17, 18, 19], [20, 21, 22, 23]];
    T3.damSegs.forEach(function (c, q) { F.forEach(function (f, m) { di.set([q * 24 + f[0], q * 24 + f[1], q * 24 + f[2], q * 24 + f[0], q * 24 + f[2], q * 24 + f[3]], q * 36 + m * 6); }); });
    dg.setIndex(new THREE.BufferAttribute(di, 1));
    T3.dam = new THREE.Mesh(dg, new THREE.MeshLambertMaterial({ color: 0xbdbbb2, side: THREE.DoubleSide })); grp.add(T3.dam);
    T3.vexBuilt = 0;
  }
  var c3 = [0, 0, 0];
  function update3D() {
    if (!T3 || !T3.terr) return;
    var nx = sim.nx, ny = sim.ny, N = sim.N, dx = sim.dx, z = sim.z, h = sim.h, kind = sim.kind, vex = +$('v_vex').value, mode = $('v_mode').value, V = VARS[mode];
    var terr = !!V.terr, tcol = terr ? terrColors(mode) : null;
    var zmid = sim.zmid, X0 = sim.Lx / 2, Z0 = sim.Ly / 2;
    var tp = T3.terr.geometry.attributes.position.array, tcn = T3.terr.geometry.attributes.aCol.array, tgA = T3.terr.geometry.attributes.aG.array;
    var S = +$('v_rock').value, Stex = S > 0 ? S : 5, ci = contourInfo(mode), gs = [0, 0];
    var wp = T3.water.geometry.attributes.position.array, wn = T3.water.geometry.attributes.normal.array, wd = T3.water.geometry.attributes.aD.array, wcl = T3.water.geometry.attributes.aC.array;
    var i, j, c, p;
    for (j = 0; j < ny; j++) for (i = 0; i < nx; i++) {
      c = j * nx + i; p = c * 3;
      var zz = kind[c] === 2 ? damGround(c) : (sim.zOut && !sim.basin[c] ? sim.zOut[c] : z[c]);   // DEM の流域の外は元の標高で見せる
      tp[p] = (i + 0.5) * dx - X0; tp[p + 1] = (zz - zmid) * vex; tp[p + 2] = Z0 - (j + 0.5) * dx;
      if (mode === 'real') { ground(kind[c] === 2 ? nonDam(c) : c, c3); }
      else if (terr) { c3[0] = tcol[p]; c3[1] = tcol[p + 1]; c3[2] = tcol[p + 2]; }
      else { var v = V.water ? NaN : value(mode, c); if (v === v) V.ramp(norm(mode, v), c3); else { c3[0] = c3[1] = kind[c] === 1 ? 120 : 150; c3[2] = c3[0] + 3; } }
      tcn[p] = c3[0] / 255; tcn[p + 1] = c3[1] / 255; tcn[p + 2] = c3[2] / 255;
      var g4 = c * 4;
      if (kind[c] !== 2 && gstat(c, gs)) { tgA[g4] = coarse(gs[0]); tgA[g4 + 1] = gs[1] * Stex; } else { tgA[g4] = 0; tgA[g4 + 1] = 0.02; }
      tgA[g4 + 2] = ci ? (ci.rel ? zz - trend((i + 0.5) * dx) : zz) : 0;
      // 水面
      var hh = h[c], eta;
      if (hh > 0.006) eta = z[c] + hh;
      else {
        eta = 1e30;
        if (i > 0 && h[c - 1] > 0.006) eta = Math.min(eta, z[c - 1] + h[c - 1]); if (i < nx - 1 && h[c + 1] > 0.006) eta = Math.min(eta, z[c + 1] + h[c + 1]);
        if (j > 0 && h[c - nx] > 0.006) eta = Math.min(eta, z[c - nx] + h[c - nx]); if (j < ny - 1 && h[c + nx] > 0.006) eta = Math.min(eta, z[c + nx] + h[c + nx]);
        // 水のないセルの水面は、描いている地面（DEM の流域の外は元の標高）より下に置く。計算の z（流域の外は高い壁）を使うと、
        // 濡れたセルとの三角形が壁の高さまで伸びて、水の柱のように見えるため
        eta = eta > 1e29 ? zz - 0.4 : Math.min(eta, zz - 0.02); hh = 0;
      }
      wp[p] = tp[p]; wp[p + 1] = (eta - zmid) * vex + 0.01; wp[p + 2] = tp[p + 2];
      var uu = hh > 0.01 ? sim.qx[c] / hh : 0, vv = hh > 0.01 ? sim.qy[c] / hh : 0, fr = hh > 0.01 ? Math.hypot(uu, vv) / Math.sqrt(G * hh) : 0;
      wd[c * 4] = hh; wd[c * 4 + 1] = kind[c] === 2 ? 0 : Math.max(0, Math.min(1, (fr - 0.95) * 0.9)); wd[c * 4 + 2] = uu; wd[c * 4 + 3] = -vv;
      if (mode !== 'real' && !terr && V.water) { var wv = value(mode, c); if (wv === wv) V.ramp(norm(mode, wv), c3); else { c3[0] = c3[1] = c3[2] = 200; } wcl[p] = c3[0] / 255; wcl[p + 1] = c3[1] / 255; wcl[p + 2] = c3[2] / 255; }
      else { wcl[p] = 0.35; wcl[p + 1] = 0.6; wcl[p + 2] = 0.85; }
    }
    for (j = 0; j < ny; j++) for (i = 0; i < nx; i++) {
      c = j * nx + i; p = c * 3;
      var ex = (wp[(i < nx - 1 ? c + 1 : c) * 3 + 1] - wp[(i > 0 ? c - 1 : c) * 3 + 1]) / (2 * dx), ez = (wp[(j > 0 ? c - nx : c) * 3 + 1] - wp[(j < ny - 1 ? c + nx : c) * 3 + 1]) / (2 * dx);
      var il = 1 / Math.sqrt(ex * ex + ez * ez + 1); wn[p] = -ex * il; wn[p + 1] = il; wn[p + 2] = -ez * il;
    }
    ['position', 'aCol', 'aG'].forEach(function (a) { T3.terr.geometry.attributes[a].needsUpdate = true; });
    T3.terr.geometry.computeBoundingSphere(); T3.terr.geometry.computeBoundingBox();
    T3.terr.geometry.computeVertexNormals();
    ['position', 'normal', 'aD', 'aC'].forEach(function (a) { T3.water.geometry.attributes[a].needsUpdate = true; });
    T3.tmat.uniforms.uReal.value = mode === 'real' ? 1 : 0; T3.tmat.uniforms.uCI.value = ci ? ci.ci : 0; T3.tmat.uniforms.uFlat.value = terr ? 0.7 : 0;
    placeRocks(mode === 'real' ? S : 0, tp, zmid, vex);
    placeTrees(mode === 'real' && $('v_tree').checked, tp, vex);
    var wc = waterCol(); T3.wmat.uniforms.uSh.value.setRGB(wc.sh[0] / 255, wc.sh[1] / 255, wc.sh[2] / 255); T3.wmat.uniforms.uDp.value.setRGB(wc.dp[0] / 255, wc.dp[1] / 255, wc.dp[2] / 255);
    T3.wmat.uniforms.uSci.value = mode === 'real' || terr ? 0 : (V.water ? 1 : 2);
    if (T3.vexBuilt !== vex) {
      T3.vexBuilt = vex;
      var base = (sim.zlo - zmid - Math.max(1.5, sim.P.erodDepth + 1.5)) * vex, sp = T3.skirt.geometry.attributes.position.array;
      T3.ring.forEach(function (c, q) { sp[q * 6] = tp[c * 3]; sp[q * 6 + 1] = tp[c * 3 + 1]; sp[q * 6 + 2] = tp[c * 3 + 2]; sp[q * 6 + 3] = tp[c * 3]; sp[q * 6 + 4] = base; sp[q * 6 + 5] = tp[c * 3 + 2]; });
      T3.skirt.geometry.attributes.position.needsUpdate = true; T3.skirt.geometry.computeVertexNormals();
      var dp = T3.dam.geometry.attributes.position.array, hw = dx * 0.5;
      var X0d = sim.Lx / 2, Z0d = sim.Ly / 2;
      var fillSeg = function (arr, sg, q) {
        var D = sg.D, g = sg.g, top = (g.top - zmid) * vex, bot = ((g.bot != null ? g.bot : D.zlow - sim.P.erodDepth - 0.5) - zmid) * vex, ht = g.steel ? Math.min(dx * 0.25, 0.4) : dx * (g.spill ? 0.2 : 0.42);
        var cor = [[g.s0, -ht], [g.s0, ht], [g.s1, ht], [g.s1, -ht]].map(function (u) {
          var x = D.cx + D.nx * u[0] + D.tx * u[1], y = D.cy + D.ny * u[0] + D.ty * u[1]; return [x - X0d, Z0d - y];
        });
        var P8 = [[cor[0][0], bot, cor[0][1]], [cor[1][0], bot, cor[1][1]], [cor[2][0], bot, cor[2][1]], [cor[3][0], bot, cor[3][1]],
          [cor[0][0], top, cor[0][1]], [cor[1][0], top, cor[1][1]], [cor[2][0], top, cor[2][1]], [cor[3][0], top, cor[3][1]]];
        var Fq = [[4, 7, 6, 5], [0, 1, 2, 3], [0, 4, 5, 1], [3, 2, 6, 7], [0, 3, 7, 4], [1, 5, 6, 2]];
        Fq.forEach(function (f, m) { f.forEach(function (vi, n) { arr.set(P8[vi], (q * 24 + m * 4 + n) * 3); }); });
      };
      if (T3.steelSegs.length) { var spa = T3.steel.geometry.attributes.position.array; T3.steelSegs.forEach(function (sg, q) { fillSeg(spa, sg, q); }); T3.steel.geometry.attributes.position.needsUpdate = true; T3.steel.geometry.computeVertexNormals(); }
      T3.damSegs.forEach(function (sg, q) {
        var D = sg.D, g = sg.g, top = (g.top - zmid) * vex, bot = (D.zlow - zmid - sim.P.erodDepth - 0.5) * vex, ht = dx * (g.spill ? 0.2 : 0.42);
        var cor = [[g.s0, -ht], [g.s0, ht], [g.s1, ht], [g.s1, -ht]].map(function (u) {   // (堰堤軸方向 s, 流下方向 a)
          var x = D.cx + D.nx * u[0] + D.tx * u[1], y = D.cy + D.ny * u[0] + D.ty * u[1]; return [x - X0d, Z0d - y];
        });
        var P8 = [[cor[0][0], bot, cor[0][1]], [cor[1][0], bot, cor[1][1]], [cor[2][0], bot, cor[2][1]], [cor[3][0], bot, cor[3][1]],
          [cor[0][0], top, cor[0][1]], [cor[1][0], top, cor[1][1]], [cor[2][0], top, cor[2][1]], [cor[3][0], top, cor[3][1]]];
        var Fq = [[4, 7, 6, 5], [0, 1, 2, 3], [0, 4, 5, 1], [3, 2, 6, 7], [0, 3, 7, 4], [1, 5, 6, 2]];
        Fq.forEach(function (f, m) { f.forEach(function (vi, n) { dp.set(P8[vi], (q * 24 + m * 4 + n) * 3); }); });
      });
      T3.dam.geometry.attributes.position.needsUpdate = true; T3.dam.geometry.computeVertexNormals();
    }
  }
  // 谷壁の森林を置く（表示だけ。地形が変わらない限り置き直さない）
  function placeTrees(on, tp, vex) {
    var C = T3.treesC, Bm = T3.treesB;
    if (!on) { C.visible = Bm.visible = false; return; }
    C.visible = Bm.visible = true;
    var nx = sim.nx, ny = sim.ny, dx = sim.dx, kind = sim.kind, z = sim.z, z0 = sim.z0, bare = sim.bare;
    var key = sim.N + '|' + vex + '|' + sim.nstep;
    if (T3.treeSim === sim && T3.treeVex === vex && (sim.nstep - T3.treeStep) < 400) return;
    T3.treeSim = sim; T3.treeVex = vex; T3.treeStep = sim.nstep;
    var X0 = sim.Lx / 2, Z0 = sim.Ly / 2, dens = Math.min(1, dx * dx / 16), nc = 0, nb = 0, ac = C.instanceMatrix.array, ab = Bm.instanceMatrix.array, cc = C.instanceColor.array, cb = Bm.instanceColor.array;
    for (var j = 1; j < ny - 1; j++) for (var i = 1; i < nx - 1; i++) {
      var c = j * nx + i;
      var bankV = sim.bankC && sim.bankC[c] && Math.abs(z[c] - z0[c]) < 0.05;
      if (!bankV && (kind[c] !== 1 || (bare && bare[c]) || z[c] - z0[c] > 0.02)) continue;
      // 谷底のすぐ脇（岸の際）は疎らにする
      var nearFloor = kind[c - 1] === 0 || kind[c + 1] === 0 || kind[c - nx] === 0 || kind[c + nx] === 0;
      var nt = dens * (nearFloor ? 0.45 : 1), cnt = Math.floor(nt + hash(c * 17 + 5));
      for (var m = 0; m < cnt; m++) {
        var key2 = c * 29 + m * 7, fx = i + hash(key2 + 0.11) - 0.5, fy = j + hash(key2 + 0.57) - 0.5;
        var i0 = Math.min(nx - 2, Math.max(0, fx | 0)), j0 = Math.min(ny - 2, Math.max(0, fy | 0)), tx = fx - i0, ty = fy - j0, b00 = (j0 * nx + i0) * 3;
        var y = (tp[b00 + 1] * (1 - tx) + tp[b00 + 4] * tx) * (1 - ty) + (tp[b00 + nx * 3 + 1] * (1 - tx) + tp[b00 + nx * 3 + 4] * tx) * ty;
        var X = (fx + 0.5) * dx - X0, Zp = Z0 - (fy + 0.5) * dx, con = hash(key2 + 0.33) < 0.55, tone = 0.75 + 0.4 * hash(key2 + 0.77);
        if (con && nc < T3.TMAX) {
          var hgt = (7 + 6 * hash(key2 + 0.91)), rad = hgt * (0.26 + 0.06 * hash(key2 + 0.21)), o = nc * 16;
          ac.set([rad, 0, 0, 0, 0, hgt, 0, 0, 0, 0, rad, 0, X, y - 0.3, Zp, 1], o);
          cc[nc * 3] = 0.13 * tone; cc[nc * 3 + 1] = 0.27 * tone; cc[nc * 3 + 2] = 0.17 * tone; nc++;
        } else if (!con && nb < T3.TMAX) {
          var rb = 2.2 + 1.6 * hash(key2 + 0.43), o2 = nb * 16, yw = hash(key2 + 0.63) * 6.28, cs = Math.cos(yw), sn = Math.sin(yw);
          ab.set([rb * cs, 0, -rb * sn, 0, 0, rb * 1.1, 0, 0, rb * sn, 0, rb * cs, 0, X, y - 0.2, Zp, 1], o2);
          var yl = hash(key2 + 0.19);
          cb[nb * 3] = (0.24 + 0.08 * yl) * tone; cb[nb * 3 + 1] = (0.40 + 0.06 * yl) * tone; cb[nb * 3 + 2] = 0.19 * tone; nb++;
        }
      }
    }
    C.count = nc; Bm.count = nb;
    C.instanceMatrix.needsUpdate = Bm.instanceMatrix.needsUpdate = true; C.instanceColor.needsUpdate = Bm.instanceColor.needsUpdate = true;
  }
  // 粗い粒径階を石として配置する（S：表示倍率、0 で表示しない）
  function placeRocks(S, tp, zmid, vex) {
    var R = T3.rocks; if (!(S > 0)) { R.count = 0; R.visible = false; return; }
    R.visible = true;
    var nx = sim.nx, ny = sim.ny, dx = sim.dx, K = sim.K, A = sim.A, d = sim.d, kind = sim.kind, h = sim.h;
    var area = dx * dx, cover = 0.4, cls = [], k, c, i, j;
    for (k = 0; k < K; k++) if (d[k] >= 0.015 && d[k] * S >= 0.06) cls.push(k);
    if (!cls.length) { R.count = 0; return; }
    var tot = 0;
    for (c = 0; c < sim.N; c++) {
      if (kind[c] === 2) continue; var ck = c * K, s = 0; for (k = 0; k < K; k++) s += A[ck + k]; if (s < 1e-6) continue;
      for (var q = 0; q < cls.length; q++) { k = cls[q]; var r = d[k] / 2 * S; tot += Math.min(8, A[ck + k] / s * cover * area / (Math.PI * r * r)); }
    }
    var nbo = 0; if (sim.fixedC) for (c = 0; c < sim.N; c++) if (sim.fixedC[c]) nbo++;
    var RM = Math.max(0, T3.RMAX - nbo);
    var thin = tot > RM ? RM / tot : 1, n = 0, arr = R.instanceMatrix.array, col = R.instanceColor.array, RV = T3.RV;
    var X0 = sim.Lx / 2, Z0 = sim.Ly / 2;
    for (j = 0; j < ny && n < RM; j++) for (i = 0; i < nx && n < RM; i++) {
      c = j * nx + i; if (kind[c] === 2) continue;
      var ck2 = c * K, s2 = 0; for (k = 0; k < K; k++) s2 += A[ck2 + k]; if (s2 < 1e-6) continue;
      for (var q2 = 0; q2 < cls.length; q2++) {
        k = cls[q2]; var rr = d[k] / 2 * S, ne = Math.min(8, A[ck2 + k] / s2 * cover * area / (Math.PI * rr * rr)) * thin;
        var cnt = Math.floor(ne + hash(c * 31 + k * 7));
        for (var m = 0; m < cnt && n < RM; m++) {
          var key = c * 53 + k * 11 + m * 3, ux = hash(key + 0.13), uz = hash(key + 0.71);
          var fx = i + ux - 0.5, fy = j + uz - 0.5; if (fx < 0) fx = 0; if (fy < 0) fy = 0; if (fx > nx - 1) fx = nx - 1; if (fy > ny - 1) fy = ny - 1;
          var i0 = Math.min(nx - 2, fx | 0), j0 = Math.min(ny - 2, fy | 0), tx = fx - i0, ty = fy - j0, b00 = (j0 * nx + i0) * 3;
          var y = (tp[b00 + 1] * (1 - tx) + tp[b00 + 4] * tx) * (1 - ty) + (tp[b00 + nx * 3 + 1] * (1 - tx) + tp[b00 + nx * 3 + 4] * tx) * ty;
          var sz = rr * (0.7 + 0.6 * hash(key + 0.37)), v = RV[(key * 7) & 63], e = v.e, o = n * 16;
          arr[o] = e[0] * sz; arr[o + 1] = e[1] * sz; arr[o + 2] = e[2] * sz; arr[o + 3] = 0;
          arr[o + 4] = e[4] * sz; arr[o + 5] = e[5] * sz; arr[o + 6] = e[6] * sz; arr[o + 7] = 0;
          arr[o + 8] = e[8] * sz; arr[o + 9] = e[9] * sz; arr[o + 10] = e[10] * sz; arr[o + 11] = 0;
          arr[o + 12] = (fx + 0.5) * dx - X0; arr[o + 13] = y + sz * v.sy * 0.35; arr[o + 14] = Z0 - (fy + 0.5) * dx; arr[o + 15] = 1;
          var tone = 0.5 + 0.38 * hash(key + 0.91), wet = h[c] > 0.01 ? 0.72 : 1, w3 = n * 3, hue = hash(key + 0.53) - 0.5;
          col[w3] = tone * (0.97 + 0.08 * hue) * wet; col[w3 + 1] = tone * 0.94 * wet; col[w3 + 2] = tone * (0.88 - 0.08 * hue) * wet;
          n++;
        }
      }
    }
    // 巨礫（動かない礫）：セルごとに大きな石を 1 つ
    if (nbo) for (c = 0; c < sim.N && n < T3.RMAX; c++) {
      if (!sim.fixedC[c]) continue;
      var bi = c % nx, bj = (c - bi) / nx, bo = (bj * nx + bi) * 3, bsz = dx * (0.55 + 0.15 * hash(c + 0.3)), bv = RV[(c * 13) & 63], be = bv.e, bo16 = n * 16;
      arr[bo16] = be[0] * bsz; arr[bo16 + 1] = be[1] * bsz * 0.8; arr[bo16 + 2] = be[2] * bsz; arr[bo16 + 3] = 0;
      arr[bo16 + 4] = be[4] * bsz; arr[bo16 + 5] = be[5] * bsz * 0.8; arr[bo16 + 6] = be[6] * bsz; arr[bo16 + 7] = 0;
      arr[bo16 + 8] = be[8] * bsz; arr[bo16 + 9] = be[9] * bsz * 0.8; arr[bo16 + 10] = be[10] * bsz; arr[bo16 + 11] = 0;
      arr[bo16 + 12] = (bi + 0.5) * dx - X0; arr[bo16 + 13] = tp[bo + 1] + bsz * 0.25; arr[bo16 + 14] = Z0 - (bj + 0.5) * dx; arr[bo16 + 15] = 1;
      var bt = 0.42 + 0.12 * hash(c + 0.8); col[n * 3] = bt; col[n * 3 + 1] = bt * 0.97; col[n * 3 + 2] = bt * 0.92;
      n++;
    }
    R.count = n; R.instanceMatrix.needsUpdate = true; R.instanceColor.needsUpdate = true;
  }
  function render3D() {
    if (!T3 || !T3.terr) return; var cv = $('c3d'), w = cv.clientWidth, hgt = Math.round(Math.max(260, Math.min(520, w * 0.46)));
    if (cv._w !== w || cv._h !== hgt) { cv._w = w; cv._h = hgt; T3.r.setSize(w, hgt, false); cv.style.height = hgt + 'px'; T3.cam.aspect = w / hgt; T3.cam.updateProjectionMatrix(); }
    // 目標の視点へ滑らかに近づける
    if (T3.gdist) {
      var now = performance.now(), dtf = Math.min(0.1, ((now - (T3.lastF || now)) / 1000) || 0.016), a = 1 - Math.exp(-dtf * 14); T3.lastF = now;
      var dth = T3.gth - T3.th; dth = Math.atan2(Math.sin(dth), Math.cos(dth));
      T3.th += dth * a; T3.ph += (T3.gph - T3.ph) * a; T3.dist *= Math.exp(Math.log(T3.gdist / T3.dist) * a); T3.tgt.lerp(T3.gtgt, a);
    }
    var d = T3.dist, cp = Math.cos(T3.ph);
    T3.cam.position.set(T3.tgt.x + d * cp * Math.cos(T3.th), T3.tgt.y + d * Math.sin(T3.ph), T3.tgt.z + d * cp * Math.sin(T3.th)); T3.cam.lookAt(T3.tgt);
    var near = Math.max(0.05, d * 0.004); if (Math.abs(T3.cam.near - near) > near * 0.2) { T3.cam.near = near; T3.cam.far = d * 30 + 4000; T3.cam.updateProjectionMatrix(); }
    T3.wmat.uniforms.uTime.value = (performance.now() - T3.t0) / 1000;
    T3.r.render(T3.scene, T3.cam);
  }

  // ------------------------------------------------------------ 実行制御
  var spd = 0, lastViz = 0, lastCh = 0, dirty = true;
  function refresh(all) {
    if (!sim) return; computeRange($('v_mode').value);
    update3D(); drawPlan(); if (all) { redrawCharts(); drawStats(); snapInfo(); }
  }
  function apply() {
    running = false; $('run').textContent = '▶ 計算開始';
    try { P = readForm(); sim = new S2.Sim(P); } catch (e) {
      var em = '地形を切り替えられませんでした：' + e.message + '　表示はいまも前の地形のままです。';
      $('gridinfo').innerHTML = ''; var sp = document.createElement('span'); sp.style.color = 'var(--err, #b3261e)'; sp.style.fontWeight = '600'; sp.textContent = em; $('gridinfo').appendChild(sp);
      if ($('g_src').value === 'csv') $('csv_info').textContent = e.message; if ($('g_src').value === 'survey') $('sv_info').textContent = e.message;
      if (sim) refresh(true); return;
    }
    DP = sim.demMode && window.SaboDEMUI && window.SaboDEMUI.paths ? window.SaboDEMUI.paths() : null; PI = 0; PS = null; xsS = null; pathOptions();
    if (pathMode()) { var PS0 = pathSamples(); xsS = Math.round(PS0.L * 0.55); $('xs_x').max = PS0.L; $('xs_x').step = sim.dx; $('xs_x').value = xsS; }
    else { $('xs_x').max = sim.Lx; $('xs_x').step = sim.dx; if (xsX == null || xsX > sim.Lx) xsX = Math.round(sim.Lx * 0.55); $('xs_x').value = xsX; }
    $('gridinfo').textContent = '格子 ' + sim.nx + ' × ' + sim.ny + '（' + sim.N.toLocaleString() + ' セル、幅 ' + sim.dx + ' m）　粒径階 ' + sim.K + '　堰堤 ' + sim.dams.length + ' 基' + (sim.fromCSV && !sim.demMode ? '　CSV 地形（平均勾配 1/' + (1 / sim.I).toFixed(0) + '）' : '') + (sim.I > 0.1 + 1e-9 && sim.P.transport !== 'debris' ? '　注意：勾配 1/' + (1 / sim.I).toFixed(0) + ' は掃流砂の式の適用範囲（1/10 程度より緩い）を超えます。' : '') + (sim.warn.length ? '　注意：' + sim.warn.join(' ') : '');
    if (P.prefill > 0) $('gridinfo').textContent += '　はじめから満砂（比 ' + P.prefill + '）';
    if (sim.demMode) {
      var dmS = window.SaboDEMUI.state.dom;
      $('gridinfo').textContent += '　DEM の流域：' + ['回転なし（出口は東向き）', '時計回りに 90° 回して', '180° 回して', '反時計回りに 90° 回して'][dmS ? dmS.rot : 0] + '、出口を右端に置きました。供給点 ' + (P.sources || []).length + '、雨の横流入 ' + (P.latIn ? P.latIn.pts.length + ' セル' : (P.rain && P.rain.on ? '格子の雨' : 'なし')) + '、斜面崩壊 ' + sim.slides.length + '、天然ダム ' + sim.ldamList.length + '。';
    }
    // 文献の式の適用範囲から外れる設定の注意
    var nag = (P.sec && (P.sec.bed || P.sec.mom)) || (P.bank && P.bank.collapse === 'hasegawa' && P.bank.erodH > 0);
    if (nag && sim.I > 0.01 + 1e-9) $('gridinfo').textContent += '　注意：長田ら（1996）の2次流・崩落モデルの検証は勾配 1/100〜1/466 の砂の水路です。勾配 1/' + (1 / sim.I).toFixed(0) + ' では定性的な目安として見てください。';
    if (P.transport === 'noneq') $('gridinfo').textContent += '　侵食・堆積速度式の係数は水路実験の値です。';
    damRangeNote();
    var src = $('g_src').value, tkey = src + '|' + (src === 'csv' ? textHash($('csv_text').value) + $('csv_cs').value + $('csv_tr').checked + $('csv_flip').checked + $('csv_I').value
      : src === 'survey' ? textHash($('sv_xs').value + '#' + $('sv_long').value) + $('sv_align').value + $('sv_wall').value
      : JSON.stringify([P.L, P.B, P.slope, P.slopeTable, P.crossSlope, P.widthTable, P.bend, P.valley, P.seed, P.sideGrad])) + '|' + sim.dx;
    if (apply.tkey && apply.tkey !== tkey && ref) { ref = null; $('saveref').textContent = '比較用に保存'; }   // 別の地形の「比較」を重ねない
    apply.tkey = tkey; $('file_note').hidden = true;
    var dimKey = sim.Lx + 'x' + sim.Ly;
    build3D(); if (!apply.done || apply.dim !== dimKey) { home3D(true); apply.done = true; } apply.dim = dimKey;
    refresh(true);
  }
  function loop(now) {
    requestAnimationFrame(loop);
    if (!sim) return;
    if (running && sim.phase !== 'done') {
      var t0 = performance.now(), n = 0;
      while (performance.now() - t0 < 18 && sim.phase !== 'done') { sim.step(); n++; }
      xsRecord();
      spd = spd * 0.9 + 0.1 * (n / ((performance.now() - t0) / 1000));
      dirty = true;
      if (sim.phase === 'done') { running = false; $('run').textContent = '計算終了'; refresh(true); }
    }
    if (dirty && now - lastViz > 140) { lastViz = now; computeRange($('v_mode').value); update3D(); drawPlan(); }
    if (dirty && now - lastCh > 450) { lastCh = now; redrawCharts(); drawStats(); snapInfo(); dirty = running; }
    render3D();
  }
  $('run').onclick = function () {
    if (!sim || sim.phase === 'done') return;
    running = !running; this.textContent = running ? '⏸ 一時停止' : '▶ 再開'; dirty = true;
  };
  $('reset').onclick = apply; $('apply').onclick = apply; $('saveref').onclick = saveRef;
  $('dem_json').onclick = function () {
    try { var pj = readForm(); saveText('sabo2d_params.json', JSON.stringify(pj), 'dem_msg', true); }
    catch (e) { $('dem_msg').hidden = false; $('dem_msg').textContent = e.message; }
  };
  $('h_mode').onchange = modeSwitch;
  $('s_tcm').onchange = function () { if (this.value === 'iwagaki') $('s_hid').value = 'egiazaroff'; tcSwitch(); };
  // CSV を読む：UTF-8（BOM の有無を問わない）で読めなければ Shift_JIS（Excel の「CSV（コンマ区切り）」）として読み直す
  function readTextFile(f, ok, ng) {
    var rd = new FileReader();
    rd.onload = function () {
      var buf = new Uint8Array(rd.result), t;
      try { t = new TextDecoder('utf-8', { fatal: true }).decode(buf); }
      catch (e) { try { t = new TextDecoder('shift_jis').decode(buf); } catch (e2) { t = new TextDecoder('utf-8').decode(buf); } }
      ok(t.replace(/^\uFEFF/, ''));
    };
    rd.onerror = ng; rd.readAsArrayBuffer(f);
  }
  $('g_src').onchange = function () {
    geoSwitch(); if ($('g_src').value === 'csv') { $('csv_msg').hidden = true; csvInfo(); matchDx(); } if ($('g_src').value === 'survey') surveyInfo();
    if ($('g_src').value !== 'param' && $('d_fill').checked) { $('d_fill').checked = false; fillSwitch(); }   // モデル渓流の「満砂」を読み込んだ地形に持ち込まない
    damRangeNote();
  };
  // 堰堤の位置が、これから使う地形の長さを超えていないか
  function damRangeNote() {
    var L = curL(), out = [];
    Array.prototype.forEach.call($('damtab').tBodies[0].rows, function (tr) { var x = parseFloat(tr.querySelector('input').value); if (isFinite(x) && (x < 0 || x > L)) out.push(x); });
    $('dam_note').hidden = !out.length;
    $('dam_note').textContent = out.length ? '注意：x = ' + out.join('、') + ' m の堰堤は、地形の範囲（0〜' + Math.round(L) + ' m）の外にあるので置きません。' : '';
  }
  ['sv_xs', 'sv_long'].forEach(function (id) {
    var t = 0; $(id).addEventListener('input', function () { clearTimeout(t); t = setTimeout(surveyInfo, 300); });
    $(id + '_file').onchange = function () { var f = this.files && this.files[0], inp = this; if (!f) return; readTextFile(f, function (t) { routeFile(t, f.name, id === 'sv_xs' ? 'xs' : 'long'); inp.value = ''; }, function () { $('sv_info').textContent = 'ファイルを読めませんでした。'; }); };
  });
  ['sv_align', 'sv_wall', 'c_dx'].forEach(function (id) { $(id).addEventListener('change', function () { if ($('g_src').value === 'survey') surveyInfo(); }); });
  $('sv_clr').onclick = function () { $('sv_xs').value = ''; $('sv_long').value = ''; FILES.xs = FILES.long = ''; surveyInfo(); };
  $('sv_test').onclick = function () { FILES.xs = FILES.long = ''; $('sv_xs').value = TEST_XS; $('sv_long').value = TEST_LONG; $('sv_align').value = 'thalweg'; setDams([]); surveyInfo(); };
  $('sv_grid').onclick = function () {
    var r = surveyInfo(); if (r.err) return;
    $('csv_text').value = r.rows.map(function (row) { return row.map(function (v) { return v.toFixed(3); }).join(','); }).join('\n');
    $('csv_cs').value = r.cs; $('csv_tr').checked = false; $('csv_flip').checked = false; csvCache.key = '';
    $('sv_info').textContent += '　格子にした地形を「格子の CSV を読み込む」の欄に入れました。そこからコピーできます。';
  };
  $('pz_in').onclick = function () { zoomBy('prof', 0.6); redrawCharts(); };
  $('pz_out').onclick = function () { zoomBy('prof', 1 / 0.6); redrawCharts(); };
  $('pz_all').onclick = function () { var z = ZOOM.prof; if (z) { z.x0 = z.f0; z.x1 = z.f1; } redrawCharts(); };
  $('b_mode').onchange = geoSwitch;
  $('b_add').onclick = function () { $('btab').tBodies[0].appendChild(bRow([Math.round(num('g_L', 400) / 2), 0])); };
  $('b_ex').onclick = function () { var L = num('g_L', 400); setB([[0, 0], [Math.round(L * 0.2), 0], [Math.round(L * 0.4), 25], [Math.round(L * 0.6), -25], [Math.round(L * 0.8), 0], [L, 0]]); };
  $('csv_file').onchange = function () {
    var f = this.files && this.files[0]; if (!f) return;
    var inp = this; readTextFile(f, function (t) { routeFile(t, f.name, 'grid'); inp.value = ''; }, function () { $('csv_info').textContent = 'ファイルを読めませんでした。'; });
  };
  var csvT = 0; ['csv_text', 'csv_cs', 'csv_tr', 'csv_flip'].forEach(function (id) { $(id).addEventListener('input', function () { clearTimeout(csvT); csvT = setTimeout(csvInfo, 300); }); $(id).addEventListener('change', csvInfo); });
  // 計算後の河床を「格子の CSV」欄に入れ、堰堤をいまの位置・向き・天端標高で固定する
  function bedToCSV() {
    if (!sim) return null;
    var B = bedCSV(), nd = fixDams();
    $('csv_text').value = B.text; $('csv_cs').value = sim.dx; $('csv_tr').checked = false; $('csv_flip').checked = false;
    $('csv_fh').value = B.floorH; $('csv_I').value = ''; csvCache.key = '';
    if ($('d_fill').checked) { $('d_fill').checked = false; fillSwitch(); }
    var opts = Array.prototype.map.call($('c_dx').options, function (o) { return +o.value; }), k = opts.indexOf(sim.dx); if (k >= 0) $('c_dx').selectedIndex = k;
    return { nd: nd, floorH: B.floorH };
  }
  $('csv_cont').onclick = function () {
    var r = bedToCSV(); if (!r) return;
    $('g_src').value = 'csv'; geoSwitch(); apply();
    $('csv_msg').hidden = false;
    $('csv_msg').textContent = 'いまの河床から計算をやり直せるようにしました。堰堤 ' + r.nd + ' 基は、前と同じ位置・向き・天端標高で置いています（表の「向き・中心・天端標高・堤長」）。粒度は初期の河床材料に戻ります。「▶ 計算開始」で続きを計算します。';
  };
  $('csv_now').onclick = function () {
    var r = bedToCSV(); if (!r) return;
    var msg = 'いまの河床を「格子の CSV」欄に入れ、堰堤 ' + r.nd + ' 基を同じ位置・向き・天端標高で表に固定しました。「地形の作り方」を「格子の CSV を読み込む」にして「適用して最初から」を押すと、この河床から計算できます。粒度は初期の河床材料に戻ります。';
    $('csv_info').textContent = msg; $('csv_msg').textContent = msg; $('csv_msg').hidden = false;
  };
  $('csv_save').onclick = function () {
    if (!sim) return;
    saveText('河床_' + Math.round(sim.t / 60) + '分.csv', bedCSV().text);
  };
  $('dam_fix').onclick = function () { var n = fixDams(); $('dam_note').hidden = false; $('dam_note').style.color = ''; $('dam_note').textContent = n ? 'いまの堰堤 ' + n + ' 基の位置・向き・天端標高・堤長を表に書き込みました。値を変えて「適用して最初から」で動かせます。' : '堰堤がありません。'; };
  $('d_adv').onchange = advSwitch;
  $('csv_copy').onclick = function () {
    var t = $('csv_text'); if (!t.value) { $('csv_now').click(); }
    var done = function () { $('csv_copy').textContent = 'コピーしました'; setTimeout(function () { $('csv_copy').textContent = 'CSV をコピー'; }, 1500); };
    try { navigator.clipboard.writeText(t.value).then(done, function () { t.hidden = false; t.select(); }); } catch (e) { t.select(); }
  };
  $('xs_x').oninput = function () { if (pathMode()) xsS = +this.value; else xsX = +this.value; drawXsec(); drawPlan(); drawProfile(); };
  $('l_path').onchange = function () { PI = +this.value; PS = null; var S0 = pathSamples(); xsS = Math.round(S0.L * 0.55); $('xs_x').max = S0.L; $('xs_x').value = xsS; drawXsec(); drawPlan(); drawProfile(); };
  $('dam_add').onclick = function () { var L = num('g_L', 400); $('damtab').tBodies[0].appendChild(damRow({ x: Math.round(L / 2), H: 1.2, b: 0 })); };
  $('w_add').onclick = function () { var L = num('g_L', 400), B = num('g_B', 50), n = $('wtab').tBodies[0].rows.length; $('wtab').tBodies[0].appendChild(wRow(n ? [Math.round(L / 2), B] : [0, B])); };
  $('w_ex').onclick = function () { var L = num('g_L', 400), B = num('g_B', 50); setW([[0, B], [Math.round(L * 0.35), B], [Math.round(L * 0.5), Math.round(B * 0.4)], [Math.round(L * 0.6), Math.round(B * 0.4)], [Math.round(L * 0.75), Math.round(B * 1.4)], [L, Math.round(B * 1.4)]]); };
  $('w_clr').onclick = function () { setW([]); };
  $('i_add').onclick = function () { var L = num('g_L', 400), I = num('g_I', 20), n = $('itab').tBodies[0].rows.length; $('itab').tBodies[0].appendChild(iRow([Math.round(L * (n ? 0.6 : 0.4) / 10) * 10, Math.round(I * (n ? 1.5 : 2))])); };
  $('i_ex1').onclick = function () { var L = num('g_L', 400), I = num('g_I', 20); setI([[Math.round(L * 0.45 / 10) * 10, Math.round(I * 3)]]); };
  $('i_ex2').onclick = function () { var L = num('g_L', 400), I = num('g_I', 20); setI([[Math.round(L * 0.3 / 10) * 10, Math.round(I * 2)], [Math.round(L * 0.65 / 10) * 10, Math.round(I * 4)]]); };
  $('i_clr').onclick = function () { setI([]); };
  $('t_add').onclick = function () { var L = num('g_L', 400), B = num('g_B', 50); $('ttab').tBodies[0].appendChild(tRow({ x: Math.round(L / 2 / 10) * 10, side: 'L', B: Math.max(4, Math.round(B / 5)), I: 6, r: 0.3, rs: 1, lag: 0 })); };
  $('t_ex').onclick = function () { var L = num('g_L', 400), B = num('g_B', 50), db = $('h_mode').value === 'debris'; setT([{ x: Math.round(L * 0.45 / 10) * 10, side: 'L', B: Math.max(4, Math.round(B / 5)), I: db ? 3 : 5, r: 0.4, rs: db ? 0.4 : 1.5, lag: db ? 60 : 0, ang: 60 }]); };
  $('t_clr').onclick = function () { setT([]); };
  // 左の条件設定欄の幅：右端をドラッグで変える（ダブルクリックで元に戻す）。幅はこのブラウザに覚える
  (function () {
    var app = document.querySelector('.app'), side = $('side'), grip = $('side_grip'), KEY = 'sabo2d.sideW';
    if (!app || !grip) return;
    function setW(w) { w = Math.max(280, Math.min(w, Math.max(320, window.innerWidth - 420))); app.style.setProperty('--sw', w + 'px'); return w; }
    try { var w0 = parseFloat(localStorage.getItem(KEY)); if (w0 > 0) setW(w0); } catch (e) {}
    var x0 = 0, wS = 0, wNow = 0;
    function mv(e) { wNow = setW(wS + e.clientX - x0); }
    function up() { document.removeEventListener('pointermove', mv); document.removeEventListener('pointerup', up); document.body.classList.remove('resizing'); try { localStorage.setItem(KEY, String(wNow)); } catch (e) {} window.dispatchEvent(new Event('resize')); }
    grip.addEventListener('pointerdown', function (e) { e.preventDefault(); x0 = e.clientX; wS = wNow = side.getBoundingClientRect().width; document.body.classList.add('resizing'); document.addEventListener('pointermove', mv); document.addEventListener('pointerup', up); });
    grip.addEventListener('dblclick', function () { app.style.removeProperty('--sw'); try { localStorage.removeItem(KEY); } catch (e) {} window.dispatchEvent(new Event('resize')); });
    grip.addEventListener('keydown', function (e) { if (e.key !== 'ArrowLeft' && e.key !== 'ArrowRight') return; e.preventDefault(); wNow = setW(side.getBoundingClientRect().width + (e.key === 'ArrowRight' ? 30 : -30)); try { localStorage.setItem(KEY, String(wNow)); } catch (er) {} window.dispatchEvent(new Event('resize')); });
  })();
  $('gw_on').onchange = gwSwitch;
  $('gr_add').onclick = function () { if ($('grtab').tBodies[0].rows.length < 6) $('grtab').tBodies[0].appendChild(grRow({ d: 0.3, bed: 0, feed: 0 })); grSummary(); };
  $('gr_pre').onchange = function () { var v = GR_PRE[this.value]; if (!v) return; $('gr_dm').value = v[0]; $('gr_sg').value = v[1]; $('gr_k').value = v[2]; $('gr_fr').value = v[3]; $('gr_gen').onclick(); };
  $('gr_gen').onclick = function () {
    var g = makeGrains(+$('gr_dm').value || 5, +$('gr_sg').value || 1, +$('gr_k').value || 5, +$('gr_fr').value || 1), tb = $('grtab').tBodies[0]; tb.innerHTML = '';
    // 割合は 0.1 % 単位で、合計がちょうど 100 % になるよう最大の階で調整
    ['bed', 'feed'].forEach(function (key) {
      var r = g.map(function (q) { return Math.round(q[key] * 1000) / 1000; }), sm = r.reduce(function (a, b) { return a + b; }, 0), im = 0;
      r.forEach(function (v, k) { if (v > r[im]) im = k; }); r[im] = Math.round((r[im] + 1 - sm) * 1000) / 1000;
      g.forEach(function (q, k) { q[key] = r[k]; });
    });
    g.forEach(function (q) { tb.appendChild(grRow({ d: +(q.d).toPrecision(3), bed: q.bed, feed: q.feed })); });
    grSummary();
  };
  $('grtab').closest('details').addEventListener('toggle', grSummary);
  $('pre_low').onclick = function () { var L = curL(), a = []; for (var x = 60; x < L - 50; x += 45) a.push({ x: x, H: 1.2, b: 0 }); setDams(a); };
  $('pre_high').onclick = function () { var L = curL(), B = num('g_B', 50); setDams([{ x: Math.round(L * 0.82), H: 8, b: Math.round(B * 0.4) }]); };
  $('pre_none').onclick = function () { setDams([]); };
  ['A', 'B', 'S'].forEach(function (k) { $('pre_' + k).onclick = function () { applyPreset(k); apply(); }; });
  $('b_ex2').onclick = function () { var L = num('g_L', 400), B = num('g_B', 50); setB([[0, 0, B], [Math.round(L * 0.25), 15, B], [Math.round(L * 0.42), 0, Math.round(B * 0.45)], [Math.round(L * 0.5), -10, Math.round(B * 0.45)], [Math.round(L * 0.7), -20, Math.round(B * 1.6)], [L, 0, B]]); };
  $('f_add').onclick = function () { var L = num('g_L', 400); $('ftab').tBodies[0].appendChild(fRow([Math.round(L / 2), 1, 30, 3])); };
  $('f_clr').onclick = function () { setF([]); };
  $('d_fill').onchange = fillSwitch;
  $('s_tr').onchange = trSwitch;
  $('damtab').addEventListener('input', function () { presetDams = false; damRangeNote(); }); $('damtab').addEventListener('click', function (e) { if (e.target.tagName === 'BUTTON') presetDams = false; setTimeout(damRangeNote, 0); });
  ['pre_low', 'pre_high', 'pre_none', 'dam_add'].forEach(function (id) { $(id).addEventListener('click', function () { presetDams = false; setTimeout(damRangeNote, 0); }); });
  $('gen_go').addEventListener('click', function () { presetDams = false; });
  $('gen_go').onclick = function () { var a = [], n = num('gen_n', 7, 1, 60), x0 = num('gen_x0', 60), s = num('gen_s', 45, 4), hh = num('gen_h', 1.2, 0.1), b = num('gen_b', 0, 0); for (var i = 0; i < n; i++) a.push({ x: x0 + s * i, H: hh, b: b }); setDams(a); };
  ['v_vex', 'v_turb', 'p_core', 'v_cont', 'v_rock', 'v_tree'].forEach(function (id) { $(id).onchange = function () { refresh(false); }; });
  $('v_mode').onchange = function () { var m = this.value; if (m === 'bed' || m === 'bedrel') $('v_cont').checked = true; refresh(false); };
  ['l_det', 'l_bed', 'l_var'].forEach(function (id) { $(id).onchange = redrawCharts; });
  window.addEventListener('resize', function () { dirty = true; refresh(true); });
  if (window.matchMedia) window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', function () { refresh(true); });

  fillForm(); advSwitch(); applyPreset('A'); geoSwitch(); init3D(); apply();
  requestAnimationFrame(loop);
  window.__sabo = { get sim() { return sim; }, start: function () { if (!running) $('run').click(); }, view: function (th, ph, dist, tx, tz) { if (!T3) return; T3.th = T3.gth = th; T3.ph = T3.gph = ph; T3.dist = T3.gdist = dist; T3.tgt.set(tx || 0, 0, tz || 0); T3.gtgt.copy(T3.tgt); }, get T3() { return T3; } };
})();
