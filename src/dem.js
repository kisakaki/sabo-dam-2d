/* DEM の読み込みと流域の解析（ブラウザと Node の両方で動く）
 *  - 読み込み：ESRI ASCII grid（.asc）、XYZ 点群（.xyz/.txt/.csv）、基盤地図情報 数値標高モデル（JPGIS GML .xml、複数可）、GeoTIFF（非圧縮・LZW・Deflate）
 *  - 窪地の埋め：Priority-Flood + ε（Barnes et al. 2014）
 *  - 流向：D8（最急勾配の 8 方向）、集水面積、谷筋（集水面積がしきい値以上）、谷頭（上流に谷筋がない谷筋のセル）
 *  - 流域：DEM の端（または欠測）へ出る出口ごとの集水域。谷筋の任意の点を出口にした小流域
 *  - 計算領域：流域を切り出し、出口が右（下流端）に来るように 90° 単位で回し、出口から右端まで水路を掘る
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.SaboDEM = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';
  // D8 の 8 方向（東から反時計回り）。行は北が 0（下に行くほど南）
  var DC = [1, 1, 0, -1, -1, -1, 0, 1], DR = [0, -1, -1, -1, 0, 1, 1, 1];
  var DL = [1, Math.SQRT2, 1, Math.SQRT2, 1, Math.SQRT2, 1, Math.SQRT2];

  function mkDEM(nx, ny, csx, csy, x0, yTop, z, src) {
    return { nx: nx, ny: ny, csx: csx, csy: csy, cs: csx, x0: x0, yTop: yTop, z: z, src: src || '' };
  }

  // ---------------------------------------------------------------- 読み込み
  function parseASC(text) {
    var lines = text.split(/\r?\n/), H = {}, i = 0;
    for (; i < lines.length; i++) {
      var m = lines[i].trim().match(/^([A-Za-z_]+)\s+([-+0-9.eE]+)$/);
      if (!m) break;
      H[m[1].toLowerCase()] = parseFloat(m[2]);
    }
    var nx = H.ncols | 0, ny = H.nrows | 0;
    if (!(nx > 1 && ny > 1)) throw new Error('ESRI ASCII grid の見出し（ncols, nrows）が読めません');
    var csx = H.cellsize || H.dx, csy = H.cellsize || H.dy || csx;
    var x0 = isFinite(H.xllcorner) ? H.xllcorner : (H.xllcenter || 0) - csx / 2;
    var yll = isFinite(H.yllcorner) ? H.yllcorner : (H.yllcenter || 0) - csy / 2;
    var nd = isFinite(H.nodata_value) ? H.nodata_value : -9999;
    var z = new Float64Array(nx * ny), k = 0;
    for (; i < lines.length && k < nx * ny; i++) {
      var t = lines[i].trim(); if (!t) continue;
      var a = t.split(/[\s,]+/);
      for (var q = 0; q < a.length && k < nx * ny; q++) { var v = parseFloat(a[q]); z[k++] = (!isFinite(v) || v === nd || v <= -9998) ? NaN : v; }
    }
    if (k < nx * ny) throw new Error('ESRI ASCII grid の値が足りません（' + k + ' / ' + nx * ny + '）');
    return mkDEM(nx, ny, csx, csy, x0, yll + ny * csy, z, 'asc');
  }

  function parseXYZ(text) {
    var pts = [], lines = text.split(/\r?\n/);
    for (var i = 0; i < lines.length; i++) {
      var a = lines[i].trim().split(/[\s,;]+/); if (a.length < 3) continue;
      var x = parseFloat(a[0]), y = parseFloat(a[1]), zz = parseFloat(a[2]);
      if (isFinite(x) && isFinite(y) && isFinite(zz)) pts.push([x, y, zz]);
    }
    if (pts.length < 9) throw new Error('XYZ の点が足りません');
    function step(vals) {
      vals.sort(function (p, q) { return p - q; });
      var u = [vals[0]];
      for (var i2 = 1; i2 < vals.length; i2++) if (vals[i2] - u[u.length - 1] > 1e-6) u.push(vals[i2]);
      var ds = []; for (i2 = 1; i2 < u.length; i2++) ds.push(u[i2] - u[i2 - 1]);
      ds.sort(function (p, q) { return p - q; });
      return { min: u[0], max: u[u.length - 1], d: ds.length ? ds[Math.floor(ds.length / 2)] : 1 };
    }
    var sx = step(pts.map(function (p) { return p[0]; })), sy = step(pts.map(function (p) { return p[1]; }));
    var cs = Math.min(sx.d, sy.d);
    var nx = Math.round((sx.max - sx.min) / cs) + 1, ny = Math.round((sy.max - sy.min) / cs) + 1;
    if (nx * ny > 4e7) throw new Error('XYZ の格子が大きすぎます（' + nx + '×' + ny + '）');
    var z = new Float64Array(nx * ny).fill(NaN);
    for (i = 0; i < pts.length; i++) {
      var c = Math.round((pts[i][0] - sx.min) / cs), r = Math.round((sy.max - pts[i][1]) / cs);
      if (c >= 0 && c < nx && r >= 0 && r < ny) z[r * nx + c] = pts[i][2];
    }
    return mkDEM(nx, ny, cs, cs, sx.min - cs / 2, sy.max + cs / 2, z, 'xyz');
  }

  // 基盤地図情報 数値標高モデル（JPGIS GML）。緯度経度の格子。複数のタイルをつなげる
  function parseGSIXMLOne(text) {
    function tag(name) { var m = text.match(new RegExp('<gml:' + name + '[^>]*>([^<]*)</gml:' + name + '>')); return m ? m[1].trim() : null; }
    var lo = tag('lowerCorner'), up = tag('upperCorner'), hi = tag('high'), sp = tag('startPoint');
    if (!lo || !up || !hi) throw new Error('基盤地図情報の XML（lowerCorner, upperCorner, high）が読めません');
    lo = lo.split(/\s+/).map(Number); up = up.split(/\s+/).map(Number); hi = hi.split(/\s+/).map(Number);
    var nx = hi[0] + 1, ny = hi[1] + 1, st = sp ? sp.split(/\s+/).map(Number) : [0, 0];
    var m = text.match(/<gml:tupleList>([\s\S]*?)<\/gml:tupleList>/);
    var z = new Float64Array(nx * ny).fill(NaN), k = st[1] * nx + st[0];
    var rows = m ? m[1].split(/\r?\n/) : [];   // tupleList が無いタイル（全部欠測）は欠測のまま
    for (var i = 0; i < rows.length && k < nx * ny; i++) {
      var t = rows[i].trim(); if (!t) continue;
      var p = t.lastIndexOf(','), v = parseFloat(p >= 0 ? t.slice(p + 1) : t);
      z[k++] = (!isFinite(v) || v <= -9998) ? NaN : v;
    }
    return { lat0: lo[0], lon0: lo[1], lat1: up[0], lon1: up[1], nx: nx, ny: ny, z: z };
  }
  function parseGSIXML(texts) {
    if (!Array.isArray(texts)) texts = [texts];
    // 読めないタイルは飛ばす（数を D.skipped に記録）
    var T = [], bad = 0, msg = '';
    texts.forEach(function (t) { try { T.push(parseGSIXMLOne(t)); } catch (e) { bad++; msg = e.message; } });
    if (!T.length) throw new Error(msg || '基盤地図情報の XML が読めません');
    // いちばん細かいタイルの間隔に合わせる（5m と 10m が混じるときなど）
    var dlat = Infinity, dlon = Infinity;
    T.forEach(function (t) { dlat = Math.min(dlat, (t.lat1 - t.lat0) / t.ny); dlon = Math.min(dlon, (t.lon1 - t.lon0) / t.nx); });
    var la0 = Infinity, la1 = -Infinity, lo0 = Infinity, lo1 = -Infinity;
    T.forEach(function (t) { la0 = Math.min(la0, t.lat0); la1 = Math.max(la1, t.lat1); lo0 = Math.min(lo0, t.lon0); lo1 = Math.max(lo1, t.lon1); });
    var NX = Math.round((lo1 - lo0) / dlon), NY = Math.round((la1 - la0) / dlat);
    if (NX * NY > 4e7) throw new Error('つなげた格子が大きすぎます（' + NX + '×' + NY + '）。範囲が離れたタイルが混じっていないか確かめてください');
    var Z = new Float64Array(NX * NY).fill(NaN);
    // 細かいタイルを先に置き、粗いタイルは欠測の所だけ埋める
    T.sort(function (a, b) { return (a.lat1 - a.lat0) / a.ny - (b.lat1 - b.lat0) / b.ny; });
    T.forEach(function (t) {
      var tl = (t.lat1 - t.lat0) / t.ny, tn = (t.lon1 - t.lon0) / t.nx, same = Math.abs(tl - dlat) < 1e-3 * dlat && Math.abs(tn - dlon) < 1e-3 * dlon;
      var R0 = Math.round((la1 - t.lat1) / dlat), C0 = Math.round((t.lon0 - lo0) / dlon);
      var R1 = Math.round((la1 - t.lat0) / dlat), C1 = Math.round((t.lon1 - lo0) / dlon);
      for (var R = Math.max(0, R0); R < Math.min(NY, R1); R++) for (var C = Math.max(0, C0); C < Math.min(NX, C1); C++) {
        var k = R * NX + C; if (!same && isFinite(Z[k])) continue;
        var r = same ? R - R0 : Math.floor((la1 - (R + 0.5) * dlat - t.lat0) / tl), c = same ? C - C0 : Math.floor((lo0 + (C + 0.5) * dlon - t.lon0) / tn);
        r = same ? r : t.ny - 1 - r;
        if (r < 0 || r >= t.ny || c < 0 || c >= t.nx) continue;
        var v = t.z[r * t.nx + c]; if (isFinite(v)) Z[k] = v;
      }
    });
    // 緯度経度 → 局所の平面（m）。正方形の格子に直す（小さい方の間隔）
    var latm = (la0 + la1) / 2 * Math.PI / 180;
    var my = 111132.92 - 559.82 * Math.cos(2 * latm) + 1.175 * Math.cos(4 * latm);
    var mx = 111412.84 * Math.cos(latm) - 93.5 * Math.cos(3 * latm);
    var csx = dlon * mx, csy = dlat * my;
    var D = mkDEM(NX, NY, csx, csy, lo0 * mx, la1 * my, Z, 'gsi');
    D.geo = { lat0: la0, lat1: la1, lon0: lo0, lon1: lo1 };
    var E = squareCells(D); E.skipped = bad; return E;
  }

  // 欠測の穴埋め：DEM の縁につながらない欠測（内側の穴）を、外側から 1 周ずつ、周り 8 セルの値の平均で埋める（同時に更新）。
  // 基盤地図情報の 5m メッシュは測量のない所が欠測になり、谷筋が途切れて水が流れなくなるため。縁につながる欠測は DEM の外として残す
  function fillHoles(D) {
    var nx = D.nx, ny = D.ny, N = nx * ny, z = D.z, ext = new Uint8Array(N), st = [], c, i, j, n = 0;
    for (c = 0; c < N; c++) { i = c % nx; j = (c - i) / nx; if (z[c] !== z[c] && (i === 0 || j === 0 || i === nx - 1 || j === ny - 1)) { ext[c] = 1; st.push(c); } }
    while (st.length) {
      c = st.pop(); i = c % nx; j = (c - i) / nx;
      for (var dj = -1; dj <= 1; dj++) for (var di = -1; di <= 1; di++) {
        var ii = i + di, jj = j + dj; if (ii < 0 || jj < 0 || ii >= nx || jj >= ny) continue;
        var o = jj * nx + ii; if (!ext[o] && z[o] !== z[o]) { ext[o] = 1; st.push(o); }
      }
    }
    var hole = []; for (c = 0; c < N; c++) if (z[c] !== z[c] && !ext[c]) hole.push(c);
    if (!hole.length) return 0;
    var Z = Float64Array.from(z), left = hole, it = 0;
    while (left.length && it++ < 5000) {
      var put = [], vals = [], rest = [];
      for (var q = 0; q < left.length; q++) {
        c = left[q]; i = c % nx; j = (c - i) / nx; var s = 0, m = 0;
        for (var dj2 = -1; dj2 <= 1; dj2++) for (var di2 = -1; di2 <= 1; di2++) {
          if (!di2 && !dj2) continue; var i2 = i + di2, j2 = j + dj2; if (i2 < 0 || j2 < 0 || i2 >= nx || j2 >= ny) continue;
          var v = Z[j2 * nx + i2]; if (v === v) { s += v; m++; }
        }
        if (m) { put.push(c); vals.push(s / m); } else rest.push(c);
      }
      if (!put.length) break;
      for (q = 0; q < put.length; q++) Z[put[q]] = vals[q];
      n += put.length; left = rest;
    }
    D.z = Z; D.filledHoles = n; return n;
  }
  // 平滑化：幅 w の箱形の平均を縦・横に 3 回（ガウス分布で標準偏差 sigma [m] に近い）。欠測は除いて平均し、欠測のまま残す
  function smooth(D, sigma) {
    var nx = D.nx, ny = D.ny, sc = sigma / D.cs;
    if (!(sc > 0.3)) return D;
    var rr = Math.max(1, Math.round((Math.sqrt(4 * sc * sc + 1) - 1) / 2)), z = Float64Array.from(D.z), t = new Float64Array(nx * ny);
    var bs = new Float64Array(Math.max(nx, ny) + 1), bn = new Float64Array(Math.max(nx, ny) + 1);
    function pass(src, dst, n, m, stride, step) {   // n 本の線、それぞれ m 個。線の先頭は stride*k、点の間隔 step
      for (var k = 0; k < n; k++) {
        var o = k * stride, i;
        bs[0] = 0; bn[0] = 0;
        for (i = 0; i < m; i++) { var v = src[o + i * step], f = v === v; bs[i + 1] = bs[i] + (f ? v : 0); bn[i + 1] = bn[i] + (f ? 1 : 0); }
        for (i = 0; i < m; i++) {
          var v0 = src[o + i * step]; if (v0 !== v0) { dst[o + i * step] = NaN; continue; }
          var a = Math.max(0, i - rr), b = Math.min(m, i + rr + 1), cnt = bn[b] - bn[a];
          dst[o + i * step] = cnt > 0 ? (bs[b] - bs[a]) / cnt : v0;
        }
      }
    }
    for (var it = 0; it < 3; it++) { pass(z, t, ny, nx, nx, 1); pass(t, z, nx, ny, 1, nx); }
    var E = mkDEM(nx, ny, D.cs, D.cs, D.x0, D.yTop, z, D.src); E.cs = D.cs; E.geo = D.geo; E.smooth = sigma; return E;
  }
  // k × k セルを平均して粗くする（大きすぎる DEM 用。欠測は除いて平均）
  function coarsen(D, k) {
    var nx = Math.floor(D.nx / k), ny = Math.floor(D.ny / k), z = new Float64Array(nx * ny);
    for (var r = 0; r < ny; r++) for (var c = 0; c < nx; c++) {
      var s = 0, n = 0;
      for (var a = 0; a < k; a++) for (var b = 0; b < k; b++) { var v = D.z[(r * k + a) * D.nx + c * k + b]; if (isFinite(v)) { s += v; n++; } }
      z[r * nx + c] = n ? s / n : NaN;
    }
    var E = mkDEM(nx, ny, D.cs * k, D.cs * k, D.x0, D.yTop, z, D.src); E.cs = D.cs * k; E.geo = D.geo; return E;
  }
  // 格子の間隔が縦横で違うときは、小さい方の間隔の正方形格子に直す（双一次補間）
  function squareCells(D) {
    if (Math.abs(D.csx - D.csy) < 1e-6 * D.csx) { D.cs = D.csx; return D; }
    var cs = Math.min(D.csx, D.csy), W = D.nx * D.csx, Hh = D.ny * D.csy;
    var nx = Math.max(2, Math.round(W / cs)), ny = Math.max(2, Math.round(Hh / cs));
    var z = new Float64Array(nx * ny);
    for (var r = 0; r < ny; r++) for (var c = 0; c < nx; c++) z[r * nx + c] = sampleBil(D, (c + 0.5) * cs / D.csx - 0.5, (r + 0.5) * cs / D.csy - 0.5);
    var E = mkDEM(nx, ny, cs, cs, D.x0, D.yTop, z, D.src); E.geo = D.geo; return E;
  }
  function sampleBil(D, fx, fy) {
    var nx = D.nx, ny = D.ny;
    if (fx < 0) fx = 0; if (fy < 0) fy = 0; if (fx > nx - 1) fx = nx - 1; if (fy > ny - 1) fy = ny - 1;
    var c0 = Math.min(nx - 2, Math.floor(fx)), r0 = Math.min(ny - 2, Math.floor(fy)), tx = fx - c0, ty = fy - r0;
    var a = D.z[r0 * nx + c0], b = D.z[r0 * nx + c0 + 1], c = D.z[(r0 + 1) * nx + c0], d = D.z[(r0 + 1) * nx + c0 + 1];
    if (!(isFinite(a) && isFinite(b) && isFinite(c) && isFinite(d))) {   // 欠測が混じるときは最も近い点
      var v = D.z[Math.round(fy) * nx + Math.round(fx)]; return v;
    }
    return (a * (1 - tx) + b * tx) * (1 - ty) + (c * (1 - tx) + d * tx) * ty;
  }

  // ---- GeoTIFF（最低限）：ストリップとタイル、非圧縮・LZW・Deflate、予測子 2・3、整数と浮動小数点
  // inflate(Uint8Array) は Deflate（zlib）の展開関数。ブラウザでは非同期（DecompressionStream）なので、parseTIFF は Promise を返す
  function parseTIFF(buf, inflate) {
    var dv = new DataView(buf instanceof ArrayBuffer ? buf : buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength));
    var u8 = new Uint8Array(dv.buffer);
    var le = dv.getUint16(0) === 0x4949;
    if (dv.getUint16(2, le) !== 42) throw new Error('GeoTIFF ではありません（BigTIFF は未対応）');
    var off = dv.getUint32(4, le), n = dv.getUint16(off, le), T = {};
    var SZ = { 1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 6: 1, 7: 1, 8: 2, 9: 4, 10: 8, 11: 4, 12: 8, 16: 8 };
    function rd(type, p) {
      switch (type) {
        case 1: case 7: return dv.getUint8(p); case 6: return dv.getInt8(p);
        case 3: return dv.getUint16(p, le); case 8: return dv.getInt16(p, le);
        case 4: return dv.getUint32(p, le); case 9: return dv.getInt32(p, le);
        case 5: return dv.getUint32(p, le) / dv.getUint32(p + 4, le); case 10: return dv.getInt32(p, le) / dv.getInt32(p + 4, le);
        case 11: return dv.getFloat32(p, le); case 12: return dv.getFloat64(p, le);
        case 16: return dv.getUint32(p, le) + dv.getUint32(p + 4, le) * 4294967296;
      }
      return 0;
    }
    for (var i = 0; i < n; i++) {
      var e = off + 2 + i * 12, tg = dv.getUint16(e, le), ty = dv.getUint16(e + 2, le), cnt = dv.getUint32(e + 4, le);
      var sz = (SZ[ty] || 1) * cnt, p = sz <= 4 ? e + 8 : dv.getUint32(e + 8, le);
      if (ty === 2) { var s = ''; for (var q = 0; q < cnt; q++) { var ch = u8[p + q]; if (ch) s += String.fromCharCode(ch); } T[tg] = s; }
      else { var arr = []; for (q = 0; q < cnt; q++) arr.push(rd(ty, p + q * (SZ[ty] || 1))); T[tg] = arr; }
    }
    var W = T[256][0], H = T[257][0], bps = (T[258] || [8])[0], cmp = (T[259] || [1])[0], spp = (T[277] || [1])[0];
    var pred = (T[317] || [1])[0], fmt = (T[339] || [1])[0], bpp = bps / 8;
    if (spp !== 1) throw new Error('GeoTIFF：1 バンドの標高だけに対応しています（' + spp + ' バンド）');
    var tiled = !!T[322], tw = tiled ? T[322][0] : W, th = tiled ? T[323][0] : (T[278] ? T[278][0] : H);
    var offs = tiled ? T[324] : T[273], cnts = tiled ? T[325] : T[279];
    var nd = T[42113] ? parseFloat(T[42113]) : NaN;
    var scale = T[33550] || [1, 1], tie = T[33922] || [0, 0, 0, 0, 0, 0];
    var z = new Float64Array(W * H);
    var nTx = Math.ceil(W / tw);
    function readVal(b, p) {
      var d = new DataView(b.buffer, b.byteOffset, b.byteLength);
      if (fmt === 3) return bps === 64 ? d.getFloat64(p, le) : d.getFloat32(p, le);
      if (fmt === 2) return bps === 8 ? d.getInt8(p) : bps === 16 ? d.getInt16(p, le) : d.getInt32(p, le);
      return bps === 8 ? d.getUint8(p) : bps === 16 ? d.getUint16(p, le) : d.getUint32(p, le);
    }
    function place(k, raw) {
      var rowsIn = tiled ? th : Math.min(th, H - k * th);
      var b = raw;
      if (pred === 2) {          // 横方向の差分（整数）
        var d = new DataView(b.buffer, b.byteOffset, b.byteLength);
        for (var r = 0; r < rowsIn; r++) for (var x = 1; x < tw; x++) {
          var p0 = (r * tw + x) * bpp, p1 = p0 - bpp;
          if (bps === 8) b[p0] = (b[p0] + b[p1]) & 255;
          else if (bps === 16) d.setUint16(p0, (d.getUint16(p0, le) + d.getUint16(p1, le)) & 65535, le);
          else d.setUint32(p0, (d.getUint32(p0, le) + d.getUint32(p1, le)) >>> 0, le);
        }
      } else if (pred === 3) {   // 浮動小数点の予測子（バイトの並べ替え＋差分）
        var rowB = tw * bpp, out = new Uint8Array(b.length);
        for (r = 0; r < rowsIn; r++) {
          var o = r * rowB;
          for (x = 1; x < rowB; x++) b[o + x] = (b[o + x] + b[o + x - 1]) & 255;
          for (x = 0; x < tw; x++) for (var bb = 0; bb < bpp; bb++) {
            var src = b[o + bb * tw + x];
            out[o + x * bpp + (le ? bpp - 1 - bb : bb)] = src;
          }
        }
        b = out;
      }
      var tx0 = tiled ? (k % nTx) * tw : 0, ty0 = tiled ? Math.floor(k / nTx) * th : k * th;
      for (var rr = 0; rr < rowsIn; rr++) {
        var Y = ty0 + rr; if (Y >= H) break;
        for (var xx = 0; xx < tw; xx++) {
          var X = tx0 + xx; if (X >= W) break;
          var v = readVal(b, (rr * tw + xx) * bpp);
          z[Y * W + X] = (!isFinite(v) || v === nd || v <= -9998 || v >= 1e30) ? NaN : v;
        }
      }
    }
    var jobs = [];
    for (var k = 0; k < offs.length; k++) {
      var seg = u8.subarray(offs[k], offs[k] + cnts[k]);
      if (cmp === 1) jobs.push(new Uint8Array(seg));
      else if (cmp === 5) jobs.push(lzw(seg));
      else if (cmp === 8 || cmp === 32946) { if (!inflate) throw new Error('Deflate の展開ができません'); jobs.push(inflate(seg)); }
      else throw new Error('GeoTIFF の圧縮形式（' + cmp + '）には対応していません。非圧縮・LZW・Deflate にしてください');
    }
    function finish(raws) {
      for (var k2 = 0; k2 < raws.length; k2++) place(k2, raws[k2]);
      var csx = scale[0], csy = scale[1] || scale[0];
      var x0 = tie[3] - tie[0] * csx, yTop = tie[4] + tie[1] * csy;
      return squareCells(mkDEM(W, H, csx, csy, x0, yTop, z, 'tif'));
    }
    if (jobs.some(function (j) { return j && typeof j.then === 'function'; })) return Promise.all(jobs).then(finish);
    return finish(jobs);
  }
  function lzw(src) {
    var out = [], dict = [], bitPos = 0, width = 9, prev = null, i;
    function reset() { dict = []; for (i = 0; i < 256; i++) dict[i] = [i]; dict[256] = null; dict[257] = null; width = 9; }
    function next() {
      var v = 0;
      for (var b = 0; b < width; b++) { var byte = src[(bitPos + b) >> 3]; v = (v << 1) | ((byte >> (7 - ((bitPos + b) & 7))) & 1); }
      bitPos += width; return v;
    }
    reset();
    while ((bitPos + width) <= src.length * 8) {
      var code = next();
      if (code === 257) break;
      if (code === 256) { reset(); prev = null; continue; }
      var entry;
      if (code < dict.length && dict[code]) entry = dict[code];
      else if (prev) entry = prev.concat([prev[0]]);
      else break;
      for (i = 0; i < entry.length; i++) out.push(entry[i]);
      if (prev) dict.push(prev.concat([entry[0]]));
      prev = entry;
      if (dict.length + 1 >= (1 << width) && width < 12) width++;
    }
    return new Uint8Array(out);
  }

  // ---------------------------------------------------------------- 解析
  // 二分ヒープ（標高の小さい順）
  function Heap(N) { this.k = new Int32Array(N); this.v = new Float64Array(N); this.n = 0; }
  Heap.prototype.push = function (key, val) {
    var i = this.n++, k = this.k, v = this.v;
    while (i > 0) { var p = (i - 1) >> 1; if (v[p] <= val) break; k[i] = k[p]; v[i] = v[p]; i = p; }
    k[i] = key; v[i] = val;
  };
  Heap.prototype.pop = function () {
    var k = this.k, v = this.v, top = k[0], n = --this.n, lk = k[n], lv = v[n], i = 0;
    while (true) {
      var l = 2 * i + 1; if (l >= n) break;
      var r = l + 1, m = (r < n && v[r] < v[l]) ? r : l;
      if (v[m] >= lv) break;
      k[i] = k[m]; v[i] = v[m]; i = m;
    }
    k[i] = lk; v[i] = lv; return top;
  };

  // 窪地の埋め（Priority-Flood + ε）。欠測（NaN）と DEM の端が出口
  function fill(D, epsIn) {
    var nx = D.nx, ny = D.ny, N = nx * ny, z = D.z, zf = new Float64Array(N), done = new Uint8Array(N);
    var H = new Heap(N), c, r, q, eps = epsIn > 0 ? epsIn : 1e-5;
    for (c = 0; c < N; c++) {
      zf[c] = z[c];
      if (!isFinite(z[c])) { done[c] = 1; continue; }
      var x = c % nx, y = (c / nx) | 0, edge = x === 0 || y === 0 || x === nx - 1 || y === ny - 1;
      if (!edge) for (q = 0; q < 8; q++) { var n2 = (y + DR[q]) * nx + x + DC[q]; if (!isFinite(z[n2])) { edge = true; break; } }
      if (edge) { done[c] = 1; H.push(c, z[c]); }
    }
    while (H.n) {
      c = H.pop(); var cx = c % nx, cy = (c / nx) | 0;
      for (q = 0; q < 8; q++) {
        var xx = cx + DC[q], yy = cy + DR[q]; if (xx < 0 || yy < 0 || xx >= nx || yy >= ny) continue;
        var o = yy * nx + xx; if (done[o]) continue;
        done[o] = 1;
        if (zf[o] <= zf[c]) zf[o] = zf[c] + eps * DL[q];
        H.push(o, zf[o]);
      }
    }
    return zf;
  }

  // D8 の流向。-1：DEM の端か欠測へ出る出口、-2：欠測
  function flowDir(D, zf) {
    var nx = D.nx, ny = D.ny, N = nx * ny, dir = new Int8Array(N), csx = D.cs;
    for (var c = 0; c < N; c++) {
      if (!isFinite(zf[c])) { dir[c] = -2; continue; }
      var x = c % nx, y = (c / nx) | 0, best = -1, bs = 0, out = false;
      for (var q = 0; q < 8; q++) {
        var xx = x + DC[q], yy = y + DR[q];
        if (xx < 0 || yy < 0 || xx >= nx || yy >= ny) { out = true; continue; }
        var o = yy * nx + xx; if (!isFinite(zf[o])) { out = true; continue; }
        var s = (zf[c] - zf[o]) / (DL[q] * csx);
        if (s > bs) { bs = s; best = q; }
      }
      dir[c] = best >= 0 ? best : (out ? -1 : -1);
    }
    return dir;
  }
  function downOf(D, dir, c) {
    var q = dir[c]; if (q < 0) return -1;
    var x = c % D.nx + DC[q], y = ((c / D.nx) | 0) + DR[q];
    return (x < 0 || y < 0 || x >= D.nx || y >= D.ny) ? -1 : y * D.nx + x;
  }
  // 標高の高い順に並べる（流れの上流から下流への順）
  function order(D, zf) {
    var N = D.nx * D.ny, idx = [], c;
    for (c = 0; c < N; c++) if (isFinite(zf[c])) idx.push(c);
    idx.sort(function (a, b) { return zf[b] - zf[a]; });
    return idx;
  }
  function accumulate(D, dir, ord) {
    var acc = new Float64Array(D.nx * D.ny), i;
    for (i = 0; i < ord.length; i++) acc[ord[i]] = 1;
    for (i = 0; i < ord.length; i++) { var c = ord[i], o = downOf(D, dir, c); if (o >= 0) acc[o] += acc[c]; }
    return acc;
  }

  // 全体の解析。thrArea：谷筋とみなす集水面積 [m2]
  function analyze(D, opt) {
    opt = opt || {};
    var nx = D.nx, ny = D.ny, N = nx * ny, a1 = D.cs * D.cs;
    var zf = fill(D), dir = flowDir(D, zf), ord = order(D, zf), acc = accumulate(D, dir, ord);
    var thr = Math.max(2, Math.round((opt.thrArea > 0 ? opt.thrArea : 20000) / a1));
    var stream = new Uint8Array(N), c, q;
    for (c = 0; c < N; c++) if (acc[c] >= thr && dir[c] !== -2) stream[c] = 1;
    // 谷頭：谷筋のうち、上流側に谷筋のセルがないもの
    var hasUp = new Uint8Array(N);
    for (c = 0; c < N; c++) if (stream[c]) { var o = downOf(D, dir, c); if (o >= 0) hasUp[o] = 1; }
    var heads = [];
    for (c = 0; c < N; c++) if (stream[c] && !hasUp[c]) heads.push(c);
    // 流域：出口（流向が -1 のセル）ごとに番号。下流から上流へ順に伝える
    var bid = new Int32Array(N).fill(-1), basins = [];
    for (var i = ord.length - 1; i >= 0; i--) {
      c = ord[i];
      var dn = downOf(D, dir, c);
      if (dn < 0) { bid[c] = basins.length; basins.push({ id: basins.length, outlet: c, n: 0 }); }
      else bid[c] = bid[dn];
    }
    for (c = 0; c < N; c++) if (bid[c] >= 0) basins[bid[c]].n++;
    basins.forEach(function (b) { b.area = b.n * a1; b.acc = acc[b.outlet]; });
    // Strahler 次数（谷筋だけ）
    var so = new Uint8Array(N), mx = new Uint8Array(N), cntMx = new Uint8Array(N);
    for (i = 0; i < ord.length; i++) {
      c = ord[i]; if (!stream[c]) continue;
      if (!so[c]) so[c] = mx[c] > 0 ? (cntMx[c] >= 2 ? mx[c] + 1 : mx[c]) : 1;
      var d2 = downOf(D, dir, c);
      if (d2 >= 0 && stream[d2]) { if (so[c] > mx[d2]) { mx[d2] = so[c]; cntMx[d2] = 1; } else if (so[c] === mx[d2]) cntMx[d2]++; }
    }
    return { zf: zf, dir: dir, ord: ord, acc: acc, stream: stream, heads: heads, bid: bid, basins: basins, thr: thr, strahler: so };
  }

  // 出口 c0 の上流の範囲（流域）
  function upstream(D, A, c0) {
    var N = D.nx * D.ny, m = new Uint8Array(N); m[c0] = 1;
    var o = A.ord, idx = -1, i;
    for (i = o.length - 1; i >= 0; i--) if (o[i] === c0) { idx = i; break; }
    for (i = idx - 1; i >= 0; i--) { var c = o[i], dn = downOf(D, A.dir, c); if (dn >= 0 && m[dn]) m[c] = 1; }
    return m;
  }

  // 谷筋のセルの局所勾配（流向の方向、下流側 n セルの平均）
  function localSlope(D, A, c, n) {
    var z0 = A.zf[c], L = 0, cc = c, k;
    for (k = 0; k < (n || 3); k++) { var dn = downOf(D, A.dir, cc); if (dn < 0) break; L += DL[A.dir[cc]] * D.cs; cc = dn; }
    return L > 0 ? Math.max(1e-3, (z0 - A.zf[cc]) / L) : 1e-2;
  }

  // 斜面の長さ（谷筋までの流下距離）と、谷筋のセルごとの直接の集水面積（斜面から直接入る分）
  function hillslopes(D, A, mask) {
    var N = D.nx * D.ny, dist = new Float64Array(N), recv = new Int32Array(N).fill(-1), i, c;
    // 下流から上流へ：谷筋なら距離 0・受け手は自分、そうでなければ下流の距離＋1 区間
    for (i = A.ord.length - 1; i >= 0; i--) {
      c = A.ord[i]; if (mask && !mask[c]) continue;
      if (A.stream[c]) { dist[c] = 0; recv[c] = c; continue; }
      var dn = downOf(D, A.dir, c);
      if (dn < 0 || (mask && !mask[dn])) { dist[c] = 0; recv[c] = -1; continue; }
      dist[c] = dist[dn] + DL[A.dir[c]] * D.cs; recv[c] = recv[dn];
    }
    return { dist: dist, recv: recv };
  }

  // ---------------------------------------------------------------- 計算領域の用意
  // basinMask の範囲を切り出し、出口の流向が右（+x）になるよう回す。出口から右端まで水路を掘る。
  // 返り値：rows（北が上の 2 次元配列、ソルバーの格子 CSV と同じ並び）、basin、chan（谷筋）、セルの対応
  function prepareDomain(D, A, mask, outlet, opt) {
    opt = opt || {};
    var nx = D.nx, ny = D.ny, c, x, y;
    // 計算に使う標高：窪地を埋めた標高（既定）。谷筋が斜めにつながる所は、横か縦の隣のセルも谷筋にして
    // 4 近傍でつながるようにする（浅水流の計算は辺を通してしか水が流れないため）
    var zs = opt.fill === false ? D.z : A.zf, st = new Uint8Array(A.stream), zz2 = new Float64Array(zs);
    for (c = 0; c < nx * ny; c++) {
      if (!st[c] || !mask[c]) continue;
      var q0 = A.dir[c]; if (q0 < 0 || !(q0 % 2)) continue;
      var dn = downOf(D, A.dir, c); if (dn < 0) continue;
      var cx0 = c % nx, cy0 = (c / nx) | 0, a1 = cy0 * nx + cx0 + DC[q0], b1 = (cy0 + DR[q0]) * nx + cx0;
      var pick = zz2[a1] <= zz2[b1] ? a1 : b1;
      var zm = 0.5 * (zz2[c] + zz2[dn]);
      if (zz2[pick] > zm) zz2[pick] = zm;
      if (mask[pick]) st[pick] = 1;
    }
    var x0 = nx, x1 = -1, y0 = ny, y1 = -1;
    for (c = 0; c < nx * ny; c++) if (mask[c]) { x = c % nx; y = (c / nx) | 0; if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y; }
    var pad = 2; x0 = Math.max(0, x0 - pad); y0 = Math.max(0, y0 - pad); x1 = Math.min(nx - 1, x1 + pad); y1 = Math.min(ny - 1, y1 + pad);
    var W = x1 - x0 + 1, H = y1 - y0 + 1;
    // 出口の流向 → 回転（時計回りに rot × 90°）
    var q = A.dir[outlet], ox = DC[q < 0 ? 0 : q], oy = DR[q < 0 ? 0 : q];
    if (q < 0 || downOf(D, A.dir, outlet) < 0) {   // DEM の端（または欠測）へ出る出口：いちばん近い端の外向き（斜めに出る端のセルも）
      var cx = outlet % nx, cy = (outlet / nx) | 0, dl = [nx - 1 - cx, cy, cx, ny - 1 - cy], km = 0;
      for (var k0 = 1; k0 < 4; k0++) if (dl[k0] < dl[km]) km = k0;
      ox = [1, 0, -1, 0][km]; oy = [0, -1, 0, 1][km];
    }
    var rot;   // 東→0、北→1、西→2、南→3（時計回りの回数）
    if (Math.abs(ox) >= Math.abs(oy)) rot = ox > 0 ? 0 : 2; else rot = oy < 0 ? 1 : 3;
    // 元の (col, row) → 回転後 (i, j)（j は北が 0）
    var RW = (rot % 2) ? H : W, RH = (rot % 2) ? W : H;
    function fwd(col, row) {
      var a = col - x0, b = row - y0;
      if (rot === 0) return [a, b];
      if (rot === 1) return [H - 1 - b, a];        // 時計回り 90°：北が東へ
      if (rot === 2) return [W - 1 - a, H - 1 - b];
      return [b, W - 1 - a];
    }
    function inv(i, j) {
      if (rot === 0) return [x0 + i, y0 + j];
      if (rot === 1) return [x0 + j, y0 + H - 1 - i];
      if (rot === 2) return [x0 + W - 1 - i, y0 + H - 1 - j];
      return [x0 + W - 1 - j, y0 + i];
    }
    // 右に足す列：3 列か、計算格子の 4 セル分の長い方（出口の先に緩衝の水路をとり、その先で水を抜くため）
    var nAdd = Math.max(3, Math.ceil(4 * ((opt && opt.dx) || D.cs) / D.cs));
    var po = fwd(outlet % nx, (outlet / nx) | 0), extra = RW - 1 - po[0] + nAdd;   // 出口から右端まで
    var NW = RW + nAdd, rows = [], bas = [], chn = [], src = [];
    var zmaxB = -Infinity;
    for (c = 0; c < nx * ny; c++) if (mask[c] && D.z[c] > zmaxB) zmaxB = D.z[c];
    for (var j = 0; j < RH; j++) {
      var r = [], b = [], s = [], m = [];
      for (var i = 0; i < NW; i++) {
        if (i < RW) {
          var o = inv(i, j), oc = o[1] * nx + o[0];
          var zz = zz2[oc]; if (!isFinite(zz)) zz = zmaxB;
          r.push(zz); b.push(mask[oc] ? 1 : 0); s.push(mask[oc] && st[oc] ? 1 : 0); m.push(oc);
        } else { r.push(r[RW - 1]); b.push(0); s.push(0); m.push(-1); }   // 右に足す列は端の列の標高のまま（流域の外）
      }
      rows.push(r); bas.push(b); chn.push(s); src.push(m);
    }
    // 出口から右端まで水路（出口の勾配で下げる）
    var Iout = Math.max(0.01, localSlope(D, A, outlet, 4)), zo = zz2[outlet], cor = [];
    for (j = 0; j < RH; j++) { var cr = []; for (var q0 = 0; q0 < NW; q0++) cr.push(0); cor.push(cr); }
    // 幅は 5 セル（中央が谷筋）。1 セル幅だと、格子を粗くしたとき両側の山と混ざって通路がふさがるため
    for (i = po[0] + 1; i < NW; i++) {
      var dd = (i - po[0]) * D.cs, zc = zo - Iout * dd;
      for (var dj = -2; dj <= 2; dj++) {
        var jj = po[1] + dj; if (jj < 0 || jj >= RH) continue;
        if (dj && bas[jj][i] && rows[jj][i] < zc) continue;
        if (!bas[jj][i]) cor[jj][i] = i - po[0];   // 出口の先の水路（流域の外だったセル）。値は出口からの列数
        rows[jj][i] = zc; bas[jj][i] = 1; if (!dj) chn[jj][i] = 1;
      }
    }
    return {
      rows: rows, cs: D.cs, basin: bas, chan: chn, corr: cor, cellOf: src, outXY: [(po[0] + 0.5) * D.cs, (RH - po[1] - 0.5) * D.cs], rot: rot, x0: x0, y0: y0, W: W, H: H, nx: NW, ny: RH,
      outlet: po, Iout: Iout, zTop: zmaxB, extra: extra,
      // DEM のセル → 計算領域の座標（m、x は右、y は上が正。ソルバーの座標）
      toXY: function (cell) { var p = fwd(cell % nx, (cell / nx) | 0); return [(p[0] + 0.5) * D.cs, (RH - p[1] - 0.5) * D.cs]; },
      // DEM の流向の単位ベクトル（計算領域の向き）
      dirXY: function (cell) {
        var qq = A.dir[cell]; if (qq < 0) return [1, 0];
        var vx = DC[qq], vy = -DR[qq], L = Math.sqrt(vx * vx + vy * vy); vx /= L; vy /= L;   // vy：北が正
        for (var k = 0; k < rot; k++) { var t = vx; vx = vy; vy = -t; }   // 時計回りに 90°
        return [vx, vy];
      }
    };
  }

  // ---------------------------------------------------------------- 例の DEM（架空の山地流域）
  // 乱数の起伏に、流れの集中で谷を刻む簡単な地形発達（stream power の侵食）を繰り返して谷を作る
  function synthDEM(opt) {
    opt = opt || {};
    var nx = opt.nx || 160, ny = opt.ny || 120, cs = opt.cs || 10, seed = opt.seed || 3, N = nx * ny;
    var s = seed >>> 0; function rnd() { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; }
    var g = [], gn = 9; for (var q = 0; q < (gn + 2) * (gn + 2) * 4; q++) g.push(rnd() * 2 - 1);
    function noise(x, y, o) {
      var xi = Math.floor(x), yi = Math.floor(y), tx = x - xi, ty = y - yi;
      function G(a, b) { return g[(((a % (gn + 1)) + (gn + 1)) % (gn + 1) + (((b % (gn + 1)) + (gn + 1)) % (gn + 1)) * (gn + 2) + o * (gn + 2) * (gn + 2)) % g.length]; }
      var sx = tx * tx * (3 - 2 * tx), sy = ty * ty * (3 - 2 * ty);
      return (G(xi, yi) * (1 - sx) + G(xi + 1, yi) * sx) * (1 - sy) + (G(xi, yi + 1) * (1 - sx) + G(xi + 1, yi + 1) * sx) * sy;
    }
    var z = new Float64Array(N);
    for (var r = 0; r < ny; r++) for (var c = 0; c < nx; c++) {
      var X = c / nx, Y = r / ny;
      // 東（右）へ下がる斜面＋南北の中央がやや低い＋起伏
      var base = 420 * (1 - X) + 90 * Math.pow(Math.abs(Y - 0.5) * 2, 1.5) + 40;
      var nz = 34 * noise(X * 4, Y * 4, 0) + 16 * noise(X * 8 + 3, Y * 8 + 1, 1) + 7 * noise(X * 16 + 5, Y * 16 + 7, 2);
      z[r * nx + c] = base + nz;
    }
    var D = mkDEM(nx, ny, cs, cs, 0, ny * cs, z, 'synth');
    // 地形発達：dz = -K A^0.5 S（谷を刻む）＋拡散（尾根を丸める）
    for (var it = 0; it < (opt.iter || 40); it++) {
      var zf = fill(D), dir = flowDir(D, zf), ord = order(D, zf), acc = accumulate(D, dir, ord), zn = new Float64Array(N);
      for (c = 0; c < N; c++) {
        var dn = downOf(D, dir, c), S = dn >= 0 ? Math.max(0, (zf[c] - zf[dn]) / (DL[dir[c]] * cs)) : 0;
        var x = c % nx, y = (c / nx) | 0, lap = 0, nn = 0;
        for (q = 0; q < 8; q += 2) { var xx = x + DC[q], yy = y + DR[q]; if (xx >= 0 && yy >= 0 && xx < nx && yy < ny) { lap += z[yy * nx + xx] - z[c]; nn++; } }
        zn[c] = z[c] - 0.9 * Math.sqrt(acc[c] * cs * cs) * S * 0.06 + 0.08 * lap;
      }
      for (c = 0; c < N; c++) { var xc = c % nx; z[c] = xc === nx - 1 ? z[c] : zn[c]; }
    }
    // 窪地は 1 % の勾配の斜面で埋める（例の地形に大きな平らな湖ができないように）
    var zf2 = fill(D, 0.01 * cs);
    for (c = 0; c < N; c++) z[c] = zf2[c];
    return D;
  }

  return {
    parseASC: parseASC, parseXYZ: parseXYZ, parseGSIXML: parseGSIXML, parseTIFF: parseTIFF, lzw: lzw,
    fill: fill, flowDir: flowDir, accumulate: accumulate, order: order, analyze: analyze, upstream: upstream,
    downOf: downOf, localSlope: localSlope, hillslopes: hillslopes, prepareDomain: prepareDomain, synthDEM: synthDEM,
    squareCells: squareCells, coarsen: coarsen, smooth: smooth, fillHoles: fillHoles, DC: DC, DR: DR, DL: DL
  };
});
