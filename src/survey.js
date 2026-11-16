/*
 * survey.js — 横断・縦断の測量値から格子地形（行列）を作る。tools/xs2grid.py と同じ方法。
 *   横断：x, s, z（断面位置、左岸からの距離、標高）
 *   縦断：x, z（最深河床高、省略可）
 * 断面の間は、最深部の左右それぞれを幅で正規化して比例配分で補間する。
 * 戻り値の rows[0] が左岸側、左の列が上流。
 */
(function (root) {
  'use strict';
  function readTable(text, ncol) {
    var rows = [];
    String(text || '').split(/\r?\n/).forEach(function (ln) {
      ln = ln.trim(); if (!ln || ln.charAt(0) === '#') return;
      var v = ln.split(/[,\t; ]+/).filter(function (t) { return t !== ''; }).map(Number);
      if (v.length < ncol || v.slice(0, ncol).some(function (x) { return !isFinite(x); })) return;
      rows.push(v.slice(0, ncol));
    });
    return rows;
  }
  function interp1(xs, ys, x) {
    var n = xs.length; if (x <= xs[0]) return ys[0]; if (x >= xs[n - 1]) return ys[n - 1];
    var lo = 0, hi = n - 1; while (hi - lo > 1) { var m = (lo + hi) >> 1; if (xs[m] <= x) lo = m; else hi = m; }
    var t = (x - xs[lo]) / (xs[hi] - xs[lo]); return ys[lo] * (1 - t) + ys[hi] * t;
  }
  function buildSections(rows, align) {
    var by = {}, keys = [];
    rows.forEach(function (r) { var k = Math.round(r[0] * 1e6) / 1e6; if (!by[k]) { by[k] = []; keys.push(k); } by[k].push([r[1], r[2]]); });
    keys.sort(function (a, b) { return a - b; });
    var secs = [];
    for (var q = 0; q < keys.length; q++) {
      var pts = by[keys[q]].sort(function (a, b) { return a[0] - b[0]; });
      if (pts.length < 3) return { err: 'x = ' + keys[q] + ' m の断面の点が 3 点未満です。' };
      var s = pts.map(function (p) { return p[0]; }), z = pts.map(function (p) { return p[1]; });
      var zmin = Math.min.apply(null, z), ss = 0, ns = 0;
      for (var i = 0; i < z.length; i++) if (z[i] <= zmin + 1e-9) { ss += s[i]; ns++; }
      secs.push({ x: keys[q], s: s, zr: z.map(function (v) { return v - zmin; }), zmin: zmin, st: ss / ns, sl: s[0], sr: s[s.length - 1] });
    }
    if (secs.length < 2) return { err: '断面が 2 つ以上必要です（いまは ' + secs.length + '）。' };
    var U;
    if (align === 'left') secs.forEach(function (c) { c.off = -c.sl; });
    else if (align === 'center') { U = Math.max.apply(null, secs.map(function (c) { return (c.sr - c.sl) / 2; })); secs.forEach(function (c) { c.off = U - (c.sl + c.sr) / 2; }); }
    else { U = Math.max.apply(null, secs.map(function (c) { return c.st - c.sl; })); secs.forEach(function (c) { c.off = U - c.st; }); }
    secs.forEach(function (c) { c.uL = c.sl + c.off; c.uT = c.st + c.off; c.uR = c.sr + c.off; });
    return { secs: secs };
  }
  function secZ(c, side, e) {
    var s = side < 0 ? c.sl + (c.st - c.sl) * e : c.st + (c.sr - c.st) * e;
    return interp1(c.s, c.zr, s);
  }
  function makeGrid(secs, longp, dx, wall) {
    var xs = secs.map(function (c) { return c.x; }), lx = longp ? longp.map(function (p) { return p[0]; }) : null, lz = longp ? longp.map(function (p) { return p[1]; }) : null;
    var all = xs.concat(lx || []), xa = Math.min.apply(null, all), xb = Math.max.apply(null, all);
    var L = xb - xa, W = Math.max.apply(null, secs.map(function (c) { return c.uR; }));
    var nx = Math.max(10, Math.round(L / dx)), ny = Math.max(3, Math.round(W / dx) + 1), rows = [], i, j;
    for (j = 0; j < ny; j++) rows.push(new Array(nx));
    for (i = 0; i < nx; i++) {
      var x = xa + (i + 0.5) * dx, a, b, t;
      if (x <= xs[0]) { a = b = secs[0]; t = 0; }
      else if (x >= xs[xs.length - 1]) { a = b = secs[secs.length - 1]; t = 0; }
      else { var k = 1; while (xs[k] < x) k++; a = secs[k - 1]; b = secs[k]; t = (x - a.x) / (b.x - a.x); }
      var uL = a.uL * (1 - t) + b.uL * t, uT = a.uT * (1 - t) + b.uT * t, uR = a.uR * (1 - t) + b.uR * t;
      var zt = longp ? interp1(lx, lz, x) : a.zmin * (1 - t) + b.zmin * t;
      var zL = secZ(a, -1, 0) * (1 - t) + secZ(b, -1, 0) * t, zR = secZ(a, 1, 1) * (1 - t) + secZ(b, 1, 1) * t;
      for (j = 0; j < ny; j++) {
        var u = (j + 0.5) * dx, zr, e;
        if (u < uL) zr = zL + wall * (uL - u);
        else if (u > uR) zr = zR + wall * (u - uR);
        else if (u <= uT) { e = uT > uL ? (u - uL) / (uT - uL) : 1; zr = secZ(a, -1, e) * (1 - t) + secZ(b, -1, e) * t; }
        else { e = uR > uT ? (u - uT) / (uR - uT) : 0; zr = secZ(a, 1, e) * (1 - t) + secZ(b, 1, e) * t; }
        rows[j][i] = zt + zr;
      }
    }
    return { rows: rows, L: L, W: W, nx: nx, ny: ny };
  }
  // 横断を与えないときの断面（谷壁＋中央へ傾く谷底）
  function templateRows(B, cross, side, Hw, x0, x1) {
    var wl = Math.max(4, Hw / Math.max(0.2, side)), out = [];
    [x0, x1].forEach(function (x) {
      var hb = cross * B / 2;
      out.push([x, 0, hb + Hw], [x, wl / 2, hb + Hw / 2], [x, wl, hb], [x, wl + B / 4, hb / 2], [x, wl + B / 2, 0], [x, wl + 3 * B / 4, hb / 2], [x, wl + B, hb], [x, wl + B + wl / 2, hb + Hw / 2], [x, wl + B + wl, hb + Hw]);
    });
    return out;
  }
  // 入口：テキストから格子地形を作る
  function toRows(o) {
    var xsRows = readTable(o.xsText, 3), longp = o.longText ? readTable(o.longText, 2) : [];
    longp.sort(function (a, b) { return a[0] - b[0]; });
    if (!longp.length) longp = null;
    var usedTemplate = false;
    if (!xsRows.length) {
      if (!longp || longp.length < 2) return { err: '横断の値か、縦断の値（2 点以上）が必要です。' };
      xsRows = templateRows(o.B || 50, o.cross || 0.02, o.side || 0.8, o.Hw || 6, longp[0][0], longp[longp.length - 1][0]);
      usedTemplate = true;
    }
    var bs = buildSections(xsRows, o.align || 'thalweg'); if (bs.err) return bs;
    var g = makeGrid(bs.secs, longp, o.dx || 2, o.wall == null ? 1 : o.wall);
    g.cs = o.dx || 2; g.nsec = bs.secs.length; g.usedTemplate = usedTemplate; g.hasLong = !!longp;
    return g;
  }
  var api = { toRows: toRows, readTable: readTable };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.SaboSurvey = api;
})(typeof window !== 'undefined' ? window : globalThis);
