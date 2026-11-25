"""sabo2d_dem.py — DEM の読み込み・流域の解析・雨の流出（ブラウザ版 dem.js・hydro.js と同じ手順）

  読み込み：ESRI ASCII grid（.asc）、XYZ 点群、基盤地図情報 数値標高モデル（JPGIS GML .xml、複数可）、GeoTIFF（非圧縮・LZW・Deflate）
  解析    ：窪地の埋め（Priority-Flood + ε、Barnes et al. 2014）、D8 の流向、集水面積、谷筋、谷頭、流域、Strahler 次数
  計算領域：流域を切り出し、出口が右（下流端）に来るように 90° 単位で回し、出口から右端まで水路を掘る
  雨      ：損失（流出係数 / f1・Rsa / 浸透能）と、合理式（時間面積法・土研式の到達時間）か貯留関数法で、谷筋のセルへの横流入
"""
import heapq
import math
import re
import struct
import zlib

import numpy as np

DC = [1, 1, 0, -1, -1, -1, 0, 1]
DR = [0, -1, -1, -1, 0, 1, 1, 1]
DL = [1, math.sqrt(2), 1, math.sqrt(2), 1, math.sqrt(2), 1, math.sqrt(2)]


class DEM:
    def __init__(self, nx, ny, csx, csy, x0, yTop, z, src=''):
        self.nx, self.ny, self.csx, self.csy, self.cs = nx, ny, csx, csy, csx
        self.x0, self.yTop, self.z, self.src = x0, yTop, np.asarray(z, dtype=float).ravel(), src
        self.geo = None


# ---------------------------------------------------------------- 読み込み
def parse_asc(text):
    lines = text.splitlines(); H = {}; i = 0
    for i, ln in enumerate(lines):
        m = re.match(r'^\s*([A-Za-z_]+)\s+([-+0-9.eE]+)\s*$', ln)
        if not m:
            break
        H[m.group(1).lower()] = float(m.group(2))
    nx, ny = int(H['ncols']), int(H['nrows'])
    csx = H.get('cellsize') or H.get('dx'); csy = H.get('cellsize') or H.get('dy') or csx
    x0 = H['xllcorner'] if 'xllcorner' in H else H.get('xllcenter', 0.0) - csx / 2
    yll = H['yllcorner'] if 'yllcorner' in H else H.get('yllcenter', 0.0) - csy / 2
    nd = H.get('nodata_value', -9999.0)
    vals = np.array(' '.join(lines[i:]).replace(',', ' ').split(), dtype=float)[:nx * ny]
    if vals.size < nx * ny:
        raise ValueError('ESRI ASCII grid の値が足りません')
    vals[(vals == nd) | (vals <= -9998)] = np.nan
    return _square(DEM(nx, ny, csx, csy, x0, yll + ny * csy, vals, 'asc'))


def parse_xyz(text):
    a = []
    for ln in text.splitlines():
        p = re.split(r'[\s,;]+', ln.strip())
        if len(p) >= 3:
            try:
                a.append([float(p[0]), float(p[1]), float(p[2])])
            except ValueError:
                pass
    a = np.array(a)

    def step(v):
        u = np.unique(np.round(v, 6)); d = np.diff(u)
        return u[0], u[-1], float(np.median(d)) if d.size else 1.0
    x0, x1, dxs = step(a[:, 0]); y0, y1, dys = step(a[:, 1]); cs = min(dxs, dys)
    nx = int(round((x1 - x0) / cs)) + 1; ny = int(round((y1 - y0) / cs)) + 1
    z = np.full(nx * ny, np.nan)
    c = np.round((a[:, 0] - x0) / cs).astype(int); r = np.round((y1 - a[:, 1]) / cs).astype(int)
    ok = (c >= 0) & (c < nx) & (r >= 0) & (r < ny)
    z[r[ok] * nx + c[ok]] = a[ok, 2]
    return DEM(nx, ny, cs, cs, x0 - cs / 2, y1 + cs / 2, z, 'xyz')


def _gsi_one(text):
    def tag(n):
        m = re.search(r'<gml:' + n + r'[^>]*>([^<]*)</gml:' + n + '>', text)
        return m.group(1).strip() if m else None
    lo = [float(v) for v in tag('lowerCorner').split()]; up = [float(v) for v in tag('upperCorner').split()]
    hi = [int(v) for v in tag('high').split()]; sp = tag('startPoint'); st = [int(v) for v in sp.split()] if sp else [0, 0]
    nx, ny = hi[0] + 1, hi[1] + 1
    mb = re.search(r'<gml:tupleList>([\s\S]*?)</gml:tupleList>', text); body = mb.group(1) if mb else ''
    z = np.full(nx * ny, np.nan); k = st[1] * nx + st[0]
    for ln in body.splitlines():
        ln = ln.strip()
        if not ln:
            continue
        v = float(ln.rsplit(',', 1)[-1])
        if k < nx * ny:
            z[k] = np.nan if v <= -9998 else v
        k += 1
    return dict(lat0=lo[0], lon0=lo[1], lat1=up[0], lon1=up[1], nx=nx, ny=ny, z=z)


def parse_gsi_xml(texts):
    if isinstance(texts, str):
        texts = [texts]
    T = []; bad = 0; msg = ''
    for t in texts:   # 読めないタイルは飛ばす
        try:
            T.append(_gsi_one(t))
        except Exception as e:
            bad += 1; msg = str(e)
    if not T:
        raise ValueError(msg or '基盤地図情報の XML が読めません')
    # いちばん細かいタイルの間隔に合わせ、細かいタイルを先に置く（粗いタイルは欠測の所だけ埋める）
    res = lambda t: ((t['lat1'] - t['lat0']) / t['ny'], (t['lon1'] - t['lon0']) / t['nx'])
    dlat = min(res(t)[0] for t in T); dlon = min(res(t)[1] for t in T)
    T.sort(key=lambda t: res(t)[0])
    la0 = min(t['lat0'] for t in T); la1 = max(t['lat1'] for t in T); lo0 = min(t['lon0'] for t in T); lo1 = max(t['lon1'] for t in T)
    NX = int(round((lo1 - lo0) / dlon)); NY = int(round((la1 - la0) / dlat))
    Z = np.full((NY, NX), np.nan)
    for t in T:
        tl, tn = res(t); zz = t['z'].reshape(t['ny'], t['nx'])
        R0 = int(round((la1 - t['lat1']) / dlat)); C0 = int(round((t['lon0'] - lo0) / dlon))
        R1 = int(round((la1 - t['lat0']) / dlat)); C1 = int(round((t['lon1'] - lo0) / dlon))
        if abs(tl - dlat) < 1e-3 * dlat and abs(tn - dlon) < 1e-3 * dlon:
            m = np.isfinite(zz); sub = Z[R0:R0 + t['ny'], C0:C0 + t['nx']]; sub[m[:sub.shape[0], :sub.shape[1]]] = zz[:sub.shape[0], :sub.shape[1]][m[:sub.shape[0], :sub.shape[1]]]
        else:
            Rs = np.arange(max(0, R0), min(NY, R1)); Cs = np.arange(max(0, C0), min(NX, C1))
            rr = t['ny'] - 1 - np.floor((la1 - (Rs + 0.5) * dlat - t['lat0']) / tl).astype(int)
            cc = np.floor((lo0 + (Cs + 0.5) * dlon - t['lon0']) / tn).astype(int)
            ok_r = (rr >= 0) & (rr < t['ny']); ok_c = (cc >= 0) & (cc < t['nx'])
            Rs, rr = Rs[ok_r], rr[ok_r]; Cs, cc = Cs[ok_c], cc[ok_c]
            v = zz[np.ix_(rr, cc)]; sub = Z[np.ix_(Rs, Cs)]; put = ~np.isfinite(sub) & np.isfinite(v); sub[put] = v[put]; Z[np.ix_(Rs, Cs)] = sub
    latm = (la0 + la1) / 2 * math.pi / 180
    my = 111132.92 - 559.82 * math.cos(2 * latm) + 1.175 * math.cos(4 * latm)
    mx = 111412.84 * math.cos(latm) - 93.5 * math.cos(3 * latm)
    D = DEM(NX, NY, dlon * mx, dlat * my, lo0 * mx, la1 * my, Z, 'gsi')
    D.geo = dict(lat0=la0, lat1=la1, lon0=lo0, lon1=lo1)
    E = _square(D); E.skipped = bad
    return E


def _lzw(src):
    out = bytearray(); bitPos = 0; width = 9; prev = None
    dic = [bytes([i]) for i in range(256)] + [b'', b'']
    nbits = len(src) * 8
    while bitPos + width <= nbits:
        v = 0
        for b in range(width):
            byte = src[(bitPos + b) >> 3]; v = (v << 1) | ((byte >> (7 - ((bitPos + b) & 7))) & 1)
        bitPos += width
        if v == 257:
            break
        if v == 256:
            dic = [bytes([i]) for i in range(256)] + [b'', b'']; width = 9; prev = None
            continue
        if v < len(dic) and (v < 256 or v > 257):
            entry = dic[v]
        elif prev is not None:
            entry = prev + prev[:1]
        else:
            break
        out += entry
        if prev is not None:
            dic.append(prev + entry[:1])
        prev = entry
        if len(dic) + 1 >= (1 << width) and width < 12:
            width += 1
    return bytes(out)


def parse_tiff(buf):
    le = buf[:2] == b'II'; e = '<' if le else '>'
    if struct.unpack(e + 'H', buf[2:4])[0] != 42:
        raise ValueError('GeoTIFF ではありません（BigTIFF は未対応）')
    off = struct.unpack(e + 'I', buf[4:8])[0]; n = struct.unpack(e + 'H', buf[off:off + 2])[0]
    SZ = {1: ('B', 1), 2: ('c', 1), 3: ('H', 2), 4: ('I', 4), 5: ('II', 8), 6: ('b', 1), 7: ('B', 1), 8: ('h', 2), 9: ('i', 4), 11: ('f', 4), 12: ('d', 8), 16: ('Q', 8)}
    T = {}
    for i in range(n):
        p = off + 2 + i * 12
        tg, ty, cnt = struct.unpack(e + 'HHI', buf[p:p + 8])
        f, s = SZ.get(ty, ('B', 1)); size = s * cnt
        q = p + 8 if size <= 4 else struct.unpack(e + 'I', buf[p + 8:p + 12])[0]
        if ty == 2:
            T[tg] = buf[q:q + cnt].split(b'\0')[0].decode('latin1')
        elif ty == 5:
            v = struct.unpack(e + 'I' * (2 * cnt), buf[q:q + 8 * cnt]); T[tg] = [v[2 * k] / v[2 * k + 1] for k in range(cnt)]
        else:
            T[tg] = list(struct.unpack(e + f * cnt, buf[q:q + size]))
    W, H = T[256][0], T[257][0]; bps = T.get(258, [8])[0]; cmp = T.get(259, [1])[0]; spp = T.get(277, [1])[0]
    pred = T.get(317, [1])[0]; fmt = T.get(339, [1])[0]; bpp = bps // 8
    if spp != 1:
        raise ValueError('1 バンドの標高だけに対応しています')
    tiled = 322 in T; tw = T[322][0] if tiled else W; th = T[323][0] if tiled else T.get(278, [H])[0]
    offs = T[324] if tiled else T[273]; cnts = T[325] if tiled else T[279]
    nd = float(T[42113]) if 42113 in T else np.nan
    scale = T.get(33550, [1, 1]); tie = T.get(33922, [0, 0, 0, 0, 0, 0])
    dt = {3: {32: 'f4', 64: 'f8'}, 2: {8: 'i1', 16: 'i2', 32: 'i4'}, 1: {8: 'u1', 16: 'u2', 32: 'u4'}}[fmt][bps]
    dt = np.dtype(dt).newbyteorder(e)
    z = np.full((H, W), np.nan); nTx = -(-W // tw)
    for k, (o, c) in enumerate(zip(offs, cnts)):
        raw = bytes(buf[o:o + c])
        raw = raw if cmp == 1 else (_lzw(raw) if cmp == 5 else zlib.decompress(raw))
        rowsIn = th if tiled else min(th, H - k * th)
        b = bytearray(raw)
        if pred == 2:
            a = np.frombuffer(bytes(b), dtype=dt)[:rowsIn * tw].reshape(rowsIn, tw).copy()
            a = np.cumsum(a, axis=1, dtype=a.dtype)
            b = bytearray(a.astype(dt).tobytes())
        elif pred == 3:
            rowB = tw * bpp; out = bytearray(len(b))
            for r in range(rowsIn):
                o2 = r * rowB
                row = np.frombuffer(bytes(b[o2:o2 + rowB]), dtype=np.uint8).astype(np.int64)
                row = (np.cumsum(row) & 255).astype(np.uint8)
                planes = row.reshape(bpp, tw)   # MSB が先
                px = planes.T[:, ::-1] if le else planes.T
                out[o2:o2 + rowB] = px.astype(np.uint8).tobytes()
            b = out
        a = np.frombuffer(bytes(b), dtype=dt)[:rowsIn * tw].reshape(rowsIn, tw).astype(float)
        tx0 = (k % nTx) * tw if tiled else 0; ty0 = (k // nTx) * th if tiled else k * th
        h_ = min(rowsIn, H - ty0); w_ = min(tw, W - tx0)
        z[ty0:ty0 + h_, tx0:tx0 + w_] = a[:h_, :w_]
    bad = ~np.isfinite(z) | (z <= -9998) | (z >= 1e30)
    if np.isfinite(nd):
        bad |= z == nd
    z[bad] = np.nan
    csx = scale[0]; csy = scale[1] or scale[0]
    return _square(DEM(W, H, csx, csy, tie[3] - tie[0] * csx, tie[4] + tie[1] * csy, z, 'tif'))


def fill_holes(D):
    """欠測の穴埋め：DEM の縁につながらない欠測を、外側から 1 周ずつ、周り 8 セルの値の平均で埋める（dem.js の fillHoles と同じ）"""
    nx, ny = D.nx, D.ny
    z = np.array(D.z, dtype=float).reshape(ny, nx); nan = np.isnan(z)
    if not nan.any():
        return 0
    try:
        from scipy import ndimage
        lab, nl = ndimage.label(nan, structure=np.ones((3, 3)))
        edge = set(np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]]))) - {0}
        hole = nan & ~np.isin(lab, list(edge))
    except ImportError:   # scipy がないときは縁からたどる
        ext = np.zeros_like(nan); st = [(j, i) for j in range(ny) for i in (0, nx - 1) if nan[j, i]] + [(j, i) for i in range(nx) for j in (0, ny - 1) if nan[j, i]]
        for p in st:
            ext[p] = True
        while st:
            j, i = st.pop()
            for dj in (-1, 0, 1):
                for di in (-1, 0, 1):
                    jj, ii = j + dj, i + di
                    if 0 <= jj < ny and 0 <= ii < nx and nan[jj, ii] and not ext[jj, ii]:
                        ext[jj, ii] = True; st.append((jj, ii))
        hole = nan & ~ext
    n = 0; it = 0
    while hole.any() and it < 5000:
        it += 1
        f = np.isfinite(z); zp = np.pad(np.where(f, z, 0.0), 1); fp = np.pad(f, 1).astype(float)
        s = np.zeros_like(z); m = np.zeros_like(z)
        for dj in (-1, 0, 1):
            for di in (-1, 0, 1):
                if dj == 0 and di == 0:
                    continue
                s += zp[1 + dj:1 + dj + ny, 1 + di:1 + di + nx]; m += fp[1 + dj:1 + dj + ny, 1 + di:1 + di + nx]
        put = hole & (m > 0)
        if not put.any():
            break
        z[put] = s[put] / m[put]; hole &= ~put; n += int(put.sum())
    D.z = z.reshape(-1); D.filledHoles = n
    return n


def smooth(D, sigma):
    """平滑化：幅 w の箱形の平均を縦・横に 3 回（標準偏差 sigma [m] のガウス分布に近い）。欠測は除いて平均し、欠測のまま残す（dem.js の smooth と同じ）"""
    sc = sigma / D.cs
    if not (sc > 0.3):
        return D
    rr = max(1, int(math.floor((math.sqrt(4 * sc * sc + 1) - 1) / 2 + 0.5)))
    z = np.array(D.z, dtype=float).reshape(D.ny, D.nx)

    def box(a, axis):
        a = np.moveaxis(a, axis, -1); m = a.shape[-1]; f = np.isfinite(a)
        cs = np.concatenate([np.zeros(a.shape[:-1] + (1,)), np.cumsum(np.where(f, a, 0.0), axis=-1)], axis=-1)
        cn = np.concatenate([np.zeros(a.shape[:-1] + (1,)), np.cumsum(f, axis=-1)], axis=-1)
        i = np.arange(m); lo = np.maximum(0, i - rr); hi = np.minimum(m, i + rr + 1)
        n = cn[..., hi] - cn[..., lo]; out = np.where(n > 0, (cs[..., hi] - cs[..., lo]) / np.maximum(n, 1), a)
        out = np.where(f, out, np.nan)
        return np.moveaxis(out, -1, axis)
    for _ in range(3):
        z = box(z, 1); z = box(z, 0)
    import copy
    E = copy.copy(D); E.z = z.reshape(-1); E.smooth = sigma
    return E


def _expand(paths):
    """ZIP（中の ZIP も）を展開して (名前, バイト列) の並びにする"""
    import zipfile, io
    out = []

    def add(name, data, depth):
        if name.lower().endswith('.zip') and depth < 3:
            with zipfile.ZipFile(io.BytesIO(data)) as z:
                for n in z.namelist():
                    if not n.endswith('/'):
                        add(n, z.read(n), depth + 1)
        else:
            out.append((name.split('/')[-1], data))
    for p in paths:
        add(p, open(p, 'rb').read(), 0)
    return out


def load(paths):
    """DEM のファイル（1 つ、基盤地図情報の XML の並び、またはそれらの ZIP）を読む。
    XML を優先し、標高の XML でないものは除く"""
    if isinstance(paths, str):
        paths = [paths]
    fs = _expand(paths)
    xml = [(n, d) for n, d in fs if n.lower().endswith('.xml') and (b'tupleList' in d[:200000] or b'ElevationModel' in d[:4000] or 'dem' in n.lower())]
    if xml:
        return parse_gsi_xml([d.decode('utf-8', errors='replace') for n, d in xml])
    tif = [(n, d) for n, d in fs if n.lower().endswith(('.tif', '.tiff'))]
    if tif:
        return parse_tiff(tif[0][1])
    oth = [(n, d) for n, d in fs if n.lower().endswith(('.asc', '.txt', '.xyz', '.csv', '.dat'))]
    if not oth:
        raise ValueError('DEM のファイルが見つかりません：' + ', '.join(n for n, d in fs[:5]))
    t = oth[0][1].decode('utf-8', errors='replace')
    return parse_asc(t) if re.match(r'^\s*ncols', t, re.I) else parse_xyz(t)


def _square(D):
    if abs(D.csx - D.csy) < 1e-6 * D.csx:
        D.cs = D.csx
        return D
    cs = min(D.csx, D.csy); W = D.nx * D.csx; H = D.ny * D.csy
    nx = max(2, int(round(W / cs))); ny = max(2, int(round(H / cs)))
    Z = D.z.reshape(D.ny, D.nx)
    fy, fx = np.mgrid[0:ny, 0:nx]
    fx = np.clip((fx + 0.5) * cs / D.csx - 0.5, 0, D.nx - 1); fy = np.clip((fy + 0.5) * cs / D.csy - 0.5, 0, D.ny - 1)
    c0 = np.minimum(D.nx - 2, np.floor(fx).astype(int)); r0 = np.minimum(D.ny - 2, np.floor(fy).astype(int)); tx = fx - c0; ty = fy - r0
    v = (Z[r0, c0] * (1 - tx) + Z[r0, c0 + 1] * tx) * (1 - ty) + (Z[r0 + 1, c0] * (1 - tx) + Z[r0 + 1, c0 + 1] * tx) * ty
    near = Z[np.round(fy).astype(int), np.round(fx).astype(int)]
    v = np.where(np.isfinite(v), v, near)
    E = DEM(nx, ny, cs, cs, D.x0, D.yTop, v, D.src); E.geo = D.geo
    return E


# ---------------------------------------------------------------- 解析
def fill(D, eps=1e-5):
    nx, ny = D.nx, D.ny; N = nx * ny; z = D.z; zf = z.copy(); done = np.zeros(N, bool)
    H = []
    Zg = z.reshape(ny, nx); fin = np.isfinite(Zg)
    edge = np.zeros((ny, nx), bool); edge[0] = edge[-1] = True; edge[:, 0] = edge[:, -1] = True
    pad = np.pad(fin, 1, constant_values=False)
    for q in range(8):
        edge |= ~pad[1 + DR[q]:1 + DR[q] + ny, 1 + DC[q]:1 + DC[q] + nx]
    edge &= fin
    done[~fin.ravel()] = True
    for c in np.nonzero(edge.ravel())[0]:
        done[c] = True; H.append((z[c], int(c)))
    heapq.heapify(H)
    while H:
        v, c = heapq.heappop(H); cx = c % nx; cy = c // nx
        for q in range(8):
            xx = cx + DC[q]; yy = cy + DR[q]
            if xx < 0 or yy < 0 or xx >= nx or yy >= ny:
                continue
            o = yy * nx + xx
            if done[o]:
                continue
            done[o] = True
            if zf[o] <= zf[c]:
                zf[o] = zf[c] + eps * DL[q]
            heapq.heappush(H, (zf[o], o))
    return zf


def flow_dir(D, zf):
    nx, ny = D.nx, D.ny; Z = zf.reshape(ny, nx); best = np.full((ny, nx), -1, np.int8); bs = np.zeros((ny, nx))
    pad = np.pad(Z, 1, constant_values=np.nan)
    for q in range(8):
        o = pad[1 + DR[q]:1 + DR[q] + ny, 1 + DC[q]:1 + DC[q] + nx]
        s = (Z - o) / (DL[q] * D.cs)
        m = np.isfinite(s) & (s > bs)
        best[m] = q; bs[m] = s[m]
    best[~np.isfinite(Z)] = -2
    return best.ravel()


def down_of(D, d, c):
    q = d[c]
    if q < 0:
        return -1
    x = c % D.nx + DC[q]; y = c // D.nx + DR[q]
    return -1 if (x < 0 or y < 0 or x >= D.nx or y >= D.ny) else y * D.nx + x


class Analysis:
    pass


def analyze(D, thrArea=20000.0):
    nx, ny = D.nx, D.ny; N = nx * ny; a1 = D.cs * D.cs
    zf = fill(D); d = flow_dir(D, zf)
    fin = np.isfinite(zf)
    ordr = np.nonzero(fin)[0]; ordr = ordr[np.argsort(-zf[ordr], kind='stable')]
    dn = np.full(N, -1, np.int64)
    for c in ordr:
        dn[c] = down_of(D, d, c)
    acc = np.zeros(N); acc[ordr] = 1.0
    for c in ordr:
        if dn[c] >= 0:
            acc[dn[c]] += acc[c]
    thr = max(2, int(round(thrArea / a1)))
    stream = (acc >= thr) & (d != -2)
    hasUp = np.zeros(N, bool)
    for c in np.nonzero(stream)[0]:
        if dn[c] >= 0:
            hasUp[dn[c]] = True
    heads = [int(c) for c in np.nonzero(stream & ~hasUp)[0]]
    bid = np.full(N, -1, np.int64); basins = []
    for c in ordr[::-1]:
        if dn[c] < 0:
            bid[c] = len(basins); basins.append(dict(id=len(basins), outlet=int(c), n=0))
        else:
            bid[c] = bid[dn[c]]
    cnt = np.bincount(bid[bid >= 0], minlength=len(basins))
    for b in basins:
        b['n'] = int(cnt[b['id']]); b['area'] = b['n'] * a1; b['acc'] = acc[b['outlet']]
    A = Analysis()
    A.zf, A.dir, A.ord, A.acc, A.stream, A.heads, A.bid, A.basins, A.thr, A.dn = zf, d, ordr, acc, stream, heads, bid, basins, thr, dn
    return A


def upstream(D, A, c0):
    m = np.zeros(D.nx * D.ny, bool); m[c0] = True
    idx = int(np.nonzero(A.ord == c0)[0][0])
    for c in A.ord[idx - 1::-1] if idx > 0 else []:
        if A.dn[c] >= 0 and m[A.dn[c]]:
            m[c] = True
    return m


def local_slope(D, A, c, n=3):
    z0 = A.zf[c]; L = 0.0; cc = c
    for _ in range(n):
        dn = A.dn[cc]
        if dn < 0:
            break
        L += DL[A.dir[cc]] * D.cs; cc = dn
    return max(1e-3, (z0 - A.zf[cc]) / L) if L > 0 else 1e-2


def prepare_domain(D, A, mask, outlet, fill_dem=True):
    """流域を切り出し、出口の流向が右（+x）になるよう回して、出口から右端まで水路を掘る（dem.js の prepareDomain と同じ）"""
    nx, ny = D.nx, D.ny; N = nx * ny
    zz2 = (A.zf if fill_dem else D.z).copy(); st = A.stream.copy()
    for c in range(N):
        if not st[c] or not mask[c]:
            continue
        q0 = A.dir[c]
        if q0 < 0 or not (q0 % 2):
            continue
        dn = A.dn[c]
        if dn < 0:
            continue
        cx0 = c % nx; cy0 = c // nx; a1 = cy0 * nx + cx0 + DC[q0]; b1 = (cy0 + DR[q0]) * nx + cx0
        pick = a1 if zz2[a1] <= zz2[b1] else b1
        zm = 0.5 * (zz2[c] + zz2[dn])
        if zz2[pick] > zm:
            zz2[pick] = zm
        if mask[pick]:
            st[pick] = True
    ys, xs = np.divmod(np.nonzero(mask)[0], nx)
    x0 = max(0, xs.min() - 2); y0 = max(0, ys.min() - 2); x1 = min(nx - 1, xs.max() + 2); y1 = min(ny - 1, ys.max() + 2)
    W = x1 - x0 + 1; H = y1 - y0 + 1
    q = A.dir[outlet]
    if q >= 0 and down_of(D, A.dir, outlet) >= 0:
        ox, oy = DC[q], DR[q]
    else:
        cx = outlet % nx; cy = outlet // nx; dl = [nx - 1 - cx, cy, cx, ny - 1 - cy]; km = int(np.argmin(dl))
        ox = [1, 0, -1, 0][km]; oy = [0, -1, 0, 1][km]
    rot = (0 if ox > 0 else 2) if abs(ox) >= abs(oy) else (1 if oy < 0 else 3)
    RW, RH = (H, W) if rot % 2 else (W, H)

    def fwd(col, row):
        a, b = col - x0, row - y0
        return [(a, b), (H - 1 - b, a), (W - 1 - a, H - 1 - b), (b, W - 1 - a)][rot]

    def inv(i, j):
        return [(x0 + i, y0 + j), (x0 + j, y0 + H - 1 - i), (x0 + W - 1 - i, y0 + H - 1 - j), (x0 + W - 1 - j, y0 + i)][rot]
    po = fwd(outlet % nx, outlet // nx); NW = RW + 3
    zmaxB = float(np.nanmax(np.where(mask, D.z, np.nan)))
    rows = np.zeros((RH, NW)); bas = np.zeros((RH, NW), np.uint8); chn = np.zeros((RH, NW), np.uint8); cellOf = np.full((RH, NW), -1, np.int64)
    for j in range(RH):
        for i in range(RW):
            o = inv(i, j); oc = o[1] * nx + o[0]
            zz = zz2[oc]
            rows[j, i] = zz if np.isfinite(zz) else zmaxB
            bas[j, i] = 1 if mask[oc] else 0; chn[j, i] = 1 if (mask[oc] and st[oc]) else 0; cellOf[j, i] = oc
        rows[j, RW:] = rows[j, RW - 1]
    Iout = max(0.01, local_slope(D, A, outlet, 4)); zo = zz2[outlet]
    cor = np.zeros_like(bas, dtype=np.uint8)
    # 幅は 5 セル（中央が谷筋）。1 セル幅だと、格子を粗くしたとき両側の山と混ざって通路がふさがるため
    for i in range(po[0] + 1, NW):
        zc = zo - Iout * (i - po[0]) * D.cs
        for dj in range(-2, 3):
            jj = po[1] + dj
            if jj < 0 or jj >= RH:
                continue
            if dj and bas[jj, i] and rows[jj, i] < zc:
                continue
            if not bas[jj, i]:
                cor[jj, i] = 1   # 出口の先の水路（流域の外だったセル）
            rows[jj, i] = zc; bas[jj, i] = 1
            if not dj:
                chn[jj, i] = 1

    class Dom:
        pass
    dom = Dom()
    dom.rows, dom.cs, dom.basin, dom.chan, dom.cellOf, dom.rot, dom.nx, dom.ny, dom.outlet, dom.Iout = rows, D.cs, bas, chn, cellOf, rot, NW, RH, po, Iout
    dom.corr = cor

    def toXY(cell):
        p = fwd(cell % nx, cell // nx)
        return [(p[0] + 0.5) * D.cs, (RH - p[1] - 0.5) * D.cs]

    def dirXY(cell):
        qq = A.dir[cell]
        if qq < 0:
            return [1.0, 0.0]
        vx, vy = DC[qq], -DR[qq]; L = math.hypot(vx, vy); vx /= L; vy /= L
        for _ in range(rot):
            vx, vy = vy, -vx
        return [vx, vy]
    dom.toXY, dom.dirXY = toXY, dirXY
    return dom


def terrain_params(dom):
    return {'rows': dom.rows.tolist(), 'cs': dom.cs, 'basin': dom.basin.tolist(), 'chan': dom.chan.tolist(), 'corr': dom.corr.tolist(), 'Iout': dom.Iout}


def synth_dem(nx=160, ny=120, cs=10.0, seed=3, iters=40):
    """例の DEM（架空の山地流域）。乱数の起伏に、stream power の侵食で谷を刻む（dem.js の synthDEM と同じ考え方。値は同じではない）"""
    rng = np.random.default_rng(seed)
    Y, X = np.mgrid[0:ny, 0:nx]; X = X / nx; Y = Y / ny
    z = 420 * (1 - X) + 90 * np.abs(Y - 0.5) ** 1.5 * 2 ** 1.5 + 40
    for amp, f in ((34, 4), (16, 8), (7, 16)):
        g = rng.standard_normal((f + 2, f + 2))
        xi = X * f; yi = Y * f; i0 = xi.astype(int); j0 = yi.astype(int); tx = xi - i0; ty = yi - j0
        sx = tx * tx * (3 - 2 * tx); sy = ty * ty * (3 - 2 * ty)
        z = z + amp * ((g[j0, i0] * (1 - sx) + g[j0, i0 + 1] * sx) * (1 - sy) + (g[j0 + 1, i0] * (1 - sx) + g[j0 + 1, i0 + 1] * sx) * sy)
    D = DEM(nx, ny, cs, cs, 0.0, ny * cs, z, 'synth')
    for _ in range(iters):
        A = analyze(D, 1e12)
        S = np.zeros(nx * ny)
        for c in range(nx * ny):
            dn = A.dn[c]
            if dn >= 0:
                S[c] = max(0.0, (A.zf[c] - A.zf[dn]) / (DL[A.dir[c]] * cs))
        Z = D.z.reshape(ny, nx); P_ = np.pad(Z, 1, mode='edge')
        lap = P_[1:-1, 2:] + P_[1:-1, :-2] + P_[2:, 1:-1] + P_[:-2, 1:-1] - 4 * Z
        zn = D.z - 0.9 * np.sqrt(A.acc * cs * cs) * S * 0.06 + 0.08 * lap.ravel()
        zn = zn.reshape(ny, nx); zn[:, -1] = Z[:, -1]
        D.z = zn.ravel()
    D.z = fill(D, 0.01 * cs)   # 窪地は 1 % の勾配の斜面で埋める
    return D


# ---------------------------------------------------------------- 雨の流出（hydro.js と同じ）
def rain_uniform(tab, t):
    if not tab:
        return 0.0
    tm = t / 60; v = 0.0
    for r in tab:
        if tm >= r[0] - 1e-9:
            v = r[1]
        else:
            break
    if len(tab) > 1 and tm >= tab[-1][0] + (tab[-1][0] - tab[-2][0]):
        v = 0.0
    return max(0.0, v)


def loss_step(L, st, r, dt):
    dth = dt / 3600
    if L.get('type') == 'f1rsa':
        st['R'] = st.get('R', 0.0) + r * dth
        R0 = L.get('R0') or 0.0
        f = 0.0 if st['R'] < R0 else ((L.get('f1') if L.get('f1') is not None else 0.5) if st['R'] < R0 + (L.get('Rsa') or 0.0) else 1.0)
        return f * r
    if L.get('type') == 'ks':
        Smax = max(0.0, (L.get('D') or 1) * ((L.get('n') or 0.45) - (L.get('th0') or 0.25)) * 1000); ks = (L.get('ks') or 1e-5) * 3.6e6
        S = st.get('S', 0.0)
        if S >= Smax:
            return r
        inf = min(r, ks); add = inf * dth
        if S + add > Smax:
            over = (S + add - Smax) / dth; st['S'] = Smax
            return r - inf + over
        st['S'] = S + add
        return r - inf
    return (L.get('f') if L.get('f') is not None else 0.8) * r


def tc_pwri(L, S):
    """土研式（自然流域）の洪水到達時間 [s]：T = 1.67×10⁻³ (L/√S)^0.7 [h]"""
    return 3600 * 1.67e-3 * (max(1.0, L) / math.sqrt(max(1e-3, S))) ** 0.7


def lateral_inflow(D, A, mask, table, loss, route, dt=60.0, tEnd=21600.0):
    """一様な雨（表 [[分, mm/h]]）で、谷筋のセルごとの横流入の時系列を求める。返り値：cells、Q（[セル][時刻] m3/s）、dt"""
    N = D.nx * D.ny; a1 = D.cs * D.cs; nT = max(2, int(math.ceil(tEnd / dt)) + 1)
    dist = np.zeros(N); recv = np.full(N, -1, np.int64)
    for c in A.ord[::-1]:
        if not mask[c]:
            continue
        if A.stream[c]:
            recv[c] = c; continue
        dn = A.dn[c]
        if dn < 0 or not mask[dn]:
            continue
        dist[c] = dist[dn] + DL[A.dir[c]] * D.cs; recv[c] = recv[dn]
    cells = [int(c) for c in range(N) if mask[c] and A.stream[c]]
    idx = {c: i for i, c in enumerate(cells)}; n = np.zeros(len(cells)); Lm = np.zeros(len(cells)); Ss = np.zeros(len(cells))
    for c in range(N):
        if not mask[c] or recv[c] < 0 or recv[c] not in idx:
            continue
        i = idx[recv[c]]; n[i] += 1; Lm[i] = max(Lm[i], dist[c])
        dn = A.dn[c]; Ss[i] += max(0.0, (A.zf[c] - A.zf[dn]) / (DL[A.dir[c]] * D.cs)) if dn >= 0 else 0.01
    Tc = np.array([route['tc'] * 60 if (route.get('tc') or 0) > 0 else max(300, min(7200, tc_pwri(Lm[i] + D.cs / 2, Ss[i] / n[i] if n[i] else 0.1))) for i in range(len(cells))])
    st = {}; re_ = np.zeros(nT)
    for k in range(nT):
        re_[k] = loss_step(loss, st, rain_uniform(table, k * dt), dt)
    cum = np.cumsum(re_)
    Q = np.zeros((len(cells), nT)); K = route.get('K') or 10.0; P = route.get('P') or 0.6; lag = int(round((route.get('Tl') or 0) * 60 / dt))
    for i in range(len(cells)):
        if route.get('type') == 'storage':
            s = 0.0; q = np.zeros(nT)
            for k in range(nT):
                dth = dt / 3600; q0 = (max(0.0, s) / K) ** (1 / P); s1 = s + dth * (re_[k] - q0); q1 = (max(0.0, s1) / K) ** (1 / P)
                s = max(0.0, s + 0.5 * dth * (2 * re_[k] - q0 - q1)); q[k] = (s / K) ** (1 / P)
        else:
            nb = max(1, int(round(Tc[i] / dt))); q = (cum - np.concatenate([np.zeros(nb), cum[:-nb]])[:nT]) / nb
        Qi = q * n[i] * a1 / 3.6e6
        Q[i, lag:] = Qi[:nT - lag] if lag else Qi
    return cells, Q, dt
