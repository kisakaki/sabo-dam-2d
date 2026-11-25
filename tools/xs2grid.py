#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
xs2grid.py — 横断測量・縦断測量の値から、河床変動シミュレータに読み込める
格子地形 CSV（行列形式）を作る。標準ライブラリだけで動く。

使い方
  テスト形状を書き出す
    python3 xs2grid.py --make-test
      → test_cross_sections.csv, test_longitudinal.csv を書き出し、続けて test_grid.csv を作る

  測量値から格子地形を作る
    python3 xs2grid.py --xs 横断.csv [--long 縦断.csv] [--dx 2] [--align thalweg] [--wall 1.0] --out grid.csv

入力の形式（1行目などの見出し行と、# で始まる行は読み飛ばす。区切りはカンマ・タブ・空白）
  文字コードは UTF-8（BOM 付きも可）か Shift_JIS。出力する CSV は BOM 付き UTF-8 で、Excel で開いても文字化けしない。
  横断 CSV：x, s, z
    x：断面の位置（流下方向の距離）[m]
    s：左岸からの距離 [m]（右岸へ向かって増える）
    z：標高 [m]
    同じ x の行が1つの断面になる。断面は2つ以上必要。
  縦断 CSV：x, z（省略可）
    最深河床高の縦断形。与えると各断面の形はそのままに、最深部の高さをこの縦断形に合わせる。
    省略すると、最深部の高さは断面の間を直線でつなぐ。

断面の間の補間（比例配分）
  1. 各断面を「最深部より左岸側」と「右岸側」の2つに分け、それぞれの幅で 0〜1 に正規化する。
  2. 隣り合う2断面の同じ正規化位置どうしを、x 方向に直線で補間する。
     左岸端・最深部・右岸端の位置も x 方向に直線で補間する。
     これで最深部（澪筋）は断面の間でもつながる。
  3. 形は最深部からの比高で補間し、最深部の高さ（縦断形）に足す。
  断面の端より外側は、端の標高から勾配 --wall で立ち上がる斜面とする（谷壁）。

断面の横方向のそろえ方（--align）
  left    ：左岸端をそろえる
  thalweg ：最深部をそろえる（既定）
  center  ：断面の中央をそろえる

出力（シミュレータの「CSV を読み込む」の行列形式）
  1行目が左岸側、左の列が上流。格子間隔は --dx。セルの中心の標高を書く。
  x は最上流の断面（または縦断の始点）を 0 とする。
"""
import argparse
import math
import sys


def read_lines(path):
    """UTF-8（BOM 付きも可）で読み、だめなら Shift_JIS（cp932）で読む。Excel で保存した CSV にも対応する。"""
    for enc in ('utf-8-sig', 'cp932'):
        try:
            with open(path, encoding=enc) as f:
                return f.read().splitlines()
        except UnicodeDecodeError:
            continue
    with open(path, encoding='utf-8', errors='replace') as f:
        return f.read().splitlines()


def read_table(path, ncol):
    rows = []
    if True:
        for ln in read_lines(path):
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            tok = [t for t in ln.replace('\t', ',').replace(';', ',').replace(' ', ',').split(',') if t != '']
            try:
                v = [float(t) for t in tok]
            except ValueError:
                continue          # 見出し行
            if len(v) < ncol:
                continue
            rows.append(v[:ncol])
    return rows


def interp1(xs, ys, x):
    """区分線形補間（範囲外は端の値）"""
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        m = (lo + hi) // 2
        if xs[m] <= x:
            lo = m
        else:
            hi = m
    t = (x - xs[lo]) / (xs[hi] - xs[lo])
    return ys[lo] * (1 - t) + ys[hi] * t


def build_sections(xs_rows, align):
    by = {}
    for x, s, z in xs_rows:
        by.setdefault(round(x, 6), []).append((s, z))
    secs = []
    for x in sorted(by):
        pts = sorted(by[x])
        if len(pts) < 3:
            raise SystemExit('x = %g の断面の点が 3 点未満です。' % x)
        s = [p[0] for p in pts]
        z = [p[1] for p in pts]
        zmin = min(z)
        idx = [i for i, v in enumerate(z) if v <= zmin + 1e-9]
        st = sum(s[i] for i in idx) / len(idx)          # 最深部の位置（同じ高さが続けば中央）
        secs.append({'x': x, 's': s, 'zr': [v - zmin for v in z], 'zmin': zmin, 'st': st, 'sl': s[0], 'sr': s[-1]})
    if len(secs) < 2:
        raise SystemExit('断面が 2 つ以上必要です。')
    # そろえ方：u は左岸側の基準線からの距離
    if align == 'left':
        for c in secs:
            c['off'] = -c['sl']
    elif align == 'center':
        U = max((c['sr'] - c['sl']) / 2 for c in secs)
        for c in secs:
            c['off'] = U - (c['sl'] + c['sr']) / 2
    else:
        U = max(c['st'] - c['sl'] for c in secs)
        for c in secs:
            c['off'] = U - c['st']
    for c in secs:
        c['uL'] = c['sl'] + c['off']
        c['uT'] = c['st'] + c['off']
        c['uR'] = c['sr'] + c['off']
    return secs


def section_z(c, eta_side, eta):
    """断面 c の、最深部の左（eta_side=-1）または右（+1）側の正規化位置 eta（0〜1）での比高"""
    if eta_side < 0:
        s = c['sl'] + (c['st'] - c['sl']) * eta
    else:
        s = c['st'] + (c['sr'] - c['st']) * eta
    return interp1(c['s'], c['zr'], s)


def make_grid(secs, longp, dx, wall):
    xs_ = [c['x'] for c in secs]
    xa = min(xs_ + ([p[0] for p in longp] if longp else []))
    xb = max(xs_ + ([p[0] for p in longp] if longp else []))
    L = xb - xa
    W = max(c['uR'] for c in secs)
    nx = max(10, int(round(L / dx)))
    ny = max(3, int(round(W / dx)) + 1)
    lx = [p[0] for p in longp] if longp else None
    lz = [p[1] for p in longp] if longp else None
    grid = [[0.0] * nx for _ in range(ny)]
    for i in range(nx):
        x = xa + (i + 0.5) * dx
        # 両側の断面と重み
        if x <= xs_[0]:
            a = b = secs[0]; t = 0.0
        elif x >= xs_[-1]:
            a = b = secs[-1]; t = 0.0
        else:
            k = 1
            while xs_[k] < x:
                k += 1
            a, b = secs[k - 1], secs[k]
            t = (x - a['x']) / (b['x'] - a['x'])
        uL = a['uL'] * (1 - t) + b['uL'] * t
        uT = a['uT'] * (1 - t) + b['uT'] * t
        uR = a['uR'] * (1 - t) + b['uR'] * t
        zt = interp1(lx, lz, x) if longp else a['zmin'] * (1 - t) + b['zmin'] * t
        zL = section_z(a, -1, 0.0) * (1 - t) + section_z(b, -1, 0.0) * t
        zR = section_z(a, 1, 1.0) * (1 - t) + section_z(b, 1, 1.0) * t
        for j in range(ny):
            u = (j + 0.5) * dx                     # j = 0 が左岸側（1行目）
            if u < uL:
                zr = zL + wall * (uL - u)
            elif u > uR:
                zr = zR + wall * (u - uR)
            elif u <= uT:
                e = (u - uL) / (uT - uL) if uT > uL else 1.0
                zr = section_z(a, -1, e) * (1 - t) + section_z(b, -1, e) * t
            else:
                e = (u - uT) / (uR - uT) if uR > uT else 0.0
                zr = section_z(a, 1, e) * (1 - t) + section_z(b, 1, e) * t
            grid[j][i] = zt + zr
    return grid, L, W


def write_grid(path, grid):
    with open(path, 'w', encoding='utf-8-sig', newline='\r\n') as f:   # BOM 付き UTF-8：Excel でも文字化けしない
        for row in grid:
            f.write(','.join('%.3f' % v for v in row) + '\n')


def make_test(prefix='test_'):
    """テスト形状：区間長 400 m の山地河川。縦断勾配 1/15 → 1/25（堆積区間）→ 1/15。
    断面は 50 m ごと。低水路（台形）＋高水敷（段丘）＋谷壁の複断面で、低水路の位置と幅が変わる。"""
    def thalweg(x):
        if x <= 150:
            z = 25.0 - x / 15.0
        elif x <= 300:
            z = 15.0 - (x - 150) / 25.0
        else:
            z = 9.0 - (x - 300) / 15.0
        return z - 0.3 * math.sin(2 * math.pi * x / 60.0)    # 瀬と淵
    # 断面ごとの形：(x, 谷底幅, 低水路幅, 低水路の中心（左岸から）, 低水路深, 高水敷の横断勾配, 谷壁高)
    spec = [(0, 50, 22, 25, 1.5, 0.02, 7), (50, 52, 24, 20, 1.6, 0.02, 7), (100, 48, 20, 30, 1.8, 0.02, 7),
            (150, 56, 26, 24, 1.4, 0.02, 7), (200, 66, 32, 30, 1.2, 0.015, 7), (250, 72, 36, 40, 1.0, 0.015, 7),
            (300, 64, 30, 34, 1.2, 0.015, 7), (350, 52, 22, 22, 1.6, 0.02, 7), (400, 50, 22, 26, 1.6, 0.02, 7)]
    xs_lines = ['# テスト横断形状　x：断面の位置 [m]（上流端から）　s：左岸からの距離 [m]　z：標高 [m]', 'x,s,z']
    for (x, Bv, Bc, cc, hc, cs, Hw) in spec:
        zt = thalweg(x)
        wl = 14.0                                  # 谷壁の水平幅（片側）
        pts = []
        # 左岸の谷壁 → 左の高水敷 → 低水路 → 右の高水敷 → 右岸の谷壁
        left_floor, right_floor = wl, wl + Bv
        c0, c1 = wl + cc - Bc / 2, wl + cc + Bc / 2
        for s in [0, 4, 8, 12]:
            pts.append((s, zt + hc + cs * (c0 - left_floor) + Hw * (1 - s / wl)))
        s = left_floor
        while s < c0 - 1e-9:
            pts.append((s, zt + hc + cs * (c0 - s)))
            s += 4
        pts.append((c0, zt + hc))
        bank = min(3.0, Bc / 6)
        pts.append((c0 + bank, zt + 0.15))
        s = c0 + bank + 3
        while s < c1 - bank - 1e-9:
            pts.append((s, zt + 0.15 * (1 - math.sin(math.pi * (s - c0) / Bc) ** 2)))
            s += 3
        pts.append((cc + wl, zt))                    # 最深部
        pts.append((c1 - bank, zt + 0.15))
        pts.append((c1, zt + hc))
        s = c1 + 4
        while s < right_floor - 1e-9:
            pts.append((s, zt + hc + cs * (s - c1)))
            s += 4
        zr = zt + hc + cs * (right_floor - c1)
        for k, s in enumerate([right_floor, right_floor + 4, right_floor + 8, right_floor + 12, right_floor + wl]):
            pts.append((s, zr + Hw * (s - right_floor) / wl))
        pts = sorted(set((round(a, 2), round(b, 3)) for a, b in pts))
        for s_, z_ in pts:
            xs_lines.append('%g,%g,%.3f' % (x, s_, z_))
    lo_lines = ['# テスト縦断形状　x：上流端からの距離 [m]　z：最深河床高 [m]', 'x,z']
    for k in range(0, 41):
        x = k * 10
        lo_lines.append('%g,%.3f' % (x, thalweg(x)))
    with open(prefix + 'cross_sections.csv', 'w', encoding='utf-8-sig', newline='\r\n') as f:
        f.write('\n'.join(xs_lines) + '\n')
    with open(prefix + 'longitudinal.csv', 'w', encoding='utf-8-sig', newline='\r\n') as f:
        f.write('\n'.join(lo_lines) + '\n')
    return prefix + 'cross_sections.csv', prefix + 'longitudinal.csv'


def main():
    ap = argparse.ArgumentParser(description='横断・縦断の測量値から格子地形 CSV を作る')
    ap.add_argument('--xs', help='横断 CSV（x, s, z）')
    ap.add_argument('--long', help='縦断 CSV（x, z）。省略可')
    ap.add_argument('--dx', type=float, default=2.0, help='格子間隔 [m]（既定 2）')
    ap.add_argument('--align', choices=['thalweg', 'left', 'center'], default='thalweg', help='断面の横方向のそろえ方')
    ap.add_argument('--wall', type=float, default=1.0, help='断面の端より外側の谷壁の勾配（既定 1.0）')
    ap.add_argument('--out', default='grid.csv', help='出力ファイル')
    ap.add_argument('--make-test', action='store_true', help='テスト形状の CSV を書き出し、格子地形も作る')
    a = ap.parse_args()
    if a.make_test:
        fx, fl = make_test()
        a.xs, a.long = fx, fl
        if a.out == 'grid.csv':
            a.out = 'test_grid.csv'
        print('書き出し：', fx, fl)
    if not a.xs:
        ap.error('--xs を指定してください（または --make-test）')
    secs = build_sections(read_table(a.xs, 3), a.align)
    longp = sorted(read_table(a.long, 2)) if a.long else None
    grid, L, W = make_grid(secs, longp, a.dx, a.wall)
    write_grid(a.out, grid)
    print('書き出し：%s（%d 行 × %d 列、格子間隔 %g m、区間長 %.0f m、幅 %.0f m、断面 %d）'
          % (a.out, len(grid), len(grid[0]), a.dx, L, W, len(secs)))


if __name__ == '__main__':
    main()
