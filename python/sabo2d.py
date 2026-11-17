#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sabo2d.py — 連続砂防堰堤 平面2次元河床変動計算（Python 版）

ブラウザ版（連続堰堤 河床変動シミュレータ）と同じ計算を、NumPy で書き直したもの。
Jupyter Notebook からも、コマンドラインからも使える。

支配方程式
  流れ  平面2次元浅水流方程式（保存形）
          ∂h/∂t + ∂(hu)/∂x + ∂(hv)/∂y = 0
          ∂(hu)/∂t + ∂(hu²+gh²/2)/∂x + ∂(huv)/∂y = -gh ∂z/∂x - g n² u√(u²+v²)/h^(1/3)
          ∂(hv)/∂t + ∂(huv)/∂x + ∂(hv²+gh²/2)/∂y = -gh ∂z/∂y - g n² v√(u²+v²)/h^(1/3)
  河床  粒径階 k ごとの流砂の連続式（Exner）  (1-λ) ∂z_k/∂t + ∂q_bx,k/∂x + ∂q_by,k/∂y = 0
離散化
  流れ  セル中心の有限体積法。HLL 近似リーマン解法＋静水圧再構築（Audusse ら 2004、Chen・Noelle 2017）。
        h と水位 η を minmod 制限付きで区分線形に再構築。時間積分は陽的 Euler、摩擦は半陰的。
  流砂  芦田・道上式（1972）。限界掃流力は芦田・高橋・水山（1978）式(24) τ*c = 0.04×10^(1.72 i)、
        岩垣式（1956）、一定値から選ぶ。粒径別は芦田・高橋・水山の結果（u*ci = u*cm）か修正 Egiazaroff 式。
        流砂の向きは流速方向＋河床勾配による重力効果（長谷川 1981 型の係数）。
  粒度  交換層（平野 1971）＋多層の堆積層。
  側岸  安息角を超えた斜面の崩落。
  堰堤  1 セル幅の固定床。天端・水通し・袖部を標高で与え、河道の中心線（両岸）に直角に置く。

ブラウザ版との違い
  ・安息角を超えた斜面の崩落は、セルを市松に分けて同時に処理する（ブラウザ版は 1 面ずつ順に処理）。
    そのため結果は完全には一致しないが、体積は同じく保存する。
  ・それ以外の式・係数・境界条件はブラウザ版と同じ。地形と初期の乱数も同じ値になる。

使い方（コマンドライン）
  python sabo2d.py --preset A --out result            # モデル渓流 A（床固工群）を計算して図を保存
  python sabo2d.py --preset B --minutes 60            # 出水の継続時間を 60 分に
  python sabo2d.py --grid test_grid.csv --cs 2 --dams "100:3:0,200:3:0"
  python sabo2d.py --xs test_cross_sections.csv --long test_longitudinal.csv --dams "150:4:12"
  python sabo2d.py --help

使い方（Python / Jupyter）
  import sabo2d as sb
  P = sb.preset('A')                 # 既定値＋モデル渓流 A
  P['hydro']['Qpeak'] = 80
  sim = sb.Sim(P)
  sb.run(sim)                        # 出水が終わるまで計算
  sb.plot_all(sim)

  # 測量値（横断 x,s,z と縦断 x,z）や格子 CSV を読むとき
  T = sb.read_survey('test_cross_sections.csv', 'test_longitudinal.csv', dx=2.0)   # または sb.read_grid('test_grid.csv', cs=2.0)
  P = sb.use_terrain(sb.default_params(), T)
  P['dams'] = [{'x': 150, 'H': 4.0, 'b': 12.0}]

  Jupyter Notebook 版（sabo2d.ipynb）は、このファイルと同じコードをセルに分けたもの。

必要なもの：Python 3.8 以上、NumPy、Matplotlib（図を描くとき）。
"""
import copy
import math
import sys
import time

import numpy as np

G = 9.81
LN19 = math.log(19.0)
M32 = 0xFFFFFFFF


# =============================================================================
# 条件（既定値）
# =============================================================================
def default_params():
    """計算条件の既定値。ブラウザ版の defaultParams と同じ。"""
    return {
        # ---- 河道地形
        'L': 400.0,            # 計算区間長 [m]
        'B': 50.0,             # 谷底（河床）幅 [m]
        'widthTable': [],      # 川幅の変化点 [[x, 幅], ...]（線形補間。空なら B で一様）
        'dx': 2.0,             # 格子幅 [m]
        'slope': 0.05,         # 縦断勾配 I
        'crossSlope': 0.02,    # 横断勾配（谷の中央・澪筋へ向かう）
        'perturb': 0.05,       # 初期河床の擾乱振幅 [m]
        'seed': 1,             # 乱数の種
        'erodDepth': 3.0,      # 侵食可能深 [m]
        'sideCells': 5,        # 谷壁のセル数（片側）
        'sideGrad': 0.8,       # 谷壁の勾配
        # 平面形：{'mode': 'none'|'sine'|'table', 'amp', 'wl', 'x0', 'x1', 'table': [[x, 横ずれ, 川幅(省略可)], ...]}
        'bend': {'mode': 'none'},
        # 渓流の作り込み（数値で指定する地形のみ。0 のとき使わない）
        'valley': {'chD': 0.0, 'chW': 0.0, 'chWl': 0.0, 'relief': 0.0, 'wallH': 0.0, 'wallRough': 0.0,
                   'fans': []},   # fans: [[x, 岸(+1 左岸 / -1 右岸), 半径, 高さ], ...]
        # 読み込む地形（None のとき数値で指定）：{'rows': 2次元リスト（1行目が左岸側、左の列が上流）, 'cs': 格子間隔}
        'terrain': None,
        'floorH': 3.0,         # 読み込んだ地形で、各横断の最低点から谷底とみなす高さ [m]
        'csvSlope': 0.0,       # >0 のとき読み込んだ地形の縦断勾配をこの値に傾け直す
        # ---- 流れ
        'manning': 0.050,
        'cfl': 0.4,
        'order': 2,            # 1: 区分一定, 2: 区分線形
        # ---- 土砂
        'grains': [            # d: 粒径 [m], bed: 河床材料の割合, feed: 供給土砂の割合
            {'d': 0.005, 'bed': 0.20, 'feed': 0.30},
            {'d': 0.020, 'bed': 0.30, 'feed': 0.35},
            {'d': 0.060, 'bed': 0.30, 'feed': 0.25},
            {'d': 0.150, 'bed': 0.20, 'feed': 0.10},
        ],
        'porosity': 0.4,
        'sg': 2.65,
        'tcMethod': 'atm',     # 'atm' 芦田・高橋・水山 | 'iwagaki' 岩垣 | 'const' 一定値 tauCm
        'hiding': 'atm',       # 'atm' | 'egiazaroff'
        'tauCm': 0.05,
        'reposeDeg': 35.0,     # 流砂の勾配項の上限角 [deg]
        # ---- 側岸侵食
        'bank': {'erodH': 0.0,       # 谷壁のうち谷底から erodH [m] までを侵食できる側岸にする（0 で固定床）
                 'reposeWet': 30.0,  # 水中の安息角 [deg]（高橋・里深・奥村 1996 は 30°）
                 'reposeDry': 50.0,  # 陸上の安息角 [deg]
                 'dryEro': 0.0,      # 水際の濡れたセルの侵食のうち、隣の乾いた側岸から削る割合（0〜1）
                 'collapse': 'repose'},  # 'repose' 安息角を超えた分を少しずつ崩す | 'hasegawa' 長谷川の崩落モデル（長田ら 1996 式17）
        # ---- 2次流（湾曲部のらせん流）。長田・細田・村本・Rahman（1996）
        'sec': {'bed': False,    # 河床付近の流れを曲率の内側へ向け、掃流砂の向きを変える（式11・14：u_nb/u_sb = N* h/r）
                'mom': False,    # 2次流による運動量の分散項（式19〜23、Kalkwijk・de Vriend）
                'Nstar': 11.5,   # 式14の係数 N*
                'hrMax': 0.05},  # h/r の上限（弱い曲率を前提とする式のため）
        # ---- 流砂の計算法：'eq' 平衡流砂量（芦田・道上式）| 'noneq' 侵食・堆積速度式（高橋・里深 1999）
        'transport': 'eq',
        'noneq': {'de': 0.6, 'dd': 0.6, 'hEdge': 0.0, 'fe': 2.0, 'fd': 10.0},
        # ---- 巨礫（動かない礫）を水みちの両岸に置く（高橋・里深・奥村 1996）
        'boulders': {'L': 0.0, 'mode': 'parallel', 'size': 1},
        'activeLayer': 0.0,    # 交換層厚（0 のとき最大粒径）
        'layerTh': 0.25,
        'morfac': 4.0,         # 河床変動の加速係数
        'gridDamp': 0.3,
        'gravK': 1.0,
        # ---- ハイドログラフ
        'hydro': {'mode': 'tri', 'Qbase': 10.0, 'Qpeak': 60.0, 'tPeak': 3600.0, 'duration': 10800.0,
                  'Cs': 0.005, 'table': []},   # table: [[t(s), Q, Qs(省略可)], ...]
        # ---- 堰堤
        'dams': [],            # [{'x': 位置, 'H': 堤高, 'b': 水通し幅（0 で全幅）}]
        'wingHeight': 1.5,
        'prefill': 0.0,        # >0 のとき堰堤の上流をはじめから満砂にする（堆砂勾配 = 元の勾配 × prefill）
    }


# モデル渓流（写真の河道を参考にした架空の地形）。ブラウザ版のボタンと同じ値
PRESETS = {
    # 床固工群が続く、谷底の広い渓流（早月川の床固工群を参考）
    'A': dict(L=600, B=60, I=18, S=3, grad=1.0, wing=1.5, dx=2.5, fill=0.5,
              bend=[[0, 0, 52], [70, 6, 58], [150, -6, 66], [230, -14, 70], [310, -4, 62], [380, 10, 74],
                    [460, 16, 68], [530, 6, 60], [600, 0, 56]],
              V=dict(chD=0.7, chW=14, chWl=0, relief=0.4, wallH=16, wallRough=3,
                     fans=[[120, -1, 34, 3.2], [340, 1, 40, 3.5], [505, -1, 28, 2.5]]),
              dams=[{'x': x, 'H': 3.0, 'b': 36.0} for x in (75, 135, 195, 255, 315, 375, 435, 495, 555)]),
    # 狭窄部に堰堤を置き、その上流に広い堆砂域がある渓流（南股入の堰堤群を参考）
    'B': dict(L=660, B=44, I=28, S=2, grad=1.0, wing=1.5, dx=2.5, fill=0.5,
              bend=[[0, 0, 40], [70, 10, 38], [140, -6, 42], [200, -2, 30], [235, 0, 24], [270, 4, 36],
                    [340, 18, 84], [410, 24, 110], [480, 18, 112], [540, 8, 70], [575, 2, 34], [610, -2, 40],
                    [660, -6, 42]],
              V=dict(chD=0.8, chW=12, chWl=0, relief=0.35, wallH=18, wallRough=3,
                     fans=[[110, 1, 30, 3], [300, -1, 36, 3.5], [440, 1, 46, 4]]),
              dams=[{'x': 235, 'H': 5.0, 'b': 16.0}, {'x': 577, 'H': 9.0, 'b': 22.0}]),
    # 単純な直線水路（以前の既定）
    'S': dict(L=400, B=50, I=20, S=2, grad=0.8, wing=1.5, dx=2.0, fill=0.0, bend=None, V={},
              dams=[{'x': x, 'H': 1.2, 'b': 0.0} for x in (60, 105, 150, 195, 240, 285, 330)]),
}


def preset(name='A', P=None):
    """既定値（または P）にモデル渓流 name の地形と堰堤を入れた条件を返す。"""
    P = copy.deepcopy(P) if P is not None else default_params()
    R = PRESETS[name]
    P['L'] = float(R['L']); P['B'] = float(R['B']); P['slope'] = 1.0 / R['I']; P['crossSlope'] = R['S'] / 100.0
    P['sideGrad'] = R['grad']; P['wingHeight'] = R['wing']; P['dx'] = R['dx']; P['widthTable'] = []
    P['bend'] = {'mode': 'table', 'table': [list(r) for r in R['bend']]} if R['bend'] else {'mode': 'none'}
    V = {'chD': 0.0, 'chW': 0.0, 'chWl': 0.0, 'relief': 0.0, 'wallH': 0.0, 'wallRough': 0.0, 'fans': []}
    V.update(copy.deepcopy(R['V']))
    P['valley'] = V
    P['dams'] = copy.deepcopy(R['dams'])
    P['prefill'] = R['fill']
    P['terrain'] = None
    return P


# =============================================================================
# 乱数（ブラウザ版と同じ値を出す）
# =============================================================================
def _imul(a, b):
    return ((a & M32) * (b & M32)) & M32


def rng(seed):
    """mulberry32（ブラウザ版と同じ列）"""
    st = [((int(seed) & M32) + 0x6D2B79F5) & M32]

    def f():
        a = (st[0] + 0x6D2B79F5) & M32
        st[0] = a
        t = _imul(a ^ (a >> 15), 1 | a)
        t = ((t + _imul(t ^ (t >> 7), 61 | t)) & M32) ^ t
        return ((t ^ (t >> 14)) & M32) / 4294967296.0
    return f


def _hashI(a, b, s):
    """格子点の乱数 -1〜1（NumPy 配列可）"""
    a = np.asarray(a, dtype=np.int64) & M32
    b = np.asarray(b, dtype=np.int64) & M32
    s = np.asarray(s, dtype=np.int64) & M32
    h = ((a * 374761393) & M32) + ((b * 668265263) & M32) + ((s * 2147483647) & M32)
    h &= M32
    h = ((h ^ (h >> 13)) * 1274126177) & M32
    return ((h ^ (h >> 16)) & M32) / 4294967296.0 * 2.0 - 1.0


def vn2(x, y, s):
    """滑らかな乱数場（値ノイズ）-1〜1"""
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    i = np.floor(x); j = np.floor(y)
    fx = x - i; fy = y - j
    fx = fx * fx * (3 - 2 * fx); fy = fy * fy * (3 - 2 * fy)
    i = i.astype(np.int64); j = j.astype(np.int64)
    a = _hashI(i, j, s); b = _hashI(i + 1, j, s); c = _hashI(i, j + 1, s); d = _hashI(i + 1, j + 1, s)
    return (a * (1 - fx) + b * fx) * (1 - fy) + (c * (1 - fx) + d * fx) * fy


def vn1(x, s):
    return 0.7 * vn2(x, 0.5, s) + 0.3 * vn2(np.asarray(x) * 2.7 + 11.3, 3.5, s + 17)


def iwagaki_uc2(dm, sg1):
    """岩垣式（1956）の限界摩擦速度の2乗 [m²/s²]。原式は cm・s 単位、σ/ρ=2.65 の値なので (σ/ρ−1)/1.65 で補正"""
    d = np.asarray(dm, dtype=float) * 100
    u = np.where(d >= 0.303, 80.9 * d,
        np.where(d >= 0.118, 134.6 * d ** (31 / 22),
        np.where(d >= 0.0565, 55.0 * d,
        np.where(d >= 0.0065, 8.41 * d ** (11 / 32), 226 * d))))
    return u * 1e-4 * sg1 / 1.65


def _interp_tab(tab, t, col):
    n = len(tab)
    if t <= tab[0][0]:
        return tab[0][col]
    if t >= tab[n - 1][0]:
        return tab[n - 1][col]
    for i in range(1, n):
        if t <= tab[i][0]:
            w = (t - tab[i - 1][0]) / max(1e-9, tab[i][0] - tab[i - 1][0])
            return tab[i - 1][col] * (1 - w) + tab[i][col] * w
    return tab[n - 1][col]


def _centerline(P, nx, dx):
    """中心線の横ずれ yo(x)（下流を向いて左が正）"""
    yo = np.zeros(nx)
    Bd = P.get('bend') or {}
    x = (np.arange(nx) + 0.5) * dx
    if Bd.get('mode') == 'sine' and Bd.get('amp', 0) > 0 and Bd.get('wl', 0) > 0:
        x0 = Bd.get('x0', 0) or 0
        x1 = Bd['x1'] if Bd.get('x1', 0) > x0 else 1e30
        tp = Bd['wl'] / 2
        a = np.minimum(1, np.minimum((x - x0) / tp, (x1 - x) / tp))
        a = a * a * (3 - 2 * a)
        inside = (x > x0) & (x < x1)
        yo = np.where(inside, Bd['amp'] * a * np.sin(2 * np.pi * (x - x0) / Bd['wl']), 0.0)
    elif Bd.get('mode') == 'table' and Bd.get('table'):
        T = sorted([r for r in Bd['table'] if np.isfinite(r[0]) and np.isfinite(r[1])], key=lambda r: r[0])
        n = len(T)
        for i in range(nx):
            xi = x[i]
            if n == 1 or xi <= T[0][0]:
                yo[i] = T[0][1]; continue
            if xi >= T[n - 1][0]:
                yo[i] = T[n - 1][1]; continue
            k = 1
            while T[k][0] < xi:
                k += 1
            p0 = T[max(0, k - 2)]; p1 = T[k - 1]; p2 = T[k]; p3 = T[min(n - 1, k + 1)]
            hs = p2[0] - p1[0]; t = (xi - p1[0]) / hs
            m1 = (p2[1] - p0[1]) / (p2[0] - p0[0]) * hs if k - 1 > 0 else (p2[1] - p1[1])
            m2 = (p3[1] - p1[1]) / (p3[0] - p1[0]) * hs if k < n - 1 else (p2[1] - p1[1])
            t2 = t * t; t3 = t2 * t
            yo[i] = (2 * t3 - 3 * t2 + 1) * p1[1] + (t3 - 2 * t2 + t) * m1 + (-2 * t3 + 3 * t2) * p2[1] + (t3 - t2) * m2
    return yo


def _sample_terrain(T, X, Y, Ly):
    """読み込んだ地形を (x, y) で双一次補間（rows[0] が左岸側 = y の大きい側）"""
    R = np.asarray(T['rows'], dtype=float); nr, nc = R.shape; cs = T['cs']
    fx = np.clip(X / cs - 0.5, 0, nc - 1); fy = np.clip((Ly - Y) / cs - 0.5, 0, nr - 1)
    i0 = np.minimum(nc - 2, np.floor(fx)).astype(int) if nc >= 2 else np.zeros_like(fx, dtype=int)
    j0 = np.minimum(nr - 2, np.floor(fy)).astype(int) if nr >= 2 else np.zeros_like(fy, dtype=int)
    tx = fx - i0; ty = fy - j0
    i1 = np.minimum(nc - 1, i0 + 1); j1 = np.minimum(nr - 1, j0 + 1)
    a = R[j0, i0] * (1 - tx) + R[j0, i1] * tx
    b = R[j1, i0] * (1 - tx) + R[j1, i1] * tx
    return a * (1 - ty) + b * ty


def _minmod(a, b):
    return np.where(a * b <= 0, 0.0, np.where(np.abs(a) < np.abs(b), a, b))


def _jround(x):
    """JavaScript の Math.round（.5 は切り上げ）"""
    return int(math.floor(x + 0.5))


# =============================================================================
# 計算本体
# =============================================================================
class Sim:
    """平面2次元河床変動計算。配列は (ny, nx)。j は横断方向（j が大きいほど左岸側）、i は流下方向。"""

    def __init__(self, P):
        self.P = P
        self._init()

    # ------------------------------------------------------------------ 初期化
    def _init(self):
        P = self.P
        dx = self.dx = float(P['dx'])
        T = P.get('terrain')
        if T is not None and (len(T['rows']) < 3 or len(T['rows'][0]) < 3):
            T = None
        self.fromCSV = T is not None
        dams_in = P.get('dams') or []
        maxH = max([D['H'] for D in dams_in] + [0])
        self.warn = []
        VL = P.get('valley') or {}
        wallMax = max(VL.get('wallH', 0) if VL.get('wallH', 0) > 0 else 6, maxH + P['wingHeight'] + 3)
        self.hL = self.hR = None
        if T is not None:
            nr, nc = len(T['rows']), len(T['rows'][0])
            nx = max(20, _jround(nc * T['cs'] / dx))
            ny = max(6, _jround(nr * T['cs'] / dx))
            ns = 0
            Bx = np.zeros(nx); yc = np.zeros(nx); cosT = np.ones(nx)
        else:
            nx = max(20, _jround(P['L'] / dx))
            bd = P.get('bend') or {}
            bdT = bd.get('table') if bd.get('mode') == 'table' and bd.get('table') else []
            wt = [[r[0], r[2]] for r in bdT if len(r) > 2 and r[2] is not None and np.isfinite(r[0]) and r[2] > 0]
            if not wt:
                wt = [list(r) for r in (P.get('widthTable') or []) if np.isfinite(r[0]) and r[1] > 0]
            wt.sort(key=lambda r: r[0])
            Bx = np.zeros(nx)
            for i in range(nx):
                bw = (_interp_tab(wt, (i + 0.5) * dx, 1) if len(wt) > 1 else wt[0][1]) if wt else P['B']
                Bx[i] = max(2, _jround(bw / dx)) * dx
            rough = max(0.0, VL.get('wallRough', 0) or 0)
            xr = (np.arange(nx) + 0.5) * dx / 28
            jl = rough * vn1(xr, P['seed'] * 7 + 1) if rough > 0 else 0.0
            jr = rough * vn1(xr, P['seed'] * 7 + 2) if rough > 0 else 0.0
            hL = np.maximum(dx, Bx / 2 + jl); hR = np.maximum(dx, Bx / 2 + jr)
            if rough > 0:
                Bx = hL + hR
            self.hL, self.hR = hL, hR
            yo = _centerline(P, nx, dx)
            ip = np.minimum(nx - 1, np.arange(nx) + 1); im = np.maximum(0, np.arange(nx) - 1)
            sl = (yo[ip] - yo[im]) / ((ip - im) * dx)
            cosT = 1 / np.sqrt(1 + sl * sl)
            smax = np.abs(sl).max()
            if smax > 1.2:
                self.warn.append('中心線の傾きが %.0f° に達します。50° を超える曲がりは精度が落ちます。' % math.degrees(math.atan(smax)))
            ymin = (yo - hR / cosT).min(); ymax = (yo + hL / cosT).max()
            ns = max(1, int(P['sideCells']))
            if VL.get('wallH', 0) > 0:
                ns = max(ns, math.ceil(wallMax / max(0.2, P['sideGrad']) / dx) + 1)
            ny = _jround((ymax - ymin) / dx) + 2 * ns
            yc = ns * dx + (yo - ymin)
        self.nx, self.ny, self.ns = nx, ny, ns
        N = self.N = nx * ny
        self.Lx = nx * dx; self.Ly = ny * dx
        self.Bx, self.yc, self.cosT = Bx, yc, cosT

        # 粒径階
        K = self.K = len(P['grains'])
        self.d = np.array([g['d'] for g in P['grains']], dtype=float)
        sb = sum(max(0, g['bed']) for g in P['grains']) or 1.0
        sf = sum(max(0, g['feed']) for g in P['grains']) or 1.0
        self.f0 = np.array([max(0, g['bed']) / sb for g in P['grains']])
        self.ff = np.array([max(0, g['feed']) / sf for g in P['grains']])
        self.dm0 = float((self.f0 * self.d).sum()); self.dmax = float(self.d.max())
        self.La = P['activeLayer'] if P['activeLayer'] > 0 else max(self.dmax, 0.05)
        self.Ls = max(0.05, P['layerTh'])
        self.tanPhi = math.tan(math.radians(P['reposeDeg']))
        self.sg1 = P['sg'] - 1
        self.n2 = P['manning'] ** 2

        X = (np.arange(nx) + 0.5) * dx
        Y = (np.arange(ny) + 0.5) * dx
        XX, YY = np.meshgrid(X, Y)
        z = np.zeros((ny, nx)); zfix = np.zeros((ny, nx))
        kind = np.zeros((ny, nx), dtype=np.uint8); fl = np.zeros((ny, nx), dtype=bool)
        self.bare = np.zeros((ny, nx), dtype=bool)

        # 初期河床の擾乱（ブラウザ版と同じ乱数列を 3×3 で平滑化）
        rnd = rng(P['seed'])
        tmp = np.array([rnd() * 2 - 1 for _ in range(N)]).reshape(ny, nx)
        pad = np.pad(tmp, 1); cnt = np.pad(np.ones_like(tmp), 1)
        s = sum(pad[1 + a:1 + a + ny, 1 + b:1 + b + nx] for a in (-1, 0, 1) for b in (-1, 0, 1))
        n_ = sum(cnt[1 + a:1 + a + ny, 1 + b:1 + b + nx] for a in (-1, 0, 1) for b in (-1, 0, 1))
        noise = s / n_ * 1.8
        fixedEnd = np.zeros((ny, nx), dtype=bool); fixedEnd[:, :2] = True; fixedEnd[:, nx - 2:] = True

        colmin = np.zeros(nx)
        self.shaped = False
        if T is not None:
            floorH = P['floorH'] if P['floorH'] > 0 else 3
            z = _sample_terrain(T, XX, YY, self.Ly)
            if P.get('csvSlope', 0) > 0:
                mn0 = z.min(axis=0)
                b0 = np.polyfit(X, mn0, 1)[0]
                z = z + (-P['csvSlope'] - b0) * (XX - X.mean())
            colmin = z.min(axis=0)
            yc[:] = (np.argmin(z, axis=0) + 0.5) * dx
            bb = np.polyfit(X, colmin, 1)[0]
            self.trB = bb; self.trA = float(np.mean(colmin) - bb * X.mean())
            self.I = max(1e-3, -bb)
            if bb > 0:
                self.warn.append('読み込んだ地形は下流ほど高くなっています。流下方向を確認してください。')
            M = T.get('mask')
            if M is not None and np.shape(M) == np.shape(T['rows']):
                # 谷底の印（計算後の河床から続けるとき）。初期の擾乱も足さない
                M = np.asarray(M); cs = T['cs']
                fy = np.clip(np.floor((self.Ly - YY) / cs).astype(int), 0, M.shape[0] - 1)
                fx = np.clip(np.floor(XX / cs).astype(int), 0, M.shape[1] - 1)
                fl = M[fy, fx] == 1
            else:
                M = None
                fl = z <= colmin[None, :] + floorH
            kind = np.where(fl, 0, 1).astype(np.uint8)
            if M is None:
                z = np.where(fl & ~fixedEnd, z + P['perturb'] * noise, z)
            zfix = np.where(fl, np.where(fixedEnd, z, z - P['erodDepth']), z)
            Bx[:] = np.maximum(1, fl.sum(axis=0)) * dx
        else:
            I0 = self.I = P['slope']
            self.trA = self.Lx * I0; self.trB = -I0
            chD = max(0.0, VL.get('chD', 0) or 0); relief = max(0.0, VL.get('relief', 0) or 0)
            wR = max(0.0, VL.get('wallRough', 0) or 0)
            fans = [f for f in (VL.get('fans') or []) if np.isfinite(f[0]) and f[2] > 0 and f[3] > 0]
            self.shaped = chD > 0 or relief > 0 or len(fans) > 0 or wR > 0
            hL, hR = self.hL, self.hR
            ych = np.zeros(nx); meanB = Bx.mean()
            chW = VL['chW'] if VL.get('chW', 0) > 0 else 0.3 * meanB
            if chD > 0:
                r2 = rng(P['seed'] * 31 + 7); ph1 = r2() * 6.283; ph2 = r2() * 6.283
                wlc = VL['chWl'] if VL.get('chWl', 0) > 0 else max(60, 4 * meanB)
                sv = 0.72 * np.sin(2 * np.pi * X / wlc + ph1) + 0.28 * np.sin(2 * np.pi * X * 2.37 / wlc + ph2)
                lo = -hR + chW / 2 + dx; hi = hL - chW / 2 - dx
                ych = np.where(hi > lo, (lo + hi) / 2 + (hi - lo) / 2 * 0.85 * sv, (hL - hR) / 2)
                wsm = max(1, _jround(chW / dx))
                for _ in range(2):
                    cs_ = np.concatenate([[0.0], np.cumsum(ych)])
                    a_ = np.maximum(0, np.arange(nx) - wsm); b_ = np.minimum(nx - 1, np.arange(nx) + wsm)
                    ych = (cs_[b_ + 1] - cs_[a_]) / (b_ - a_ + 1)
            self.ych = ych
            cs0 = P['crossSlope']
            zc = (self.Lx - X) * I0

            def plane(o):          # 谷底の面（澪筋へ向かう横断勾配）。o は (ny, nx)
                return zc[None, :] + cs0 * np.abs(o - ych[None, :])
            offs = (YY - yc[None, :]) * cosT[None, :]            # 中心線からの法線方向距離（左が正）
            floor = (offs < hL[None, :]) & (offs > -hR[None, :])
            zz = plane(offs)
            if chD > 0:
                uu = np.abs(offs - ych[None, :]) / (chW / 2)
                zz = zz - chD * np.where(uu < 0.55, 1.0, np.where(uu < 1, 0.5 * (1 + np.cos(np.pi * (uu - 0.55) / 0.45)), 0.0))
            if relief > 0:
                zz = zz + relief * (0.7 * vn2(XX / 16, YY / 9, P['seed'] * 13 + 3) + 0.3 * vn2(XX / 6, YY / 6, P['seed'] * 13 + 4))
            zfloor = zz + np.where(fixedEnd, 0.0, P['perturb'] * noise)
            left = offs >= hL[None, :]
            dw = np.where(left, offs - hL[None, :], -hR[None, :] - offs)
            eo = np.where(left, hL[None, :], -hR[None, :])
            zwall = plane(eo) + np.minimum(P['sideGrad'] * dw, wallMax)
            if wR > 0:
                zwall = zwall + wR * 0.6 * vn2(XX / 11, YY / 11, P['seed'] * 13 + 5) * np.minimum(1, dw / (3 * dx))
            z = np.where(floor, zfloor, zwall)
            fl = floor.copy()
            kind = np.where(floor, 0, 1).astype(np.uint8)
            # 支渓の扇状地（崖錐）
            for Fa in fans:
                fi = max(0, min(nx - 1, _jround(Fa[0] / dx - 0.5))); sd = 1 if Fa[1] >= 0 else -1
                eo2 = hL[fi] if sd > 0 else -hR[fi]
                ax_ = (fi + 0.5) * dx; ay_ = yc[fi] + eo2 / cosT[fi]; RR = Fa[2]; HH = Fa[3]
                rr = np.hypot(XX - ax_, YY - ay_)
                oc = np.clip(offs, -hR[None, :], hL[None, :])
                zf = plane(oc) + HH * np.clip(1 - rr / RR, 0, None) ** 1.5
                hit = (rr < RR) & (zf > z)
                self.bare |= hit & (kind == 1)
                z = np.where(hit, zf, z)
            zfix = np.where(kind == 0, np.where(fixedEnd, z, z - P['erodDepth']), z)
            colmin = zc.copy()
            if self.shaped:
                colmin = np.where(fl.any(axis=0), np.where(fl, z, np.inf).min(axis=0), colmin)
        I = self.I

        # ---- 堰堤（河道の中心線に直角＝両岸に直角）
        self.dams = []
        noSl = np.zeros((ny, nx), dtype=bool)
        if T is not None:
            win = max(3, _jround(15 / dx * 2))
            cs_ = np.concatenate([[0.0], np.cumsum(yc)])
            a_ = np.maximum(0, np.arange(nx) - win); b_ = np.minimum(nx - 1, np.arange(nx) + win)
            ycs = (cs_[b_ + 1] - cs_[a_]) / (b_ - a_ + 1)
        else:
            ycs = yc.copy()
        self.ycs = ycs
        sideG = 0.3 if T is not None else max(0.2, P['sideGrad'])
        for D in sorted(dams_in, key=lambda D: D['x']):
            if not D['H'] > 0:
                continue
            if D['x'] < 0 or D['x'] > self.Lx:
                self.warn.append('x = %g m の堰堤は、地形の範囲（0〜%.0f m）の外なので置きませんでした。' % (D['x'], self.Lx))
                continue
            idd = min(max(_jround(D['x'] / dx - 0.5), 4), nx - 5)
            # 中心 y・回転角・天端標高・堤長・水通しの中心は、与えられていればその値を使う
            fin = lambda k: D.get(k) is not None and np.isfinite(D.get(k))
            cx = (idd + 0.5) * dx; cy = D['y'] if fin('y') else ycs[idd]
            i1 = min(nx - 1, idd + 1); i0 = max(0, idd - 1)
            th = math.radians(D['ang']) if fin('ang') else math.atan2(ycs[i1] - ycs[i0], (i1 - i0) * dx)
            tx, ty = math.cos(th), math.sin(th); nxv, nyv = -ty, tx
            zbase = colmin[idd] if T is not None else (self.Lx - cx) * I
            crest = D['zc'] if fin('zc') else zbase + D['H']; wingTop = crest + P['wingHeight']
            zbase = crest - D['H']
            b = D.get('b', 0) or 0
            spFix = fin('so') and b > 0
            bn = (b if spFix else min(b, Bx[idd])) if b > 0 else Bx[idd]
            hw0 = max(self.hL[idd], self.hR[idd]) if self.hL is not None else Bx[idd] / 2
            ext = D['len'] / 2 if (D.get('len') or 0) > 0 else hw0 + (D['H'] + P['wingHeight']) / sideG + 2 * dx
            thick = (abs(tx) + abs(ty)) * dx / 2 + 1e-9
            px = XX - cx; py = YY - cy
            al = px * tx + py * ty; sl_ = px * nxv + py * nyv
            band = (np.abs(al) <= thick) & (np.abs(sl_) <= ext)
            band[:, :1] = False; band[:, nx - 1:] = False
            iA = max(1, math.floor((cx - ext - dx) / dx)); iB = min(nx - 2, math.ceil((cx + ext + dx) / dx))
            band[:, :iA] = False; band[:, iB + 1:] = False
            inSp = (np.abs(sl_ - D['so']) <= bn / 2 + 1e-9) if spFix else (fl & ((b <= 0) | (np.abs(sl_) <= bn / 2 + 1e-9)))
            top = np.where(inSp, crest, wingTop)
            cells = band & ((top > z) | (inSp if spFix else False))   # 固定した水通しは、地盤が天端より高くても置く
            if not cells.any() or (cells & (kind == 2)).any():
                self.warn.append('x = %g m の堰堤は、ほかの堰堤と重なるので置きませんでした。' % D['x'])
                continue
            z = np.where(cells, top, z); zfix = np.where(cells, top, zfix); kind = np.where(cells, 2, kind).astype(np.uint8)
            nb = cells.copy()
            nb[:, 1:] |= cells[:, :-1]; nb[:, :-1] |= cells[:, 1:]; nb[1:, :] |= cells[:-1, :]; nb[:-1, :] |= cells[1:, :]
            noSl |= nb
            spc = cells & inSp
            sv_ = sl_[cells]; sp_ = sl_[spc]
            spMin, spMax = (sp_.min(), sp_.max()) if sp_.size else (0.0, 0.0)
            self.dams.append(dict(i=idd, x=cx, cx=cx, cy=cy, th=th, tx=tx, ty=ty, nx=nxv, ny=nyv, ext=ext, thick=thick,
                                  H=D['H'], b=bn, crest=crest, wingTop=wingTop, zbase=zbase, zlow=float(z[cells].min()),
                                  cells=cells, spill=spc, sMin=sv_.min() - dx / 2, sMax=sv_.max() + dx / 2,
                                  spMin=spMin - dx / 2, spMax=spMax + dx / 2))
        nd = len(self.dams)
        region = np.full((ny, nx), nd, dtype=np.int16)
        done = np.zeros((ny, nx), dtype=bool)
        for q, Dq in enumerate(self.dams):
            ppx = XX - Dq['cx']; ppy = YY - Dq['cy']
            slq = ppx * Dq['nx'] + ppy * Dq['ny']; alq = ppx * Dq['tx'] + ppy * Dq['ty']
            up = np.where(np.abs(slq) <= Dq['ext'], alq <= Dq['thick'], ppx <= 0) & ~done
            region[up] = q; done |= up
        self.region = region
        # 堰堤の上流をはじめから満砂にする
        if P.get('prefill', 0) > 0 and nd:
            relF = VL.get('relief', 0) if (T is None and VL.get('relief', 0) > 0) else 0
            zt = np.full((ny, nx), -np.inf)
            for q, Dq in enumerate(self.dams):
                m = (region == q) & (kind != 2)
                upd = -((XX - Dq['cx']) * Dq['tx'] + (YY - Dq['cy']) * Dq['ty'])
                tq = Dq['crest'] + upd * I * P['prefill']
                if relF > 0:
                    tq = tq + relF * 0.6 * vn2(XX / 14, YY / 8, P['seed'] * 13 + 9) * np.minimum(1, upd / 6)
                zt = np.where(m & (upd >= 0), tq, zt)
            raise_ = z < zt
            wall2floor = raise_ & (kind == 1)
            kind = np.where(wall2floor, 0, kind).astype(np.uint8); fl = fl | wall2floor
            z = np.where(raise_, zt, z)
            colmin = np.where((kind == 0).any(axis=0), np.where(kind == 0, z, np.inf).min(axis=0), colmin)
            self.shaped = True
        # ---- 侵食できる側岸：谷壁のうち、谷底の高さから erodH までの部分
        BK = P.get('bank') or {}; eH = BK.get('erodH', 0) or 0
        self.bankC = np.zeros((ny, nx), dtype=bool)
        if eH > 0:
            f0m = kind == 0
            with np.errstate(invalid='ignore'):
                fmax = np.where(f0m.any(axis=0), np.where(f0m, z, -np.inf).max(axis=0), np.nan)
                fmin = np.where(f0m.any(axis=0), np.where(f0m, z, np.inf).min(axis=0), np.nan)
            conv = (kind == 1) & np.isfinite(fmax)[None, :] & (z <= fmax[None, :] + eH)
            kind = np.where(conv, 0, kind).astype(np.uint8)
            self.bankC = conv
            zfb = np.where(fixedEnd, z, np.minimum(z, fmin[None, :] - P['erodDepth']))
            zfix = np.where(conv, zfb, zfix)
        # ---- 巨礫
        BL = P.get('boulders') or {}; bL = BL.get('L', 0) or 0
        self.fixedC = np.zeros((ny, nx), dtype=bool); self.nBoulder = 0
        if bL > 0:
            bsz = max(1, min(4, int(BL.get('size', 1) or 1)))
            useCh = (T is None) and (VL.get('chD', 0) or 0) > 0 and getattr(self, 'ych', None) is not None
            if useCh:
                meanB_ = Bx.mean(); chW2 = (VL['chW'] if VL.get('chW', 0) > 0 else 0.3 * meanB_) / 2

            def edge_at(i_, side):
                if useCh:
                    yt = yc[i_] + (self.ych[i_] + side * chW2) / cosT[i_]; jj = _jround(yt / dx - 0.5)
                    return jj if (0 <= jj < ny and kind[jj, i_] == 0) else -1
                col = np.nonzero((kind[:, i_] == 0) & fl[:, i_])[0]
                if not col.size:
                    return -1
                return int(col.max() if side > 0 else col.min())
            for side in (-1, 1):
                x0b = bL / 2 + (bL / 2 if (BL.get('mode') == 'zigzag' and side < 0) else 0)
                xb = x0b
                while xb < self.Lx - 3 * dx:
                    ib = _jround(xb / dx - 0.5)
                    if 3 <= ib <= nx - 4:
                        jb = edge_at(ib, side)
                        if jb >= 0:
                            for di in range(bsz):
                                for dj in range(bsz):
                                    ii = ib + di - (bsz >> 1); jj2 = jb + side * dj
                                    if 1 <= ii < nx - 1 and 0 <= jj2 < ny and kind[jj2, ii] == 0:
                                        self.fixedC[jj2, ii] = True; zfix[jj2, ii] = z[jj2, ii]
                    xb += bL
            self.nBoulder = int(self.fixedC.sum())
        self.z, self.zfix, self.kind, self.fl, self.noSl, self.colmin = z, zfix, kind, fl, noSl, colmin
        self.zmid = float(z[fl].mean()) if fl.any() else 0.0
        self.zlo = float(z.min()); self.zhi = float(z.max())
        self.Bf = Bx[0]
        self.z0 = z.copy()

        # ---- 土砂の層構造（A：交換層内の粒径別体積、SF：堆積層の粒径別割合）
        La, Ls = self.La, self.Ls
        maxH = max([D['H'] for D in self.dams] + [0])
        thMax = float((z - zfix).max())
        NL = self.NL = min(100, max(math.ceil((P['erodDepth'] + maxH + 3) / Ls) + 2, math.ceil(thMax / Ls) + 2))
        zf = z.ravel(); zx = zfix.ravel()
        dl = zf - zx
        self.A = np.zeros((N, K))
        a = np.where(dl > 1e-9, np.minimum(dl, La), 0.0)
        self.A[:] = a[:, None] * self.f0[None, :]
        rem = np.where(dl > 1e-9, dl - a, 0.0)
        nl = np.where(rem > 1e-9, np.minimum(NL, np.ceil(rem / Ls - 1e-9)), 0).astype(np.int32)
        self.nl = nl
        self.topTh = np.where(nl > 0, rem - (nl - 1) * Ls, 0.0)
        self.SF = np.zeros((N, NL, K), dtype=np.float32)
        lay = np.arange(NL)[None, :] < nl[:, None]
        self.SF[lay] = self.f0.astype(np.float32)
        self.tcs = np.full((ny, nx), np.nan)
        self.Cc = np.zeros((N, K))          # 侵食・堆積速度式：流れの中の粒径別土砂量 C_k h（実質体積の層厚） [m]
        self.dzs = np.zeros((ny, nx))       # 1 回の河床変動計算での河床変化
        self.kap = np.zeros((ny, nx))       # 流線の曲率 1/r（左曲がりが正） [1/m]
        self.nHase = 0                      # 長谷川モデルの崩落回数

        # ---- 初期の水
        h = np.zeros((ny, nx)); qx = np.zeros((ny, nx)); qy = np.zeros((ny, nx))
        Q0 = self.Qat(0)
        if Q0 > 0.01:
            q0 = Q0 / Bx
            hn = (q0 * P['manning'] / math.sqrt(max(I, 1e-4))) ** 0.6; un = q0 / hn
            if T is not None or self.shaped:
                hw = colmin[None, :] + hn[None, :] - z
                m = (kind == 0) & (hw > 0)
                h = np.where(m, hw, 0.0); qx = np.where(m, hw * un[None, :], 0.0)
            else:
                m = kind == 0
                h = np.where(m, hn[None, :], 0.0); qx = np.where(m, q0[None, :], 0.0)
        if nd:   # 堰堤の上流を天端まで湛水させておく
            cr = np.array([D['crest'] for D in self.dams] + [-np.inf])[region]
            m = (region < nd) & (kind != 2) & (z < cr) & (cr - z > h)
            h = np.where(m, cr - z, h); qx = np.where(m, 0.0, qx); qy = np.where(m, 0.0, qy)
        self.h, self.qx, self.qy = h, qx, qy

        # ---- 時刻・収支
        self.t = 0.0; self.th = 0.0; self.nstep = 0; self.dt = 0.01
        self.phase = 'spinup' if Q0 > 0.01 else 'run'
        self.spinT = 0.0
        self.Qin = Q0; self.Qout = 0.0; self.QsIn = 0.0; self.QsOut = 0.0
        self.Vwin = 0.0; self.Vwout = 0.0; self.Vw0 = self.water_volume()
        self.sedIn = 0.0; self.sedOut = 0.0; self.sedOutK = np.zeros(K)
        self.nLimW = 0; self.nLimS = 0
        self.series = {'t': [], 'Qin': [], 'Qout': [], 'QsIn': [], 'QsOut': [], 'stor': []}
        self._nextSample = 0.0
        self.Uref = max(0.5, (Q0 / self.Bf) / ((Q0 / self.Bf) * P['manning'] / math.sqrt(max(I, 1e-4))) ** 0.6) if Q0 > 0.01 else 1.0
        self.FxM = np.zeros((ny, nx + 1)); self.GyM = np.zeros((ny + 1, nx))
        self.qin = np.zeros(ny); self.uin = np.zeros(ny)

    # ------------------------------------------------------------ ハイドログラフ
    def duration(self):
        H = self.P['hydro']
        if H['mode'] == 'table' and len(H.get('table') or []) > 1:
            return H['table'][-1][0]
        return H['duration']

    def Qat(self, t):
        H = self.P['hydro']
        if H['mode'] == 'table' and len(H.get('table') or []) > 1:
            return max(0.0, _interp_tab(H['table'], t, 1))
        if H['mode'] == 'const':
            return H['Qpeak']
        if t <= 0:
            return H['Qbase']
        if t < H['tPeak']:
            return H['Qbase'] + (H['Qpeak'] - H['Qbase']) * t / H['tPeak']
        if t < H['duration']:
            return H['Qpeak'] + (H['Qbase'] - H['Qpeak']) * (t - H['tPeak']) / max(1e-9, H['duration'] - H['tPeak'])
        return H['Qbase']

    def QsAt(self, t):
        """供給土砂量（実質体積） [m3/s]"""
        H = self.P['hydro']
        tab = H.get('table') or []
        if H['mode'] == 'table' and len(tab) > 1 and len(tab[0]) > 2 and np.isfinite(tab[0][2]):
            return max(0.0, _interp_tab(tab, t, 2))
        return H['Cs'] * self.Qat(t)

    # ------------------------------------------------------------ 1 ステップ
    def step(self):
        P = self.P; dx = self.dx; nx, ny = self.nx, self.ny
        h, qx, qy, z = self.h, self.qx, self.qy, self.z
        HDRY, HS, UMAX, g = 1e-3, 0.01, 30.0, G
        hg = 0.5 * g
        second = P['order'] != 1

        # ---- 流入境界（上流端の谷底セルに流量を配分。水深の 5/3 乗で重み付け）
        Q = self.Qat(self.t); self.Qin = Q
        sqI = math.sqrt(max(self.I, 1e-4)); nman = P['manning']
        w = np.where(self.fl[:, 0], (h[:, 0] + 0.05) ** (5 / 3), 0.0)
        qin = Q * w / w.sum() / dx if w.sum() > 0 else np.zeros(ny)
        uin = np.where(qin > 0, qin / np.maximum(qin * nman / sqI, 1e-30) ** 0.6, 0.0)
        self.qin, self.uin = qin, uin

        # ---- 2次流：流線の曲率（前のステップの面の流束から求める）
        SEC = P.get('sec') or {}
        if SEC.get('bed') or SEC.get('mom'):
            self._curvature()

        # ---- 原始変数
        wet = h > HDRY
        hs = np.where(wet, h, 1.0)
        u = np.clip(np.where(wet, qx / hs, 0.0), -UMAX, UMAX)
        v = np.clip(np.where(wet, qy / hs, 0.0), -UMAX, UMAX)
        qx = np.where(wet, u * h, 0.0); qy = np.where(wet, v * h, 0.0)
        eta = z + h
        dhx = np.zeros_like(h); dex = np.zeros_like(h); dhy = np.zeros_like(h); dey = np.zeros_like(h)
        if second:
            ok = (h > HS) & ~self.noSl
            okx = ok[:, 1:-1] & (h[:, :-2] > HS) & (h[:, 2:] > HS)
            dhx[:, 1:-1] = np.where(okx, _minmod(h[:, 1:-1] - h[:, :-2], h[:, 2:] - h[:, 1:-1]), 0.0)
            dex[:, 1:-1] = np.where(okx, _minmod(eta[:, 1:-1] - eta[:, :-2], eta[:, 2:] - eta[:, 1:-1]), 0.0)
            oke = ok[:, -1] & (h[:, -2] > HS)            # 下流端：初期勾配で続く仮想セル
            dex[:, -1] = np.where(oke, _minmod(eta[:, -1] - eta[:, -2], -self.I * dx), 0.0)
            oky = ok[1:-1, :] & (h[:-2, :] > HS) & (h[2:, :] > HS)
            dhy[1:-1, :] = np.where(oky, _minmod(h[1:-1, :] - h[:-2, :], h[2:, :] - h[1:-1, :]), 0.0)
            dey[1:-1, :] = np.where(oky, _minmod(eta[1:-1, :] - eta[:-2, :], eta[2:, :] - eta[1:-1, :]), 0.0)

        # ---- 数値流束（HLL + 静水圧再構築）。axis=1 が x、axis=0 が y
        def hll(hL, hR, eL, eR, uL, uR, tL, tR):
            zL = eL - hL; zR = eR - hR
            zs = np.minimum(np.maximum(zL, zR), np.minimum(eL, eR))
            hLs = np.minimum(eL - zs, hL); hRs = np.minimum(eR - zs, hR)
            dryb = (hL < 1e-9) & (hR < 1e-9)
            drys = (hLs < 1e-9) & (hRs < 1e-9)
            cL = np.sqrt(g * np.maximum(hLs, 0)); cR = np.sqrt(g * np.maximum(hRs, 0))
            Ld = hLs < 1e-9; Rd = hRs < 1e-9
            SL = np.where(Ld, uR - 2 * cR, np.where(Rd, uL - cL, np.minimum(uL - cL, uR - cR)))
            SR = np.where(Ld, uR + cR, np.where(Rd, uL + 2 * cL, np.maximum(uL + cL, uR + cR)))
            mL = hLs * uL; mR = hRs * uR
            fL1 = mL * uL + hg * hLs * hLs; fR1 = mR * uR + hg * hRs * hRs
            den = SR - SL; den = np.where(np.abs(den) < 1e-12, 1e-12, den)
            F0h = (SR * mL - SL * mR + SL * SR * (hRs - hLs)) / den
            F1h = (SR * fL1 - SL * fR1 + SL * SR * (mR - mL)) / den
            F0 = np.where(SL >= 0, mL, np.where(SR <= 0, mR, F0h))
            F1 = np.where(SL >= 0, fL1, np.where(SR <= 0, fR1, F1h))
            F0 = np.where(drys | dryb, 0.0, F0); F1 = np.where(drys | dryb, 0.0, F1)
            sa = np.where(drys | dryb, 0.0, np.maximum(-SL, SR))
            Ft = F0 * np.where(F0 >= 0, tL, tR)
            CL = np.where(dryb, 0.0, hg * (hL + hLs) * (zs - zL))
            CR = np.where(dryb, 0.0, hg * (hR + hRs) * (zs - zR))
            return F0, F1, Ft, CL, CR, (sa.max() if sa.size else 0.0)

        FxM = np.zeros((ny, nx + 1)); FxU = np.zeros_like(FxM); FxV = np.zeros_like(FxM)
        FxCL = np.zeros_like(FxM); FxCR = np.zeros_like(FxM)
        F0, F1, Ft, CL, CR, s1 = hll(h[:, :-1] + 0.5 * dhx[:, :-1], h[:, 1:] - 0.5 * dhx[:, 1:],
                                      eta[:, :-1] + 0.5 * dex[:, :-1], eta[:, 1:] - 0.5 * dex[:, 1:],
                                      u[:, :-1], u[:, 1:], v[:, :-1], v[:, 1:])
        FxM[:, 1:nx] = F0; FxU[:, 1:nx] = F1; FxV[:, 1:nx] = Ft; FxCL[:, 1:nx] = CL; FxCR[:, 1:nx] = CR
        FxCR[:, 0] = hg * h[:, 0] ** 2                       # 上流端：壁（流入は生成項）
        hb = h[:, -1]; ub = u[:, -1]                         # 下流端：自由流出
        out = (hb > HDRY) & (ub > 0)
        FxM[:, nx] = np.where(out, hb * ub, 0.0)
        FxU[:, nx] = np.where(out, hb * ub * ub + hg * hb * hb, hg * hb * hb)
        FxV[:, nx] = FxM[:, nx] * v[:, -1]
        s2 = (ub + np.sqrt(g * np.maximum(hb, 0)))[out].max() if out.any() else 0.0
        GyM = np.zeros((ny + 1, nx)); GyU = np.zeros_like(GyM); GyV = np.zeros_like(GyM)
        GyCL = np.zeros_like(GyM); GyCR = np.zeros_like(GyM)
        F0, F1, Ft, CL, CR, s3 = hll(h[:-1, :] + 0.5 * dhy[:-1, :], h[1:, :] - 0.5 * dhy[1:, :],
                                      eta[:-1, :] + 0.5 * dey[:-1, :], eta[1:, :] - 0.5 * dey[1:, :],
                                      v[:-1, :], v[1:, :], u[:-1, :], u[1:, :])
        GyM[1:ny] = F0; GyV[1:ny] = F1; GyU[1:ny] = Ft; GyCL[1:ny] = CL; GyCR[1:ny] = CR
        GyCR[0] = hg * h[0] ** 2; GyCL[ny] = hg * h[-1] ** 2   # 側方：壁
        smax = max(1e-6, s1, s2, s3)

        # ---- 時間刻み（CFL）
        dt = min(P['cfl'] * dx / smax, 0.5)
        if self.dt > 0:
            dt = min(dt, 1.5 * self.dt)
        self.dt = dt
        cf = dt / dx

        # ---- 水深の非負性を保つ流出制限
        outv = cf * (np.maximum(FxM[:, 1:], 0) - np.minimum(FxM[:, :-1], 0) + np.maximum(GyM[1:], 0) - np.minimum(GyM[:-1], 0))
        sc = np.where(outv > h, np.where(outv > 0, h / np.where(outv > 0, outv, 1), 0.0), 1.0)
        if (sc < 1).any():
            self.nLimW += 1
            r = np.ones((ny, nx + 1))
            r[:, 1:] = np.where(FxM[:, 1:] > 0, sc, np.concatenate([sc[:, 1:], np.ones((ny, 1))], axis=1))
            r[:, 0] = 1.0
            FxM *= r; FxU *= r; FxV *= r
            ry = np.ones((ny + 1, nx))
            ry[1:ny] = np.where(GyM[1:ny] > 0, sc[:-1], sc[1:])
            GyM *= ry; GyU *= ry; GyV *= ry

        # ---- 更新（勾配項・流入・摩擦）
        hn = h - cf * (FxM[:, 1:] - FxM[:, :-1] + GyM[1:] - GyM[:-1])
        qxn = qx - cf * ((FxU[:, 1:] + FxCL[:, 1:]) - (FxU[:, :-1] + FxCR[:, :-1]) + GyU[1:] - GyU[:-1])
        qyn = qy - cf * (FxV[:, 1:] - FxV[:, :-1] + (GyV[1:] + GyCL[1:]) - (GyV[:-1] + GyCR[:-1]))
        if second:
            pos = h > 0
            qxn = qxn - np.where(pos, cf * g * h * (dex - dhx), 0.0)
            qyn = qyn - np.where(pos, cf * g * h * (dey - dhy), 0.0)
        if SEC.get('mom'):
            # 2次流による運動量の分散項（面の値は両側のセル中心の平均。どちらかが乾いていれば 0）
            Txx, Txy, msk = self._sec_stress()
            def face(T):
                E = np.zeros_like(T); W = np.zeros_like(T); Nf = np.zeros_like(T); Sf = np.zeros_like(T)
                E[:, :-1] = np.where(msk[:, 1:], 0.5 * (T[:, :-1] + T[:, 1:]), 0.0)
                W[:, 1:] = np.where(msk[:, :-1], 0.5 * (T[:, 1:] + T[:, :-1]), 0.0)
                Nf[:-1] = np.where(msk[1:], 0.5 * (T[:-1] + T[1:]), 0.0)
                Sf[1:] = np.where(msk[:-1], 0.5 * (T[1:] + T[:-1]), 0.0)
                return E, W, Nf, Sf
            tE, tW, tN, tS = face(Txx); sE, sW, sN, sS = face(Txy)
            qxn = qxn - np.where(msk, cf * ((tE - tW) + (sN - sS)), 0.0)
            qyn = qyn - np.where(msk, cf * ((sE - sW) - (tN - tS)), 0.0)
        hn[:, 0] += np.where(qin > 0, cf * qin, 0.0)
        qxn[:, 0] += np.where(qin > 0, cf * qin * uin, 0.0)
        hn = np.maximum(hn, 0.0)
        wetn = hn > HDRY
        hs = np.where(wetn, hn, 1.0)
        Vv = np.sqrt(qxn * qxn + qyn * qyn)
        fr = 0.5 * (1 + np.sqrt(1 + 4 * dt * g * self.n2 * Vv / (hs * hs * np.cbrt(hs))))
        self.h = hn
        self.qx = np.where(wetn, qxn / fr, 0.0)
        self.qy = np.where(wetn, qyn / fr, 0.0)
        self.u, self.v = u, v
        self.FxM, self.GyM = FxM, GyM
        self.Qout = FxM[:, nx].sum() * dx
        self.Vwin += Q * dt; self.Vwout += self.Qout * dt
        self.th += dt; self.nstep += 1

        # ---- 助走計算（水だけを流して流れを立ち上げる）
        if self.phase == 'spinup':
            self.spinT += dt
            tmin = 0.5 * self.Lx / self.Uref; tmax = 3 * self.Lx / self.Uref
            if (self.spinT > tmin and abs(self.Qout - Q) < 0.08 * Q) or self.spinT > tmax:
                self.phase = 'run'; self.Vwin = 0.0; self.Vwout = 0.0; self.Vw0 = self.water_volume()
            return
        if self.phase == 'done':
            return

        # ---- 河床変動
        dtm = dt * max(1.0, P['morfac'])
        self._morph(dtm)
        self.t += dtm
        if self.t >= self._nextSample:
            S = self.series
            S['t'].append(self.t); S['Qin'].append(Q); S['Qout'].append(self.Qout)
            S['QsIn'].append(self.QsIn); S['QsOut'].append(self.QsOut); S['stor'].append(self.storage())
            self._nextSample = self.t + self.duration() / 400
        if self.t >= self.duration():
            self.phase = 'done'

    # ------------------------------------------------------------ 限界掃流力
    def _flow_slope(self, ux, uy):
        """流れの向きの河床勾配 i（下り勾配を正、0〜0.2）。堰堤のセルは除く"""
        z, kind, dx = self.z, self.kind, self.dx
        nd = kind != 2
        zl = np.empty_like(z); zl[:, 1:] = z[:, :-1]; zl[:, 0] = z[:, 0]
        zr = np.empty_like(z); zr[:, :-1] = z[:, 1:]; zr[:, -1] = z[:, -1]
        al = np.zeros_like(nd); al[:, 1:] = nd[:, :-1]
        ar = np.zeros_like(nd); ar[:, :-1] = nd[:, 1:]
        gx = np.where(al & ar, (zr - zl) / (2 * dx), np.where(al, (z - zl) / dx, np.where(ar, (zr - z) / dx, -self.I)))
        zd = np.empty_like(z); zd[1:] = z[:-1]; zd[0] = z[0]
        zu = np.empty_like(z); zu[:-1] = z[1:]; zu[-1] = z[-1]
        ad = np.zeros_like(nd); ad[1:] = nd[:-1]
        au = np.zeros_like(nd); au[:-1] = nd[1:]
        gy = np.where(ad & au, (zu - zd) / (2 * dx), np.where(ad, (z - zd) / dx, np.where(au, (zu - z) / dx, 0.0)))
        V = np.sqrt(ux * ux + uy * uy)
        sl = np.where(V > 1e-9, -(gx * ux + gy * uy) / np.where(V > 1e-9, V, 1), 0.0)
        return np.clip(sl, 0.0, 0.2)

    def _crit_us2(self, dm, ux, uy):
        """平均粒径 dm の限界摩擦速度の2乗 u*cm² [m²/s²]"""
        P = self.P
        if P['tcMethod'] == 'iwagaki':
            return iwagaki_uc2(dm, self.sg1)
        if P['tcMethod'] == 'const':
            return P['tauCm'] * self.sg1 * G * dm
        return 0.04 * 10 ** (1.72 * self._flow_slope(ux, uy)) * self.sg1 * G * dm   # 芦田・高橋・水山（1978）式(24)

    # ------------------------------------------------------------ 河床変動 1 ステップ
    def _morph(self, dtm):
        self.dzs = np.zeros((self.ny, self.nx))
        if self.P.get('transport') == 'noneq':
            nsub = max(1, math.ceil(dtm / max(self.dt, 1e-6) - 1e-9))
            for _ in range(nsub):
                self._noneq_step(dtm / nsub)
        else:
            self._morph_eq(dtm)
        self._dry_bank_erosion()
        self._hasegawa_collapse()
        self._collapse()

    def _morph_eq(self, dtm):
        """平衡流砂量（芦田・道上式）による河床変動"""
        P = self.P; nx, ny, K, dx = self.nx, self.ny, self.K, self.dx
        h, z, zfix = self.h, self.z, self.zfix
        A = self.A.reshape(ny, nx, K); d = self.d
        g = G; sg1 = self.sg1; La = self.La; tanPhi = self.tanPhi
        HSED = 0.02
        invMu = (1.0 if P.get('gravK') is None else P['gravK']) / math.sqrt(1.0 * 0.45)
        cGrid = 0.3 if P.get('gridDamp') is None else P['gridDamp']
        FxM, GyM = self.FxM, self.GyM

        # ---- 粒径別掃流砂量（芦田・道上式）。流速は面の質量流束の平均から求める
        As = A.sum(axis=2)
        hs = np.where(h >= HSED, h, 1.0)
        uu = 0.5 * (FxM[:, :-1] + FxM[:, 1:]) / hs
        vv = 0.5 * (GyM[:-1] + GyM[1:]) / hs
        V2 = uu * uu + vv * vv
        base = (h >= HSED) & (As >= 1e-6) & (V2 >= 1e-4) & ~self.fixedC     # 巨礫のセルは動かない
        V = np.sqrt(np.where(base, V2, 1.0))
        us2 = g * self.n2 * V2 / np.cbrt(hs)
        dm = np.where(base, (A * d).sum(axis=2) / np.where(As > 0, As, 1), self.dm0)
        taum = us2 / (sg1 * g * dm)
        arg = hs / (dm * (1 + 2 * taum))
        den = 6.0 + 5.75 * np.log10(np.maximum(arg, 0.2))
        ue = V / np.maximum(den, 2.0)
        ue2 = ue * ue
        lim = ue2 > us2
        ue2 = np.where(lim, us2, ue2); ue = np.where(lim, np.sqrt(us2), ue)
        uss = np.sqrt(np.where(base, us2, 1.0))
        cap = 17 * ue2 * ue / (sg1 * g) / La
        uscm2 = self._crit_us2(dm, uu, vv)
        self.tcs = np.where(base, uscm2 / (sg1 * g * dm), np.nan)
        epsMin = np.minimum(cGrid * dx / hs, 4.0)
        qbx = np.zeros((ny, nx, K)); qby = np.zeros((ny, nx, K)); gk = np.zeros((ny, nx, K))
        act = np.zeros((ny, nx), dtype=bool)
        hideEg = P['hiding'] == 'egiazaroff'
        # 2次流で河床付近の流れが曲率の内側へ向く（長田ら 1996 式14：u_nb/u_sb = N* h/r。n は流れの左向き）
        SEC = P.get('sec') or {}
        secD = (SEC.get('Nstar') or 11.5) * np.clip(hs * self.kap, -(SEC.get('hrMax') or 0.05), SEC.get('hrMax') or 0.05) if SEC.get('bed') else 0.0
        for k in range(K):
            rr = d[k] / dm
            if hideEg:
                lg = LN19 / np.log(19 * np.maximum(rr, 1e-9))
                usc2 = np.where(rr >= 0.4, uscm2 * lg * lg * rr, 0.85 * uscm2)
            else:
                usc2 = uscm2
            m = base & (us2 > usc2) & (A[:, :, k] > 0)
            usc = np.sqrt(usc2)
            qb = np.where(m, cap * A[:, :, k] * (1 - usc2 / np.where(m, us2, 1)) * (1 - usc / uss), 0.0)
            qbx[:, :, k] = qb * (uu - secD * vv) / V; qby[:, :, k] = qb * (vv + secD * uu) / V
            epsg = np.maximum(usc / uss * invMu, epsMin)
            gk[:, :, k] = qb * epsg
            act |= m
        qbx[~act] = 0; qby[~act] = 0; gk[~act] = 0

        # ---- 面の流砂量（上流化＋勾配による重力項）
        SFx = np.zeros((ny, nx + 1, K)); SFy = np.zeros((ny + 1, nx, K))
        aL = act[:, :-1]; aR = act[:, 1:]
        S = np.clip((z[:, 1:] - z[:, :-1]) / dx, -tanPhi, tanPhi)
        bLR = (z[:, 1:] - zfix[:, 1:] < 1e-6) & (z[:, 1:] - z[:, :-1] > 0.5 * h[:, :-1])   # 固定床への段上がりは通さない
        bRL = (z[:, :-1] - zfix[:, :-1] < 1e-6) & (z[:, :-1] - z[:, 1:] > 0.5 * h[:, 1:])
        qL = qbx[:, :-1]; qR = qbx[:, 1:]
        fa = np.where((aL & ~bLR)[..., None] & (qL > 0), qL, 0) + np.where((aR & ~bRL)[..., None] & (qR < 0), qR, 0)
        gg = np.where(aL[..., None], gk[:, :-1], 0) + np.where(aR[..., None], gk[:, 1:], 0)
        SFx[:, 1:nx] = np.where((self.fixedC[:, :-1] | self.fixedC[:, 1:])[..., None], 0.0, fa - 0.5 * gg * S[..., None])
        SFx[:, nx] = np.where(act[:, -1][:, None] & (qbx[:, -1] > 0), qbx[:, -1], 0)       # 下流端：流出
        aL = act[:-1]; aR = act[1:]
        S = np.clip((z[1:] - z[:-1]) / dx, -tanPhi, tanPhi)
        bLR = (z[1:] - zfix[1:] < 1e-6) & (z[1:] - z[:-1] > 0.5 * h[:-1])
        bRL = (z[:-1] - zfix[:-1] < 1e-6) & (z[:-1] - z[1:] > 0.5 * h[1:])
        qL = qby[:-1]; qR = qby[1:]
        fa = np.where((aL & ~bLR)[..., None] & (qL > 0), qL, 0) + np.where((aR & ~bRL)[..., None] & (qR < 0), qR, 0)
        gg = np.where(aL[..., None], gk[:-1], 0) + np.where(aR[..., None], gk[1:], 0)
        SFy[1:ny] = np.where((self.fixedC[:-1] | self.fixedC[1:])[..., None], 0.0, fa - 0.5 * gg * S[..., None])

        # ---- 粒径別の持ち出し制限（交換層内の量を超えない）
        cf = dtm / ((1 - P['porosity']) * dx)
        outv = cf * (np.maximum(SFx[:, 1:], 0) - np.minimum(SFx[:, :-1], 0) + np.maximum(SFy[1:], 0) - np.minimum(SFy[:-1], 0))
        scS = np.where(outv > A, np.where(outv > 0, A / np.where(outv > 0, outv, 1), 0.0), 1.0)
        if (scS < 1).any():
            self.nLimS += 1
            r = np.ones_like(SFx)
            r[:, 1:] = np.where(SFx[:, 1:] > 0, scS, np.concatenate([scS[:, 1:], np.ones((ny, 1, K))], axis=1))
            SFx *= r
            ry = np.ones_like(SFy)
            ry[1:ny] = np.where(SFy[1:ny] > 0, scS[:-1], scS[1:])
            SFy *= ry

        # ---- 上流端の給砂
        Qs = self.QsAt(self.t) if self.phase == 'run' else 0.0
        self.QsIn = Qs
        Qw = self.Qin
        feedCoef = Qs / Qw * dtm / ((1 - P['porosity']) * dx) if (Qs > 0 and Qw > 0) else 0.0
        self.sedIn += Qs * dtm

        # ---- Exner 式
        out_k = SFx[:, nx].sum(axis=0)
        self.sedOutK += out_k * dx * dtm
        dA = -cf * (SFx[:, 1:] - SFx[:, :-1] + SFy[1:] - SFy[:-1])
        if feedCoef > 0:
            dA[:, 0, :] += np.where(self.qin > 0, feedCoef * self.qin, 0.0)[:, None] * self.ff[None, :]
        An = A + dA
        neg = An < 0
        dA = np.where(neg, -A, dA); An = np.where(neg, 0.0, An)
        dz = dA.sum(axis=2)
        self.A = An.reshape(-1, K)
        self.z = z + dz
        self.dzs += dz
        chg = np.flatnonzero(dz.ravel() != 0)
        if chg.size:
            self._fix_layers(chg)
        self.QsOut = out_k.sum() * dx
        self.sedOut += self.QsOut * dtm

    def _noneq_step(self, dts):
        """侵食・堆積速度式による非平衡の河床変動（高橋・里深 1999、高橋・中川・里深・奥村・安本 1998）
        i_k = f_bk δe (C∞−C)/C* |u*−u*ck|（C∞ > C かつ u* > u*ck）、f_0k δd (C∞−C)/C* |u*−u*ck|（C∞ ≤ C かつ u* < u*ck）
        ∂(C_k h)/∂t + ∇·(C_k h u') = C* i_k、∂z/∂t = −Σ i_k。u' は長谷川（1983）の横断方向流砂で向きを補正した流速"""
        P = self.P; NE = P.get('noneq') or {}; nx, ny, K, dx = self.nx, self.ny, self.K, self.dx
        h = self.h; z = self.z; d = self.d; g = G; sg1 = self.sg1; Cst = 1 - P['porosity']; HSED = 0.02
        A = self.A.reshape(ny, nx, K); Cc = self.Cc.reshape(ny, nx, K)
        de0 = NE.get('de') or 0.6; dd0 = NE.get('dd') or 0.6; hE = NE.get('hEdge') or 2.5 * self.dm0
        fe = NE.get('fe') or 1.0; fd = NE.get('fd') or 1.0; gK = (1.0 if P.get('gravK') is None else P['gravK']) / math.sqrt(1.0 * 0.45)
        dz = np.zeros((ny, nx))
        # 水が引いたセル：流れの中の土砂をその場に堆積させる
        dry = h < HSED
        if (dry[..., None] & (Cc > 0)).any():
            dep = np.where(dry[..., None], Cc / Cst, 0.0)
            A = A + dep; dz += dep.sum(axis=2); Cc = np.where(dry[..., None], 0.0, Cc)
        wet = ~dry
        hs = np.where(wet, h, 1.0)
        uu = 0.5 * (self.FxM[:, :-1] + self.FxM[:, 1:]) / hs; vv = 0.5 * (self.GyM[:-1] + self.GyM[1:]) / hs
        V = np.sqrt(uu * uu + vv * vv)
        As = A.sum(axis=2); hasA = As > 1e-9
        fb = np.where(hasA[..., None], A / np.where(hasA, As, 1)[..., None], 0.0)
        dm = np.where(hasA, (fb * d).sum(axis=2), self.dm0)
        order = np.argsort(d); cum = np.cumsum(fb[..., order], axis=2)
        ks = np.where(hasA, d[order][np.minimum(K - 1, (cum < 0.9 - 1e-9).sum(axis=2))], d.max())
        den = np.maximum(6.0 + 5.75 * np.log10(np.maximum(hs / ks, 0.2)), 2.0)
        us = V / den; us2 = us * us
        uscm2 = self._crit_us2(dm, uu, vv)
        self.tcs = np.where(wet, uscm2 / (sg1 * g * dm), np.nan)
        # 輸送速度（横断方向の流砂による向きの補正）
        b1 = np.arctan2(vv, uu); nxn = -np.sin(b1); nyn = np.cos(b1)
        gy, gx = np.gradient(z, dx)
        dzdn = gx * nxn + gy * nyn
        Kc = np.where(us2 > 1e-12, np.sqrt(np.minimum(1.0, uscm2 / np.where(us2 > 1e-12, us2, 1))) * gK, 0.0)
        SEC = P.get('sec') or {}       # 2次流による河床付近の流向の偏り（長田ら 1996 式11・14）
        sd0 = (SEC.get('Nstar') or 11.5) * np.clip(hs * self.kap, -(SEC.get('hrMax') or 0.05), SEC.get('hrMax') or 0.05) if SEC.get('bed') else 0.0
        b2 = np.clip(np.arctan(sd0 - Kc * dzdn), -1.0, 1.0)
        mv_ = wet & (V > 1e-6)
        tu = np.where(mv_, V * np.cos(b1 + b2), 0.0); tv = np.where(mv_, V * np.sin(b1 + b2), 0.0)
        # 侵食・堆積速度
        q = hs * V; Ct = Cc.sum(axis=2) / hs
        hideEg = P['hiding'] == 'egiazaroff'
        usck2 = np.zeros((ny, nx, K)); Cinf = np.zeros((ny, nx))
        for k in range(K):
            rr = d[k] / dm
            if hideEg:
                lg = LN19 / np.log(19 * np.maximum(rr, 1e-9))
                usck2[..., k] = np.where(rr >= 0.4, uscm2 * lg * lg * rr, 0.85 * uscm2)
            else:
                usck2[..., k] = uscm2
            okk = wet & hasA & (q > 1e-9) & (us2 > usck2[..., k])
            tk = us2 / (sg1 * g * d[k]); tck = usck2[..., k] / (sg1 * g * d[k])
            qbu = np.where(okk, 17 * math.sqrt(sg1 * g * d[k] ** 3) * np.maximum(tk, 0) ** 1.5 * (1 - tck / np.where(okk, tk, 1)) * (1 - np.sqrt(usck2[..., k]) / np.where(okk, us, 1)), 0.0)
            Cinf += fb[..., k] * qbu / np.where(q > 1e-9, q, 1)
        edge = h < hE
        de = np.where(edge, de0 * fe, de0); ddd = np.where(edge, dd0 * fd, dd0)
        diffC = Cinf - Ct
        usck = np.sqrt(usck2)
        f0 = np.where((Ct > 0)[..., None], Cc / np.where(Ct > 0, Ct * hs, 1)[..., None], 0.0)
        ero = (diffC > 0)[..., None] & (us[..., None] > usck) & hasA[..., None]
        depo = (diffC <= 0)[..., None] & (Ct > 0)[..., None] & (us[..., None] < usck)
        rate = np.where(ero, fb * de[..., None] * diffC[..., None] / Cst * (us[..., None] - usck),
               np.where(depo, f0 * ddd[..., None] * diffC[..., None] / Cst * (usck - us[..., None]), 0.0))
        act = wet & ~self.fixedC
        dzk = np.where(act[..., None], -rate * dts, 0.0)
        dzk = np.where(dzk < 0, np.maximum(dzk, -A), np.minimum(dzk, Cc / Cst))
        A = A + dzk; Cc = Cc - dzk * Cst; dz += dzk.sum(axis=2)
        self.A = A.reshape(-1, K); self.z = z + dz; self.dzs += dz
        chg = np.flatnonzero(dz.ravel() != 0)
        if chg.size:
            self._fix_layers(chg)
        # 流れの中の土砂の移流（風上差分）
        Qs = self.QsAt(self.t) if self.phase == 'run' else 0.0
        self.QsIn = Qs; self.sedIn += Qs * dts
        Cin = Qs / self.Qin if (Qs > 0 and self.Qin > 0) else 0.0
        cf = dts / dx
        wl = wet[:, :-1]; wr = wet[:, 1:]
        uf = np.where(wl & wr, 0.5 * (tu[:, :-1] + tu[:, 1:]), np.where(wl, np.maximum(0, tu[:, :-1]), np.where(wr, np.minimum(0, tu[:, 1:]), 0.0)))
        SFx = np.zeros((ny, nx + 1, K)); SFy = np.zeros((ny + 1, nx, K))
        SFx[:, 1:nx] = uf[..., None] * np.where((uf > 0)[..., None], Cc[:, :-1], Cc[:, 1:])
        SFx[:, nx] = np.where(wet[:, -1], np.maximum(0, tu[:, -1]), 0.0)[:, None] * Cc[:, -1]
        wl = wet[:-1]; wr = wet[1:]
        vf = np.where(wl & wr, 0.5 * (tv[:-1] + tv[1:]), np.where(wl, np.maximum(0, tv[:-1]), np.where(wr, np.minimum(0, tv[1:]), 0.0)))
        SFy[1:ny] = vf[..., None] * np.where((vf > 0)[..., None], Cc[:-1], Cc[1:])
        outv = cf * (np.maximum(SFx[:, 1:], 0) - np.minimum(SFx[:, :-1], 0) + np.maximum(SFy[1:], 0) - np.minimum(SFy[:-1], 0))
        scS = np.where(outv > Cc + 1e-15, np.where(outv > 0, Cc / np.where(outv > 0, outv, 1), 0.0), 1.0)
        if (scS < 1).any():
            r = np.ones_like(SFx)
            r[:, 1:] = np.where(SFx[:, 1:] > 0, scS, np.concatenate([scS[:, 1:], np.ones((ny, 1, K))], axis=1))
            SFx *= r
            ry = np.ones_like(SFy); ry[1:ny] = np.where(SFy[1:ny] > 0, scS[:-1], scS[1:]); SFy *= ry
        out_k = SFx[:, nx].sum(axis=0)
        self.sedOutK += out_k * dx * dts
        Cn = Cc - cf * (SFx[:, 1:] - SFx[:, :-1] + SFy[1:] - SFy[:-1])
        if Cin > 0:
            Cn[:, 0, :] += (cf * Cin * np.maximum(self.qin, 0))[:, None] * self.ff[None, :]
        self.Cc = np.maximum(Cn, 0.0).reshape(-1, K)
        self.QsOut = out_k.sum() * dx
        self.sedOut += self.QsOut * dts

    # ------------------------------------------------------------ 2次流（長田・細田・村本・Rahman 1996）
    def _curvature(self):
        """流線の曲率 κ = (u(u v_x + v v_y) − v(u u_x + v u_y)) / V³（左に曲がるとき正）。
        流速は前のステップの面の流束の平均。濡れたセルだけで 3×3 平均を 2 回とり、|h κ| を hrMax で頭打ちにする。"""
        nx, ny, dx, h = self.nx, self.ny, self.dx, self.h
        HS = 0.05; hrMax = (self.P.get('sec') or {}).get('hrMax') or 0.05
        w = h > HS; hs = np.where(w, h, 1.0)
        uf = np.where(w, 0.5 * (self.FxM[:, :-1] + self.FxM[:, 1:]) / hs, 0.0)
        vf = np.where(w, 0.5 * (self.GyM[:-1] + self.GyM[1:]) / hs, 0.0)
        self.uf, self.vf = uf, vf
        tmp = np.zeros((ny, nx))
        ok = w[1:-1, 1:-1] & w[1:-1, :-2] & w[1:-1, 2:] & w[:-2, 1:-1] & w[2:, 1:-1]
        uu = uf[1:-1, 1:-1]; vv = vf[1:-1, 1:-1]; V2 = uu * uu + vv * vv
        ok &= V2 >= 0.01
        ux = (uf[1:-1, 2:] - uf[1:-1, :-2]) / (2 * dx); uy = (uf[2:, 1:-1] - uf[:-2, 1:-1]) / (2 * dx)
        vx = (vf[1:-1, 2:] - vf[1:-1, :-2]) / (2 * dx); vy = (vf[2:, 1:-1] - vf[:-2, 1:-1]) / (2 * dx)
        V2s = np.where(ok, V2, 1.0)
        tmp[1:-1, 1:-1] = np.where(ok, (uu * (uu * vx + vv * vy) - vv * (uu * ux + vv * uy)) / (V2s * np.sqrt(V2s)), 0.0)
        wp = np.pad(w, 1)
        for _ in range(2):
            sp = np.pad(tmp, 1); sm = np.zeros((ny, nx)); nn = np.zeros((ny, nx))
            for dj in (-1, 0, 1):
                for di in (-1, 0, 1):
                    m = wp[1 + dj:1 + dj + ny, 1 + di:1 + di + nx]
                    sm = sm + np.where(m, sp[1 + dj:1 + dj + ny, 1 + di:1 + di + nx], 0.0)
                    nn = nn + m
            tmp = np.where(w, sm / np.maximum(nn, 1), 0.0)
        lim = hrMax / np.maximum(h, 1e-3)
        self.kap = np.clip(tmp, -lim, lim)

    def _sec_stress(self):
        """2次流による運動量の分散項（Kalkwijk・de Vriend 1980、長田ら 1996 式19〜23）。
        Txx = −Tyy = 2 k_sn h² κ u v、Txy = −k_sn h² κ (u² − v²)、k_sn = 5a − 15.6a² + 37.5a³、a = √g/(κ_K C) ≤ 0.4"""
        h = self.h; HS = 0.05
        msk = h > HS; hs = np.where(msk, h, 1.0)
        a = np.minimum(math.sqrt(G) * self.P['manning'] / (0.4 * hs ** (1 / 6)), 0.4)
        ksn = 5 * a - 15.6 * a * a + 37.5 * a * a * a
        wv = np.where(msk, ksn * hs * hs * self.kap, 0.0)
        uu, vv = self.uf, self.vf
        return 2 * wv * uu * vv, -wv * (uu * uu - vv * vv), msk

    def _hasegawa_collapse(self):
        """長谷川の側岸崩落モデル（長谷川 1983、長田ら 1996 式17）。乾いた側岸セルの位置の水面から水中安息角で下る線を引き、
        線より上の土塊 A と線より下の洗掘の不足 B を横断方向（4 方向）に数え、B ≥ A なら A を崩して B を埋め戻す（ΔBc = 1 セル）"""
        BK = self.P.get('bank') or {}
        if BK.get('collapse') != 'hasegawa':
            return
        nx, ny, dx = self.nx, self.ny, self.dx
        HW = 0.01; MAXM = 40
        self.z = np.ascontiguousarray(self.z)
        tb = math.tan(math.radians(BK.get('reposeWet') or self.P['reposeDeg'])) * dx
        z = self.z.ravel(); h = self.h.ravel(); kind = self.kind.ravel(); fxC = self.fixedC.ravel()
        okc = (h >= HW) & (kind == 0) & ~fxC            # 水に浸かった移動床
        dry = (h < HW) & (kind == 0) & ~fxC
        ok2 = okc.reshape(ny, nx); nbw = np.zeros((ny, nx), dtype=bool)
        nbw[:, :-1] |= ok2[:, 1:]; nbw[:, 1:] |= ok2[:, :-1]; nbw[:-1] |= ok2[1:]; nbw[1:] |= ok2[:-1]
        cand = np.flatnonzero(dry & nbw.ravel())
        DI = (1, -1, 0, 0); DJ = (0, 0, 1, -1)
        for c in cand:            # ブラウザ版と同じく 1 セルずつ順に処理する
            i0 = c % nx; j0 = c // nx
            for q in range(4):
                i1 = i0 + DI[q]; j1 = j0 + DJ[q]
                if i1 < 0 or i1 >= nx or j1 < 0 or j1 >= ny:
                    continue
                w1 = j1 * nx + i1
                if not okc[w1]:
                    continue
                H = z[w1] + h[w1]
                if z[c] <= H:
                    continue
                exC = [c]; exV = [z[c] - H]; Aex = z[c] - H; dfC = []; dfV = []; Bdf = 0.0; stage = 0
                for m in range(1, MAXM + 1):
                    ii = i0 + m * DI[q]; jj = j0 + m * DJ[q]
                    if ii < 0 or ii >= nx or jj < 0 or jj >= ny:
                        break
                    cc = jj * nx + ii
                    if not okc[cc]:
                        break
                    dzl = z[cc] - (H - m * tb)
                    if dzl > 0:
                        if stage == 0:
                            exC.append(cc); exV.append(dzl); Aex += dzl
                        else:
                            break
                    elif dzl < 0:
                        stage = 1; dfC.append(cc); dfV.append(-dzl); Bdf -= dzl
                    elif stage == 1:
                        break
                if not dfC or Bdf < Aex:
                    continue
                for e in range(len(exC)):
                    for f in range(len(dfC)):
                        self._move_sed(exC[e], dfC[f], exV[e] * dfV[f] / Bdf)
                self.nHase += 1
                break

    def _move_sed(self, a, b, dz):
        """セル a の上から厚さ dz の土砂を（堆積層の粒度で）セル b へ移す。動かせた厚さを返す"""
        A = self.A; z = self.z.ravel(); K = self.K
        left = dz; it = 0
        while left > 1e-9 and it < 400:
            it += 1
            av = 0.0
            for k in range(K):
                av += A[a, k]
            if av < 1e-9:
                break
            take = left if left < av else av
            frc = take / av
            for k in range(K):
                t = A[a, k] * frc; A[a, k] -= t; A[b, k] += t
            z[a] -= take; z[b] += take; left -= take
            self._fix_layers(np.array([a])); self._fix_layers(np.array([b]))
        return dz - left

    def _dry_bank_erosion(self):
        """水際の側岸侵食：濡れたセルが侵食されたとき、その一部（割合 dryEro）を隣の乾いた側岸セルから削る
        （Delft3D の dry-cell erosion と同じ考え方）"""
        th = (self.P.get('bank') or {}).get('dryEro', 0) or 0
        if th <= 0:
            return
        nx, ny, K = self.nx, self.ny, self.K
        HW = 0.01
        z = self.z.ravel(); h = self.h.ravel(); kind = self.kind.ravel(); fxC = self.fixedC.ravel(); A = self.A
        cand = np.flatnonzero((self.dzs.ravel() < 0) & (h >= HW) & (kind == 0) & ~fxC)
        for c in cand:            # ブラウザ版と同じく 1 セルずつ順に処理する
            i = c % nx
            nb = []
            for n in (c - 1 if i > 0 else -1, c + 1 if i < nx - 1 else -1, c - nx, c + nx):
                if 0 <= n < self.N and h[n] < HW and kind[n] == 0 and not fxC[n] and z[n] > z[c]:
                    nb.append(n)
            if not nb:
                continue
            want = th * (-self.dzs.ravel()[c]) / len(nb)
            for n in nb:
                av = A[n].sum()
                if av < 1e-6:
                    continue
                mv = min(want, 0.5 * av, 0.5 * (z[n] - z[c]))
                if mv < 1e-9:
                    continue
                t = A[n] * (mv / av); A[n] -= t; A[c] += t
                z[n] -= mv; z[c] += mv
                self._fix_layers(np.array([n, c]))
        self.z = z.reshape(ny, nx)

    def sediment_in_water(self):
        """流れの中にある土砂の量（実質体積） [m3]"""
        return float(self.Cc.sum() * self.dx * self.dx)

    def _collapse(self):
        """安息角を超えた斜面の崩落。水に浸かった斜面は水中安息角、乾いた斜面は陸上安息角で判定する。
        ブラウザ版と同じ順序で処理する（x 方向の面を左から順に、y 方向の面を下から順に。同じ列・行の面は互いに独立なのでまとめて計算）"""
        nx, ny, K, dx = self.nx, self.ny, self.K, self.dx
        BK = self.P.get('bank') or {}
        limW = math.tan(math.radians(BK.get('reposeWet') or self.P['reposeDeg'])) * dx
        limD = math.tan(math.radians(BK.get('reposeDry') or self.P['reposeDeg'])) * dx
        HW = 0.01
        hase = BK.get('collapse') == 'hasegawa'
        self.z = np.ascontiguousarray(self.z)
        z = self.z.ravel(); h = self.h.ravel(); A = self.A
        fxC = self.fixedC.ravel(); k0 = self.kind.ravel() == 0
        wet = h > HW
        def faces(ca, cb):
            """面（セル ca と cb の間）の組を一度に処理する。動かしたら True"""
            lim = np.where(wet[ca] | wet[cb], limW, limD)
            diff = z[ca] - z[cb]
            hiA = diff > lim; hiB = -diff > lim
            ok = (hiA | hiB) & ~(fxC[ca] | fxC[cb])
            if hase:   # 長谷川モデルでは、水際の乾いた側岸（高い側が乾き、低い側が濡れている）は崩落モデルに任せる
                kk = k0[ca] & k0[cb]
                ok &= ~(kk & ((hiA & ~wet[ca] & wet[cb]) | (hiB & wet[ca] & ~wet[cb])))
            if not ok.any():
                return False
            ca = ca[ok]; cb = cb[ok]; hiA = hiA[ok]; dd = np.abs(diff[ok]); lim = lim[ok]
            hi = np.where(hiA, ca, cb); lo = np.where(hiA, cb, ca)
            av = np.zeros(hi.size)
            for k in range(K):
                av = av + A[hi, k]
            m = av >= 1e-6
            mv = np.minimum(0.25 * (dd - lim), 0.5 * av)
            m &= mv >= 1e-6
            if not m.any():
                return False
            hi = hi[m]; lo = lo[m]; mv = mv[m]; frc = mv / av[m]
            t = A[hi] * frc[:, None]
            A[hi] -= t; A[lo] += t
            z[hi] -= mv; z[lo] += mv
            self._fix_layers(hi); self._fix_layers(lo)
            return True
        jx = np.arange(ny) * nx
        for i in range(nx - 1):          # x 方向の面
            faces(jx + i, jx + i + 1)
        ix = np.arange(nx)
        for j in range(ny - 1):          # y 方向の面
            faces(j * nx + ix, (j + 1) * nx + ix)

    def _fix_layers(self, cells):
        """交換層厚を La に保つよう、堆積層との間で土砂を出し入れする（cells：平坦化した番号の配列）"""
        A = self.A; SF = self.SF; nl = self.nl; th = self.topTh
        La, Ls, NL = self.La, self.Ls, self.NL
        As = A[cells].sum(axis=1)
        # 堆積：余剰分を交換層の粒度のまま堆積層へ
        dep = cells[As > La + 1e-12]
        if dep.size:
            Ad = A[dep]; s = Ad.sum(axis=1); e = s - La; fr = e / s
            n0 = nl[dep].copy(); t0 = th[dep].copy()
            new = n0 == 0
            n0[new] = 1; t0[new] = 0.0
            tk = Ad * fr[:, None]
            A[dep] = Ad - tk
            tn = t0 + e
            lay = n0 - 1
            SF[dep, lay] = ((SF[dep, lay] * t0[:, None] + tk) / tn[:, None]).astype(np.float32)
            t0 = tn
            grow = (t0 > Ls) & (n0 < NL)
            if grow.any():
                g_ = dep[grow]
                SF[g_, n0[grow]] = SF[g_, n0[grow] - 1]
                n0[grow] += 1; t0[grow] -= Ls
            nl[dep] = n0; th[dep] = t0
        # 侵食：不足分を堆積層最上層の粒度で補う
        ero = cells[(As < La - 1e-12) & (nl[cells] > 0)]
        if ero.size:
            deficit = La - A[ero].sum(axis=1)
            for _ in range(NL + 1):
                m = (deficit > 1e-12) & (nl[ero] > 0)
                if not m.any():
                    break
                c = ero[m]
                take = np.minimum(deficit[m], th[c])
                A[c] += take[:, None] * SF[c, nl[c] - 1]
                th[c] -= take; deficit[m] -= take
                empty = th[c] <= 1e-9
                ce = c[empty]
                nl[ce] -= 1
                th[ce] = np.where(nl[ce] > 0, Ls, 0.0)

    # ------------------------------------------------------------ 集計
    def water_volume(self):
        return float(self.h.sum() * self.dx * self.dx)

    def storage(self):
        """計算区間内の堆積土砂量の増分（空隙込み） [m3]"""
        return float((self.z - self.z0).sum() * self.dx * self.dx)

    def storage_by_dam(self):
        """堰堤ごとの堆砂量（各堰堤の上流側区間。最後の要素は最下流堰堤より下流） [m3]"""
        nd = len(self.dams)
        return np.bincount(self.region.ravel(), weights=(self.z - self.z0).ravel(), minlength=nd + 1)[:nd + 1] * self.dx ** 2

    def surface_dm(self):
        """表層（交換層）の平均粒径 [m]（(ny, nx) の配列）"""
        A = self.A.reshape(self.ny, self.nx, self.K); s = A.sum(axis=2)
        return np.where(s > 1e-9, (A * self.d).sum(axis=2) / np.where(s > 1e-9, s, 1), np.nan)

    def profile(self):
        """縦断分布（横断方向に集約）"""
        fl = self.fl & (self.kind != 2)
        z, z0, h = self.z, self.z0, self.h
        x = (np.arange(self.nx) + 0.5) * self.dx
        with np.errstate(invalid='ignore', divide='ignore'):
            cnt = fl.sum(axis=0)
            zmin = np.where(cnt > 0, np.where(fl, z, np.inf).min(axis=0), np.nan)
            z0min = np.where(cnt > 0, np.where(fl, z0, np.inf).min(axis=0), np.nan)
            zmean = np.where(cnt > 0, np.where(fl, z, 0).sum(axis=0) / cnt, np.nan)
            z0mean = np.where(cnt > 0, np.where(fl, z0, 0).sum(axis=0) / cnt, np.nan)
            w = fl & (h > 0.03)
            hw = np.where(w, h, 0).sum(axis=0)
            ws = np.where(hw > 0, np.where(w, (z + h) * h, 0).sum(axis=0) / hw, np.nan)
            wet = w.sum(axis=0) * self.dx
            dms = self.surface_dm()
            dm = np.where(w.sum(axis=0) > 0, np.nansum(np.where(w, dms, 0), axis=0) / np.maximum(1, (w & np.isfinite(dms)).sum(axis=0)), np.nan)
        return dict(x=x, zmin=zmin, zmean=zmean, z0min=z0min, z0mean=z0mean, ws=ws, wet=wet, dm=dm)

    def dams_fixed(self):
        """いまの堰堤を、位置・回転角・中心 y・天端標高・堤長・水通しの中心まで固定した諸元にする（続きの計算で同じ場所に置くため）"""
        out = []
        for D in self.dams:
            out.append({'x': round(D['x'], 3), 'H': D['H'], 'b': round(max(0.0, D['spMax'] - D['spMin']), 3),
                        'ang': round(math.degrees(D['th']), 4), 'y': round(float(D['cy']), 4), 'zc': round(float(D['crest']), 3),
                        'len': round(2 * max(abs(D['sMin']), abs(D['sMax'])), 2), 'so': round((D['spMin'] + D['spMax']) / 2, 4)})
        return out

    def continue_params(self):
        """いまの河床を初期地形にし、堰堤を同じ位置・向き・天端標高で置いた条件を返す（粒度は初期の河床材料に戻る）"""
        P = copy.deepcopy(self.P)
        zz = self.z.copy()
        dam = self.kind == 2
        if dam.any():
            pad = np.pad(np.where(dam, np.nan, zz), 1, constant_values=np.nan)
            nb = np.stack([pad[1:-1, :-2], pad[1:-1, 2:], pad[:-2, 1:-1], pad[2:, 1:-1]])
            with np.errstate(invalid='ignore'):
                mean = np.nanmean(nb, axis=0)
            zz = np.where(dam, np.minimum(np.where(np.isfinite(mean), mean, zz), self.z - 0.01), zz)   # 堰堤の下は天端より少し低く
        fl = self.kind == 0
        with np.errstate(invalid='ignore'):
            r = np.where(fl, self.z, np.nan)
            fh = float(np.nanmax(np.nanmax(r, axis=0) - np.nanmin(r, axis=0)))
        mask = ((self.kind == 0) | ((self.kind == 2) & self.fl)).astype(int)
        P['terrain'] = {'rows': zz[::-1].tolist(), 'cs': self.dx, 'mask': mask[::-1].tolist()}
        P['floorH'] = min(20.0, max(0.5, math.ceil((fh + 0.05) * 10) / 10)); P['csvSlope'] = 0.0; P['prefill'] = 0.0
        P['dams'] = self.dams_fixed()
        return P

    def bed_csv(self, path):
        """いまの河床を格子 CSV（行列形式。1行目が左岸側）に書き出す。ブラウザ版でも読み込める"""
        zz = self.z.copy()
        dam = self.kind == 2
        if dam.any():       # 堰堤のセルは、周りの地盤の高さにする
            pad = np.pad(np.where(dam, np.nan, zz), 1, constant_values=np.nan)
            nb = np.stack([pad[1:-1, :-2], pad[1:-1, 2:], pad[:-2, 1:-1], pad[2:, 1:-1]])
            with np.errstate(invalid='ignore'):
                mean = np.nanmean(nb, axis=0)
            zz = np.where(dam, np.minimum(np.where(np.isfinite(mean), mean, zz), self.z - 0.01), zz)   # 堰堤の下は天端より少し低く
        fl = self.kind == 0
        with np.errstate(invalid='ignore'):
            rng_ = np.where(fl, self.z, np.nan)
            fh = float(np.nanmax(np.nanmax(rng_, axis=0) - np.nanmin(rng_, axis=0)))
        fh = min(20.0, max(0.5, math.ceil((fh + 0.05) * 10) / 10))
        with open(path, 'w', encoding='utf-8-sig', newline='\r\n') as f:
            f.write('# 計算後の河床（出水の時刻 %.1f 分）　1行目が左岸側、左の列が上流　格子間隔 %g m\n' % (self.t / 60, self.dx))
            f.write('# floorH=%g\n' % fh)
            for d in self.dams_fixed():
                f.write('# dam ' + ' '.join('%s=%s' % (k, repr(float(d[k]))) for k in ('x', 'H', 'b', 'ang', 'y', 'zc', 'len', 'so')) + '\n')
            mask = (self.kind == 0) | ((self.kind == 2) & self.fl)       # 谷底の印（1 が谷底）
            for j in range(self.ny - 1, -1, -1):
                f.write('# floor ' + ''.join('1' if v else '0' for v in mask[j]) + '\n')
            for j in range(self.ny - 1, -1, -1):
                f.write(','.join('%.3f' % v for v in zz[j]) + '\n')


# =============================================================================
# 計算を進める
# =============================================================================
def run(sim, until=None, max_steps=10 ** 8, report=True, every=30.0, callback=None):
    """出水が終わるまで（または sim.t >= until [s] まで）計算する。
    callback(sim) を与えると、約 every 秒（計算機の時間）ごとに呼ぶ（途中の図を描くときなど）。"""
    t0 = time.time(); last = t0; n = 0
    while sim.phase != 'done' and n < max_steps and (until is None or sim.t < until):
        sim.step(); n += 1
        now = time.time()
        if now - last > every:
            last = now
            if report:
                print('  %s  出水の時刻 %6.0f / %.0f s  ステップ %d  Δt %.3f s  流出 %.1f m3/s  （%.0f 秒経過）'
                      % ('助走' if sim.phase == 'spinup' else '計算', sim.t, sim.duration(), sim.nstep, sim.dt, sim.Qout, now - t0))
            if callback:
                callback(sim)
    if report:
        sb = sim.storage_by_dam()
        print('終了：ステップ %d、%.0f 秒。供給土砂 %.0f m3、流出土砂 %.0f m3（実質体積）、区間内の堆積 %+.0f m3（空隙込み）'
              % (sim.nstep, time.time() - t0, sim.sedIn, sim.sedOut, sim.storage()))
        print('  土砂の体積収支の誤差 %.1e m3、水の体積収支の誤差 %.1e m3'
              % ((sim.sedIn - sim.sedOut) / (1 - sim.P['porosity']) - sim.storage(),
                 sim.Vw0 + sim.Vwin - sim.Vwout - sim.water_volume()))
    return sim


def dam_table(sim):
    """堰堤ごとの堆砂量の表（行のリスト）"""
    sb = sim.storage_by_dam(); rows = []
    for q, D in enumerate(sim.dams):
        rows.append(('堰堤 %d の上流' % (q + 1), D['x'], D['H'], D['b'], sb[q]))
    rows.append(('最下流の堰堤より下流', float('nan'), float('nan'), float('nan'), sb[len(sim.dams)]))
    return rows


def print_dam_table(sim):
    print('%-16s %8s %6s %8s %12s' % ('区間', 'x [m]', '堤高', '水通し', '堆砂量 [m3]'))
    for r in dam_table(sim):
        if r[1] == r[1]:
            print('%-16s %8.0f %6.1f %8.0f %12.0f' % r)
        else:
            print('%-16s %8s %6s %8s %12.0f' % (r[0], '', '', '', r[4]))


# =============================================================================
# 地形ファイルの読み込み
# =============================================================================
def read_lines(path):
    """UTF-8（BOM 付きも可）で読み、だめなら Shift_JIS（cp932）で読む"""
    for enc in ('utf-8-sig', 'cp932'):
        try:
            with open(path, encoding=enc) as f:
                return f.read().splitlines()
        except UnicodeDecodeError:
            continue
    with open(path, encoding='utf-8', errors='replace') as f:
        return f.read().splitlines()


def _numeric_rows(lines):
    rows = []
    for ln in lines:
        ln = ln.strip()
        if not ln or ln.startswith('#'):
            continue
        tok = [t for t in ln.replace('\t', ',').replace(';', ',').replace(' ', ',').split(',') if t != '']
        try:
            rows.append([float(t) for t in tok])
        except ValueError:
            continue      # 見出し行
    return rows


def read_grid(path, cs=2.0, transpose=False, flip=False):
    """格子 CSV を読む。行列形式（1行目が左岸側、左の列が上流）か x, y, z 形式（y は左岸側が大きい）。
    戻り値 {'rows': 2次元配列, 'cs': 格子間隔}"""
    lines = read_lines(path)
    data = _numeric_rows(lines)
    if len(data) < 3:
        raise ValueError('数値の行が 3 行以上必要です。')
    dams, floorH, mask = [], None, []
    import re as _re
    for ln in lines:      # bed_csv・ブラウザ版が書いた堰堤の諸元（# dam x=… H=… …）
        m = _re.match(r'\s*#\s*dam\s+(.*)$', ln, _re.I)
        if m:
            d = {k: float(v) for k, v in _re.findall(r'(\w+)=(-?[\d.eE+-]+)', m.group(1))}
            if 'x' in d and d.get('H', 0) > 0:
                dams.append(d)
        f = _re.match(r'\s*#\s*floorH=(-?[\d.]+)', ln, _re.I)
        if f:
            floorH = float(f.group(1))
        fm = _re.match(r'\s*#\s*floor\s+([01]+)\s*$', ln)
        if fm:
            mask.append([int(ch) for ch in fm.group(1)])
    rows = None
    if len(data) >= 9 and all(len(r) == 3 for r in data):
        arr = np.array(data)
        ux = np.unique(np.round(arr[:, 0], 6)); uy = np.unique(np.round(arr[:, 1], 6))
        if len(ux) >= 3 and len(uy) >= 3 and len(ux) * len(uy) <= len(data) * 1.6 + 4:
            csx = float(np.median(np.diff(ux))); csy = float(np.median(np.diff(uy)))
            nc = _jround((ux[-1] - ux[0]) / csx) + 1; nr = _jround((uy[-1] - uy[0]) / csy) + 1
            g = np.full((nr, nc), np.nan)
            ii = np.round((arr[:, 0] - ux[0]) / csx).astype(int); jj = nr - 1 - np.round((arr[:, 1] - uy[0]) / csy).astype(int)
            ok = (ii >= 0) & (ii < nc) & (jj >= 0) & (jj < nr)
            g[jj[ok], ii[ok]] = arr[ok, 2]
            for _ in range(200):      # 抜けを周りの平均で埋める
                miss = np.isnan(g)
                if not miss.any():
                    break
                pad = np.pad(g, 1, constant_values=np.nan)
                nb = np.stack([pad[1:-1, :-2], pad[1:-1, 2:], pad[:-2, 1:-1], pad[2:, 1:-1]])
                with np.errstate(invalid='ignore'):
                    m = np.nanmean(nb, axis=0)
                g = np.where(miss & np.isfinite(m), m, g)
            if abs(csy - csx) > 0.05 * csx:     # y 方向の間隔を x にそろえる
                nr2 = max(3, _jround((nr - 1) * csy / csx) + 1)
                f = np.arange(nr2) * (nr - 1) / (nr2 - 1); a = np.floor(f).astype(int); t = f - a; b = np.minimum(nr - 1, a + 1)
                g = g[a] * (1 - t[:, None]) + g[b] * t[:, None]
            rows = g; cs = csx
    if rows is None:
        L = max(set(len(r) for r in data), key=lambda n: sum(1 for r in data if len(r) == n))
        rows = np.array([r[:L] for r in data if len(r) >= L])
    rows = np.asarray(rows, dtype=float)
    if transpose:
        rows = rows.T
    if flip:
        rows = rows[:, ::-1]
    if rows.shape[1] < 10 or rows.shape[0] < 3:
        raise ValueError('流下方向に 10 点以上、横断方向に 3 点以上が必要です（いまは %d × %d）。横断・縦断の CSV なら read_survey を使ってください。' % (rows.shape[1], rows.shape[0]))
    mask = np.array(mask) if mask and not transpose and not flip and len(set(map(len, mask))) == 1 else None
    if mask is not None and mask.shape != rows.shape:
        mask = None
    return {'rows': rows, 'cs': float(cs), 'dams': dams, 'floorH': floorH, 'mask': mask}


def survey_to_grid(xs_rows, long_rows=None, dx=2.0, align='thalweg', wall=1.0):
    """横断（x, s, z）と縦断（x, z）の測量値から格子地形を作る（xs2grid.py・ブラウザ版と同じ方法）。
    断面の間は、最深部の左右それぞれを幅で正規化して比例配分で補間する。"""
    by = {}
    for r in xs_rows:
        by.setdefault(round(r[0], 6), []).append((r[1], r[2]))
    secs = []
    for x in sorted(by):
        pts = sorted(by[x])
        if len(pts) < 3:
            raise ValueError('x = %g m の断面の点が 3 点未満です。' % x)
        s = np.array([p[0] for p in pts]); z = np.array([p[1] for p in pts]); zmin = z.min()
        st = s[z <= zmin + 1e-9].mean()
        secs.append(dict(x=x, s=s, zr=z - zmin, zmin=zmin, st=st, sl=s[0], sr=s[-1]))
    if len(secs) < 2:
        raise ValueError('断面が 2 つ以上必要です。')
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
        c['uL'] = c['sl'] + c['off']; c['uT'] = c['st'] + c['off']; c['uR'] = c['sr'] + c['off']

    def secz(c, side, e):
        s = c['sl'] + (c['st'] - c['sl']) * e if side < 0 else c['st'] + (c['sr'] - c['st']) * e
        return np.interp(s, c['s'], c['zr'])
    longp = sorted(long_rows) if long_rows else None
    xs_ = [c['x'] for c in secs]
    allx = xs_ + ([p[0] for p in longp] if longp else [])
    xa, xb = min(allx), max(allx)
    L = xb - xa; W = max(c['uR'] for c in secs)
    nx = max(10, _jround(L / dx)); ny = max(3, _jround(W / dx) + 1)
    grid = np.zeros((ny, nx))
    u = (np.arange(ny) + 0.5) * dx
    for i in range(nx):
        x = xa + (i + 0.5) * dx
        if x <= xs_[0]:
            a = b = secs[0]; t = 0.0
        elif x >= xs_[-1]:
            a = b = secs[-1]; t = 0.0
        else:
            k = 1
            while xs_[k] < x:
                k += 1
            a, b = secs[k - 1], secs[k]; t = (x - a['x']) / (b['x'] - a['x'])
        uL = a['uL'] * (1 - t) + b['uL'] * t; uT = a['uT'] * (1 - t) + b['uT'] * t; uR = a['uR'] * (1 - t) + b['uR'] * t
        zt = np.interp(x, [p[0] for p in longp], [p[1] for p in longp]) if longp else a['zmin'] * (1 - t) + b['zmin'] * t
        zL = secz(a, -1, 0.0) * (1 - t) + secz(b, -1, 0.0) * t; zR = secz(a, 1, 1.0) * (1 - t) + secz(b, 1, 1.0) * t
        eL = np.where(uT > uL, (u - uL) / max(uT - uL, 1e-12), 1.0); eR = np.where(uR > uT, (u - uT) / max(uR - uT, 1e-12), 0.0)
        zl = secz(a, -1, np.clip(eL, 0, 1)) * (1 - t) + secz(b, -1, np.clip(eL, 0, 1)) * t
        zr_ = secz(a, 1, np.clip(eR, 0, 1)) * (1 - t) + secz(b, 1, np.clip(eR, 0, 1)) * t
        zrel = np.where(u < uL, zL + wall * (uL - u), np.where(u > uR, zR + wall * (u - uR), np.where(u <= uT, zl, zr_)))
        grid[:, i] = zt + zrel
    return {'rows': grid, 'cs': float(dx), 'nsec': len(secs), 'L': L, 'W': W}


def read_survey(xs_path, long_path=None, dx=2.0, align='thalweg', wall=1.0):
    """横断 CSV（x, s, z）と縦断 CSV（x, z）を読んで格子地形にする"""
    xs = [r[:3] for r in _numeric_rows(read_lines(xs_path)) if len(r) >= 3]
    lo = [r[:2] for r in _numeric_rows(read_lines(long_path)) if len(r) >= 2] if long_path else None
    if not xs:
        raise ValueError('%s に x, s, z の行がありません。' % xs_path)
    return survey_to_grid(xs, lo, dx, align, wall)


def classify_csv(path):
    """CSV の中身から種類を見分ける：'xs'（横断 x,s,z）／'long'（縦断 x,z）／'grid'（格子）"""
    data = _numeric_rows(read_lines(path))
    if len(data) < 2:
        return 'none'
    lens = [len(r) for r in data]; m = max(set(lens), key=lens.count)
    if m == 2:
        return 'long'
    if m >= 4:
        return 'grid'
    try:
        arr = np.array([r for r in data if len(r) == 3])
        ux = np.unique(np.round(arr[:, 0], 6)); uy = np.unique(np.round(arr[:, 1], 6))
        if len(arr) >= 9 and len(ux) >= 3 and len(uy) >= 3 and len(ux) * len(uy) <= len(arr) * 1.6 + 4:
            return 'grid'
    except Exception:
        pass
    return 'xs'


def use_terrain(P, T, floorH=None, csvSlope=0.0):
    """読み込んだ地形 T（read_grid / read_survey の戻り値）を条件 P に入れる。満砂の設定は外す。
    T に堰堤の諸元（bed_csv で書いた CSV）があれば、その堰堤も入れる"""
    P['terrain'] = {'rows': np.asarray(T['rows'], dtype=float).tolist(), 'cs': T['cs']}
    if T.get('mask') is not None:
        P['terrain']['mask'] = np.asarray(T['mask']).tolist()
    P['floorH'] = floorH if floorH is not None else (T.get('floorH') or 3.0)
    P['csvSlope'] = csvSlope; P['prefill'] = 0.0
    if T.get('dams'):
        P['dams'] = [dict(d) for d in T['dams']]
    return P


# =============================================================================
# 図（Matplotlib）
# =============================================================================
def setup_japanese_font():
    """日本語の文字が出るフォントを探して使う（見つからなければ何もしない）"""
    import matplotlib
    from matplotlib import font_manager
    names = {f.name for f in font_manager.fontManager.ttflist}
    for cand in ('Noto Sans CJK JP', 'Noto Sans JP', 'IPAexGothic', 'IPAGothic', 'Yu Gothic', 'Meiryo',
                 'Hiragino Sans', 'Hiragino Kaku Gothic ProN', 'MS Gothic', 'TakaoGothic', 'BIZ UDGothic'):
        if cand in names:
            matplotlib.rcParams['font.family'] = cand
            matplotlib.rcParams['axes.unicode_minus'] = False
            return cand
    print('注意：日本語のフォントが見つかりません。図の文字が □ になるときは、IPAexGothic などを入れてください。')
    return None


def _hillshade(z, dx, vex=2.0):
    gy, gx = np.gradient(z * vex, dx)
    nxv, nyv, nzv = -gx, -gy, np.ones_like(z)
    nrm = np.sqrt(nxv ** 2 + nyv ** 2 + nzv ** 2)
    L = np.array([-0.45, 0.4, 0.8]); L = L / np.linalg.norm(L)
    return np.clip((nxv * L[0] + nyv * L[1] + nzv * L[2]) / nrm, 0, 1)


def ground_rgb(sim):
    """地表の色（写実）：森林・裸地・礫（表層の粒径で濃淡）・堰堤。(ny, nx, 3)、0〜1"""
    rgb = np.zeros((sim.ny, sim.nx, 3))
    dm = sim.surface_dm()
    f = np.clip(np.log(np.nan_to_num(dm, nan=0.02) / 0.004) / math.log(0.15 / 0.004), 0, 1)
    rgb[..., 0] = (208 - 72 * f) / 255; rgb[..., 1] = (194 - 64 * f) / 255; rgb[..., 2] = (166 - 40 * f) / 255
    rnd = (_hashI(np.arange(sim.N).reshape(sim.ny, sim.nx), 7, 3) + 1) / 2
    wall = (sim.kind == 1) & (sim.z - sim.z0 < 0.02)
    forest = wall & ~sim.bare
    rgb[forest] = np.stack([(46 + 30 * rnd) / 255, (78 + 22 * rnd) / 255, (42 + 12 * rnd) / 255], axis=-1)[forest]
    rgb[wall & sim.bare] = np.array([184, 176, 160]) / 255
    bankV = sim.bankC & (np.abs(sim.z - sim.z0) < 0.05)          # まだ削れていない側岸（段丘の森林）
    rgb[bankV] = np.stack([(52 + 30 * rnd) / 255, (84 + 22 * rnd) / 255, (44 + 12 * rnd) / 255], axis=-1)[bankV]
    rgb[sim.fixedC] = np.array([128, 126, 120]) / 255              # 動かない巨礫
    rgb[sim.kind == 2] = np.array([190, 188, 180]) / 255
    return rgb


def plan_image(sim, vex=2.0):
    """平面図の画像（陰影＋水）。(ny, nx, 3)"""
    rgb = ground_rgb(sim) * (0.55 + 0.6 * _hillshade(sim.z, sim.dx, vex))[..., None]
    h = sim.h
    a = np.clip(1 - np.exp(-h * 3.2), 0, 0.88)[..., None] * (h > 0.01)[..., None]
    t = (1 - np.exp(-h * 1.4))[..., None]
    wc = np.array([96, 172, 170]) / 255 * (1 - t) + np.array([22, 78, 104]) / 255 * t
    return np.clip(rgb * (1 - a) + wc * a, 0, 1)


def plot_plan(sim, ax=None, var='real', cores=True, title=None):
    """平面図。var: 'real'（写実）, 'dz'（河床変動量）, 'dm'（表層の平均粒径）, 'depth'（水深）, 'bed'（標高）"""
    import matplotlib.pyplot as plt
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 12 * sim.Ly / sim.Lx + 1))
    ext = [0, sim.Lx, 0, sim.Ly]
    if var == 'real':
        ax.imshow(plan_image(sim), origin='lower', extent=ext, interpolation='bilinear')
    else:
        if var == 'dz':
            data = np.where((sim.kind == 0) | (np.abs(sim.z - sim.z0) > 1e-3), sim.z - sim.z0, np.nan)
            m = max(0.5, np.nanpercentile(np.abs(data), 99)); kw = dict(cmap='RdBu_r', vmin=-m, vmax=m); lab = '河床変動量 [m]（赤：堆積）'
        elif var == 'dm':
            data = sim.surface_dm() * 1000; kw = dict(cmap='YlOrBr'); lab = '表層の平均粒径 [mm]'
        elif var == 'depth':
            data = np.where(sim.h > 0.01, sim.h, np.nan); kw = dict(cmap='Blues'); lab = '水深 [m]'
        else:
            data = sim.z; kw = dict(cmap='terrain'); lab = '標高 [m]'
        ax.imshow(_hillshade(sim.z, sim.dx), origin='lower', extent=ext, cmap='gray', vmin=-0.6, vmax=1.1)
        im = ax.imshow(data, origin='lower', extent=ext, alpha=0.85, interpolation='nearest', **kw)
        plt.colorbar(im, ax=ax, label=lab, fraction=0.025, pad=0.01)
    for D in sim.dams:
        for s0, s1, sp in ((D['sMin'], D['spMin'], False), (D['spMin'], D['spMax'], True), (D['spMax'], D['sMax'], False)):
            if s1 - s0 > 1e-6:
                ax.plot([D['cx'] + D['nx'] * s0, D['cx'] + D['nx'] * s1], [D['cy'] + D['ny'] * s0, D['cy'] + D['ny'] * s1],
                        color='white' if sp else '#222', lw=3 if sp else 2.5, solid_capstyle='butt')
    if cores and var == 'real':
        for xs_, ys_, w in flow_cores(sim):
            ax.plot(xs_, ys_, color='white', lw=w, alpha=0.9)
    ax.set_xlim(0, sim.Lx); ax.set_ylim(0, sim.Ly); ax.set_aspect('equal')
    ax.set_xlabel('流下方向の距離 x [m]（流れは左→右）'); ax.set_ylabel('y [m]（上が左岸）')
    ax.set_title(title or ('平面図　出水の時刻 %.0f 分' % (sim.t / 60)))
    return ax


def flow_cores(sim, min_len=6):
    """流心線：各横断で単位幅流量の極大を求め、隣の横断の極大どうしを結ぶ。[(x, y, 線の太さ), ...]"""
    q = np.where((sim.kind != 2) & (sim.h > 0.03), np.hypot(sim.qx, sim.qy), 0.0)
    if q.max() <= 0:
        return []
    k = np.array([1, 2, 3, 2, 1], float); k /= k.sum()
    qs = np.apply_along_axis(lambda c: np.convolve(c, k, mode='same'), 0, q)
    qs = np.apply_along_axis(lambda r: np.convolve(r, np.array([1, 2, 1]) / 4, mode='same'), 1, qs)
    thr = 0.12 * qs.max(); dx = sim.dx
    tracks = []; open_ = []
    for i in range(sim.nx):
        c = qs[:, i]
        pk = [j for j in range(1, sim.ny - 1) if c[j] >= thr and c[j] >= c[j - 1] and c[j] > c[j + 1]]
        new_open = []
        for j in pk:
            best = None
            for tr in open_:
                if abs(tr['j'][-1] - j) <= 3 and (best is None or abs(tr['j'][-1] - j) < abs(best['j'][-1] - j)):
                    best = tr
            if best is not None and best['i'][-1] == i - 1 and best not in new_open:
                best['i'].append(i); best['j'].append(j); best['q'].append(c[j]); new_open.append(best)
            else:
                tr = {'i': [i], 'j': [j], 'q': [c[j]]}
                if best is not None:          # 分岐：前の点から始める
                    tr['i'].insert(0, best['i'][-1]); tr['j'].insert(0, best['j'][-1]); tr['q'].insert(0, best['q'][-1])
                tracks.append(tr); new_open.append(tr)
        open_ = new_open
    out = []
    qmax = qs.max()
    for tr in tracks:
        if len(tr['i']) < min_len:
            continue
        out.append(((np.array(tr['i']) + 0.5) * dx, (np.array(tr['j']) + 0.5) * dx, 0.6 + 2.2 * float(np.mean(tr['q'])) / qmax))
    return out


def plot_profile(sim, ax=None, use_min=False, detrend=False):
    """縦断図：初期河床・いまの河床・水面・堰堤"""
    import matplotlib.pyplot as plt
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 4))
    p = sim.profile(); x = p['x']
    tr = (sim.trA + sim.trB * x) if detrend else 0 * x
    z0 = p['z0min'] if use_min else p['z0mean']; z1 = p['zmin'] if use_min else p['zmean']
    ax.plot(x, z0 - tr, color='#888', lw=1.2, label='初期河床')
    ax.fill_between(x, np.nanmin(z1 - tr) - 2, z1 - tr, color='#c9b38f', alpha=0.6)
    ax.plot(x, z1 - tr, color='#7a5a32', lw=1.6, label='いまの河床（%s）' % ('最深' if use_min else '横断平均'))
    ax.plot(x, p['ws'] - tr, color='#2a7fb8', lw=1.2, label='水面')
    for D in sim.dams:
        xd = D['cx']; zb = D['zlow'] - (sim.trA + sim.trB * xd if detrend else 0)
        ax.plot([xd, xd], [zb, D['crest'] - (sim.trA + sim.trB * xd if detrend else 0)], color='#333', lw=3)
    ax.set_xlabel('x [m]'); ax.set_ylabel('標高 [m]' + ('（初期勾配を差し引いた値）' if detrend else ''))
    ax.set_title('縦断図'); ax.legend(loc='best', fontsize=9); ax.grid(alpha=0.3)
    return ax


def plot_xsec(sim, x, ax=None):
    """横断図（下流を向いて見る。左が左岸）"""
    import matplotlib.pyplot as plt
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 3.5))
    i = int(np.clip(round(x / sim.dx - 0.5), 0, sim.nx - 1))
    s = sim.Ly - (np.arange(sim.ny) + 0.5) * sim.dx          # 左岸からの距離
    zz = sim.z[:, i]; z0 = sim.z0[:, i]; ws = np.where(sim.h[:, i] > 0.01, zz + sim.h[:, i], np.nan)
    ax.fill_between(s, zz.min() - 2, zz, color='#d8c7a8')
    ax.plot(s, z0, color='#888', lw=1, label='初期河床'); ax.plot(s, zz, color='#7a5a32', lw=1.6, label='河床')
    ax.plot(s, ws, color='#2a7fb8', lw=1.4, label='水面')
    ax.set_xlabel('左岸からの距離 [m]'); ax.set_ylabel('標高 [m]'); ax.set_title('横断図（x = %.0f m）' % ((i + 0.5) * sim.dx))
    ax.legend(fontsize=9); ax.grid(alpha=0.3)
    return ax


def plot_3d(sim, ax=None, vex=2.0, stride=None, elev=30, azim=-72):
    """鳥瞰図（Matplotlib の 3D 面。重いので格子を間引く）"""
    import matplotlib.pyplot as plt
    if ax is None:
        fig = plt.figure(figsize=(13, 6)); ax = fig.add_subplot(111, projection='3d')
        fig.subplots_adjust(left=0, right=1, bottom=0, top=0.95)
    st = stride or max(1, int(max(sim.nx, sim.ny) / 220))
    X = (np.arange(sim.nx) + 0.5) * sim.dx; Y = (np.arange(sim.ny) + 0.5) * sim.dx
    XX, YY = np.meshgrid(X, Y)
    zw = np.where(sim.h > 0.01, sim.z + sim.h, sim.z)
    col = plan_image(sim, vex)
    ax.plot_surface(XX[::st, ::st], YY[::st, ::st], (zw[::st, ::st] - sim.zmid) * vex, facecolors=col[::st, ::st],
                    rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)
    try:
        ax.set_box_aspect((sim.Lx, sim.Ly, max(1.0, (sim.zhi - sim.zlo) * vex)), zoom=1.0)
    except TypeError:          # 古い Matplotlib
        ax.set_box_aspect((sim.Lx, sim.Ly, max(1.0, (sim.zhi - sim.zlo) * vex)))
    ax.view_init(elev=elev, azim=azim); ax.set_axis_off()
    ax.set_title('鳥瞰図（高さ %g 倍）' % vex)
    return ax


def plot_series(sim, ax=None):
    """流量と流砂量の時間変化"""
    import matplotlib.pyplot as plt
    if ax is None:
        fig, ax = plt.subplots(1, 2, figsize=(12, 3.5))
    S = sim.series; t = np.array(S['t']) / 60
    ax[0].plot(t, S['Qin'], label='流入', color='#2a7fb8'); ax[0].plot(t, S['Qout'], label='流出', color='#d0702a')
    ax[0].set_xlabel('出水の時刻 [分]'); ax[0].set_ylabel('流量 [m³/s]'); ax[0].legend(); ax[0].grid(alpha=0.3)
    ax[1].plot(t, S['QsIn'], label='供給', color='#2a7fb8'); ax[1].plot(t, S['QsOut'], label='流出', color='#d0702a')
    ax[1].set_xlabel('出水の時刻 [分]'); ax[1].set_ylabel('流砂量 [m³/s]（実質体積）'); ax[1].legend(); ax[1].grid(alpha=0.3)
    return ax


def plot_all(sim, path=None, x_section=None):
    """平面（写実・河床変動量）、縦断、横断、時間変化をまとめて描く。path を与えると保存する"""
    import matplotlib.pyplot as plt
    setup_japanese_font()
    asp = sim.Ly / sim.Lx
    fig = plt.figure(figsize=(14, 14 * asp * 2 + 9))
    gs = fig.add_gridspec(5, 2, height_ratios=[asp * 14, asp * 14, 4, 3.5, 3.5])
    plot_plan(sim, fig.add_subplot(gs[0, :]))
    plot_plan(sim, fig.add_subplot(gs[1, :]), var='dz', title='河床変動量')
    plot_profile(sim, fig.add_subplot(gs[2, :]))
    plot_xsec(sim, x_section if x_section is not None else sim.Lx * 0.55, fig.add_subplot(gs[3, 0]))
    ax = fig.add_subplot(gs[3, 1]); ax.axis('off')
    rows = dam_table(sim)
    txt = '\n'.join('%-14s  堆砂量 %+8.0f m³' % (r[0], r[4]) for r in rows)
    ax.text(0, 1, '堰堤ごとの堆砂量（空隙込み）\n' + txt, va='top', family=plt.rcParams['font.family'], fontsize=9)
    plot_series(sim, [fig.add_subplot(gs[4, 0]), fig.add_subplot(gs[4, 1])])
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=110)
    return fig


# =============================================================================
# コマンドライン
# =============================================================================
def _parse_dams(s):
    out = []
    for tok in s.split(','):
        tok = tok.strip()
        if not tok:
            continue
        a = [float(v) for v in tok.split(':')]
        out.append({'x': a[0], 'H': a[1], 'b': a[2] if len(a) > 2 else 0.0})
    return out


def main(argv=None):
    import argparse
    import os
    ap = argparse.ArgumentParser(description='連続砂防堰堤 平面2次元河床変動計算（Python 版）')
    ap.add_argument('--preset', choices=['A', 'B', 'S'], default='A', help='モデル渓流 A：床固工群、B：狭窄部と堆砂域、S：直線水路')
    ap.add_argument('--grid', help='格子 CSV（行列形式か x, y, z）')
    ap.add_argument('--cs', type=float, default=2.0, help='格子 CSV の格子間隔 [m]（行列形式のとき）')
    ap.add_argument('--xs', help='横断 CSV（x, s, z）')
    ap.add_argument('--long', help='縦断 CSV（x, z）')
    ap.add_argument('--align', choices=['thalweg', 'left', 'center'], default='thalweg')
    ap.add_argument('--floor', type=float, help='読み込んだ地形で谷底とみなす高さ [m]（既定 3。bed_csv の CSV はファイルの値）')
    ap.add_argument('--dams', help='堰堤 "x:堤高:水通し幅,..."（水通し幅 0 で全幅）')
    ap.add_argument('--fill', type=float, help='はじめから満砂にする（堆砂勾配と元の勾配の比。0 で空）')
    ap.add_argument('--dx', type=float, help='格子幅 [m]')
    ap.add_argument('--qpeak', type=float, help='ピーク流量 [m3/s]')
    ap.add_argument('--qbase', type=float, help='基底流量 [m3/s]')
    ap.add_argument('--minutes', type=float, help='出水の継続時間 [分]（ピークは 1/3 の時刻）')
    ap.add_argument('--cs-feed', type=float, help='供給土砂濃度 Qs/Q [%%]')
    ap.add_argument('--morfac', type=float, help='河床変動の加速係数')
    ap.add_argument('--noneq', action='store_true', help='非平衡の流砂モデル（侵食・堆積速度式。高橋・里深 1999）を使う')
    ap.add_argument('--bank-h', type=float, help='侵食される側岸の高さ [m]（谷底より上の壁をこの高さまで掘れるようにする。0 で側岸は固定）')
    ap.add_argument('--dry-ero', type=float, help='水際の乾いた側岸を削る割合（0〜1。Delft3D の dry-cell erosion）')
    ap.add_argument('--repose-wet', type=float, help='水中の安息角 [度]（既定 30）')
    ap.add_argument('--repose-dry', type=float, help='陸上の安息角 [度]（既定 50）')
    ap.add_argument('--collapse', choices=['repose', 'hasegawa'], help='側岸の崩れ方（安息角で少しずつ／長谷川の崩落モデル。長田ら 1996）')
    ap.add_argument('--secondary', action='store_true', help='2次流で流砂の向きを内岸へ曲げる（長田ら 1996 式14）')
    ap.add_argument('--secondary-mom', action='store_true', help='2次流による運動量の分散項を加える（Kalkwijk・de Vriend 1980）')
    ap.add_argument('--boulders', type=float, help='岸に置く動かない巨礫の間隔 [m]（高橋・里深・奥村 1996。0 で置かない）')
    ap.add_argument('--boulder-mode', choices=['parallel', 'zigzag'], help='巨礫の並べ方（両岸で向かい合わせ／千鳥）')
    ap.add_argument('--out', default='result', help='図と CSV を書き出すフォルダ')
    a = ap.parse_args(argv)
    P = preset(a.preset)
    if a.dx:
        P['dx'] = a.dx
    if a.grid or a.xs:
        if a.grid:
            T = read_grid(a.grid, a.cs)
        else:
            T = read_survey(a.xs, a.long, P['dx'], a.align)
        P['dams'] = []
        use_terrain(P, T, a.floor if a.floor is not None else None)
    if a.dams is not None:
        P['dams'] = _parse_dams(a.dams)
    if a.fill is not None:
        P['prefill'] = a.fill
    H = P['hydro']
    if a.qpeak:
        H['Qpeak'] = a.qpeak
    if a.qbase:
        H['Qbase'] = a.qbase
    if a.minutes:
        H['duration'] = a.minutes * 60; H['tPeak'] = a.minutes * 20
    if getattr(a, 'cs_feed', None) is not None:
        H['Cs'] = a.cs_feed / 100
    if a.morfac:
        P['morfac'] = a.morfac
    if a.noneq:
        P['transport'] = 'noneq'
    B = P.setdefault('bank', {})
    if a.bank_h is not None:
        B['erodH'] = a.bank_h
    if a.dry_ero is not None:
        B['dryEro'] = a.dry_ero
    if a.repose_wet is not None:
        B['reposeWet'] = a.repose_wet
    if a.repose_dry is not None:
        B['reposeDry'] = a.repose_dry
    if a.collapse:
        B['collapse'] = a.collapse
    if a.secondary or a.secondary_mom:
        SC = P.setdefault('sec', {}); SC['bed'] = bool(a.secondary); SC['mom'] = bool(a.secondary_mom)
    if a.boulders is not None:
        P.setdefault('boulders', {})['L'] = a.boulders
    if a.boulder_mode:
        P.setdefault('boulders', {})['mode'] = a.boulder_mode
    sim = Sim(P)
    print('格子 %d × %d（%d セル、幅 %g m）、堰堤 %d 基、%s' % (sim.nx, sim.ny, sim.N, sim.dx, len(sim.dams),
          '読み込んだ地形' if sim.fromCSV else 'モデル渓流 ' + a.preset))
    for w in sim.warn:
        print('注意：' + w)
    run(sim)
    print_dam_table(sim)
    os.makedirs(a.out, exist_ok=True)
    plot_all(sim, os.path.join(a.out, 'result.png'))
    sim.bed_csv(os.path.join(a.out, 'bed_after.csv'))
    print('書き出し：%s/result.png（図）、%s/bed_after.csv（計算後の河床。ブラウザ版の「格子の CSV」で読める）' % (a.out, a.out))


if __name__ == '__main__':
    main()
