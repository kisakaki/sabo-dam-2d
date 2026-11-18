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
import os
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
        # 支川：{'x': 合流点, 'side': 'L'|'R', 'B': 幅, 'I': 勾配 1/n, 'len': 計算区間に入れる長さ, 'r': 流量比, 'rs': 土砂量比, 'lag': 遅れ [s]}
        'tribs': [],
        # 伏流（河床への浸透と地下水の流れ。Darcy 則・Dupuit 近似）。k = 0 なら粒度から Kozeny–Carman 式
        'gw': {'on': False, 'k': 0.0, 'kMax': 0.05, 's0': 0.0},
        # ---- DEM の流域（地形を DEM から切り出したとき）
        'dem': {'soilD': 1.0, 'chW': 6.0, 'ks': 1e-5, 's0': 0.3, 'duration': 3600.0},
        'sources': [],      # 点の流入 [{x, y, ux, uy, I, hydro}]
        'latIn': None,      # 雨の横流入 {dt, pts: [[x, y, ux, uy, I]], Q: [[...]]}
        'rain': {'on': False, 'method': 'lumped', 'mode': 'uniform', 'table': [], 'frames': [], 'loss': {'type': 'f', 'f': 0.8}},
        'slides': [],       # 斜面崩壊
        'slideChk': 1.0,
        'ldams': [],        # 天然ダム
        'ldam': {'k': 0.0, 'fail': 'auto', 'slipLevel': 0.9, 'slipDeg': 20.0, 'thDep': 10.0, 'phi': 35.0, 'c': 0.0, 'gs': 20.0, 'Ts': 0.5, 's0': 0.3},
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
        'hiding': 'egiazaroff',  # 'egiazaroff'（選択的な輸送） | 'atm'（全粒径が同じ割合で動く）
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
        'transport': 'eq',     # 'eq' 平衡流砂量 | 'noneq' 侵食・堆積速度式 | 'debris' 土石流（高橋の式）
        # 土石流（高橋 2007。Kanako 2D と同じ形の式）
        'debris': {'de': 0.0007, 'dd': 0.05, 'tanPhi': 0.7, 'rho': 1000.0, 'Vmin': 0.1, 'manningFloor': True},
        # 透過型堰堤の共通設定（部材の幅、閉塞の判定 l ≤ beta·d95、閉塞する濃度、通過粒径 d < l/alpha、閉塞後の開口率）
        'openDam': {'w': 0.5, 'beta': 1.5, 'Cblk': 0.2, 'alpha': 1.0, 'phiBlock': 0.05},
        'noneq': {'de': 0.6, 'dd': 0.6, 'hEdge': 0.0, 'fe': 2.0, 'fd': 10.0},
        # ---- 巨礫（動かない礫）を水みちの両岸に置く（高橋・里深・奥村 1996）
        'boulders': {'L': 0.0, 'mode': 'parallel', 'size': 1},
        'activeLayer': 0.0,    # 交換層厚（0 のとき最大粒径）
        'layerTh': 0.25,
        'morfac': 4.0,         # 河床変動の加速係数
        'snapEvery': 600.0,    # 時間ごとの結果を記録する間隔 [s]（出水の時刻）。0 で記録しない
        'gridDamp': 0.3,
        'gravK': 1.0,
        # ---- ハイドログラフ
        'hydro': {'mode': 'tri', 'Qbase': 10.0, 'Qpeak': 60.0, 'tPeak': 3600.0, 'duration': 10800.0,
                  'Cs': 0.005, 'table': [],    # table: [[t(s), Q, Qs(省略可)], ...]
                  # mode 'debris'：基底流量に三角形の土石流の波（t0 から立ち上がり tRise、減衰 tFall、濃度 Cd）を重ねる。Iup：流入点より上流の勾配 1/n
                  't0': 60.0, 'tRise': 30.0, 'tFall': 120.0, 'Cd': 0.4, 'Iup': 0.0},
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


def _wave_qh(H, t):
    a = t - (H.get('t0') or 0.0); r = max(1.0, H.get('tRise') or 30.0); f = max(1.0, H.get('tFall') or 120.0)
    if a <= 0 or a >= r + f:
        return 0.0
    return (H['Qpeak'] - H['Qbase']) * (a / r if a < r else 1 - (a - r) / f)


def QatH(H, t):
    """ハイドロ H（上流端と同じ形の設定）の時刻 t の流量 [m3/s]"""
    if H['mode'] == 'debris':
        return H['Qbase'] + _wave_qh(H, t)
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


def QsAtH(H, t):
    """供給土砂量（実質体積） [m3/s]"""
    tab = H.get('table') or []
    if H['mode'] == 'table' and len(tab) > 1 and len(tab[0]) > 2 and np.isfinite(tab[0][2]):
        return max(0.0, _interp_tab(tab, t, 2))
    if H['mode'] == 'debris':
        return (H.get('Cd') or 0.0) * _wave_qh(H, t) + (H.get('Cs') or 0.0) * H['Qbase']
    return (H.get('Cs') or 0.0) * QatH(H, t)


def _rain_uni(tab, t):
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


def _rain_frame(frames, t):
    if not frames:
        return None
    tm = t / 60; f = None
    for fr in frames:
        if tm >= fr['t'] - 1e-9:
            f = fr
        else:
            break
    if f is not None and len(frames) > 1:
        dtf = frames[-1]['t'] - frames[-2]['t']
        if tm >= frames[-1]['t'] + dtf:
            f = None
    return f


def _loss_cell(L, Sa, Ra, c, r, dt):
    """損失（有効雨量）：流出係数 f、一次流出率 f1・飽和雨量 Rsa・初期損失 R0、浸透能 ks と土層の空き容量"""
    dth = dt / 3600
    if L.get('type') == 'f1rsa':
        Ra[c] += r * dth
        R0 = L.get('R0') or 0.0
        f = 0.0 if Ra[c] < R0 else ((L['f1'] if L.get('f1') is not None else 0.5) if Ra[c] < R0 + (L.get('Rsa') or 0.0) else 1.0)
        return f * r
    if L.get('type') == 'ks':
        Smax = max(0.0, (L.get('D') or 1.0) * ((L.get('n') or 0.45) - (L.get('th0') or 0.25)) * 1000); ks = (L.get('ks') or 1e-5) * 3.6e6
        if Sa[c] >= Smax:
            return r
        inf = min(r, ks); add = inf * dth
        if Sa[c] + add > Smax:
            over = (Sa[c] + add - Smax) / dth; Sa[c] = Smax
            return r - inf + over
        Sa[c] += add
        return r - inf
    return (L['f'] if L.get('f') is not None else 0.8) * r


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
        self.tribs = []
        if T is not None and P.get('tribs'):
            self.warn.append('支川は、数値で指定した地形でだけ使えます（CSV の地形では使いません）。')
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
            nsL = nsR = ns         # 支川があれば、その側の谷壁の外側を支川の長さだけ広げる
            for tb in (P.get('tribs') or []):
                nn = math.ceil(max(10, tb.get('len') or 40) / dx) + 1
                if tb.get('side') == 'R':
                    nsR = max(nsR, nn)
                else:
                    nsL = max(nsL, nn)
            ny = _jround((ymax - ymin) / dx) + nsL + nsR
            yc = nsR * dx + (yo - ymin)
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
        self.Ix = np.zeros(nx)          # 各列の基準の縦断勾配（勾配の表があれば区間ごと）
        self.zbAt = None
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
            self.Ix[:] = self.I
            trA_, trB_ = self.trA, self.trB
            self.zbAt = lambda x: trA_ + trB_ * x
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
            if T.get('basin') is not None and np.shape(T['basin']) == np.shape(T['rows']):
                self._setup_basin(T, z, zfix, kind, fl, colmin, Bx)
        else:
            # 縦断：勾配の表（位置 x から先の勾配 1/n）があれば、区間ごとの勾配をつないだ折れ線にする
            st = sorted([list(r) for r in (P.get('slopeTable') or []) if np.isfinite(r[0]) and r[1] > 0], key=lambda r: r[0])
            Lx_ = self.Lx
            if st and st[0][0] > 0:
                st.insert(0, [0.0, 1.0 / P['slope']])
            brk = [max(0.0, min(Lx_, r[0])) for r in st] + [Lx_]

            def zbAt(x):
                if not st:
                    return (Lx_ - x) * P['slope']
                xa = max(0.0, min(Lx_, x)); zz_ = 0.0
                for q2 in range(len(st)):
                    a0, b0 = brk[q2], brk[q2 + 1]
                    if b0 <= xa:
                        continue
                    zz_ += (b0 - max(a0, xa)) / st[q2][1]
                return zz_

            def segI(x):
                v = P['slope']
                for r in st:
                    if x >= r[0] - 1e-9:
                        v = 1.0 / r[1]
                return v
            self.zbAt = zbAt
            I0 = self.I = max(1e-4, zbAt(0) / Lx_) if st else P['slope']
            self.trA = self.Lx * I0; self.trB = -I0
            self.Ix[:] = [segI(x) for x in X] if st else P['slope']
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
            zc = np.array([zbAt(x) for x in X]) if st else (self.Lx - X) * I0

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
            # 支川：谷壁を切り込んだ水路（勾配 I）を、計算区間の側方の端から本川の谷底の縁までつなぐ
            TR = [tb for tb in (P.get('tribs') or []) if np.isfinite(tb.get('x', np.nan)) and (tb.get('B') or 0) > 0 and (tb.get('I') or 0) > 0]
            for tb in TR:
                # 合流角 ang：本川の流下方向と支川の流下方向のなす角（90° で直角、小さいほど本川の下流向きに斜めに入る）
                sdT = -1 if tb.get('side') == 'R' else 1; It = 1.0 / tb['I']
                angD = tb.get('ang'); angD = 90.0 if angD is None or not np.isfinite(angD) else max(30.0, min(150.0, float(angD)))
                aR = angD * math.pi / 180; sA = 1.0 if angD == 90 else math.sin(aR); cA = 0.0 if angD == 90 else math.cos(aR)
                cot = cA / sA; hw = tb['B'] / 2; xj = tb['x']
                jbT = ny - 1 if sdT > 0 else 0
                inC = []
                zE = np.array([zc[i] + cs0 * abs((hL[i] if sdT > 0 else -hR[i]) - ych[i]) for i in range(nx)])
                for i in range(2, nx - 2):
                    offT = (Y - yc[i]) * cosT[i]
                    dist = offT - hL[i] if sdT > 0 else -hR[i] - offT
                    xr = (i + 0.5) * dx - xj
                    pp = xr * sA + dist * cA                    # 支川の中心線からの距離（中心線に直角）
                    iw = np.clip(np.floor(((i + 0.5) * dx + dist * cot) / dx), 0, nx - 1).astype(int)
                    zb = zE[iw] + dist / sA * It
                    ch = (dist > 0) & (pp >= -hw) & (pp < hw)
                    z[ch, i] = zb[ch]; kind[ch, i] = 0; fl[ch, i] = True; self.bare[ch, i] = False
                    if ch[jbT]:
                        inC.append(jbT * nx + i)
                    # 支川の谷：支川の河床より低い谷壁を、水路から離れるほど高くなる V 字の谷（側壁の勾配 1:2）に盛る
                    lat = np.abs(pp) - hw
                    ztb = zb + np.minimum(1.5, dist * 0.5) + 0.5 * lat
                    m = (dist > 0) & ~ch & (kind[:, i] == 1) & (lat < 30) & (z[:, i] < ztb)
                    z[m, i] = ztb[m]
                if not inC:
                    self.warn.append('支川（x = %g m）が計算区間の側方の端に届きません。合流角か位置を変えてください。' % tb['x'])
                if inC:
                    self.tribs.append(dict(x=tb['x'], side=sdT, I=It, B=tb['B'], ang=angD, sA=sA, cA=cA, r=tb.get('r') or 0.0,
                                           rs=tb['rs'] if (tb.get('rs') is not None and tb['rs'] >= 0) else (tb.get('r') or 0.0),
                                           lag=tb.get('lag') or 0.0, cells=np.array(inC), q=0.0, v=0.0, Q=0.0, Qs=0.0))
            zfix = np.where(kind == 0, np.where(fixedEnd, z, z - P['erodDepth']), z)
            colmin = zc.copy()
            if self.shaped:
                colmin = np.where(fl.any(axis=0), np.where(fl, z, np.inf).min(axis=0), colmin)
        I = self.I

        # ---- 堰堤（河道の中心線に直角＝両岸に直角）
        self.dams = []
        self.openC = np.zeros((ny, nx), dtype=np.int16)    # 透過型堰堤の開口のセル（堰堤の番号 + 1）
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
            zbase = colmin[idd] if T is not None else self.zbAt(cx)
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
            zlow_ = float(z[cells].min())
            # ---- 透過型：水通しの中に開口を開ける（開口のセルは開口の底の高さの移動床）
            OD = P.get('openDam') or {}; typ = D.get('type') or 'closed'; op = None
            if typ != 'closed':
                lw = D['l'] if (D.get('l') or 0) > 0 else 1.0; wm = OD.get('w') or 0.5; nO = max(1, _jround(D.get('n') or 1))
                zo = zbase + max(0.0, min(D['H'] - 0.1, D.get('hb') or 0.0)); sc0 = 0.5 * (spMin + spMax)
                if typ == 'slit':
                    nO = 1; wz = lw
                elif typ == 'mslit':
                    wz = nO * lw + (nO - 1) * wm
                else:
                    wz = spMax - spMin + dx
                ocm = spc & (np.abs(sl_ - sc0) <= wz / 2 + 1e-9)
                if not ocm.any() and spc.any():      # 開口が格子より狭いときは、水通しの中心に最も近いセルを 1 つ開ける
                    fi = np.flatnonzero(spc); dd_ = np.abs(sl_.ravel()[fi] - sc0)
                    ocm = np.zeros_like(spc); ocm.ravel()[fi[int(np.argmin(dd_))]] = True
                if ocm.any():
                    cover = ocm.sum() * dx / max(1, _jround(thick * 2 / dx))
                    if typ in ('slit', 'mslit'):
                        phiW = min(1.0, nO * lw / max(cover, 1e-9)); phiS = phiW
                    else:
                        fr1 = lw / (lw + wm); phiW = fr1 * fr1 if typ == 'grid' else fr1; phiS = fr1
                    z = np.where(ocm, zo, z); zfix = np.where(ocm, zo, zfix); kind = np.where(ocm, 0, kind).astype(np.uint8)
                    fl = fl | ocm; self.openC[ocm] = len(self.dams) + 1
                    cells = cells & ~ocm; spc = spc & ~ocm
                    zs0 = sc0 - wz / 2 if typ in ('slit', 'mslit') else spMin - dx / 2
                    zs1 = sc0 + wz / 2 if typ in ('slit', 'mslit') else spMax + dx / 2
                    op = dict(type=typ, l=lw, n=nO, w=wm, zo=zo, cells=np.flatnonzero(ocm), phiW=phiW, phiS=phiS, blocked=False,
                              tBlock=None, beta=D['beta'] if (D.get('beta') or 0) > 0 else (OD.get('beta') or 1.5), s0=zs0, s1=zs1, hb=zo - zbase,
                              Cblk=None, d95=None)
            self.dams.append(dict(i=idd, x=cx, cx=cx, cy=cy, th=th, tx=tx, ty=ty, nx=nxv, ny=nyv, ext=ext, thick=thick,
                                  H=D['H'], b=bn, crest=crest, wingTop=wingTop, zbase=zbase, zlow=zlow_,
                                  cells=cells, spill=spc, sMin=sv_.min() - dx / 2, sMax=sv_.max() + dx / 2,
                                  spMin=spMin - dx / 2, spMax=spMax + dx / 2, type=typ if op else 'closed', op=op, hb=D.get('hb') or 0.0))
        nd = len(self.dams)
        self.hasOpen = any(Dd['op'] is not None for Dd in self.dams)
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
                tq = Dq['crest'] + upd * self.Ix[max(0, Dq['i'] - 1)] * P['prefill']
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
        dl = np.where(kind.ravel() == 3, 0.0, zf - zx)   # 斜面の土層（kind 3）は流れで侵食しない（層を持たない）
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
        if not hasattr(self, 'demMode'):
            self.demMode = False
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
            cr = np.array([(-np.inf if D['op'] is not None else D['crest']) for D in self.dams] + [-np.inf])[region]   # 透過型の上流は湛水させない
            m = (region < nd) & (kind != 2) & (self.openC == 0) & (z < cr) & (cr - z > h)
            h = np.where(m, cr - z, h); qx = np.where(m, 0.0, qx); qy = np.where(m, 0.0, qy)
        self.h, self.qx, self.qy = h, qx, qy

        # 支川の流入セル
        self.tribIn = np.zeros(N, dtype=np.int16)
        for q_, tb in enumerate(self.tribs):
            self.tribIn[tb['cells']] = q_ + 1
        self._setup_points()
        # 伏流：河床（堆積層）の中の水 W [m]。水に浸かったセルは飽和、乾いたセルは飽和度 s0
        GW0 = P.get('gw') or {}
        self.W = np.zeros((ny, nx)); self.kgw = np.zeros((ny, nx)); self.gwOut = 0.0
        self.por = (self.kind == 0) | (self.kind == 3)   # 伏流の層があるセル（河床と斜面の土層）
        if GW0.get('on'):
            Tc = np.maximum(0.0, self.z - self.zfix)
            DM0 = P.get('dem') or {}
            s0h = DM0['s0'] if (DM0.get('s0') is not None and DM0['s0'] >= 0) else (GW0.get('s0') or 0.0)
            sat = np.where(self.h > 1e-3, 1.0, np.where(self.kind == 3, min(1.0, max(0.0, s0h)), min(1.0, max(0.0, GW0.get('s0') or 0.0))))
            self.W = np.where(self.por, P['porosity'] * Tc * sat, 0.0)
        self.W0 = float(self.W.sum()) * dx * dx
        self.voidAir = 0.0   # 伏流なしの土石流で、堆積物の間隙を満たす水が足りなかった体積（空気）
        self.hOld = np.zeros((ny, nx)); self.dfC = np.zeros((ny, nx)); self.dfD = np.zeros((ny, nx))
        self.dzs = np.zeros((ny, nx))
        self.rainVol = 0.0; self.rainLoss = 0.0; self.slideLog = []
        self._setup_slides()
        # ---- 時刻・収支
        self.t = 0.0; self.th = 0.0; self.nstep = 0; self.dt = 0.01
        self.phase = 'spinup' if Q0 > 0.01 else 'run'
        self.spinT = 0.0
        self.Qin = Q0; self.Qout = 0.0; self.QsIn = 0.0; self.QsOut = 0.0
        self.Vwin = 0.0; self.Vwout = 0.0; self.Vw0 = self.water_volume()
        self.sedIn = 0.0; self.sedOut = 0.0; self.sedOutK = np.zeros(K)
        self.nLimW = 0; self.nLimS = 0
        self.series = {'t': [], 'Qin': [], 'Qout': [], 'QsIn': [], 'QsOut': [], 'stor': []}
        self.snaps = []; self.snapFull = False; self._nextSnap = 0.0   # 時間ごとの結果（記録は出水の開始から）
        self._nextSample = 0.0
        self.Uref = max(0.5, (Q0 / self.Bf) / ((Q0 / self.Bf) * P['manning'] / math.sqrt(max(I, 1e-4))) ** 0.6) if Q0 > 0.01 else 1.0
        self.FxM = np.zeros((ny, nx + 1)); self.GyM = np.zeros((ny + 1, nx))
        self.qin = np.zeros(ny); self.uin = np.zeros(ny)
        # 天然ダムを置いたあとの地形を初期の地形とする
        self.z0 = self.z.copy(); self.Vw0 = self.water_volume(); self.W0 = self.ground_water()

    # ------------------------------------------------------------ ハイドログラフ
    def duration(self):
        if getattr(self, 'demMode', False) and ((self.P.get('dem') or {}).get('duration') or 0) > 0:
            return float(self.P['dem']['duration'])
        H = self.P['hydro']
        if H['mode'] == 'table' and len(H.get('table') or []) > 1:
            return H['table'][-1][0]
        return H['duration']

    def waveQ(self, t):
        """土石流の一波（三角形）の流量 [m3/s]（混合物）"""
        H = self.P['hydro']; a = t - (H.get('t0') or 0.0); r = max(1.0, H.get('tRise') or 30.0); f = max(1.0, H.get('tFall') or 120.0)
        if a <= 0 or a >= r + f:
            return 0.0
        return (H['Qpeak'] - H['Qbase']) * (a / r if a < r else 1 - (a - r) / f)

    def Qat(self, t):
        return 0.0 if getattr(self, 'demMode', False) else QatH(self.P['hydro'], t)

    def QsAt(self, t):
        """供給土砂量（実質体積） [m3/s]"""
        return 0.0 if getattr(self, 'demMode', False) else QsAtH(self.P['hydro'], t)

    # ------------------------------------------------------------ 1 ステップ
    def step(self):
        P = self.P; dx = self.dx; nx, ny = self.nx, self.ny
        h, qx, qy, z = self.h, self.qx, self.qy, self.z
        HDRY, HS, UMAX, g = 1e-3, 0.01, 30.0, G
        hg = 0.5 * g
        second = P['order'] != 1

        # ---- 流入境界（上流端の谷底セルに流量を配分。水深の 5/3 乗で重み付け）
        Q = self.Qat(self.t); self.Qin = Q
        nman = P['manning']
        w = np.where(self.fl[:, 0], (h[:, 0] + 0.05) ** (5 / 3), 0.0)
        qin = Q * w / w.sum() / dx if w.sum() > 0 else np.zeros(ny)
        dbIn = P.get('transport') == 'debris'
        Cin0 = min(0.9 * (1 - P['porosity']), self.QsAt(self.t) / Q) if (dbIn and Q > 0) else 0.0
        Hd = P['hydro']; Iu = 1.0 / Hd['Iup'] if (Hd.get('Iup') or 0) > 0 else self.Ix[0]
        if dbIn:
            uin = np.array([qj / self.debris_depth(qj, Cin0, Iu) if qj > 0 else 0.0 for qj in qin])
        else:
            uin = np.where(qin > 0, qin / np.maximum(qin * nman / math.sqrt(max(Iu, 1e-4)), 1e-30) ** 0.6, 0.0)
        self.qin, self.uin = qin, uin
        # 支川の流入（本川の流量の r 倍、lag 秒遅れ）
        QtSum = 0.0
        for tb in self.tribs:
            tl = max(0.0, self.t - tb['lag']); Qt = tb['r'] * self.Qat(tl)
            tb['Q'] = Qt; tb['q'] = Qt / (len(tb['cells']) * dx); QtSum += Qt
            Ct = min(0.9 * (1 - P['porosity']), tb['rs'] * self.QsAt(tl) / Qt) if (dbIn and Qt > 0 and self.phase == 'run') else 0.0
            qw = tb['q'] * tb['sA']   # 支川の単位幅流量（斜めに入ると、端のセルの列は幅 B/sin になる）
            tb['v'] = (qw / self.debris_depth(qw, Ct, tb['I']) if dbIn else qw / (qw * nman / math.sqrt(tb['I'])) ** 0.6) if qw > 0 else 0.0

        # 点の流入（谷頭の供給点・雨の横流入）
        QpSum = self._pin_step() if self.pin else 0.0
        rainH = None

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
            dex[:, -1] = np.where(oke, _minmod(eta[:, -1] - eta[:, -2], -self.Ix[-1] * dx), 0.0)
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
        if self.rainOn:
            rainH = self._rain_step(dt)   # 格子に降る雨（このステップの有効雨量の水深 [m]）
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

        # ---- 土石流：濃度と代表粒径（抵抗則に使う）
        dbr = P.get('transport') == 'debris'
        if dbr:
            self.hOld = h.copy(); self._debris_params()
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
        for tb in self.tribs:
            jj_, ii_ = np.divmod(tb['cells'], nx)
            hn[jj_, ii_] += cf * tb['q']; qyn[jj_, ii_] -= cf * tb['q'] * tb['v'] * tb['sA'] * tb['side']; qxn[jj_, ii_] += cf * tb['q'] * tb['v'] * tb['cA']
        if self.pin:
            pj, pi = np.divmod(self.pinCells, nx)
            hn[pj, pi] += dt * self.pH[self.pinCells]; qxn[pj, pi] += dt * self.pMx[self.pinCells]; qyn[pj, pi] += dt * self.pMy[self.pinCells]
        if rainH is not None:
            hn = np.where(rainH > 0, hn + rainH, hn)
        hn = np.maximum(hn, 0.0)
        wetn = hn > HDRY
        hs = np.where(wetn, hn, 1.0)
        Vv = np.sqrt(qxn * qxn + qyn * qyn)
        frX = 4 * dt * g * self.n2 * Vv / (hs * hs * np.cbrt(hs))
        if dbr:
            # 高橋の抵抗則：石礫型（C ≥ 0.4C*）Ψ = d²/{8 h² (C + (1−C)ρ/σ)((C*/C)^(1/3) − 1)²}、乱流型 Ψ = d²/(0.49 h²)
            DB = P.get('debris') or {}; Cst_ = 1 - P['porosity']; rhoS = 1000 * P['sg']; rhoW = DB.get('rho') or 1000.0
            Cd_ = np.minimum(self.dfC, 0.98 * Cst_); dd_ = self.dfD; use = self.dfC >= 0.01
            Cds = np.where(use, Cd_, 0.5 * Cst_)
            X_ = np.cbrt(Cst_ / Cds) - 1
            psi = np.where(Cds >= 0.4 * Cst_, dd_ * dd_ / (8 * hs * hs * (Cds + (1 - Cds) * rhoW / rhoS) * X_ * X_), dd_ * dd_ / (0.49 * hs * hs))
            frD = 4 * dt * psi * Vv / (hs * hs)
            frX = np.where(use, np.maximum(frX, frD) if DB.get('manningFloor', True) is not False else frD, frX)
        if self.hasOpen:
            # 透過型堰堤の開口：多孔板の損失係数（Idelchik）ζ = {(1 + 0.707√(1−φ) − φ)/φ}²。天端より上を流れる分は損失なし
            phiB = (P.get('openDam') or {}).get('phiBlock') or 0.05
            oc = self.openC > 0
            for q_, Dq in enumerate(self.dams):
                op = Dq['op']
                if op is None:
                    continue
                jj_, ii_ = np.divmod(op['cells'], nx)
                ph = phiB if op['blocked'] else op['phiW']
                hq = hs[jj_, ii_]; Hs_ = z[jj_, ii_] + hq
                hsub = np.maximum(0.0, np.minimum(Hs_, Dq['crest']) - z[jj_, ii_]); hov = np.maximum(0.0, Hs_ - Dq['crest'])
                pe = (ph * hsub + hov) / hq
                pe2 = np.maximum(pe, 1e-3); zt_ = (1 + 0.707 * np.sqrt(np.maximum(0.0, 1 - pe2)) - pe2) / pe2
                add = np.where(pe < 1, 4 * dt * (zt_ * zt_ / (2 * dx * hq)) * Vv[jj_, ii_], 0.0)
                frX[jj_, ii_] = frX[jj_, ii_] + add
        fr = 0.5 * (1 + np.sqrt(1 + frX))
        self.h = hn
        self.qx = np.where(wetn, qxn / fr, 0.0)
        self.qy = np.where(wetn, qyn / fr, 0.0)
        self.u, self.v = u, v
        self.FxM, self.GyM = FxM, GyM
        self.Qout = FxM[:, nx].sum() * dx
        self.Vwin += (Q + QtSum + QpSum) * dt + (self._rainV if rainH is not None else 0.0); self.Vwout += self.Qout * dt
        self.th += dt; self.nstep += 1
        gwOn_ = bool((P.get('gw') or {}).get('on'))

        # ---- 助走計算（水だけを流して流れを立ち上げる）
        if self.phase == 'spinup':
            self.spinT += dt
            tmin = 0.5 * self.Lx / self.Uref; tmax = 3 * self.Lx / self.Uref
            if (self.spinT > tmin and abs(self.Qout - Q) < 0.08 * Q) or self.spinT > tmax:
                self.phase = 'run'; self.Vwin = 0.0; self.Vwout = 0.0; self.Vw0 = self.water_volume()
                self.W0 = self.ground_water(); self.gwOut = 0.0; self.voidAir = 0.0
                self.snapshot()
            if gwOn_:
                self._gw_step(dt)
            return
        if self.phase == 'done':
            return

        # ---- 河床変動
        if self.t == 0 and not self.snaps:
            self.snapshot()             # 出水の開始時の状態
        dtm = dt * (1.0 if P.get('transport') == 'debris' else max(1.0, P['morfac']))   # 土石流は加速しない
        self._morph(dtm)
        if gwOn_:
            self._gw_step(dt)        # 伏流は流れの時間刻みで（加速しない）
        if self.slides or self.depOn:
            self._slide_step(dt)
        self.t += dtm
        if self.t >= self._nextSample:
            S = self.series
            S['t'].append(self.t); S['Qin'].append(Q); S['Qout'].append(self.Qout)
            S['QsIn'].append(self.QsIn); S['QsOut'].append(self.QsOut); S['stor'].append(self.storage())
            self._nextSample = self.t + self.duration() / 400
        if self.t >= self.duration():
            self.phase = 'done'
        ev = float(P.get('snapEvery') or 0)
        if ev > 0 and (self.t >= self._nextSnap - 1e-9 or (self.phase == 'done' and self.snaps and self.snaps[-1]['t'] < self.t)):
            self.snapshot()

    # ------------------------------------------------------------ 時間ごとの結果の記録
    SNAP_MAX = 25e6     # 記録する値の総数の上限（約 100 MB）

    def snapshot(self):
        """河床変動量・水深・流速・表層の平均粒径を 4 バイトの実数で残す（河床高は z0 + 河床変動量で戻す）"""
        ev = float(self.P.get('snapEvery') or 0)
        if not ev > 0:
            return
        self._nextSnap = (math.floor(self.t / ev + 1e-9) + 1) * ev
        if self.snaps and abs(self.snaps[-1]['t'] - self.t) < 1e-9:
            return
        if (len(self.snaps) + 1) * self.N * 4 > self.SNAP_MAX:
            self.snapFull = True
            return
        h = self.h; w = h > 1e-3; hs = np.where(w, h, 1.0)
        u = 0.5 * (self.FxM[:, :-1] + self.FxM[:, 1:]) / hs; v = 0.5 * (self.GyM[:-1] + self.GyM[1:]) / hs
        dm = np.where(self.kind == 0, self.surface_dm(), np.nan)
        self.snaps.append({'t': float(self.t), 'dz': (self.z - self.z0).astype(np.float32), 'h': h.astype(np.float32),
                           'V': np.where(w, np.sqrt(u * u + v * v), 0.0).astype(np.float32), 'dm': dm.astype(np.float32),
                           'Qin': float(self.Qin), 'Qout': float(self.Qout), 'QsIn': float(self.QsIn), 'QsOut': float(self.QsOut),
                           'sedIn': float(self.sedIn), 'sedOut': float(self.sedOut), 'stor': float(self.storage()), 'byDam': self.storage_by_dam().copy()})

    def set_snap_every(self, sec):
        """記録の間隔を途中で変える（次の記録は新しい間隔の倍数の時刻）"""
        self.P['snapEvery'] = sec if sec > 0 else 0.0
        self._nextSnap = (math.floor(self.t / sec + 1e-9) + 1) * sec if sec > 0 else float('inf')
        if sec > 0 and self.phase != 'spinup' and not self.snaps:
            self.snapshot()

    # ------------------------------------------------------------ 限界掃流力
    def _flow_slope(self, ux, uy):
        """流れの向きの河床勾配 i（下り勾配を正、0〜0.2）。堰堤のセルは除く"""
        z, kind, dx = self.z, self.kind, self.dx
        nd = kind != 2
        zl = np.empty_like(z); zl[:, 1:] = z[:, :-1]; zl[:, 0] = z[:, 0]
        zr = np.empty_like(z); zr[:, :-1] = z[:, 1:]; zr[:, -1] = z[:, -1]
        al = np.zeros_like(nd); al[:, 1:] = nd[:, :-1]
        ar = np.zeros_like(nd); ar[:, :-1] = nd[:, 1:]
        gx = np.where(al & ar, (zr - zl) / (2 * dx), np.where(al, (z - zl) / dx, np.where(ar, (zr - z) / dx, -self.Ix[None, :])))
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
        if self.P.get('transport') == 'debris':
            self._morph_debris(dtm)
        elif self.P.get('transport') == 'noneq':
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

        if self.hasOpen:
            self._check_block(qb=(qbx, qby)); self._open_filter(SFx, SFy)
        # ---- 粒径別の持ち出し制限（交換層内の量を超えない）
        cf = dtm / ((1 - P['porosity']) * dx)
        outv = cf * (np.maximum(SFx[:, 1:], 0) - np.minimum(SFx[:, :-1], 0) + np.maximum(SFy[1:], 0) - np.minimum(SFy[:-1], 0))
        with np.errstate(over='ignore', divide='ignore', invalid='ignore'):   # 選ばれない側の割り算の警告を出さない
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
        self._trib_feed(dtm, True)

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
        if self.hasOpen:
            self._check_block(Cc=Cc); self._open_filter(SFx, SFy)
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
        self._trib_feed(dts, False)

    # ------------------------------------------------------------ 透過型堰堤（開口を通る土砂と閉塞）
    def _open_filter(self, SFx, SFy):
        """開口のセルの面を通る流砂量に通過率 f_k を掛ける（閉塞後 0、d_k ≥ l/α の礫も 0、それ以外は開口率 φs）"""
        nx, K, d = self.nx, self.K, self.d
        al = (self.P.get('openDam') or {}).get('alpha') or 1.0
        for Dq in self.dams:
            op = Dq['op']
            if op is None:
                continue
            fK = np.array([0.0 if op['blocked'] else (op['phiS'] if d[k] * al < op['l'] else 0.0) for k in range(K)])
            jj, ii = np.divmod(op['cells'], nx)
            fx = set(); fy = set()
            for j_, i_ in zip(jj, ii):
                fx.add((j_, i_)); fx.add((j_, i_ + 1)); fy.add((j_, i_)); fy.add((j_ + 1, i_))
            for (j_, i_) in fx:
                SFx[j_, i_, :] *= fK
            for (j_, i_) in fy:
                SFy[j_, i_, :] *= fK

    def _hold_volume(self, dX, dY, cf):
        """土石流：開口で止めた土砂の体積を、混合物の水深として上流側に残す（下流側へは水だけが通る）"""
        h = self.h; ny, nx = h.shape
        for j, i in zip(*np.nonzero(dX)):
            if i < 1 or i >= nx:
                continue
            t = cf * dX[j, i]
            t = min(t, h[j, i]) if t > 0 else -min(-t, h[j, i - 1])
            h[j, i - 1] += t; h[j, i] -= t
        for j, i in zip(*np.nonzero(dY)):
            if j < 1 or j >= ny:
                continue
            t = cf * dY[j, i]
            t = min(t, h[j, i]) if t > 0 else -min(-t, h[j - 1, i])
            h[j - 1, i] += t; h[j, i] -= t

    def _check_block(self, qb=None, Cc=None):
        """閉塞の判定：開口のすぐ上流の流れの土砂濃度（流砂量 ÷ 流量）が Cblk 以上で、流れている土砂の d95 に対し l ≤ β d95 なら閉塞"""
        nx, ny, K, dx = self.nx, self.ny, self.K, self.dx
        h = self.h.ravel(); kind = self.kind.ravel(); openC = self.openC.ravel()
        FxM = self.FxM.ravel(); GyM = self.GyM.ravel()
        OD = self.P.get('openDam') or {}; Cb = OD.get('Cblk') or 0.2
        ord_ = np.argsort(self.d, kind='stable'); d = self.d
        if qb is not None:
            qbx = qb[0].reshape(-1, K); qby = qb[1].reshape(-1, K)
        else:
            Ccf = Cc.reshape(-1, K)
        for Dq in self.dams:
            op = Dq['op']
            if op is None or op['blocked']:
                continue
            sS = 0.0; sQ = 0.0; tk = [0.0] * K
            for c0 in op['cells']:
                i0 = c0 % nx; j0 = c0 // nx
                for nb in (c0 - 1 if i0 > 0 else -1, c0 + 1 if i0 < nx - 1 else -1, c0 - nx if j0 > 0 else -1, c0 + nx if j0 < ny - 1 else -1):
                    if nb < 0 or openC[nb] or kind[nb] != 0 or h[nb] < 0.02:
                        continue
                    ib = nb % nx; jb = nb // nx
                    if ((ib + 0.5) * dx - Dq['cx']) * Dq['tx'] + ((jb + 0.5) * dx - Dq['cy']) * Dq['ty'] >= 0:
                        continue
                    fw = jb * (nx + 1) + ib
                    uu = 0.5 * (FxM[fw] + FxM[fw + 1]) / h[nb]; vv = 0.5 * (GyM[nb] + GyM[nb + nx]) / h[nb]; V = math.sqrt(uu * uu + vv * vv)
                    sQ += h[nb] * V
                    for k in range(K):
                        qk = Ccf[nb, k] * V if qb is None else math.sqrt(qbx[nb, k] * qbx[nb, k] + qby[nb, k] * qby[nb, k])
                        tk[k] += qk; sS += qk
            if not sQ > 1e-9 or sS / sQ < Cb:
                continue
            cum = 0.0; d95 = d[ord_[K - 1]]
            for k in range(K):
                cum += tk[ord_[k]]
                if cum >= 0.95 * sS:
                    d95 = d[ord_[k]]; break
            if op['l'] <= op['beta'] * d95:
                op['blocked'] = True; op['tBlock'] = self.t; op['Cblk'] = sS / sQ; op['d95'] = d95

    # ------------------------------------------------------------ 支川からの給砂
    def _trib_feed(self, dtm, bulk):
        """本川の供給土砂量の rs 倍（lag 秒遅れ）。bulk：河床に置く（平衡流砂量）、そうでなければ流れの中へ"""
        if self.phase == 'run' and self.pin:
            self._pin_feed(dtm, bulk)
        if not self.tribs or self.phase != 'run':
            return
        K, dx, lam = self.K, self.dx, 1 - self.P['porosity']
        for tb in self.tribs:
            Qs = tb['rs'] * self.QsAt(max(0.0, self.t - tb['lag']))
            if self.P.get('transport') == 'debris' and Qs > 0.9 * lam * tb['Q']:
                Qs = 0.9 * lam * tb['Q']       # 混合物の濃度は 0.9C* まで
            tb['Qs'] = Qs
            if not Qs > 0:
                continue
            self.sedIn += Qs * dtm
            per = Qs * dtm / (len(tb['cells']) * dx * dx)
            if bulk:
                zf = self.z.ravel(); dzf = self.dzs.ravel()
                for c in tb['cells']:
                    self.A[c] += per / lam * self.ff; zf[c] += per / lam; dzf[c] += per / lam
                self._fix_layers(np.array(tb['cells']))
            else:
                for c in tb['cells']:
                    self.Cc[c] += per * self.ff

    # ======================================================================================
    # DEM の流域（sabo2d_dem.prepare_domain で切り出した地形）
    #  流域の外：高い壁。谷筋（幅 chW に広げる）：侵食できる河床（kind 0）。斜面：kind 3（流れでは侵食しない土層）
    # ======================================================================================
    def _setup_basin(self, T, z, zfix, kind, fl, colmin, Bx):
        P = self.P; nx, ny, dx = self.nx, self.ny, self.dx; N = nx * ny
        DM = P.get('dem') or {}
        soilD = DM.get('soilD') if (DM.get('soilD') or 0) > 0 else 1.0
        chR = max(0, _jround(((DM.get('chW') if (DM.get('chW') or 0) > 0 else dx) / dx - 1) / 2))
        cs = T['cs']; rows = np.asarray(T['rows'], dtype=float); nr, nc = rows.shape
        jj, ii = np.mgrid[0:ny, 0:nx]
        r = np.clip(np.floor((self.Ly - (jj + 0.5) * dx) / cs).astype(int), 0, nr - 1)
        q = np.clip(np.floor((ii + 0.5) * dx / cs).astype(int), 0, nc - 1)
        self.tIdx = (r * nc + q).ravel()
        bas = (np.asarray(T['basin'])[r, q] != 0)
        ch = (np.asarray(T['chan'])[r, q] != 0) if T.get('chan') is not None else np.zeros((ny, nx), bool)
        XX = (ii + 0.5) * dx; YY = (jj + 0.5) * dx
        z[:] = _sample_terrain(T, XX, YY, self.Ly)
        zTop = float(z[bas].max()) if bas.any() else float(z.max())
        chD = np.zeros((ny, nx), bool)
        for dj in range(-chR, chR + 1):
            for di in range(-chR, chR + 1):
                sh = np.zeros((ny, nx), bool)
                ys0, ys1 = max(0, -dj), ny - max(0, dj); xs0, xs1 = max(0, -di), nx - max(0, di)
                sh[ys0 + dj:ys1 + dj, xs0 + di:xs1 + di] = ch[ys0:ys1, xs0:xs1]
                chD |= sh
        chD &= bas
        wallZ = zTop + 30
        self.zOut = np.where(~bas, z, 0.0).astype(np.float32).ravel()
        fixedEnd = np.zeros((ny, nx), bool); fixedEnd[:, :2] = True; fixedEnd[:, nx - 2:] = True
        kind[:] = np.where(~bas, 1, np.where(chD, 0, 3)).astype(np.uint8)
        fl[:] = bas & chD
        zfix[:] = np.where(~bas, wallZ, np.where(chD, np.where(fixedEnd, z, z - P['erodDepth']), z - soilD))
        z[:] = np.where(~bas, wallZ, z)
        Bx[:] = np.maximum(1, fl.sum(axis=0)) * dx
        self.basin = bas.ravel().astype(np.uint8)
        self.demMode = True; self.zTop = zTop; self.wallZ = wallZ
        Io = T.get('Iout') if (T.get('Iout') or 0) > 0 else 0.01
        self.I = Io; self.Ix[:] = Io
        zb = np.where(bas, z, np.inf); cm = zb.min(axis=0)
        colmin[:] = np.where(np.isfinite(cm), cm, wallZ)
        return Io

    def _cell_at(self, x, y):
        i, j = int(math.floor(x / self.dx)), int(math.floor(y / self.dx))
        if i < 0 or j < 0 or i >= self.nx or j >= self.ny:
            return -1
        return j * self.nx + i

    def _snap(self, c):
        kind = self.kind.ravel()
        if c < 0 or kind[c] != 1:
            return c
        nx, ny = self.nx, self.ny; ci, cj = c % nx, c // nx; best, bd = -1, 1e9
        for dj in range(-3, 4):
            for di in range(-3, 4):
                ii, jj = ci + di, cj + dj
                if ii < 0 or jj < 0 or ii >= nx or jj >= ny:
                    continue
                o = jj * nx + ii
                if kind[o] == 1 or kind[o] == 2:
                    continue
                dd = di * di + dj * dj
                if dd < bd:
                    bd, best = dd, o
        return best if best >= 0 else c

    def _setup_points(self):
        """点の流入：谷頭などの供給点（それぞれのハイドロ）と、雨の横流入（谷筋のセルごとの流量の表）"""
        P = self.P; N = self.nx * self.ny
        self.pin = []
        for S in P.get('sources') or []:
            c = self._snap(self._cell_at(S['x'], S['y']))
            if c < 0:
                continue
            ux = S.get('ux', 1.0); uy = S.get('uy', 0.0); L = math.sqrt(ux * ux + uy * uy) or 1.0
            self.pin.append(dict(c=c, ux=ux / L, uy=uy / L, I=max(1e-3, S.get('I') or 0.05), H=S.get('hydro'), tab=None, Q=0.0, Qs=0.0, v=0.0))
        LI = P.get('latIn')
        if LI and LI.get('pts') and LI.get('Q'):
            for k, pt in enumerate(LI['pts']):
                c = self._snap(self._cell_at(pt[0], pt[1]))
                if c < 0:
                    continue
                self.pin.append(dict(c=c, ux=pt[2], uy=pt[3], I=max(1e-3, pt[4] or 0.05), H=None, tab=np.asarray(LI['Q'][k], dtype=float),
                                     dtab=LI['dt'] if (LI.get('dt') or 0) > 0 else 60.0, Q=0.0, Qs=0.0, v=0.0, lat=True))
        self.pinCells = np.array(sorted(set(p['c'] for p in self.pin)), dtype=np.int64)
        self.pH = np.zeros(N); self.pMx = np.zeros(N); self.pMy = np.zeros(N)
        R = P.get('rain') or {}
        self.rainOn = bool(R.get('on') and R.get('method') == 'grid' and getattr(self, 'demMode', False))
        if self.rainOn:
            self.rS = np.zeros(N); self.rR = np.zeros(N)

    def _pin_step(self):
        P = self.P; t = self.t; dx = self.dx; dbIn = P.get('transport') == 'debris'; Cmax = 0.9 * (1 - P['porosity']); nman = P['manning']
        tot = 0.0
        for p in self.pin:
            c0 = p['c']; self.pH[c0] = 0.0; self.pMx[c0] = 0.0; self.pMy[c0] = 0.0
        for p in self.pin:
            if p['tab'] is not None:
                f = t / p['dtab']; k0 = int(math.floor(f)); w = f - k0; nT = len(p['tab'])
                Q = p['tab'][nT - 1] if k0 >= nT - 1 else p['tab'][k0] * (1 - w) + p['tab'][k0 + 1] * w; Qs = 0.0
            else:
                Q = QatH(p['H'], t); Qs = QsAtH(p['H'], t)
            if not Q > 0:
                p['Q'] = 0.0; p['Qs'] = 0.0; p['v'] = 0.0
                continue
            if dbIn and Qs > Cmax * Q:
                Qs = Cmax * Q
            p['Q'] = Q; p['Qs'] = Qs if Qs > 0 else 0.0; tot += Q
            q = Q / dx
            p['v'] = q / self.debris_depth(q, p['Qs'] / Q, p['I']) if dbIn else q / (q * nman / math.sqrt(p['I'])) ** 0.6
            a = Q / (dx * dx)
            self.pH[p['c']] += a; self.pMx[p['c']] += a * p['v'] * p['ux']; self.pMy[p['c']] += a * p['v'] * p['uy']
        return tot

    def _pin_feed(self, dtm, bulk):
        if not self.pin:
            return
        K, dx, lam = self.K, self.dx, 1 - self.P['porosity']; add = 0.0
        zf = self.z.ravel(); dzf = self.dzs.ravel()
        for p in self.pin:
            if not p['Qs'] > 0:
                continue
            per = p['Qs'] * dtm / (dx * dx); c = p['c']; add += p['Qs'] * dtm
            if bulk:
                self.A[c] += per / lam * self.ff; zf[c] += per / lam; dzf[c] += per / lam; self._fix_layers(np.array([c]))
            else:
                self.Cc[c] += per * self.ff
        self.sedIn += add; self.QsIn += add / max(1e-9, dtm)

    def rain_at_cell(self, c, t):
        R = self.P.get('rain') or {}
        if R.get('mode') == 'grid' and R.get('frames'):
            fr = _rain_frame(R['frames'], t)
            return max(0.0, fr['r'][self.tIdx[c]] or 0.0) if (fr is not None and getattr(self, 'tIdx', None) is not None) else 0.0
        return _rain_uni(R.get('table'), t)

    def _rain_step(self, dt):
        P = self.P; R = P['rain']; N = self.nx * self.ny; bas = self.basin; t = self.t; a = self.dx * self.dx
        L = R.get('loss') or {'type': 'f', 'f': 1.0}; gwOn = bool((P.get('gw') or {}).get('on'))
        fr = _rain_frame(R['frames'], t) if (R.get('mode') == 'grid' and R.get('frames')) else None
        ru = 0.0 if fr is not None else _rain_uni(R.get('table'), t)
        rr = np.zeros(N); vol = 0.0; tot = 0.0
        frr = np.asarray(fr['r'], dtype=float) if fr is not None else None
        for c in range(N):
            if not bas[c]:
                continue
            r = max(0.0, frr[self.tIdx[c]]) if frr is not None else ru
            if not r > 0:
                continue
            e = r if gwOn else _loss_cell(L, self.rS, self.rR, c, r, dt)
            rr[c] = e / 3.6e6 * dt; vol += rr[c]; tot += r / 3.6e6 * dt
        self._rainV = vol * a; self.rainVol += tot * a; self.rainLoss += (tot - vol) * a
        return rr.reshape(self.ny, self.nx)

    # ======================================================================================
    # 斜面崩壊・天然ダム（ブラウザ版と同じ順で 1 セルずつ処理する）
    # ======================================================================================
    def _setup_slides(self):
        P = self.P; nx, ny, dx = self.nx, self.ny, self.dx; N = nx * ny
        z = self.z.ravel(); kind = self.kind.ravel(); bas = getattr(self, 'basin', None)
        self.slides = []; self.dep = np.zeros(N, dtype=np.uint8); self.depOn = False; self.ldamList = []; self._chkT = 0.0
        for q, S in enumerate(P.get('slides') or []):
            R = max(dx, S.get('R') or 20.0); D = max(0.1, S.get('D') or 2.0); deep = S.get('type') == 'deep'
            cells = []; zs = []
            i0 = max(0, int(math.floor((S['x'] - R) / dx))); i1 = min(nx - 1, int(math.floor((S['x'] + R) / dx)))
            j0 = max(0, int(math.floor((S['y'] - R) / dx))); j1 = min(ny - 1, int(math.floor((S['y'] + R) / dx)))
            zt = -1e30; sx = 0; sy = 0
            for j in range(j0, j1 + 1):
                for i in range(i0, i1 + 1):
                    c = j * nx + i
                    if kind[c] == 2 or (bas is not None and not bas[c]):
                        continue
                    rx = (i + 0.5) * dx - S['x']; ry = (j + 0.5) * dx - S['y']; rr = math.sqrt(rx * rx + ry * ry) / R
                    if rr > 1:
                        continue
                    th = D * math.sqrt(max(0.0, 1 - rr * rr)) if deep else D
                    if th < 0.05:
                        continue
                    cells.append(c); zs.append(z[c] - th)
                    if z[c] > zt:
                        zt = z[c]
                    sx += i; sy += j
            if not cells:
                continue
            ne = S['ne'] if (S.get('ne') or 0) > 0 else (0.05 if deep else 0.15)
            m0 = S['m0'] if (S.get('m0') is not None and S['m0'] >= 0) else 0.2
            self.slides.append(dict(id=q, S=S, cells=cells, zslip=np.array(zs), deep=deep, D=D, R=R, ne=ne, fired=False, done=np.zeros(len(cells), np.uint8),
                                    Sw=m0 * D * ne * 1000, Fs=float('nan'), zTop=zt, ci=_jround(sx / len(cells)), cj=_jround(sy / len(cells)), m=m0, fsc=None))
        for L in P.get('ldams') or []:
            self._place_ldam(L)

    def _surf_slope(self, c):
        nx, ny, dx = self.nx, self.ny, self.dx; z = self.z.ravel(); i = c % nx; j = c // nx
        a = c - 1 if i > 0 else c; b = c + 1 if i < nx - 1 else c; s = c - nx if j > 0 else c; n = c + nx if j < ny - 1 else c
        gx = 0.0 if b == a else (z[b] - z[a]) / ((b - a) * dx)
        gy = 0.0 if n == s else (z[n] - z[s]) / (((n - s) // nx) * dx)
        return math.sqrt(gx * gx + gy * gy)

    def _slide_hw(self, sl, k):
        P = self.P; c = sl['cells'][k]; z = self.z.ravel(); T = z[c] - sl['zslip'][k]
        if not sl['deep'] and (P.get('gw') or {}).get('on') and self.por.ravel()[c]:
            n = P['porosity']; b = self.W.ravel()[c] / n; base = self.zfix.ravel()[c]
            return max(0.0, min(T, base + b - sl['zslip'][k]))
        return max(0.0, min(T, sl['Sw'] / 1000 / sl['ne']))

    def _slide_fs(self, sl):
        S = sl['S']; gt = S['gt'] if (S.get('gt') or 0) > 0 else 18.0; gs = S['gs'] if (S.get('gs') or 0) > 0 else 20.0; gw = 9.81
        coh = S['c'] if (S.get('c') is not None and S['c'] >= 0) else (30.0 if sl['deep'] else 5.0)
        tp = math.tan((S['phi'] if (S.get('phi') or 0) > 0 else (30.0 if sl['deep'] else 35.0)) * math.pi / 180)
        z = self.z.ravel(); dx = self.dx; nx, ny = self.nx, self.ny; n = len(sl['cells'])
        if not sl['deep']:
            if sl['fsc'] is None:
                sl['fsc'] = np.zeros(n)
            fsc = sl['fsc']; mn = 1e9; msum = 0.0; cnt = 0
            for k in range(n):
                if sl['done'][k]:
                    fsc[k] = 9.0; continue
                c = sl['cells'][k]; D = z[c] - sl['zslip'][k]
                if D <= 0.01:
                    fsc[k] = 9.0; continue
                hw = self._slide_hw(sl, k); tb = self._surf_slope(c); b = math.atan(tb); cb = math.cos(b); sb = math.sin(b)
                res = coh + (gt * (D - hw) + (gs - gw) * hw) * cb * cb * tp; drv = (gt * (D - hw) + gs * hw) * sb * cb
                fsc[k] = res / drv if drv > 1e-9 else 9.0
                if fsc[k] < mn:
                    mn = fsc[k]
                msum += hw / D; cnt += 1
            sl['m'] = msum / cnt if cnt else 0.0
            return mn
        gxs = 0.0; gys = 0.0
        for k in range(n):
            c = sl['cells'][k]; i = c % nx; j = c // nx
            if 0 < i < nx - 1:
                gxs += (z[c + 1] - z[c - 1]) / (2 * dx)
            if 0 < j < ny - 1:
                gys += (z[c + nx] - z[c - nx]) / (2 * dx)
        L = math.sqrt(gxs * gxs + gys * gys)
        if L < 1e-9:
            return 9.0
        ux = -gxs / L; uy = -gys / L; num = 0.0; den = 0.0; a1 = dx * dx; ms = 0.0; Nn = nx * ny
        zsl = {}
        for k in range(n):
            zsl[sl['cells'][k]] = sl['zslip'][k]
        for k in range(n):
            c = sl['cells'][k]; T = z[c] - sl['zslip'][k]
            if T <= 0.01:
                continue
            i2 = c % nx; j2 = c // nx
            sxp = zsl[c + 1] if (i2 < nx - 1 and (c + 1) in zsl) else z[min(c + 1, Nn - 1)]
            sxm = zsl[c - 1] if (i2 > 0 and (c - 1) in zsl) else z[max(0, c - 1)]
            syp = zsl[c + nx] if (j2 < ny - 1 and (c + nx) in zsl) else (z[c + nx] if j2 < ny - 1 else zsl[c])
            sym = zsl[c - nx] if (j2 > 0 and (c - nx) in zsl) else (z[c - nx] if j2 > 0 else zsl[c])
            al = math.atan(-((sxp - sxm) / (2 * dx) * ux + (syp - sym) / (2 * dx) * uy))
            hw2 = self._slide_hw(sl, k); Wt = (gt * (T - hw2) + gs * hw2) * a1; Ab = a1 / max(0.2, math.cos(al)); u = gw * hw2 * math.cos(al) * math.cos(al)
            num += coh * Ab + max(0.0, Wt * math.cos(al) - u * Ab) * tp; den += Wt * math.sin(al); ms += hw2 / T
        sl['m'] = ms / n if n else 0.0; sl['ux'] = ux; sl['uy'] = uy
        return num / den if den > 1e-9 else 9.0

    def _slide_step(self, dt):
        P = self.P; t = self.t
        for sl in self.slides:
            S = sl['S']
            if sl['fired']:
                continue
            Tk = (S['Tk'] if (S.get('Tk') or 0) > 0 else (48.0 if sl['deep'] else 6.0)) * 3600
            r = self.rain_at_cell(sl['cells'][0], t) if (P.get('rain') or {}).get('on') else 0.0
            sl['Sw'] += (r / 3600 - sl['Sw'] / Tk) * dt
            if sl['Sw'] > sl['D'] * sl['ne'] * 1000:
                sl['Sw'] = sl['D'] * sl['ne'] * 1000
        self._chkT += dt
        chk = P['slideChk'] if (P.get('slideChk') or 0) > 0 else 1.0
        if self._chkT < chk:
            return
        self._chkT = 0.0
        for sl in self.slides:
            S2 = sl['S']
            if sl['fired']:
                continue
            Fs = self._slide_fs(sl); sl['Fs'] = Fs
            if S2.get('trigger') == 'time':
                if t >= (S2.get('t') or 0):
                    self._release_slide(sl, None)
            elif sl['deep']:
                if Fs < 1:
                    self._release_slide(sl, None)
            elif Fs < 1:
                lst = [k for k in range(len(sl['cells'])) if not sl['done'][k] and sl['fsc'][k] < 1]
                if lst:
                    self._release_slide(sl, lst)
        if self.depOn and (P.get('gw') or {}).get('on'):
            self._progressive()
        self._ldam_slip()

    def _take_sed(self, c, dz, out):
        K = self.K; left = dz; it = 0; lam = 1 - self.P['porosity']
        while left > 1e-9 and it < 400:
            it += 1
            av = float(self.A[c].sum())
            if av < 1e-9:
                break
            take = left if left < av else av; frc = take / av
            tt = self.A[c] * frc; self.A[c] -= tt
            if out is not None:
                out += tt * lam
            left -= take; self._fix_layers(np.array([c]))
        return dz - left

    def _add_deposit(self, c, dz, fr):
        kind = self.kind.ravel(); z = self.z.ravel(); zfix = self.zfix.ravel(); fl = self.fl.ravel()
        if kind[c] != 0:
            kind[c] = 0; zfix[c] = z[c]; fl[c] = True
            if getattr(self, 'fixedC', None) is not None:
                self.fixedC.ravel()[c] = False
        self.A[c] += dz * fr
        z[c] += dz; self.dzs.ravel()[c] += dz; self._fix_layers(np.array([c]))
        self.por.ravel()[c] = True; self.dep[c] = 1; self.depOn = True

    def _release_slide(self, sl, lst):
        P = self.P; S = sl['S']; n = P['porosity']; Cst = 1 - n; K = self.K; dx = self.dx; a1 = dx * dx
        z = self.z.ravel(); h = self.h.ravel(); zfix = self.zfix.ravel(); kind = self.kind.ravel(); W = self.W.ravel(); por = self.por.ravel()
        gwOn = bool((P.get('gw') or {}).get('on'))
        idx = lst if lst is not None else list(range(len(sl['cells'])))
        self._slide_fs(sl)
        mode = S.get('move') or 'auto'
        if mode == 'auto':
            mode = 'fluid' if (sl['m'] >= (S['mc'] if (S.get('mc') or 0) > 0 else 0.8) and P.get('transport') == 'debris') else 'block'
        if mode == 'fluid' and P.get('transport') != 'debris':
            mode = 'block'
        CL = S['CL'] if (S.get('CL') or 0) > 0 else 0.5; Sr0 = S['Sr0'] if (S.get('Sr0') is not None and S['Sr0'] >= 0) else 0.5
        vol = 0.0; wat = 0.0; sed = np.zeros(K); added = 0.0; sedTot = 0.0
        for k in idx:
            if sl['done'][k]:
                continue
            sl['done'][k] = 1
            c = sl['cells'][k]; T = z[c] - sl['zslip'][k]
            if T <= 1e-6:
                continue
            hw = self._slide_hw(sl, k)
            if gwOn and por[c] and not sl['deep']:
                wsl = max(0.0, min(W[c], W[c] - n * max(0.0, sl['zslip'][k] - zfix[c]))); W[c] -= wsl
            elif gwOn:
                wsl = n * (hw + Sr0 * max(0.0, T - hw)); self.Vwin += wsl * a1
            else:
                wsl = n * T if mode == 'fluid' else 0.0
            frc = np.zeros(K); got = 0.0
            if kind[c] == 0:
                got = self._take_sed(c, min(T, max(0.0, z[c] - zfix[c])), frc)
            frc += (T - got) * Cst * self.f0
            z[c] = sl['zslip'][k]; self.dzs.ravel()[c] -= T
            if zfix[c] > z[c]:
                zfix[c] = z[c]
            if kind[c] == 3:
                zfix[c] = z[c]
            if gwOn and por[c]:
                cap = n * max(0.0, z[c] - zfix[c])
                if W[c] > cap:
                    h[c] += W[c] - cap; W[c] = cap
            if mode == 'fluid':
                ss = float(frc.sum()); self.Cc[c] += frc
                h[c] += ss + wsl
                Sx = float(self.Cc[c].sum())
                if Sx > CL * h[c]:
                    ad = Sx / CL - h[c]; h[c] += ad; added += ad * a1
                sedTot += ss * a1
            else:
                sed += frc * a1; vol += T * a1; wat += wsl * a1
        self.Vwin += added
        if mode == 'block' and vol > 0:
            ss2 = float(sed.sum())
            frc2 = sed / ss2 if ss2 > 0 else self.f0.copy()
            self._runout(sl, ss2 / Cst, wat, frc2)
        if lst is None or sl['done'].all():
            sl['fired'] = True
        self.slideLog.append(dict(t=self.t, id=sl['id'], mode=mode, vol=vol if mode == 'block' else sedTot / Cst, Fs=sl['Fs'], m=sl['m'], water=added))

    def _runout(self, sl, V, Wv, frc):
        S = sl['S']; nx, ny, dx = self.nx, self.ny, self.dx; z = self.z.ravel(); kind = self.kind.ravel(); bas = getattr(self, 'basin', None)
        if (S.get('mu') or 0) > 0:
            mu = math.tan(S['mu'] * math.pi / 180)
        else:
            mu = max(0.05, min(0.9, 10 ** (-0.15666 * math.log(max(1.0, V)) / math.log(10) + 0.62419)))
        c = sl['cj'] * nx + sl['ci']; L = 0.0; it = 0; top = sl['zTop']
        DCX = [1, 1, 0, -1, -1, -1, 0, 1]; DCY = [0, 1, 1, 1, 0, -1, -1, -1]
        while it < 20000:
            it += 1
            ci = c % nx; cj = c // nx; best = -1; bs = 0.0; bl = 0.0
            for q in range(8):
                ii = ci + DCX[q]; jj = cj + DCY[q]
                if ii < 0 or jj < 0 or ii >= nx or jj >= ny:
                    continue
                o = jj * nx + ii
                if kind[o] == 2 or (bas is not None and not bas[o]):
                    continue
                dl = (math.sqrt(2) if q % 2 else 1.0) * dx; s = (z[c] - z[o]) / dl
                if s > bs:
                    bs, best, bl = s, o, dl
            if best < 0:
                break
            L += bl; c = best
            if z[c] >= top - mu * (L + sl['R']):
                break
        sl['stop'] = c; sl['mu'] = mu; sl['L'] = L
        self._deposit_cone(c, V, math.tan((S['thDep'] if (S.get('thDep') or 0) > 0 else 10.0) * math.pi / 180), Wv, frc)

    def _deposit_cone(self, c0, V, tanT, Wv, frc):
        nx, ny, dx = self.nx, self.ny, self.dx; z = self.z.ravel(); kind = self.kind.ravel(); a1 = dx * dx; bas = getattr(self, 'basin', None)
        ci = c0 % nx; cj = c0 // nx; z0 = z[c0]

        def vol(H, apply=False):
            R = H / tanT; rr = int(math.ceil(R / dx)) + 1; s = 0.0; cells = [] if apply else None
            for j in range(max(0, cj - rr), min(ny - 1, cj + rr) + 1):
                for i in range(max(0, ci - rr), min(nx - 1, ci + rr) + 1):
                    c = j * nx + i
                    if kind[c] == 2 or (bas is not None and not bas[c]):
                        continue
                    r = math.sqrt((i - ci) * (i - ci) + (j - cj) * (j - cj)) * dx; zt = z0 + H - tanT * r
                    if zt > z[c]:
                        s += (zt - z[c]) * a1
                        if apply:
                            cells.append((c, zt - z[c]))
            return cells if apply else s
        lo = 0.0; hi = max(1.0, V ** (1.0 / 3.0)); it = 0
        while it < 60 and vol(hi) < V:
            hi *= 2; it += 1
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if vol(mid) < V:
                lo = mid
            else:
                hi = mid
        cells = vol(hi, True); tot = sum(cq[1] for cq in cells)
        self._deposit_cells(cells, V / (tot * a1) if tot > 0 else 0.0, Wv, frc)

    def _deposit_cells(self, cells, sc, Wv, frc):
        P = self.P; n = P['porosity']; a1 = self.dx * self.dx; gwOn = bool((P.get('gw') or {}).get('on'))
        z = self.z.ravel(); zfix = self.zfix.ravel(); W = self.W.ravel(); h = self.h.ravel()
        fb = np.array(frc, dtype=float)
        tot = 0.0
        for cq in cells:
            tot += cq[1] * sc
        for c, th in cells:
            t = th * sc
            if not t > 0:
                continue
            self._add_deposit(c, t, fb)
            w = Wv / a1 * t / tot if tot > 0 else 0.0
            if gwOn:
                cap = n * max(0.0, z[c] - zfix[c]); room = max(0.0, cap - W[c]); ww = min(w, room); W[c] += ww; h[c] += w - ww
            else:
                h[c] += w
            self.Vdep = getattr(self, 'Vdep', 0.0) + t * a1

    def _place_ldam(self, L):
        nx, ny, dx = self.nx, self.ny, self.dx; z = self.z.ravel(); kind = self.kind.ravel(); bas = getattr(self, 'basin', None)
        c0 = self._snap(self._cell_at(L['x'], L['y']))
        if c0 < 0:
            return
        ux = L.get('ux', 1.0); uy = L.get('uy', 0.0); nn = math.sqrt(ux * ux + uy * uy) or 1.0; ux /= nn; uy /= nn
        H = max(0.5, L.get('H') or 10.0); Wc = max(0.0, L.get('Wc') or 0.0)
        tu = math.tan((L['up'] if (L.get('up') or 0) > 0 else 30.0) * math.pi / 180); td = math.tan((L['dn'] if (L.get('dn') or 0) > 0 else 25.0) * math.pi / 180)
        lat = L['lat'] if (L.get('lat') or 0) > 0 else max(50.0, 8 * H)
        x0 = (c0 % nx + 0.5) * dx; y0 = (c0 // nx + 0.5) * dx; zc = z[c0]
        ext = max(H / tu, H / td) + Wc + lat; cells = []
        for j in range(max(0, int(math.floor((y0 - ext) / dx))), min(ny - 1, int(math.floor((y0 + ext) / dx))) + 1):
            for i in range(max(0, int(math.floor((x0 - ext) / dx))), min(nx - 1, int(math.floor((x0 + ext) / dx))) + 1):
                c = j * nx + i
                if kind[c] == 2 or (bas is not None and not bas[c]):
                    continue
                px = (i + 0.5) * dx - x0; py = (j + 0.5) * dx - y0; s = px * ux + py * uy; w = -px * uy + py * ux
                if abs(w) > lat:
                    continue
                zt = zc + H - ((s - Wc / 2) * td if s > Wc / 2 else ((-Wc / 2 - s) * tu if s < -Wc / 2 else 0.0))
                if zt > z[c] + 1e-3:
                    cells.append((c, zt - z[c]))
        LD = self.P.get('ldam') or {}; K = self.K; fr = self.f0.copy()
        gr = LD.get('grains')
        if gr and len(gr) == K:
            sm = sum(max(0.0, v) for v in gr)
            fr = np.array([max(0.0, v) / sm for v in gr]) if sm > 0 else self.f0.copy()
        n = self.P['porosity']; gwOn = bool((self.P.get('gw') or {}).get('on')); s0 = LD['s0'] if (LD.get('s0') is not None and LD['s0'] >= 0) else 0.3
        V = 0.0; W = self.W.ravel()
        for c, t in cells:
            self._add_deposit(c, t, fr); V += t * dx * dx
            if gwOn:
                W[c] += n * t * s0
        self.ldamList.append(dict(c0=c0, x0=x0, y0=y0, zc=zc, H=H, ux=ux, uy=uy, crest=zc + H, lat=lat, V=V, cells=[c for c, _ in cells], slipped=False, L=L))

    def _progressive(self):
        """進行性崩壊：湧水の出ている（地表まで飽和した）天然ダムの法面で、浸透流が斜面に平行な無限長斜面の安全率が 1 を下回れば、
        下の隣のセルへ限界勾配 tanβc = (γ'/γsat) tanφ になるまで土を移す"""
        P = self.P; LD = P.get('ldam') or {}; nx, ny, dx = self.nx, self.ny, self.dx; N = nx * ny; n = P['porosity']
        z = self.z.ravel(); zfix = self.zfix.ravel(); h = self.h.ravel(); W = self.W.ravel(); kind = self.kind.ravel(); por = self.por.ravel(); fl = self.fl.ravel()
        bas = getattr(self, 'basin', None)
        gs = LD['gs'] if (LD.get('gs') or 0) > 0 else 20.0; gw = 9.81; tp = math.tan((LD['phi'] if (LD.get('phi') or 0) > 0 else 35.0) * math.pi / 180)
        coh = LD['c'] if (LD.get('c') or 0) > 0 else 0.0; Ts = LD['Ts'] if (LD.get('Ts') or 0) > 0 else 0.5
        tanC = (gs - gw) / gs * tp; dep = self.dep; cnt = 0; vol = 0.0
        for c in np.nonzero(dep)[0]:
            c = int(c)
            if not dep[c]:
                continue
            Tq = z[c] - zfix[c]
            if Tq < 1e-3:
                continue
            if W[c] < 0.98 * n * Tq or h[c] > 0.1:
                continue
            i = c % nx; j = c // nx; best = -1; bs = 0.0
            for o in (c + 1 if i < nx - 1 else -1, c - 1 if i > 0 else -1, c + nx if j < ny - 1 else -1, c - nx if j > 0 else -1):
                if o < 0 or kind[o] == 2 or (bas is not None and not bas[o]):
                    continue
                s = (z[c] - z[o]) / dx
                if s > bs:
                    bs, best = s, o
            if best < 0:
                continue
            b = math.atan(bs); Fs = coh / (gs * Ts * math.sin(b) * math.cos(b)) + (gs - gw) / gs * tp / bs
            if Fs >= 1:
                continue
            dl = min((z[c] - z[best] - tanC * dx) / 2, Tq)
            if not dl > 1e-4:
                continue
            got = self._move_sed(c, best, dl)
            if not got > 0:
                continue
            self.dzs.ravel()[c] -= got; self.dzs.ravel()[best] += got
            if kind[best] != 0:
                kind[best] = 0; zfix[best] = z[best] - got; fl[best] = True
            por[best] = True; dep[best] = 1
            wm = min(W[c], n * got); W[c] -= wm
            cap = n * max(0.0, z[best] - zfix[best]); room = max(0.0, cap - W[best]); ww = min(wm, room); W[best] += ww; h[best] += wm - ww
            cnt += 1; vol += got * dx * dx
        if cnt:
            self.progCount = getattr(self, 'progCount', 0) + cnt; self.progVol = getattr(self, 'progVol', 0.0) + vol

    def _ldam_slip(self):
        LD = self.P.get('ldam') or {}
        if LD.get('fail') != 'slip' or not self.ldamList:
            return
        nx, ny, dx = self.nx, self.ny, self.dx; z = self.z.ravel(); h = self.h.ravel(); zfix = self.zfix.ravel(); a1 = dx * dx; K = self.K
        for m, D in enumerate(self.ldamList):
            if D['slipped']:
                continue
            wl = -1e30; sTo = -1e30
            for c in D['cells']:
                px = (c % nx + 0.5) * dx - D['x0']; py = (c // nx + 0.5) * dx - D['y0']; s = px * D['ux'] + py * D['uy']
                if z[c] - zfix[c] > 0.05 and s > sTo:
                    sTo = s
            r0 = int(math.ceil((D['H'] * 4 + D['lat']) / dx)); ci = D['c0'] % nx; cj = D['c0'] // nx
            for j in range(max(0, cj - r0), min(ny - 1, cj + r0) + 1):
                for i in range(max(0, ci - r0), min(nx - 1, ci + r0) + 1):
                    c = j * nx + i; px2 = (i + 0.5) * dx - D['x0']; py2 = (j + 0.5) * dx - D['y0']; s2 = px2 * D['ux'] + py2 * D['uy']
                    if s2 < 0 and h[c] > 0.05 and z[c] + h[c] > wl:
                        wl = z[c] + h[c]
            lev = D['zc'] + LD['slipLevel'] * D['H'] if (LD.get('slipLevel') or 0) > 0 else D['crest'] - 0.05 * D['H']
            if wl < lev:
                continue
            D['slipped'] = True
            al = math.tan((LD['slipDeg'] if (LD.get('slipDeg') or 0) > 0 else 20.0) * math.pi / 180)
            tanD = math.tan((LD['thDep'] if (LD.get('thDep') or 0) > 0 else 10.0) * math.pi / 180)
            sed = np.zeros(K); Vb = 0.0; Wm = 0.0; n = self.P['porosity']; ztoe = D['zc']; W = self.W.ravel()
            for c in D['cells']:
                px3 = (c % nx + 0.5) * dx - D['x0']; py3 = (c // nx + 0.5) * dx - D['y0']; s3 = px3 * D['ux'] + py3 * D['uy']
                zp = ztoe + (sTo - s3) * al; base = max(zfix[c], zp); ex = z[c] - base
                if ex <= 1e-3:
                    continue
                got = self._take_sed(c, ex, sed); z[c] -= got; self.dzs.ravel()[c] -= got; Vb += got * a1
                cap = n * max(0.0, z[c] - zfix[c])
                if W[c] > cap:
                    Wm += (W[c] - cap) * a1; W[c] = cap
            if Vb <= 0:
                continue
            ss = float(sed.sum()); fr = sed / ss if ss > 0 else self.f0.copy()
            self._deposit_ramp(D, sTo, Vb, tanD, Wm, fr)
            self.slideLog.append(dict(t=self.t, id='ldam%d' % m, mode='slip', vol=Vb, water=0.0))

    def _deposit_ramp(self, D, s0, V, tanT, Wv, frc):
        nx, ny, dx = self.nx, self.ny, self.dx; z = self.z.ravel(); kind = self.kind.ravel(); a1 = dx * dx; bas = getattr(self, 'basin', None)
        r0 = int(math.ceil((D['H'] * 12 + D['lat']) / dx)); ci = D['c0'] % nx; cj = D['c0'] // nx

        def vol(Z, apply=False):
            s = 0.0; cells = [] if apply else None
            for j in range(max(0, cj - r0), min(ny - 1, cj + r0) + 1):
                for i in range(max(0, ci - r0), min(nx - 1, ci + r0) + 1):
                    c = j * nx + i
                    if kind[c] == 2 or (bas is not None and not bas[c]):
                        continue
                    px = (i + 0.5) * dx - D['x0']; py = (j + 0.5) * dx - D['y0']; sv = px * D['ux'] + py * D['uy']; w = -px * D['uy'] + py * D['ux']
                    if sv < s0 - dx or abs(w) > D['lat']:
                        continue
                    zt = Z - tanT * max(0.0, sv - s0)
                    if zt > z[c]:
                        s += (zt - z[c]) * a1
                        if apply:
                            cells.append((c, zt - z[c]))
            return cells if apply else s
        lo = D['zc'] - 50; hi = D['crest'] + 50; it = 0
        while it < 60 and vol(hi) < V:
            hi += 50; it += 1
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if vol(mid) < V:
                lo = mid
            else:
                hi = mid
        cells = vol(hi, True); tot = sum(cq[1] for cq in cells)
        self._deposit_cells(cells, V / (tot * a1) if tot > 0 else 0.0, Wv, frc)

    # ------------------------------------------------------------ 伏流（河床への浸透・地下水の流れ・湧出）
    def balance(self):
        """土砂と水の体積収支の誤差 [m3]（土砂は実質体積。流れの中の土砂・伏流の水も数える）"""
        P = self.P; lam = 1 - P['porosity']; a2 = self.dx * self.dx; tr = P.get('transport')
        sw = self.sediment_in_water() if tr in ('noneq', 'debris') else 0.0
        eb = (self.sedIn - self.sedOut) - lam * self.storage() - sw
        hs = float(self.h.sum()) * a2
        if tr == 'debris' and (P.get('gw') or {}).get('on'):
            ew = (hs - sw + self.ground_water()) - (self.Vw0 + self.W0 + (self.Vwin - self.sedIn) - (self.Vwout - self.sedOut) - self.gwOut)
        elif tr == 'debris':
            ew = (hs + self.storage() - self.voidAir) - (self.Vw0 + self.Vwin - self.Vwout)
        else:
            ew = (hs + self.ground_water()) - (self.Vw0 + self.W0 + self.Vwin - self.Vwout - self.gwOut)
        return eb, ew

    def ground_water(self):
        return float(self.W.sum()) * self.dx * self.dx

    def _gw_step(self, dt):
        """浸透：水に浸かったセルで、堆積層の空隙が残っていれば単位動水勾配で k の速さ。地下水：飽和厚 b = W/n、水頭 zfix + b の差で
        Darcy 流（Dupuit 近似）。湧出：W が空隙量 n(z − zfix) を超えた分は地表へ。k は指定がなければ Kozeny–Carman 式"""
        P = self.P; GW = P.get('gw') or {}; nx, ny, dx, K = self.nx, self.ny, self.dx, self.K
        n = P['porosity']; Cst = 1 - n; HD = 5e-4; kMax = GW.get('kMax') or 0.05
        h = self.h; z = self.z; zfix = self.zfix; W = self.W
        mov = self.por
        A = self.A.reshape(ny, nx, K)
        DMk = (P.get('dem') or {}).get('ks') or 0.0; LDk = (P.get('ldam') or {}).get('k') or 0.0
        if (GW.get('k') or 0) > 0:
            kk = np.full((ny, nx), float(GW['k']))
        else:
            kc0 = G / 1e-6 * n * n * n / ((1 - n) * (1 - n)) / 180
            sA = np.zeros((ny, nx)); sh = np.zeros((ny, nx))
            for k in range(K):
                sA = sA + A[:, :, k]; sh = sh + A[:, :, k] / self.d[k]
            dh = np.where(sh > 0, sA / np.where(sh > 0, sh, 1.0), self.dm0)
            kk = kc0 * dh * dh
        if LDk > 0:
            kk = np.where(self.dep.reshape(ny, nx) > 0, LDk, kk)
        kk = np.where(self.kind == 3, DMk if DMk > 0 else 1e-5, kk)
        kg = np.where(mov, np.minimum(kk, kMax), 0.0)
        self.kgw = kg
        Wm = n * np.maximum(0.0, z - zfix)
        ex = mov & (W > Wm)
        h = np.where(ex, h + (W - Wm), h); W = np.where(ex, Wm, W)
        hh = h
        av = hh - HD
        if P.get('transport') == 'debris':
            S = self.Cc.reshape(ny, nx, K).sum(axis=2)
            av = np.minimum(av, hh - S / (0.98 * Cst))
        f = np.minimum(np.minimum(kg * dt, Wm - W), av)
        inf = mov & (hh > HD) & (W < Wm) & (f > 0)
        f = np.where(inf, f, 0.0)
        r = np.where(inf, (hh - f) / np.where(hh > 0, hh, 1.0), 1.0)
        self.qx = np.where(inf, self.qx * r, self.qx); self.qy = np.where(inf, self.qy * r, self.qy)
        h = np.where(inf, hh - f, hh); W = W + f
        b = np.where(mov, W / n, 0.0); bmax = float(b.max()) if b.size else 0.0; kmx = float(kg.max()) if kg.size else 0.0
        dtg = 0.2 * n * dx * dx / (kmx * bmax) if kmx * bmax > 0 else dt
        ns = min(50, max(1, math.ceil(dt / dtg))); dts = dt / ns
        Ib = self.Ix[nx - 1]
        for _ in range(ns):
            Hg = zfix + W / n
            dW = np.zeros((ny, nx))
            # x 方向の面（c と c+1）
            okx = mov[:, :-1] & mov[:, 1:]
            dH = Hg[:, :-1] - Hg[:, 1:]
            upW = np.where(dH > 0, W[:, :-1], W[:, 1:])
            ka, kb = kg[:, :-1], kg[:, 1:]
            kf = np.where((ka > 0) & (kb > 0), 2 * ka * kb / np.where(ka + kb > 0, ka + kb, 1.0), 0.0)
            volx = kf * (upW / n) * dH / dx * dts / dx
            capx = 0.25 * upW; volx = np.clip(volx, -capx, capx); volx = np.where(okx & (dH != 0), volx, 0.0)
            oky = mov[:-1] & mov[1:]
            dHy = Hg[:-1] - Hg[1:]
            upWy = np.where(dHy > 0, W[:-1], W[1:])
            ka, kb = kg[:-1], kg[1:]
            kfy = np.where((ka > 0) & (kb > 0), 2 * ka * kb / np.where(ka + kb > 0, ka + kb, 1.0), 0.0)
            voly = kfy * (upWy / n) * dHy / dx * dts / dx
            capy = 0.25 * upWy; voly = np.clip(voly, -capy, capy); voly = np.where(oky & (dHy != 0), voly, 0.0)
            # ブラウザ版と同じ順に足す（南の面、西の面、東の面、北の面、下流端）
            dW[1:, :] = dW[1:, :] + voly
            dW[:, 1:] = dW[:, 1:] + volx
            dW[:, :-1] = dW[:, :-1] - volx
            dW[:-1, :] = dW[:-1, :] - voly
            vo = np.where(mov[:, -1], kg[:, -1] * (W[:, -1] / n) * Ib * dts / dx, 0.0)
            vo = np.minimum(vo, 0.25 * W[:, -1])
            dW[:, -1] = dW[:, -1] - vo
            for j_ in range(ny):
                if mov[j_, -1]:
                    self.gwOut += vo[j_] * dx * dx
            W = W + dW
        Wm2 = n * np.maximum(0.0, z - zfix)
        ex2 = mov & (W > Wm2)
        h = np.where(ex2, h + (W - Wm2), h); W = np.where(ex2, Wm2, W)
        self.h = h; self.W = W

    # ------------------------------------------------------------ 土石流（高橋 2007。Kanako 2D と同じ形の式）
    def cinf_t(self, t, h, d):
        """平衡濃度 C∞（tanθ は水面勾配）。配列でも可"""
        P = self.P; DB = P.get('debris') or {}; Cst = 1 - P['porosity']; tp = DB.get('tanPhi') or 0.7
        rho = DB.get('rho') or 1000.0; sig = 1000 * P['sg']
        t = np.asarray(t, dtype=float); h = np.asarray(h, dtype=float); d = np.asarray(d, dtype=float)
        tt = np.where(t > 0, np.minimum(t, tp - 1e-9), 1e-9)
        c1 = rho * tt / ((sig - rho) * (tp - tt))
        sgr = sig / rho
        ts = h * tt / ((sgr - 1) * np.where(d > 0, d, 1.0)); tc = 0.04 * 10 ** (1.72 * tt)
        st = sig * tt / (sig - rho); a2 = 2 * (0.425 - st) / (1 - st)
        tsS = np.where(ts > 0, ts, 1.0); a2S = np.where(a2 > 0, a2, 0.0)
        Cb = (1 + 5 * tt) * tt / (sgr - 1) * (1 - a2S * tc / tsS) * (1 - np.sqrt(np.maximum(0.0, a2S * tc / tsS)))
        Cb = np.where((ts > tc) & (a2 > 0), Cb, 0.0)
        C = np.where(tt > 0.138, c1, np.where(tt > 0.03, 6.7 * c1 * c1, Cb))
        C = np.where(C < 0, 0.0, np.where(C > 0.9 * Cst, 0.9 * Cst, C))
        C = np.where(t >= tp, 0.9 * Cst, C)
        return np.where(t > 0, C, 0.0)

    def _debris_params(self):
        """各セルの流れの土砂濃度と代表粒径（流れている土砂の平均粒径。なければ河床表層の値）"""
        nx, ny, K = self.nx, self.ny, self.K; Cst = 1 - self.P['porosity']
        Cc = self.Cc.reshape(ny, nx, K); S = np.zeros((ny, nx)); Sd = np.zeros((ny, nx))
        for k in range(K):
            S = S + Cc[:, :, k]; Sd = Sd + Cc[:, :, k] * self.d[k]
        h = self.h
        self.dfC = np.where(h > 1e-3, np.minimum(0.98 * Cst, S / np.where(h > 1e-3, h, 1.0)), 0.0)
        dm = self.surface_dm()
        self.dfD = np.where(S > 1e-9, Sd / np.where(S > 1e-9, S, 1.0), np.where(np.isfinite(dm), dm, self.dm0))

    def debris_depth(self, q, C, I):
        """流入する土石流の等流水深（単位幅流量 q、濃度 C、勾配 I）"""
        P = self.P; DB = P.get('debris') or {}; Cst = 1 - P['porosity']; rhoS = 1000 * P['sg']; rhoW = DB.get('rho') or 1000.0
        nman = P['manning']; sn = I / math.sqrt(1 + I * I); mfl = DB.get('manningFloor', True) is not False
        dd = 0.0
        for k in range(self.K):
            dd += self.ff[k] * self.d[k]
        Cd = min(C, 0.98 * Cst)

        def vel(hh):
            Vm = hh ** (2 / 3) * math.sqrt(max(I, 1e-4)) / nman
            if Cd < 0.01:
                return Vm
            if Cd >= 0.4 * Cst:
                X = (Cst / Cd) ** (1 / 3) - 1; psi = dd * dd / (8 * hh * hh * (Cd + (1 - Cd) * rhoW / rhoS) * X * X)
            else:
                psi = dd * dd / (0.49 * hh * hh)
            Vd = math.sqrt(G * hh * sn / psi)
            return min(Vd, Vm) if mfl else Vd
        lo, hi = 1e-4, 50.0
        for _ in range(60):
            mid = math.sqrt(lo * hi)
            if mid * vel(mid) > q:
                hi = mid
            else:
                lo = mid
        return math.sqrt(lo * hi)

    def _morph_debris(self, dt):
        """混合物の流束で土砂を運び（風上）、高橋の侵食・堆積速度式で河床と交換する"""
        P = self.P; DB = P.get('debris') or {}; nx, ny, K, dx = self.nx, self.ny, self.K, self.dx
        Cst = 1 - P['porosity']; HSED = 0.02
        de = DB.get('de') or 0.0007; dd = DB.get('dd') or 0.05; Vmin = DB.get('Vmin') if DB.get('Vmin') is not None else 0.1
        gwOn = bool((P.get('gw') or {}).get('on'))
        h = self.h; hO = self.hOld; FxM = self.FxM; GyM = self.GyM
        Cc = self.Cc.reshape(ny, nx, K).copy()
        # ---- 1. 移流（面の混合物流束 × 風上セルの濃度）
        SFx = np.zeros((ny, nx + 1, K)); SFy = np.zeros((ny + 1, nx, K))
        F = FxM[:, 1:nx]
        cuL = F > 0
        huF = np.where(cuL, hO[:, :-1], hO[:, 1:])
        CcU = np.where(cuL[..., None], Cc[:, :-1], Cc[:, 1:])
        SFx[:, 1:nx] = np.where((huF >= 1e-6)[..., None], F[..., None] * CcU / np.where(huF >= 1e-6, huF, 1.0)[..., None], 0.0)
        Fe = FxM[:, nx]; okE = (Fe > 0) & (hO[:, -1] >= 1e-6)
        SFx[:, nx] = np.where(okE[:, None], Fe[:, None] * Cc[:, -1] / np.where(okE, hO[:, -1], 1.0)[:, None], 0.0)
        Gf = GyM[1:ny]
        cuS = Gf > 0
        hu = np.where(cuS, hO[:-1], hO[1:])
        CcV = np.where(cuS[..., None], Cc[:-1], Cc[1:])
        SFy[1:ny] = np.where((hu >= 1e-6)[..., None], Gf[..., None] * CcV / np.where(hu >= 1e-6, hu, 1.0)[..., None], 0.0)
        cf = dt / dx
        if self.hasOpen:
            sx0 = SFx.sum(axis=2); sy0 = SFy.sum(axis=2)
            self._check_block(Cc=Cc); self._open_filter(SFx, SFy)
            self._hold_volume(sx0 - SFx.sum(axis=2), sy0 - SFy.sum(axis=2), cf)
        outv = cf * (np.maximum(SFx[:, 1:], 0) - np.minimum(SFx[:, :-1], 0) + np.maximum(SFy[1:], 0) - np.minimum(SFy[:-1], 0))
        scS = np.where(outv > Cc + 1e-15, np.where(outv > 0, Cc / np.where(outv > 0, outv, 1), 0.0), 1.0)
        if (scS < 1).any():
            r = np.ones_like(SFx)
            r[:, 1:] = np.where(SFx[:, 1:] > 0, scS, np.concatenate([scS[:, 1:], np.ones((ny, 1, K))], axis=1))
            SFx *= r
            ry = np.ones_like(SFy); ry[1:ny] = np.where(SFy[1:ny] > 0, scS[:-1], scS[1:]); SFy *= ry
        Qs = self.QsAt(self.t) if self.phase == 'run' else 0.0; Qw = self.Qin
        self.QsIn = Qs; self.sedIn += Qs * dt
        Cin = Qs / Qw if (Qs > 0 and Qw > 0) else 0.0
        out_k = SFx[:, nx].sum(axis=0)
        self.sedOutK += out_k * dx * dt
        Cn = Cc - cf * (SFx[:, 1:] - SFx[:, :-1] + SFy[1:] - SFy[:-1])
        if Cin > 0:
            Cn[:, 0, :] += np.where(self.qin > 0, cf * Cin * self.qin, 0.0)[:, None] * self.ff[None, :]
        Cc = np.maximum(Cn, 0.0)
        self.QsOut = out_k.sum() * dx; self.sedOut += out_k.sum() * dx * dt
        self.Cc = Cc.reshape(-1, K)
        self._trib_feed(dt, False)
        Cc = self.Cc.reshape(ny, nx, K)
        # ---- 2. 河床との交換
        z = self.z; zfix = self.zfix; A = self.A.reshape(ny, nx, K); W = self.W
        eta = z + h
        S = np.zeros((ny, nx))
        for k in range(K):
            S = S + Cc[:, :, k]
        dry = h < HSED
        # 乾いたら、流れの中の土砂はその場に堆積させる
        drp = dry & (S > 0)
        dz0 = np.where(drp, S / Cst, 0.0)
        Anew = A + np.where(drp[..., None], Cc / Cst, 0.0)
        Ccn = np.where(drp[..., None], 0.0, Cc)
        zn = z + dz0; dzs = dz0.copy()
        if gwOn:
            pw0 = np.minimum((1 - Cst) * dz0, np.maximum(0.0, h - S))
            Wn = W + np.where(drp, pw0, 0.0); hn = np.where(drp, np.maximum(0.0, h - S - pw0), h)
        else:
            Wn = W.copy(); hn = np.where(drp, h - dz0, h)
        act = ~dry & ~self.fixedC
        hs = np.where(act, h, 1.0)
        fw = FxM; uu = 0.5 * (fw[:, :-1] + fw[:, 1:]) / hs; vv = 0.5 * (GyM[:-1] + GyM[1:]) / hs; V = np.sqrt(uu * uu + vv * vv)
        wet = ~dry
        wl = np.zeros_like(wet); wl[:, 1:] = wet[:, :-1]; wr = np.zeros_like(wet); wr[:, :-1] = wet[:, 1:]
        ws = np.zeros_like(wet); ws[1:] = wet[:-1]; wn = np.zeros_like(wet); wn[:-1] = wet[1:]
        eL = np.zeros_like(eta); eL[:, 1:] = eta[:, :-1]; eR = np.zeros_like(eta); eR[:, :-1] = eta[:, 1:]
        eS = np.zeros_like(eta); eS[1:] = eta[:-1]; eN = np.zeros_like(eta); eN[:-1] = eta[1:]
        # 水面勾配は流れの下流側の片側差分（中心差分では 1 セルおきの凹凸が見えず、堆積面が市松模様になる）
        gxF = np.where(wr, (eR - eta) / dx, np.where(wl, (eta - eL) / dx, -self.Ix[None, :]))
        gxB = np.where(wl, (eta - eL) / dx, np.where(wr, (eR - eta) / dx, -self.Ix[None, :]))
        gx = np.where(uu >= 0, gxF, gxB)
        gyF = np.where(wn, (eN - eta) / dx, np.where(ws, (eta - eS) / dx, 0.0))
        gyB = np.where(ws, (eta - eS) / dx, np.where(wn, (eN - eta) / dx, 0.0))
        gy = np.where(vv >= 0, gyF, gyB)
        tq = np.where(V > 1e-6, -(gx * uu + gy * vv) / np.where(V > 1e-6, V, 1.0), 0.0)
        C = S / hs
        As = np.zeros((ny, nx)); dmB = np.zeros((ny, nx))
        for k in range(K):
            As = As + A[:, :, k]; dmB = dmB + A[:, :, k] * self.d[k]
        dRep = np.where(C > 1e-6, self.dfD, np.where(As > 1e-9, dmB / np.where(As > 1e-9, As, 1.0), self.dm0))
        dRep = np.where(dRep > 0, dRep, self.dm0)
        Ci = self.cinf_t(tq, hs, dRep)
        # 侵食
        ero = act & (C < Ci) & (As >= 1e-9)
        dB = dmB / np.where(As > 1e-9, As, 1.0)
        dzE = de * (Ci - C) / (Cst - Ci) * hs * V / np.where(dB > 0, dB, 1.0) * dt
        eroMode = DB.get('erosion') or 'sat'
        if eroMode != 'sat':
            # 不飽和河床の侵食速度式（高橋・里深）：i/√(gh) = Kt sin^{3/2}θ {1 − (σ−ρ)/ρ C (tanφ/tanθ − 1)}^{1/2} (tanφ/tanθ − 1)(C∞ − C) h/d
            Kt = DB.get('Kt') or 0.0019; tanPd = DB.get('tanPhi') or 0.7; sgm = P['sg'] * 1000 / (DB.get('rho') or 1000.0)
            Wd = self.W if gwOn else np.zeros((ny, nx))
            unsat = np.full((ny, nx), eroMode == 'unsat') | ((eroMode == 'auto') & gwOn & self.por & (Wd < 0.9 * (1 - Cst) * np.maximum(0.0, z - zfix)))
            use = unsat & (tq > 1e-4) & (tq < tanPd)
            tqs = np.where(use, tq, 0.5)
            tr_ = tanPd / tqs - 1; br_ = 1 - (sgm - 1) * C * tr_; sn_ = tqs / np.sqrt(1 + tqs * tqs)
            dzU = np.where(br_ > 0, Kt * np.sqrt(G * hs) * np.power(sn_, 1.5) * np.sqrt(np.maximum(br_, 0.0)) * tr_ * (Ci - C) * hs / np.where(dB > 0, dB, 1.0) * dt, 0.0)
            dzE = np.where(use, dzU, dzE)
        dzE = np.minimum(dzE, (Ci * hs - S) / (Cst - Ci)); dzE = np.minimum(dzE, As)
        dzE = np.minimum(dzE, np.maximum(0.0, z - zfix))
        ero &= dzE > 0
        dzE = np.where(ero, dzE, 0.0)
        frac = np.where(ero[..., None], A / np.where(As > 0, As, 1.0)[..., None], 0.0)
        e = dzE[..., None] * frac
        # 堆積
        dep = act & (C > Ci) & (S > 0) & ~ero
        dzD = dd * (C - Ci) / Cst * np.maximum(V, Vmin) * dt
        dzD = np.minimum(dzD, (S - Ci * hs) / (Cst - Ci)); dzD = np.minimum(dzD, S / Cst)
        dep &= dzD > 0
        dzD = np.where(dep, dzD, 0.0)
        gq = dzD[..., None] * np.where(dep[..., None], Cc / np.where(S > 0, S, 1.0)[..., None], 0.0)
        Anew = np.where(ero[..., None], A - e, np.where(dep[..., None], A + gq, Anew))
        Ccn = np.where(ero[..., None], Cc + Cst * e, np.where(dep[..., None], np.maximum(0.0, Cc - Cst * gq), Ccn))
        if gwOn:
            Tq = z - zfix; bq = W / (1 - Cst)
            win = (1 - Cst) * np.maximum(0.0, np.minimum(dzE, bq - (Tq - dzE)))
            pw = np.minimum((1 - Cst) * dzD, np.maximum(0.0, h - S))
            Wn = np.where(ero, W - win, np.where(dep, W + pw, Wn))
            hn = np.where(ero, h + Cst * dzE + win, np.where(dep, h - Cst * dzD - pw, hn))
        else:
            hn = np.where(ero, h + dzE, np.where(dep, h - dzD, hn))
            ng = hn < 0
            if ng.any():
                self.voidAir -= float(hn[ng].sum()) * dx * dx; hn = np.where(ng, 0.0, hn)
        zn = np.where(ero, z - dzE, np.where(dep, z + dzD, zn))
        dzs = np.where(ero, -dzE, np.where(dep, dzD, dzs))
        self.A = Anew.reshape(-1, K); self.Cc = Ccn.reshape(-1, K)
        self.z = zn; self.h = hn; self.W = Wn; self.dzs = self.dzs + dzs
        chg = np.flatnonzero((drp | ero | dep).ravel())
        if chg.size:
            self._fix_layers(chg)

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
            if D.get('op') is not None:
                op = D['op']; out[-1].update(type=op['type'], l=op['l'], n=op['n'])
                if op['hb'] > 0:
                    out[-1]['hb'] = round(float(op['hb']), 3)
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
                extra = ''
                if d.get('type') and d['type'] != 'closed':
                    extra = ' ty=%d l=%s n=%d' % (_DCODE[d['type']], repr(float(d['l'])), int(d['n'])) + (' hb=%s' % repr(float(d['hb'])) if d.get('hb') else '')
                f.write('# dam ' + ' '.join('%s=%s' % (k, repr(float(d[k]))) for k in ('x', 'H', 'b', 'ang', 'y', 'zc', 'len', 'so')) + extra + '\n')
            mask = (self.kind == 0) | ((self.kind == 2) & self.fl)       # 谷底の印（1 が谷底）
            for j in range(self.ny - 1, -1, -1):
                f.write('# floor ' + ''.join('1' if v else '0' for v in mask[j]) + '\n')
            for j in range(self.ny - 1, -1, -1):
                f.write(','.join('%.3f' % v for v in zz[j]) + '\n')

    # ------------------------------------------------------------ 時間ごとの結果の書き出し
    def _snap_label(self, t):
        ev = float(self.P.get('snapEvery') or 600); w = len(str(int(round(self.duration() / 60))))
        if abs(ev / 60 - round(ev / 60)) < 1e-9:
            return str(int(round(t / 60))).zfill(w)
        return '%.1f' % (t / 60)

    def snapshots_long_csv(self, path):
        """時間ごとの結果を縦長の CSV に書く（1 行が 1 つの時刻・1 つの移動床セル）"""
        nx, ny, dx = self.nx, self.ny, self.dx
        jj, ii = np.nonzero(self.kind == 0)
        order = np.lexsort((jj, ii)); jj, ii = jj[order], ii[order]          # x の順、同じ x では y の順
        with open(path, 'w', encoding='utf-8-sig', newline='\r\n') as f:
            f.write('# 時間ごとの結果（記録の間隔 %g 分、格子間隔 %g m）。移動床のセルだけ\n' % ((self.P.get('snapEvery') or 0) / 60, dx))
            f.write('# x：上流端からの流下方向の距離。y：右岸側の端からの距離（左岸側が正）。河床変動量は計算開始時の河床からの変化（＋が堆積）。水深 0 は乾いたセル\n')
            f.write('時刻_分,x_m,y_m,河床高_m,河床変動量_m,水深_m,水位_m,流速_m/s,表層平均粒径_mm\n')
            for sn in self.snaps:
                tm = '%.2f' % (sn['t'] / 60)
                dz = sn['dz'][jj, ii]; h = sn['h'][jj, ii]; V = sn['V'][jj, ii]; dm = sn['dm'][jj, ii]
                z = self.z0[jj, ii] + dz
                for k in range(len(jj)):
                    hh = float(h[k]) if h[k] > 1e-3 else 0.0
                    f.write('%s,%.2f,%.2f,%.3f,%.4f,%.3f,%s,%.3f,%s\n' % (tm, (ii[k] + 0.5) * dx, (jj[k] + 0.5) * dx, z[k], dz[k], hh,
                            ('%.3f' % (z[k] + hh)) if hh > 0 else '', V[k] if hh > 0 else 0.0, ('%.2f' % (dm[k] * 1000)) if np.isfinite(dm[k]) else ''))

    def write_snapshots(self, folder):
        """時間ごとの結果をフォルダに書く（ブラウザ版の zip と同じ構成）"""
        os.makedirs(folder, exist_ok=True)
        for sub in ('河床高', '河床変動量', '水深'):
            os.makedirs(os.path.join(folder, sub), exist_ok=True)
        S = self.snaps; nd = len(self.dams)
        with open(os.path.join(folder, '集計_時系列.csv'), 'w', encoding='utf-8-sig', newline='\r\n') as f:
            f.write('# 時系列の集計。土砂量は実質体積、堆砂量は空隙込み（河床の上昇・低下の体積）\n')
            f.write(','.join(['時刻_分', '流入流量_m3/s', '流出流量_m3/s', '流入土砂量_m3/s', '流出土砂量_m3/s', '累計流入土砂_m3', '累計流出土砂_m3', '堆砂量合計_m3']
                             + ['堰堤%dの上流_m3' % (q + 1) for q in range(nd)] + ['最下流堰堤より下流_m3']) + '\n')
            for sn in S:
                f.write(','.join(['%.2f' % (sn['t'] / 60), '%.3f' % sn['Qin'], '%.3f' % sn['Qout'], '%.4e' % sn['QsIn'], '%.4e' % sn['QsOut'],
                                  '%.2f' % sn['sedIn'], '%.2f' % sn['sedOut'], '%.2f' % sn['stor']] + ['%.2f' % v for v in sn['byDam'][:nd + 1]]) + '\n')
        fl = self.fl & (self.kind == 0)
        with open(os.path.join(folder, '縦断_時刻別.csv'), 'w', encoding='utf-8-sig', newline='\r\n') as f:
            f.write('# 時刻別の縦断。谷底（移動床）のセルの最深値と平均値。堰堤の列は空欄\n')
            f.write(','.join(['x_m'] + ['最深河床高_%s分' % self._snap_label(sn['t']) for sn in S] + ['平均河床高_%s分' % self._snap_label(sn['t']) for sn in S]) + '\n')
            n = fl.sum(axis=0)
            zs = [np.where(fl, self.z0 + sn['dz'], np.nan) for sn in S]
            with np.errstate(invalid='ignore'):
                mins = [np.nanmin(np.where(fl, z, np.inf), axis=0) for z in zs]; means = [np.nansum(z, axis=0) / np.maximum(n, 1) for z in zs]
            for i in range(self.nx):
                row = ['%.2f' % ((i + 0.5) * self.dx)] + [('%.3f' % m[i]) if n[i] else '' for m in mins] + [('%.3f' % m[i]) if n[i] else '' for m in means]
                f.write(','.join(row) + '\n')
        zsave, tsave = self.z, self.t
        try:
            for sn in S:
                lb = self._snap_label(sn['t']); tm = sn['t'] / 60
                self.z = self.z0 + sn['dz']; self.t = sn['t']
                self.bed_csv(os.path.join(folder, '河床高', '河床高_%s分.csv' % lb))
                for name, arr, fmt, ttl in (('河床変動量', sn['dz'], '%.4f', '河床変動量 [m]（出水の時刻 %.1f 分、＋が堆積）' % tm),
                                            ('水深', np.where(sn['h'] > 1e-3, sn['h'], 0.0), '%.3f', '水深 [m]（出水の時刻 %.1f 分）' % tm)):
                    with open(os.path.join(folder, name, '%s_%s分.csv' % (name, lb)), 'w', encoding='utf-8-sig', newline='\r\n') as f:
                        f.write('# %s　1行目が左岸側、左の列が上流　格子間隔 %g m\n' % (ttl, self.dx))
                        for j in range(self.ny - 1, -1, -1):
                            f.write(','.join(fmt % v for v in arr[j]) + '\n')
        finally:
            self.z, self.t = zsave, tsave


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
        eb, ew = sim.balance()
        print('  土砂の体積収支の誤差 %.1e m3、水の体積収支の誤差 %.1e m3' % (eb, ew))
        for L_ in getattr(sim, 'slideLog', []):
            print('  崩壊 %s：%.0f 秒、%s、%.0f m3（安全率 %.2f、地下水の割合 m %.2f）' % (L_['id'], L_['t'], {'fluid': '土石流化', 'block': '土塊', 'slip': 'すべり破壊'}.get(L_['mode'], L_['mode']), L_['vol'], L_.get('Fs', float('nan')), L_.get('m', float('nan'))))
        if getattr(sim, 'rainVol', 0) > 0:
            print('  格子に降った雨 %.0f m3（損失 %.0f m3）' % (sim.rainVol, sim.rainLoss))
        if getattr(sim, 'progCount', 0):
            print('  天然ダムの進行性崩壊：%d 回、%.0f m3' % (sim.progCount, sim.progVol))
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
            if 'ty' in d and 0 <= int(d['ty']) < len(_DNAME):
                d['type'] = _DNAME[int(d.pop('ty'))]
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
_DNAME = ['closed', 'slit', 'mslit', 'grid', 'sslit']          # 堰堤の形式（CSV の ty= の番号）
_DCODE = {k: i for i, k in enumerate(_DNAME)}


def _parse_dams(s):
    out = []
    for tok in s.split(','):
        tok = tok.strip()
        if not tok:
            continue
        a = tok.split(':')
        D = {'x': float(a[0]), 'H': float(a[1]), 'b': float(a[2]) if len(a) > 2 else 0.0}
        if len(a) > 3 and a[3] and a[3] != 'closed':     # 形式：slit / mslit / grid / sslit、開口幅 l、開口数 n
            D['type'] = a[3]; D['l'] = float(a[4]) if len(a) > 4 else 1.0; D['n'] = int(float(a[5])) if len(a) > 5 else 1
        out.append(D)
    return out


def _parse_pairs(s):
    return [[float(v) for v in tok.split(':')] for tok in s.split(',') if tok.strip()]


def _parse_tribs(s):
    out = []
    for tok in s.split(','):
        a = tok.strip().split(':')
        if len(a) < 4:
            continue
        out.append({'x': float(a[0]), 'side': 'R' if a[1].upper().startswith('R') else 'L', 'B': float(a[2]), 'I': float(a[3]), 'len': 40.0,
                    'r': float(a[4]) if len(a) > 4 else 0.3, 'rs': float(a[5]) if len(a) > 5 else 1.0, 'lag': float(a[6]) if len(a) > 6 else 0.0,
                    'ang': float(a[7]) if len(a) > 7 else 90.0})
    return out


def dem_params(P, a):
    """DEM から流域を切り出し、供給点・雨・崩壊・天然ダムを計算条件に入れる（ブラウザ版の「DEM から流域を選ぶ」と同じ）"""
    import sabo2d_dem as SD
    D = SD.synth_dem() if a.dem[0] == 'synth' else SD.load(a.dem)
    A = SD.analyze(D, a.dem_thr * 1e6)

    def cell_of(X, Y):
        c = int(np.clip(math.floor((X - D.x0) / D.cs), 0, D.nx - 1)); r = int(np.clip(math.floor((D.yTop - Y) / D.cs), 0, D.ny - 1))
        return r * D.nx + c

    def near_stream(c, rr=4):
        x, y = c % D.nx, c // D.nx; best, bd = -1, 1e18
        for dy in range(-rr, rr + 1):
            for dx_ in range(-rr, rr + 1):
                xx, yy = x + dx_, y + dy
                if 0 <= xx < D.nx and 0 <= yy < D.ny and A.stream[yy * D.nx + xx] and dx_ * dx_ + dy * dy < bd:
                    bd, best = dx_ * dx_ + dy * dy, yy * D.nx + xx
        return best
    if a.dem_outlet:
        X, Y = [float(v) for v in a.dem_outlet.split(',')]
        out = near_stream(cell_of(X, Y))
        mask = SD.upstream(D, A, out)
    else:
        b = max(A.basins, key=lambda b: b['area']); out = b['outlet']; mask = A.bid == b['id']
    dom = SD.prepare_domain(D, A, mask, out)
    P['terrain'] = SD.terrain_params(dom); P['csvSlope'] = 0; P['dams'] = []; P['tribs'] = []; P['prefill'] = 0
    P['dx'] = max(P['dx'], D.cs)
    dur = a.dem_minutes * 60
    P['dem'] = {'soilD': 1.0, 'chW': 10.0, 'ks': 1e-5, 's0': 0.3, 'duration': dur}
    heads = sorted([c for c in A.heads if mask[c]], key=lambda c: -A.acc[c])[:max(0, a.dem_heads)]
    db = P.get('transport') == 'debris'
    srcs = []
    for c in heads:
        xy = dom.toXY(c); u = dom.dirXY(c)
        H = ({'mode': 'debris', 'Qbase': 0.2, 'Qpeak': a.dem_q, 't0': 30.0, 'tRise': 60.0, 'tFall': 300.0, 'Cd': a.dem_c, 'Cs': 0.0} if db else
             {'mode': 'table', 'table': [[0, 0.2, 0.0], [30, 0.2, 0.0], [90, a.dem_q, a.dem_q * 0.01], [390, 0.2, 0.0], [1e9, 0.2, 0.0]]})
        srcs.append({'x': xy[0], 'y': xy[1], 'ux': u[0], 'uy': u[1], 'I': SD.local_slope(D, A, c, 3), 'hydro': H})
    P['sources'] = srcs
    if a.rain:
        tab = _parse_pairs(a.rain)
        lp = a.rain_loss.split(':')
        loss = {'type': 'f', 'f': float(lp[1]) if len(lp) > 1 else 0.8}
        if lp[0] == 'f1rsa':
            loss = {'type': 'f1rsa', 'f1': float(lp[1]), 'Rsa': float(lp[2]) if len(lp) > 2 else 100.0}
        elif lp[0] == 'ks':
            loss = {'type': 'ks', 'ks': float(lp[1]), 'D': float(lp[2]) if len(lp) > 2 else 1.0, 'n': 0.45, 'th0': 0.25}
        rp = a.rain_route.split(':')
        route = {'type': 'storage', 'K': float(rp[1]), 'P': float(rp[2]) if len(rp) > 2 else 0.6} if rp[0] == 'storage' else {'type': 'rational', 'tc': float(rp[1]) if len(rp) > 1 else 0.0}
        if a.rain_method == 'grid':
            P['rain'] = {'on': True, 'method': 'grid', 'mode': 'uniform', 'table': tab, 'frames': [], 'loss': loss}
        else:
            cells, Q, dtq = SD.lateral_inflow(D, A, mask, tab, loss, route, 60.0, dur)
            P['latIn'] = {'dt': dtq, 'pts': [dom.toXY(c) + dom.dirXY(c) + [SD.local_slope(D, A, c, 3)] for c in cells], 'Q': Q.tolist()}
            P['rain'] = {'on': True, 'method': 'lumped', 'mode': 'uniform', 'table': tab, 'frames': [], 'loss': loss}
    P['slides'] = []
    for s_ in a.slide or []:
        f = s_.split(',')
        c = cell_of(float(f[0]), float(f[1])); xy = dom.toXY(c)
        deep = len(f) > 4 and f[4] == 'deep'
        P['slides'].append({'x': xy[0], 'y': xy[1], 'R': float(f[2]) if len(f) > 2 else 25.0, 'D': float(f[3]) if len(f) > 3 else 1.5, 'type': 'deep' if deep else 'shallow',
                            'trigger': f[5] if len(f) > 5 else 'time', 't': float(f[6]) if len(f) > 6 else 30.0, 'move': f[7] if len(f) > 7 else 'auto',
                            'phi': 30.0 if deep else 35.0, 'c': 30.0 if deep else 3.0, 'm0': 0.9, 'thDep': 10.0, 'CL': 0.5})
    P['ldams'] = []
    for s_ in a.ldam or []:
        f = s_.split(','); c = near_stream(cell_of(float(f[0]), float(f[1])))
        if c < 0:
            continue
        xy = dom.toXY(c); u = dom.dirXY(c)
        P['ldams'].append({'x': xy[0], 'y': xy[1], 'ux': u[0], 'uy': u[1], 'H': float(f[2]) if len(f) > 2 else 15.0, 'Wc': float(f[3]) if len(f) > 3 else 10.0, 'up': 30.0, 'dn': 25.0})
    if P['ldams']:
        P['gw']['on'] = True; P['ldam'].update(k=a.ldam_k, fail=a.ldam_fail); P['debris']['erosion'] = 'auto'
    if P['slides'] or P['ldams']:
        P['transport'] = 'debris'
    print('DEM：%d × %d（格子 %.2f m）、流域 %.3f km²、谷頭 %d 点（供給点 %d）、回転 %d × 90°' % (D.nx, D.ny, D.cs, mask.sum() * D.cs ** 2 / 1e6, sum(1 for c in A.heads if mask[c]), len(srcs), dom.rot))


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
    ap.add_argument('--dams', help='堰堤 "x:堤高:水通し幅[:形式:開口幅:開口数],..."（水通し幅 0 で全幅。形式は slit / mslit / grid / sslit）')
    ap.add_argument('--slope-table', help='勾配の屈折 "x:n,x:n,..."（その位置から下流を勾配 1/n にする）')
    ap.add_argument('--tribs', help='支川 "x:L|R:幅:n:流量比:土砂量比:遅れ秒:合流角度,..."（勾配は 1/n、合流角は度で 90 が直角）')
    ap.add_argument('--gw', action='store_true', help='伏流（河床への浸透と地下水の流れ）を計算する')
    ap.add_argument('--gw-k', type=float, help='透水係数 [m/s]（省略時は粒度から Kozeny–Carman 式）')
    ap.add_argument('--debris', action='store_true', help='土石流の一波を流す（高橋の式。流砂モデルも土石流になる）')
    ap.add_argument('--debris-peak', type=float, default=100.0, help='土石流のピーク流量 [m3/s]（混合物）')
    ap.add_argument('--debris-c', type=float, default=0.4, help='土石流の土砂濃度（容積）')
    ap.add_argument('--debris-rise', type=float, default=30.0, help='立ち上がり時間 [s]')
    ap.add_argument('--debris-fall', type=float, default=120.0, help='減衰時間 [s]')
    ap.add_argument('--debris-iup', type=float, default=0.0, help='流入点より上流の勾配 1/n（0 で上流端の勾配）')
    ap.add_argument('--debris-minutes', type=float, default=10.0, help='土石流の計算時間 [分]')
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
    ap.add_argument('--snap-every', type=float, help='時間ごとの結果を記録する間隔 [分]（出水の時刻。既定 10、0 で記録しない）')
    ap.add_argument('--params', help='ブラウザ版で保存した計算条件の JSON（DEM の流域・供給点・雨・崩壊もそのまま再現）')
    ap.add_argument('--dem', nargs='+', help='DEM のファイル（.asc / .xyz / .tif / 基盤地図情報 .xml、複数可）。"synth" で例の DEM')
    ap.add_argument('--dem-thr', type=float, default=0.02, help='谷筋とみなす集水面積 [km2]')
    ap.add_argument('--dem-outlet', help='流域の出口 "X,Y"（DEM の座標。省略で最も大きい流域）')
    ap.add_argument('--dem-heads', type=int, default=3, help='供給点にする谷頭の数（集水面積の大きい順）')
    ap.add_argument('--dem-q', type=float, default=5.0, help='供給点のピーク流量 [m3/s]（土石流なら混合物）')
    ap.add_argument('--dem-c', type=float, default=0.3, help='供給点の土石流の濃度（--debris と一緒に）')
    ap.add_argument('--dem-minutes', type=float, default=30.0, help='DEM の流域の計算時間 [分]')
    ap.add_argument('--rain', help='一様な雨 "分:mm/h,分:mm/h,..."（その時刻から次の時刻までの強度）')
    ap.add_argument('--rain-method', choices=['lumped', 'grid'], default='lumped', help='谷への集め方（斜面で流出計算 / 格子に直接降らせる）')
    ap.add_argument('--rain-loss', default='f:0.8', help='損失 "f:0.8" / "f1rsa:f1:Rsa" / "ks:ks:D"')
    ap.add_argument('--rain-route', default='rational', help='斜面の流出 "rational[:到達時間 分]" / "storage:K:P"')
    ap.add_argument('--slide', action='append', help='斜面崩壊 "X,Y,半径,深さ,shallow|deep,fs|time,時刻s,auto|fluid|block"（DEM の座標）')
    ap.add_argument('--ldam', action='append', help='天然ダム "X,Y,高さ[,天端幅]"（DEM の座標。谷筋に置く）')
    ap.add_argument('--ldam-k', type=float, default=1e-3, help='天然ダムの堤体の透水係数 [m/s]')
    ap.add_argument('--ldam-fail', choices=['auto', 'slip'], default='auto', help='決壊の形（透水性で決まる / すべり破壊）')
    ap.add_argument('--out', default='result', help='図と CSV を書き出すフォルダ')
    a = ap.parse_args(argv)
    P = preset(a.preset)
    if a.params:
        import json
        P = default_params(); P.update(json.load(open(a.params, encoding='utf-8')))
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
    if a.slope_table:
        P['slopeTable'] = _parse_pairs(a.slope_table)
    if a.tribs:
        P['tribs'] = _parse_tribs(a.tribs)
    if a.gw:
        P['gw'] = {'on': True, 'k': a.gw_k or 0.0, 'kMax': 0.05, 's0': 0.0}
    if a.fill is not None:
        P['prefill'] = a.fill
    H = P['hydro']
    if a.qpeak:
        H['Qpeak'] = a.qpeak
    if a.qbase:
        H['Qbase'] = a.qbase
    if a.minutes:
        H['duration'] = a.minutes * 60; H['tPeak'] = a.minutes * 20
    if a.debris:
        H.update(mode='debris', Qbase=a.qbase or 1.0, Qpeak=a.debris_peak, t0=30.0, tRise=a.debris_rise, tFall=a.debris_fall,
                 Cd=a.debris_c, Iup=a.debris_iup, Cs=0.0, duration=a.debris_minutes * 60)
        P['transport'] = 'debris'
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
    if a.snap_every is not None:
        P['snapEvery'] = a.snap_every * 60
    if a.secondary or a.secondary_mom:
        SC = P.setdefault('sec', {}); SC['bed'] = bool(a.secondary); SC['mom'] = bool(a.secondary_mom)
    if a.boulders is not None:
        P.setdefault('boulders', {})['L'] = a.boulders
    if a.boulder_mode:
        P.setdefault('boulders', {})['mode'] = a.boulder_mode
    if a.dem:
        dem_params(P, a)
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
    if sim.snaps:
        sim.snapshots_long_csv(os.path.join(a.out, 'snapshots_long.csv'))
        sim.write_snapshots(os.path.join(a.out, 'snapshots'))
        print('時間ごとの結果（%d 回）：%s/snapshots_long.csv（縦長）、%s/snapshots/（時刻ごとの行列・縦断・集計）' % (len(sim.snaps), a.out, a.out))


if __name__ == '__main__':
    main()
