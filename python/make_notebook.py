# -*- coding: utf-8 -*-
"""sabo2d.py から Jupyter Notebook（sabo2d.ipynb）を作る"""
import re
import nbformat as nbf

src = open('sabo2d.py', encoding='utf-8').read()
src = src.split('\n# =============================================================================\n# コマンドライン')[0]
parts = re.split(r'\n# =============================================================================\n# (.+)\n# =============================================================================\n', src)
head = parts[0]; secs = list(zip(parts[1::2], parts[2::2]))
doc = head.split('"""')[1]
imports = head.split('"""', 2)[2].strip()

nb = nbf.v4.new_notebook()
C = nb.cells
md = lambda t: C.append(nbf.v4.new_markdown_cell(t))
code = lambda t: C.append(nbf.v4.new_code_cell(t.strip('\n')))

md("""# 連続砂防堰堤 平面2次元河床変動計算（Jupyter Notebook 版）

ブラウザ版の「河床変動シミュレータ」と同じ計算を、Python（NumPy）で書いたものです。
同じ条件なら、地形・初期の乱数・流れの計算はブラウザ版と同じ値になります（確認した差は 10⁻¹⁰ m 程度）。

**使い方**

1. 上から順にセルを実行します（メニューの「Run All」でも構いません）。
2. 「計算条件」のセルで、モデル渓流・ハイドログラフ・堰堤を変えます。
3. 「計算」のセルを実行すると、出水が終わるまで計算します。モデル渓流 A で 1 時間の出水なら、2 分ほどかかります。
4. その下のセルで、平面図・河床変動量・縦断図・横断図・鳥瞰図・時間変化を描きます。

測量値（横断・縦断）や格子 CSV を読み込む例は、最後のほうにあります。

必要なもの：Python 3.8 以上、NumPy、Matplotlib。

<details><summary>モデルの説明（クリックで開く）</summary>

```
""" + doc.strip() + """
```
</details>""")
code('%matplotlib inline\n' + imports + '\nimport matplotlib.pyplot as plt')
titles = {'条件（既定値）': '## 1. 計算条件の既定値とモデル渓流\n\n`default_params()` が既定値、`preset(\'A\')` などがモデル渓流です。',
          '乱数（ブラウザ版と同じ値を出す）': '## 2. 乱数と補助の関数\n\nブラウザ版と同じ乱数列を出すので、同じ地形ができます。',
          '計算本体': '## 3. 計算本体（Sim）\n\n初期化（地形・堰堤・満砂・土砂の層）、流れ（HLL＋静水圧再構築）、流砂（芦田・道上式）、交換層、斜面の崩落。',
          '計算を進める': '## 4. 計算を進める関数と堆砂量の表',
          '地形ファイルの読み込み': '## 5. 地形ファイルの読み込み\n\n格子 CSV（`read_grid`）、横断・縦断の測量値（`read_survey`）、ファイルの種類の判定（`classify_csv`）。文字コードは UTF-8 と Shift_JIS のどちらでも読めます。',
          '図（Matplotlib）': '## 6. 図'}
for t, body in secs:
    md(titles.get(t, '## ' + t))
    code(body)

md("""---
## 7. 計算条件

`preset('A')`：床固工群の渓流（早月川を参考にした架空の地形）　`preset('B')`：狭窄部と広い堆砂域（南股入を参考）　`preset('S')`：単純な直線水路

主な項目（単位は m、s、m³/s）

| 項目 | 意味 |
| --- | --- |
| `P['hydro']['Qbase']`, `['Qpeak']` | 基底流量、ピーク流量 |
| `P['hydro']['tPeak']`, `['duration']` | ピーク時刻、継続時間 [s] |
| `P['hydro']['Cs']` | 供給土砂濃度 Qs/Q（実質体積） |
| `P['dams']` | 堰堤 `[{'x': 位置, 'H': 堤高, 'b': 水通し幅（0 で全幅）}, ...]` |
| 堰堤の任意の項目 | `'ang'` 回転角 [°]（x 軸に直角な向きから反時計回り）、`'y'` 中心 y [m]、`'zc'` 天端標高 [m]、`'len'` 堤長 [m]。与えなければ自動 |
| `P['prefill']` | 0 より大きいと、堰堤の上流をはじめから満砂にする（堆砂勾配と元の勾配の比） |
| `P['bend']['table']` | 平面形 `[[x, 中心線の横ずれ, 川幅], ...]`（川幅の列で狭窄・拡幅） |
| `P['valley']` | 澪筋の深さ・幅、谷底の起伏、谷壁の高さ、岸の出入り、支渓の扇状地 |
| `P['grains']` | 粒径階 `[{'d': 粒径, 'bed': 河床の割合, 'feed': 供給の割合}, ...]` |
| `P['tcMethod']` | 限界掃流力：`'atm'` 芦田・高橋・水山、`'iwagaki'` 岩垣、`'const'` 一定値 |
| `P['morfac']` | 河床変動の加速係数 |""")
code("""setup_japanese_font()

P = preset('A')                       # 'A'・'B'・'S'
P['hydro'].update(Qbase=10, Qpeak=60, tPeak=20 * 60, duration=60 * 60, Cs=0.005)
# 堰堤を変えるとき（例）：
# P['dams'] = [{'x': x, 'H': 3.0, 'b': 36.0} for x in range(75, 600, 60)]
# P['prefill'] = 0.0                  # 空の堰堤から始める

sim = Sim(P)
print('格子 %d × %d（%d セル、幅 %g m）、堰堤 %d 基' % (sim.nx, sim.ny, sim.N, sim.dx, len(sim.dams)))
for w in sim.warn:
    print('注意：' + w)
plot_plan(sim, title='初期の地形')
plt.show()""")
md('## 8. 計算\n\n出水が終わるまで計算します。30 秒ごとに途中の状態を表示します。')
code('run(sim)\nprint_dam_table(sim)')
md('## 9. 結果の図')
code("""plot_plan(sim); plt.show()
plot_plan(sim, var='dz', title='河床変動量'); plt.show()
plot_plan(sim, var='dm', title='表層の平均粒径'); plt.show()""")
code("""plot_profile(sim); plt.show()
plot_profile(sim, detrend=True); plt.show()        # 初期勾配を差し引いて、堆積・侵食を見やすく
plot_xsec(sim, x=300); plt.show()                  # x を変えると別の横断""")
code("""plot_3d(sim, vex=2, elev=35, azim=-60); plt.show()
plot_series(sim); plt.show()""")
code("""# まとめて 1 枚に描いて保存する。計算後の河床は CSV に書き出すと、ブラウザ版の「格子の CSV」で読めます
fig = plot_all(sim, 'result.png')
sim.bed_csv('bed_after.csv')""")
md("""## 10. 計算後の河床から続ける

`sim.continue_params()` は、いまの河床を初期地形にし、堰堤を**前と同じ位置・向き・天端標高**で置いた条件を返します（谷底の範囲もそのまま。粒度は初期の河床材料に戻ります）。
`sim.bed_csv('ファイル名.csv')` で書き出した CSV にも堰堤の諸元と谷底の範囲が入っていて、`read_grid` → `use_terrain` で読むと同じ状態に戻ります。ブラウザ版でも読めます。""")
code("""P3 = sim.continue_params()
# 堰堤を手で動かすとき（例）：1 基目を 10° 回し、堤長を 70 m にする
# P3['dams'][0]['ang'] += 10; P3['dams'][0]['len'] = 70
P3['hydro'].update(Qbase=10, Qpeak=60, tPeak=20 * 60, duration=60 * 60)
sim3 = Sim(P3)
print('堰堤：', [(round(D['x']), round(np.degrees(D['th']), 1), round(D['crest'], 2)) for D in sim3.dams])
run(sim3)
print_dam_table(sim3)
plot_plan(sim3, var='dz', title='2 回目の出水での河床変動量'); plt.show()""")
md("""## 11. 測量値・格子 CSV を読み込む例

配布したテスト CSV（`test_cross_sections.csv`、`test_longitudinal.csv`、`test_grid.csv`）をこのノートと同じフォルダに置いて実行します。
`classify_csv` で、ファイルの種類（横断・縦断・格子）を確かめられます。""")
code("""for f in ('test_cross_sections.csv', 'test_longitudinal.csv', 'test_grid.csv'):
    try:
        print(f, '→', classify_csv(f))
    except FileNotFoundError:
        print(f, '：ファイルがありません')""")
code("""# 横断・縦断の測量値から地形を作る（格子 CSV なら T = read_grid('test_grid.csv', cs=2.0)）
T = read_survey('test_cross_sections.csv', 'test_longitudinal.csv', dx=2.0, align='thalweg')
P2 = use_terrain(default_params(), T, floorH=3.0)     # 谷底とみなす高さ floorH [m]
P2['dx'] = 2.0
P2['dams'] = [{'x': 150, 'H': 3.0, 'b': 0.0}, {'x': 300, 'H': 3.0, 'b': 0.0}]   # 位置は地形の最上流を 0 とした距離
P2['hydro'].update(Qbase=10, Qpeak=50, tPeak=15 * 60, duration=45 * 60)
sim2 = Sim(P2)
print('格子 %d × %d、堰堤 %d 基' % (sim2.nx, sim2.ny, len(sim2.dams)), sim2.warn)
plot_plan(sim2, title='読み込んだ地形'); plt.show()""")
code("""run(sim2)
print_dam_table(sim2)
plot_plan(sim2); plt.show()
plot_profile(sim2, detrend=True); plt.show()""")
md("""## 12. 側岸侵食の比べ方

谷底いっぱいに流れる狭い水路（幅 10 m、勾配 1/40、流量 25 m³/s を 30 分）で、側岸を固定したときと、段丘の側岸を削れるようにしたときを比べます。

- `P['bank']['erodH']`：谷底より上の壁を、この高さまで侵食できるようにする [m]
- `P['bank']['reposeWet']`, `['reposeDry']`：水中・陸上の安息角 [度]（水中 30°：高橋ら 1996、陸上は XBeach・Delft3D の wetslp／dryslp と同じ考え方）
- `P['bank']['dryEro']`：濡れたセルが削れたとき、隣の乾いた側岸を同時に削る割合（Delft3D の dry-cell erosion）
- `P['transport'] = 'noneq'`：侵食・堆積速度式による非平衡の流砂（高橋・里深 1999、高橋ら 1998）
- `P['boulders']`：岸に動かない巨礫を間隔 `L` [m] で置く（高橋・里深・奥村 1996）""")
code("""def narrow_channel(bank, transport='eq'):
    P = default_params()
    P.update(L=240, B=10, slope=1 / 40, crossSlope=0.02, dx=1.0, sideCells=8, sideGrad=1.2, dams=[], morfac=4, transport=transport)
    P['valley'] = {'chD': 0, 'chW': 0, 'chWl': 0, 'relief': 0.1, 'wallH': 3, 'wallRough': 0, 'fans': []}
    P['bank'] = bank
    P['hydro'] = {'mode': 'const', 'Qbase': 25, 'Qpeak': 25, 'tPeak': 0, 'duration': 1800, 'Cs': 0.0, 'table': []}
    return P

cases = {'側岸は固定': narrow_channel({'erodH': 0, 'reposeWet': 30, 'reposeDry': 50, 'dryEro': 0}),
         '側岸を侵食（高さ 2 m）': narrow_channel({'erodH': 2, 'reposeWet': 30, 'reposeDry': 50, 'dryEro': 0.5})}
sims = {}
for name, P in cases.items():
    s = Sim(P)
    while s.phase != 'done':
        s.step()
    i0, i1 = int(60 / s.dx), int(200 / s.dx)
    width = ((s.h[:, i0:i1] > 0.05).sum(axis=0) * s.dx).mean()
    vb = ((s.z0 - s.z) * (s.bankC & (s.z < s.z0))).sum() * s.dx ** 2
    print('%s：水面幅 %.1f m、側岸の侵食 %.0f m³' % (name, width, vb))
    sims[name] = s""")
code("""fig, axs = plt.subplots(2, 1, figsize=(12, 5))
for ax, (name, s) in zip(axs, sims.items()):
    plot_plan(s, ax=ax, var='dz', title=name)
plt.tight_layout(); plt.show()""")
md("""## 13. 湾曲部の2次流と長谷川の崩落モデル（長田ら 1996）

蛇行した水路（幅 12 m、勾配 1/100、20 m³/s を 30 分）で、2次流を入れない場合と入れた場合を比べます。

- `P['sec']['bed']`：河床付近の流れを曲率の内側へ向け、掃流砂の向きを変える（u_nb/u_sb = N* h/r、N* = 11.5）
- `P['sec']['mom']`：2次流による運動量の分散項（Kalkwijk・de Vriend 1980）。流速の最大が外岸側へ寄る
- `P['bank']['collapse'] = 'hasegawa'`：水際の側岸は根元の洗掘が溜まってから一度に崩れる（長谷川 1983）

**適用条件**：長田らの検証は、勾配 1/100〜1/466 の砂（0.7〜1.4 mm）の水路実験で、掃流砂が主で変化がゆるやかな場合です。弱い曲率（h/r ≪ 1）を前提とするため、h/r を `hrMax`（既定 0.05）で頭打ちにしています。急勾配の渓流では定性的な目安として見てください。""")
code("""def meander(sec, collapse='repose', erodH=0.0):
    P = default_params()
    P.update(L=360, B=12, slope=1 / 100, crossSlope=0, dx=1.0, sideCells=6, sideGrad=1.5, dams=[], morfac=4)
    P['valley'] = {'chD': 0, 'chW': 0, 'chWl': 0, 'relief': 0.0, 'wallH': 3, 'wallRough': 0, 'fans': []}
    P['bend'] = {'mode': 'sine', 'amp': 12, 'wl': 120, 'x0': 60, 'x1': 300}
    P['grains'] = [{'d': 0.01, 'bed': 0.5, 'feed': 0.5}, {'d': 0.03, 'bed': 0.5, 'feed': 0.5}]
    P['hydro'] = {'mode': 'const', 'Qbase': 20, 'Qpeak': 20, 'tPeak': 0, 'duration': 1800, 'Cs': 0.002, 'table': []}
    P['sec'] = sec
    P['bank'] = {'erodH': erodH, 'reposeWet': 30, 'reposeDry': 50, 'dryEro': 0.0, 'collapse': collapse}
    return P

cases = {'2次流なし': meander({'bed': False, 'mom': False}),
         '2次流あり（流砂の向き＋分散項）': meander({'bed': True, 'mom': True}),
         '2次流あり＋長谷川の崩落（側岸 2 m）': meander({'bed': True, 'mom': True}, 'hasegawa', 2.0)}
msims = {}
for name, P in cases.items():
    s = Sim(P)
    while s.phase != 'done':
        s.step()
    msims[name] = s
    print(name, '：長谷川モデルの崩落', s.nHase, '回')""")
code("""fig, axs = plt.subplots(3, 1, figsize=(12, 7.5))
for ax, (name, s) in zip(axs, msims.items()):
    plot_plan(s, ax=ax, var='dz', title=name)
plt.tight_layout(); plt.show()""")
md("""## 14. 時間ごとの結果を取り出す

計算中、出水の時刻で `P['snapEvery']` 秒ごと（既定 600 秒 = 10 分）に、各セルの河床変動量・水深・流速・表層の平均粒径を `sim.snaps` に残します（出水の開始時と終了時も残します）。途中で間隔を変えるときは `sim.set_snap_every(秒)`、記録しないときは `P['snapEvery'] = 0` です。

- `sim.snapshots_long_csv(パス)`：縦長の CSV（1 行が 1 つの時刻・1 つの移動床セル。列は 時刻・x・y・河床高・河床変動量・水深・水位・流速・表層平均粒径）
- `sim.write_snapshots(フォルダ)`：時刻ごとの行列 CSV（河床高・河床変動量・水深）、時刻別の縦断、時系列の集計。河床高の CSV はブラウザ版の「格子の CSV」でも読める

コマンドラインでは `--snap-every 分` で間隔を指定し、`result/snapshots_long.csv` と `result/snapshots/` に書き出します。""")
code("""s = msims['2次流あり（流砂の向き＋分散項）']
print('記録した時刻 [分]：', [round(q['t'] / 60, 1) for q in s.snaps])
s.snapshots_long_csv('snapshots_long.csv')
s.write_snapshots('snapshots')
print('書き出し：snapshots_long.csv、snapshots/')""")
code("""# 時刻ごとの河床変動量の縦断（谷底セルの平均）
fl = s.fl & (s.kind == 0)
x = (np.arange(s.nx) + 0.5) * s.dx
fig, ax = plt.subplots(figsize=(10, 3.5))
for q in s.snaps[1:]:
    dz = np.where(fl, q['dz'], np.nan)
    ax.plot(x, np.nanmean(dz, axis=0), label='%.0f 分' % (q['t'] / 60))
ax.axhline(0, color='0.6', lw=0.8)
ax.set_xlabel('流下方向の距離 x [m]'); ax.set_ylabel('平均河床変動量 [m]'); ax.legend(ncol=4, fontsize=8)
plt.tight_layout(); plt.show()""")
md("""## 15. 土石流の一波・勾配の屈折・透過型堰堤・支川・伏流

勾配 1/4 から 1/12 へ屈折する渓流に、土石流の一波（ピーク 120 m³/s、濃度 0.35）を流します。左岸の急な支川（勾配 1/3）からも、30 秒遅れて土石流が入ります。屈折点の下流に、不透過型と鋼製格子型（純間隔 0.6 m）の堰堤を置いて比べます。

- `P['hydro']['mode'] = 'debris'`：基底流量に三角形の土石流の波（`t0`, `tRise`, `tFall`, 濃度 `Cd`、流入点より上流の勾配 `Iup`）を重ねる。流砂モデルは `P['transport'] = 'debris'`（高橋の式：侵食・堆積速度式、平衡濃度、石礫型・乱流型の抵抗則）
- `P['slopeTable'] = [[x, n], ...]`：位置 x から下流を勾配 1/n にする（勾配の屈折）
- 堰堤の `type`：`'slit'`（スリット）、`'mslit'`（複数スリット、開口数 `n`）、`'grid'`（鋼製格子型）、`'sslit'`（鋼製スリット型）。開口幅 `l`。共通設定は `P['openDam']`（部材の幅 `w`、閉塞の判定 l ≤ `beta`·d95、閉塞が起きる濃度 `Cblk`）
- `P['tribs']`：支川（合流点 x、岸、幅、勾配 1/n、流量比、土砂量比、遅れ、合流角 `'ang'` [°]（90 で直角））
- `P['gw']`：伏流（浸透と地下水の流れ。透水係数は粒度から Kozeny–Carman 式）

**適用条件**：土石流の式は代表粒径 1 つの式で、粒径別の分級は扱いません。閉塞の判定（β = 1.5 は土石流対策技術指針（案）2000、濃度のしきい値は本プログラムの仮定）は簡易なモデルです。""")
code("""def debris_case(dam_type):
    P = default_params()
    P.update(L=400, B=20, slope=1 / 4, dx=2.0, transport='debris')
    P['slopeTable'] = [[220, 12]]
    P['valley'] = {'chD': 0, 'chW': 0, 'chWl': 0, 'relief': 0.0, 'wallH': 6, 'wallRough': 0, 'fans': []}
    P['grains'] = [{'d': 0.01, 'bed': 0.3, 'feed': 0.3}, {'d': 0.05, 'bed': 0.4, 'feed': 0.4}, {'d': 0.2, 'bed': 0.2, 'feed': 0.2}, {'d': 0.5, 'bed': 0.1, 'feed': 0.1}]
    P['hydro'] = {'mode': 'debris', 'Qbase': 1, 'Qpeak': 120, 't0': 30, 'tRise': 20, 'tFall': 100, 'Cd': 0.35, 'Cs': 0, 'duration': 600, 'Iup': 4, 'table': []}
    P['dams'] = [{'x': 260, 'H': 5, 'b': 0} if dam_type == 'closed' else {'x': 260, 'H': 5, 'b': 0, 'type': dam_type, 'l': 0.6}]
    P['tribs'] = [{'x': 150, 'side': 'L', 'B': 8, 'I': 3, 'len': 40, 'r': 0.4, 'rs': 0.4, 'lag': 30, 'ang': 60}]
    P['gw'] = {'on': True, 'k': 0, 'kMax': 0.05, 's0': 0}
    return P

dsims = {}
for name, typ in [('不透過型', 'closed'), ('鋼製格子型（純間隔 0.6 m）', 'grid')]:
    s = Sim(debris_case(typ))
    while s.phase != 'done':
        s.step()
    op = s.dams[0]['op']
    eb, ew = s.balance()
    print('%s：流出土砂 %.0f m³、閉塞 %s、土砂の収支の誤差 %.1e m³、水 %.1e m³' % (name, s.sedOut, ('%.0f 秒' % op['tBlock']) if op and op['blocked'] else 'なし', eb, ew))
    dsims[name] = s""")
code("""fig, axs = plt.subplots(2, 1, figsize=(12, 6))
for ax, (name, s) in zip(axs, dsims.items()):
    plot_plan(s, ax=ax, var='dz', title=name)
plt.tight_layout(); plt.show()""")
dem_src = open('sabo2d_dem.py', encoding='utf-8').read()
md("""## 16. DEM から流域を選ぶ・雨・斜面崩壊・天然ダム

DEM を読み込み（`SD.load('ファイル')`。ESRI ASCII grid、XYZ、GeoTIFF、基盤地図情報の XML）、窪地を埋めて流向・集水面積・谷筋・流域・谷頭を求めます。
ここでは例の DEM（`SD.synth_dem()`、架空の山地流域）で、最も大きい流域を切り出し、

- 谷頭 3 つから土石流の一波を流す（渓床侵食型）
- 一様な雨を降らせ、斜面から谷筋へ合理式（土研式の到達時間）で横流入させる
- 斜面の 1 か所を 30 秒に崩して流動化させる（崩壊型）

を一度に計算します。天然ダム（`P['ldams']`）は谷筋の点に置きます（堤体の透水係数 `P['ldam']['k']` で越流侵食と進行性崩壊が分かれます）。
ブラウザ版で作った計算条件は「計算条件を JSON で保存」で保存し、`P = default_params(); P.update(json.load(open('sabo2d_params.json')))` で同じ計算ができます。

**適用条件**：合理式・貯留関数法は流域スケールの経験的な扱いで、斜面ごとの到達時間は土研式（自然流域、A < 50 km²）を斜面に当てた近似です。
表層崩壊の安全率は無限長斜面、深層は Hovland の 3 次元簡便法、土塊の移動は瞬時移動（エネルギー線）の簡易モデルです。流動化の判定（地下水の割合 m ≥ 0.8）と流動化濃度 0.5 は仮定です。""")
code("import types\nSD = types.ModuleType('sabo2d_dem')\nexec(" + repr(dem_src) + ", SD.__dict__)")
code("""D = SD.synth_dem()
A = SD.analyze(D, thrArea=20000)
b = max(A.basins, key=lambda b: b['area']); mask = A.bid == b['id']
dom = SD.prepare_domain(D, A, mask, b['outlet'])
heads = sorted([c for c in A.heads if mask[c]], key=lambda c: -A.acc[c])[:3]
print('DEM %d × %d、流域 %.2f km²、谷頭 %d 点' % (D.nx, D.ny, mask.sum() * D.cs ** 2 / 1e6, sum(1 for c in A.heads if mask[c])))
fig, ax = plt.subplots(figsize=(9, 6))
Z = D.z.reshape(D.ny, D.nx)
ax.imshow(Z, cmap='terrain'); ax.contour(mask.reshape(D.ny, D.nx), [0.5], colors='r')
ys, xs = np.divmod(np.nonzero(A.stream)[0], D.nx); ax.plot(xs, ys, ',', color='b')
hy, hx = np.divmod(np.array(heads), D.nx); ax.plot(hx, hy, 'o', color='orange')
ax.set_title('流域（赤）、谷筋（青）、供給点（橙）'); plt.show()""")
code("""P = default_params()
P['terrain'] = SD.terrain_params(dom); P['dx'] = D.cs; P['csvSlope'] = 0; P['dams'] = []; P['tribs'] = []
P['transport'] = 'debris'; P['snapEvery'] = 0
P['dem'] = {'soilD': 1.0, 'chW': 10.0, 'ks': 1e-5, 's0': 0.3, 'duration': 300.0}
P['sources'] = [{'x': dom.toXY(c)[0], 'y': dom.toXY(c)[1], 'ux': dom.dirXY(c)[0], 'uy': dom.dirXY(c)[1], 'I': SD.local_slope(D, A, c, 3),
                 'hydro': {'mode': 'debris', 'Qbase': 0.5, 'Qpeak': 10, 't0': 10, 'tRise': 30, 'tFall': 120, 'Cd': 0.3, 'Cs': 0}} for c in heads]
rain = [[0, 50], [10, 80], [20, 0]]
cells, Q, dtq = SD.lateral_inflow(D, A, mask, rain, {'type': 'f', 'f': 0.8}, {'type': 'rational', 'tc': 0}, 60.0, 300.0)
P['latIn'] = {'dt': dtq, 'pts': [dom.toXY(c) + dom.dirXY(c) + [SD.local_slope(D, A, c, 3)] for c in cells], 'Q': Q.tolist()}
c1 = heads[1]; x1, y1 = dom.toXY(c1)
P['slides'] = [{'x': x1 + 15, 'y': y1 + 15, 'R': 35, 'D': 2, 'type': 'shallow', 'trigger': 'time', 't': 30, 'move': 'fluid'}]
s = Sim(P)
while s.phase != 'done':
    s.step()
eb, ew = s.balance()
print('流出土砂 %.0f m³、崩壊 %s、土砂の収支の誤差 %.1e m³、水 %.1e m³' % (s.sedOut, [(round(L['t']), L['mode'], round(L['vol'])) for L in s.slideLog], eb, ew))
fig, axs = plt.subplots(1, 2, figsize=(14, 5))
plot_plan(s, ax=axs[0], var='depth', title='水深'); plot_plan(s, ax=axs[1], var='dz', title='河床変動量')
plt.tight_layout(); plt.show()""")
nb.metadata['kernelspec'] = {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'}
nb.metadata['language_info'] = {'name': 'python'}
nbf.write(nb, 'sabo2d.ipynb')
print('sabo2d.ipynb: %d cells' % len(nb.cells))
