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

ブラウザ版の「連続堰堤 河床変動シミュレータ」と同じ計算を、Python（NumPy）で書いたものです。
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
nb.metadata['kernelspec'] = {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'}
nb.metadata['language_info'] = {'name': 'python'}
nbf.write(nb, 'sabo2d.ipynb')
print('sabo2d.ipynb: %d cells' % len(nb.cells))
