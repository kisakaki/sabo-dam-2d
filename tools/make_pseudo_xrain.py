"""疑似 XRAIN（格子の雨）の CSV を作る。

本物の XRAIN（国土交通省の XバンドMPレーダ雨量計）は 250 m メッシュ・1 分間隔の雨量強度 [mm/h] です。
ここでは、その形をまねた架空の雨を作ります（観測値ではありません）。

雨の形：一様な層状の雨（背景）＋ 北東へ動く対流セル（ガウス分布、強さは時間で山形）。
出力：1 行に「時刻 [分], 経度, 緯度, 雨量強度 [mm/h]」。シミュレータの「格子の雨（XRAIN などを変換した CSV）」で読めます。

使い方：
  python tools/make_pseudo_xrain.py                      # 基盤地図情報 5337-06（553706）の範囲、120 分、5 分おき
  python tools/make_pseudo_xrain.py --minutes 180 --step 1 --peak 100 --out rain.csv
"""
import argparse, math

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--lat0', type=float, default=36.6667, help='南端の緯度（既定：2 次メッシュ 553706）')
ap.add_argument('--lat1', type=float, default=36.75, help='北端の緯度')
ap.add_argument('--lon0', type=float, default=137.75, help='西端の経度')
ap.add_argument('--lon1', type=float, default=137.875, help='東端の経度')
ap.add_argument('--mesh', type=float, default=250.0, help='メッシュの大きさ [m]（XRAIN は 250 m）')
ap.add_argument('--minutes', type=float, default=120.0, help='期間 [分]')
ap.add_argument('--step', type=float, default=5.0, help='時間の間隔 [分]（XRAIN は 1 分。ファイルが大きくなるので既定 5 分）')
ap.add_argument('--base', type=float, default=8.0, help='背景の雨の強さ [mm/h]')
ap.add_argument('--peak', type=float, default=80.0, help='対流セルの中心の最大の強さ [mm/h]')
ap.add_argument('--radius', type=float, default=2500.0, help='対流セルの大きさ（ガウス分布の標準偏差）[m]')
ap.add_argument('--speed', type=float, default=15.0, help='対流セルの移動速度 [km/h]（北東へ）')
ap.add_argument('--tpeak', type=float, default=50.0, help='対流セルが最も強くなる時刻 [分]')
ap.add_argument('--out', default='examples/pseudo_xrain_553706.csv')
a = ap.parse_args()

latm = math.radians((a.lat0 + a.lat1) / 2)
my = 111132.92 - 559.82 * math.cos(2 * latm) + 1.175 * math.cos(4 * latm)    # 緯度 1 度の長さ [m]
mx = 111412.84 * math.cos(latm) - 93.5 * math.cos(3 * latm)                 # 経度 1 度の長さ [m]
dlat, dlon = a.mesh / my, a.mesh / mx
W, H = (a.lon1 - a.lon0) * mx, (a.lat1 - a.lat0) * my
nx, ny = int(round(W / a.mesh)), int(round(H / a.mesh))
# 対流セル：南西の外から入り、北東へ。中心が範囲の中央を tpeak に通る
v = a.speed * 1000 / 60                                     # m/分
ux = uy = math.sqrt(0.5)
n = 0
with open(a.out, 'w', encoding='utf-8') as f:
    f.write('# 疑似 XRAIN（架空の雨。観測値ではない）。%g m メッシュ、%g 分おき。列：時刻[分],経度,緯度,雨量強度[mm/h]\n' % (a.mesh, a.step))
    f.write('time_min,lon,lat,mm_h\n')
    t = 0.0
    while t <= a.minutes + 1e-9:
        cx = W / 2 + ux * v * (t - a.tpeak); cy = H / 2 + uy * v * (t - a.tpeak)
        amp = a.peak * math.exp(-((t - a.tpeak) / (0.45 * a.minutes)) ** 2 * 4)   # 時間の山形
        bg = a.base * (1 if 10 <= t <= a.minutes - 10 else 0.3)
        for j in range(ny):
            for i in range(nx):
                x, y = (i + 0.5) * a.mesh, (j + 0.5) * a.mesh
                r = bg + amp * math.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2 * a.radius ** 2))
                f.write('%g,%.6f,%.6f,%.1f\n' % (t, a.lon0 + (i + 0.5) * dlon, a.lat0 + (j + 0.5) * dlat, r))
                n += 1
        t += a.step
print('%s：%d 行（%d × %d メッシュ × %d 時刻）' % (a.out, n, nx, ny, int(a.minutes / a.step) + 1))
