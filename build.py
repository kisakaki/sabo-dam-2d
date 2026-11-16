# -*- coding: utf-8 -*-
"""src/ の部品から公開用の HTML を作る。
  index.html        … three.js を中に含む 1 ファイル版（GitHub Pages とオフライン用）
  dist/sabo2d-cdn.html … three.js を CDN から読む軽い版
  dist/sabo2d-solver.js … 計算部分だけ（Node.js で画面なしに計算するとき）
three.js は `npm install` で node_modules に入れておく。"""
import json
import os

page = open('src/page.html', encoding='utf-8').read()
sol = open('src/sabo2d-solver.js', encoding='utf-8').read()
ui = open('src/ui.js', encoding='utf-8').read()
sv = open('src/survey.js', encoding='utf-8').read()
txs = open('tools/test_cross_sections.csv', encoding='utf-8-sig').read()
tlo = open('tools/test_longitudinal.csv', encoding='utf-8-sig').read()
ui = ui.replace('/*TESTCSV*/', 'var TEST_XS = ' + json.dumps(txs, ensure_ascii=False) + ', TEST_LONG = ' + json.dumps(tlo, ensure_ascii=False) + ';')
ui = sv + '\n' + ui
body = page.replace('/*SOLVER*/', sol).replace('/*UI*/', ui)
os.makedirs('dist', exist_ok=True)
cdn = body.replace('<!--THREE-->', '<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>')
head = '''<!DOCTYPE html>
<html lang="ja">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="description" content="連続砂防堰堤（低ダム群・高ダム）の効果を、平面2次元の河床変動計算で調べるブラウザ用シミュレータ。混合粒径・網状流路・砂州を扱う。">
<meta property="og:type" content="website">
<meta property="og:url" content="https://kisakaki.github.io/sabo-dam-2d/">
<meta property="og:title" content="連続堰堤 河床変動シミュレータ">
<meta property="og:description" content="連続砂防堰堤の効果を、平面2次元の河床変動計算で調べるブラウザ用シミュレータ。">
<meta property="og:image" content="https://kisakaki.github.io/sabo-dam-2d/og.png">
<meta property="og:locale" content="ja_JP">
<meta name="twitter:card" content="summary_large_image">
<meta name="twitter:image" content="https://kisakaki.github.io/sabo-dam-2d/og.png">
</head>
<body style="margin:0">
'''
open('dist/sabo2d-cdn.html', 'w', encoding='utf-8').write(head + cdn + '\n</body>\n</html>\n')
three_path = 'node_modules/three/build/three.min.js'
if os.path.exists(three_path):
    three = open(three_path, encoding='utf-8').read()
    sa = body.replace('<!--THREE-->', '<script>\n' + three.replace('</script>', '<\\/script>') + '\n</script>')
    open('index.html', 'w', encoding='utf-8').write(head + sa + '\n</body>\n</html>\n')
    print('index.html', len(sa))
else:
    print('three.js がありません。`npm install` のあとでもう一度実行すると index.html を作ります。')
open('dist/sabo2d-solver.js', 'w', encoding='utf-8').write(sol)
print('dist/sabo2d-cdn.html', len(cdn))
