// 取扱説明書の画面を撮る（Playwright）。
//   node tools/manual_shots.js <シミュレータの html> <基盤地図情報 DEM5A の ZIP> [出力先]
// 例：node tools/manual_shots.js dist/sabo2d.html FG-GML-553706-DEM5A-20251113.zip docs/manual/img
const { chromium } = require('playwright');
const path = require('path'), fs = require('fs');
const HTML = path.resolve(process.argv[2] || 'dist/sabo2d.html'), ZIP = process.argv[3], OUT = path.resolve(process.argv[4] || 'docs/manual/img');
const ONLY = process.env.ONLY ? process.env.ONLY.split(',') : null;
const META = { slidePoly: [[538.6, 1215.4], [544.2, 1217.5], [550.4, 1214.6], [553.1, 1209.2], [551.3, 1203.8], [542.2, 1205.9], [537.8, 1210.7]], outlet: [560, 784] };   // 1.64 km² の小流域（平滑化 10 m・欠測の穴埋めのあと）
(async () => {
  const b = await chromium.launch({ args: ['--use-gl=swiftshader', '--enable-webgl', '--ignore-gpu-blocklist'] });
  const pg = await b.newPage({ viewport: { width: 1600, height: 1000 } });
  // 設定欄はタブに分かれているので、操作する前にその要素のタブを開く。「適用して最初から」は隠しボタン
  for (const m of ['fill', 'click', 'selectOption', 'check', 'uncheck', 'setInputFiles', 'dispatchEvent']) {
    const f = pg[m].bind(pg);
    pg[m] = async (s, ...a) => {
      if (typeof s === 'string') await pg.evaluate((q) => { try { const e = document.querySelector(q); if (e && window.SaboTab) window.SaboTab(e); } catch (er) {} }, s).catch(() => {});
      if (m === 'click' && s === '#apply') return pg.evaluate(() => document.getElementById('apply').click());
      return f(s, ...a);
    };
  }
  const errs = []; pg.on('pageerror', e => errs.push(e.message));
  await pg.addInitScript(() => { window.__sims = []; const iv = setInterval(() => { const S = window.Sabo2D; if (S && S.Sim && !S.Sim.__w) { const O = S.Sim; const W = function (P) { const s = new O(P); window.__sims.push(s); return s; }; W.prototype = O.prototype; W.__w = 1; S.Sim = W; clearInterval(iv); } }, 5); });
  const log = (...a) => console.log(new Date().toISOString().slice(11, 19), ...a);
  const W = (ms) => pg.waitForTimeout(ms);
  const shot = async (sel, name) => { const e = typeof sel === 'string' ? await pg.$(sel) : sel; if (!e) { log('見つからない', sel); return; } await e.scrollIntoViewIfNeeded(); await W(250); await e.screenshot({ path: path.join(OUT, name) }); log('撮影', name); };
  const ZO = +(process.env.ZO || 3); const shot3d = async (name) => { const dem = await pg.evaluate(() => !!(window.SaboDEMUI && window.SaboDEMUI.state.A && document.getElementById('g_src').value === 'dem')); if (dem) { if (process.env.VEX) await sel('v_vex', process.env.VEX); await pg.click('#v_home'); await W(300); if (process.env.V) { await sel('v_view', process.env.V); await W(600); } for (let k = 0; k < ZO; k++) { await pg.click('#v_zo'); await W(150); } await W(500); } await shot(await card('c3d'), name); };
  const card = async (id) => (await pg.evaluateHandle((id) => document.getElementById(id).closest('.card') || document.getElementById(id), id)).asElement();
  const det = async (title, open = true) => pg.evaluate(([t, o]) => { for (const d of document.querySelectorAll('details')) { const s = d.querySelector('summary'); if (s && s.textContent.trim().startsWith(t)) { d.open = o; if (o && window.SaboTab) window.SaboTab(d); if (o) { let p = d.parentElement; while (p) { if (p.tagName === 'DETAILS') p.open = true; p = p.parentElement; } } return true; } } return false; }, [title, open]);
  const detEl = async (title) => (await pg.evaluateHandle((t) => { for (const d of document.querySelectorAll('details')) { const s = d.querySelector('summary'); if (s && s.textContent.trim().startsWith(t)) return d; } return null; }, title)).asElement();
  const closeAll = () => pg.evaluate(() => { for (const d of document.querySelectorAll('#side details')) d.open = false; });
  const sel = (id, v) => pg.evaluate(([id, v]) => { const e = document.getElementById(id); e.value = v; e.dispatchEvent(new Event('change')); e.dispatchEvent(new Event('input')); }, [id, v]);
  // 時刻 T [s] まで計算する。描画を止めて計算だけ進め（速い）、最後に一度描き直す
  const runTo = async (T, maxMs = 1800000) => {
    const t0 = Date.now();
    for (;;) {
      const r = await pg.evaluate((T) => { const s = window.__sims[window.__sims.length - 1]; const w = performance.now(); while (performance.now() - w < 2000 && s.t < T && s.phase !== 'done') s.step(); return [s.t, s.phase]; }, T);
      if (r[0] >= T || r[1] === 'done' || Date.now() - t0 > maxMs) break;
    }
    await pg.click('#run'); await W(400); const ph = await pg.evaluate(() => window.__sims[window.__sims.length - 1].phase); if (ph !== 'done') await pg.click('#run');
    await W(2500);
  };
  const want = (k) => !ONLY || ONLY.includes(k);
  // DEM のセル → 地図の画面座標（必要なら表示を寄せる）
  const demXY = (x, y, sc) => pg.evaluate(([x, y, sc]) => {
    const U = window.SaboDEMUI, S = U.state; if (sc) { S.cam = { cx: x + 0.5, cy: y + 0.5, sc: sc }; U.draw(); }
    const v = S.view, r = document.getElementById('demmap').getBoundingClientRect(); return [r.left + (x + 0.5 - v.x0) * v.sc, r.top + (y + 0.5 - v.y0) * v.sc];
  }, [x, y, sc || 0]);
  const mode = (m) => sel('dm_mode', m);

  await pg.goto('file://' + HTML); await W(2500);

  // ================================================================ A. 標準の画面（数値で作る河道）
  if (want('A')) try {
    await shot('body', 'a00_overview.png');
    await closeAll();
    for (const [t, f] of [['ハイドロ', 'a01_hydro.png'], ['河道地形', 'a02_terrain.png'], ['堰堤の諸元', 'a03_dams.png'], ['透過型', 'a04_open.png'], ['流砂モデル', 'a06_sedmodel.png'], ['計算条件', 'a07_calc.png']]) { await det(t); await shot(await detEl(t), f); await det(t, false); }
    // 粒径：ひな形・材料 B・区間
    await det('粒径と河床材料');
    await sel('gr_pre', 'sg'); await W(200);
    await pg.click('#gr_madd'); await sel('gr_to', '1'); await pg.fill('#gr_dm', '60'); await pg.fill('#gr_sg', '2'); await pg.click('#gr_gen'); await W(200);
    await pg.click('#bz_add'); await pg.evaluate(() => { const i = document.querySelectorAll('#bztab tbody tr')[0].querySelectorAll('input'); i[0].value = 120; i[1].value = 220; });
    await shot(await detEl('粒径と河床材料'), 'a05_grain.png');
    await det('粒径と河床材料', false);
    // 計算（既定の条件、20 分）
    await pg.click('#apply'); await W(1500);
    await runTo(1200);
    await shot(await card('c3d'), 'a10_3d.png');
    await sel('v_mode', 'dz'); await W(800); await shot(await card('c3d'), 'a11_3d_dz.png');
    await sel('v_mode', 'dm'); await W(800); await shot(await card('c3d'), 'a12_3d_dm.png');
    await sel('v_mode', 'red'); await W(800); await shot(await card('c3d'), 'a13_3d_red.png');
    await sel('v_mode', 'real'); await W(500);
    await shot(await card('plan'), 'a14_plan.png');
    await shot(await card('xsec'), 'a15_xsec.png');
    await shot(await card('prof'), 'a16_prof.png');
    await shot(await card('tsq'), 'a17_tsq.png'); await shot(await card('tss'), 'a18_tss.png');
    await shot(await card('tiles'), 'a19_stats.png');
  } catch (e) { log('A でエラー', e.message); }
  if (want('A2')) try {
    // 土石流の流入（波形：土石流の一波、流砂モデル：土石流）
    await det('ハイドロ'); await sel('h_mode', 'debris'); await W(300); await shot(await detEl('ハイドロ'), 'a20_hydro_debris.png'); await det('ハイドロ', false);
    await pg.click('#apply'); await W(1500); await runTo(600);
    await shot(await card('c3d'), 'a21_debris_3d.png'); await shot(await card('prof'), 'a22_debris_prof.png'); await shot(await card('tss'), 'a23_debris_tss.png');
    await sel('h_mode', 'tri');
  } catch (e) { log('A2 でエラー', e.message); }

  // ================================================================ B. DEM（基盤地図情報 DEM5A）
  if (ZIP && want('B')) try {
    await sel('g_src', 'dem'); await W(500); await closeAll(); await det('河道地形');
    await shot(await detEl('河道地形'), 'b00_dem_panel.png');
    await pg.setInputFiles('#dem_file', [ZIP]);
    await pg.waitForFunction(() => !/読み込み中/.test(document.getElementById('dem_info').textContent) && window.SaboDEMUI.state.A, null, { timeout: 120000, polling: 300 });
    await W(800);
    await pg.click('#dm_zall'); await W(600);
    await shot('#demcard', 'b01_dem_all.png');
    await sel('dm_shade', 'red'); await W(800); await shot('#demcard', 'b02_dem_red.png');
    await sel('dm_shade', 'slope'); await W(800); await shot('#demcard', 'b03_dem_slope.png');
    await sel('dm_shade', 'hs'); await W(500);
    // 出口を選ぶ（谷筋のセルをクリック）
    await mode('outlet'); await W(200);
    await (await pg.$('#demmap')).scrollIntoViewIfNeeded(); let p = await demXY(META.outlet[0], META.outlet[1], 3); await pg.mouse.click(p[0], p[1]); await W(1500);
    await pg.click('#dm_zbas'); await W(800);
    await shot('#demcard', 'b04_dem_basin.png');
    // 供給点とハイドロ：土石流の一波（ピーク 20 m³/s、濃度 0.3）
    await det('土石流の発生形態'); await det('谷頭の供給点');
    await pg.evaluate(() => { document.querySelectorAll('#srctab tbody tr').forEach(tr => { const s = tr.querySelector('select'); s.value = 'debris'; s.dispatchEvent(new Event('change')); const set = (k, v) => { const i = tr.querySelector('input[data-k="' + k + '"]'); i.value = v; i.dispatchEvent(new Event('change')); }; set('Qb', 0.5); set('Qp', 20); set('t0', 30); set('tr', 60); set('tf', 300); set('C', 0.3); }); });
    await shot(await detEl('谷頭の供給点'), 'b05_sources.png');
    await pg.click('#apply'); await W(2500);
    await runTo(900);
    await shot3d('b10_debris_3d.png');
    await sel('v_mode', 'depth'); await W(800); await shot3d('b11_debris_depth.png'); await sel('v_mode', 'real');
    await shot(await card('plan'), 'b12_debris_plan.png');
    await shot(await card('prof'), 'b13_debris_prof.png');
    await shot(await card('xsec'), 'b14_debris_xsec.png');
    await shot(await card('tsq'), 'b15_debris_tsq.png'); await shot(await card('tss'), 'b16_debris_tss.png');
    // 支川の線へ切り替え
    await sel('l_path', '1'); await W(800); await shot(await card('prof'), 'b17_prof_trib.png'); await sel('l_path', '0');
    await det('谷頭の供給点', false);
  } catch (e) { log('B でエラー', e.message); }

  // ---- 雨：一様のハイエトグラフ（集中化）と疑似 XRAIN（格子に降らせる）
  if (ZIP && want('R')) try {
    await det('土石流の発生形態'); await sel('gen_type', 'bed'); await W(300);
    await det('雨（谷へ集まる水）');
    await pg.evaluate(() => { const c = document.getElementById('rn_on'); if (!c.checked) c.click(); }); await W(300);
    await pg.click('#rn_calc'); await W(2500);
    await shot(await detEl('雨（谷へ集まる水）'), 'b20_rain_panel.png');
    await shot('#rn_chartw', 'b21_rain_chart.png');
    // 疑似 XRAIN
    await sel('rn_mode', 'grid'); await W(300);
    await pg.setInputFiles('#rn_file', [path.resolve(__dirname, '../examples/pseudo_xrain_553706.csv')]); await W(4000);
    await pg.click('#rn_calc'); await W(3000);   // 集め方が「斜面で流出計算」のうちに図を出す
    await shot('#rn_chartw', 'b23_xrain_chart.png');
    await sel('rn_meth', 'grid'); await W(300);
    await shot(await detEl('雨（谷へ集まる水）'), 'b22_xrain_panel.png');
    await shot('#demcard', 'b24_xrain_map.png');
    await pg.click('#apply'); await W(2500);
    await runTo(1800);
    await shot(await card('plan'), 'b25_xrain_plan.png'); await shot(await card('tsq'), 'b26_xrain_tsq.png'); await shot3d('b27_xrain_3d.png');
    await pg.evaluate(() => { const c = document.getElementById('rn_on'); if (c.checked) c.click(); }); await det('雨（谷へ集まる水）', false);
  } catch (e) { log('R でエラー', e.message); }

  // ---- 斜面崩壊（多角形・表層・流動化 → 土石流化）と深層崩壊（土塊）
  if (ZIP && want('S')) try {
    await det('土石流の発生形態'); await sel('gen_type', 'slide'); await W(300);
    await mode('slide'); await sel('dm_shape', 'poly'); await W(200);
    // 1 番目の供給点の近くの急な斜面に、流下方向へ伸びた 7 角形
    const P0 = await pg.evaluate(() => {
      const S = window.SaboDEMUI.state, D = S.D, A = S.A, M = window.SaboDEM, h0 = S.src[0].c, hx = h0 % D.nx, hy = (h0 / D.nx) | 0; let best = -1, bs = 0;
      for (let dy = -60; dy <= 60; dy++) for (let dx = -60; dx <= 60; dx++) { const x = hx + dx, y = hy + dy, c = y * D.nx + x; if (!S.mask[c] || A.stream[c]) continue; const sl = M.localSlope(D, A, c, 2), d = Math.hypot(dx, dy); if (d > 15 && d < 60 && sl > bs && sl < 1.2) { bs = sl; best = c; } }
      let cc = best; for (let k = 0; k < 10; k++) { const d = M.downOf(D, A.dir, cc); if (d < 0) break; cc = d; }
      const bx = best % D.nx, by = (best / D.nx) | 0; let ax = cc % D.nx - bx, ay = ((cc / D.nx) | 0) - by; const L = Math.hypot(ax, ay) || 1; ax /= L; ay /= L; const px = -ay, py = ax;
      return [[-4, 0], [-1, 5], [5, 6], [11, 3], [12, -2], [4, -6], [-2, -4]].map(([a, b]) => [bx + a * ax + b * px, by + a * ay + b * py]);
    });
    let q = await demXY(P0[0][0], P0[0][1], 6);
    await (await pg.$('#demmap')).scrollIntoViewIfNeeded();
    for (let k = 0; k < P0.length; k++) { q = await demXY(P0[k][0], P0[k][1]); if (k < P0.length - 1) await pg.mouse.click(q[0], q[1]); else await pg.mouse.dblclick(q[0], q[1]); await W(120); }
    await W(600);
    await det('斜面崩壊'); await pg.evaluate(() => { const r = document.querySelector('#sltab tbody tr'); if (!r) return; const s = r.querySelectorAll('select'); s.forEach(x => { if ([...x.options].some(o => o.value === 'time')) { x.value = 'time'; x.dispatchEvent(new Event('change')); } if ([...x.options].some(o => o.value === 'fluid')) { x.value = 'fluid'; x.dispatchEvent(new Event('change')); } }); const t = r.querySelector('input[data-k="t"]'); if (t) { t.value = 60; t.dispatchEvent(new Event('change')); } });
    await shot('#demcard', 'b30_slide_map.png');
    await shot(await detEl('斜面崩壊'), 'b31_slide_panel.png');
    await pg.click('#apply'); await W(2500); await runTo(600);
    await shot3d('b32_slide_3d.png'); await sel('v_mode', 'dz'); await W(800); await shot3d('b33_slide_dz.png'); await sel('v_mode', 'real');
    await shot(await card('plan'), 'b34_slide_plan.png');
    // 深層崩壊・土塊
    await pg.evaluate(() => { const r = document.querySelector('#sltab tbody tr'); r.querySelectorAll('select').forEach(x => { if ([...x.options].some(o => o.value === 'deep')) { x.value = 'deep'; x.dispatchEvent(new Event('change')); } }); });
    await pg.evaluate(() => { const r = document.querySelector('#sltab tbody tr'); r.querySelectorAll('select').forEach(x => { if ([...x.options].some(o => o.value === 'block')) { x.value = 'block'; x.dispatchEvent(new Event('change')); } }); const d = r.querySelector('input[data-k="D"]'); if (d) { d.value = 8; d.dispatchEvent(new Event('change')); } });
    await pg.click('#apply'); await W(2500); await runTo(300);
    await sel('v_mode', 'dz'); await W(800); await shot3d('b35_deep_dz.png'); await sel('v_mode', 'real');
    await pg.click('#sl_clr').catch(() => {}); await det('斜面崩壊', false);
  } catch (e) { log('S でエラー', e.message); }

  // ---- 天然ダム（例を置く → 越流、透水係数を大きく → 進行性）
  if (ZIP && want('L')) try {
    await det('土石流の発生形態'); await sel('gen_type', 'ldam'); await W(300);
    await det('天然ダム'); await pg.click('#ld_ex'); await W(1500);
    await shot('#demcard', 'b40_ldam_map.png'); await shot(await detEl('天然ダム'), 'b41_ldam_panel.png');
    await pg.click('#apply'); await W(2500); await runTo(1200);
    await shot3d('b42_ldam_3d.png'); await shot(await card('prof'), 'b43_ldam_prof.png'); await shot(await card('tsq'), 'b44_ldam_tsq.png');
    await pg.fill('#ld_k', '0.03'); await pg.dispatchEvent('#ld_k', 'change'); await pg.click('#apply'); await W(2500); await runTo(1200);
    await shot(await card('prof'), 'b45_ldam_prog_prof.png'); await shot(await card('tsq'), 'b46_ldam_prog_tsq.png');
    await pg.fill('#ld_k', '0'); await pg.dispatchEvent('#ld_k', 'change'); await pg.click('#ld_clr'); await det('天然ダム', false);
  } catch (e) { log('L でエラー', e.message); }

  // ---- DEM の砂防堰堤と河床材料の範囲
  if (ZIP && want('D')) try {
    await sel('gen_type', 'bed'); await W(200);
    const cells = await pg.evaluate(() => { const S = window.SaboDEMUI.state, D = S.D, A = S.A; let c = S.outlet, p = [c]; for (let k = 0; k < 2000; k++) { const x = c % D.nx, y = (c / D.nx) | 0; let bu = -1, ba = 0; for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) { const o = (y + dy) * D.nx + x + dx; if ((dx || dy) && o >= 0 && S.mask[o] && window.SaboDEM.downOf(D, A.dir, o) === c && A.acc[o] > ba) { ba = A.acc[o]; bu = o; } } if (bu < 0) break; c = bu; p.push(c); } return [0.2, 0.35].map(f => { const q = p[Math.floor(p.length * f)]; return [q % D.nx, (q / D.nx) | 0]; }); });
    await mode('dam');
    await (await pg.$('#demmap')).scrollIntoViewIfNeeded();
    for (const [x, y] of cells) { const q = await demXY(x, y, 2.5); await pg.mouse.click(q[0], q[1]); await W(300); }
    await det('砂防堰堤（DEM）'); await pg.evaluate(() => { const r = document.querySelectorAll('#dmtab tbody tr'); if (r[1]) { const s = r[1].querySelector('select'); s.value = 'slit'; s.dispatchEvent(new Event('change')); } });
    await pg.click('#dm_zbas'); await W(500);
    await shot('#demcard', 'b50_dam_map.png'); await shot(await detEl('砂防堰堤（DEM）'), 'b51_dam_panel.png');
    await pg.click('#apply'); await W(2500); await runTo(900);
    await shot(await card('prof'), 'b52_dam_prof.png'); await shot3d('b53_dam_3d.png'); await shot(await card('damres'), 'b54_dam_table.png');
    await pg.click('#dmt_clr');
    // 河床材料の範囲（多角形）
    await det('粒径と河床材料');
    await mode('mat'); await W(200);
    const c0 = cells[0]; const sq = [[-25, -25], [25, -25], [25, 25], [-25, 25]];
    await (await pg.$('#demmap')).scrollIntoViewIfNeeded();
    for (let k = 0; k < 4; k++) { const q = await demXY(c0[0] + sq[k][0], c0[1] + sq[k][1], k ? 0 : 1.5); if (k < 3) await pg.mouse.click(q[0], q[1]); else await pg.mouse.dblclick(q[0], q[1]); await W(150); }
    await W(500); await shot('#demcard', 'b55_mat_map.png');
    await pg.click('#apply'); await W(2500);
    await sel('v_mode', 'dm'); await W(800); await shot3d('b56_mat_3d.png'); await sel('v_mode', 'real');
  } catch (e) { log('D でエラー', e.message); }

  // ---- 疑似 XRAIN を長め（90 分）に。DEM と流域の選び直しから
  if (ZIP && want('X')) try {
    if (!(await pg.evaluate(() => window.SaboDEMUI && window.SaboDEMUI.state.A))) {
      await sel('g_src', 'dem'); await W(500); await det('河道地形');
      await pg.setInputFiles('#dem_file', [ZIP]);
      await pg.waitForFunction(() => !/読み込み中/.test(document.getElementById('dem_info').textContent) && window.SaboDEMUI.state.A, null, { timeout: 120000, polling: 300 });
      await mode('outlet'); await W(200); await (await pg.$('#demmap')).scrollIntoViewIfNeeded();
      const p = await demXY(META.outlet[0], META.outlet[1], 3); await pg.mouse.click(p[0], p[1]); await W(1500);
    }
    await pg.evaluate(() => { const S = window.SaboDEMUI.state; S.src = []; S.slides = []; S.ldams = []; S.dams = []; S.matz = []; window.SaboDEMUI.draw(); });
    await det('河道地形'); await pg.fill('#dem_dur', '90'); await pg.dispatchEvent('#dem_dur', 'change');
    await det('土石流の発生形態'); await det('雨（谷へ集まる水）');
    await pg.evaluate(() => { const c = document.getElementById('rn_on'); if (!c.checked) c.click(); }); await W(300);
    await sel('rn_mode', 'grid'); await W(300);
    await pg.setInputFiles('#rn_file', [path.resolve(__dirname, '../examples/pseudo_xrain_553706.csv')]); await W(4000);
    await sel('rn_meth', 'lumped'); await W(300); await pg.click('#rn_calc'); await W(3000);
    await shot('#rn_chartw', 'b23_xrain_chart.png');
    await sel('rn_meth', 'grid'); await W(300);
    await pg.click('#apply'); await W(2500);
    await runTo(3600);
    await sel('v_mode', 'depth'); await W(800); await shot3d('b27_xrain_3d.png'); await sel('v_mode', 'real');
    await shot(await card('plan'), 'b25_xrain_plan.png'); await shot(await card('tsq'), 'b26_xrain_tsq.png');
    await pg.evaluate(() => { const c = document.getElementById('rn_on'); if (c.checked) c.click(); }); await pg.fill('#dem_dur', '30');
  } catch (e) { log('X でエラー', e.message); }

  log('ページのエラー', errs.length ? errs : 'なし');
  await b.close();
})();
