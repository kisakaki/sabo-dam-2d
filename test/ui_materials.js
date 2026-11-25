const { chromium } = require('playwright');
(async () => {
  const b = await chromium.launch(); const pg = await b.newPage({ viewport: { width: 1500, height: 1100 } });
  const errs = []; pg.on('pageerror', e => errs.push(e.message + ' ' + (e.stack || '').split('\n')[1]));
  await pg.addInitScript(() => { window.__sims = []; const iv = setInterval(() => { const S = window.Sabo2D; if (S && S.Sim && !S.Sim.__w) { const O = S.Sim; const W = function (P) { const s = new O(P); window.__sims.push(s); return s; }; W.prototype = O.prototype; W.__w = 1; S.Sim = W; clearInterval(iv); } }, 5); });
  await pg.goto('file://' + require('path').resolve(__dirname, '../dist/' + (require('fs').existsSync(require('path').resolve(__dirname, '../dist/sabo2d.html')) ? 'sabo2d.html' : 'sabo2d-cdn.html')) + ''); await pg.waitForTimeout(1500);
  await pg.evaluate(() => { document.getElementById('grtab').closest('details').open = true; });
  console.log('start cols', await pg.evaluate(() => document.querySelectorAll('#grtab thead th').length));
  // 粒径階 12
  await pg.selectOption('#gr_to', 'A'); await pg.fill('#gr_dm', '20'); await pg.fill('#gr_sg', '4'); await pg.selectOption('#gr_k', '12'); await pg.click('#gr_gen');
  // 材料を 4 つ（A, B, C, D）に
  for (let i = 0; i < 3; i++) await pg.click('#gr_madd');
  for (const [m, dm] of [['1', '5'], ['2', '60'], ['3', '150']]) { await pg.selectOption('#gr_to', m); await pg.fill('#gr_dm', dm); await pg.fill('#gr_sg', '2'); await pg.click('#gr_gen'); }
  console.log('cols', await pg.evaluate(() => Array.from(document.querySelectorAll('#grtab thead th')).map(t => t.textContent).join('|')), 'rows', await pg.evaluate(() => document.querySelectorAll('#grtab tbody tr').length));
  console.log(await pg.evaluate(() => document.getElementById('gr_sum').textContent.slice(0, 400)));
  await pg.click('#bz_add'); await pg.click('#bz_add');
  await pg.evaluate(() => { const r = document.querySelectorAll('#bztab tbody tr'); r[1].querySelector('select').value = '3'; const i = r[1].querySelectorAll('input'); i[0].value = 200; i[1].value = 300; });
  const d = await pg.evaluateHandle(() => document.getElementById('grtab').closest('details')); await d.asElement().screenshot({ path: '/tmp/m12.png' });
  await pg.click('#apply'); await pg.waitForTimeout(1200);
  console.log('normal', await pg.evaluate(() => { const s = window.__sims[window.__sims.length - 1]; const n = [0, 0, 0, 0]; for (let c = 0; c < s.N; c++) n[s.matIdx[c]]++; return JSON.stringify({ K: s.K, mats: n, nMat: s.mats.length }); }));
  await pg.click('#run'); await pg.waitForTimeout(3000); await pg.click('#run');
  // DEM：材料 D の多角形
  await pg.evaluate(() => { const g = document.getElementById('g_src'); g.value = 'dem'; g.dispatchEvent(new Event('change')); });
  await pg.click('#dem_sample'); await pg.waitForTimeout(1200);
  await pg.evaluate(() => { const m = document.getElementById('dm_mode'); m.value = 'mat'; m.dispatchEvent(new Event('change')); });
  console.log('dm_matsel', await pg.evaluate(() => Array.from(document.querySelectorAll('#dm_matsel option')).map(o => o.textContent).join(',')));
  await pg.selectOption('#dm_matsel', '3');
  const cv = await pg.$('#demmap'); await cv.scrollIntoViewIfNeeded(); const bb = await cv.boundingBox();
  const pts = [[0.35, 0.35], [0.6, 0.35], [0.6, 0.6], [0.35, 0.6]].map(p => [bb.x + bb.width * p[0], bb.y + bb.height * p[1]]);
  for (const p of pts.slice(0, 3)) { await pg.mouse.click(p[0], p[1]); await pg.waitForTimeout(80); }
  await pg.mouse.dblclick(pts[3][0], pts[3][1]); await pg.waitForTimeout(300);
  await pg.click('#apply'); await pg.waitForTimeout(1500);
  console.log('dem', await pg.evaluate(() => { const s = window.__sims[window.__sims.length - 1]; const n = [0, 0, 0, 0]; for (let c = 0; c < s.N; c++) n[s.matIdx[c]]++; return JSON.stringify({ K: s.K, mats: n }); }));
  await pg.click('#run'); await pg.waitForTimeout(3000); await pg.click('#run');
  console.log('errors', errs); await b.close();
})();
