const { chromium } = require('playwright');
(async () => {
  const b = await chromium.launch(); const pg = await b.newPage({ viewport: { width: 1500, height: 1000 } });
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
  await pg.goto('file://' + require('path').resolve(__dirname, '../index.html')); await pg.waitForTimeout(1500);
  await pg.evaluate(() => { const g = document.getElementById('g_src'); g.value = 'dem'; g.dispatchEvent(new Event('change')); const t = document.getElementById('gen_type'); t.value = 'slide'; t.dispatchEvent(new Event('change')); });
  await pg.click('#dem_sample'); await pg.waitForTimeout(1500);
  await pg.evaluate(() => { const s = document.getElementById('dm_shape'); s.value = 'poly'; s.dispatchEvent(new Event('change')); for (const d of document.querySelectorAll('#dem_more details')) d.open = true; });
  const cv = await pg.$('#demmap'); await cv.scrollIntoViewIfNeeded();
  // 選んだ流域の中で谷筋でないセルを中心に 5 角形の頂点を作る
  const pts = await pg.evaluate(() => {
    const S = window.SaboDEMUI.state, D = S.D, v = S.view, r = document.getElementById('demmap').getBoundingClientRect();
    let c0 = -1; for (let c = 0; c < S.mask.length; c++) if (S.mask[c] && !S.A.stream[c] && (c % D.nx) > 50 && (c % D.nx) < 110) { const y = (c / D.nx) | 0; if (y > 40 && y < 80) { c0 = c; break; } }
    const cx = c0 % D.nx + 6, cy = ((c0 / D.nx) | 0) + 6, out = [];
    for (let k = 0; k < 5; k++) { const a = k * 2 * Math.PI / 5, rr = k % 2 ? 4 : 7; out.push([r.left + (cx + rr * Math.cos(a) + 0.5 - v.x0) * v.sc, r.top + (cy + rr * Math.sin(a) + 0.5 - v.y0) * v.sc]); }
    return out;
  });
  for (const p of pts.slice(0, 4)) { await pg.mouse.click(p[0], p[1]); await pg.waitForTimeout(80); }
  await pg.mouse.dblclick(pts[4][0], pts[4][1]); await pg.waitForTimeout(300);
  console.log(await pg.evaluate(() => { const S = window.SaboDEMUI.state; return JSON.stringify({ n: S.slides.length, poly: S.slides[0] && S.slides[0].poly.length, R: S.slides[0] && S.slides[0].R, msg: document.getElementById('dem_sel').textContent, cell: document.querySelector('#sltab tbody td:nth-child(3)').textContent }); }));
  await (await pg.$('#demcard')).screenshot({ path: '/home/claude/work/gsi/ui_poly.png' });
  await pg.click('#apply'); await pg.waitForTimeout(1500);
  console.log(await pg.evaluate(() => document.getElementById('gridinfo').textContent.slice(-80)));
  await pg.click('#run'); await pg.waitForTimeout(6000); await pg.click('#run');
  console.log('errors', errs); await b.close();
})();
