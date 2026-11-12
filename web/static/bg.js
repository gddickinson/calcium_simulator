// Subtle animated background: a few slow, fading Ca2+ "puffs" (expanding rings). Off under prefers-reduced-motion.
export function startBackground(canvas) {
  if (!canvas || window.matchMedia('(prefers-reduced-motion: reduce)').matches) return;
  const ctx = canvas.getContext('2d');
  const puffs = [];
  let w = 0, h = 0, last = 0;
  const resize = () => { w = canvas.width = innerWidth; h = canvas.height = innerHeight; };
  addEventListener('resize', resize); resize();
  const frame = t => {
    if (t - last > 50) {
      last = t;
      if (puffs.length < 6 && Math.random() < 0.04) puffs.push({ x: Math.random() * w, y: Math.random() * h, r: 2, a: 0.8 });
      ctx.clearRect(0, 0, w, h);
      for (const p of puffs) {
        p.r += 1.2; p.a *= 0.975;
        const g = ctx.createRadialGradient(p.x, p.y, 0, p.x, p.y, p.r);
        g.addColorStop(0, `rgba(255,170,40,${p.a})`); g.addColorStop(1, 'rgba(255,60,20,0)');
        ctx.fillStyle = g; ctx.beginPath(); ctx.arc(p.x, p.y, p.r, 0, 7); ctx.fill();
      }
      for (let i = puffs.length - 1; i >= 0; i--) if (puffs[i].a < 0.02) puffs.splice(i, 1);
    }
    requestAnimationFrame(frame);
  };
  requestAnimationFrame(frame);
}
