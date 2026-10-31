/* Validación de la FASE 3: electrón, dinámica relativista y dirección arbitraria de B
 * Uso:  node test_fase3.js          (termina con código 1 si alguna comprobación falla)
 * R1  radio y frecuencia relativistas del electrón (V = 0)        R2  sentido de giro (e⁻ vs p⁺)
 * R3  límite clásico: a baja energía relativista ≈ clásico         R4  conservación de la energía
 * B1  B hacia −z = reflejo de B hacia +z                           B2  B ∥ v: sin fuerza
 * B3  B en +y: círculo en el plano xz                              B4  B inclinado: hélice (radio y velocidad paralela)
 * E1  electrón a 1 keV: queda atrapado en el hueco                 E2  electrón a 2-5 MeV: oscila, nunca se acelera
 * E3  fase: el electrón también gana energía con φ = π/2
 * P1  protón relativista: se atasca con f_E = f_c y sale con f_E ≈ 0.995 f_c                         */
const P = require('./cyclotron_physics.js');
const { theory, PARTICLES, E_CHARGE: e, C_LIGHT: c } = P;
const el = PARTICLES.electron, pr = PARTICLES.proton;
let fallos = 0;
const check = (nombre, ok, detalle) => { if (!ok) fallos++; console.log(`  [${ok ? 'OK   ' : 'FALLO'}] ${nombre}${detalle ? '  →  ' + detalle : ''}`); };
const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
const B = 1;

// Cicloide sin campo E: partícula libre en B (V = 0, D enormes para que nunca salga)
function libre(pt, K0, extra, pasos, spt) {
  const fc = theory.frequency(pt.q, pt.m, B);
  const sim = P.createCyclotron(Object.assign({ q: pt.q, m: pt.m, B, V: 0, gap: 0.01, R: 1e3, K0, y0: 0, dt: (1 / fc) / (spt || 400) }, extra));
  const tray = [sim.state()];
  for (let i = 0; i < pasos; i++) { sim.step(); tray.push(sim.state()); }
  return { sim, tray };
}

console.log('=== R1: electrón de 1 MeV en B = 1 T, dinámica relativista (V = 0) ===');
{
  const K0 = 1e6 * e, g = theory.gamma(el.m, K0), rTeo = theory.radiusRel(el.q, el.m, B, K0);
  const spt = 400, Tteo = g * theory.period(el.q, el.m, B);
  const fc = theory.frequency(el.q, el.m, B), n = Math.round(Tteo * fc * spt * 5);   // 5 vueltas
  const { tray } = libre(el, K0, { relativistic: true }, n, spt);
  const ys = tray.map((s) => s.y), diam = Math.max(...ys) - Math.min(...ys);
  check('radio = p/(eB)', rel(diam / 2, rTeo) < 2e-3, `simulado ${(diam / 2 * 1e3).toFixed(4)} mm, teoría ${(rTeo * 1e3).toFixed(4)} mm (clásico daría ${(theory.radius(el.q, el.m, B, c) * 1e3).toFixed(2)} mm con v = c)`);
  // período: ángulo acumulado alrededor del centro (0, +r) para q < 0
  let ang = 0, prev = Math.atan2(tray[0].y - rTeo, tray[0].x);
  for (let i = 1; i < tray.length; i++) { const a = Math.atan2(tray[i].y - rTeo, tray[i].x); let d = a - prev; if (d > Math.PI) d -= 2 * Math.PI; if (d < -Math.PI) d += 2 * Math.PI; ang += d; prev = a; }
  const omega = Math.abs(ang) / tray[tray.length - 1].t;
  check('ω = eB/(γm)', rel(omega, 2 * Math.PI / Tteo) < 2e-3, `γ = ${g.toFixed(4)}; ω simulada ${omega.toExponential(4)}, teoría ${(2 * Math.PI / Tteo).toExponential(4)} rad/s`);
  check('v < c siempre', tray.every((s) => s.beta < 1), `β = ${tray[0].beta.toFixed(5)}`);
}

console.log('\n=== R2: sentido de giro con B hacia +z (visto desde arriba) ===');
{
  const Lz = (s) => s.x * s.vy - s.y * s.vx;                       // >0 antihorario, <0 horario
  const cy = (pt, K0, r) => libre(pt, K0, { relativistic: true }, 80).tray;
  const t1 = libre(pr, 1e5 * e, {}, 80).tray, t2 = libre(el, 1e6 * e, { relativistic: true }, 80).tray;
  const giro = (t) => { let a = 0, p = Math.atan2(t[0].vy, t[0].vx); for (let i = 1; i < t.length; i++) { const q = Math.atan2(t[i].vy, t[i].vx); let d = q - p; if (d > Math.PI) d -= 2 * Math.PI; if (d < -Math.PI) d += 2 * Math.PI; a += d; p = q; } return a; };
  check('protón gira en sentido horario', giro(t1) < 0, `Δ(ángulo de v) = ${giro(t1).toFixed(3)} rad`);
  check('electrón gira en sentido antihorario', giro(t2) > 0, `Δ(ángulo de v) = ${giro(t2).toFixed(3)} rad`);
}

console.log('\n=== R3: límite clásico (protón de 1 keV, 3 vueltas): relativista ≈ clásico ===');
{
  const fc = theory.frequency(pr.q, pr.m, B), n = 3 * 400;
  const a = libre(pr, 1e3 * e, {}, n).tray, b = libre(pr, 1e3 * e, { relativistic: true }, n).tray;
  const r = theory.radius(pr.q, pr.m, B, theory.speedFromK(pr.m, 1e3 * e));
  let maxd = 0; for (let i = 0; i < a.length; i++) maxd = Math.max(maxd, Math.hypot(a[i].x - b[i].x, a[i].y - b[i].y));
  check('trayectorias casi idénticas', maxd / r < 1e-4, `diferencia máx = ${(maxd / r).toExponential(2)} × radio`);
}

console.log('\n=== R4: conservación de la energía con V = 0 (electrón 1 MeV, 20 vueltas, relativista) ===');
{
  const { tray } = libre(el, 1e6 * e, { relativistic: true }, 20 * 400);
  const K0 = tray[0].K, dev = Math.max(...tray.map((s) => Math.abs(s.K - K0) / K0));
  check('K constante', dev < 1e-9, `desviación máx ${dev.toExponential(2)}`);
}

console.log('\n=== B1: B hacia −z es el reflejo de B hacia +z (protón de 100 keV, V = 0) ===');
{
  const a = libre(pr, 1e5 * e, { Bdir: [0, 0, 1] }, 800).tray, b = libre(pr, 1e5 * e, { Bdir: [0, 0, -1] }, 800).tray;
  let d = 0; for (let i = 0; i < a.length; i++) d = Math.max(d, Math.abs(a[i].x - b[i].x), Math.abs(a[i].y + b[i].y));
  const r = theory.radius(pr.q, pr.m, B, theory.speedFromK(pr.m, 1e5 * e));
  check('x igual, y con signo opuesto', d / r < 1e-12, `diferencia ${(d / r).toExponential(2)} × radio`);
}

console.log('\n=== B2: B paralelo a la velocidad (B hacia +x): sin fuerza magnética ===');
{
  const { tray } = libre(pr, 1e5 * e, { Bdir: [1, 0, 0] }, 400);
  const u = tray[0], f = tray[tray.length - 1];
  check('velocidad constante', Math.abs(f.vx - u.vx) / u.vx < 1e-12 && Math.abs(f.vy) < 1e-6 && Math.abs(f.vz) < 1e-6, `vy = ${f.vy.toExponential(1)}, vz = ${f.vz.toExponential(1)}`);
  check('trayectoria recta', Math.abs(f.y) < 1e-12 && Math.abs(f.z) < 1e-12, `y = ${f.y.toExponential(1)}, z = ${f.z.toExponential(1)}`);
}

console.log('\n=== B3: B hacia +y, v hacia +x: círculo en el plano xz ===');
{
  const K0 = 1e5 * e, r = theory.radius(pr.q, pr.m, B, theory.speedFromK(pr.m, K0));
  const { tray } = libre(pr, K0, { Bdir: [0, 1, 0] }, 400);
  const zs = tray.map((s) => s.z), xs = tray.map((s) => s.x), ymax = Math.max(...tray.map((s) => Math.abs(s.y)));
  check('diámetro en z = 2r', rel(Math.max(...zs) - Math.min(...zs), 2 * r) < 2e-3, `${((Math.max(...zs) - Math.min(...zs)) * 1e3).toFixed(3)} mm vs ${(2 * r * 1e3).toFixed(3)} mm`);
  check('diámetro en x = 2r', rel(Math.max(...xs) - Math.min(...xs), 2 * r) < 2e-3);
  check('nunca sale del plano xz', ymax < 1e-12, `|y| máx = ${ymax.toExponential(1)}`);
}

console.log('\n=== B4: B inclinado θ = 30° desde +z (azimut 0): hélice ===');
{
  const th = Math.PI / 6, bdir = [Math.sin(th), 0, Math.cos(th)], K0 = 1e5 * e, v0 = theory.speedFromK(pr.m, K0);
  const vpar = v0 * Math.sin(th), vperp = v0 * Math.cos(th);             // v = v0·x̂, b̂ = (sinθ, 0, cosθ)
  const rTeo = pr.m * Math.sqrt(v0 * v0 - vpar * vpar) / (pr.q * B);
  const { sim, tray } = libre(pr, K0, { Bdir: bdir }, 3 * 400);
  const par = tray.map((s) => s.vx * bdir[0] + s.vy * bdir[1] + s.vz * bdir[2]);
  check('v∥ constante', Math.max(...par.map((v) => Math.abs(v - vpar))) / vpar < 1e-9, `v∥/v = ${(par[10] / v0).toFixed(5)} (teoría ${Math.sin(th).toFixed(5)})`);
  check('K constante', Math.max(...tray.map((s) => Math.abs(s.K - tray[0].K))) / tray[0].K < 1e-9);
  check('radio de la hélice = m·v⊥/(qB)', rel(tray[5].rCurv, rTeo) < 1e-6, `${(tray[5].rCurv * 1e3).toFixed(4)} mm vs ${(rTeo * 1e3).toFixed(4)} mm`);
  // avance a lo largo de B después de t: d∥ = v∥·t
  const f = tray[tray.length - 1], d = f.x * bdir[0] + f.y * bdir[1] + f.z * bdir[2];
  check('avance a lo largo de B = v∥·t', rel(d, vpar * f.t) < 1e-3, `${(d * 1e3).toFixed(3)} mm vs ${(vpar * f.t * 1e3).toFixed(3)} mm`);
}

console.log('\n=== E1: electrón de 1 keV en el ciclotrón (hueco de 1 cm): no llega a las D ===');
{
  const fc = theory.frequency(el.q, el.m, B);
  const sim = P.createCyclotron({ q: el.q, m: el.m, B, V: 50e3, relativistic: true, halfHeight: 0.05 });
  const s = sim.runUntil(300 / fc);
  const rTeo = theory.radiusRel(el.q, el.m, B, s.Kmax);
  check('0 cruces completos entre D en 300 vueltas', s.crossings === 0, `K máx ${(s.Kmax / e / 1e3).toFixed(1)} keV; radio de giro a esa energía ${(rTeo * 1e3).toFixed(2)} mm < medio hueco (5 mm)`);
}

console.log('\n=== E2: electrón inyectado a 2 y 5 MeV (f_E = f_c/γ₀): oscila y nunca se acelera ===');
{
  const fc = theory.frequency(el.q, el.m, B);
  for (const KMeV of [2, 5]) for (const V of [50e3, 100e3]) {
    const K0 = KMeV * 1e6 * e, g0 = theory.gamma(el.m, K0);
    const sim = P.createCyclotron({ q: el.q, m: el.m, B, V, relativistic: true, halfHeight: 0.05, K0, fE: fc / g0 });
    const s = sim.runUntil(300 * g0 / fc);
    check(`K₀ = ${KMeV} MeV, V = ${V / 1e3} kV: no sale y K máx < 1.5·K₀`, !s.exited && s.Kmax < 1.5 * K0, `vueltas ${s.turns}, K máx = ${(s.Kmax / K0).toFixed(2)} × K₀ (K máx teórica ${(theory.kineticMaxRel(el.q, el.m, B, 0.5) / e / 1e6).toFixed(0)} MeV)`);
  }
}

console.log('\n=== E3: electrón de 10 MeV, signo de ΔK en el primer cruce: φ = π/2 acelera, φ = 3π/2 frena ===');
{
  const K0 = 10e6 * e, g0 = theory.gamma(el.m, K0), fc = theory.frequency(el.q, el.m, B), qV = e * 100e3;
  const primer = (phi) => { const sim = P.createCyclotron({ q: el.q, m: el.m, B, V: 100e3, relativistic: true, K0, fE: fc / g0, phase: phi, inject: 'dee', y0: 0 }); sim.runUntil(3 * g0 / fc); return sim.crossings[0]; };
  const a = primer(Math.PI / 2), b = primer(1.5 * Math.PI);
  check('φ = π/2: ΔK > 0 (≈ +qV)', a.dK > 0.8 * qV, `ΔK = ${(a.dK / qV).toFixed(3)} qV`);
  check('φ = 3π/2: ΔK < 0 (≈ −qV)', b.dK < -0.8 * qV, `ΔK = ${(b.dK / qV).toFixed(3)} qV`);
}

console.log('\n=== P1: protón relativista (V = 50 kV, 400 vueltas) ===');
{
  const fc = theory.frequency(pr.q, pr.m, B), go = (r, relat) => P.createCyclotron({ q: pr.q, m: pr.m, B, V: 50e3, relativistic: relat, fE: r * fc, halfHeight: 0.05 }).runUntil(400 / fc);
  const a = go(1, true), b = go(0.995, true), cl = go(1, false);
  check('clásico, f_E = f_c: sale (referencia)', cl.exited, `K = ${(cl.K / e / 1e6).toFixed(2)} MeV`);
  check('relativista, f_E = f_c: se atasca', !a.exited && a.Kmax / e / 1e6 < 10, `K máx ${(a.Kmax / e / 1e6).toFixed(2)} MeV`);
  check('relativista, f_E = 0.995 f_c: sale', b.exited, `K = ${(b.K / e / 1e6).toFixed(2)} MeV, ${b.turns} vueltas`);
}

console.log('\n=== Ciclotrón con B hacia −z (protón, V = 50 kV): mismo resultado que +z ===');
{
  const fc = theory.frequency(pr.q, pr.m, B), go = (d) => P.createCyclotron({ q: pr.q, m: pr.m, B, Bdir: d, V: 50e3, halfHeight: 0.05 }).runUntil(400 / fc);
  const a = go([0, 0, 1]), b = go([0, 0, -1]);
  check('sale con la misma energía y vueltas', a.exited && b.exited && rel(a.K, b.K) < 5e-3 && a.turns === b.turns, `+z: ${(a.K / e / 1e6).toFixed(3)} MeV / ${a.turns} v.   −z: ${(b.K / e / 1e6).toFixed(3)} MeV / ${b.turns} v.`);
}

console.log(fallos ? `\n${fallos} comprobación(es) FALLARON` : '\nTodas las comprobaciones pasaron');
process.exit(fallos ? 1 : 0);
