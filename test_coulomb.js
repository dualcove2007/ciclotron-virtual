/* Validación de la carga externa (Coulomb) en el ciclotrón
 * Uso:  node test_coulomb.js
 * C1: aceleración inicial vs. k·q·Q/(m·d²)   C2: signo (repulsión / atracción)
 * C3: conservación de K + q·Φ con V = 0 (solo B y la carga externa) */
const P = require('./cyclotron_physics.js');
const { theory, PARTICLES, E_CHARGE: e } = P;
const pt = PARTICLES.proton, k = 8.9875517923e9, EPS2 = 0.005 ** 2;
const B = 1, T = theory.period(pt.q, pt.m, B), dt = T / 200;
const base = { q: pt.q, m: pt.m, B, V: 0, gap: 0.01, R: 10, y0: 0, dt };
// Boris es leapfrog: la velocidad guardada está en t - dt/2, así que las comparaciones usan esos tiempos.

console.log('=== C1 + C2: aceleración inicial y signo (protón en reposo en el origen, V = 0, B = 0 para aislar a Coulomb) ===');
console.log('Q (µC) | posición Q  | a simulada (m/s²) | a teórica (m/s²) | error   | sentido');
for (const [Q, qx, qy] of [[10e-6, 0.5, 0], [-10e-6, 0.5, 0], [10e-6, 0, 0.5], [30e-6, 0.3, 0.4]]) {
  const sim = P.createCyclotron({ ...base, B: 0, K0: 0, Q, qx, qy });
  const n = 20; for (let i = 0; i < n; i++) sim.step();
  const s = sim.state(), t = (n - 0.5) * dt, v = Math.hypot(s.vx, s.vy);   // v está en t = (n - 1/2)·dt
  const d = Math.hypot(qx, qy), a_th = (k * Math.abs(Q) * pt.q) / (pt.m * (d * d + EPS2) ** 1.5) * d;
  // sentido: producto punto de la velocidad con el vector (origen - Q): >0 = se aleja de Q
  const alejando = (s.vx * -qx + s.vy * -qy) > 0;
  console.log(`${(Q * 1e6).toFixed(0).padStart(6)} | (${qx}, ${qy})`.padEnd(24) +
    ` | ${(v / t).toExponential(4).padStart(17)} | ${a_th.toExponential(4).padStart(16)} | ${(100 * (v / t - a_th) / a_th).toFixed(3).padStart(6)} % | ${alejando ? 'se aleja (repulsión)' : 'se acerca (atracción)'}`);
}

console.log('\n=== C3: conservación de K + qΦ (protón 1 MeV, B = 1 T, V = 0, 20 vueltas) ===');
for (const Q of [20e-6, -20e-6]) {
  const qx = 0.3, qy = -0.1;
  const sim = P.createCyclotron({ ...base, K0: 1e6 * e, Q, qx, qy });
  const phi = (x, y) => (k * Q * pt.q) / Math.sqrt((x - qx) ** 2 + (y - qy) ** 2 + EPS2);
  // K (con v en t - dt/2) se compara con la posición a mitad de paso: (x_ahora + x_antes)/2
  let s = sim.state(), U0 = s.K + phi(s.x, s.y), maxDev = 0, dmin = Infinity, xp = s.x, yp = s.y;
  for (let i = 0; i < 20 * 200; i++) {
    sim.step(); s = sim.state();
    maxDev = Math.max(maxDev, Math.abs(s.K + phi((s.x + xp) / 2, (s.y + yp) / 2) - U0) / s.K);
    xp = s.x; yp = s.y;
    dmin = Math.min(dmin, Math.hypot(s.x - qx, s.y - qy));
  }
  console.log(`Q = ${(Q * 1e6).toFixed(0).padStart(4)} µC | desviación máx de K + qΦ: ${(100 * maxDev).toFixed(4)} % | distancia mínima a la carga: ${(dmin * 100).toFixed(1)} cm`);
}

// ---------------------------------------------------------------------------------------------
console.log('\n=== C4: deriva E×B del centro de giro (V = 0, protón 1 MeV, B = 1 T, 20 vueltas) ===');
console.log('Teoría: v_d = E×B/B²  con E evaluado en el centro de giro inicial (0, -r).');
console.log('Q (µC) | posición Q   | v_d medida (m/s)        | v_d teórica (m/s)       | error módulo | error ángulo');
{
  const e1 = P.E_CHARGE, K0 = 1e6 * e1, v0 = theory.speedFromK(pt.m, K0), r = theory.radius(pt.q, pt.m, B, v0);
  const nT = 20, spt = 400, n = nT * spt;
  for (const [Q, qx, qy] of [[1e-6, 0.75, 0], [-1e-6, 0.75, 0], [1e-6, 0, 0.75], [2e-6, -0.6, 0.5]]) {
    const sim = P.createCyclotron({ ...base, dt: T / spt, K0, Q, qx, qy });
    let sx = 0, sy = 0;
    for (let i = 0; i < n; i++) { sim.step(); const s = sim.state(); sx += s.x; sy += s.y; }
    const tt = n * T / spt;                    // el promedio de posición de un centro que deriva a u·t vale u·t/2
    const ux = (sx / n) / (tt / 2), uy = (sy / n + r) / (tt / 2);
    const cx = 0, cy = -r, dx = cx - qx, dy = cy - qy, d2 = dx * dx + dy * dy + EPS2;
    const f = k * Q / (d2 * Math.sqrt(d2)), Ex = f * dx, Ey = f * dy;          // E en el centro de giro
    const tx = Ey / B, ty = -Ex / B;                                            // E×B/B² con B en +z
    const em = (Math.hypot(ux, uy) - Math.hypot(tx, ty)) / Math.hypot(tx, ty);
    const ang = (Math.atan2(uy, ux) - Math.atan2(ty, tx)) * 180 / Math.PI;
    console.log(`${(Q * 1e6).toFixed(0).padStart(6)} | (${qx}, ${qy})`.padEnd(24) +
      ` | (${ux.toExponential(2)}, ${uy.toExponential(2)}) | (${tx.toExponential(2)}, ${ty.toExponential(2)}) | ${(100 * em).toFixed(1).padStart(6)} %    | ${ang.toFixed(1)}°`);
  }
}

console.log('\n=== C5: efecto del blindaje de las D (protón, V = 50 kV, carga en (0.75, 0)) ===');
console.log('Q (µC) | sin blindaje: K/Kteo, vueltas | con blindaje: K/Kteo, vueltas | sin carga: K/Kteo');
{
  const Kth = theory.kineticMax(pt.q, pt.m, B, 0.5);
  const go = (Q, shield) => P.createCyclotron({ q: pt.q, m: pt.m, B, V: 50e3, gap: 0.01, R: 0.5, Q, qx: 0.75, qy: 0, shield }).runUntil(3e-4);
  const r0 = go(0, false);
  for (const Q of [0.1e-6, 1e-6, 3e-6, -1e-6]) {
    const a = go(Q, false), b = go(Q, true);
    console.log(`${(Q * 1e6).toFixed(1).padStart(6)} | ${(100 * a.K / Kth).toFixed(2).padStart(8)} %, ${a.turns.toFixed(1).padStart(6)}      | ${(100 * b.K / Kth).toFixed(2).padStart(8)} %, ${b.turns.toFixed(1).padStart(6)}       | ${(100 * r0.K / Kth).toFixed(2)} %`);
  }
}
