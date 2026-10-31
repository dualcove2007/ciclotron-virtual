/* Validación de la FASE 1: partícula en B uniforme (E = 0)
 * Uso:  node test_fase1.js
 * Pruebas: V1 radio, V2 período, V8 Boris vs RK4, V9 sensibilidad a dt */
const P = require('./cyclotron_physics.js');
const { theory, PARTICLES, E_CHARGE } = P;
const cross = P.vec.cross, norm = P.vec.norm;

const B = 1.0;                 // tesla
const K0_MEV = 1.0;            // energía inicial de prueba (MeV)
const pct = (x) => (100 * x).toFixed(4) + ' %';

/* Corre `turns` vueltas y mide: error de radio (máx desviación respecto al círculo ideal),
 * error de período (por ángulo acumulado) y deriva relativa de energía. */
function runCircle(pt, integrator, dtFactor, turns) {
  const { q, m } = pt;
  const T = theory.period(q, m, B);
  const dt = T / dtFactor;
  const v0 = theory.speedFromK(m, K0_MEV * 1e6 * E_CHARGE);
  const r = theory.radius(q, m, B, v0);

  const sim = P.createSimulation({ q, m, B, v0: [v0, 0, 0], dt, integrator });
  // centro teórico del círculo: c = x + m/(qB²) (v × B)
  const vxB = cross([v0, 0, 0], [0, 0, B]);
  const c = [0 + (m / (q * B * B)) * vxB[0], 0 + (m / (q * B * B)) * vxB[1], 0];

  const K_ini = sim.state().K;
  let maxRadErr = 0, angle = 0, prevAng = Math.atan2(0 - c[1], 0 - c[0]);
  let tPeriod = null;
  const total = Math.round(turns * dtFactor);
  for (let i = 0; i < total; i++) {
    sim.step();
    const s = sim.state();
    const dx = s.x[0] - c[0], dy = s.x[1] - c[1];
    maxRadErr = Math.max(maxRadErr, Math.abs(Math.hypot(dx, dy) - r) / r);
    let a = Math.atan2(dy, dx), d = a - prevAng;
    if (d > Math.PI) d -= 2 * Math.PI; if (d < -Math.PI) d += 2 * Math.PI;
    const prevAngle = angle; angle += d; prevAng = a;
    if (tPeriod === null && Math.abs(angle) >= 2 * Math.PI) {   // interpolación lineal
      const f = (2 * Math.PI - Math.abs(prevAngle)) / (Math.abs(angle) - Math.abs(prevAngle));
      tPeriod = s.t - sim.params.dt + f * sim.params.dt;
    }
  }
  const K_fin = sim.state().K;
  return { T, r, maxRadErr, periodErr: tPeriod ? Math.abs(tPeriod - T) / T : NaN, driftE: Math.abs(K_fin - K_ini) / K_ini };
}

console.log('=== V1 + V2: radio y período (Boris, dt = T/200, 10 vueltas, B = 1 T, K0 = 1 MeV) ===');
console.log('Partícula        | r teórico (cm) | error radio máx | T teórico (ns) | error período');
for (const k of Object.keys(PARTICLES).filter((k) => k !== 'electron')) {
  const pt = PARTICLES[k], res = runCircle(pt, 'boris', 200, 10);
  console.log(`${pt.name.padEnd(16)} | ${(res.r * 100).toFixed(4).padStart(14)} | ${pct(res.maxRadErr).padStart(15)} | ${(res.T * 1e9).toFixed(4).padStart(14)} | ${pct(res.periodErr)}`);
}

console.log('\n=== V8: Boris vs RK4, deriva de energía tras 1000 vueltas (protón, dt = T/50) ===');
for (const ig of ['boris', 'rk4']) {
  const res = runCircle(PARTICLES.proton, ig, 50, 1000);
  console.log(`${ig.padEnd(6)} | deriva relativa de energía: ${res.driftE.toExponential(3)} | error radio máx: ${pct(res.maxRadErr)}`);
}

console.log('\n=== V9: sensibilidad a dt (protón, Boris, 10 vueltas) ===');
console.log('dt      | error radio máx | error período');
for (const f of [25, 50, 100, 200, 400]) {
  const res = runCircle(PARTICLES.proton, 'boris', f, 10);
  console.log(`T/${String(f).padEnd(5)} | ${pct(res.maxRadErr).padStart(15)} | ${pct(res.periodErr)}`);
}
