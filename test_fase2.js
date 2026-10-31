/* Validación de la FASE 2: ciclotrón (hueco + campo E alternante)
 * Uso:  node test_fase2.js
 * Pruebas: V0 consistencia con la fase 1, V3 ganancia por cruce, V4 K vs cruces,
 *          V5 energía final, V6 barrido de frecuencia, V7 barrido de fase */
const P = require('./cyclotron_physics.js');
const { theory, PARTICLES, E_CHARGE: e } = P;

const B = 1.0, V = 50e3, GAP = 0.01, R = 0.5, K0 = 1e3 * e;
const fmtP = (x, d = 2) => (100 * x).toFixed(d) + ' %';
const MeV = (K) => (K / e / 1e6).toFixed(3);

function make(pt, extra) {
  return P.createCyclotron(Object.assign({ q: pt.q, m: pt.m, B, V, gap: GAP, R, K0 }, extra));
}
function run(pt, extra, maxTurns) {
  const sim = make(pt, extra);
  const st = sim.runUntil(maxTurns / theory.frequency(pt.q, pt.m, B));
  return { sim, st };
}
function linreg(xs, ys) {
  const n = xs.length, mx = xs.reduce((a, b) => a + b) / n, my = ys.reduce((a, b) => a + b) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (xs[i] - mx) * (ys[i] - my); sxx += (xs[i] - mx) ** 2; syy += (ys[i] - my) ** 2; }
  const slope = sxy / sxx;
  return { slope, intercept: my - slope * mx, r2: (sxy * sxy) / (sxx * syy) };
}

const proton = PARTICLES.proton;
const qV = proton.q * V;
const Kth = (pt) => theory.kineticMax(pt.q, pt.m, B, R);

console.log(`Configuración: B = ${B} T, V = ${V / 1e3} kV, hueco = ${GAP * 100} cm, R = ${R} m, K0 = 1 keV, modo ideal, f_E = f_c\n`);

// ---------- V0 ----------
{
  const T = theory.period(proton.q, proton.m, B), dt = T / 200, v0 = theory.speedFromK(proton.m, K0);
  const ref = P.createSimulation({ q: proton.q, m: proton.m, B, v0: [v0, 0, 0], dt });
  const cyc = P.createCyclotron({ q: proton.q, m: proton.m, B, V: 0, gap: GAP, R: 1e3, K0, y0: 0, dt });
  let maxDiff = 0;
  for (let i = 0; i < 2000; i++) {
    ref.step(); cyc.step();
    const a = ref.state().x, b = cyc.state();
    maxDiff = Math.max(maxDiff, Math.hypot(a[0] - b.x, a[1] - b.y));
  }
  const r = theory.radius(proton.q, proton.m, B, v0);
  console.log('=== V0: el ciclotrón con V = 0 reproduce la fase 1 (2000 pasos) ===');
  console.log(`Diferencia máxima de posición: ${(maxDiff / r).toExponential(2)} (relativa al radio)\n`);
}

// ---------- V3 y V4 ----------
const base = run(proton, {}, 300);
const cr = base.sim.crossings;
{
  console.log('=== V3: ganancia de energía por cruce del hueco (protón) ===');
  console.log('Cruce #  | ΔK / qV');
  cr.slice(0, 6).forEach((c, i) => console.log(`${String(i + 1).padStart(7)}  | ${(c.dK / qV).toFixed(4)}`));
  const rest = cr.slice(6).map((c) => c.dK / qV);
  const mean = rest.reduce((a, b) => a + b) / rest.length;
  console.log(`Cruces 7 a ${cr.length}: media ${mean.toFixed(4)}, mín ${Math.min(...rest).toFixed(4)}, máx ${Math.max(...rest).toFixed(4)}`);
  console.log(`Rebotes (entra y vuelve al mismo lado): ${base.st.bounces}\n`);

  const xs = cr.map((_, i) => i + 1), ys = cr.map((c) => c.K_out);
  const fit = linreg(xs, ys);
  console.log('=== V4: energía vs. número de cruces (ajuste lineal) ===');
  console.log(`Pendiente: ${(fit.slope / qV).toFixed(5)} qV por cruce (teoría: 1)   R² = ${fit.r2.toFixed(8)}`);
  console.log(`Ordenada: ${(fit.intercept / e / 1e3).toFixed(2)} keV (teoría K0 + qV/2 = ${((K0 + qV / 2) / e / 1e3).toFixed(2)} keV)`);
  const Nturns = base.st.turns;
  console.log(`Relación K = 2·N·qV: N = ${Nturns} vueltas → 2NqV = ${MeV(2 * Nturns * qV)} MeV, K simulada = ${MeV(base.st.K)} MeV\n`);
}

// ---------- V5 ----------
console.log('=== V5: energía final vs. K = q²B²R²/(2m) ===');
console.log('Partícula        | K simulada (MeV) | K teórica (MeV) | error   | vueltas | f_c (MHz)');
for (const k of Object.keys(PARTICLES).filter((k) => k !== 'electron')) {
  const pt = PARTICLES[k], { sim, st } = run(pt, {}, 400);
  console.log(`${pt.name.padEnd(16)} | ${MeV(st.K).padStart(16)} | ${MeV(Kth(pt)).padStart(15)} | ${fmtP((st.K - Kth(pt)) / Kth(pt)).padStart(7)} | ${String(st.turns).padStart(7)} | ${(sim.fc / 1e6).toFixed(3)}`);
}
console.log('Nota: K simulada es la energía de la última semicircunferencia, así que puede diferir ~0,4 % por cuantización (ΔK por cruce / K).\n');

// ---------- V6 ----------
const y0res = base.sim.y0;                  // mismo punto de inyección en todo el barrido
const fc = theory.frequency(proton.q, proton.m, B);
console.log('=== V6: barrido de frecuencia (protón, hasta 300 vueltas) ===');
console.log('f_E / f_c | sale | K máx alcanzada / K teórica | vueltas');
for (const r of [0.95, 0.98, 0.99, 0.995, 0.998, 0.999, 1.0, 1.001, 1.002, 1.005, 1.01, 1.02, 1.05]) {
  const { st } = run(proton, { fE: r * fc, y0: y0res }, 300);
  console.log(`${r.toFixed(3).padStart(9)} | ${(st.exited ? 'sí' : 'no').padEnd(4)} | ${fmtP(st.Kmax / Kth(proton), 1).padStart(27)} | ${st.turns}`);
}
console.log('La resonancia es estrecha: la fase se desliza 2π·(f_E − f_c)·t, y con ~120 vueltas (7,6 µs) basta un desajuste de ~0,2 % para perderla.\n');

// ---------- V7 ----------
{
  const KD = 10e3 * e;                                           // 10 keV al inyectar
  const ref = make(proton, { K0: KD, inject: 'dee', phase: Math.PI / 2 });
  const y0dee = ref.y0;
  console.log('=== V7: barrido de la fase φ, partícula inyectada en la D (10 keV, f_E = f_c, hasta 300 vueltas) ===');
  console.log('φ        | ideal: sale | rebotes | ideal: K final / K teór. | sinus.: sale | sinus.: K final / K teór. | sinus.: vueltas');
  for (let i = 0; i < 16; i++) {
    const phi = (i * Math.PI) / 8;
    const a = run(proton, { K0: KD, inject: 'dee', phase: phi, y0: y0dee, mode: 'ideal' }, 300).st;
    const b = run(proton, { K0: KD, inject: 'dee', phase: phi, y0: y0dee, mode: 'sine' }, 300).st;
    const lab = `${(i / 8).toFixed(3)} π`;
    console.log(`${lab.padEnd(8)} | ${(a.exited ? 'sí' : 'no').padEnd(11)} | ${String(a.bounces).padStart(7)} | ${fmtP(a.Kmax / Kth(proton), 1).padStart(23)} | ${(b.exited ? 'sí' : 'no').padEnd(12)} | ${fmtP(b.Kmax / Kth(proton), 1).padStart(24)} | ${b.turns}`);
  }
  console.log('Convención: E_x = (V/g)·s(t), s = sign(sen(2πft+φ)) (ideal) o sen(2πft+φ) (sinusoidal). La partícula llega al hueco con E_x < 0 si φ ∈ (0, π).');
  console.log('φ ∈ (0, π): gana energía desde el primer cruce (0 rebotes). φ ∈ (π, 2π): el campo la frena y la rebota en el primer cruce; el retraso de media');
  console.log('vuelta la deja de nuevo en fase y termina saliendo (autocorrección). Sale con K algo menor (~89 %) porque el punto de inyección queda descentrado.');
  console.log('En sinusoidal la ganancia por cruce es ~qV·sen(φ): las vueltas necesarias son mínimas en φ = π/2 y crecen hacia los bordes.\n');

  console.log('--- Contraste: inyección en el CENTRO del hueco (modo ideal, y0 recalibrado para cada φ) ---');
  console.log('φ        | sale | K final / K teórica');
  for (const f of [0.25, 0.5, 1.25, 1.5]) {
    const { st } = run(proton, { phase: f * Math.PI, inject: 'center', mode: 'ideal' }, 300);
    console.log(`${(f.toFixed(2) + ' π').padEnd(8)} | ${(st.exited ? 'sí' : 'no').padEnd(4)} | ${fmtP(st.K / Kth(proton), 2)}`);
  }
  console.log('Con inyección central, el primer empujón decide hacia qué D va la partícula: se autosincroniza y φ casi no importa.');
}
