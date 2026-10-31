/* Exporta los datos de validación (fases 1 y 2) a CSV en ./datos
 * Uso:  node export_datos.js   (luego: python graficas.py) */
const fs = require('fs');
const P = require('./cyclotron_physics.js');
const { theory, PARTICLES, E_CHARGE: e } = P;

const B = 1.0, V = 50e3, GAP = 0.01, R = 0.5, K0 = 1e3 * e;
const proton = PARTICLES.proton;
const qV = proton.q * V;
const Kth = (pt) => theory.kineticMax(pt.q, pt.m, B, R);
const fc = theory.frequency(proton.q, proton.m, B);

fs.mkdirSync('datos', { recursive: true });
const save = (name, header, rows) => {
  fs.writeFileSync(`datos/${name}`, [header, ...rows.map((r) => r.join(','))].join('\n') + '\n');
  console.log(`datos/${name}  (${rows.length} filas)`);
};
const make = (pt, extra) => P.createCyclotron(Object.assign({ q: pt.q, m: pt.m, B, V, gap: GAP, R, K0 }, extra));
const run = (pt, extra, maxTurns) => {
  const sim = make(pt, extra);
  return { sim, st: sim.runUntil(maxTurns / theory.frequency(pt.q, pt.m, B)) };
};

// 1) K vs cruces (protón) + ganancia por cruce
{
  const { sim } = run(proton, {}, 300);
  save('k_vs_cruces.csv', 'cruce,K_MeV,dK_sobre_qV',
    sim.crossings.map((c, i) => [i + 1, (c.K_out / e / 1e6).toFixed(6), (c.dK / qV).toFixed(6)]));
}

// 2) Trayectoria del protón (submuestreada)
{
  const sim = make(proton, {});
  const rows = [];
  let i = 0;
  while (!sim.state().exited && sim.state().t < 300 / fc) {
    sim.step();
    if (i++ % 40 === 0) { const s = sim.state(); rows.push([s.t.toExponential(5), s.x.toFixed(6), s.y.toFixed(6), (s.K / e / 1e6).toFixed(6)]); }
  }
  save('trayectoria_proton.csv', 't_s,x_m,y_m,K_MeV', rows);
}

// 3) Energía final por partícula
{
  const rows = Object.keys(PARTICLES).filter((k) => k !== 'electron').map((k) => {
    const pt = PARTICLES[k], { sim, st } = run(pt, {}, 400);
    const Ks = st.K / e / 1e6, Kt = Kth(pt) / e / 1e6;
    return [pt.name, Ks.toFixed(4), Kt.toFixed(4), (100 * (Ks - Kt) / Kt).toFixed(4), st.turns, (sim.fc / 1e6).toFixed(4)];
  });
  save('energia_final.csv', 'particula,K_sim_MeV,K_teorica_MeV,error_pct,vueltas,fc_MHz', rows);
}

// 4) Barrido de frecuencia (protón)
{
  const y0 = run(proton, {}, 300).sim.y0, rows = [];
  for (let r = 0.95; r <= 1.0501; r += 0.0025) {
    const { st } = run(proton, { fE: r * fc, y0 }, 300);
    rows.push([r.toFixed(4), st.exited ? 1 : 0, (st.Kmax / Kth(proton)).toFixed(5), st.turns]);
  }
  for (const r of [0.998, 0.999, 1.001, 1.002]) {          // puntos finos cerca de la resonancia
    const { st } = run(proton, { fE: r * fc, y0 }, 300);
    rows.push([r.toFixed(4), st.exited ? 1 : 0, (st.Kmax / Kth(proton)).toFixed(5), st.turns]);
  }
  rows.sort((a, b) => a[0] - b[0]);
  save('barrido_frecuencia.csv', 'fE_sobre_fc,sale,Kmax_sobre_Kteo,vueltas', rows);
}

// 5) Barrido de fase (inyección en la D, 10 keV)
{
  const KD = 10e3 * e, y0 = make(proton, { K0: KD, inject: 'dee', phase: Math.PI / 2 }).y0, rows = [];
  for (let i = 0; i < 32; i++) {
    const phi = (i * Math.PI) / 16;
    const a = run(proton, { K0: KD, inject: 'dee', phase: phi, y0, mode: 'ideal' }, 300).st;
    const b = run(proton, { K0: KD, inject: 'dee', phase: phi, y0, mode: 'sine' }, 300).st;
    rows.push([(i / 16).toFixed(4), (a.Kmax / Kth(proton)).toFixed(5), a.bounces, (b.Kmax / Kth(proton)).toFixed(5), b.turns]);
  }
  save('barrido_fase.csv', 'phi_sobre_pi,ideal_K_sobre_Kteo,ideal_rebotes,sinus_K_sobre_Kteo,sinus_vueltas', rows);
}

// 6) Boris vs RK4: deriva de energía (dt = T/50, 1000 vueltas, cada 10 vueltas)
{
  const T = theory.period(proton.q, proton.m, B), dt = T / 50;
  const v0 = theory.speedFromK(proton.m, 1e6 * e);
  const sims = ['boris', 'rk4'].map((ig) => P.createSimulation({ q: proton.q, m: proton.m, B, v0: [v0, 0, 0], dt, integrator: ig }));
  const K_ini = sims.map((s) => s.state().K), rows = [];
  for (let n = 1; n <= 1000; n++) {
    for (let i = 0; i < 50; i++) sims.forEach((s) => s.step());
    if (n % 10 === 0) rows.push([n, ...sims.map((s, j) => Math.abs(s.state().K - K_ini[j]) / K_ini[j])].map((x, j) => (j ? x.toExponential(4) : x)));
  }
  save('boris_vs_rk4.csv', 'vueltas,deriva_boris,deriva_rk4', rows);
}

// 7) Sensibilidad a dt (error de radio máx, 10 vueltas, K0 = 1 MeV)
{
  const v0 = theory.speedFromK(proton.m, 1e6 * e), r = theory.radius(proton.q, proton.m, B, v0);
  const T = theory.period(proton.q, proton.m, B), cx = 0, cy = -r;   // centro: q>0, B en +z, v en +x
  const rows = [];
  for (const f of [25, 50, 100, 200, 400]) {
    const sim = P.createSimulation({ q: proton.q, m: proton.m, B, v0: [v0, 0, 0], dt: T / f });
    let maxErr = 0;
    for (let i = 0; i < 10 * f; i++) { sim.step(); const s = sim.state(); maxErr = Math.max(maxErr, Math.abs(Math.hypot(s.x[0] - cx, s.x[1] - cy) - r) / r); }
    rows.push([f, (100 * maxErr).toFixed(5)]);
  }
  save('sensibilidad_dt.csv', 'pasos_por_periodo,error_radio_max_pct', rows);
}

// ---------------------------------------------------------------------------------------------
// 8) CARGA EXTERNA: K final vs Q (con y sin blindaje de las D), carga en (0.75, 0)
const Kth_p = Kth(proton);
const conCarga = (Q, qx, qy, shield) => P.createCyclotron({ q: proton.q, m: proton.m, B, V, gap: GAP, R, K0, Q, qx, qy, shield });
{
  const rows = [];
  for (const Qn of [-3000, -1000, -300, -100, -30, -10, -3, -1, 0, 1, 3, 10, 30, 100, 300, 1000, 3000]) {
    const a = conCarga(Qn * 1e-9, 0.75, 0, false).runUntil(3e-4), b = conCarga(Qn * 1e-9, 0.75, 0, true).runUntil(3e-4);
    rows.push([Qn, (a.K / Kth_p).toFixed(5), a.turns, (b.K / Kth_p).toFixed(5), b.turns]);
  }
  save('carga_vs_Q.csv', 'Q_nC,sin_blindaje_K_sobre_Kteo,sin_blindaje_vueltas,con_blindaje_K_sobre_Kteo,con_blindaje_vueltas', rows);
}
// 9) K final vs distancia de la carga (Q = +1 µC sobre el eje x)
{
  const rows = [];
  for (const d of [0.55, 0.6, 0.7, 0.8, 1.0, 1.25, 1.5, 2.0, 3.0, 5.0]) {
    const a = conCarga(1e-6, d, 0, false).runUntil(3e-4), b = conCarga(1e-6, d, 0, true).runUntil(3e-4);
    rows.push([d, (a.K / Kth_p).toFixed(5), (b.K / Kth_p).toFixed(5)]);
  }
  save('carga_vs_distancia.csv', 'd_m,sin_blindaje_K_sobre_Kteo,con_blindaje_K_sobre_Kteo', rows);
}
// 10) Deriva E×B medida vs teórica (V = 0, protón 1 MeV)
{
  const K = 1e6 * e, Tc = 1 / fc, spt = 400, nT = 20, n = nT * spt, kC = 8.9875517923e9, EPS2 = 0.005 ** 2;
  const r = theory.radius(proton.q, proton.m, B, theory.speedFromK(proton.m, K));
  const rows = [];
  for (const [Q, qx, qy] of [[1e-6, 0.75, 0], [-1e-6, 0.75, 0], [1e-6, 0, 0.75], [2e-6, -0.6, 0.5], [0.5e-6, 0.5, 0.5], [3e-6, 1.2, -0.8], [-2e-6, 0, -1.0]]) {
    const sim = P.createCyclotron({ q: proton.q, m: proton.m, B, V: 0, gap: GAP, R: 10, y0: 0, K0: K, dt: Tc / spt, Q, qx, qy });
    let sx = 0, sy = 0;
    for (let i = 0; i < n; i++) { sim.step(); const s = sim.state(); sx += s.x; sy += s.y; }
    const tt = n * Tc / spt, ux = (sx / n) / (tt / 2), uy = (sy / n + r) / (tt / 2);
    const dx = -qx, dy = -r - qy, d2 = dx * dx + dy * dy + EPS2, f = kC * Q / (d2 * Math.sqrt(d2));
    const tx = (f * dy) / B, ty = -(f * dx) / B;
    rows.push([Q * 1e6, qx, qy, Math.hypot(ux, uy).toFixed(2), Math.hypot(tx, ty).toFixed(2)]);
  }
  save('deriva_ExB.csv', 'Q_uC,qx_m,qy_m,vd_medida_ms,vd_teorica_ms', rows);
}
// 11) Trayectorias: sin carga, con Q = 1 µC, y con Q = 1 µC y blindaje
{
  const rows = [];
  for (const [nombre, Q, sh] of [['sin_carga', 0, false], ['Q1uC_sin_blindaje', 1e-6, false], ['Q1uC_con_blindaje', 1e-6, true]]) {
    const sim = conCarga(Q, 0.75, 0, sh);
    let i = 0;
    while (!sim.state().exited && sim.state().t < 300 / fc) {
      sim.step();
      if (i++ % 40 === 0) { const s = sim.state(); rows.push([nombre, s.x.toFixed(6), s.y.toFixed(6), (s.K / e / 1e6).toFixed(6)]); }
    }
    const s = sim.state(); rows.push([nombre, s.x.toFixed(6), s.y.toFixed(6), (s.K / e / 1e6).toFixed(6)]);
  }
  save('trayectorias_carga.csv', 'caso,x_m,y_m,K_MeV', rows);
}
