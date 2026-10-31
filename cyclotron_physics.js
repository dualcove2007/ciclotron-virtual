/* ==========================================================
 * Ciclotrón virtual 3D - CAPA DE FÍSICA (sin gráficos)
 * Fase 1: partícula en campo magnético uniforme (E = 0)
 * Funciona en el navegador (window.CyclotronPhysics) y en Node.
 * Unidades: SI. B uniforme a lo largo de +z.
 * ========================================================== */
(function (root) {
  'use strict';

  const E_CHARGE = 1.602176634e-19;      // C
  const PARTICLES = {
    proton:   { name: 'Protón',         q: E_CHARGE,     m: 1.67262192e-27 },
    deuteron: { name: 'Deuterón',       q: E_CHARGE,     m: 3.3435838e-27  },
    alpha:    { name: 'Partícula alfa', q: 2 * E_CHARGE, m: 6.6446573e-27  },
    electron: { name: 'Electrón',       q: -E_CHARGE,    m: 9.1093837015e-31 },
  };
  const C_LIGHT = 299792458;             // m/s

  // ---------- utilidades vectoriales (arrays [x,y,z]) ----------
  const cross = (a, b) => [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const norm = (a) => Math.sqrt(dot(a, a));

  // ---------- fórmulas teóricas ----------
  const theory = {
    radius:    (q, m, B, v) => (m * v) / (Math.abs(q) * B),
    period:    (q, m, B)    => (2 * Math.PI * m) / (Math.abs(q) * B),
    frequency: (q, m, B)    => (Math.abs(q) * B) / (2 * Math.PI * m),
    omega:     (q, m, B)    => (Math.abs(q) * B) / m,
    kineticMax:(q, m, B, R) => (q * q * B * B * R * R) / (2 * m),   // K = q²B²R²/(2m)
    speedFromK:(m, K)       => Math.sqrt((2 * K) / m),
    // ---- versión relativista (p = γmv, K = (γ-1)mc²) ----
    gamma:     (m, K)       => 1 + K / (m * C_LIGHT * C_LIGHT),
    speedFromKRel: (m, K)   => { const g = 1 + K / (m * C_LIGHT * C_LIGHT); return C_LIGHT * Math.sqrt(1 - 1 / (g * g)); },
    // radio de giro en B con momento p:  r = p/(|q|B)
    radiusRel: (q, m, B, K) => { const mc2 = m * C_LIGHT * C_LIGHT, pc = Math.sqrt(K * (K + 2 * mc2)); return pc / (C_LIGHT * Math.abs(q) * B); },
    // energía cinética cuando el radio de giro llega a R:  pc = |q|BRc,  K = sqrt((pc)² + (mc²)²) - mc²
    kineticMaxRel: (q, m, B, R) => { const mc2 = m * C_LIGHT * C_LIGHT, pc = Math.abs(q) * B * R * C_LIGHT; return Math.sqrt(pc * pc + mc2 * mc2) - mc2; },
    // frecuencia de giro con la masa relativista: f = qB/(2πγm)
    frequencyRel: (q, m, B, K) => (Math.abs(q) * B) / (2 * Math.PI * m * (1 + K / (m * C_LIGHT * C_LIGHT))),
  };

  // ---------- paso de Boris (no relativista) ----------
  // v_nuevo = Boris(v, E, B, q/m, dt). Funciona también con dt negativo
  // (sirve para el medio paso hacia atrás de la inicialización leapfrog).
  function borisVelocity(v, E, Bvec, qm, dt) {
    const h = 0.5 * qm * dt;
    const vm = [v[0] + h * E[0], v[1] + h * E[1], v[2] + h * E[2]];          // medio impulso E
    const t  = [h * Bvec[0], h * Bvec[1], h * Bvec[2]];
    const t2 = dot(t, t);
    const s  = [2 * t[0] / (1 + t2), 2 * t[1] / (1 + t2), 2 * t[2] / (1 + t2)];
    const vp = cross(vm, t); vp[0] += vm[0]; vp[1] += vm[1]; vp[2] += vm[2];  // v' = v- + v- x t
    const vs = cross(vp, s); vs[0] += vm[0]; vs[1] += vm[1]; vs[2] += vm[2];  // v+ = v- + v' x s
    return [vs[0] + h * E[0], vs[1] + h * E[1], vs[2] + h * E[2]];            // medio impulso E
  }

  // ---------- paso RK4 (solo para comparar con Boris) ----------
  function rk4Step(x, v, E, Bvec, qm, dt) {
    const acc = (vel) => {
      const c = cross(vel, Bvec);
      return [qm * (E[0] + c[0]), qm * (E[1] + c[1]), qm * (E[2] + c[2])];
    };
    const add = (a, b, k) => [a[0] + k * b[0], a[1] + k * b[1], a[2] + k * b[2]];
    const k1v = acc(v),                 k1x = v;
    const k2v = acc(add(v, k1v, dt / 2)), k2x = add(v, k1v, dt / 2);
    const k3v = acc(add(v, k2v, dt / 2)), k3x = add(v, k2v, dt / 2);
    const k4v = acc(add(v, k3v, dt)),     k4x = add(v, k3v, dt);
    const vn = [0, 0, 0], xn = [0, 0, 0];
    for (let i = 0; i < 3; i++) {
      vn[i] = v[i] + (dt / 6) * (k1v[i] + 2 * k2v[i] + 2 * k3v[i] + k4v[i]);
      xn[i] = x[i] + (dt / 6) * (k1x[i] + 2 * k2x[i] + 2 * k3x[i] + k4x[i]);
    }
    return { x: xn, v: vn };
  }

  // ---------- simulación ----------
  // params: { q, m, B (T), x0:[..], v0:[..], dt, integrator:'boris'|'rk4' }
  function createSimulation(params) {
    const p = Object.assign({ x0: [0, 0, 0], v0: [0, 0, 0], integrator: 'boris', B: 1 }, params);
    const qm = p.q / p.m;
    const Bvec = [0, 0, p.B];
    const Zero = [0, 0, 0];
    const state = { t: 0, x: p.x0.slice(), v: p.v0.slice(), steps: 0 };

    // Boris es "leapfrog": v vive en t - dt/2. Se retrocede medio paso al iniciar.
    if (p.integrator === 'boris') state.v = borisVelocity(state.v, Zero, Bvec, qm, -p.dt / 2);

    const kinetic = () => 0.5 * p.m * dot(state.v, state.v);

    function step() {
      if (p.integrator === 'boris') {
        state.v = borisVelocity(state.v, Zero, Bvec, qm, p.dt);
        state.x = [state.x[0] + state.v[0] * p.dt, state.x[1] + state.v[1] * p.dt, state.x[2] + state.v[2] * p.dt];
      } else {
        const r = rk4Step(state.x, state.v, Zero, Bvec, qm, p.dt);
        state.x = r.x; state.v = r.v;
      }
      state.t += p.dt; state.steps++;
    }

    return {
      params: p,
      step,
      state: () => ({ t: state.t, x: state.x.slice(), v: state.v.slice(), K: kinetic(), speed: norm(state.v), steps: state.steps }),
    };
  }

  // ==========================================================
  // FASE 2: ciclotrón = campo B uniforme + hueco con campo E alternante
  // Geometría: hueco de ancho `gap` centrado en x = 0 (E solo ahí, E = Ex).
  //            Cada D es un semicírculo de radio R centrado en el borde del hueco (±gap/2, 0).
  // FASE 3 (ampliación):
  //   - Carga q de cualquier signo (electrón). Para q < 0 la fase se desplaza π internamente, de modo que
  //     φ = π/2 siga siendo "la fase buena" para cualquier partícula (equivale a invertir la polaridad de las D).
  //   - Dirección de B arbitraria: `Bdir` = [bx, by, bz] (se normaliza). El movimiento deja de ser plano (hélice).
  //   - Dinámica relativista opcional (`relativistic: true`): paso de Boris con u = γv.
  //   - `halfHeight`: semialtura de las D (m). Si |z| la supera, la partícula "sale" (choca con la tapa).
  // ==========================================================
  function suggestDt(q, m, B, gap, R, samplesPerTurn, samplesPerGap, rel) {
    const T = theory.period(q, m, B);
    const Kmax = rel ? theory.kineticMaxRel(q, m, B, R) : theory.kineticMax(q, m, B, R);
    const vmax = rel ? theory.speedFromKRel(m, Kmax) : theory.speedFromK(m, Kmax);
    // dt debe resolver tanto la vuelta como el cruce del hueco a máxima velocidad
    return Math.min(T / samplesPerTurn, gap / (samplesPerGap * vmax));
  }

  // params: { q, m, B, Bdir, V, gap, R, fE, phase, mode:'ideal'|'sine', K0, y0, dt, relativistic, halfHeight, ... }
  function createCyclotron(params) {
    const p = Object.assign({
      B: 1, Bdir: [0, 0, 1], V: 50e3, gap: 0.01, R: 0.5, fE: null, phase: Math.PI / 2, mode: 'ideal',
      K0: 1e3 * E_CHARGE, y0: 'auto', inject: 'center', samplesPerTurn: 200, samplesPerGap: 20, dt: null,
      relativistic: false, halfHeight: Infinity,
      Q: 0, qx: 0.75, qy: 0,        // carga externa fija Q (C) en (qx, qy) metros; Q = 0 la desactiva
      shield: false,               // true = las D (metal) blindan el campo de Q: dentro de ellas Eext = 0
    }, params);
    const { q, m, B } = p;
    const rel = !!p.relativistic;
    const qm = q / m;
    const fc = theory.frequency(q, m, B);
    if (p.fE == null) p.fE = fc;
    if (p.dt == null) p.dt = suggestDt(q, m, B, p.gap, p.R, p.samplesPerTurn, p.samplesPerGap, rel);
    const dt = p.dt, hg = p.gap / 2, EV = p.V / p.gap, R2 = p.R * p.R, zMax = p.halfHeight;
    const omegaE = 2 * Math.PI * p.fE;
    const ideal = p.mode === 'ideal';
    const phaseEff = p.phase + (q < 0 ? Math.PI : 0);   // ver nota de la fase 3
    const C2 = C_LIGHT * C_LIGHT;

    // dirección de B (unitaria) y vector B
    const bn = norm(p.Bdir);
    const bhat = bn > 0 ? [p.Bdir[0] / bn, p.Bdir[1] / bn, p.Bdir[2] / bn] : [0, 0, 1];
    const Bx = B * bhat[0], By = B * bhat[1], Bz = B * bhat[2];
    const alongZ = Math.abs(bhat[2]) > 1 - 1e-9;       // B paralelo (o antiparalelo) a z: el movimiento es plano

    // Inyección (y0): el radio crece más rápido al inicio, así que el espiral NO queda centrado
    // si se inyecta en y = 0. Con y0 = 'auto' se hace una pre-corrida (sin salida) con y0 = 0,
    // se mide el centro del espiral a media energía y se inyecta desplazado para centrarlo
    // (equivale a centrar el haz en una máquina real). Trasladar y0 traslada todo el espiral.
    // Solo tiene sentido si el movimiento es plano (B ∥ z); si no, y0 = 0.
    let y0;
    if (p.y0 === 'auto') {
      if (alongZ) {
        const pre = createCyclotron(Object.assign({}, p, { y0: 0, R: 10 * p.R, dt, fE: p.fE, Q: 0, halfHeight: Infinity }));
        const Kmax0 = rel ? theory.kineticMaxRel(q, m, B, p.R) : theory.kineticMax(q, m, B, p.R);
        const Ktarget = 0.5 * Kmax0;
        while (pre.state().K < Ktarget && pre.state().t < 400 / fc) pre.step();
        const c = pre.crossings, n = c.length;
        y0 = (n >= 4) ? -(c[n - 3].y + 2 * c[n - 2].y + c[n - 1].y) / 4 : 0;
      } else y0 = 0;
    } else y0 = p.y0;

    // estado: posición en t, velocidad en t - dt/2 (esquema leapfrog de Boris)
    // inject: 'center' = nace en el centro del hueco; 'dee' = nace en el borde de la D derecha
    // (llega al hueco en un instante fijo, así la fase del campo decide si gana o pierde energía)
    let t = 0, x = (p.inject === 'dee') ? hg + 1e-9 : 0, y = y0, z = 0;
    let vx = rel ? theory.speedFromKRel(m, p.K0) : theory.speedFromK(m, p.K0), vy = 0, vz = 0;

    // Campo de Coulomb de la carga externa (con suavizado de 5 mm para evitar la singularidad).
    // La carga está en el plano z = 0; si la partícula se sale del plano (B inclinado) se usa la distancia 3D.
    const K_COULOMB = 8.9875517923e9, EPS2 = 0.005 * 0.005;
    let Qe = p.Q, qxe = p.qx, qye = p.qy, shield = !!p.shield;
    // Dentro de una D (semicírculo de radio R, fuera del hueco) el conductor blinda el campo externo.
    // Es un modelo simple: no calcula cargas inducidas, solo anula Eext dentro de las D.
    const inDee = (xx, yy) => { const d = Math.abs(xx) - hg; return d > 0 && d * d + yy * yy < R2; };
    function Eext(xx, yy, zz) {
      if (Qe === 0) return [0, 0, 0];
      if (shield && inDee(xx, yy)) return [0, 0, 0];
      const dx = xx - qxe, dy = yy - qye, d2 = dx * dx + dy * dy + zz * zz + EPS2;
      const f = K_COULOMB * Qe / (d2 * Math.sqrt(d2));
      return [f * dx, f * dy, f * zz];
    }
    function Efield(xx, yy, zz, tt) {                 // devuelve Ex en (x, y, z, t)
      if (Math.abs(xx) < hg && Math.abs(yy) < p.R && Math.abs(zz) < zMax) {
        const s = Math.sin(omegaE * tt + phaseEff);
        return EV * (ideal ? (s >= 0 ? 1 : -1) : s);
      }
      return 0;
    }
    // Paso de Boris con B en dirección arbitraria. En modo relativista se trabaja con u = γv:
    //   u- = u + (q dt/2m)E ;  t = (q dt/2m)B/γ- ;  rotación ;  u+ + (q dt/2m)E
    // (con rel = false y B ∥ z coincide exactamente con la versión escalar de la fase 2).
    function kick(Ex, dtl, Ey, Ez) {
      Ey = Ey || 0; Ez = Ez || 0;
      const h = 0.5 * qm * dtl;
      let ux = vx, uy = vy, uz = vz;
      if (rel) { const g0 = 1 / Math.sqrt(1 - (vx * vx + vy * vy + vz * vz) / C2); ux *= g0; uy *= g0; uz *= g0; }
      const umx = ux + h * Ex, umy = uy + h * Ey, umz = uz + h * Ez;          // medio impulso E
      const gm = rel ? Math.sqrt(1 + (umx * umx + umy * umy + umz * umz) / C2) : 1;
      const tx = h * Bx / gm, ty = h * By / gm, tz = h * Bz / gm;
      const t2 = tx * tx + ty * ty + tz * tz, f = 2 / (1 + t2);
      const upx = umx + (umy * tz - umz * ty), upy = umy + (umz * tx - umx * tz), upz = umz + (umx * ty - umy * tx);   // u' = u- + u- x t
      const sx = f * tx, sy = f * ty, sz = f * tz;
      const unx = umx + (upy * sz - upz * sy) + h * Ex;                                // u+ = u- + u' x s, y medio impulso E
      const uny = umy + (upz * sx - upx * sz) + h * Ey;
      const unz = umz + (upx * sy - upy * sx) + h * Ez;
      if (rel) { const g1 = Math.sqrt(1 + (unx * unx + uny * uny + unz * unz) / C2); vx = unx / g1; vy = uny / g1; vz = unz / g1; }
      else { vx = unx; vy = uny; vz = unz; }
    }
    { const e0 = Eext(x, y, z); kick(Efield(x, y, z, 0) + e0[0], -dt / 2, e0[1], e0[2]); }   // medio paso hacia atrás (incluye la carga externa)

    const crossings = [];                             // pasos completos de una D a la otra
    let nBounce = 0, exited = false, exitInfo = null;
    const gammaNow = () => rel ? 1 / Math.sqrt(1 - (vx * vx + vy * vy + vz * vz) / C2) : 1;
    const K = () => {
      const v2 = vx * vx + vy * vy + vz * vz;
      return rel ? (gammaNow() - 1) * m * C2 : 0.5 * m * v2;
    };
    // radio de giro (de la hélice): r = γ m v⊥ / (|q| B), con v⊥ = componente de v perpendicular a B
    const rCurv = () => {
      const v2 = vx * vx + vy * vy + vz * vz, vpar = vx * bhat[0] + vy * bhat[1] + vz * bhat[2];
      return gammaNow() * m * Math.sqrt(Math.max(0, v2 - vpar * vpar)) / (Math.abs(q) * B);
    };
    let inGap = Math.abs(x) < hg, entrySide = 0, K_in = K();
    let Kmax = K_in;

    function step() {
      if (exited) return;
      const xPrev = x;
      const ee = Eext(x, y, z);
      kick(Efield(x, y, z, t) + ee[0], dt, ee[1], ee[2]);
      x += vx * dt; y += vy * dt; z += vz * dt; t += dt;
      const k = K(); if (k > Kmax) Kmax = k;

      const nowInGap = Math.abs(x) < hg;
      if (!inGap && nowInGap) {                       // entra al hueco
        entrySide = Math.sign(xPrev); K_in = k;
      } else if (inGap && !nowInGap) {                // sale del hueco
        const side = Math.sign(x);
        if (entrySide === 0) { /* salida inicial desde el centro: no cuenta como cruce completo */ }
        else if (side !== entrySide) crossings.push({ t, K_in, K_out: k, dK: k - K_in, side, y });
        else nBounce++;
        entrySide = 0;
      }
      inGap = nowInGap;

      const dxD = Math.abs(x) - hg;                  // distancia al borde del hueco
      if ((dxD > 0 && dxD * dxD + y * y >= R2) || Math.abs(z) >= zMax) {
        exited = true;
        exitInfo = { t, K: k, rCurv: rCurv(), x, y, z, porTapa: Math.abs(z) >= zMax };
      }
    }

    return {
      params: p, fc, dt, y0, crossings, bhat,
      step,
      setExternal(Q, qx, qy) { Qe = Q; qxe = qx; qye = qy; },   // se puede mover la carga en vivo
      setShield(on) { shield = !!on; },
      runUntil(maxTime, maxSteps) {                   // útil para validar sin render
        let n = 0; maxSteps = maxSteps || Infinity;
        while (!exited && t < maxTime && n++ < maxSteps) step();
        return this.state();
      },
      state: () => {
        const g = gammaNow(), sp = Math.sqrt(vx * vx + vy * vy + vz * vz);
        return {
          t, x, y, z, vx, vy, vz, K: K(), Kmax, rCurv: rCurv(),
          gamma: g, beta: sp / C_LIGHT, speed: sp,
          fcNow: fc / g,                               // frecuencia de giro instantánea (masa relativista)
          crossings: crossings.length, turns: crossings.length / 2, bounces: nBounce,
          exited, exitInfo, Ex: Efield(x, y, z, t), fc, fE: p.fE,
        };
      },
    };
  }

  const api = { E_CHARGE, C_LIGHT, PARTICLES, theory, borisVelocity, rk4Step, createSimulation,
                createCyclotron, suggestDt, vec: { cross, dot, norm } };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.CyclotronPhysics = api;
})(typeof window !== 'undefined' ? window : globalThis);
