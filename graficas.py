"""Genera las gráficas de validación a partir de ./datos (creados con export_datos.js).
Uso:  python graficas.py      (requiere: pip install matplotlib numpy)"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("graficas", exist_ok=True)
plt.rcParams.update({"figure.dpi": 130, "axes.grid": True, "grid.alpha": 0.3, "font.size": 10})


def leer(nombre):
    with open(f"datos/{nombre}", newline="", encoding="utf-8") as f:
        filas = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) if k != "particula" else r[k] for r in filas]) if k != "particula"
            else [r[k] for r in filas] for k in filas[0]}


def guardar(fig, nombre):
    fig.tight_layout()
    fig.savefig(f"graficas/{nombre}.png")
    plt.close(fig)
    print(f"graficas/{nombre}.png")


# 1) Espiral del protón, coloreada por energía
d = leer("trayectoria_proton.csv")
fig, ax = plt.subplots(figsize=(6, 6))
sc = ax.scatter(d["x_m"], d["y_m"], c=d["K_MeV"], s=2, cmap="coolwarm")
ax.add_patch(plt.Circle((0, 0), 0.5, fill=False, ls="--", color="gray"))
ax.axvspan(-0.005, 0.005, color="gold", alpha=0.4, label="hueco")
ax.set_aspect("equal"); ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)")
ax.set_title("Trayectoria del protón (B = 1 T, V = 50 kV)")
fig.colorbar(sc, label="K (MeV)"); ax.legend(loc="upper right")
guardar(fig, "01_espiral")

# 2) K vs cruces y ganancia por cruce (V3, V4)
d = leer("k_vs_cruces.csv")
n = np.arange(1, len(d["K_MeV"]) + 1)
pend, ord_ = np.polyfit(n, d["K_MeV"], 1)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4))
a1.plot(n, d["K_MeV"], ".", ms=3, label="simulación")
a1.plot(n, pend * n + ord_, "r-", lw=1, label=f"ajuste: {pend * 1e3:.2f} keV/cruce (teoría 50)")
a1.set_xlabel("cruces del hueco"); a1.set_ylabel("K (MeV)"); a1.set_title("V4: energía vs. cruces"); a1.legend()
a2.plot(n, d["dK_sobre_qV"], ".", ms=3); a2.axhline(1, color="r", lw=1)
a2.set_xlabel("cruce"); a2.set_ylabel("ΔK / qV"); a2.set_title("V3: ganancia por cruce")
guardar(fig, "02_energia_vs_cruces")

# 3) Barrido de frecuencia (V6)
d = leer("barrido_frecuencia.csv")
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(d["fE_sobre_fc"], 100 * d["Kmax_sobre_Kteo"], "o-", ms=3)
ax.axvline(1, color="r", ls="--", lw=1)
ax.set_xlabel("f_E / f_c"); ax.set_ylabel("K máx alcanzada / K teórica (%)"); ax.set_title("V6: resonancia")
guardar(fig, "03_resonancia")

# 4) Barrido de fase (V7)
d = leer("barrido_fase.csv")
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(d["phi_sobre_pi"], 100 * d["ideal_K_sobre_Kteo"], "o-", ms=3, label="campo ideal")
ax.plot(d["phi_sobre_pi"], 100 * d["sinus_K_sobre_Kteo"], "s-", ms=3, label="campo sinusoidal")
ax.set_xlabel("φ / π"); ax.set_ylabel("K máx / K teórica (%)"); ax.set_title("V7: fase inicial (inyección en la D)"); ax.legend()
guardar(fig, "04_fase")

# 5) Boris vs RK4 (V8)
d = leer("boris_vs_rk4.csv")
fig, ax = plt.subplots(figsize=(7, 4))
ax.semilogy(d["vueltas"], d["deriva_boris"], label="Boris")
ax.semilogy(d["vueltas"], d["deriva_rk4"], label="RK4")
ax.set_xlabel("vueltas"); ax.set_ylabel("deriva relativa de energía"); ax.set_title("V8: conservación de energía (dt = T/50)"); ax.legend()
guardar(fig, "05_boris_vs_rk4")

# 6) Sensibilidad a dt (V9)
d = leer("sensibilidad_dt.csv")
fig, ax = plt.subplots(figsize=(6, 4))
ax.loglog(d["pasos_por_periodo"], d["error_radio_max_pct"], "o-")
ax.set_xlabel("pasos por período (T/N)"); ax.set_ylabel("error de radio máx (%)"); ax.set_title("V9: sensibilidad a dt (orden 2)")
guardar(fig, "06_sensibilidad_dt")

# ---------------------------------------------------------------------------------------------
# CARGA EXTERNA
# 7) K final vs Q, con y sin blindaje de las D
d = leer("carga_vs_Q.csv")
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(d["Q_nC"], 100 * d["sin_blindaje_K_sobre_Kteo"], "o-", ms=4, label="sin blindaje (campo de Coulomb en todo el plano)")
ax.plot(d["Q_nC"], 100 * d["con_blindaje_K_sobre_Kteo"], "s-", ms=4, label="con blindaje (Eext = 0 dentro de las D)")
ax.set_xscale("symlog", linthresh=10)
ax.set_xlabel("Q (nC), carga en (0.75, 0) m"); ax.set_ylabel("K final / K teórica (%)")
ax.set_title("C: energía final con una carga externa cerca"); ax.legend(fontsize=8)
guardar(fig, "07_carga_vs_Q")

# 8) K final vs distancia
d = leer("carga_vs_distancia.csv")
fig, ax = plt.subplots(figsize=(7, 4))
ax.semilogx(d["d_m"], 100 * d["sin_blindaje_K_sobre_Kteo"], "o-", ms=4, label="sin blindaje")
ax.semilogx(d["d_m"], 100 * d["con_blindaje_K_sobre_Kteo"], "s-", ms=4, label="con blindaje")
ax.set_xlabel("distancia de la carga al centro (m), Q = +1 µC"); ax.set_ylabel("K final / K teórica (%)")
ax.set_title("C: efecto vs. distancia"); ax.legend()
guardar(fig, "08_carga_vs_distancia")

# 9) Deriva E×B: simulada vs teórica
d = leer("deriva_ExB.csv")
fig, ax = plt.subplots(figsize=(5, 5))
m = max(d["vd_teorica_ms"].max(), d["vd_medida_ms"].max()) * 1.05
ax.plot([0, m], [0, m], "r--", lw=1, label="y = x (teoría)")
ax.plot(d["vd_teorica_ms"], d["vd_medida_ms"], "o")
ax.set_aspect("equal"); ax.set_xlabel("v_d teórica = E/B (m/s)"); ax.set_ylabel("v_d medida en la simulación (m/s)")
ax.set_title("C4: deriva E×B del centro de giro"); ax.legend()
guardar(fig, "09_deriva_ExB")

# 10) Trayectorias con y sin carga (un panel por caso para que no se tapen)
with open("datos/trayectorias_carga.csv", newline="", encoding="utf-8") as f:
    filas = list(csv.DictReader(f))
casos = [("sin_carga", "sin carga", "gray"), ("Q1uC_sin_blindaje", "Q = 1 µC, sin blindaje", "tab:red"),
         ("Q1uC_con_blindaje", "Q = 1 µC, con blindaje", "tab:green")]
fig, axs = plt.subplots(1, 3, figsize=(13, 4.6), sharex=True, sharey=True)
for ax, (clave, etq, col) in zip(axs, casos):
    xs = [float(r["x_m"]) for r in filas if r["caso"] == clave]; ys = [float(r["y_m"]) for r in filas if r["caso"] == clave]
    ax.plot(xs, ys, color=col, lw=0.35)
    ax.add_patch(plt.Circle((0, 0), 0.5, fill=False, ls="--", color="k", lw=0.8))
    ax.plot([0.75], [0], "o", color="tab:red" if clave != "sin_carga" else "none", ms=7)
    ax.set_aspect("equal"); ax.set_xlabel("x (m)"); ax.set_title(etq, fontsize=10)
    ax.set_xlim(-0.6, 0.85); ax.set_ylim(-0.6, 0.6)
axs[0].set_ylabel("y (m)")
fig.suptitle("Espiral del protón con una carga externa en (0.75, 0): deriva E×B y blindaje")
guardar(fig, "10_trayectorias_carga")
