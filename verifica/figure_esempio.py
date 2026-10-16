"""Genera le figure di esempio in verifica/figure/ (controllo visivo dei grafici di gasdyn)."""
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RADICE = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RADICE))
import gasdyn as gd  # noqa: E402

OUT = RADICE / "verifica" / "figure"
OUT.mkdir(parents=True, exist_ok=True)
aria = gd.Gas(1.4)


def salva(fig, nome):
    fig.savefig(OUT / nome, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(OUT / nome)


# 1. Sod: x–t, p–u, profili
sod = gd.riemann(gd.State(aria, 1.0, 0.0, 1.0), gd.State(aria, 0.125, 0.0, 0.1))
fig, ax = plt.subplots(figsize=(6, 4.5))
gd.grafico_xt(sod, 0.25, (-0.5, 0.5), ax=ax, titolo=f"Sod: {sod.configurazione}")
salva(fig, "sod_xt.png")
fig, ax = plt.subplots(figsize=(6, 4.5))
gd.grafico_pu(sod, ax=ax, titolo="Sod: piano p–u")
salva(fig, "sod_pu.png")
salva(gd.grafico_profili(sod, 0.2, (-0.5, 0.5)), "sod_profili.png")

# 2. Ventaglio transonico e vuoto
tr = gd.riemann(gd.State(aria, 1.0, 0.75, 1.0), gd.State(aria, 0.125, 0.0, 0.1))  # Toro, test 1
fig, ax = plt.subplots(figsize=(6, 4.5))
gd.grafico_xt(tr, 0.3, (-0.8, 0.8), ax=ax, titolo=f"ventaglio transonico: {tr.configurazione}")
salva(fig, "transonico_xt.png")
vu = gd.riemann(gd.State(aria, 1.0, -7.0, 1.0), gd.State(aria, 1.0, 7.0, 1.0))
fig, ax = plt.subplots(figsize=(6, 4.5))
gd.grafico_xt(vu, 0.1, (-1, 1), ax=ax, titolo=f"vuoto: {vu.configurazione}")
salva(fig, "vuoto_xt.png")

# 3. B51-E01: incontro di due correnti, piano p–θ e odografo
sotto = gd.State.da_M(gd.Gas(1.4, 287.0), 6.0, p=0.01, T=260.0, nome="inf")
sopra = gd.State.da_M(gd.Gas(1.4, 287.0), 4.0, p=0.50, T=350.0, nome="sup")
r = gd.riemann_stazionario(sopra, sotto)
fig, ax = plt.subplots(figsize=(6, 4.5))
gd.grafico_p_theta([(sopra, "C+"), (sotto, "C-")],
                   [sopra, sotto, r.onda_sopra.valle.con(nome="sup'=inf'")],
                   ax=ax, log=True, p_min_rel=0.02, titolo="B51-E01: piano p–θ")
salva(fig, "b51_p_theta.png")
fig, ax = plt.subplots(figsize=(5.5, 5))
gd.grafico_odografo([(sopra, "C+"), (sotto, "C-")], [sopra, sotto], ax=ax, titolo="B51-E01: odografo")
salva(fig, "b51_odografo.png")

# 4. Riflessione regolare M = 2.8, θw = 10°
st1 = gd.State.da_M(aria, 2.8, p=1.0, a=1.0, nome="1")
rr = gd.riflessione_regolare(st1, 10.0)
fig, ax = plt.subplots(figsize=(6, 4.5))
gd.grafico_p_theta([(st1, "entrambe"), (rr["2"], "C+")], [rr["1"], rr["2"], rr["3"]], ax=ax,
                   titolo="riflessione regolare M1 = 2.8, θw = 10°")
salva(fig, "rr_p_theta.png")

# 5. Diagramma θ–M1 di transizione RR/MR
import numpy as np  # noqa: E402

d = gd.diagramma_rr_mr(np.linspace(1.5, 6, 30), aria)
fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(d["M"], d["theta_D"], "k", label="θ_D (distacco)")
ax.plot(d["M"], d["theta_N"], "k--", label="θ_N (von Neumann)")
ax.set_xlabel("M1")
ax.set_ylabel("θw [°]")
ax.legend()
ax.grid(alpha=0.25)
salva(fig, "diagramma_rr_mr.png")

# 6. Ugello di lunghezza minima Me = 2.5
u = gd.ugello_minima_lunghezza(2.5, aria, n=20)
fig, ax = plt.subplots(figsize=(7, 3))
ax.plot(u["x_parete"], u["y_parete"], "k", lw=2)
for (k, j), (x, y, th, nu) in u["nodi"].items():
    ax.plot(x, y, ".", color="0.5", ms=2)
ax.set_aspect("equal")
ax.set_title(f"ugello di lunghezza minima Me = 2.5: L/h* = {u['L_su_h_gola']:.3f}")
salva(fig, "ugello_minima_lunghezza.png")

# 7. Pistone con legge parabolica (B49-E01) e formazione dell'urto
q = gd.State(aria, 1.0, 0.0, 1.0)
P = gd.pistone_legge(lambda t: 0.74 * t**2, q, 1.0)
fig, ax = plt.subplots(figsize=(6, 4.5))
for tp in np.linspace(0, 1, 15):
    x0, t0, lam = P.caratteristica(tp)
    tt = np.array([t0, 1.2])
    ax.plot(x0 + lam * (tt - t0), tt, color="0.5", lw=0.7)
tt = np.linspace(0, 1.2, 100)
ax.plot(np.where(tt < 1, 0.74 * tt**2, 0.74 + 1.48 * (tt - 1)), tt, **gd._STILE["pistone"])
ax.plot(P.x_bu, P.t_bu, "ro")
ax.set_xlim(0, 2.5)
ax.set_ylim(0, 1.2)
ax.set_xlabel("x")
ax.set_ylabel("t")
ax.set_title(f"B49-E01: t_BU = {P.t_bu:.5f}, x_BU = {P.x_bu:.5f}")
salva(fig, "pistone_b49.png")
