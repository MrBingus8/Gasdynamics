"""Prova del 14/09/2021, esercizio 3 (B65-E03): perché la variazione di pressione totale
attraverso un urto obliquo stazionario dipende solo da M_n1 = M1 sin(sigma), se p0 "contiene"
tutta l'energia della corrente (anche la parte cinetica tangenziale).

Quesito teorico: i numeri servono solo come esempio e come verifica delle affermazioni.
Riproduce tutti i numeri e le figure dello svolgimento.
Uso:  python calcoli.py   (dalla cartella dello svolgimento)

Esempio scelto (la traccia non dà numeri): aria, gamma = 1,4, M1 = 5 (come l'es. 2 della stessa
prova) e sigma tale che M_n1 = 2, confrontato con l'urto normale a M = 2 e a M = 5.
Fuori modulo (scritto qui riusando gasdyn): solo il cambio di riferimento galileiano, cioè la
costruzione degli stati visti dall'osservatore che scorre lungo l'urto con V_t, e i grafici
(schema delle velocità, piano T-s, perdite). Tutte le relazioni d'urto vengono da gasdyn.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

QUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(QUI, "..", ".."))
import gasdyn as gd  # noqa: E402

FIG = os.path.join(QUI, "figure")
os.makedirs(FIG, exist_ok=True)

aria = gd.Gas(1.4)
g, d = aria.gamma, aria.delta
k = g / (g - 1)               # esponente gamma/(gamma-1)

# =============================================================================
# 1. ESEMPIO: M1 = 5, M_n1 = 2 (riferimento del laboratorio)
# =============================================================================
M1, Mn1 = 5.0, 2.0
sigma = np.degrees(np.arcsin(Mn1 / M1))
s1 = gd.State.da_M(aria, M1, 0.0, p=1.0, a=1.0, nome="1")     # p1 = 1, a1 = 1
ob = gd.urto_obliquo(s1, sigma=sigma)                          # F §12
s2 = ob.valle
print("=== Urto obliquo, riferimento del laboratorio (F §12) ===")
print(f"M1 = {M1}, M_n1 = {Mn1}: sigma = asin(Mn1/M1) = {sigma:.3f} deg")
print(f"theta = {ob.deflessione:.3f} deg, M2 = {s2.M:.4f}, Mn1 (modulo) = {ob.Mn1:.4f}")
print(f"p2/p1 = {ob.p2_p1:.4f}, T2/T1 = {s2.a**2 / s1.a**2:.4f}")
print(f"p01/p1 = {s1.p0 / s1.p:.2f}, p02/p01 = {s2.p0 / s1.p0:.5f} (modulo: {ob.p02_p01:.5f})")
print(f"a01 = {s1.a0:.6f}, a02 = {s2.a0:.6f}  -> T0 conservata")

# Componenti di velocità rispetto all'urto (la tangenziale si conserva)
sig = np.radians(sigma)
th2 = np.radians(s2.theta)
Vn1, Vt1 = s1.V * np.sin(sig), s1.V * np.cos(sig)
Vn2, Vt2 = s2.V * np.sin(sig - th2), s2.V * np.cos(sig - th2)
print(f"V_t1 = {Vt1:.6f}, V_t2 = {Vt2:.6f}  ([V_t] = {Vt2 - Vt1:.1e})")
print(f"V_n1 = {Vn1:.4f}, V_n2 = {Vn2:.4f}, V_n1/V_n2 = {Vn1 / Vn2:.4f} = rho2/rho1")

# =============================================================================
# 2. CAMBIO DI RIFERIMENTO: osservatore che scorre lungo l'urto con V_t
#    (fuori modulo: stati costruiti a mano con la sola componente normale)
# =============================================================================
n1 = gd.State.da_M(aria, Vn1 / s1.a, 0.0, p=s1.p, a=s1.a, nome="1n")
nr = gd.urto_normale(Mn1, aria)                                 # F §7 / §12
n2 = gd.urto_normale_stato(n1)
n2 = n2.valle if hasattr(n2, "valle") else n2
print("\n=== Riferimento che scorre lungo l'urto con V_t: urto normale a M_n1 ===")
print(f"p2/p1 = {n2.p / n1.p:.4f} (lab {s2.p / s1.p:.4f}); a2 = {n2.a:.6f} (lab {s2.a:.6f});"
      f" rho2 = {n2.rho:.6f} (lab {s2.rho:.6f})")
print(f"p0n1/p1 = {n1.p0 / n1.p:.4f}, p0n2/p0n1 = {n2.p0 / n1.p0:.5f} (lab {s2.p0 / s1.p0:.5f})")
assert abs(n2.p - s2.p) < 1e-9 and abs(n2.rho - s2.rho) < 1e-9
assert abs(n2.p0 / n1.p0 - s2.p0 / s1.p0) < 1e-9

# Fattore comune: p0 = p0n * (T0/T0n)^k, uguale ai due lati
C1 = (s1.a0**2 / n1.a0**2)**k
C2 = (s2.a0**2 / n2.a0**2)**k
print(f"fattore (T0/T0n)^(g/(g-1)): lato 1 = {C1:.4f}, lato 2 = {C2:.4f}")
print(f"T0/T1 = {s1.a0**2:.4f}, T0n/T1 = {n1.a0**2:.4f}, (T0 - T0n)/T1 = d*(Vt/a1)^2 = {d * Vt1**2:.4f}")
assert abs(C1 - C2) < 1e-9

# Entropia: p02/p01 = exp(-ds/R) perché T02 = T01
ds = gd.salto_entropia(s1, s2)
print(f"\nDelta s/R = {ds['ds_R']:.5f}; exp(-Delta s/R) = {np.exp(-ds['ds_R']):.5f}; p02/p01 = {s2.p0 / s1.p0:.5f}")
assert abs(np.exp(-ds["ds_R"]) - s2.p0 / s1.p0) < 1e-9

# =============================================================================
# 3. TABELLA DI CONFRONTO: normale M=2 / obliquo M=5, Mn=2 / normale M=5
# =============================================================================
print("\n=== Confronto ===")
print(" caso                 Mn1   p2/p1   T2/T1   p01/p1   p02/p01   (p01-p02)/p1   ds/R")
casi = [("normale M=2", 2.0, 2.0), ("obliquo M=5, Mn=2", 5.0, 2.0), ("normale M=5", 5.0, 5.0)]
tab = []
for nome, M, Mn in casi:
    r = gd.urto_normale(Mn, aria)
    p01 = gd.isentropiche(M, aria)["p0_p"]
    loss = p01 * (1 - r["p02_p01"])
    tab.append((nome, Mn, r["p2_p1"], r["T2_T1"], p01, r["p02_p01"], loss, r["ds_R"]))
    print(f" {nome:20s} {Mn:4.1f} {r['p2_p1']:7.3f} {r['T2_T1']:7.4f} {p01:8.2f} {r['p02_p01']:9.4f}"
          f" {loss:13.2f} {r['ds_R']:7.4f}")

# Limiti: Mn1 -> 1 (onda di Mach) e Mn1 -> infinito
for Mn in (1.0001, 1.05, 1.2):
    r = gd.urto_normale(Mn, aria)
    print(f"Mn1 = {Mn}: p02/p01 = {r['p02_p01']:.6f}, ds/R = {r['ds_R']:.2e}, "
          f"ds/R / (Mn^2-1)^3 = {r['ds_R'] / (Mn**2 - 1)**3:.4f}")
print(f"coefficiente teorico 2g/(3(g+1)^2) = {2 * g / (3 * (g + 1)**2):.4f}")

# Applicazione: rampa di 11 gradi a M1 = 5 (dati dell'es. 2 della stessa prova)
o11 = gd.urto_obliquo(s1, theta=11.0)
print(f"\nRampa 11 deg a M1 = 5: sigma = {o11.sigma:.2f}, Mn1 = {o11.Mn1:.3f}, p02/p01 = {o11.p02_p01:.4f}")

# Due urti obliqui in serie vs un urto normale (presa d'aria)
o_a = gd.urto_obliquo(s1, theta=11.0)
o_b = gd.urto_obliquo(o_a.valle, theta=11.0)
r_n = gd.urto_normale(o_b.valle.M, aria)
tot = o_a.p02_p01 * o_b.p02_p01 * r_n["p02_p01"]
print(f"Due rampe da 11 deg + urto normale finale: M dopo le rampe = {o_b.valle.M:.3f}, "
      f"p0 recuperata = {tot:.4f} (urto normale diretto a M=5: {gd.urto_normale(5.0, aria)['p02_p01']:.4f})")
print(f"   singoli: {o_a.p02_p01:.4f}, {o_b.p02_p01:.4f}, normale finale {r_n['p02_p01']:.4f}")

# =============================================================================
# 4. FIGURE
# =============================================================================
plt.rcParams.update({"font.size": 10})

# --- Fig. 1: schema delle velocità nei due riferimenti -------------------------
fig, axs = plt.subplots(1, 2, figsize=(9.6, 4.2))
sg = np.radians(40.0)                         # inclinazione solo illustrativa
n_hat = np.array([np.sin(sg), -np.cos(sg)])    # normale all'urto, verso valle
t_hat = np.array([np.cos(sg), np.sin(sg)])     # tangente all'urto
vn1, vt, vn2 = 1.0, 1.15, 0.42                # proporzioni illustrative
V1 = vn1 * n_hat + vt * t_hat
V2 = vn2 * n_hat + vt * t_hat
ax = axs[0]
L = 2.2
ax.plot([-L * t_hat[0], L * t_hat[0]], [-L * t_hat[1], L * t_hat[1]], "r-", lw=2.5)
ax.text(1.55 * t_hat[0] + 0.08, 1.55 * t_hat[1] - 0.05, "urto", color="r")
A = np.array([-0.15, 0.0]) - 1.15 * n_hat - 0.6 * t_hat
P1 = A - V1 * 0.0
ax.annotate("", xy=P1 + V1, xytext=P1, arrowprops=dict(arrowstyle="-|>", lw=2, color="k"))
ax.text(*(P1 + 0.45 * V1 + np.array([-0.15, 0.12])), r"$\mathbf{V}_1$")
ax.annotate("", xy=P1 + vn1 * n_hat, xytext=P1, arrowprops=dict(arrowstyle="-|>", lw=1.3, color="tab:blue"))
ax.annotate("", xy=P1 + V1, xytext=P1 + vn1 * n_hat,
            arrowprops=dict(arrowstyle="-|>", lw=1.3, color="tab:green"))
ax.text(*(P1 + 0.5 * vn1 * n_hat + np.array([-0.45, -0.12])), r"$V_{n1}$", color="tab:blue")
ax.text(*(P1 + vn1 * n_hat + 0.5 * vt * t_hat + np.array([0.12, -0.2])), r"$V_t$", color="tab:green")
P2 = np.array([0.35, 0.0]) + 0.25 * n_hat - 0.1 * t_hat
ax.annotate("", xy=P2 + V2, xytext=P2, arrowprops=dict(arrowstyle="-|>", lw=2, color="k"))
ax.text(*(P2 + 0.55 * V2 + np.array([-0.32, 0.08])), r"$\mathbf{V}_2$")
ax.annotate("", xy=P2 + vn2 * n_hat, xytext=P2, arrowprops=dict(arrowstyle="-|>", lw=1.3, color="tab:blue"))
ax.annotate("", xy=P2 + V2, xytext=P2 + vn2 * n_hat,
            arrowprops=dict(arrowstyle="-|>", lw=1.3, color="tab:green"))
ax.text(*(P2 + vn2 * n_hat + np.array([0.02, -0.22])), r"$V_{n2}$", color="tab:blue")
ax.text(*(P2 + vn2 * n_hat + 0.5 * vt * t_hat + np.array([0.1, -0.1])), r"$V_t$", color="tab:green")
ax.text(-1.6, 0.9, "1", fontsize=14, weight="bold")
ax.text(1.2, -0.9, "2", fontsize=14, weight="bold")
ax.set_title("(a) riferimento del laboratorio:\nurto obliquo")
ax.set_aspect("equal"); ax.set_xlim(-2.0, 2.0); ax.set_ylim(-1.6, 1.8); ax.axis("off")

ax = axs[1]
ax.plot([0, 0], [-1.5, 1.5], "r-", lw=2.5)
ax.text(0.06, 1.3, "urto", color="r")
ax.annotate("", xy=(-0.25, 0.3), xytext=(-0.25 - vn1 * 1.2, 0.3),
            arrowprops=dict(arrowstyle="-|>", lw=2, color="tab:blue"))
ax.text(-1.2, 0.42, r"$V_{n1}$", color="tab:blue")
ax.annotate("", xy=(0.25 + vn2 * 1.2, 0.3), xytext=(0.25, 0.3),
            arrowprops=dict(arrowstyle="-|>", lw=2, color="tab:blue"))
ax.text(0.3, 0.42, r"$V_{n2}$", color="tab:blue")
ax.annotate("", xy=(0.15, -0.5), xytext=(0.15, -1.2),
            arrowprops=dict(arrowstyle="-|>", lw=1.6, color="tab:green"))
ax.text(0.25, -0.95, "l'osservatore scorre\nlungo l'urto con $V_t$", fontsize=9, color="tab:green")
ax.text(-1.4, -0.4, "1", fontsize=14, weight="bold")
ax.text(1.0, 0.8, "2", fontsize=14, weight="bold")
ax.text(-1.55, -1.3, r"$M_{n1}=M_1\sin\sigma$", fontsize=10)
ax.set_title("(b) riferimento che scorre lungo l'urto:\nurto normale a $M_{n1}$")
ax.set_aspect("equal"); ax.set_xlim(-1.7, 1.7); ax.set_ylim(-1.6, 1.8); ax.axis("off")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "riferimenti.pdf"))
plt.close(fig)

# --- Fig. 2: piano T-s con stati statici e di ristagno nei due riferimenti -----
T1 = 1.0
T2 = s2.a**2 / s1.a**2
T0 = s1.a0**2 / s1.a**2
T0n = n1.a0**2 / n1.a**2
dsR = ds["ds_R"]


def isobara(p_rel, sR):
    """T/T1 lungo l'isobara p/p1 = p_rel nel piano T-s: T = p^((g-1)/g) exp(s/cp)."""
    return p_rel**(1 / k) * np.exp(sR / k)


sR = np.linspace(-0.12, dsR + 0.15, 200)
fig, ax = plt.subplots(figsize=(7.2, 5.2))
p01, p02 = s1.p0, s2.p0
p0n1, p0n2 = n1.p0, n2.p0
for p_rel, lab, st in [(1.0, "$p_1$", "-"), (s2.p, "$p_2$", "-"),
                       (p0n1, "$p_{0n,1}$", "--"), (p0n2, "$p_{0n,2}$", "--"),
                       (p01, "$p_{0,1}$", ":"), (p02, "$p_{0,2}$", ":")]:
    ax.plot(sR, isobara(p_rel, sR), "0.6", ls=st, lw=1)
    ax.text(sR[-1] + 0.005, isobara(p_rel, sR[-1]), lab, fontsize=9, va="center")
pts = {"1": (0, T1), "2": (dsR, T2), "0n,1": (0, T0n), "0n,2": (dsR, T0n), "0,1": (0, T0), "0,2": (dsR, T0)}
for nm, (x, y) in pts.items():
    ax.plot(x, y, "ko", ms=5)
    ax.text(x + (-0.045 if x == 0 else 0.012), y + (0.12 if nm != "2" else -0.3), nm, fontsize=10)
ax.annotate("", xy=(0, T0), xytext=(0, T1), arrowprops=dict(arrowstyle="-|>", color="tab:blue"))
ax.annotate("", xy=(dsR, T0), xytext=(dsR, T2), arrowprops=dict(arrowstyle="-|>", color="tab:blue"))
ax.annotate("", xy=(dsR, T2), xytext=(0, T1), arrowprops=dict(arrowstyle="-|>", color="r", lw=1.6))
ax.text(dsR / 2 - 0.04, 1.45, "urto", color="r")
ax.annotate("", xy=(dsR + 0.07, T0), xytext=(dsR + 0.07, T0n),
            arrowprops=dict(arrowstyle="<|-|>", color="tab:green"))
ax.text(dsR + 0.08, (T0 + T0n) / 2, r"$\frac{V_t^2}{2c_p}$", color="tab:green", fontsize=12)
ax.annotate("", xy=(-0.07, T0), xytext=(-0.07, T0n),
            arrowprops=dict(arrowstyle="<|-|>", color="tab:green"))
ax.text(-0.115, (T0 + T0n) / 2, r"$\frac{V_t^2}{2c_p}$", color="tab:green", fontsize=12)
ax.annotate("", xy=(dsR, 0.75), xytext=(0, 0.75), arrowprops=dict(arrowstyle="<|-|>", color="k"))
ax.text(dsR / 2 - 0.03, 0.8, r"$\Delta s$")
ax.set_xlabel(r"$(s-s_1)/R$")
ax.set_ylabel(r"$T/T_1$")
ax.set_xlim(-0.14, dsR + 0.22); ax.set_ylim(0.6, 6.8)
ax.set_title(r"Piano $T$–$s$: $M_1=5$, $M_{n1}=2$ (aria)")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "piano_Ts.pdf"))
plt.close(fig)

# --- Fig. 3: p02/p01 in funzione di sigma e di M_n1 -----------------------------
fig, axs = plt.subplots(1, 2, figsize=(9.6, 4.0))
Mlist = [2.0, 3.0, 5.0, 10.0]
for Mx in Mlist:
    mu = np.degrees(np.arcsin(1 / Mx))
    sg_v = np.linspace(mu + 1e-6, 90.0, 300)
    Mn_v = Mx * np.sin(np.radians(sg_v))
    r_v = np.array([gd.urto_normale(m, aria)["p02_p01"] if m > 1 else 1.0 for m in Mn_v])
    axs[0].plot(sg_v, r_v, label=f"$M_1={Mx:g}$")
    axs[1].plot(Mn_v, r_v, lw=4 if Mx == 10.0 else 1.6, alpha=0.6 if Mx == 10.0 else 1.0,
                label=f"$M_1={Mx:g}$")
axs[0].plot(sigma, s2.p0 / s1.p0, "ko")
axs[0].annotate(r"$M_1=5,\ M_{n1}=2$", (sigma, s2.p0 / s1.p0), (35, 0.85), fontsize=8,
                arrowprops=dict(arrowstyle="-", lw=0.6))
axs[0].set_xlabel(r"$\sigma$ [°]"); axs[0].set_ylabel(r"$p_{0,2}/p_{0,1}$")
axs[0].set_title("(a) in funzione dell'angolo d'urto")
axs[0].set_xlim(0, 90); axs[0].set_ylim(0, 1.02); axs[0].grid(alpha=0.3); axs[0].legend(fontsize=8)
axs[1].plot(2.0, s2.p0 / s1.p0, "ko")
axs[1].set_xlabel(r"$M_{n1}=M_1\sin\sigma$"); axs[1].set_ylabel(r"$p_{0,2}/p_{0,1}$")
axs[1].set_title(r"(b) in funzione di $M_{n1}$: una sola curva")
axs[1].set_xlim(1, 10); axs[1].set_ylim(0, 1.02); axs[1].grid(alpha=0.3); axs[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "rapporto_p0.pdf"))
plt.close(fig)

# --- Fig. 4: perdita assoluta (p01 - p02)/p1 in funzione di M_n1 ---------------
fig, axs = plt.subplots(1, 2, figsize=(9.6, 4.0))
for Mx in Mlist:
    Mn_v = np.linspace(1.0001, Mx, 300)
    p01_v = gd.isentropiche(Mx, aria)["p0_p"]
    loss = np.array([p01_v * (1 - gd.urto_normale(m, aria)["p02_p01"]) for m in Mn_v])
    for ax in axs:
        ax.plot(Mn_v, loss, label=f"$M_1={Mx:g}$")
for ax in axs:
    ax.plot([2.0, 2.0], [tab[0][6], tab[1][6]], "ko", ms=4)
    ax.set_xlabel(r"$M_{n1}$"); ax.set_ylabel(r"$(p_{0,1}-p_{0,2})/p_1$")
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)
axs[0].set_title("(a) scala lineare"); axs[0].set_xlim(1, 5); axs[0].set_ylim(0, 600)
axs[1].set_yscale("log"); axs[1].set_title("(b) scala logaritmica"); axs[1].set_xlim(1, 10)
axs[1].set_ylim(1e-3, 1e5)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "perdita_assoluta.pdf"))
plt.close(fig)

print("\nFigure salvate in", FIG)
