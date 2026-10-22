"""Prova del 13/02/2023, esercizio 2 (B70-E02): condizione al contorno d'uscita
((γ+1)/2) M^(4/(γ+1)) + M^(-2(γ-1)/(γ+1)) = K, con M = u/a nella sezione d'uscita, e suo
effetto su un'onda semplice acustica che arriva dall'interno del tubo.

Riproduce tutti i numeri e le figure dello svolgimento.
Uso:  python calcoli.py   (dalla cartella dello svolgimento)

Parti fuori modulo (scritte qui riusando gasdyn): studio della funzione f(M) della condizione al
contorno (minimo, radici) e coefficiente di riflessione; gli stati 0, 1, 2 dell'esempio sono
costruiti con gd.onda_semplice (F §4).
Il caso numerico (aria, K = 1,8, compressione incidente p1/p0 = 1,10) è solo illustrativo:
la traccia non dà numeri.
"""
import os
import sys

import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

QUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(QUI, "..", ".."))
import gasdyn as gd  # noqa: E402

FIG = os.path.join(QUI, "figure")
os.makedirs(FIG, exist_ok=True)

aria = gd.Gas(1.4)            # adimensionale: a0 = 1
g, d = aria.gamma, aria.delta
n_esp = 4 / (g + 1)           # esponente di M nel primo termine
m_esp = 2 * (g - 1) / (g + 1)  # esponente di 1/M nel secondo termine


def f(M):
    """Primo membro della condizione al contorno, funzione del solo M = u/a."""
    return (g + 1) / 2 * M**n_esp + M**(-m_esp)


def df(M):
    """f'(M) = 2 M^(-m-1) [M^2 - (γ-1)/(γ+1)]   (n + m = 2)."""
    return 2 * M**(-m_esp - 1) * (M**2 - (g - 1) / (g + 1))


def rapporto(M):
    """ΔR1/ΔR2 al contorno = Δa_rifl/Δa_inc = (1 - δM)/(1 + δM) (M costante, F §3)."""
    return (1 - d * M) / (1 + d * M)


def rapporto_gradienti(M):
    """a_x,rifl / a_x,inc = -(1+M)(1-δM)/((1-M)(1+δM)) (forma differenziale)."""
    return -(1 + M) * (1 - d * M) / ((1 - M) * (1 + d * M))


# =============================================================================
# 1. STUDIO DI f(M)
# =============================================================================
print(f"esponenti: n = 4/(γ+1) = {n_esp:.4f}, m = 2(γ-1)/(γ+1) = {m_esp:.4f}, n + m = {n_esp + m_esp:.4f}")
Mstar = np.sqrt((g - 1) / (g + 1))
fmin = f(Mstar)
fmin_formula = (g + 1) / 2 * ((g + 1) / (g - 1))**((g - 1) / (g + 1))
f1 = f(1.0)
print(f"minimo in M* = sqrt((γ-1)/(γ+1)) = {Mstar:.4f}, f(M*) = {fmin:.4f} (formula {fmin_formula:.4f})")
print(f"f(1) = (γ+3)/2 = {f1:.4f}")
# controllo numerico della derivata
for Mt in (0.2, 0.6, 0.9):
    h = 1e-6
    assert abs((f(Mt + h) - f(Mt - h)) / (2 * h) - df(Mt)) < 1e-6
assert abs(fmin - fmin_formula) < 1e-12

K = 1.8
Ma = brentq(lambda M: f(M) - K, 1e-8, Mstar)
Mb = brentq(lambda M: f(M) - K, Mstar, 1.0)
print(f"\nK = {K}: radici subsoniche M_a = {Ma:.4f}, M_b = {Mb:.4f}")
for Mx in (Ma, Mb):
    print(f"   M = {Mx:.3f}: Δa_r/Δa_i = {rapporto(Mx):.4f}, a_x,r/a_x,i = {rapporto_gradienti(Mx):.4f}")

print("\nTabella del rapporto (aria):")
for Mx in (0.0, 0.2, Mstar, 0.6, 0.8, 1.0):
    gr = rapporto_gradienti(Mx) if Mx < 1 else float("-inf")
    print(f"   M = {Mx:.3f}: Δa_r/Δa_i = {rapporto(Mx):.4f}, a_x,r/a_x,i = {gr:.3f}")
print(f"rapporto nullo solo per M = 1/δ = {1/d:.1f} (non subsonico)")

# =============================================================================
# 2. ESEMPIO: M_e = M_b, compressione incidente p1/p0 = 1,10
# =============================================================================
Me = Mb
s0 = gd.State.da_pa(aria, 1.0, 1.0, u=Me, nome="0")
s1 = gd.onda_semplice(s0, "u+a", p=1.10, nome="1")          # R1 uniforme
a2 = s1.R2 / (1 + d * Me)                                    # R2 uniforme + u2 = Me a2
s2 = gd.onda_semplice(s1, "u-a", a=a2, nome="2")             # R2 uniforme
assert abs(s2.u / s2.a - Me) < 1e-12
assert abs(s0.R1 - rapporto(Me) * s0.R2) < 1e-12           # stato 0 compatibile col contorno
print("\n=== Esempio (stato 0 di riferimento: p = 1, a = 1) ===")
for s in (s0, s1, s2):
    print(f"stato {s.nome}: a = {s.a:.4f}, u = {s.u:.4f}, p = {s.p:.4f}, M = {s.u/s.a:.4f}, "
          f"R1 = {s.R1:.4f}, R2 = {s.R2:.4f}, u-a = {s.u - s.a:.4f}, u+a = {s.u + s.a:.4f}")
dai, dar = s1.a - s0.a, s2.a - s1.a
dpi, dpr = s1.p - s0.p, s2.p - s1.p
print(f"Δa_i = {dai:.4f}, Δa_r = {dar:.4f}, rapporto = {dar/dai:.4f} (formula {rapporto(Me):.4f})")
print(f"Δp_i = {dpi:.4f}, Δp_r = {dpr:.4f}, rapporto = {dpr/dpi:.4f} (acustico {rapporto(Me):.4f})")
print(f"Δu_i = {s1.u - s0.u:.4f}, Δu al contorno = {s2.u - s0.u:.4f}, "
      f"formula 2δM/(1+δM) Δu_i = {2*d*Me/(1+d*Me)*(s1.u - s0.u):.4f}")
# confronto: uscita a pressione imposta (p = p0 => a = a0, F §9)
s2_aperto = gd.onda_semplice(s1, "u-a", a=s0.a, nome="2'")  # p imposta = p0
print(f"confronto, p imposta (estremo aperto): a = {s2_aperto.a:.4f}, u = {s2_aperto.u:.4f}, "
      f"Δa_r/Δa_i = {(s2_aperto.a - s1.a)/dai:.3f}")

# =============================================================================
# 3. FIGURE
# =============================================================================
plt.rcParams.update({"font.size": 10})

# --- Figura 1: f(M) e rapporto di riflessione ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.9))
M = np.linspace(0.03, 1.15, 400)
ax1.plot(M, f(M), "k", lw=1.6)
ax1.axhspan(1.0, fmin, color="0.88")
ax1.text(0.86, 1.22, "1 < K < f(M*):\nnessuno stato d'uscita", fontsize=8.5, ha="center")
ax1.axhline(K, color="C0", ls="--", lw=1)
ax1.plot([Ma, Mb], [K, K], "o", color="C0")
ax1.annotate(f"$M_a$ = {Ma:.2f}", (Ma, K), xytext=(Ma + 0.03, K + 0.22), fontsize=9, color="C0")
ax1.annotate(f"$M_b$ = {Mb:.2f}", (Mb, K), xytext=(Mb - 0.15, K + 0.22), fontsize=9, color="C0")
ax1.text(0.04, K + 0.04, f"K = {K}", color="C0", fontsize=9)
ax1.plot(Mstar, fmin, "s", color="C3")
ax1.annotate(f"M* = {Mstar:.3f}\nf = {fmin:.2f}", (Mstar, fmin), xytext=(Mstar - 0.30, fmin - 0.30),
             fontsize=8.5, color="C3")
ax1.plot(1, f1, "^", color="k")
ax1.annotate(f"M = 1, f = {f1:.1f}", (1, f1), xytext=(0.80, f1 + 0.15), fontsize=8.5)
ax1.axvline(1, color="0.5", lw=0.8, ls=":")
ax1.set_xlim(0, 1.15)
ax1.set_ylim(1.0, 3.2)
ax1.set_xlabel("M = u/a nella sezione d'uscita")
ax1.set_ylabel("f(M)")
ax1.set_title("(a) primo membro della condizione, γ = 1,4")

Mr = np.linspace(0, 0.95, 300)
ax2.plot(Mr, rapporto(Mr), "C0", lw=1.6, label=r"$\Delta a_r/\Delta a_i=(1-\delta M)/(1+\delta M)$")
ax2.plot(Mr, -rapporto_gradienti(Mr), "C2", lw=1.3, ls="-.",
         label=r"$-a_{x,r}/a_{x,i}$ (gradienti)")
ax2.axhline(1, color="k", lw=0.9, ls="--")
ax2.axhline(-1, color="C3", lw=0.9, ls="--")
ax2.axhline(0, color="0.5", lw=0.6)
ax2.text(0.02, 1.05, "parete (u = 0)", fontsize=8.5)
ax2.text(0.02, -0.93, "estremo aperto (p = cost.)", fontsize=8.5, color="C3")
ax2.text(0.02, 0.05, "contorno non riflettente", fontsize=8.5, color="0.4")
ax2.set_xlim(0, 1)
ax2.set_ylim(-1.2, 3.0)
ax2.set_xlabel(r"$M_e$")
ax2.set_ylabel("rapporto riflessa / incidente")
ax2.set_title("(b) intensità dell'onda riflessa")
ax2.legend(fontsize=8, loc="upper center")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fM_rapporto.pdf"))
plt.close(fig)

# --- Figura 2: piano u-a ---
fig, ax = plt.subplots(figsize=(6.2, 4.6))
ax.plot([0, 1.0], [0, 1.0 / Me], "C0", lw=1.6, label=r"contorno: $u = M_e\,a$")
ug = np.linspace(0.6, 1.0, 50)
ax.plot(ug, s0.R1 + d * ug, "k", lw=1, ls="--", label=r"onda $u+a$: $a-\delta u=R_{1,0}$")
ax.plot(ug, s1.R2 - d * ug, "C2", lw=1, ls="-.", label=r"onda $u-a$: $a+\delta u=R_{2,1}$")
ax.axhline(s0.a, color="C3", lw=0.9, ls=":", label=r"estremo aperto: $p=p_0$ ($a=a_0$)")
for s, off in ((s0, (6, -12)), (s1, (6, 2)), (s2, (-16, 6))):
    ax.plot(s.u, s.a, "ko", ms=5)
    ax.annotate(s.nome, (s.u, s.a), xytext=off, textcoords="offset points", fontsize=11)
ax.plot(s2_aperto.u, s2_aperto.a, "o", color="C3", ms=5)
ax.annotate("2'", (s2_aperto.u, s2_aperto.a), xytext=(4, -13), textcoords="offset points",
            fontsize=10, color="C3")
ax.set_xlim(0.6, 0.95)
ax.set_ylim(0.95, 1.07)
ax.set_xlabel("u / a0")
ax.set_ylabel("a / a0")
ax.set_title(f"Piano u–a nella sezione d'uscita (M_e = {Me:.2f}, p1/p0 = 1,10)")
ax.legend(fontsize=8, loc="lower right")
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "piano_ua.pdf"))
plt.close(fig)

# --- Figura 3: piano x-t qualitativo (caratteristiche semirette con le pendenze degli stati uniformi) ---
# onda acustica (p1/p0 = 1,02), così le caratteristiche non si incontrano nel tratto disegnato
s1 = gd.onda_semplice(s0, "u+a", p=1.02, nome="1")
s2 = gd.onda_semplice(s1, "u-a", a=s1.R2 / (1 + d * Me), nome="2")
L = 1.0
fig, ax = plt.subplots(figsize=(6.4, 4.6))
x_head0, x_tail0 = 0.0, -0.45      # fronti dell'onda incidente a t = 0 (onda che arriva da sinistra)
lam_h, lam_t = s0.u + s0.a, s1.u + s1.a
tA = (L - x_head0) / lam_h
tB = (L - x_tail0) / lam_t
n = 6
for k in range(n):
    x0 = x_head0 + (x_tail0 - x_head0) * k / (n - 1)
    lam = lam_h + (lam_t - lam_h) * k / (n - 1)
    tL = (L - x0) / lam
    t0 = max(0.0, -x0 / lam)
    ax.plot([x0 + lam * t0, L], [t0, tL], "k", lw=1.2 if k in (0, n - 1) else 0.6)
# onda riflessa u-a: dal contorno tra tA e tB; testa con u1-a1, coda con u2-a2
mu_h, mu_t = s1.u - s1.a, s2.u - s2.a
tmax = 4.6
for k in range(n):
    tb = tA + (tB - tA) * k / (n - 1)
    lam = mu_h + (mu_t - mu_h) * k / (n - 1)
    x_end = L + lam * (tmax - tb)
    xe = max(0.0, x_end)
    te = tb + (xe - L) / lam
    ax.plot([L, xe], [tb, te], "C2", lw=1.2 if k in (0, n - 1) else 0.6)
# zona d'interazione (due famiglie) vicino all'uscita
tI = (lam_t * tB - mu_h * tA) / (lam_t - mu_h)   # testa riflessa ∩ coda incidente
xI = L + mu_h * (tI - tA)
ax.fill([L, L, xI], [tA, tB, tI], color="0.75", zorder=0)
ax.annotate("zona a 2 famiglie", (xI + 0.01, 0.5 * (tA + tB)), xytext=(0.55, 1.35), fontsize=8,
            arrowprops=dict(arrowstyle="->", lw=0.6))

ax.axvline(L, color="C0", lw=2.5)
ax.text(L + 0.015, 0.05, "uscita\n$M=M_e$", color="C0", fontsize=9)
ax.text(0.55, 0.12, "0", fontsize=13)
ax.text(0.15, 1.6, "1", fontsize=13)
ax.text(0.80, 3.6, "2", fontsize=13)
ax.text(0.22, 0.62, "$u+a$ (compressione)", fontsize=8.5)
ax.text(0.05, 3.1, "$u-a$ riflessa\n(compressione)", fontsize=8.5, color="C2")
ax.set_xlim(0, L + 0.17)
ax.set_ylim(0, tmax)
ax.set_xlabel("x / L")
ax.set_ylabel(r"$a_0 t / L$")
ax.set_title("Piano x–t (qualitativo)")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "xt.pdf"))
plt.close(fig)
print(f"\nx-t: tA = {tA:.3f}, tB = {tB:.3f}, pendenze riflessa testa {mu_h:.4f}, coda {mu_t:.4f}")
print("figure salvate in", FIG)
