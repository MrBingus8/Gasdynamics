"""Prova del 14/09/2021, esercizio 1 (B65-E01): problema 1D instazionario con dati iniziali
isentropici, espansione centrata in x = 0 e compressione lineare in 0 <= x <= 1.

Riproduce tutti i numeri e le figure dello svolgimento.
Uso:  python calcoli.py   (dalla cartella dello svolgimento)

Dati: x <= 0: p = 1, u = 0, a = 1,522;  0 <= x <= 1: u = 2,61 - 1,856 x, a = 1 + 0,4995 x;
x >= 1: u = 0,754, a = 1,4995.  gamma = 1,4 assunto (unico valore per cui a + delta u e'
continuo in x = 0, come richiede un'onda semplice u-a).

Parti fuori modulo (scritte qui riusando gasdyn):
- invarianti R1, R2 dei dati iniziali e velocita' caratteristiche in funzione degli invarianti;
- stima del punto di separazione delle due famiglie con pendenze medie (F §9, reticolo);
- profili a t* + eps (ventaglio con F §5, onda u+a non centrata stimata per interpolazione);
- curve qualitative di urto e contatto dopo t* nel piano x-t (solo andamento);
- controllo con un reticolo delle caratteristiche completo della zona a due famiglie
  (solo per la sintesi in chat: NON compare nel PDF).
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
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3})

aria = gd.Gas(1.4)
g, d = aria.gamma, aria.delta
esp = 2 * g / (g - 1)          # 2γ/(γ−1) = 7


def riga(titolo):
    print("\n" + "=" * 78 + "\n" + titolo + "\n" + "=" * 78)


# =============================================================================
# 1. DATI E INVARIANTI (F §3, §4)
# =============================================================================
riga("1. Dati iniziali e invarianti di Riemann")
u_lin = lambda x: 2.61 - 1.856 * x
a_lin = lambda x: 1 + 0.4995 * x

S1 = gd.State.da_pa(aria, 1.0, 1.522, u=0.0).con(nome="1")
S2 = gd.isentropica(S1, a=1.0).con(u=2.61, nome="2")          # x = 0+
S4 = gd.isentropica(S1, a=1.4995).con(u=0.754, nome="4")      # x >= 1
for s in (S1, S2, S4):
    print(f"stato {s.nome}: p = {s.p:.4f}, u = {s.u:.4f}, a = {s.a:.4f}, M = {s.M:.4f}, "
          f"R1 = a-δu = {s.R1:.4f}, R2 = a+δu = {s.R2:.4f}, s = {s.s:.5f}")
print(f"continuita' in x = 1: u = {u_lin(1):.4f}, a = {a_lin(1):.4f}")
print(f"γ tale che a+δu sia continuo in x = 0: δ = (1.522-1)/2.61 = {(1.522 - 1) / 2.61:.5f}")
print(f"γ tale che a+δu sia uniforme in [0,1]: δ = 0.4995/1.856 = {0.4995 / 1.856:.5f}  (diverso)")
# In [0,1]: R1 = 0.478 + 0.8707 x, R2 = 1.522 + 0.1283 x
dR1 = 0.4995 + d * 1.856
dR2 = 0.4995 - d * 1.856
print(f"in [0,1]: R1 = {S2.R1:.4f} + {dR1:.4f} x,  R2 = {S2.R2:.4f} + {dR2:.4f} x  -> entrambe le famiglie")
print(f"          u-a = {u_lin(0) - a_lin(0):.4f} - {1.856 + 0.4995:.4f} x,  "
      f"u+a = {u_lin(0) + a_lin(0):.4f} - {1.856 - 0.4995:.4f} x")


def R1i(x):
    return S2.R1 + dR1 * x


def R2i(x):
    return S2.R2 + dR2 * x


def lam1(r1, r2):
    """u - a in funzione degli invarianti: u = (R2-R1)/(2δ), a = (R1+R2)/2 (F §3)."""
    return (r2 - r1) / (2 * d) - (r1 + r2) / 2


def lam2(r1, r2):
    return (r2 - r1) / (2 * d) + (r1 + r2) / 2


print(f"con δ = 0.2: u-a = 2R2 - 3R1 -> {lam1(0.0, 1.0):.3f} R2 {lam1(1.0, 0.0):+.3f} R1;  "
      f"u+a = 3R2 - 2R1 -> {lam2(0.0, 1.0):.3f} R2 {lam2(1.0, 0.0):+.3f} R1")

# =============================================================================
# 2. ESPANSIONE CENTRATA IN x = 0 (F §5)
# =============================================================================
riga("2. Discontinuita' in x = 0: espansione centrata u-a")
E = gd.onda_verso(S1, "u-a", u=S2.u)
print(E)
print(f"testa u1-a1 = {S1.u - S1.a:.4f}, coda u2-a2 = {S2.u - S2.a:.4f}, transonica = {E.transonico}")
xi_son = 0.0
print(f"punto sonico: ξ = 0, u = a = R2/(1+δ) = {S1.R2 / (1 + d):.4f}, p = {(S1.R2 / (1 + d) / S1.a)**esp:.4f}")

# =============================================================================
# 3. FORMAZIONE DELL'URTO (F §6)
# =============================================================================
riga("3. Formazione dell'urto dai dati iniziali")
fu = gd.formazione_urto_ci(u_lin, a_lin, aria, (0.0, 1.0))
for fam in ("u-a", "u+a"):
    r = fu[fam]
    print(f"{fam}: v' = {r['dv_min']:.4f}, t_BU = {r['t_bu']:.5f}, x_BU = {r['x_bu']:.4f}, "
          f"onda semplice = {r['onda_semplice']}")
t_s = 0.4246                                   # istante della traccia
x_s = (S2.u - S2.a) * t_s                      # tutte le u-a convergono: x = x0 + (1.61 - 2.3555 x0) t
print(f"compressione centrata: t* = 1/2.3555 = {1 / 2.3555:.5f} (traccia 0.4246), x* = 1.61 t* = {x_s:.4f}")
x0s = np.linspace(0, 1, 5)
print("verifica: x(t*) delle u-a da x0 =", np.round(x0s + (u_lin(x0s) - a_lin(x0s)) * t_s, 4))
# con le pendenze dopo la separazione (R2 = 1.522): u-a = 1.61 - 2.612 x0
gr_dopo = 3 * dR1
print(f"pendenze dopo la separazione: u-a = 2R2,1 - 3R1(x0) -> gradiente {gr_dopo:.4f}, "
      f"t = {1 / gr_dopo:.4f}, x = {(S2.u - S2.a) / gr_dopo:.4f}")

# =============================================================================
# 4. STATO 3 TRA LE FAMIGLIE SEPARATE (F §3)
# =============================================================================
riga("4. Stato 3: R2 dello stato 1 e R1 dello stato 4")
S3 = gd.punto_caratteristiche(S1, S4).con(nome="3")
print(f"a3 = (R2,1 + R1,4)/2 = {S3.a:.5f}, u3 = (R2,1 - R1,4)/(2δ) = {S3.u:.5f}, "
      f"p3 = {S3.p:.5f}, M3 = {S3.M:.4f}")
print(f"u3-a3 = {S3.u - S3.a:.4f}, u3+a3 = {S3.u + S3.a:.4f}, u4-a4 = {S4.u - S4.a:.4f}, "
      f"u4+a4 = {S4.u + S4.a:.4f}")

# separazione delle famiglie: u+a da (0,0) e u-a da (1,0), pendenze medie fra iniziale e finale
m2 = 0.5 * ((S2.u + S2.a) + (S3.u + S3.a))
m1 = 0.5 * ((S4.u - S4.a) + (S3.u - S3.a))
tS = 1 / (m2 - m1)
xS = m2 * tS
print(f"separazione (pendenze medie {m2:.4f}, {m1:.4f}): t_S = {tS:.4f}, x_S = {xS:.4f}"
      f"   [pendenze iniziali: t = {1 / ((S2.u + S2.a) - (S4.u - S4.a)):.4f}]")
coda_ts = xS + (S3.u + S3.a) * (t_s - tS)
testa_ts = 1 + (S4.u + S4.a) * t_s
print(f"onda u+a a t*: coda ≈ {coda_ts:.4f}, testa = {testa_ts:.4f}")
print(f"onda u+a dopo la separazione: u+a = 3R2(x0) - 2R1,4 = {lam2(S4.R1, S2.R2):.4f} + "
      f"{3 * dR2:.4f} x0  -> crescente: espansione")

# =============================================================================
# 5. PROBLEMA DI RIEMANN NEL CENTRO DI COMPRESSIONE (F §8)
# =============================================================================
riga("5. Riemann in (x*, t*) fra 2 (sinistra) e 3 (destra)")
RP = gd.riemann(S2, S3, x0=x_s, t0=t_s)
print(RP)
print("iterazioni (k, u^k, p2^k, p3^k):")
for it in RP.iterazioni:
    print("   ", tuple(round(float(v), 6) for v in it))
S5 = RP.stati["2"].con(nome="5")
S6 = RP.stati["3"].con(nome="6")
urto = RP.onda_sx
print(f"u(0) primo tentativo isentropico = {RP.iterazioni[0][1]:.5f}")
print(f"urto: W = {urto.w:.4f}, Ms = {urto.Ms:.4f}, p5/p2 = {urto.p2_p1:.4f}, "
      f"rho5/rho2 = {urto.rho2_rho1:.4f}")
print(f"stato 5: p = {S5.p:.4f}, u = {S5.u:.4f}, a = {S5.a:.4f}, rho = {S5.rho:.4f}, M = {S5.M:.4f}, "
      f"s = {S5.s:.4f}")
print(f"stato 6: p = {S6.p:.4f}, u = {S6.u:.4f}, a = {S6.a:.4f}, rho = {S6.rho:.4f}, M = {S6.M:.4f}, "
      f"s = {S6.s:.4f}")
print(f"s5/s2 = {S5.s / S2.s:.4f};  Δs/R = {gd.salto_entropia(S2, S5)['ds_R']:.4f}")
print(f"espansione 6->3: coda u6+a6 = {S6.u + S6.a:.4f}, testa u3+a3 = {S3.u + S3.a:.4f}")
print(f"Lax: λ1(2) = {S2.u - S2.a:.4f} > W = {urto.w:.4f} > λ1(5) = {S5.u - S5.a:.4f}")
print(f"impedenze: rho5 a5 = {S5.rho * S5.a:.4f}, rho6 a6 = {S6.rho * S6.a:.4f}")
print(f"stesso risultato con lo stato 4 a destra: p* = {gd.riemann(S2, S4).p_star:.6f}")
print(f"residuo |1 - p2/p3| = {abs(1 - RP.stati['2'].p / RP.stati['3'].p):.2e}")
# derivazione del Mach relativo dal salto di velocità (F §7, §8)
X = (g + 1) / 4 * (S2.u - S5.u) / S2.a
print(f"X = (γ+1)/4 Δu/a2 = {X:.4f}, Mr = X + sqrt(1+X^2) = {X + np.sqrt(1 + X**2):.4f}")

# =============================================================================
# 6. SOLUZIONE ASINTOTICA: PSEUDO-RIEMANN FRA 1 E 4 (F §9)
# =============================================================================
riga("6. Soluzione asintotica: pseudo-Riemann (1, 4)")
AS = gd.interazione(S1, S4)
print(AS)
Sinf = AS.stati["2"]
print(f"p_inf = {AS.p_star:.4f}, u_inf = {AS.u_star:.4f}, a_inf = {Sinf.a:.4f}, M_inf = {Sinf.M:.4f}")
print(f"coincide con lo stato 3: Δp = {AS.p_star - S3.p:.2e}, Δu = {AS.u_star - S3.u:.2e}")
print(f"vuoto? u4-u1 = {S4.u - S1.u:.3f} < (a1+a4)/δ = {(S1.a + S4.a) / d:.2f}")
v = AS.velocita() if hasattr(AS, "velocita") else None
print("velocita' delle onde asintotiche:", v)

# =============================================================================
# 7. PROFILI A t* + eps
# =============================================================================
riga("7. Profili di p e M a t* + eps")
xt_head = (S1.u - S1.a) * t_s


def profilo_ts(x):
    """p, M a t* + eps (ventaglio esatto F §5; onda u+a con R2 interpolato linearmente: stima)."""
    p = np.empty_like(x)
    M = np.empty_like(x)
    for k, xx in enumerate(x):
        if xx <= xt_head:
            s = S1
        elif xx < x_s:
            s = E.stato_xi(xx / t_s)
        elif xx < coda_ts:
            s = S3
        elif xx < testa_ts:
            r2 = S3.R2 + (S4.R2 - S3.R2) * (xx - coda_ts) / (testa_ts - coda_ts)
            s = gd.isentropica(S1, a=(S4.R1 + r2) / 2).con(u=(r2 - S4.R1) / (2 * d))
        else:
            s = S4
        p[k], M[k] = s.p, s.u / s.a
    return p, M


def profilo_0(x):
    p = np.where(x <= 0, S1.p, np.where(x >= 1, S4.p, (a_lin(np.clip(x, 0, 1)) / S1.a)**esp))
    u = np.where(x <= 0, 0.0, np.where(x >= 1, S4.u, u_lin(np.clip(x, 0, 1))))
    a = np.where(x <= 0, S1.a, np.where(x >= 1, S4.a, a_lin(np.clip(x, 0, 1))))
    return p, u / a


xx = np.concatenate([np.linspace(-1.2, x_s, 800), [x_s + 1e-9], np.linspace(x_s + 1e-6, 2.6, 800)])
pp, MM = profilo_ts(xx)
x0g = np.concatenate([np.linspace(-1.2, 0, 50), [1e-9], np.linspace(1e-6, 2.6, 400)])
p0g, M0g = profilo_0(x0g)
print(f"testa ventaglio x = {xt_head:.4f}; urto x* = {x_s:.4f}; coda u+a ≈ {coda_ts:.4f}; "
      f"testa u+a = {testa_ts:.4f}")
print(f"salto in x*: p {S2.p:.4f} -> {S3.p:.4f}, M {S2.M:.3f} -> {S3.M:.3f}")
print(f"punto sonico a t*: x = 0")

fig, axs = plt.subplots(3, 1, figsize=(7.2, 9.0), sharex=True)
for ax, logy in ((axs[0], False), (axs[1], True)):
    ax.plot(x0g, p0g, color="0.6", ls="--", lw=1.2, label="$t=0$")
    ax.plot(xx, pp, "k", lw=1.6, label=r"$t=t^*+\varepsilon$")
    ax.set_ylabel("$p$" + (" (scala log)" if logy else ""))
    if logy:
        ax.set_yscale("log")
        ax.set_ylim(0.03, 1.5)
    else:
        ax.set_ylim(0, 1.1)
axs[0].legend(loc="lower right")
axs[2].plot(x0g, M0g, color="0.6", ls="--", lw=1.2, label="$t=0$")
axs[2].plot(xx, MM, "k", lw=1.6, label=r"$t=t^*+\varepsilon$")
axs[2].axhline(1, color="0.4", lw=0.8, ls=":")
axs[2].set_ylabel("$M=u/a$")
axs[2].set_xlabel("$x$")
axs[2].set_ylim(-0.1, 2.9)
for ax in axs:
    for xv in (xt_head, x_s, coda_ts, testa_ts):
        ax.axvline(xv, color="0.75", lw=0.7)
lab = [(xt_head - 0.35, "1"), ((xt_head + x_s) / 2, "ventaglio"), ((x_s + coda_ts) / 2, "3"),
       ((coda_ts + testa_ts) / 2, "$u+a$"), (testa_ts + 0.3, "4")]
for xl, tl in lab:
    axs[2].text(xl, 2.55, tl, ha="center", fontsize=9)
axs[2].annotate("urto (2 → 5|6 → 3)", xy=(x_s, 1.5), xytext=(x_s + 0.35, 2.0),
                arrowprops=dict(arrowstyle="->", lw=0.8), fontsize=9)
axs[2].annotate("sonico", xy=(0, 1), xytext=(-0.55, 1.6),
                arrowprops=dict(arrowstyle="->", lw=0.8), fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "profili_tstar.pdf"))
plt.close(fig)

# struttura interna del salto (variabile di similitudine ξ = (x-x*)/(t-t*))
xi = np.linspace(-1.2, 2.3, 1400)
c = RP.campiona(x_s + xi * 1.0, t_s + 1.0)
fig, axs = plt.subplots(1, 2, figsize=(9.0, 3.3))
axs[0].plot(xi, c["p"], "k")
axs[0].set_ylabel("$p$")
axs[1].plot(xi, c["M"], "k")
axs[1].set_ylabel("$M$")
for ax in axs:
    ax.set_xlabel(r"$(x-x^*)/(t-t^*)$")
    for xv, tl in ((urto.w, "urto"), (RP.u_star, "contatto")):
        ax.axvline(xv, color="0.6", lw=0.7, ls="--")
    ax.axvspan(S6.u + S6.a, S3.u + S3.a, color="0.85")
axs[0].text(-1.0, 0.1, "2", fontsize=10)
axs[0].text(-0.2, 0.62, "5", fontsize=10)
axs[0].text(0.9, 0.62, "6", fontsize=10)
axs[0].text(2.0, 0.62, "3", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "struttura_salto.pdf"))
plt.close(fig)

# =============================================================================
# 8. PIANO p-u (F §8): scala lineare e logaritmica
# =============================================================================
riga("8. Piano p-u")
pc = np.geomspace(0.03, 1.3, 600)
uL1 = gd.curva_onda(S1, "u-a", pc)          # isentropica u-a per 1 (R2 = 1.522): contiene 2 e 3
uR4 = gd.curva_onda(S4, "u+a", pc)          # isentropica u+a per 4 (R1 = 1.3487): contiene 3 e 6
uR4_is = S4.u - S4.a / d * (1 - (pc / S4.p)**aria.q)   # ramo isentropico anche per p > p4
ph = np.geomspace(S2.p, 1.3, 400)
uH2 = gd.curva_onda(S2, "u-a", ph)          # Hugoniot u-a da 2
xd = np.linspace(0, 1, 100)
pd_, ud_ = (a_lin(xd) / S1.a)**esp, u_lin(xd)
for logp, nome in ((False, "pu_lin.pdf"), (True, "pu_log.pdf")):
    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    ax.plot(uL1, pc, "b", lw=1.3, label=r"$u-a$ isentropica per 1 ($R_2=1{,}522$)")
    ax.plot(uR4_is, pc, "g", lw=1.3, label=r"$u+a$ isentropica per 4 ($R_1=1{,}349$)")
    ax.plot(uH2, ph, "r", lw=1.3, label=r"urto $u-a$ da 2 (Rankine--Hugoniot)")
    ax.plot(ud_, pd_, color="0.5", ls=":", lw=1.6, label="dati iniziali in $0<x<1$")
    for s, off in ((S1, (6, 4)), (S2, (6, -12)), (S3, (6, 4)), (S4, (6, -12)), (S5, (-30, 6))):
        ax.plot(s.u, s.p, "ko", ms=4)
        ax.annotate("5=6" if s is S5 else s.nome, (s.u, s.p), textcoords="offset points", xytext=off)
    ax.set_xlabel("$u$")
    ax.set_ylabel("$p$")
    if logp:
        ax.set_yscale("log")
        ax.set_ylim(0.03, 1.4)
    else:
        ax.set_ylim(0, 1.25)
    ax.set_xlim(-0.2, 2.9)
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, nome))
    plt.close(fig)
print("figure p-u salvate")

# =============================================================================
# 9. PIANO x-t (qualitativo, con fronti calcolati)
# =============================================================================
riga("9. Piano x-t")
fig, ax = plt.subplots(figsize=(7.6, 6.2))
tm = 1.3
tt = np.linspace(t_s, tm, 300)
xi_inf = S3.u - S3.a
tau = 0.25
# urto dopo t* (andamento qualitativo: parte con W, tende alla coda asintotica u-a = u3-a3)
vel = xi_inf + (urto.w - xi_inf) * np.sqrt(tau / (tau + tt - t_s))
xs_c = x_s + np.concatenate([[0], np.cumsum(0.5 * (vel[1:] + vel[:-1]) * np.diff(tt))])
ax.plot(xs_c, tt, "r", lw=2.2)
# ventaglio centrato in 0: ogni caratteristica si ferma sull'urto (dopo t*)
for xi_ in np.linspace(S1.u - S1.a, S2.u - S2.a, 9):
    k = np.where(xs_c <= xi_ * tt)[0]
    tc = tt[k[0]] if len(k) else tm
    ax.plot([0, xi_ * tc], [0, tc], "b", lw=0.8)
# compressione u-a dai dati in [0,1] (rette, centrata in (x*, t*))
for x0 in np.linspace(0, 1, 6):
    ax.plot([x0, x_s], [0, t_s], "b", lw=0.8, ls="-")
# u-a dallo stato 4 (pendenza u4-a4, poi u3-a3 dopo l'onda u+a) fino all'urto
for x0 in (1.4, 1.9):
    # incontro con la testa u+a  x = 1 + (u4+a4) t
    t_in = (x0 - 1) / ((S4.u + S4.a) - (S4.u - S4.a))
    x_in = x0 + (S4.u - S4.a) * t_in
    t_out = (x_in - xS + (S3.u + S3.a) * tS - (S3.u - S3.a) * t_in) / ((S3.u + S3.a) - (S3.u - S3.a))
    x_out = x_in + 0.5 * ((S4.u - S4.a) + (S3.u - S3.a)) * (t_out - t_in)
    tl = np.linspace(t_out, tm, 200)
    xl = x_out + (S3.u - S3.a) * (tl - t_out)
    j = np.where((tl >= t_s) & (xl <= np.interp(tl, tt, xs_c)))[0]
    tf = tl[j[0]] if len(j) else tm
    ax.plot([x0, x_in, x_out, x_out + (S3.u - S3.a) * (tf - t_out)], [0, t_in, t_out, tf], "b",
            lw=0.8, ls="--")
# zona a due famiglie
ax.fill([0, 1, xS], [0, 0, tS], color="0.85", zorder=0)
# u+a dai dati in [0,1]: rette fino al bordo della zona a due famiglie, poi pendenza finale
for x0 in np.linspace(0, 1, 6):
    s_in = u_lin(x0) + a_lin(x0)
    s_fin = lam2(S4.R1, R2i(x0))
    # bordo destro della zona: retta da (1,0) a S
    sb = (xS - 1) / tS
    tq = (1 - x0) / (s_in - sb) if x0 < 1 else 0.0
    tq = min(tq, tS)
    xq = x0 + s_in * tq
    ax.plot([x0, xq, xq + s_fin * (tm - tq)], [0, tq, tm], "g", lw=0.8)
# testa e coda dell'onda u+a non centrata
ax.plot([1, 1 + (S4.u + S4.a) * tm], [0, tm], "g", lw=1.4)
ax.plot([xS, xS + (S3.u + S3.a) * (tm - tS)], [tS, tm], "g", lw=1.4)
# espansione u+a centrata in (x*, t*) fra 6 e 3
for sl in np.linspace(S6.u + S6.a, S3.u + S3.a, 3):
    ax.plot([x_s, x_s + sl * (tm - t_s)], [t_s, tm], "g", lw=0.9)
# contatto (qualitativo: parte con u5, tende a u3)
velc = S3.u + (RP.u_star - S3.u) * np.sqrt(tau / (tau + tt - t_s))
xc = x_s + np.concatenate([[0], np.cumsum(0.5 * (velc[1:] + velc[:-1]) * np.diff(tt))])
ax.plot(xc, tt, "k--", lw=1.3)
ax.plot([-1.6, 3.4], [t_s, t_s], color="0.4", lw=0.6, ls=":")
ax.text(3.3, t_s + 0.015, r"$t^*$", ha="right", fontsize=9)
ax.plot(x_s, t_s, "ko", ms=4)
txt = [(-1.35, 0.25, "1"), (-0.3, 0.12, "ventaglio\n$u-a$"), (0.5, 0.05, "zona a 2\nfamiglie"),
       (2.2, 0.15, "4"), (0.88, 0.5, "3"), (1.1, 1.0, "6"), (0.25, 1.02, "5"),
       (2.6, 0.75, "onda $u+a$"), (-0.75, 1.2, "urto $u-a$")]
for xl, tl, s in txt:
    ax.text(xl, tl, s, fontsize=9, ha="center")
ax.annotate("contatto", xy=(np.interp(1.15, tt, xc), 1.15), xytext=(1.25, 1.2),
            arrowprops=dict(arrowstyle="->", lw=0.8), fontsize=9)
ax.annotate(r"$C=(x^*,t^*)$", xy=(x_s, t_s), xytext=(-0.75, 0.55),
            arrowprops=dict(arrowstyle="->", lw=0.8), fontsize=9)
ax.annotate("S", xy=(xS, tS), xytext=(xS + 0.25, tS - 0.06),
            arrowprops=dict(arrowstyle="->", lw=0.8), fontsize=9)
ax.set_xlim(-1.6, 3.4)
ax.set_ylim(0, tm)
ax.set_xlabel("$x$")
ax.set_ylabel("$t$")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "xt.pdf"))
plt.close(fig)
print("figura x-t salvata")

# =============================================================================
# 10. VERIFICHE
# =============================================================================
riga("10. Verifiche")
print(f"a+δu continuo in x=0: {abs(S1.R2 - S2.R2) < 1e-9}; u-a cresce da {S1.u - S1.a:.3f} a "
      f"{S2.u - S2.a:.3f}: espansione (ammissibile), transonica")
print(f"a > 0 ovunque: min a = {min(S2.a, S5.a, S6.a, S3.a):.3f}")
print(f"famiglie separate prima di t*: t_S = {tS:.3f} < t* = {t_s}")
print(f"urto: Ms = {urto.Ms:.3f} > 1, s5 > s2: {S5.s > S2.s}; u5 < u2 (urto a sinistra) e "
      f"u6 < u3 (espansione a destra): {S5.u < S2.u and S6.u < S3.u}")
print(f"asintotico: u* > u1 (R) e u* < u4 (R): {AS.u_star > S1.u and AS.u_star < S4.u}")
print(f"urto più lento delle u-a del ventaglio (W < 1.61) -> il ventaglio lo raggiunge e lo indebolisce")

# =============================================================================
# 11. CONTROLLO CON RETICOLO COMPLETO (solo per la sintesi, fuori dal PDF)
# =============================================================================
riga("11. Controllo: reticolo delle caratteristiche nella zona a due famiglie (non nel PDF)")
N = 400
xs = np.linspace(0, 1, N + 1)
Xg = np.full((N + 1, N + 1), np.nan)
Tg = np.full((N + 1, N + 1), np.nan)
for i in range(N + 1):
    Xg[i, i], Tg[i, i] = xs[i], 0.0
for i in range(1, N + 1):
    for j in range(i - 1, -1, -1):
        xa, ta = Xg[i, j + 1], Tg[i, j + 1]
        xb, tb = Xg[i - 1, j], Tg[i - 1, j]
        s1 = 0.5 * (lam1(R1i(xs[i]), R2i(xs[j + 1])) + lam1(R1i(xs[i]), R2i(xs[j])))
        s2 = 0.5 * (lam2(R1i(xs[i - 1]), R2i(xs[j])) + lam2(R1i(xs[i]), R2i(xs[j])))
        t = (xb - xa + s1 * ta - s2 * tb) / (s1 - s2)
        Xg[i, j], Tg[i, j] = xa + s1 * (t - ta), t
print(f"separazione vera: x_S = {Xg[N, 0]:.4f}, t_S = {Tg[N, 0]:.4f} (stima {xS:.4f}, {tS:.4f})")
print(f"coda u+a a t*: {Xg[N, 0] + (S3.u + S3.a) * (t_s - Tg[N, 0]):.4f} (stima {coda_ts:.4f})")
sl = lam1(R1i(xs), S1.R2)
xe, te = Xg[:, 0], Tg[:, 0]
tc = []
for i in range(N):
    a1, b1 = sl[i], xe[i] - sl[i] * te[i]
    a2, b2 = sl[i + 1], xe[i + 1] - sl[i + 1] * te[i + 1]
    tc.append(((b2 - b1) / (a1 - a2), i))
tc.sort()
tb, ib = tc[0]
print(f"primo incrocio di u-a (blow-up vero): t = {tb:.4f}, x = {xe[ib] + sl[ib] * (tb - te[ib]):.4f}")
tcr = (xe[-1] - sl[-1] * te[-1]) / ((S2.u - S2.a) - sl[-1])
print(f"incrocio delle u-a estreme: t = {tcr:.4f}, x = {(S2.u - S2.a) * tcr:.4f}")
