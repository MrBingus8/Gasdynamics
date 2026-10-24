"""Prova del 17/07/2023, esercizio 2 (B72-E02): urto piano in aria in quiete (p = 1, a = 1) con Mach
assoluto a valle M1 = 1,87, che incontra una brusca riduzione d'area A1/A2 = 1,5.

Riproduce tutti i numeri e le figure dello svolgimento.
Uso:  python calcoli.py   (dalla cartella dello svolgimento)

Stati: 0 quiete, 1 dietro l'urto incidente (condotto largo), 2 ingresso del condotto stretto (collegamento
stazionario isentropico 1 -> 2), 3 coda dell'espansione u-a (a sinistra del contatto), 4 dietro l'urto
trasmesso (a destra del contatto). Grandezze adimensionali: p0 = 1, a0 = 1, rho0 = gamma.

Modello del corso: F §7 (urti mobili), F §16 (legge delle aree), F §9 "Variazione d'area" (collegamento
stazionario con m, h0, s conservati, poi nuovo problema di Riemann), F §8 (solutore di Riemann).
Fuori modulo (scritto qui): formula esplicita di Ms dal Mach assoluto a valle (equazione di secondo grado in
Ms^2 ricavata da Rankine-Hugoniot), controllata contro gd.urto_mobile(..., M_valle=...).
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
M1 = 1.87          # Mach assoluto dietro l'urto incidente
A1_A2 = 1.5


def riga(s, nome):
    print(f"stato {nome}: p = {s.p:9.3f}, rho = {s.rho:7.4f}, u = {s.u:8.4f}, a = {s.a:7.4f}, "
          f"M = {s.u/s.a:6.4f}, u-a = {s.u - s.a:7.4f}, u+a = {s.u + s.a:7.4f}, s = p/rho^g = {s.s:.4f}, "
          f"p0 = {s.p0:9.2f}, a0 = {s.a0:.4f}")


# =============================================================================
# 1. URTO INCIDENTE DAL MACH ASSOLUTO A VALLE (F §7)
# =============================================================================
# u1 = w (1 - rho0/rho1)            -> u1/a0 = 2/(g+1) (Ms - 1/Ms)
# a1^2/a0^2 = T1/T0 = [2g Ms^2 - (g-1)] [(g-1) Ms^2 + 2] / [(g+1)^2 Ms^2]
# M1^2 = u1^2/a1^2 = 4 (Ms^2 - 1)^2 / {[2g Ms^2 - (g-1)] [(g-1) Ms^2 + 2]}
# => [4 - 2g(g-1) M1^2] X^2 - [8 + (4g - (g-1)^2) M1^2] X + 4 + 2(g-1) M1^2 = 0,  X = Ms^2
s0 = gd.State.da_pa(aria, 1.0, 1.0, nome="0")
cA = 4 - 2 * g * (g - 1) * M1**2
cB = 8 + (4 * g - (g - 1)**2) * M1**2
cC = 4 + 2 * (g - 1) * M1**2
X = np.array([(cB + sg * np.sqrt(cB**2 - 4 * cA * cC)) / (2 * cA) for sg in (1, -1)])
print("URTO INCIDENTE")
print(f"coefficienti: A = {cA:.5f}, B = {cB:.4f}, C = {cC:.4f}; radici Ms^2 = {X}, Ms = {np.sqrt(X)}")
Ms = np.sqrt(X[0])
print(f"radice fisica (Ms > 1): Ms = {Ms:.4f}; l'altra (Ms = {np.sqrt(X[1]):.3f} < 1) non è un urto")
print(f"limite per urto fortissimo: M1,lim = sqrt(2/(g(g-1))) = {gd.M_valle_limite(aria):.4f} (A -> 0)")
inc = gd.urto_mobile(s0, "u+a", M_valle=M1)
assert abs(inc.Ms - Ms) < 1e-6
s1 = inc.valle.con(nome="1")
r = gd.urto_normale(Ms, aria)
print(f"p1/p0 = {r['p2_p1']:.3f}, rho1/rho0 = {r['rho2_rho1']:.4f}, T1/T0 = {r['T2_T1']:.3f}, a1/a0 = {r['a2_a1']:.4f}")
print(f"w_i = Ms a0 = {inc.w:.4f}; u1 = w (1 - rho0/rho1) = {inc.w*(1 - 1/r['rho2_rho1']):.4f}")
print(f"Mach relativo a valle M_r1 = {inc.Mr_valle:.4f}")
riga(s1, "1")
print("sensibilità al dato M1 (solo sintesi):")
for MM in (1.86, 1.865, 1.87, 1.875, 1.88):
    ii = gd.urto_mobile(s0, "u+a", M_valle=MM)
    rr = gd.riemann_area(ii.valle, s0, A1_A2, 1.0)
    print(f"   M1 = {MM}: Ms = {ii.Ms:7.3f}, regime {rr['regime']}, Ms_t/Ms_i = {rr['onda_trasmessa'].Ms/ii.Ms:.4f}, "
          f"(p4-p0)/(p1-p0) = {(rr['riemann_valle'].p_star - 1)/(ii.valle.p - 1):.4f}")

# =============================================================================
# 2. SALTO D'AREA: COLLEGAMENTO STAZIONARIO ISENTROPICO 1 -> 2 (F §9, §16)
# =============================================================================
print("\nSALTO D'AREA")
AA1 = gd.area_su_area_critica(M1, aria)
AA2 = AA1 / A1_A2
print(f"A1/A* = A(M1) = {AA1:.4f} > A1/A2 = {A1_A2}: la portata passa, A2/A* = {AA2:.4f} > 1")
print(f"soglia: la contrazione non strozza se M1 >= {gd.mach_da_area(A1_A2, aria, 'supersonico'):.4f}")
print(f"margine sulla pressione totale: strozza se p0_2/p0_1 < A1/A2 / A(M1) = {A1_A2/AA1:.4f}")
M2 = gd.mach_da_area(AA2, aria, "supersonico")
f1, f2 = 1 + d * M1**2, 1 + d * M2**2
p2 = s1.p * (f1 / f2)**(g / (g - 1))
a2 = s1.a * np.sqrt(f1 / f2)
s2 = gd.isentropica(s1, p=p2).con(u=M2 * a2, nome="2")
print(f"M2 = {M2:.4f} (ramo supersonico); p2 = p1 [(1+dM1^2)/(1+dM2^2)]^(g/(g-1)) = {p2:.2f}; "
      f"a2 = a1 sqrt(...) = {a2:.4f}")
riga(s2, "2")
sv = gd.variazione_area(s1, A1_A2, 1.0)
assert abs(sv.p - s2.p) / s2.p < 1e-8
print(f"portata rho u A: {s1.rho*s1.u*A1_A2:.4f} (1), {s2.rho*s2.u:.4f} (2); h0: {s1.h0:.4f}, {s2.h0:.4f}; "
      f"s: {s1.s:.4f}, {s2.s:.4f}")
print(f"ramo subsonico (con urto): M = {gd.mach_da_area(AA2, aria, 'subsonico'):.4f} (escluso)")

# =============================================================================
# 3. PROBLEMA DI RIEMANN NEL CONDOTTO STRETTO: 2 (sinistra) | 0 (destra) (F §8)
# =============================================================================
print("\nRIEMANN 2 | 0")
sol = gd.riemann(s2, s0)
z = s0.a / s2.a * (s2.p / s0.p)**((g - 1) / (2 * g))
u0 = (z * (s2.a + d * s2.u) - (s0.a - d * s0.u)) / (d * (1 + z))
print(f"z = {z:.4f}, u(0) = {u0:.4f}  (u(0) < u2 = {s2.u:.3f}: il primo tentativo indica un urto a sinistra; la natura si decide alla convergenza)")
for k, uk, psx, pdx in sol.iterazioni:
    print(f"   k = {k}: u = {uk:9.5f}, p_sx = {psx:9.3f}, p_dx = {pdx:9.3f}, |1 - p_sx/p_dx| = {abs(1 - psx/pdx):.1e}")
print(sol)
s3 = sol.stati["2"].con(nome="3")
s4 = sol.stati["3"].con(nome="4")
riga(s3, "3")
riga(s4, "4")
tras = sol.onda_dx
print(f"stato 3: a3 = a2 - d (u* - u2) = {s2.a - d*(sol.u_star - s2.u):.4f}; p3 = p2 (a3/a2)^(2g/(g-1)) = {s3.p:.2f}")
print(f"urto trasmesso: Mr = {tras.Ms:.4f}, w_t = {tras.w:.4f}, p4/p0 = {tras.p2_p1:.3f}, rho4/rho0 = {tras.rho2_rho1:.4f}")
print(f"espansione u-a: testa u2 - a2 = {s2.u - s2.a:.4f}, coda u3 - a3 = {s3.u - s3.a:.4f}")
ra = gd.riemann_area(s1, s0, A1_A2, 1.0)
assert ra["regime"] == "supersonico" and abs(ra["riemann_valle"].p_star - sol.p_star) < 1e-6
print(f"gd.riemann_area: regime {ra['regime']}, p* = {ra['riemann_valle'].p_star:.3f} (coincide)")

# =============================================================================
# 4. VARIAZIONE D'INTENSITÀ
# =============================================================================
print("\nVARIAZIONE D'INTENSITÀ")
print(f"Ms: {inc.Ms:.3f} -> {tras.Ms:.3f}  (rapporto {tras.Ms/inc.Ms:.4f}, +{100*(tras.Ms/inc.Ms - 1):.1f} %)")
print(f"w : {inc.w:.3f} -> {tras.w:.3f}")
print(f"p a valle / p0: {s1.p:.1f} -> {s4.p:.1f}  (rapporto {s4.p/s1.p:.4f}, +{100*(s4.p/s1.p - 1):.1f} %)")
print(f"(p4 - p0)/(p1 - p0) = {(s4.p - 1)/(s1.p - 1):.4f}")
print(f"u a valle: {s1.u:.3f} -> {s4.u:.3f}")

# =============================================================================
# 5. VERIFICHE
# =============================================================================
print("\nVERIFICHE")
print(f"urto incidente: Ms = {inc.Ms:.2f} > 1, s1/s0 = {s1.s/s0.s:.3f} > 1, M1 = {s1.M:.4f} (dato)")
print(f"stato 1 supersonico: u1 - a1 = {s1.u - s1.a:.3f} > 0 (nessuna u-a risale il condotto largo)")
print(f"stato 2 supersonico: u2 - a2 = {s2.u - s2.a:.4f} > 0 (testa del ventaglio si allontana dal salto)")
print(f"espansione: u* = {sol.u_star:.3f} > u2 = {s2.u:.3f}, p3 = {s3.p:.1f} < p2 = {s2.p:.1f}")
print(f"ordine: testa {s2.u - s2.a:.3f} < coda {s3.u - s3.a:.3f} < contatto {sol.u_star:.3f} < urto {tras.w:.3f}")
print(f"urto trasmesso: Ms = {tras.Ms:.2f} > 1, s4/s0 = {s4.s/s0.s:.3f} > 1")
print(f"contatto: p3 = p4 ({abs(s3.p - s4.p)/s4.p:.1e}), u3 = u4; rho3 = {s3.rho:.3f}, rho4 = {s4.rho:.3f}, "
      f"a3 = {s3.a:.3f}, a4 = {s4.a:.3f}; s3/s0 = {s3.s/s0.s:.3f}, s4/s0 = {s4.s/s0.s:.3f}")
print(f"impedenze rho a: 1: {s1.rho*s1.a:.1f}, 2: {s2.rho*s2.a:.1f}")

# =============================================================================
# 6. FIGURE
# =============================================================================
plt.rcParams.update({"font.size": 10})
COL = dict(urto="C3", contatto="C0", esp="C2", part="0.55", car="C4")

# --- piano x-t (origine nell'arrivo dell'urto sul salto d'area) ---
fig, ax = plt.subplots(figsize=(7.6, 5.6))
tmin, tmax = -1.0, 1.0
ax.axvspan(-22, 0, color="0.94", zorder=0)
ax.axvline(0, color="k", lw=1.4)
ax.text(-21.3, 0.93, "condotto largo $A_1$", fontsize=8.5, va="top")
ax.text(0.6, 0.93, "condotto stretto $A_2$", fontsize=8.5, va="top")
ax.plot([inc.w * tmin, 0], [tmin, 0], color=COL["urto"], lw=2.2)
ax.plot([0, tras.w * tmax], [0, tmax], color=COL["urto"], lw=2.2)
ax.plot([0, sol.u_star * tmax], [0, tmax], color=COL["contatto"], lw=1.6, ls="--")
for lam in np.linspace(s2.u - s2.a, s3.u - s3.a, 6):
    ax.plot([0, lam * tmax], [0, tmax], color=COL["esp"], lw=0.9)
# traiettorie delle particelle (1: u1; 2: u2 fino alla testa; ventaglio: integrazione; 3: u*)
for xp0 in np.arange(-20, 0.1, 2.5):
    tt = [xp0 / inc.w]
    xx = [xp0]
    tt.insert(0, tmin)
    xx.insert(0, xp0)
    t, x = tt[-1], xx[-1]
    dt = 2e-4
    while t < tmax and x < 22:
        if x < 0:
            u = s1.u
        elif t <= 0:
            u = s1.u
        else:
            u = sol.stato_xi(x / t).u if x > 0 else s2.u
        x += u * dt
        t += dt
        xx.append(x)
        tt.append(t)
    ax.plot(xx, tt, color=COL["part"], lw=0.6, ls=":")
# caratteristiche u-a nello stato 1 (vanno verso destra più lente dell'urto: nessuna informazione risale)
for tc0 in (-0.85, -0.55, -0.25):
    xc0 = -21.5
    t_fin = min(tmax, tc0 + (0 - xc0) / (s1.u - s1.a))
    ax.plot([xc0, xc0 + (s1.u - s1.a) * (t_fin - tc0)], [tc0, t_fin], color=COL["car"], lw=0.9, ls="-.")
lab = dict(fontsize=11, ha="center", va="center", bbox=dict(boxstyle="circle,pad=0.15", fc="white", ec="0.5"))
ax.text(-7.0, 0.80, "1", **lab)
ax.text(-19.5, -0.92, "1", **lab)
ax.text(-3.0, -0.85, "0", **lab)
ax.text(18, 0.45, "0", **lab)
ax.text(0.9, 0.62, "2", **lab)
ax.text(10.5, 0.85, "3", **lab)
ax.text(17.2, 0.97, "4", **lab)
ax.annotate(f"urto incidente\n$w_i$ = {inc.w:.1f}", (-9, -9 / inc.w), xytext=(-20.5, -0.25), fontsize=8.5,
            color=COL["urto"], arrowprops=dict(arrowstyle="->", color=COL["urto"], lw=0.8))
ax.annotate(f"urto trasmesso\n$w_t$ = {tras.w:.1f}", (18.5, 18.5 / tras.w), xytext=(9.5, 0.25), fontsize=8.5,
            color=COL["urto"], arrowprops=dict(arrowstyle="->", color=COL["urto"], lw=0.8))
ax.annotate(f"contatto $u^*$ = {sol.u_star:.1f}", (14, 14 / sol.u_star), xytext=(12.5, 0.6), fontsize=8.5,
            color=COL["contatto"], arrowprops=dict(arrowstyle="->", color=COL["contatto"], lw=0.8))
ax.annotate(f"espansione $u-a$\ntesta {s2.u - s2.a:.2f}, coda {s3.u - s3.a:.2f}", (3.2, 0.8),
            xytext=(7.0, -0.55), fontsize=8.5, color=COL["esp"],
            arrowprops=dict(arrowstyle="->", color=COL["esp"], lw=0.8))
ax.annotate(f"caratteristiche $u-a$ nello stato 1\n($u_1-a_1$ = {s1.u - s1.a:.2f} > 0)", (-15.4, 0.33),
            xytext=(-21.3, 0.55), fontsize=8, color=COL["car"],
            arrowprops=dict(arrowstyle="->", color=COL["car"], lw=0.8))
ax.plot([], [], color=COL["part"], lw=0.8, ls=":", label="traiettorie delle particelle")
ax.legend(fontsize=8, loc="lower right")
ax.set_xlim(-22, 22)
ax.set_ylim(tmin, tmax)
ax.set_xlabel("x  (salto d'area in x = 0)")
ax.set_ylabel("t  (t = 0: arrivo dell'urto sul salto d'area)")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "xt.pdf"))
plt.close(fig)

# --- piano p-u, scala lineare e logaritmica ---
fig, axs = plt.subplots(1, 2, figsize=(11, 4.8))
MM = np.linspace(M2, M1, 100)                     # collegamento stazionario al variare dell'area
pl = s1.p * ((1 + d * M1**2) / (1 + d * MM**2))**(g / (g - 1))
ul = MM * s1.a * np.sqrt((1 + d * M1**2) / (1 + d * MM**2))
for ax, log in zip(axs, (False, True)):
    pS = np.geomspace(1.0, 1500, 500)
    ax.plot(gd.curva_onda(s0, "u+a", pS), pS, color=COL["urto"], lw=1.5, label="urti $u+a$ da 0")
    pE = np.geomspace(200, s2.p, 300)
    ax.plot(gd.curva_onda(s2, "u-a", pE), pE, color=COL["esp"], lw=1.5, label="espansione $u-a$ da 2")
    pU = np.geomspace(s2.p, 1500, 100)
    ax.plot(gd.curva_onda(s2, "u-a", pU), pU, color=COL["esp"], lw=1.0, ls=":", label="urti $u-a$ da 2")
    ax.plot(ul, pl, color="0.3", lw=1.3, ls="--", label="salto d'area (stazionario, $x=0$)")
    ax.annotate("", xy=(ul[0], pl[0]), xytext=(ul[50], pl[50]),
                arrowprops=dict(arrowstyle="->", color="0.3", lw=1.1))
    for st, nome, off in ((s0, "0", (6, 4)), (s1, "1", (7, -4)), (s2, "2", (-14, 2)), (s4, "3 = 4", (6, -12))):
        ax.plot(st.u, st.p, "ko", ms=5, zorder=5)
        ax.annotate(nome, (st.u, st.p), xytext=off, textcoords="offset points", fontsize=11)
    ax.set_xlabel("u  (unità di $a_0$)")
    ax.set_ylabel("p  (unità di $p_0$)")
    ax.set_xlim(-1, 20)
    if log:
        ax.set_yscale("log")
        ax.set_ylim(0.7, 2000)
        ax.set_title("(b) scala logaritmica")
    else:
        ax.set_ylim(0, 1400)
        ax.set_title("(a) scala lineare")
        ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3, which="both")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "pu.pdf"))
plt.close(fig)
print("\nfigure salvate in", FIG)
