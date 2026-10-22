"""Prova di settembre 2022, esercizio 2: urto piano M = 3 in aria in quiete che incontra una brusca
riduzione d'area A1/A2 = 10 in un condotto assialsimmetrico (riflessione e trasmissione).

Riproduce tutti i numeri e le figure dello svolgimento.
Uso:  python calcoli.py   (dalla cartella dello svolgimento)

Stati (etichette della traccia): 0 quiete, 1 dietro l'urto incidente, 2 dietro l'urto riflesso (condotto
largo), 3 sezione d'ingresso del condotto stretto (sonica), 4 coda dell'espansione u-a (a sinistra del
contatto), 5 dietro l'urto trasmesso (a destra del contatto). fa = (p5 - p0)/(p1 - p0).

Modello del corso: F §9 "Variazione d'area" (collegamento stazionario isentropico ṁ, h0, s) e S p.48-49;
il regime strozzato è descritto dal corso solo a parole (S p.48 punto 4) ed è risolto da gd.riemann_area
(fuori formulario, vedi gasdyn_guida §4). Qui lo stesso risultato è ricostruito anche passo per passo,
come si farebbe con la calcolatrice, e confrontato con gd.riemann_area.
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

aria = gd.Gas(1.4, 287.0)
g, d = aria.gamma, aria.delta
ATM = gd.ATM
A1_A2 = 10.0
Ms = 3.0


def riga(s, nome):
    print(f"stato {nome}: p = {s.p/ATM:8.3f} atm, rho = {s.rho:7.3f} kg/m3, u = {s.u:8.2f} m/s, "
          f"a = {s.a:7.2f} m/s, T = {s.T:7.1f} K, M = {s.u/s.a:6.4f}, p0 = {s.p0/ATM:7.2f} atm, "
          f"s = p/rho^g = {s.s:.4e}")


# =============================================================================
# 1. STATO 0 E URTO INCIDENTE (F §7)
# =============================================================================
s0 = gd.State(aria, 1.29, 0.0, ATM, nome="0")
print(f"a0 = sqrt(g p0/rho0) = {s0.a:.2f} m/s, T0 = {s0.T:.1f} K")
r = gd.urto_normale(Ms, aria)
print(f"Ms = 3: p1/p0 = {r['p2_p1']:.4f}, rho1/rho0 = {r['rho2_rho1']:.4f}, T1/T0 = {r['T2_T1']:.4f}")
inc = gd.urto_mobile(s0, "u+a", Ms=Ms)
s1 = inc.valle.con(nome="1")
W_i = inc.w
print(f"W_i = Ms a0 = {W_i:.2f} m/s; u1 = W_i (1 - rho0/rho1) = {W_i*(1 - 1/r['rho2_rho1']):.2f} m/s")
riga(s1, "1")
AA1 = gd.area_su_area_critica(s1.u / s1.a, aria)
print(f"A/A*(M1) = {AA1:.4f}  <<  A1/A2 = {A1_A2}: la corrente supersonica 1 non passa la contrazione")

# =============================================================================
# 2. URTO RIFLESSO E STRANGOLAMENTO (F §7, §16)
# =============================================================================
M2 = gd.mach_da_area(A1_A2, aria, "subsonico")
print(f"\nsezione A2 sonica => A1/A* = {A1_A2} => M2 = {M2:.5f} (ramo subsonico)")


def riflesso(u2):
    """Urto u-a nello stato 1 che porta la velocità a u2 (formula di M_rel del formulario §8)."""
    X = (g + 1) / 4 * (s1.u - u2) / s1.a
    Mr = X + np.sqrt(1 + X**2)
    p2_p1 = 1 + 2 * g / (g + 1) * (Mr**2 - 1)
    a2 = s1.a * np.sqrt((g + 1 + (g - 1) * p2_p1) / (g + 1 + (g - 1) / p2_p1))
    return Mr, p2_p1, a2


u2 = 0.0   # primo tentativo: riflessione come su parete (M2 piccolo)
print("iterazione su u2 (primo tentativo u2 = 0, riflessione su parete):")
for k in range(6):
    Mr, p2_p1, a2 = riflesso(u2)
    u2_nuovo = M2 * a2
    print(f"   k = {k}: u2 = {u2:7.2f} m/s -> M_r = {Mr:.4f}, p2/p1 = {p2_p1:.4f}, a2 = {a2:.2f} m/s, "
          f"M2 a2 = {u2_nuovo:.2f} m/s")
    if abs(u2_nuovo - u2) < 1e-6:
        break
    u2 = u2_nuovo
rifl = gd.urto_mobile(s1, "u-a", u_valle=u2)
s2 = rifl.valle.con(nome="2")
W_r = rifl.w
print(f"urto riflesso: M_r = {rifl.Ms:.4f}, W_r = u1 - M_r a1 = {s1.u - rifl.Ms*s1.a:.2f} m/s "
      f"(modulo: {W_r:.2f})")
riga(s2, "2")
parete = gd.urto_riflesso_parete(s1, "u-a")
print(f"confronto parete chiusa: M_r = {parete.Ms:.4f}, W = {parete.w:.2f} m/s, p = {parete.valle.p/ATM:.2f} atm")

# stato 3: isentropico stazionario da 2 (stessi p0, a0, s), M3 = 1
p02, a02 = s2.p0, s2.a0
p3 = p02 * (2 / (g + 1))**(g / (g - 1))
a3 = a02 * np.sqrt(2 / (g + 1))
s3 = gd.isentropica(s2, p=p3).con(u=a3, nome="3")
print(f"p02 = p2 (1 + d M2^2)^(g/(g-1)) = {p02/ATM:.3f} atm, a02 = {a02:.2f} m/s")
print(f"p3 = 0,5283 p02 = {p3/ATM:.3f} atm, a3 = u3 = 0,9129 a02 = {a3:.2f} m/s")
riga(s3, "3")
print(f"portata: rho2 u2 A1 / (rho3 u3 A2) = {s2.rho*s2.u*A1_A2/(s3.rho*s3.u):.6f}; "
      f"h0: {s2.h0:.1f} vs {s3.h0:.1f}")

# =============================================================================
# 3. PROBLEMA DI RIEMANN NEL CONDOTTO STRETTO: 3 (sinistra) | 0 (destra) (F §8)
# =============================================================================
sol = gd.riemann(s3, s0)
print("\n", sol)
z = s0.a / s3.a * (s3.p / s0.p)**((g - 1) / (2 * g))
print(f"z = {z:.4f}")
for k, uk, p2k, p3k in sol.iterazioni:
    print(f"   k = {k}: u = {uk:8.2f} m/s, p_sx = {p2k/ATM:7.3f} atm, p_dx = {p3k/ATM:7.3f} atm, "
          f"|dp|/p = {abs(p2k - p3k)/p3k:.2e}")
s4 = sol.stati["2"].con(nome="4")
s5 = sol.stati["3"].con(nome="5")
riga(s4, "4")
riga(s5, "5")
tras = sol.onda_dx
W_t = tras.w
print(f"espansione u-a: testa u3 - a3 = {s3.u - s3.a:.3e} m/s, coda u4 - a4 = {s4.u - s4.a:.2f} m/s")
print(f"contatto: u* = {sol.u_star:.2f} m/s; urto trasmesso: W_t = {W_t:.2f} m/s, Ms = {tras.Ms:.4f}")
fa = (s5.p - s0.p) / (s1.p - s0.p)
print(f"\nfa = (p5 - p0)/(p1 - p0) = ({s5.p/ATM:.3f} - 1)/({s1.p/ATM:.3f} - 1) = {fa:.4f}")
print(f"(p - p0)/(p1 - p0): stato 2 {(s2.p - s0.p)/(s1.p - s0.p):.3f}, stato 3 {(s3.p - s0.p)/(s1.p - s0.p):.3f}")
print(f"velocità di 5 rispetto a 1: u5/u1 = {s5.u/s1.u:.3f}; p5/p0 = {s5.p/s0.p:.3f}")

# controllo con la funzione del modulo
ra = gd.riemann_area(s1, s0, A1_A2, 1.0)
assert ra["regime"] == "strozzato"
assert abs(ra["stato_2p"].p - s2.p) / s2.p < 1e-6
assert abs(ra["riemann_valle"].p_star - sol.p_star) / sol.p_star < 1e-6
print(f"gd.riemann_area: regime {ra['regime']}, p* = {ra['riemann_valle'].p_star/ATM:.3f} atm (coincide)")

# =============================================================================
# 4. VERIFICHE
# =============================================================================
print("\nVERIFICHE")
print(f"urto incidente: Ms = {inc.Ms:.3f} > 1, s1/s0 = {s1.s/s0.s:.3f} > 1")
print(f"urto riflesso: M_r monte = {rifl.Mr_monte:.3f} > 1, s2/s1 = {s2.s/s1.s:.3f} > 1, W_r = {W_r:.1f} < 0 (risale nel condotto largo)")
print(f"stato 2 subsonico: M2 = {s2.u/s2.a:.4f}; u2 - a2 = {s2.u - s2.a:.1f} < 0 (le u-a dal salto d'area risalgono)")
print(f"espansione: p4 = {s4.p/ATM:.2f} < p3 = {s3.p/ATM:.2f} atm, u4 = {s4.u:.1f} > u3 = {s3.u:.1f}; "
      f"coda {s4.u - s4.a:.1f} > 0 = testa (tutta nel condotto stretto)")
print(f"ordine: coda {s4.u - s4.a:.1f} < contatto {sol.u_star:.1f} < urto {W_t:.1f}")
print(f"urto trasmesso: Ms = {tras.Ms:.3f} > 1, s5/s0 = {s5.s/s0.s:.3f} > 1")
print(f"contatto: p4 = p5 ({abs(s4.p - s5.p)/s5.p:.1e}), u4 = u5; rho4 = {s4.rho:.3f}, rho5 = {s5.rho:.3f}, "
      f"T4 = {s4.T:.0f} K, T5 = {s5.T:.0f} K")
print(f"entropia: s4/s0 = {s4.s/s0.s:.3f} (due urti), s5/s0 = {s5.s/s0.s:.3f} (un urto)")

# confronto con altri rapporti d'area (discussione)
print("\nfa al variare di A1/A2 (stesso modello):")
for AR in (1.5, 2.0, 5.0, 10.0, 100.0, 1e4):
    rr = gd.riemann_area(s1, s0, AR, 1.0)
    ps = rr["riemann_valle"].p_star
    print(f"   A1/A2 = {AR:>7g}: regime {rr['regime']}, p5 = {ps/ATM:6.2f} atm, "
          f"fa = {(ps - s0.p)/(s1.p - s0.p):.3f}, Ms trasmesso = {rr['onda_trasmessa'].Ms:.3f}")
print(f"soglia di inghiottimento: A1/A2 < A/A*(M1) = {AA1:.4f}")

# =============================================================================
# 5. FIGURE
# =============================================================================
plt.rcParams.update({"font.size": 10})
ms = 1e3  # t in ms

# --- campo a un istante: funzione (x, t) -> stato ---


def campo(x, t):
    """Restituisce (p, M) nel punto x [m] all'istante t [s] (t = 0 arrivo dell'urto in x = 0)."""
    if t <= 0:
        return (s1.p, s1.u / s1.a) if x < W_i * t else (s0.p, 0.0)
    if x < W_r * t:
        return s1.p, s1.u / s1.a
    if x < 0:
        return s2.p, s2.u / s2.a
    c = sol.campiona(np.array([x]), t)
    return float(c["p"][0]), float(c["M"][0])


# --- Figura 1: piano x-t ---
fig, ax = plt.subplots(figsize=(7.4, 5.0))
t_min, t_max = -0.5e-3, 1.0e-3
x_min, x_max = -0.6, 1.3
ax.axvspan(x_min, 0, color="0.93")
ax.text(-0.58, -0.45, "condotto largo $A_1$", fontsize=9)
ax.text(1.27, -0.45, "condotto stretto $A_2=A_1/10$", fontsize=9, ha="right")
ax.axvline(0, color="k", lw=1.2)
ax.plot([W_i * t_min, 0], [t_min * ms, 0], "r", lw=2.2)
ax.plot([0, W_r * t_max], [0, t_max * ms], "r", lw=2.2)
ax.plot([0, W_t * t_max], [0, t_max * ms], "r", lw=2.2)
ax.plot([0, sol.u_star * t_max], [0, t_max * ms], "b--", lw=1.4)
for xi in np.linspace(0, s4.u - s4.a, 6)[1:]:
    ax.plot([0, xi * t_max], [0, t_max * ms], "g", lw=0.9)
# traiettoria di una particella del condotto largo (u1 fino all'urto riflesso, poi u2 << u1)
xp, tp = -0.45, t_min
tA = (xp - s1.u * tp) / (W_r - s1.u)  # incontro con l'urto riflesso
xA = W_r * tA
ax.plot([xp, xA], [tp * ms, tA * ms], "k:", lw=1)
ax.plot([xA, xA + s2.u * (t_max - tA)], [tA * ms, t_max * ms], "k:", lw=1)
ax.text(xA - 0.06, tA * ms - 0.12, "traiettoria", fontsize=8, ha="right")
lab = dict(fontsize=12, ha="center", va="center",
           bbox=dict(boxstyle="circle,pad=0.15", fc="white", ec="0.5"))
ax.text(-0.35, -0.05, "1", **lab)
ax.text(0.15, -0.35, "0", **lab)
ax.text(-0.15, 0.75, "2", **lab)
ax.text(0.45, 0.8, "4", **lab)
ax.text(0.98, 0.95, "5", **lab)
ax.text(1.15, 0.5, "0", **lab)
ax.text(-0.5, 0.75, "1", **lab)
ax.annotate("3 (sezione sonica)", (0.0, 0.3), xytext=(0.18, 0.12), fontsize=8.5,
            arrowprops=dict(arrowstyle="->", lw=0.8))
ax.text(W_i * t_min * 0.55, t_min * ms * 0.55 + 0.04, "urto incidente", fontsize=8.5, color="r",
        rotation=60, ha="right")
ax.text(W_r * 0.6e-3 - 0.01, 0.6, "urto riflesso", fontsize=8.5, color="r", rotation=-70, ha="right")
ax.text(W_t * 0.75e-3 + 0.04, 0.75, "urto trasmesso", fontsize=8.5, color="r", rotation=40)
ax.text(sol.u_star * 0.5e-3 - 0.02, 0.56, "contatto", fontsize=8.5, color="b", rotation=46, ha="right")
ax.text((s4.u - s4.a) * 0.9e-3 - 0.06, 0.93, "espansione $u-a$", fontsize=8.5, color="g", ha="right")
ax.set_xlim(x_min, x_max)
ax.set_ylim(t_min * ms, t_max * ms)
ax.set_xlabel("x [m]  (salto d'area in x = 0)")
ax.set_ylabel("t [ms]  (t = 0: arrivo dell'urto)")
ax.set_title("Piano x–t (pendenze in scala: W_i = 995, W_r = −365, coda 302, u* = 900, W_t = 1174 m/s)",
             fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "xt.pdf"))
plt.close(fig)

# --- Figura 2: piano p-u (lineare e logaritmico) ---


def curva(stato, famiglia, p):
    return gd.curva_onda(stato, famiglia, p)


fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
for ax, log in zip(axs, (False, True)):
    pp = np.geomspace(0.5 * ATM, 60 * ATM, 400)
    # curva u+a uscente da 0 (urti a destra): urto incidente e urto trasmesso
    pS = pp[pp >= s0.p]
    ax.plot(curva(s0, "u+a", pS), pS / ATM, "r", lw=1.4, label="urti $u+a$ da 0 (NCS)")
    # curva u-a da 1 (urto riflesso)
    pR = pp[pp >= s1.p]
    ax.plot(curva(s1, "u-a", pR), pR / ATM, "m", lw=1.4, label="urto $u-a$ da 1 (riflesso)")
    # curva u-a da 3 (espansione nel condotto stretto)
    pE = np.geomspace(5 * ATM, s3.p, 200)
    ax.plot(curva(s3, "u-a", pE), pE / ATM, "g", lw=1.4, label="espansione $u-a$ da 3")
    ax.annotate("", xy=(s3.u, s3.p / ATM), xytext=(s2.u, s2.p / ATM),
                arrowprops=dict(arrowstyle="->", ls="--", color="0.3", lw=1.1))
    ax.text(0.5 * (s2.u + s3.u) + 20, 0.5 * (s2.p + s3.p) / ATM * (1.05 if log else 1.0),
            "salto d'area\n(stazionario, x = 0)", fontsize=8, color="0.3")
    for s, nome, off in ((s0, "0", (6, -12)), (s1, "1", (6, 2)), (s2, "2", (6, 2)), (s3, "3", (6, 2)),
                         (s5, "4 = 5", (6, 4))):
        ax.plot(s.u, s.p / ATM, "ko", ms=5)
        ax.annotate(nome, (s.u, s.p / ATM), xytext=off, textcoords="offset points", fontsize=11)
    ax.set_xlabel("u [m/s]")
    ax.set_ylabel("p [atm]")
    ax.set_xlim(-50, 1200)
    if log:
        ax.set_yscale("log")
        ax.set_ylim(0.8, 70)
        ax.set_title("(b) scala logaritmica")
    else:
        ax.set_ylim(0, 55)
        ax.set_title("(a) scala lineare")
        ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3, which="both")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "pu.pdf"))
plt.close(fig)

# --- Figura 3: profili di M e fa lungo il condotto a diversi istanti ---
fig, (axM, axf) = plt.subplots(2, 1, figsize=(7.4, 6.4), sharex=True)
x = np.linspace(-0.6, 1.3, 2400)
for t, col in ((-0.3e-3, "0.5"), (0.3e-3, "C0"), (0.8e-3, "C3")):
    pM = np.array([campo(xx, t) for xx in x])
    axM.plot(x, pM[:, 1], color=col, lw=1.5, label=f"t = {t*ms:+.1f} ms")
    axf.plot(x, (pM[:, 0] - s0.p) / (s1.p - s0.p), color=col, lw=1.5, label=f"t = {t*ms:+.1f} ms")
for ax in (axM, axf):
    ax.axvline(0, color="k", lw=1)
    ax.grid(alpha=0.3)
axM.axhline(1, color="0.6", lw=0.8, ls=":")
axM.set_ylabel("M = u/a")
axM.legend(fontsize=8.5, loc="center right")
axM.text(-0.55, 1.42, "1", fontsize=11)
axM.text(-0.12, 0.12, "2", fontsize=11)
axM.text(0.45, 1.58, "4", fontsize=11)
axM.text(0.85, 1.58, "5", fontsize=11)
axM.set_ylim(-0.1, 1.85)
axf.axhline(1, color="0.6", lw=0.8, ls=":")
axf.axhline(fa, color="C3", lw=0.8, ls=":")
axf.text(1.02, fa + 0.12, f"fa = {fa:.2f}", color="C3", fontsize=9)
axf.text(-0.55, 1.15, "1", fontsize=11)
axf.text(-0.12, 5.25, "2", fontsize=11)
axf.set_ylim(-0.3, 5.9)
axf.text(0.02, 2.85, "3", fontsize=11)
axf.set_ylabel(r"$(p-p_0)/(p_1-p_0)$")
axf.set_xlabel("x [m]  (salto d'area in x = 0)")
axM.set_title("Profili lungo il condotto prima e dopo l'interazione", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "profili.pdf"))
plt.close(fig)

# --- Figura 4: schema della traccia (ritaglio della foto) ---
try:
    from PIL import Image
    im = Image.open(os.path.join(FIG, "traccia.png"))
    im.crop((110, 295, 560, 435)).save(os.path.join(FIG, "schema.png"))
except Exception as e:  # pragma: no cover
    print("ritaglio non riuscito:", e)
print("\nfigure salvate in", FIG)
