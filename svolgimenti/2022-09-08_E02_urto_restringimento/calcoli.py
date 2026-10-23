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

Fuori modulo (scritto qui): misure sulle figure della traccia (foto figure/originale_figure.jpg, raddrizzata
con un'omografia sui quattro spigoli di ciascun riquadro), confronto modello/figure e variante del modello
con perdita di pressione totale nel restringimento.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image, ImageOps
from scipy.signal import find_peaks
from scipy.optimize import brentq

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
# 0. MISURE SULLE FIGURE DELLA TRACCIA (fuori modulo)
# =============================================================================
# Foto del foglio con i due grafici (M e fa contro L/A, istanti T = 40, 80, 120, 160, 200).
# Ogni riquadro è raddrizzato con un'omografia che porta i suoi quattro spigoli (trovati come intersezione
# delle rette interpolate sui lati del riquadro) in un rettangolo: x da -6 a 28 (40 px per unità),
# M da 0 a 2,4 (400 px per unità), fa in scala logaritmica da 1 a 7 (960 px per ln 7).
FOTO = ImageOps.exif_transpose(Image.open(os.path.join(FIG, "originale_figure.jpg"))).convert("L")
IM = np.asarray(FOTO).astype(float)
W_R, H_R = 1360, 960


def _lato(x0=None, y0=None, campioni=(), w=12):
    """Retta interpolata sui minimi di luminosità attorno a un lato del riquadro."""
    pts = []
    for c in campioni:
        if x0 is not None:
            seg = IM[c, x0 - w:x0 + w + 1]
            pts.append((x0 - w + int(np.argmin(seg)), c))
        else:
            seg = IM[y0 - w:y0 + w + 1, c]
            pts.append((c, y0 - w + int(np.argmin(seg))))
    p = np.array(pts, float)
    return np.polyfit(p[:, 1], p[:, 0], 1) if x0 is not None else np.polyfit(p[:, 0], p[:, 1], 1)


def spigoli(xl, xr, yt, yb):
    L = _lato(x0=xl, campioni=range(yt + 60, yb - 60, 15))
    R = _lato(x0=xr, campioni=range(yt + 60, yb - 60, 15))
    T = _lato(y0=yt, campioni=range(xl + 250, xr - 30, 15))
    B = _lato(y0=yb, campioni=range(xl + 30, xr - 30, 15))

    def inc(v, h):  # x = v0 y + v1 ; y = h0 x + h1
        y = (h[0] * v[1] + h[1]) / (1 - h[0] * v[0])
        return v[0] * y + v[1], y
    return [inc(L, T), inc(R, T), inc(R, B), inc(L, B)]   # TL, TR, BR, BL


def raddrizza(corners):
    dst = [(0, 0), (W_R, 0), (W_R, H_R), (0, H_R)]
    A = []
    for (x, y), (u, v) in zip(dst, corners):
        A.append([x, y, 1, 0, 0, 0, -u * x, -u * y, -u])
        A.append([0, 0, 0, x, y, 1, -v * x, -v * y, -v])
    Hm = np.linalg.svd(np.array(A))[2][-1].reshape(3, 3)
    c = (Hm / Hm[2, 2]).flatten()[:8]
    return np.asarray(FOTO.transform((W_R, H_R), Image.PERSPECTIVE, tuple(c), Image.BICUBIC)).astype(float)


cM = spigoli(445, 1092, 370, 883)
cF = spigoli(445, 1080, 1100, 1600)
RM, RF = raddrizza(cM), raddrizza(cF)
print("spigoli riquadro M :", [tuple(round(v, 1) for v in c) for c in cM])
print("spigoli riquadro fa:", [tuple(round(v, 1) for v in c) for c in cF])
xpx = lambda i: i / 40 - 6               # noqa: E731
ypM = lambda M: int(H_R - M * 400)       # noqa: E731
ypF = lambda f: int(H_R * (1 - np.log(f) / np.log(7)))  # noqa: E731
fondo = np.median(RM)

# urti: linee verticali nel grafico di M (minimo di luminosità per colonna fra M = 0,1 e 1,2)
colmin = np.array([np.sort(RM[ypM(1.2):ypM(0.1), i])[:40].mean() for i in range(W_R)])
scuro = np.where(colmin < fondo - 40)[0]
gruppi = []
for i in scuro:
    if gruppi and i - gruppi[-1][-1] <= 3:
        gruppi[-1].append(i)
    else:
        gruppi.append([i])
pos = [xpx(np.mean(gr)) for gr in gruppi]
rifl_x = sorted([p for p in pos if -5.6 < p < -0.8], reverse=True)          # T = 40, 80, 120
tras_x = sorted([p for p in pos if 3.0 < p < 25.0])                          # T = 40 ... 200
T_rifl, T_tras = np.array([40, 80, 120]), np.array([40, 80, 120, 160, 200])
assert len(rifl_x) == 3 and len(tras_x) == 5, (rifl_x, tras_x)
c_t, q_t = np.polyfit(T_tras, tras_x, 1)
c_r, q_r = np.polyfit(T_rifl, rifl_x, 1)
T0_t, T0_r = -q_t / c_t, -q_r / c_r
print("\nMISURE SULLE FIGURE")
print("urto riflesso  x(T) =", [round(v, 2) for v in rifl_x], f"-> {c_r:.4f} L/A per unità di T, T0 = {T0_r:.1f}")
print("urto trasmesso x(T) =", [round(v, 2) for v in tras_x], f"-> {c_t:.4f} L/A per unità di T, T0 = {T0_t:.1f}")

# contatto: piccoli massimi di M sopra il plateau (inviluppo superiore delle curve)
inv = []
for i in range(W_R):
    col = RM[ypM(1.75):ypM(1.40), i]
    dd = np.where(col < fondo - 35)[0]
    inv.append(1.75 - dd[0] / 400 if len(dd) else 1.40)
inv = np.array(inv)
picchi, _ = find_peaks(inv, height=1.515, distance=40, prominence=0.02)
cont_x = np.array([xpx(i) for i in picchi])
assert len(cont_x) == 5, cont_x
c_c, q_c = np.polyfit(T_tras[1:], cont_x[1:], 1)      # T = 40 escluso: transitorio
print("contatto (massimi di M) x(T) =", [round(v, 2) for v in cont_x], f"-> {c_c:.4f} L/A per unità di T")
print(f"rapporti: W_r/W_t = {c_r/c_t:.3f}, u*/W_t = {c_c/c_t:.3f}")


def livello(R, x0, x1, ylo, yhi, conv):
    ys = []
    for i in range(int((x0 + 6) * 40), int((x1 + 6) * 40)):
        dd = np.where(R[yhi:ylo, i] < fondo - 35)[0]
        if len(dd):
            ys.append(yhi + np.mean(dd))
    return conv(np.median(ys))


Mval = lambda y: (H_R - y) / 400               # noqa: E731
Fval = lambda y: 7 ** (1 - y / H_R)            # noqa: E731
mis = {
    "M1": livello(RM, -5.5, -2.0, ypM(1.2), ypM(1.5), Mval),
    "M45": livello(RM, 19.0, 22.3, ypM(1.35), ypM(1.53), Mval),
    "fa2_vicino": livello(RF, -1.3, -0.3, ypF(4.5), ypF(6.2), Fval),
    "fa2_lontano": livello(RF, -5.9, -5.0, ypF(4.5), ypF(6.2), Fval),
    "fa3": livello(RF, 0.05, 0.2, ypF(2.0), ypF(3.2), Fval),
    "fa5": np.median([livello(RF, a, b, ypF(1.25), ypF(1.55), Fval) for a, b in ((5, 8.5), (9, 12.5), (14, 17), (19, 22.3))]),
}
for k, v in mis.items():
    print(f"   {k} = {v:.3f}")
print("(incertezze stimate: posizioni ±0,1 L/A, velocità ±2 %, M ±0,02, fa ±0,02 sul plateau, ±0,15 nello stato 2)")

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
# 4b. CONFRONTO CON LE FIGURE E PERDITA NEL RESTRINGIMENTO (fuori modulo)
# =============================================================================


def con_perdita(eta):
    """Stesso modello con p0_3 = eta p0_2 nel restringimento (T0 conservata), sezione A2 sonica.
    Con la perdita la portata critica di A2 cala: A*_2 = eta A2, quindi A1/A*_2 = (A1/A2)/eta."""
    M2e = gd.mach_da_area(A1_A2 / eta, aria, "subsonico")
    u = 0.0
    for _ in range(50):
        r_ = gd.urto_mobile(s1, "u-a", u_valle=u)
        un = M2e * r_.valle.a
        if abs(un - u) < 1e-9:
            break
        u = un
    s2e = r_.valle
    a3e = s2e.a0 * np.sqrt(2 / (g + 1))
    s3e = gd.State.da_pa(aria, eta * s2e.p0 * (2 / (g + 1))**(g / (g - 1)), a3e, u=a3e)
    sole = gd.riemann(s3e, s0)
    return dict(fa2=(s2e.p - s0.p) / (s1.p - s0.p), fa5=(sole.p_star - s0.p) / (s1.p - s0.p),
                Wr_Wt=r_.w / sole.onda_dx.w, ustar_Wt=sole.u_star / sole.onda_dx.w,
                M4=sole.stati["2"].M, M5=sole.stati["3"].M)


eta_fit = brentq(lambda e: con_perdita(e)["fa5"] - mis["fa5"], 0.6, 1.0)
cp = con_perdita(eta_fit)
fa_parete = (parete.valle.p - s0.p) / (s1.p - s0.p)
print("\nCONFRONTO MODELLO / FIGURE")
print(f"M1: modello {s1.u/s1.a:.3f}, figura {mis['M1']:.2f}")
print(f"M4, M5: modello {s4.M:.3f}, {s5.M:.3f}; figura (plateau) {mis['M45']:.2f}")
print(f"fa2: modello {(s2.p - s0.p)/(s1.p - s0.p):.2f}, parete {fa_parete:.2f}; figura {mis['fa2_lontano']:.2f}-{mis['fa2_vicino']:.2f}")
print(f"fa3: modello {(s3.p - s0.p)/(s1.p - s0.p):.2f}; figura {mis['fa3']:.2f}")
print(f"fa5: modello {fa:.3f}; figura {mis['fa5']:.3f}")
print(f"W_r/W_t: modello {W_r/W_t:.3f}, parete {parete.w/W_t:.3f}; figura {c_r/c_t:.3f}")
print(f"u*/W_t: modello {sol.u_star/W_t:.3f}; figura {c_c/c_t:.3f}")
print(f"coda/W_t: modello {(s4.u - s4.a)/W_t:.3f}")
Mt_fig = np.sqrt(((g + 1) * (1 + mis['fa5'] * (s1.p / s0.p - 1)) + g - 1) / (2 * g))
print(f"Mach dell'urto trasmesso dalla figura (fa5 -> p5/p0 -> Mt): {Mt_fig:.2f} (modello {tras.Ms:.2f})")
print(f"perdita che riproduce fa5 della figura: p0_3/p0_2 = {eta_fit:.3f} -> "
      + ", ".join(f"{k} = {v:.3f}" for k, v in cp.items()))

# =============================================================================
# 5. FIGURE
# =============================================================================
plt.rcParams.update({"font.size": 10})
v_fig = W_t / c_t     # m/s per (L/A per unità di T): scala che fa coincidere l'urto trasmesso
T0 = T0_t


def stato_xi(v):
    """Stato del modello alla velocità v = x/t [m/s] (problema autosimile, x = 0 nel salto d'area)."""
    if v < W_r:
        return s1
    if v < 0:
        return s2
    return sol.stato_xi(v)


def profilo(x, T):
    p, M = [], []
    for xx in x:
        st = stato_xi(xx / (T - T0) * v_fig) if T > T0 else (s1 if xx < 0 else s0)
        p.append(st.p)
        M.append(st.u / st.a)
    return (np.array(p) - s0.p) / (s1.p - s0.p), np.array(M)


# --- Figura 1: figure della traccia raddrizzate, annotate, con il modello a T = 200 ---
fig, (axM, axF) = plt.subplots(2, 1, figsize=(7.6, 9.4))
axM.imshow(RM, cmap="gray", extent=(-6, 28, 0, 2.4), aspect="auto", vmin=40, vmax=230)
axF.imshow(RF, cmap="gray", extent=(-6, 28, 0, 1), aspect="auto", vmin=40, vmax=230)
x = np.linspace(-6, 28, 3000)
faT, MT = profilo(x, 200)
axM.plot(x, MT, color="C3", lw=1.3, alpha=0.85, label="modello quasi-1D, T = 200")
faT_plot = np.log(np.clip(faT, 1.0, 7)) / np.log(7)
axF.plot(x[faT >= 1.0], faT_plot[faT >= 1.0], color="C3", lw=1.3, alpha=0.85, label="modello quasi-1D, T = 200")
axF.set_yticks([np.log(v) / np.log(7) for v in (1, 1.2, 1.4, 1.6, 2, 2.5, 3, 4, 5, 6, 7)])
axF.set_yticklabels(["1", "1,2", "1,4", "1,6", "2", "2,5", "3", "4", "5", "6", "7"])
for ax in (axM, axF):
    ax.set_xlim(-6, 28)
    ax.set_xlabel("L/A")
axM.set_ylim(0, 2.4)
axF.set_ylim(0, 1)
axM.set_ylabel("M")
axF.set_ylabel("fa")
box = dict(boxstyle="round,pad=0.2", fc="white", ec="C0", alpha=0.9)
kw = dict(fontsize=8, color="C0", bbox=box, arrowprops=dict(arrowstyle="->", color="C0", lw=0.9))
axM.annotate("urto riflesso", (rifl_x[1], 0.7), xytext=(-5.5, 1.75), **kw)
axM.annotate("sezione sonica\nM = 1", (0.0, 1.0), xytext=(1.5, 0.5), **kw)
axM.annotate("espansione u−a\n(ventaglio centrato)", (1.5, 1.33), xytext=(1.0, 2.0), **kw)
axM.annotate("contatto", (cont_x[3], 1.54), xytext=(11.5, 1.9), **kw)
axM.annotate("urto trasmesso", (tras_x[4], 0.9), xytext=(19.5, 0.45), **kw)
lab = dict(fontsize=10, color="C3", ha="center", fontweight="bold")
axM.text(-5.6, 1.45, "1", **lab)
axM.text(-1.0, 0.15, "2", **lab)
axM.text(0.4, 0.75, "3", **lab)
axM.text(8.0, 1.62, "4", **lab)
axM.text(20.5, 1.62, "5", **lab)
axM.text(25.0, 0.1, "0", **lab)
yl = lambda f: np.log(f) / np.log(7)   # noqa: E731
axF.annotate("urto riflesso", (rifl_x[2], yl(3.0)), xytext=(-5.5, yl(1.9)), **kw)
axF.annotate("sezione sonica", (0.1, yl(2.6)), xytext=(2.0, yl(3.3)), **kw)
axF.annotate("espansione u−a", (1.3, yl(1.75)), xytext=(4.5, yl(2.2)), **kw)
axF.annotate("contatto: non\nvisibile in p", (cont_x[3], yl(1.38)), xytext=(10.5, yl(1.75)), **kw)
axF.annotate("urto trasmesso", (tras_x[4], yl(1.15)), xytext=(22.5, yl(1.6)), **kw)
axF.text(-3.6, yl(5.65), "2", **lab)
axF.text(0.5, yl(2.9), "3", **lab)
axF.text(8.0, yl(1.45), "4", **lab)
axF.text(20.5, yl(1.45), "5", **lab)
for xx, TT in list(zip(tras_x, T_tras)) + list(zip(rifl_x, T_rifl)):
    axF.text(xx, yl(1.07), f"T={TT}", fontsize=6.5, color="C2", ha="center",
             bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8))
axM.legend(fontsize=8, loc="upper center")
axF.legend(fontsize=8, loc="upper center")
axM.set_title("(a) numero di Mach", fontsize=10)
axF.set_title("(b) fattore di amplificazione (scala logaritmica)", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "figure_annotate.pdf"), dpi=200)
plt.close(fig)

# --- Figura 2: piano x-t nelle unità delle figure, con le posizioni misurate ---
fig, ax = plt.subplots(figsize=(7.4, 5.2))
Tm = 230
ax.axvspan(-12, 0, color="0.93")
ax.axvline(0, color="k", lw=1.2)
TT = np.array([T0, Tm])
ax.plot(c_r * (TT - T0_r), TT, "r", lw=2)
ax.plot(c_t * (TT - T0_t), TT, "r", lw=2)
ax.plot(c_c * (TT - T0) + 0 * TT, TT, "b--", lw=1.4)
ax.plot([-12, 0], [T0 - 12 / (W_i / v_fig), T0], "r", lw=2)   # incidente (pendenza del modello)
for k in range(1, 6):
    ax.plot([0, (s4.u - s4.a) / v_fig * k / 5 * (Tm - T0)], [T0, Tm], "g", lw=0.8)
ax.plot(rifl_x, T_rifl, "ko", ms=4, label="posizioni misurate sulle figure")
ax.plot(tras_x, T_tras, "ko", ms=4)
ax.plot(cont_x[1:], T_tras[1:], "o", mfc="white", mec="b", ms=4, label="massimi di M (contatto)")
for vel, nome in ((W_r, "riflesso"), (W_t, "trasmesso"), (sol.u_star, "contatto")):
    ax.plot(vel / v_fig * (TT - T0), TT, color="0.4", lw=0.8, ls=":")
ax.plot([], [], color="0.4", lw=0.8, ls=":", label="modello (stessa scala dell'urto trasmesso)")
lab = dict(fontsize=11, ha="center", va="center", bbox=dict(boxstyle="circle,pad=0.15", fc="white", ec="0.5"))
ax.text(-9.5, 130, "1", **lab)
ax.text(-3.0, 170, "2", **lab)
ax.text(5.0, 200, "4", **lab)
ax.text(19.5, 200, "5", **lab)
ax.text(18.0, 60, "0", **lab)
ax.text(-4, 20, "1", **lab)
ax.text(5, -10, "0", **lab)
ax.text(c_r * (150 - T0_r) - 0.4, 150, "urto riflesso", color="r", fontsize=8.5, ha="right", rotation=-80)
ax.text(c_t * (120 - T0_t) + 0.6, 120, "urto trasmesso", color="r", fontsize=8.5, rotation=50)
ax.text(c_c * (150 - T0) - 0.6, 150, "contatto", color="b", fontsize=8.5, ha="right", rotation=55)
ax.text(2.2, 225, "espansione u−a", color="g", fontsize=8.5, ha="center", va="top")
ax.annotate("3 (sezione sonica)", (0, 60), xytext=(2.5, 30), fontsize=8.5, arrowprops=dict(arrowstyle="->", lw=0.8))
ax.set_xlim(-12, 28)
ax.set_ylim(-30, Tm)
ax.set_xlabel("x  [L/A]   (salto d'area in x = 0)")
ax.set_ylabel("T")
ax.legend(fontsize=8, loc="lower right")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "xt.pdf"))
plt.close(fig)

# --- Figura 3: piano p-u (lineare e logaritmico) ---
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
for ax, log in zip(axs, (False, True)):
    pp = np.geomspace(0.5 * ATM, 60 * ATM, 400)
    pS = pp[pp >= s0.p]
    ax.plot(gd.curva_onda(s0, "u+a", pS), pS / ATM, "r", lw=1.4, label="urti $u+a$ da 0")
    pR = pp[pp >= s1.p]
    ax.plot(gd.curva_onda(s1, "u-a", pR), pR / ATM, "m", lw=1.4, label="urto $u-a$ da 1 (riflesso)")
    pE = np.geomspace(5 * ATM, s3.p, 200)
    ax.plot(gd.curva_onda(s3, "u-a", pE), pE / ATM, "g", lw=1.4, label="espansione $u-a$ da 3")
    ax.annotate("", xy=(s3.u, s3.p / ATM), xytext=(s2.u, s2.p / ATM),
                arrowprops=dict(arrowstyle="->", ls="--", color="0.3", lw=1.1))
    ax.text(0.5 * (s2.u + s3.u) + 20, 0.5 * (s2.p + s3.p) / ATM * (1.05 if log else 1.0),
            "salto d'area\n(stazionario, x = 0)", fontsize=8, color="0.3")
    for st, nome, off in ((s0, "0", (6, -12)), (s1, "1", (6, 2)), (s2, "2", (6, 2)), (s3, "3", (6, 2)),
                          (s5, "4 = 5", (6, 4))):
        ax.plot(st.u, st.p / ATM, "ko", ms=5)
        ax.annotate(nome, (st.u, st.p / ATM), xytext=off, textcoords="offset points", fontsize=11)
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

# --- ritagli della foto per il PDF ---
try:
    Image.open(os.path.join(FIG, "traccia.png")).crop((110, 295, 560, 435)).save(os.path.join(FIG, "schema.png"))
except Exception as e:  # pragma: no cover
    print("ritaglio non riuscito:", e)
print("\nfigure salvate in", FIG)
