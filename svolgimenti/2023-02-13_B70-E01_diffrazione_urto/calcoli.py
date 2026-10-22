"""Prova del 13/02/2023, esercizio 1 (B70-E01): diffrazione di un urto piano M_s = 4 attorno a uno spigolo di 105°.

Riproduce tutte le misure sull'immagine, tutti i numeri e le figure dello svolgimento.
Uso:  python calcoli.py   (dalla cartella dello svolgimento)

Immagine: la foto della traccia (figure/originale.jpeg) è una foto del foglio stampato; per le misure si usa
la stessa schlieren estratta, a risoluzione nativa, dal PDF della traccia in fonti/soluzioni_A/A_G_B70-E01.pdf
(figure/schlieren_traccia.jpg, 638 x 481 px, senza prospettiva). Coordinate immagine: x verso destra,
y verso il basso, in pixel. Angoli: gradi, antiorari dall'orizzontale (verso della corrente).
"""
import os
import sys
import math

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.signal import find_peaks
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

QUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(QUI, "..", ".."))
import gasdyn as gd  # noqa: E402

FIG = os.path.join(QUI, "figure")
G = gd.Gas(1.4, 287.0)
deg, rad = math.degrees, math.radians


def nu(M):
    return gd.prandtl_meyer(M, G)


def M_da_nu(v):
    return gd.mach_da_pm(v, G)


def p_p0(M):
    return 1.0 / gd.isentropiche(M, G)["p0_p"]


def a_a0(M):
    return 1.0 / gd.isentropiche(M, G)["a0_a"]


def mu(M):
    return gd.angolo_mach(M)


# =============================================================================
# 1. MISURE SULL'IMMAGINE
# =============================================================================
im = np.asarray(Image.open(os.path.join(FIG, "schlieren_traccia.jpg")).convert("L")).astype(float)
g = gaussian_filter(im, 1.2)
H, W_img = im.shape

# parete superiore: banda nera in alto
y_parete = int(np.argmax(im[:, 400] > 100))
# spigolo C: incrocio della verticale e dell'orizzontale tracciate nella traccia
xC = 150 + int(np.argmin(im[20:200, 150:190].mean(0)))
yC = 230 + int(np.argmin(im[230:255, 180:240].mean(1)))
C = np.array([xC, yC], float)

# parete inclinata: cresta luminosa sottile, fit x(y) (F: geometria, nessuna formula)
ys, xs = [], []
for y in range(yC + 12, 475, 5):
    xe = xC - (y - yC) * 0.275
    lo = int(xe - 10)
    if y == 435:            # in questa riga c'è un graffio della stampa
        continue
    xs.append(lo + np.argmax(im[y, lo:lo + 20]))
    ys.append(y)
pw = np.polyfit(ys, xs, 1)
ang_parete_vert = deg(math.atan(-pw[0]))           # inclinazione sulla verticale
deflessione_parete = 90.0 + ang_parete_vert         # rotazione della parete rispetto alla corrente
dir_parete = np.array([pw[0], 1.0]) / math.hypot(pw[0], 1.0)   # versore C -> B (verso il basso)

# piede B: massimo di luminosità lungo la parete (incrocio di DS con la parete)
s = np.arange(100, 220, 0.5)
P = C[None, :] + s[:, None] * dir_parete[None, :]
prof = map_coordinates(g, [P[:, 1], P[:, 0]], order=1)
sB_px = s[np.argmax(prof)]
B = C + sB_px * dir_parete

# testa A: cresta di DS subito sotto la parete superiore
righe = range(y_parete + 3, y_parete + 9)
xA = np.mean([480 + np.argmax(g[y, 480:625]) for y in righe])
xA_px = xA - xC
rapporto_misurato = xA_px / sB_px

# tratto rettilineo AD: cresta di DS riga per riga, scostamento dalla verticale per A
y_ds = np.arange(y_parete + 4, 200, 2)
x_ds = np.array([480 + np.argmax(g[y, 480:625]) for y in y_ds], float)
x_ds_l = np.convolve(x_ds, np.ones(5) / 5, mode="same")
scost = xA - x_ds_l
i_D = np.argmax((scost > 3.0) & (y_ds > 30))
yD_oss = y_ds[i_D]
hD_oss = (yC - yD_oss) / xA_px                       # quota di D sopra C, in unità x_A

# linea di separazione SL (tratto nero spesso disegnato): due segmenti
pts1 = []
for y in range(yC + 3, 290, 2):
    seg = im[y, 160:290]
    if seg.min() < 30:
        pts1.append((160 + np.argmin(seg), y))
pts1 = np.array(pts1, float)
theta_SL1 = deg(math.atan(np.polyfit(pts1[:, 0], pts1[:, 1], 1)[0]))
pts2 = []
for y in range(296, 368, 2):
    seg = im[y, 215:290]
    if seg.min() < 40:
        pts2.append((215 + np.argmin(seg), y))
pts2 = np.array(pts2, float)
theta_SL12 = deg(math.atan(1.0 / np.polyfit(pts2[:, 1], pts2[:, 0], 1)[0]))
fine_SL = pts2[-1]
# F: intersezione dei due segmenti
p1 = np.polyfit(pts1[:, 0], pts1[:, 1], 1)
p2i = np.polyfit(pts2[:, 1], pts2[:, 0], 1)
yF = brentq(lambda y: np.polyval(p1, np.polyval(p2i, y)) - y, 280, 320)
F = np.array([np.polyval(p2i, yF), yF])


def creste_raggio(ang, r0=20, r1=480, prom=12):
    """Creste luminose lungo il raggio uscente da C con angolo ang (gradi, antiorario)."""
    t = rad(ang)
    r = np.arange(r0, r1, 0.5)
    x = C[0] + r * math.cos(t)
    y = C[1] - r * math.sin(t)
    ok = (x > 0) & (x < W_img - 1) & (y > y_parete + 1) & (y < H - 1)
    r, x, y = r[ok], x[ok], y[ok]
    v = map_coordinates(g, [y, x], order=1)
    pk, pr = find_peaks(v, prominence=prom, distance=10)
    return r[pk], pr["prominences"]


# urto di ricompressione RC: la cresta più interna (oltre r = 150 px) fra -50° e -2°
RC_pts = []
for ang in np.arange(-50, -1, 4):
    rr, pp = creste_raggio(ang)
    sel = rr[(rr > 150) & (rr < 215)]
    if len(sel):
        r0 = sel[0]
        RC_pts.append((ang, r0, C[0] + r0 * math.cos(rad(ang)), C[1] - r0 * math.sin(rad(ang))))
RC_pts = np.array(RC_pts)
# normale a RC nel tratto inferiore e velocità normale (autosimilarità: v = (r . n)/t)
i0, i1 = 0, 4
tang = RC_pts[i1, 2:] - RC_pts[i0, 2:]
n_RC = np.array([tang[1], -tang[0]])
n_RC = -n_RC / np.linalg.norm(n_RC) if n_RC[0] < 0 else n_RC / np.linalg.norm(n_RC)
ang_n_RC = deg(math.atan2(-n_RC[1], n_RC[0]))       # direzione della normale (antioraria)
rm = 0.5 * (RC_pts[i0, 2:] + RC_pts[i1, 2:]) - C
WRC_su_WB = float(np.dot(rm, n_RC)) / sB_px

# bordo della zona scura che parte da C (testa dell'espansione): luminosità media sui raggi
ang_r = np.arange(20, 80, 1.0)
lum = []
for a_ in ang_r:
    r = np.linspace(40, 120, 60)
    x = C[0] + r * math.cos(rad(a_))
    y = C[1] - r * math.sin(rad(a_))
    lum.append(map_coordinates(gaussian_filter(im, 2), [y, x], order=1).mean())
lum = np.array(lum)
meta = 0.5 * (lum.min() + lum.max())
ang_testa_oss = float(np.interp(meta, lum, ang_r))   # lum crescente con l'angolo

print("=== MISURE SULL'IMMAGINE ===")
print(f"C = ({C[0]:.0f}, {C[1]:.0f}) px, parete superiore y = {y_parete}")
print(f"parete inclinata: {ang_parete_vert:.1f}° dalla verticale -> deflessione {deflessione_parete:.1f}° (dato 105°)")
print(f"x_A = {xA_px:.0f} px, s_B = {sB_px:.0f} px, rapporto = {rapporto_misurato:.2f} (dato 2,5)")
print(f"SL: primo tratto {theta_SL1:.1f}° (dato 44°), secondo tratto {theta_SL12:.1f}° -> theta2 = {theta_SL12 - theta_SL1:.1f}° (dato 17°)")
print(f"F = ({F[0]:.0f}, {F[1]:.0f}), fine SL = ({fine_SL[0]:.0f}, {fine_SL[1]:.0f})")
print(f"D osservato (scostamento > 3 px dalla verticale per A): quota {hD_oss:.2f} x_A sopra C")
print(f"bordo della zona scura da C (metà del salto di luminosità): {ang_testa_oss:.0f}°")
print("RC (angolo da C, distanza/s_B):", [(int(a_), round(r_ / sB_px, 2)) for a_, r_, _, _ in RC_pts])
print(f"normale a RC (tratto inferiore): {ang_n_RC:.0f}°, W_RC,n / W_B = {WRC_su_WB:.2f}")
print("incertezze: angoli ±1-2° (tratti disegnati spessi 3-4 px), distanze ±2 px su 160-400 px (≈1%);"
      " la quota di D osservata dipende dalla soglia (2-4 px => 0,30-0,40 x_A).")

# =============================================================================
# 2. STATI
# =============================================================================
print("\n=== STATI ===")
# stato 1: aria in quiete
s1 = gd.State.da_pT(G, 1e5, 293.0, u=0.0)
a1 = s1.a
# stato 2: dietro l'urto piano (F §7, urto mobile u+a con M_s = 4)
U2 = gd.urto_mobile(s1, "u+a", Ms=4.0)
s2 = U2.valle
W = U2.w
u2, a2, p2, T2 = s2.u, s2.a, s2.p, s2.T
M2 = u2 / a2
rn = gd.urto_normale(4.0, G)
print(f"1: p = {s1.p:.3g} Pa, T = {s1.T} K, a1 = {a1:.1f} m/s, rho1 = {s1.rho:.4f}")
print(f"urto piano: W = M_s a1 = {W:.0f} m/s; p2/p1 = {rn['p2_p1']:.3f}, rho2/rho1 = {rn['rho2_rho1']:.4f},"
      f" T2/T1 = {rn['T2_T1']:.4f}")
print(f"2: u2 = {u2:.0f} m/s, a2 = {a2:.1f} m/s, M2 = {M2:.4f}, p2 = {p2:.4g} Pa, T2 = {T2:.0f} K")
p02 = p2 / p_p0(M2)
a02 = a2 / a_a0(M2)
print(f"   p02 = {p02:.4g} Pa, a02 = {a02:.1f} m/s, p2/p02 = {p_p0(M2):.4f}")
Ms_lim = brentq(lambda m: gd.urto_mobile(s1, 'u+a', Ms=m).valle.u / gd.urto_mobile(s1, 'u+a', Ms=m).valle.a - 1, 1.5, 3)
print(f"   M2 = 1 per M_s = {Ms_lim:.3f}: con M_s = 4 la corrente 2 è supersonica")

# piede B: urto normale alla parete, gas 1 in quiete; autosimilarità => W_B = W / 2,5
WB = W / 2.5
MsB = WB / a1
UB = gd.urto_mobile(s1, "u+a", w=WB)
s5 = UB.valle
rB = gd.urto_normale(MsB, G)
M5 = s5.u / s5.a
print(f"\nB: W_B = W/2,5 = {WB:.0f} m/s, M_sB = {MsB:.3f}; p5/p1 = {rB['p2_p1']:.3f}, rho5/rho1 = {rB['rho2_rho1']:.3f},"
      f" T5/T1 = {rB['T2_T1']:.3f}")
print(f"5: u5 = {s5.u:.0f} m/s (lungo la parete, verso B), a5 = {s5.a:.0f} m/s, M5 = {M5:.3f}, p5 = {s5.p:.3g} Pa")
CE_CB = s5.u / WB
print(f"E: la particella che era in C all'istante 0 scorre sulla parete con u5 => CE/CB = u5/W_B = {CE_CB:.3f}")

# testa dell'espansione: cerchio acustico emesso da C a t = 0 (centro u2 t, raggio a2 t), in unità W t
xc, Rc = u2 / W, a2 / W
hD = math.sqrt(Rc**2 - (1 - xc)**2)
mu2 = mu(M2)
LT = math.sqrt(xc**2 - Rc**2)                      # distanza C-T (tangenza)
T_pt = (LT * math.cos(rad(mu2)), LT * math.sin(rad(mu2)))
print(f"\ntesta: mu2 = {mu2:.2f}°, (W - u2) = {W - u2:.0f} m/s < a2 = {a2:.0f} m/s")
print(f"   cerchio: centro {xc:.3f}, raggio {Rc:.3f} (unità W t); D a quota {hD:.3f} x_A; "
      f"T = ({T_pt[0]:.3f}, {T_pt[1]:.3f})")
print(f"   in pixel: quota D prevista = {hD * xA_px:.0f} px sopra C (osservata ~{hD_oss * xA_px:.0f})")

# espansione in C (F §11, onda semplice C+: theta + nu costante)
th1, th2 = 44.0, 17.0
nu2 = nu(M2)
nu3 = nu2 + th1
M3 = M_da_nu(nu3)
p3 = p2 * p_p0(M3) / p_p0(M2)
a3 = a2 * a_a0(M3) / a_a0(M2)
V3 = M3 * a3
mu3 = mu(M3)
print(f"\n3: nu2 = {nu2:.2f}°, nu3 = {nu3:.2f}°, M3 = {M3:.3f}, p3 = {p3:.4g} Pa (p3/p1 = {p3 / s1.p:.3f}),"
      f" a3 = {a3:.0f}, V3 = {V3:.0f} m/s, mu3 = {mu3:.2f}°")
print(f"   ventaglio in C: testa a +{mu2:.1f}°, coda a {-th1 + mu3:.1f}°; p3/p2 = {p3 / p2:.4f}")
# espansione in F
nu4 = nu3 + th2
M4 = M_da_nu(nu4)
p4 = p2 * p_p0(M4) / p_p0(M2)
a4 = a2 * a_a0(M4) / a_a0(M2)
V4 = M4 * a4
mu4 = mu(M4)
print(f"4: nu4 = {nu4:.2f}°, M4 = {M4:.3f}, p4 = {p4:.4g} Pa (p4/p3 = {p4 / p3:.4f}), a4 = {a4:.0f}, V4 = {V4:.0f},"
      f" mu4 = {mu4:.2f}°")
print(f"   ventaglio in F: testa a {-th1 + mu3:.1f}°, coda a {-(th1 + th2) + mu4:.1f}°")
# verifica con gasdyn.espansione_pm
o3 = gd.espansione_pm(gd.State.da_M(G, M2, p=p2, a=a2), th1, "C+")
o4 = gd.espansione_pm(o3.valle, th2, "C+")
print(f"   controllo espansione_pm: M3 = {o3.valle.M:.3f}, p3 = {o3.valle.p:.4g}; M4 = {o4.valle.M:.3f}, p4 = {o4.valle.p:.4g}")

# parete attaccata (ipotetica): rotazione di 105°
nu_w = nu2 + 105.0
M_w = M_da_nu(nu_w)
p_w = p2 * p_p0(M_w) / p_p0(M2)
print(f"\nparete attaccata (ipotetica): nu = {nu_w:.1f}° < nu_max = {gd.nu_max(G):.2f}°, M = {M_w:.1f}, p = {p_w:.3g} Pa")

# urto di ricompressione RC: autosimilarità, W_RC ≈ W_B (verso della corrente), urto quasi normale
WRC = WB
res_RC = {}
for nome, V, a, p in (("3", V3, a3, p3), ("4", V4, a4, p4)):
    Mrel = (V - WRC) / a
    rr = gd.urto_normale(Mrel, G)
    p6 = p * rr["p2_p1"]
    # stato a valle nel riferimento del laboratorio
    v_rel_valle = (V - WRC) / rr["rho2_rho1"]
    a6 = a * rr["a2_a1"]
    u6 = WRC + v_rel_valle
    res_RC[nome] = (Mrel, rr["p2_p1"], p6, u6 / a6, rr["M2"])
    print(f"RC su corrente {nome}: M_rel = ({V:.0f} - {WRC:.0f})/{a:.0f} = {Mrel:.3f}; p6/p{nome} = {rr['p2_p1']:.3f},"
          f" p6 = {p6:.3g} Pa; M_rel valle = {rr['M2']:.3f}; nel laboratorio V6 = {u6:.0f} m/s, M6 = {u6 / a6:.3f}")
# dove RC si annulla: (V - W_RC)/a = 1 dentro il ventaglio di C
M_RC0 = brentq(lambda M: M * a02 * a_a0(M) - WRC - a02 * a_a0(M), 1.01, 3.0)
th_RC0 = -(nu(M_RC0) - nu2)
print(f"RC si annulla dove (V - W_RC)/a = 1: M = {M_RC0:.3f}, theta = {th_RC0:.1f}°, raggio del ventaglio a "
      f"{th_RC0 + mu(M_RC0):.1f}° da C")

# stato asintotico 7: getto sottoespanso dal labbro C con contorno a p = p1
M7 = gd.mach_da_rapporto(p02 / s1.p, G, "p0/p")
nu7 = nu(M7)
th7 = nu7 - nu2
print(f"\n7: p7 = p1, p02/p1 = {p02 / s1.p:.2f} -> M7 = {M7:.3f}, nu7 = {nu7:.2f}°, theta_SL,asint = {th7:.1f}°"
      f" (n = p2/p1 = {p2 / s1.p:.1f})")

# entropie (per la discussione del contatto)
ds2 = rn["ds_R"]
dsB = rB["ds_R"]
print(f"\nDelta s / R: urto piano {ds2:.3f}, piede B {dsB:.3f}")

# =============================================================================
# 3. FIGURE
# =============================================================================
plt.rcParams.update({"font.size": 10, "savefig.dpi": 220})
img_rgb = np.asarray(Image.open(os.path.join(FIG, "schlieren_traccia.jpg")).convert("RGB"))


def px(xi, eta):
    """Coordinate autosimili (unità W t, origine C, eta verso l'alto) -> pixel."""
    return C[0] + xi * xA_px, C[1] - eta * xA_px


def semiretta(ax, P0, ang, L, **kw):
    ax.plot([P0[0], P0[0] + L * math.cos(rad(ang))], [P0[1], P0[1] - L * math.sin(rad(ang))], **kw)


def etichetta(ax, x, y, t, col="yellow", fs=11):
    ax.text(x, y, t, color=col, fontsize=fs, fontweight="bold", ha="center", va="center",
            bbox=dict(boxstyle="circle,pad=0.15", fc="black", ec=col, alpha=0.75))


# --- 3a. figura annotata ---
fig, ax = plt.subplots(figsize=(8.6, 6.6))
ax.imshow(img_rgb)
# parete inclinata e misura delle distanze
ax.plot([C[0], B[0]], [C[1], B[1]], color="lime", lw=1.6)
ax.plot([C[0], C[0] + xA_px], [y_parete + 6] * 2, color="cyan", lw=1.4)
ax.annotate("", xy=(C[0] + xA_px, y_parete + 14), xytext=(C[0], y_parete + 14),
            arrowprops=dict(arrowstyle="<->", color="cyan", lw=1.2))
ax.text(C[0] + 0.5 * xA_px, y_parete + 24, r"$x_A$", color="cyan", fontsize=11, ha="center")
ax.text(B[0] - 22, 0.5 * (C[1] + B[1]) + 25, r"$s_B$", color="lime", fontsize=11)
# testa del ventaglio di C: retta di Mach fino a T, poi arco del cerchio acustico fino a D
xT, yT = px(*T_pt)
ax.plot([C[0], xT], [C[1], yT], "--", color="orange", lw=1.6)
phT = math.atan2(T_pt[1], T_pt[0] - xc)
phD = math.atan2(hD, 1 - xc)
ph = np.linspace(phT, phD, 60)
xa, ya = px(xc + Rc * np.cos(ph), Rc * np.sin(ph))
ax.plot(xa, ya, "--", color="orange", lw=1.6)
xD, yD = px(1.0, hD)
ax.plot(xD, yD, "o", color="orange", ms=6)
ax.text(xD - 14, yD - 2, "D (calc.)", color="orange", fontsize=9, ha="right",
        bbox=dict(fc="black", ec="none", alpha=0.6, pad=1))
# coda del ventaglio di C e ventaglio in F
semiretta(ax, C, -th1 + mu3, 120, color="orange", lw=1.0, ls=":")
semiretta(ax, F, -th1 + mu3, 70, color="magenta", lw=1.0, ls=":")
semiretta(ax, F, -(th1 + th2) + mu4, 70, color="magenta", lw=1.0, ls=":")
# RC misurato
ax.plot(RC_pts[:, 2], RC_pts[:, 3], ".", color="red", ms=5)
# E previsto
Ecalc = C + CE_CB * sB_px * dir_parete
ax.plot(Ecalc[0], Ecalc[1], "s", color="cyan", ms=6)
ax.text(Ecalc[0] - 8, Ecalc[1] + 4, "E (calc.)", color="cyan", fontsize=9, ha="right",
        bbox=dict(fc="black", ec="none", alpha=0.6, pad=1))
ax.plot(B[0], B[1], "o", color="lime", ms=6)
# angoli misurati
ax.text(C[0] + 70, C[1] - 6, f"{theta_SL1:.0f}°", color="white", fontsize=9)
ax.text(F[0] + 50, F[1] - 8, f"+{theta_SL12 - theta_SL1:.0f}°", color="white", fontsize=9)
ax.text(C[0] - 12, C[1] + 30, f"{deflessione_parete:.0f}°", color="lime", fontsize=9, ha="right")
# stati
for (x, y, t) in ((60, 140, "2"), (610, 140, "1"), (215, 468, "1"), (245, 262, "3"), (350, 380, "6"),
                  (182, 392, "5"), (196, 318, "7")):
    etichetta(ax, x, y, t)
etichetta(ax, 262, 348, "4", fs=8)
ax.set_xlim(0, W_img)
ax.set_ylim(H, 0)
ax.axis("off")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "annotata.png"))
plt.close(fig)

# --- 3b. linee di corrente (riferimento del laboratorio, istante della foto) ---
# Traccia di DS e CD letta sulla figura (punti scelti a vista sulla cresta luminosa): serve solo per far
# arrivare le linee di corrente su DS con la direzione della sua normale (gas fermo davanti a DS ⇒ dietro
# DS la velocità è normale al fronte).
DS_px = np.array([(570, 10), (569, 60), (567, 125), (556, 200), (541, 258), (520, 300), (490, 345), (452, 385),
                  (405, 418), (350, 440), (290, 450), (235, 446), (189, 433), (151, 414), (127, 402)], float)


def normale_DS(Pp):
    """Normale uscente (verso il gas 1) a DS nel punto della polilinea più vicino a Pp."""
    d = np.hypot(*(DS_px - Pp).T)
    i = int(np.clip(np.argmin(d), 1, len(DS_px) - 2))
    t = DS_px[i + 1] - DS_px[i - 1]
    n = np.array([-t[1], t[0]])
    if np.dot(n, DS_px[i] - C) < 0:
        n = -n
    return DS_px[i], n / np.linalg.norm(n)


def cubica(P0, d0, P3, d3, L0, L3, n=80):
    P0, P3 = np.asarray(P0, float), np.asarray(P3, float)
    P1 = P0 + L0 * np.asarray(d0, float)
    P2 = P3 - L3 * np.asarray(d3, float)
    t = np.linspace(0, 1, n)[:, None]
    out = (1 - t)**3 * P0 + 3 * (1 - t)**2 * t * P1 + 3 * (1 - t) * t**2 * P2 + t**3 * P3
    return out[:, 0], out[:, 1]


def vers(ang):
    return np.array([math.cos(rad(ang)), -math.sin(rad(ang))])


fig, ax = plt.subplots(figsize=(8.6, 6.6))
ax.imshow(img_rgb)
col = "red"


def freccia(xs_, ys_, f=0.5, lw=1.6):
    ax.plot(xs_, ys_, color=col, lw=lw)
    i = int(f * (len(xs_) - 1))
    ax.annotate("", xy=(xs_[i + 1], ys_[i + 1]), xytext=(xs_[i - 1], ys_[i - 1]),
                arrowprops=dict(arrowstyle="-|>", color=col, lw=1.3, mutation_scale=14))


def linea(y0, x_curva, P_DS, L0, L3, f=0.6):
    """Tratto rettilineo nello stato 2 fino a x_curva, poi raccordo fino a DS con direzione normale."""
    Q, n = normale_DS(np.array(P_DS, float))
    xr = np.linspace(120, x_curva, 20)
    xb, yb = cubica((x_curva, y0), (1, 0), Q, n, L0, L3)
    xs_ = np.concatenate([xr, xb[1:]])
    ys_ = np.concatenate([np.full(20, y0), yb[1:]])
    freccia(xs_, ys_, f)


# sopra il cerchio acustico: nessuna perturbazione, rettilinea fino all'urto piano
freccia(np.linspace(120, xA - 3, 40), np.full(40, 32.0), 0.45)
# attraverso il ventaglio di C (deviazione crescente), poi verso DS
linea(110, 345, (556, 200), 70, 50)
linea(170, 240, (520, 300), 130, 50)
linea(212, 180, (452, 385), 150, 90)
# lungo SL (C-F-fine di SL), poi attraverso RC e CD fino a DS
Q, n = normale_DS(np.array((330, 445)))
xs_ = np.concatenate([np.linspace(120, C[0], 15), [F[0] + 3], [fine_SL[0] + 6]])
ys_ = np.concatenate([np.full(15, C[1] - 4), [F[1] - 3], [fine_SL[1] - 4]])
xb, yb = cubica((fine_SL[0] + 6, fine_SL[1] - 4), vers(-62), Q, n, 40, 40)
freccia(np.concatenate([xs_, xb[1:]]), np.concatenate([ys_, yb[1:]]), 0.25)
# ricircolo nella bolla fra la parete CE e SL (verso orario: giù lungo SL, su lungo la parete)
u_ax = np.array([math.cos(rad(80)), math.sin(rad(80))])
v_ax = np.array([math.sin(rad(80)), -math.cos(rad(80))])
tt = np.linspace(2 * math.pi, 0, 120)
cen = np.array([190.0, 309.0])
ell = cen[None, :] + 44 * np.cos(tt)[:, None] * u_ax + 15 * np.sin(tt)[:, None] * v_ax
ax.plot(ell[:, 0], ell[:, 1], color=col, lw=1.3)
for k in (28, 88):
    ax.annotate("", xy=ell[k + 1], xytext=ell[k - 1],
                arrowprops=dict(arrowstyle="-|>", color=col, lw=1.2, mutation_scale=12))
# lungo la parete da E verso B (gas 5, u5 parallela alla parete)
Ecalc = C + CE_CB * sB_px * dir_parete
P0 = Ecalc + 4 * dir_parete + np.array([3, 0])
P1 = B - 3 * dir_parete + np.array([3, 0])
freccia(np.linspace(P0[0], P1[0], 30), np.linspace(P0[1], P1[1], 30), 0.55)
# gas messo in moto da DS (fra CD e DS): dalla zona di E verso DS, ortogonali a DS
for Pstart, dstart, P_DS, L0, L3 in (((205, 362), vers(-82), (189, 433), 25, 25),
                                     ((262, 392), vers(-88), (268, 449), 20, 20)):
    Q, n = normale_DS(np.array(P_DS, float))
    xb, yb = cubica(Pstart, dstart, Q, n, L0, L3)
    freccia(xb, yb, 0.55, lw=1.3)
ax.set_xlim(0, W_img)
ax.set_ylim(H, 0)
ax.axis("off")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "linee_corrente.png"))
plt.close(fig)

# --- 3c. piano p-theta (espansione isentropica dallo stato 2), lineare e logaritmico ---
th_c = np.linspace(0, 105, 400)
p_c = np.array([p2 * p_p0(M_da_nu(nu2 + t)) / p_p0(M2) for t in th_c]) / 1e5
pts_pt = [("2", 0.0, p2), ("3", -th1, p3), ("4", -(th1 + th2), p4), ("7", -th7, s1.p),
          ("parete", -105.0, p_w)]
fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.3))
for ax, logy in zip(axs, (False, True)):
    ax.plot(-th_c, p_c, color="tab:blue", lw=1.6, label="espansione isentropica da 2")
    for lab, th, p in pts_pt:
        ax.plot(th, p / 1e5, "o", color="k" if lab != "parete" else "gray", ms=5)
        off = (-14, -14) if lab == "7" else ((6, 6) if lab != "parete" else (6, 8))
        ax.annotate(lab, (th, p / 1e5), textcoords="offset points", xytext=off, fontsize=10)
    ax.axhline(s1.p / 1e5, color="tab:green", ls="--", lw=1, label=r"$p_1$ (quiete)")
    ax.axhline(s5.p / 1e5, color="tab:orange", ls=":", lw=1.2, label=r"$p_5$ (piede B)")
    ax.axvline(-105.0, color="0.5", ls="-.", lw=0.8)
    ax.set_xlabel(r"$\theta$ [°]")
    ax.set_ylabel(r"$p$ [bar]")
    ax.set_xlim(-110, 5)
    ax.grid(alpha=0.3)
    if logy:
        ax.set_yscale("log")
        ax.set_ylim(2e-6, 40)
        ax.set_title("scala logaritmica")
    else:
        ax.set_ylim(0, 20)
        ax.set_title("scala lineare")
axs[0].legend(loc="upper left", fontsize=8.5)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "p_theta.pdf"))
plt.close(fig)

# --- 3d. schema autosimile (unità W t) ---
fig, ax = plt.subplots(figsize=(7.6, 5.6))
hw = (C[1] - y_parete) / xA_px
Lw = 0.9 / dir_parete[1]
xw, yw = Lw * dir_parete[0], -Lw * dir_parete[1]
ax.fill([-0.5, 0, xw, -0.5], [0, 0, yw, yw], color="0.8")
ax.plot([-0.5, 0, xw], [0, 0, yw], "k", lw=2)
ax.plot([-0.5, 1.15], [hw, hw], "k", lw=2)
ax.plot([1, 1], [hD, hw], color="tab:red", lw=2)
ax.text(1.03, 0.5 * (hD + hw) + 0.02, "urto piano", color="tab:red", va="center", fontsize=9)
ph = np.linspace(0, 2 * math.pi, 300)
ax.plot(xc + Rc * np.cos(ph), Rc * np.sin(ph), ":", color="0.5", lw=1)
ax.plot([0, T_pt[0]], [0, T_pt[1]], "--", color="tab:orange", lw=1.6)
ph = np.linspace(phT, phD, 60)
ax.plot(xc + Rc * np.cos(ph), Rc * np.sin(ph), "--", color="tab:orange", lw=1.6,
        label="testa dell'espansione (retta di Mach + cerchio acustico)")
ax.plot(xc, 0, "+", color="0.4")
ax.text(xc, -0.05, r"$u_2t$", ha="center", va="top", fontsize=9, color="0.3")
for (x_, y_, t) in ((1, hD, "D"), (T_pt[0], T_pt[1], "T"), (0, 0, "C")):
    ax.plot(x_, y_, "o", color="k", ms=4)
    ax.text(x_ + (0.02 if t != "D" else -0.07), y_ + (0.03 if t != "D" else -0.06), t, fontsize=10)
xB_s, yB_s = 0.4 * dir_parete[0], -0.4 * dir_parete[1]
xE_s, yE_s = 0.4 * CE_CB * dir_parete[0], -0.4 * CE_CB * dir_parete[1]
ax.plot(xB_s, yB_s, "o", color="tab:green")
ax.text(xB_s + 0.03, yB_s, r"B: $W_B t = 0{,}4\,Wt$", fontsize=9, color="tab:green", va="center")
ax.plot(xE_s, yE_s, "s", color="tab:cyan")
ax.text(xE_s + 0.03, yE_s, r"E: $u_5 t$", fontsize=9, color="tab:cyan", va="center")
semiretta_xy = lambda P0, ang, L, **kw: ax.plot([P0[0], P0[0] + L * math.cos(rad(ang))],
                                                [P0[1], P0[1] + L * math.sin(rad(ang))], **kw)
semiretta_xy((0, 0), -th1 + mu3, 0.45, color="tab:orange", lw=1, ls=":")
semiretta_xy((0, 0), -th1, 0.2, color="k", lw=2.2)
Fs = (0.2 * math.cos(rad(-th1)), 0.2 * math.sin(rad(-th1)))
semiretta_xy(Fs, -(th1 + th2), 0.17, color="k", lw=2.2)
ax.text(Fs[0] + 0.02, Fs[1] + 0.02, "F", fontsize=10)
ax.text(0.17, -0.07, "SL", fontsize=9)
ax.set_aspect("equal")
ax.set_xlim(-0.5, 1.25)
ax.set_ylim(-0.85, hw + 0.05)
ax.set_xlabel(r"$x/(Wt)$")
ax.set_ylabel(r"$y/(Wt)$")
ax.legend(loc="lower right", fontsize=8.5)
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "schema.pdf"))
plt.close(fig)
print("\nfigure salvate in", FIG)
