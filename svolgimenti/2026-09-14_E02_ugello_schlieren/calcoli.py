"""Prova del 14/09/2026, esercizio 2: distribuzione del numero di Mach in un ugello (schlieren).

Riproduce tutte le misure sull'immagine, tutti i numeri e le figure dello svolgimento.
Uso:  python calcoli.py   (dalla cartella dello svolgimento)

L'immagine è una foto dello schermo: prima di misurare si raddrizza la prospettiva
con l'omografia che porta i quattro spigoli della cornice nera in un rettangolo.
Coordinate immagine (raddrizzata): x verso destra, y verso il basso, in pixel.
Coordinate fisiche: X = (x - x_gola)/h*, Y = (y_asse - y)/h*, angoli antiorari dall'asse.
"""
import os
import sys
import math

import numpy as np
from PIL import Image, ImageOps
from scipy.ndimage import gaussian_filter, center_of_mass, map_coordinates
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

QUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(QUI, "..", ".."))
import gasdyn as gd  # noqa: E402

FIG = os.path.join(QUI, "figure")
aria = gd.Gas(1.4)
deg, rad = math.degrees, math.radians


def mu(M):
    return deg(math.asin(1.0 / M))


def nu(M):
    return gd.prandtl_meyer(M, aria)


def M_da_nu(v):
    return gd.mach_da_pm(v, aria)


def p_p0(M):
    return 1.0 / gd.isentropiche(M, aria)["p0_p"]


# =============================================================================
# 1. RADDRIZZAMENTO DELLA FOTO
# =============================================================================
foto = ImageOps.exif_transpose(Image.open(os.path.join(FIG, "originale.jpg")))
g0 = np.asarray(foto.convert("L")).astype(float)
# cornice nera dell'immagine d'esame: pixel scuri nella zona della figura
X0, Y0, X1, Y1 = 490, 1430, 1880, 2620
ys, xs = np.nonzero(g0[Y0:Y1, X0:X1] < 70)
xs, ys = xs + X0, ys + Y0
spigoli = []
for f in (-(xs + ys), xs - ys, xs + ys, ys - xs):          # TL, TR, BR, BL
    i = np.argmax(f)
    spigoli.append((xs[i], ys[i]))
src = np.array(spigoli, float)
Wr = 0.5 * (np.hypot(*(src[1] - src[0])) + np.hypot(*(src[2] - src[3])))
Hr = 0.5 * (np.hypot(*(src[3] - src[0])) + np.hypot(*(src[2] - src[1])))
dst = np.array([[0, 0], [Wr, 0], [Wr, Hr], [0, Hr]])


def omografia(s, d):
    A = []
    for (x, y), (u, v) in zip(s, d):
        A.append([x, y, 1, 0, 0, 0, -u * x, -u * y, -u])
        A.append([0, 0, 0, x, y, 1, -v * x, -v * y, -v])
    h = np.linalg.svd(np.array(A))[2][-1]
    return (h / h[-1]).reshape(3, 3)


Hinv = np.linalg.inv(omografia(src, dst))
rett = foto.transform((int(Wr), int(Hr)), Image.PERSPECTIVE, tuple(Hinv.flatten()[:8]), Image.BICUBIC)
rett.save(os.path.join(FIG, "raddrizzata.png"))
im = np.asarray(rett.convert("L")).astype(float)
g1 = gaussian_filter(im, 1.0)
g15 = gaussian_filter(im, 1.5)
print(f"Spigoli cornice (foto): {spigoli}; immagine raddrizzata {Wr:.0f} x {Hr:.0f} px")

# controllo: la finestra circolare deve restare circolare (arco sinistro e arco destro)
arco = []
for y in range(150, 960, 6):
    c = np.where(g15[y, 15:330] > 45)[0]
    if len(c):
        arco.append((15 + c[0], y))
arco = np.array(arco, float)
A_ = np.c_[2 * arco[:, 0], 2 * arco[:, 1], np.ones(len(arco))]
for _ in range(2):
    c_ = np.linalg.lstsq(A_, arco[:, 0] ** 2 + arco[:, 1] ** 2, rcond=None)[0]
    cxw, cyw = c_[0], c_[1]
    Rw_ = math.sqrt(c_[2] + cxw ** 2 + cyw ** 2)
    r_ = np.hypot(arco[:, 0] - cxw, arco[:, 1] - cyw) - Rw_
    k_ = np.abs(r_) < 2 * r_.std()
    arco, A_ = arco[k_], A_[k_]
x_destro = max(1180 + np.where(g15[y, 1180:1262] > 35)[0][-1] for y in range(520, 600, 10))
diam_or = x_destro - (cxw - Rw_)
print(f"Finestra: raggio dall'arco sinistro {Rw_:.0f} px, diametro orizzontale {diam_or:.0f} px "
      f"-> rapporto {diam_or / (2 * Rw_):.3f} (1 = nessuna distorsione residua)")

# =============================================================================
# 2. MISURE SULL'IMMAGINE
# =============================================================================
# --- fori che individuano la gola
lab = []
for (xa, xb, ya, yb) in [(500, 550, 95, 140), (500, 550, 965, 1010)]:
    cy_, cx_ = center_of_mass(im[ya:yb, xa:xb] > 150)
    lab.append((xa + cx_, ya + cy_))
x_gola = 0.5 * (lab[0][0] + lab[1][0])
incl_fori = deg(math.atan2(lab[1][0] - lab[0][0], lab[1][1] - lab[0][1]))
print(f"Fori: {lab[0][0]:.1f},{lab[0][1]:.1f} e {lab[1][0]:.1f},{lab[1][1]:.1f}; x_gola = {x_gola:.1f}, "
      f"retta dei fori inclinata di {incl_fori:.2f} gradi sulla verticale")

# --- pareti: bordo fra il blocco nero e la corrente
sup, inf = {}, {}
for x in range(555, 906):
    sup[x] = 40 + np.where(g1[40:600, x] > 40)[0][0]
for x in range(555, 836):
    inf[x] = 1031 - np.where(g1[600:1032, x][::-1] > 40)[0][0]
xs_s, ys_s = np.array(list(sup.keys()), float), np.array(list(sup.values()), float)
xs_i, ys_i = np.array(list(inf.keys()), float), np.array(list(inf.values()), float)
k_i = xs_i >= 558                       # la punta della cuspide inferiore è troppo sottile
xs_i, ys_i = xs_i[k_i], ys_i[k_i]
c_sup = np.polyfit(xs_s, ys_s, 3)
c_inf = np.polyfit(xs_i, ys_i, 3)
h_star = float(np.mean([inf[x] - sup[x] for x in range(555, 561)]))
y_asse = float(np.mean([(inf[x] + sup[x]) / 2 for x in range(560, 836)]))
print(f"Altezza di gola h* = {h_star:.0f} px; asse (mezzeria delle pareti) y = {y_asse:.1f} px "
      f"(scarto della mezzeria {np.ptp([(inf[x] + sup[x]) / 2 for x in range(560, 836)]):.1f} px)")
print(f"Residui dei profili cubici delle pareti: sup {np.std(ys_s - np.polyval(c_sup, xs_s)):.2f} px, "
      f"inf {np.std(ys_i - np.polyval(c_inf, xs_i)):.2f} px")


def theta_sup(x):          # inclinazione della parete superiore [gradi, antiorario]
    return -deg(math.atan(np.polyval(np.polyder(c_sup), x)))


def theta_inf(x):          # parete inferiore (negativa nel divergente)
    return -deg(math.atan(np.polyval(np.polyder(c_inf), x)))


def X(x):
    return (x - x_gola) / h_star


# --- linee rettilinee: retta che massimizza (cresta chiara) o minimizza (linea scura) la
# luminosità media lungo il tratto [xa, xb]; ricerca su inclinazione e posizione
def retta(xc, yc, ang0, xa, xb, modo, dy=12.0, dang=5.0):
    xs_ = np.arange(xa, xb, 1.0)
    best = None
    for ang in np.arange(ang0 - dang, ang0 + dang + 1e-9, 0.1):
        t = math.tan(rad(ang))
        for off in np.arange(-dy, dy + 1e-9, 0.5):
            v = map_coordinates(g15, [yc + off + t * (xs_ - xc), xs_], order=1).mean()
            s = v if modo == "max" else -v
            if best is None or s > best[0]:
                best = (s, t, yc + off - t * xc, v)
    return best[2], best[1], best[3]          # y = a + b x (immagine), luminosità media


a_s, b_s, v_s = retta(900, 353, 40.0, 725, 1080, "max")     # onda dalla parete superiore (cresta chiara)
a_i, b_i, v_i = retta(900, 710, -43.6, 720, 1080, "max")    # onda dalla parete inferiore (cresta chiara)
phi_s = -deg(math.atan(b_s))       # inclinazione fisica (y verso l'alto): negativa, onda C-
phi_i = -deg(math.atan(b_i))       # positiva, onda C+
x_P = (a_i - a_s) / (b_s - b_i)
y_P = a_s + b_s * x_P
apertura = phi_i - phi_s
print(f"Onda superiore: inclinazione {phi_s:.1f} gradi (luminosità media {v_s:.0f}); "
      f"onda inferiore: {phi_i:.1f} gradi ({v_i:.0f})")
print(f"Angolo fra le due onde {apertura:.1f} gradi; incrocio P a x = {x_P:.0f}, y = {y_P:.0f} px "
      f"({(y_asse - y_P):.0f} px sopra la mezzeria)")

# prolungamenti oltre l'incrocio (linee scure)
a_s2, b_s2, _ = retta(1160, 563, 38.0, 1110, 1205, "min", dy=10)
a_i2, b_i2, _ = retta(1155, 474, -43.0, 1110, 1200, "min", dy=10)
print(f"Oltre l'incrocio: prolungamento dell'onda superiore {-deg(math.atan(b_s2)):.1f} gradi, "
      f"dell'onda inferiore {-deg(math.atan(b_i2)):.1f} gradi")

# piedi delle onde sulle pareti (intersezione retta-profilo)
x_Ws = brentq(lambda x: a_s + b_s * x - np.polyval(c_sup, x), 600, 760)
x_Wi = brentq(lambda x: a_i + b_i * x - np.polyval(c_inf, x), 600, 760)
th_Ws, th_Wi = theta_sup(x_Ws), theta_inf(x_Wi)
print(f"Piede superiore W: x = {x_Ws:.0f} px (X = {X(x_Ws):.3f}), parete a {th_Ws:.1f} gradi")
print(f"Piede inferiore W': x = {x_Wi:.0f} px (X = {X(x_Wi):.3f}), parete a {th_Wi:.1f} gradi")

# linea quasi orizzontale con frangia irregolare che parte dal piede inferiore
a_l, b_l, v_l = retta(930, 918, 0.0, 760, 1110, "min", dy=10)
th_strato = -deg(math.atan(b_l))
print(f"Linea dal piede inferiore verso valle: inclinazione {th_strato:+.1f} gradi, "
      f"y = {a_l + b_l * x_Wi:.0f} px sotto il piede (parete a y = {np.polyval(c_inf, x_Wi):.0f})")

# fine del tratto di parete misurabile
x_fin = 835
print(f"Parete visibile fino a X = {X(x_fin):.2f} (inf) e {X(905):.2f} (sup); "
      f"inclinazione massima letta {theta_sup(900):.0f} gradi")

# --- estremo destro della finestra: bordo netto chiaro -> scuro della zona chiara, riga per riga
# (gradiente orizzontale più negativo fra x = 1170 e 1250), poi retta x = a + b y
gE = gaussian_filter(im, 2.5)
gEx = np.gradient(gE, axis=1)


def bordo(y0, y1):
    ys_ = np.arange(y0, y1 + 1, 2)
    xs_ = np.array([1170 + np.argmin(gEx[y, 1170:1250]) for y in ys_])
    b, a = np.polyfit(ys_, xs_, 1)
    return a, b, np.std(xs_ - (a + b * ys_))


a_g, b_g, r_g = bordo(430, 590)        # tratto quasi verticale vicino all'asse: gambo di Mach
a_Fs, b_Fs, r_Fs = bordo(345, 400)     # tratto che sale dal punto triplo superiore
a_Fi, b_Fi, r_Fi = bordo(625, 715)     # tratto che scende dal punto triplo inferiore


def incl_asse(b):                      # angolo acuto fra la retta x = a + b y e l'asse
    return 90.0 if abs(b) < 1e-9 else deg(math.atan(1.0 / abs(b)))


# punti tripli: incontro dei prolungamenti delle onde oltre P con il gambo
x_T = (a_g + b_g * a_i2) / (1 - b_g * b_i2)      # T, sopra: arriva l'onda partita da W'
y_T = a_i2 + b_i2 * x_T
x_Ti = (a_g + b_g * a_s2) / (1 - b_g * b_s2)     # T', sotto: arriva l'onda partita da W
y_Ti = a_s2 + b_s2 * x_Ti
sig_Fs, sig_Fi = incl_asse(b_Fs), incl_asse(b_Fi)
sig_F = 0.5 * (sig_Fs + sig_Fi)
print(f"Gambo di Mach: x = {a_g + b_g * y_asse:.0f} px (X = {X(a_g + b_g * y_asse):.2f}), inclinato di {incl_asse(b_g):.1f} "
      f"gradi sull'asse (residuo {r_g:.1f} px)")
print(f"Punti tripli: T a X = {X(x_T):.2f}, Y = {(y_asse - y_T) / h_star:+.2f}; T' a X = {X(x_Ti):.2f}, "
      f"Y = {(y_asse - y_Ti) / h_star:+.2f}")
print(f"Urti che partono dai punti tripli: da T {sig_Fs:.1f} gradi, da T' {sig_Fi:.1f} gradi sull'asse "
      f"(residui {r_Fs:.1f}, {r_Fi:.1f} px); media {sig_F:.1f}")

# =============================================================================
# 3. CALCOLO
# =============================================================================
# a) separazione: si usa la metà inferiore, dove lo strato di scorrimento è netto
th_sep = th_strato - th_Wi             # deviazione = angolo fra strato separato e tangente alla parete
sig_sep = sig_Wi = phi_i - th_Wi       # angolo fra urto e tangente alla parete
M_sep = gd.mach_da_theta_sigma(th_sep, sig_sep, aria)          # (F §12, inversa θ-σ-M)
print(f"\n--- Separazione (metà inferiore): strato a {th_sep:.1f} gradi dalla tangente alla parete, "
      f"urto a {sig_sep:.1f} gradi -> M_1 = {M_sep:.2f}")
print(f"   (metà superiore: urto a {th_Ws - phi_s:.1f} gradi dalla parete; strato non leggibile)")
s1W = gd.State.da_M(aria, M_sep, theta=th_Wi, p=p_p0(M_sep), a=1.0, nome="1")
uW = gd.urto_obliquo(s1W, sigma=sig_sep, famiglia="C+", nome="2'")
thm, sig_thm = gd.theta_max(M_sep, aria)
print(f"Urto di separazione: Mn1 = {uW.Mn1:.3f}, deviazione {uW.deflessione:.1f}, p2/p1 = {uW.p2_p1:.2f}, "
      f"M2 = {uW.valle.M:.2f}, theta2 = {uW.valle.theta:+.1f}, p1/p0 = {s1W.p:.3f}, p2/p0 = {uW.valle.p:.3f}")
print(f"   theta_max({M_sep:.2f}) = {thm:.1f}, sigma a theta_max = {sig_thm:.1f}: ramo debole")
M_PM = M_da_nu(-th_Wi)
print(f"Verifica a parete (onda semplice, nu = |theta_w| = {-th_Wi:.1f}): M = {M_PM:.2f}; "
      f"metà superiore nu = {th_Ws:.1f} -> M = {M_da_nu(th_Ws):.2f}")
p_B = uW.valle.p
print(f"Zona separata: p_B/p0 = {p_B:.3f}; p_1/p_B = {s1W.p / p_B:.2f}")

# b) asse: l'onda arriva e riparte con la stessa inclinazione -> onda di Mach
angoli_asse = [-phi_s, phi_i, deg(math.atan(b_s2)), -deg(math.atan(b_i2))]
sig_ax = float(np.mean(angoli_asse))
M_ax = 1.0 / math.sin(rad(sig_ax))
print(f"\n--- Asse: inclinazioni {', '.join(f'{v:.1f}' for v in angoli_asse)} -> media {sig_ax:.1f} gradi; "
      f"onda di Mach: M = 1/sin = {M_ax:.2f} (con 43 gradi: {1 / math.sin(rad(43)):.2f})")
M_ax = round(1.0 / math.sin(rad(43.0)), 2)      # come nello svolgimento: sigma letto ~43 gradi
s1T = gd.State.da_M(aria, M_ax, theta=0.0, p=p_p0(M_ax), a=1.0, nome="1")
print(f"Stato 1 sull'asse: M = {M_ax:.2f}, p/p0 = {s1T.p:.3f}, mu = {mu(M_ax):.1f}")

# c) punto triplo: gambo di Mach (urto normale) e urto obliquo forte
un = gd.urto_normale(M_ax, aria)
print(f"\n--- Gambo di Mach (urto normale con M1 = {M_ax:.2f}): {un}")
sig_T = round(sig_F)
uF = gd.urto_obliquo(s1T, sigma=sig_T, famiglia="C-", nome="4'")
thmT, sig_thmT = gd.theta_max(M_ax, aria)
print(f"Urto forte da T' (sigma = {sig_T} gradi): Mn1 = {uF.Mn1:.3f}, deviazione {uF.deflessione:.1f}, "
      f"p2/p1 = {uF.p2_p1:.2f}, M2 = {uF.valle.M:.2f}, p2/p0 = {uF.valle.p:.3f}")
print(f"   theta_max({M_ax:.2f}) = {thmT:.1f} a sigma = {sig_thmT:.1f}: sigma > sigma(theta_max) -> ramo forte")
for s_ in (sig_Fi, sig_Fs):
    print(f"   sensibilita': sigma = {s_:.1f} -> deviazione {gd.theta_da_sigma(M_ax, s_, aria):.1f}")
s3 = gd.State.da_M(aria, un["M2"], theta=0.0, p=s1T.p * un["p2_p1"], a=1.0, nome="3")
print(f"Stato 3 (dietro il gambo): M = {s3.M:.2f}, p/p0 = {s3.p:.3f}")

# =============================================================================
# 4. FIGURE
# =============================================================================
plt.rcParams.update({"font.size": 9})


def px(Xs, Ys):
    return x_gola + Xs * h_star, y_asse - Ys * h_star


# --- figura annotata
fig, ax = plt.subplots(figsize=(7.2, 6.0))
ax.imshow(rett, cmap="gray")
ax.set_xlim(0, Wr)
ax.set_ylim(Hr, 0)
ax.plot([x_gola, 1240], [y_asse] * 2, ls="-.", color="yellow", lw=0.8)
ax.text(30, y_asse - 10, "asse", color="yellow", fontsize=8)
ax.plot([x_gola, x_gola], [np.polyval(c_sup, 556), np.polyval(c_inf, 558)], color="cyan", lw=1.2)
ax.text(x_gola + 12, y_asse + 200, "gola\n(fori)", color="cyan", fontsize=8)
xs_p = np.linspace(555, 905, 60)
ax.plot(xs_p, np.polyval(c_sup, xs_p), color="lime", lw=1.2)
xs_p = np.linspace(558, 835, 60)
ax.plot(xs_p, np.polyval(c_inf, xs_p), color="lime", lw=1.2)
for (a, b, xa, xb) in [(a_s, b_s, x_Ws, x_P), (a_i, b_i, x_Wi, x_P)]:
    ax.plot([xa, xb], [a + b * xa, a + b * xb], color="red", lw=1.4)
ax.plot([x_P, x_T], [a_i2 + b_i2 * x_P, y_T], color="red", lw=1.0, ls="--")
ax.plot([x_P, x_Ti], [a_s2 + b_s2 * x_P, y_Ti], color="red", lw=1.0, ls="--")
ax.plot([a_g + b_g * y_T, a_g + b_g * y_Ti], [y_T, y_Ti], color="magenta", lw=1.6)
ax.plot([x_T, a_Fs + b_Fs * 335], [y_T, 335], color="magenta", lw=1.4)
ax.plot([x_Ti, a_Fi + b_Fi * 725], [y_Ti, 725], color="magenta", lw=1.4)
ax.plot([x_Wi, 1120], [a_l + b_l * x_Wi, a_l + b_l * 1120], color="orange", lw=1.2, ls="--")
for (x, y, t, dx, dy) in [(x_Ws, np.polyval(c_sup, x_Ws), "W", -10, -22), (x_Wi, np.polyval(c_inf, x_Wi), "W'", -18, 40),
                          (x_P, y_P, "P", -12, -16), (x_T, y_T, "T", -26, -6), (x_Ti, y_Ti, "T'", -32, 22)]:
    ax.text(x + dx, y + dy, t, color="white", fontsize=10, fontweight="bold")
for (x, y, t) in [(800, 560, "1"), (930, 300, "2"), (900, 830, "2'"), (1240, 470, "3"), (1244, 700, "4'"),
                  (1236, 360, "4"), (960, 115, "B"), (1000, 985, "B'")]:
    ax.text(x, y, t, color="yellow", fontsize=9 if t in ("3", "4", "4'") else 11, fontweight="bold", ha="center",
            va="center", bbox=dict(boxstyle="circle,pad=0.12", fc="black", ec="yellow", lw=0.6))
ax.axis("off")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "annotata.png"), dpi=220, bbox_inches="tight", pad_inches=0.02)
plt.close(fig)

# --- dettaglio dell'estremo destro (punto triplo)
fig, ax = plt.subplots(figsize=(3.4, 5.0))
ax.imshow(rett, cmap="gray")
ax.set_xlim(1040, 1262)
ax.set_ylim(780, 300)
ax.plot([x_P, x_T], [a_i2 + b_i2 * x_P, y_T], color="red", lw=1.0, ls="--")
ax.plot([x_P, x_Ti], [a_s2 + b_s2 * x_P, y_Ti], color="red", lw=1.0, ls="--")
ax.plot([a_g + b_g * y_T, a_g + b_g * y_Ti], [y_T, y_Ti], color="magenta", lw=1.4)
ax.plot([x_T, a_Fs + b_Fs * 335], [y_T, 335], color="magenta", lw=1.2)
ax.plot([x_Ti, a_Fi + b_Fi * 725], [y_Ti, 725], color="magenta", lw=1.2)
ax.plot([1040, 1262], [y_asse] * 2, ls="-.", color="yellow", lw=0.8)
for (x, y, t, dx, dy) in [(x_P, y_P, "P", -14, -12), (x_T, y_T, "T", -22, -4), (x_Ti, y_Ti, "T'", -28, 18)]:
    ax.text(x + dx, y + dy, t, color="white", fontsize=10, fontweight="bold")
ax.text(1150, 640, f"σ ≈ {sig_T}°", color="white", fontsize=8)
ax.axis("off")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "dettaglio_T.png"), dpi=220, bbox_inches="tight", pad_inches=0.02)
plt.close(fig)

# --- schema del campo, unità h* (pareti misurate, onde con le inclinazioni misurate)
fig, ax = plt.subplots(figsize=(7.0, 4.6))


def Yw(Xs):
    return (y_asse - np.polyval(c_sup, x_gola + Xs * h_star)) / h_star


Xg = np.linspace(X(555), X(905), 100)
ax.plot(Xg, Yw(Xg), "k", lw=1.6)
ax.plot(Xg, -Yw(Xg), "k", lw=1.6)
Xc = np.linspace(-0.25, X(555), 20)
ax.plot(Xc, 0.5 + 0.35 * Xc ** 2, "k", lw=1.6, ls=(0, (2, 2)))
ax.plot(Xc, -(0.5 + 0.35 * Xc ** 2), "k", lw=1.6, ls=(0, (2, 2)))
ax.axhline(0, color="gray", ls="-.", lw=0.7)
ax.plot([0, 0], [-0.5, 0.5], color="C0", lw=0.8, ls=":")
ys_ = np.linspace(-0.5, 0.5, 41)
ax.plot(0.04 - 0.16 * ys_ ** 2, ys_, color="C0", lw=1.2)
ax.text(0.05, 0.08, "linea\nsonica", color="C0", fontsize=7)
XW_, YW_ = 0.5 * (X(x_Ws) + X(x_Wi)), 0.5 * (Yw(X(x_Ws)) + Yw(X(x_Wi)))
XP_, XT_ = X(x_P), 0.5 * (X(x_T) + X(x_Ti))
YT_ = (XT_ - XP_) * math.tan(rad(sig_ax))
ax.plot([XW_, XP_], [YW_, 0], "r", lw=1.6)
ax.plot([XW_, XP_], [-YW_, 0], "r", lw=1.6)
ax.plot([XP_, XT_], [0, YT_], "r", lw=1.0, ls="--")
ax.plot([XP_, XT_], [0, -YT_], "r", lw=1.0, ls="--")
ax.plot([XT_, XT_], [-YT_, YT_], color="m", lw=1.8)
LF = 0.22
for sgn in (1, -1):
    ax.plot([XT_, XT_ + LF * math.cos(rad(sig_T))], [sgn * YT_, sgn * (YT_ + LF * math.sin(rad(sig_T)))], color="m", lw=1.4)
    ax.plot([XT_, XT_ + 0.12], [sgn * YT_, sgn * (YT_ + 0.02)], color="m", lw=0.8, ls=":")
    ax.plot([XW_, XT_ + 0.1], [sgn * YW_, sgn * (YW_ + 0.012 * (XT_ + 0.1 - XW_))], color="C1", lw=1.2, ls="--")
ax.text(XW_ - 0.03, YW_ + 0.03, "W")
ax.text(XW_ - 0.04, -YW_ - 0.08, "W'")
ax.text(XP_ - 0.03, 0.03, "P")
ax.text(XT_ - 0.07, YT_ + 0.01, "T")
ax.text(XT_ - 0.08, -YT_ - 0.06, "T'")
for (x, y, t) in [(0.45, 0.0, "1"), (0.62, 0.38, "2"), (0.62, -0.38, "2'"), (XT_ + 0.07, 0.0, "3"),
                  (XT_ + 0.14, 0.42, "4"), (XT_ + 0.14, -0.42, "4'"), (0.78, 0.62, "B"), (0.78, -0.62, "B'")]:
    ax.text(x, y, t, fontsize=10, ha="center", va="center", bbox=dict(boxstyle="circle,pad=0.2", fc="white", ec="k", lw=0.6))
ax.set_aspect("equal")
ax.set_xlim(-0.3, XT_ + 0.25)
ax.set_ylim(-0.95, 0.95)
ax.set_xlabel("x / h*")
ax.set_ylabel("y / h*")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "schema.pdf"))
plt.close(fig)

# --- distribuzione del numero di Mach (parete: onda semplice; asse: valori stimati)
fig, ax = plt.subplots(figsize=(6.6, 4.0))
XWi = X(x_Wi)
Xp = np.linspace(X(560), XWi, 40)
Mpar = [M_da_nu(max(-theta_inf(x_gola + z * h_star), 1e-6)) for z in Xp]
ax.plot(Xp, Mpar, "C3", lw=1.6, label="parete (onda semplice, ν = |θ$_w$|)")
ax.plot([-0.03, Xp[0]], [1.0, Mpar[0]], "C3", lw=1.0, ls=":")
ax.plot([XWi, XWi], [M_sep, uW.valle.M], "C3", lw=0.9, ls="--")
ax.plot([XWi], [M_sep], "o", color="C3", ms=4, label=f"monte dell'urto di separazione (θ–σ–M)")
ax.plot([XWi], [uW.valle.M], "s", color="C3", ms=5, label="dietro l'urto di separazione (2')")
XTm = 0.5 * (X(x_T) + X(x_Ti))
ax.plot([0.03, XP_], [1.0, M_ax], "C2", lw=1.2, ls=":", label="asse, andamento qualitativo")
ax.plot([XP_, XTm], [M_ax, M_ax], "C2", lw=1.6, label="asse (angolo di Mach dell'onda in P)")
ax.plot([XP_], [M_ax], "o", color="C2", ms=4)
ax.plot([XTm, XTm], [M_ax, s3.M], "C2", lw=0.9, ls="--")
ax.plot([XTm], [s3.M], "s", color="C2", ms=5, label="dietro il gambo di Mach (3)")
ax.plot([XTm + 0.02], [uF.valle.M], "D", color="C4", ms=5, label="dietro l'urto forte da T' (4')")
for xv, t in ((XWi, "W'"), (XP_, "P"), (XTm, "T")):
    ax.axvline(xv, color="gray", lw=0.6, ls=":")
    ax.text(xv + 0.01, 0.62, t, fontsize=8)
ax.set_xlabel("x / h*  (dalla gola)")
ax.set_ylabel("M")
ax.set_xlim(-0.05, XTm + 0.08)
ax.set_ylim(0.6, 1.8)
ax.grid(alpha=0.3)
ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(0.24, 1.0))
fig.tight_layout()
fig.savefig(os.path.join(FIG, "mach.pdf"))
plt.close(fig)

# --- piano p-theta (p/p0, scala lineare: le pressioni stanno fra 0,2 e 0,7 p0)
fig, ax = plt.subplots(figsize=(6.6, 4.4))
u = gd.polare_p_theta(s1W, famiglia="C+", p_min_rel=0.5)["urto"]
ax.plot(u["theta"], u["p"], color="C3", lw=1.3, label="polare C$^+$ di 1 a parete (W')")
u = gd.polare_p_theta(s1T, famiglia="C-", p_min_rel=0.5)["urto"]
ax.plot(u["theta"], u["p"], color="C2", lw=1.3, label="polare C$^-$ di 1 sull'asse (T')")
th_is = np.linspace(0, 16, 60)
ax.plot(-th_is, [p_p0(M_da_nu(max(t, 1e-6))) for t in th_is], color="gray", lw=1.0, ls=":",
        label="corrente attaccata a parete (ν = |θ|)")
for st, t, off in [(s1W, "1 (W')", (5, -10)), (uW.valle, "2'", (5, 4)), (s1T, "1 (asse)", (5, -10)),
                   (s3, "3", (5, 3)), (uF.valle, "4'", (-14, 4))]:
    ax.plot(st.theta, st.p, "ko", ms=3.5)
    ax.annotate(t, (st.theta, st.p), textcoords="offset points", xytext=off, fontsize=8)
ax.axhline(p_B, color="C1", lw=0.9, ls="--")
ax.text(-21, p_B + 0.008, "zona separata B'", color="C1", fontsize=8)
ax.set_xlabel("θ [°]")
ax.set_ylabel("p / p$_0$")
ax.set_xlim(-22, 12)
ax.grid(alpha=0.3)
ax.legend(fontsize=7, loc="center right", bbox_to_anchor=(1.0, 0.6))
fig.tight_layout()
fig.savefig(os.path.join(FIG, "p_theta.pdf"))
plt.close(fig)
print("\nFigure scritte in", FIG)
