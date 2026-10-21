"""B63-E01 (15/03/2021): cilindro in corrente a M = 10, schlieren.

Riproduce tutte le misure sull'immagine, tutti i numeri e le figure dello svolgimento.
Uso:  python calcoli.py   (dalla cartella dello svolgimento)
"""
import os
import sys
import math

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

QUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(QUI, "..", ".."))
import gasdyn as gd  # noqa: E402

FIG = os.path.join(QUI, "figure")
IMG = os.path.join(FIG, "originale_scansione.jpg")
aria = gd.Gas(1.4)
M_inf = 10.0

# =============================================================================
# 1. MISURE SULL'IMMAGINE (coordinate in pixel, x verso destra, y verso il basso)
# =============================================================================
rgb = np.asarray(Image.open(IMG).convert("RGB")).astype(float)
im = rgb.mean(2)
g = gaussian_filter(im, 1.5)
H, W = im.shape

# --- scala: tacche della barra 0-10-20 mm (colonne con pixel chiari nelle righe 75-92)
sub = im[75:92, 900:1200] > 180
col = np.where(sub.sum(0) > 8)[0] + 900
gruppi = np.split(col, np.where(np.diff(col) > 1)[0] + 1)
tacche = [gr.mean() for gr in gruppi]
px_mm = (tacche[-1] - tacche[0]) / 20.0
px_mm_err = 1.0 / 20.0           # +-1 px sulla lunghezza della barra


def mm(px):
    return px / px_mm


# --- cilindro: anello chiaro sovrapposto al bordo; fit ai minimi quadrati di una circonferenza
cx, cy, R = 337.0, 666.6, 171.2  # stima iniziale
for _ in range(4):
    pts = []
    for ang in np.linspace(3, 177, 300):
        a = math.radians(ang)
        best = None
        for r in np.arange(R - 8, R + 8, 0.25):
            x, y = cx + r * math.cos(a), cy - r * math.sin(a)
            if not (0 <= x < W - 1 and 0 <= y < H - 1):
                continue
            xi, yi = cx + (r - 6) * math.cos(a), cy - (r - 6) * math.sin(a)
            v, vin = im[int(round(y)), int(round(x))], im[int(round(yi)), int(round(xi))]
            if vin < 40 and (best is None or v > best[0]):
                best = (v, x, y)
        if best and best[0] > 170:
            pts.append(best)
    p = np.array(pts)
    A = np.c_[2 * p[:, 1], 2 * p[:, 2], np.ones(len(p))]
    c = np.linalg.lstsq(A, p[:, 1] ** 2 + p[:, 2] ** 2, rcond=None)[0]
    cx, cy = c[0], c[1]
    R = math.sqrt(c[2] + cx ** 2 + cy ** 2)
R_err_px = 1.5                   # semispessore dell'anello disegnato
y_asse = cy                      # asse di simmetria: retta orizzontale per il centro

# --- urto di prua sull'asse: cresta chiara dell'urto nelle righe vicine all'asse
righe = np.arange(int(y_asse) - 6, int(y_asse) + 5)
prof = im[righe, 60:180].mean(0)
x_urto_asse = 60 + np.argmax(gaussian_filter(prof, 1.0))
x_urto_err = 4.0                 # larghezza della cresta chiara
x_naso = cx - R
Delta_px = x_naso - x_urto_asse
Delta_mm = mm(Delta_px)
Delta_err = mm(math.hypot(x_urto_err, R_err_px)) + Delta_mm * px_mm_err / px_mm
R_mm = mm(R)


def traccia(xs, f, half=8, modo="max"):
    """Per ogni colonna x cerca la cresta (max) o il fronte chiaro->scuro (min del gradiente)."""
    out = []
    for x in xs:
        ys = np.arange(int(f(x) - half), int(f(x) + half) + 1)
        pr = g[ys, x]
        if modo == "max":
            y = ys[np.argmax(pr)]
        else:
            y = ys[np.argmin(np.diff(pr))] + 0.5
        out.append((x, y))
    return np.array(out, float)


def retta(pt):
    """Retta ai minimi quadrati y = a x + b; angolo sotto l'orizzontale (gradi) e scarto tipo."""
    a, b = np.polyfit(pt[:, 0], pt[:, 1], 1)
    res = pt[:, 1] - (a * pt[:, 0] + b)
    se = math.sqrt((res ** 2).sum() / (len(pt) - 2) / ((pt[:, 0] - pt[:, 0].mean()) ** 2).sum())
    return a, b, math.degrees(math.atan(a)), math.degrees(se)


# --- strato di scorrimento (contorno della bolla): primo tratto e tratto prima del collo
ss1 = traccia(range(500, 621, 4), lambda x: 565 + 0.286 * (x - 480), 10)
ss2 = traccia(range(620, 741, 4), lambda x: 565 + 0.286 * (x - 480), 10)
ss = np.vstack([ss1, ss2])
a_ss1, b_ss1, ang_ss1, _ = retta(ss1)
a_ss2, b_ss2, ang_ss2, _ = retta(ss2)
a_ss, b_ss, ang_ss, _ = retta(ss)

# --- punto di separazione: intersezione del primo tratto con la circonferenza
A_ = 1 + a_ss1 ** 2
B_ = -2 * cx + 2 * a_ss1 * (b_ss1 - cy)
C_ = cx ** 2 + (b_ss1 - cy) ** 2 - R ** 2
x_sep = (-B_ + math.sqrt(B_ ** 2 - 4 * A_ * C_)) / (2 * A_)
y_sep = a_ss1 * x_sep + b_ss1
theta_s = math.degrees(math.atan2(cy - y_sep, x_sep - cx))   # dal punto di ristagno posteriore

# --- urto di ricompressione (riattacco): fronte chiaro->scuro. L'urto si incurva (coalescenza delle
#     compressioni dello strato che si raddrizza): si misura al piede, nel tratto formato prima
#     dell'incrocio con l'urto di separazione, e dopo l'incrocio.
f_ur = lambda x: 596 - 0.175 * (x - 900) if x < 1100 else 561 - 0.22 * (x - 1100)
ur0 = traccia(range(890, 1001, 3), f_ur, 14, "fronte")
ur1 = traccia(range(1100, 1201, 3), f_ur, 14, "fronte")
ur2 = traccia(range(1262, 1344, 3), f_ur, 14, "fronte")
urt = traccia(range(890, 1341, 3), f_ur, 14, "fronte")
a_ur0, b_ur0, ang_ur0, se_ur0 = retta(ur0)
a_ur1, b_ur1, ang_ur1, se_ur1 = retta(ur1)
a_ur2, b_ur2, ang_ur2, se_ur2 = retta(ur2)
sig_ax_piede = -ang_ur0          # inclinazioni sull'asse (positive verso l'alto)
sig_ax = -ang_ur1                # tratto formato: valore adottato
sig_ax_oltre = -ang_ur2
x_radice_ur = (b_ss2 - b_ur0) / (a_ur0 - a_ss2)       # incontro con il tratto finale dello strato

# --- urto di separazione (linea quasi orizzontale sopra la bolla)
us = traccia(range(500, 881, 4), lambda x: 537, 8)
a_us, b_us, ang_us, se_us = retta(us)
x_incrocio = (b_us - b_ur1) / (a_ur1 - a_us)

# --- urto di prua: fronte chiaro->scuro lungo x, riga per riga
prua = []
xg = x_urto_asse
for y in range(int(y_asse), 140, -3):
    xs = np.arange(max(int(xg - 25), 1), min(int(xg + 25), W - 2))
    d = np.diff(g[y, xs])
    xg = xs[np.argmin(d)] + 0.5
    prua.append((xg, y))
prua = np.array(prua)
# inclinazione locale dell'urto di prua nel tratto più alto tracciato (150 < y < 260 px)
alto = prua[(prua[:, 1] > 150) & (prua[:, 1] < 260)]
a_alto = np.polyfit(alto[:, 1], alto[:, 0], 1)[0]         # dx/dy
sig_prua_alto = math.degrees(math.atan2(1.0, a_alto))      # angolo con l'orizzontale

# valori adottati nello svolgimento: misure arrotondate al grado, come con il goniometro
theta_s_mis, sig_ax_mis = theta_s, sig_ax
theta_s = float(round(theta_s))
sig_ax = float(round(sig_ax))

# incertezze (solo per la sintesi in chat, non per il PDF)
Dth = 18.0                       # inclinazione dello strato a monte del riattacco (adottata)
Dth_err = 1.5
sig_err = 1.0
theta_s_err = 3.0

# =============================================================================
# 2. CALCOLO
# =============================================================================
st1 = gd.State.da_M(aria, M_inf, theta=0.0, p=1.0, a=1.0, nome="1")
iso1 = gd.isentropiche(M_inf, aria)
un = gd.urto_normale(M_inf, aria)                 # urto normale sull'asse (F §5)
st2 = gd.urto_normale_stato(st1, nome="2")
p0p = un["p02_p1"]                               # pressione di ristagno al naso / p1
Mstar1 = st1.M_star

# a) stand-off: misura e correlazione empirica di Billig (fuori formulario)
Delta_Billig_R = 0.386 * math.exp(4.67 / M_inf ** 2)
Delta_Billig = Delta_Billig_R * R_mm

# e) Mach a monte dell'urto di riattacco: relazione theta-sigma-M invertita (F §12)
sigma_R = sig_ax + Dth
M_W = gd.mach_da_theta_sigma(Dth, sigma_R, aria)
assert abs(gd.theta_da_sigma(M_W, sigma_R, aria) - Dth) < 1e-9
sens = {}
for lab, th, sa in [("Dth-", Dth - Dth_err, sig_ax), ("Dth+", Dth + Dth_err, sig_ax),
                    ("sig-", Dth, sig_ax - sig_err), ("sig+", Dth, sig_ax + sig_err),
                    ("piede", Dth, sig_ax_piede)]:
    sens[lab] = gd.mach_da_theta_sigma(th, sa + th, aria)
M_W_min = min(gd.mach_da_theta_sigma(Dth + i * Dth_err, sig_ax + j * sig_err + Dth + i * Dth_err, aria)
              for i in (-1, 1) for j in (-1, 1))
M_W_max = max(gd.mach_da_theta_sigma(Dth + i * Dth_err, sig_ax + j * sig_err + Dth + i * Dth_err, aria)
              for i in (-1, 1) for j in (-1, 1))
stW = gd.State.da_totali(aria, M_W, p0=p0p, a0=st1.a0, theta=-Dth, nome="W")
urR = gd.urto_obliquo(stW, theta=Dth, famiglia="C+")
st6 = urR.valle.con(nome="6")
thmax_W, sigmax_W = gd.theta_max(M_W, aria)

# d) pressione nella bolla = pressione della corrente W che la lambisce
pB = stW.p
pB_range = [gd.State.da_totali(aria, M, p0=p0p, a0=st1.a0).p for M in (M_W_max, M_W_min)]

# c) separazione: compressione isentropica C+ da theta_sep a -Dth
theta_wall_sep = -(90.0 - theta_s)              # direzione della parete al punto di separazione
theta_R = -Dth - theta_wall_sep                  # deviazione della compressione alla separazione
nu_W = gd.prandtl_meyer(M_W, aria)
nu_sep = nu_W + theta_R
M_sep = gd.mach_da_pm(nu_sep, aria)
st_sep = gd.State.da_totali(aria, M_sep, p0=p0p, a0=st1.a0, theta=theta_wall_sep, nome="S")
M_sep_range = [gd.mach_da_pm(gd.prandtl_meyer(M_W_min, aria) + theta_R - theta_s_err - Dth_err, aria),
               gd.mach_da_pm(gd.prandtl_meyer(M_W_max, aria) + theta_R + theta_s_err + Dth_err, aria)]
# alternativa scartata: un solo urto obliquo di deviazione theta_R. M a valle limite per M1 -> inf
M2_lim_urto = gd._M2_obliquo(1e4, gd.sigma_da_theta(1e4, theta_R, aria), aria)

# verifica: urto che realizzasse tutta la deviazione theta_R da S (direzione della linea attesa)
sig_sepS = gd.sigma_da_theta(M_sep, theta_R, aria)
dir_sepS = theta_wall_sep + sig_sepS
M2_sepS = gd.urto_obliquo(st_sep, theta=theta_R, famiglia="C+").valle.M

# b) linea sonica sul cilindro
# (i) dalla catena: theta + nu costante lungo la parete (onda semplice C+) e attraverso la compressione
theta_star = -Dth + nu_W                         # = theta_W + nu_W = theta* + nu(1)
theta_c = 90.0 - theta_star                      # angolo al centro dal ristagno anteriore
theta_c_range = [90.0 - (-Dth + gd.prandtl_meyer(M_W_max, aria)),
                 90.0 - (-Dth + gd.prandtl_meyer(M_W_min, aria))]
theta_c_piede = 90.0 - (-Dth + gd.prandtl_meyer(sens["piede"], aria))
st3 = gd.State.da_totali(aria, 1.0, p0=p0p, a0=st1.a0, theta=theta_star, nome="3")
# (ii) regola di Busemann (F §18): pendenza della parete = deviazione sonica della polare di M_inf
th_son, sig_son = gd.sonico(M_inf, aria)
thmax1, sigmax1 = gd.theta_max(M_inf, aria)
theta_c_B = 90.0 - th_son
# conseguenze della regola di Busemann sulla catena (a ritroso)
nu_W_B = th_son + Dth
M_W_B = gd.mach_da_pm(nu_W_B, aria)
sigma_R_B = gd.sigma_da_theta(M_W_B, Dth, aria)
M_sep_B = gd.mach_da_pm(th_son - theta_wall_sep, aria)
pB_B = gd.State.da_totali(aria, M_W_B, p0=p0p, a0=st1.a0).p
# Newton modificato (Lees), solo come confronto (fuori formulario)
ps_p0 = 1.0 / gd.isentropiche(1.0, aria)["p0_p"]
cos2 = (ps_p0 - 1.0 / p0p) / (1.0 - 1.0 / p0p)
theta_c_N = math.degrees(math.acos(math.sqrt(cos2)))

# =============================================================================
# 3. STAMPA
# =============================================================================
print("=== MISURE ===")
print(f"tacche barra: {[round(t, 1) for t in tacche]} px -> {px_mm:.3f} px/mm")
print(f"cilindro: centro ({cx:.1f}, {cy:.1f}) px, R = {R:.1f} px = {R_mm:.2f} mm, D = {2*R_mm:.1f} mm")
print(f"urto sull'asse x = {x_urto_asse:.1f} px, naso x = {x_naso:.1f} px -> Delta = {Delta_px:.1f} px"
      f" = {Delta_mm:.2f} +- {Delta_err:.2f} mm, Delta/R = {Delta_mm/R_mm:.3f}")
print(f"strato di scorrimento: tratto 1 {ang_ss1:.2f} deg, tratto 2 {ang_ss2:.2f} deg, globale {ang_ss:.2f} deg")
print(f"separazione: ({x_sep:.1f}, {y_sep:.1f}) px, theta_s = {theta_s_mis:.1f} deg dal ristagno posteriore (adottato {theta_s:.0f})")
print(f"urto di riattacco sull'asse: piede {sig_ax_piede:.2f}+-{se_ur0:.2f}, formato {sig_ax_mis:.2f}+-{se_ur1:.2f} (adottato {sig_ax:.0f}),"
      f" oltre l'incrocio {sig_ax_oltre:.2f}+-{se_ur2:.2f} deg; radice sullo strato x = {x_radice_ur:.0f} px,"
      f" incrocio con l'urto di separazione x = {x_incrocio:.0f} px")
print(f"collo (strato prolungato sull'asse) x = {(y_asse - b_ss2)/a_ss2:.0f} px")
print(f"urto di separazione: {-ang_us:.2f} +- {se_us:.2f} deg sull'orizzontale (positivo verso l'alto)")
print(f"urto di prua nella parte alta: {180 - sig_prua_alto:.1f} deg sull'orizzontale; mu(M_inf) = {gd.angolo_mach(M_inf):.2f}")
print("=== CALCOLO ===")
print(f"M* inf = {Mstar1:.4f}; p0/p1 = {iso1['p0_p']:.4g}")
print(f"urto normale: p2/p1 = {un['p2_p1']:.2f}, rho2/rho1 = {un['rho2_rho1']:.4f}, M2 = {un['M2']:.4f},"
      f" p02/p01 = {un['p02_p01']:.5f}, p0'/p1 = {p0p:.2f}")
print(f"a) Delta = {Delta_mm:.2f} mm (misura); Billig Delta/R = {Delta_Billig_R:.4f} -> {Delta_Billig:.2f} mm")
print(f"e) sigma_R = {sigma_R:.2f} deg, theta = {Dth} -> M_W = {M_W:.3f} (intervallo {M_W_min:.2f}-{M_W_max:.2f});"
      f" sensibilita {dict((k, round(v, 3)) for k, v in sens.items())}")
print(f"   theta_max(M_W) = {thmax_W:.2f} (sigma {sigmax_W:.2f}); stato 6: M = {st6.M:.3f}, p6/p1 = {st6.p:.3f}, "
      f"p06/p0W = {urR.p02_p01:.4f}, Mn1 = {urR.Mn1:.3f}, p6/pW = {urR.p2_p1:.3f}")
print(f"d) p_B/p1 = {pB:.4f}  (p/p0(M_W) = {pB/p0p:.5f}); intervallo {[round(x, 3) for x in pB_range]}")
print(f"c) parete alla separazione {theta_wall_sep:.1f} deg, compressione theta_R = {theta_R:.1f} deg;"
      f" nu_W = {nu_W:.2f}, nu_sep = {nu_sep:.2f}, M_sep = {M_sep:.2f} (intervallo {M_sep_range[0]:.1f}-{M_sep_range[1]:.1f});"
      f" p_sep/p1 = {st_sep.p:.2e}; nu_max = {gd.nu_max(aria):.2f}")
print(f"   urto da S con deviazione theta_R: sigma = {sig_sepS:.1f}, direzione {dir_sepS:.1f} deg, M a valle {M2_sepS:.2f}")
print(f"   un solo urto di deviazione {theta_R:.1f} deg: M2 max (M1->inf) = {M2_lim_urto:.3f} < M_W")
print(f"b) catena: theta* = {theta_star:.2f} deg -> theta_c = {theta_c:.1f} deg dal ristagno anteriore;"
      f" intervallo {[round(x, 1) for x in theta_c_range]}; con sigma al piede {theta_c_piede:.1f}; p3/p1 = {st3.p:.2f}")
print(f"   Busemann: theta_s(M=10) = {th_son:.2f} deg (sigma {sig_son:.2f}), theta_max = {thmax1:.2f} ->"
      f" theta_c = {theta_c_B:.1f} deg; a ritroso M_W = {M_W_B:.3f}, sigma_R = {sigma_R_B:.2f}"
      f" (sull'asse {sigma_R_B - Dth:.2f}), M_sep = {M_sep_B:.2f}, p_B/p1 = {pB_B:.3f}")
print(f"   Newton modificato (confronto): theta_c = {theta_c_N:.1f} deg")
print(f"stato 3 (sonico) M* = {st3.M_star:.3f}; W: M* = {stW.M_star:.4f}, S: M* = {st_sep.M_star:.4f}, 6: M* = {st6.M_star:.4f}")

# =============================================================================
# 4. FIGURE
# =============================================================================
plt.rcParams.update({"font.size": 9})

# --- figura annotata
fig, ax = plt.subplots(figsize=(10, 5.2))
ax.imshow(rgb.astype(np.uint8), extent=(0, W, H, 0))
t = np.linspace(0, math.pi, 200)
ax.plot(cx + R * np.cos(t), cy - R * np.sin(t), color="cyan", lw=0.8)
ax.axhline(y_asse, color="yellow", lw=0.8, ls="--")
ax.text(1180, y_asse - 6, "asse di simmetria", color="yellow", fontsize=7)
ax.plot(prua[:, 0], prua[:, 1], color="red", lw=1.0)
ax.plot(ss[:, 0], ss[:, 1], ".", color="orange", ms=2)
xx = np.array([x_sep, 760]); ax.plot(xx, a_ss1 * xx + b_ss1, color="orange", lw=1.0)
ax.plot(urt[:, 0], urt[:, 1], ".", color="magenta", ms=2)
xx = np.array([1060, 1240]); ax.plot(xx, a_ur1 * xx + b_ur1, color="magenta", lw=1.2)
xx = np.array([x_radice_ur, 1010]); ax.plot(xx, a_ur0 * xx + b_ur0, color="magenta", lw=0.8, ls=":")
xx = np.array([480, 1200]); ax.plot(xx, a_us * xx + b_us, color="lime", lw=0.9, ls="--")
ax.plot([cx, x_sep], [cy, y_sep], color="cyan", lw=0.8)
ax.plot([x_urto_asse, x_naso], [y_asse - 12, y_asse - 12], color="white", lw=1.2)
ax.annotate(f"Δ ≈ {Delta_mm:.1f} mm".replace(".", ","), (x_urto_asse, y_asse - 14), xytext=(10, 560), color="white",
            fontsize=8, arrowprops=dict(arrowstyle="->", color="white", lw=0.7))
for ang, lab in [(180 - theta_c, "3 (sonico, catena)"), (180 - theta_c_B, "3' (Busemann)")]:
    a = math.radians(ang)
    ax.plot(cx + R * math.cos(a), cy - R * math.sin(a), "o", color="yellow", ms=4)
    ax.plot([cx, cx + R * math.cos(a)], [cy, cy - R * math.sin(a)], color="yellow", lw=0.6)
ax.annotate("3", (cx + R * math.cos(math.radians(180 - theta_c)), cy - R * math.sin(math.radians(180 - theta_c))),
            xytext=(255, 455), color="yellow", fontsize=8, arrowprops=dict(arrowstyle="-", color="yellow", lw=0.5))
ax.annotate("3'", (cx + R * math.cos(math.radians(180 - theta_c_B)), cy - R * math.sin(math.radians(180 - theta_c_B))),
            xytext=(215, 488), color="yellow", fontsize=8, arrowprops=dict(arrowstyle="-", color="yellow", lw=0.5))
ax.text(cx + 30, cy - 18, f"θs ≈ {theta_s:.0f}°", color="cyan", fontsize=8)
etichette = [((60, 300), "1  (M = 10)", "white"), ((170, 330), "urto di prua", "red"),
             ((235, 420), "2 (subsonico)", "white"), ((560, 470), "corrente espansa", "white"),
             ((492, 520), "S", "white"), ((610, 575), "W", "white"),
             ((560, 640), "B (bolla di ricircolo)", "white"), ((1150, 600), "6", "white"),
             ((990, 520), "urto di separazione", "lime"), ((1000, 640), "urto di riattacco", "magenta"),
             ((640, 612), "strato di scorrimento", "orange"), ((800, 655), "collo", "white")]
for (x, y), s, c in etichette:
    ax.text(x, y, s, color=c, fontsize=7.5, bbox=dict(facecolor="black", alpha=0.35, lw=0, pad=1))
ax.text(560, 625, f"{ang_ss1:.0f}°", color="orange", fontsize=7)
ax.text(1180, 585, f"{sig_ax:.0f}° (formato)", color="magenta", fontsize=7)
ax.text(840, 625, f"{sig_ax_piede:.0f}° (piede)", color="magenta", fontsize=7)
ax.set_xlim(0, W); ax.set_ylim(H, 0)
ax.axis("off")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "annotata.png"), dpi=200)
plt.close(fig)

# --- schema geometrico (semicirconferenza e angoli, come nello schizzo d'esame)
fig, ax = plt.subplots(figsize=(7, 3.4))
t = np.linspace(0, math.pi, 300)
ax.plot(np.cos(t), np.sin(t), "k", lw=1.2)
ax.plot([-1.6, 2.6], [0, 0], "k", lw=0.8)
ax.plot(-1 - Delta_mm / R_mm, 0, "r|", ms=12)
ax.text(-1 - Delta_mm / R_mm - 0.12, 0.08, "urto", color="r", fontsize=8)
# punto sonico (catena e Busemann) e tangenti
for thc, c, lab in [(theta_c, "C0", "3"), (theta_c_B, "C2", "3'")]:
    a = math.radians(180 - thc)
    P = np.array([math.cos(a), math.sin(a)])
    ax.plot([0, P[0]], [0, P[1]], color=c, lw=0.8)
    tg = np.array([math.sin(a), -math.cos(a)])           # tangente
    ax.plot([P[0] - 0.6 * tg[0], P[0] + 0.5 * tg[0]], [P[1] - 0.6 * tg[1], P[1] + 0.5 * tg[1]], color=c, lw=0.8, ls="--")
    ax.text(P[0] - 0.15, P[1] + 0.05, lab, color=c)
ax.text(-0.45, 0.12, f"θc = {theta_c:.0f}° / {theta_c_B:.0f}°", fontsize=8)
# separazione
a = math.radians(theta_s)
Ps = np.array([math.cos(a), math.sin(a)])
ax.plot([0, Ps[0]], [0, Ps[1]], "C3", lw=0.8)
ax.text(0.22, 0.1, f"θs = {theta_s:.0f}°", fontsize=8, color="C3")
tg = np.array([math.sin(a), -math.cos(a)])
ax.plot([Ps[0], Ps[0] + 0.75 * tg[0]], [Ps[1], Ps[1] + 0.75 * tg[1]], "C3", lw=0.8, ls="--")
ax.text(Ps[0] + 0.05, Ps[1] - 0.35, f"parete {theta_wall_sep:.0f}°", fontsize=8, color="C3")
# strato di scorrimento e urto di riattacco (schematici, in R)
L = (Ps[1]) / math.tan(math.radians(Dth))
Q = Ps + L * np.array([1, -math.tan(math.radians(Dth))])
ax.plot([Ps[0], Q[0]], [Ps[1], Q[1]], color="orange", lw=1.2)
ax.text(Ps[0] + 0.35, Ps[1] - 0.05, f"Δθ = {Dth:.0f}°", fontsize=8, color="orange")
ax.plot([Q[0] - 0.1, Q[0] + 1.0], [0, 1.1 * math.tan(math.radians(sig_ax))], color="m", lw=1.0)
ax.text(Q[0] + 0.3, 0.2, f"urto, {sig_ax:.0f}° sull'asse\nσ = {sigma_R:.0f}° sulla corrente W", fontsize=8, color="m")
ax.set_aspect("equal"); ax.axis("off")
ax.set_xlim(-1.7, Q[0] + 1.2); ax.set_ylim(-0.15, 1.25)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "schema.pdf"))
plt.close(fig)

# --- piano p-theta: esatto, in scala logaritmica e lineare (ampio intervallo di pressioni)
for scala in ("log", "lin"):
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    pol1 = gd.polare_p_theta(st1, "C+", n=600)
    u = pol1["urto"]; deb = u["ramo"] == "debole"
    ax.plot(u["theta"][deb], u["p"][deb], "C0", lw=1.5, label="polare di 1 (M = 10), urto di prua")
    ax.plot(u["theta"][~deb], u["p"][~deb], "C0", lw=1.5, ls=":")
    ax.plot(thmax1, pol1["punti"]["theta_max"][1], "C0|", ms=10)
    # punto sonico della polare 1 (regola di Busemann)
    st_son = gd.urto_obliquo(st1, theta=th_son, famiglia="C+")
    ax.plot(th_son, st_son.valle.p, "s", color="C2", ms=5)
    ax.annotate("punto sonico della polare\n(regola di Busemann)", (th_son, st_son.valle.p), xytext=(5, 230) if scala == "log" else (8, 110), fontsize=7.5, color="C2",
                arrowprops=dict(arrowstyle="->", color="C2", lw=0.6))
    # curva isentropica della corrente di parete (p0 = p0', theta + nu = theta*)
    Ms = np.linspace(1.0, M_sep * 1.02, 600)
    ths = theta_star - gd.prandtl_meyer(Ms, aria)
    ps = p0p / gd.isentropiche(Ms, aria)["p0_p"]
    ax.plot(ths, ps, "C1--", lw=1.3, label="parete: onda semplice C+ (p0 = p0', θ + ν = θ*)")
    ax.plot([90, 90], [p0p * 0.95, p0p * 1.05], alpha=0)
    # polare dell'urto di riattacco da W
    polW = gd.polare_p_theta(stW, "C+", n=400)
    u = polW["urto"]; deb = u["ramo"] == "debole"
    ax.plot(u["theta"][deb], u["p"][deb], "C3", lw=1.3, label="polare di W (urto di riattacco)")
    ax.plot(u["theta"][~deb], u["p"][~deb], "C3", lw=1.0, ls=":")
    ax.axhline(pB, color="0.4", lw=0.8, ls="-.")
    if scala == "log":
        ax.text(-75, pB * 1.15, f"p bolla = p_W = {pB:.2f} p1", fontsize=7.5, color="0.3")
    ax.axhline(1.0, color="0.7", lw=0.6)
    for st, nm, off in [(st1, "1", (5, 5)), (st2, "2 (normale)", (-55, 6)), (st3, "3 (M = 1)", (6, 0)),
                        (st_sep, "S", (6, -2)), (stW, "W", (-14, -12)), (st6, "6", (6, 4))]:
        ax.plot(st.theta, st.p, "ko", ms=4)
        ax.annotate(nm, (st.theta, st.p), textcoords="offset points", xytext=off, fontsize=8)
    if scala == "log":
        ax.set_yscale("log")
        ax.set_xlim(-80, 65); ax.set_ylim(2e-3, 400)
        ax.set_ylabel("p / p1  (scala logaritmica)")
        ax.legend(fontsize=7.5, loc="lower right")
    else:
        ax.set_xlim(-60, 55); ax.set_ylim(0, 125)
        ax.set_ylabel("p / p1")
        ax.legend(fontsize=7.5, loc="center left")
    ax.set_xlabel("θ [°]")
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "p_theta.pdf" if scala == "log" else "p_theta_lin.pdf"))
    plt.close(fig)

# --- piano odografo (u*, v*) = V/a* (cos theta, sin theta): a* unico (h0 uniforme)
fig, ax = plt.subplots(figsize=(7.0, 5.6))
o1 = gd.polare_odografo(st1, "C+", n=600)
ax.plot(o1["urto"]["us"], o1["urto"]["vs"], "C0", lw=1.5, label="polare di 1 (strofoide)")
oW = gd.polare_odografo(stW, "C+", n=400)
ax.plot(oW["urto"]["us"], oW["urto"]["vs"], "C3", lw=1.3, label="polare di W")
Mst = np.sqrt(2.4 * Ms ** 2 / (2 + 0.4 * Ms ** 2))
ax.plot(Mst * np.cos(np.radians(ths)), Mst * np.sin(np.radians(ths)), "C1--", lw=1.3,
        label="epicicloide di parete (θ + ν = θ*)")
tt = np.linspace(-math.pi, math.pi, 400)
ax.plot(np.cos(tt), np.sin(tt), color="0.5", lw=0.8, ls="-.", label="M* = 1")
Mlim = math.sqrt(6.0)
ax.plot(Mlim * np.cos(tt), Mlim * np.sin(tt), color="0.7", lw=0.8, ls=":", label="M* limite = √6")
mark = [(st1, "1", (5, 5)), (st2, "2", (5, 5)), (st3, "3", (5, 3)), (st_sep, "S", (6, -6)),
        (stW, "W", (6, -8)), (st6, "6", (5, 5))]
for st, nm, off in mark:
    x, y = st.M_star * math.cos(math.radians(st.theta)), st.M_star * math.sin(math.radians(st.theta))
    ax.plot(x, y, "ko", ms=4)
    ax.annotate(nm, (x, y), textcoords="offset points", xytext=off, fontsize=8)
st3B = gd.State.da_totali(aria, 1.0, p0=p0p, a0=st1.a0, theta=th_son)
ax.plot(math.cos(math.radians(th_son)), math.sin(math.radians(th_son)), "s", color="C2", ms=5)
ax.annotate("3' (Busemann)", (math.cos(math.radians(th_son)), math.sin(math.radians(th_son))),
            textcoords="offset points", xytext=(-75, 8), fontsize=7.5, color="C2")
ax.set_aspect("equal")
ax.set_xlim(-0.2, 2.7); ax.set_ylim(-2.1, 1.3)
ax.set_xlabel("u* = u/a*"); ax.set_ylabel("v* = v/a*")
ax.grid(alpha=0.25)
ax.legend(fontsize=7.5, loc="lower left")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "odografo.pdf"))
plt.close(fig)
print("figure salvate in", FIG)
