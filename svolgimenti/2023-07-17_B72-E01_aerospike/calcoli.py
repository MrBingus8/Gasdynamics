"""Prova del 17/07/2023, esercizio 1 (B72-E01): ugello aerospike piano in schlieren, p0/pa = 25.

Riproduce tutte le misure sull'immagine, tutti i numeri e le figure dello svolgimento.
Uso:  python calcoli.py   (dalla cartella dello svolgimento)

Immagine: figure/originale.jpg è la copia digitale nitida della figura della traccia (la stessa
della foto figure/traccia_foto.jpeg, estratta da fonti/soluzioni_A/A_G_B72-E01.pdf): è già dritta
(le facce verticali dei carenaggi sono verticali), quindi non serve raddrizzarla.
Coordinate immagine: x verso destra, y verso il basso, in pixel (640 x 608).
Coordinate fisiche: X = (x - x_labbro)/h, Y = (y_asse - y)/h, con h = semidistanza fra i labbri;
angoli in gradi, antiorari dall'asse (θ < 0 = verso l'asse nella metà superiore).
"""
import os
import sys
import math

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, sobel
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

QUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(QUI, "..", ".."))
import gasdyn as gd  # noqa: E402

FIG = os.path.join(QUI, "figure")
aria = gd.Gas(1.4)                 # aria (galleria a freddo), γ assunto
deg, rad = math.degrees, math.radians


def mu(M):
    return gd.angolo_mach(M)


def nu(M):
    return gd.prandtl_meyer(M, aria)            # F §11


def M_da_nu(v):
    return gd.mach_da_pm(v, aria)


def p0_p(M):
    return gd.isentropiche(M, aria)["p0_p"]     # F §8


img = Image.open(os.path.join(FIG, "originale.jpg")).convert("L")
A = np.asarray(img).astype(float)
H, W = A.shape

# =============================================================================
# 1. MISURE SULL'IMMAGINE
# =============================================================================
print("=" * 78)
print("1. MISURE SULL'IMMAGINE (pixel; angoli antiorari dall'orizzontale)")
print("=" * 78)


def fit_retta(pts):
    """Retta ai minimi quadrati y = a + b x; angolo 'matematico' (y verso l'alto) e
    incertezza dell'angolo dall'errore standard della pendenza."""
    p = np.asarray(pts, float)
    x, y = p[:, 0], p[:, 1]
    b, a = np.polyfit(x, y, 1)
    res = y - (a + b * x)
    n = len(x)
    sb = math.sqrt(res.var(ddof=2) / ((x - x.mean()) ** 2).sum()) if n > 2 else float("nan")
    ang = -deg(math.atan(b))
    dang = deg(sb / (1 + b * b))
    return a, b, ang, dang


# --- punti letti sulle griglie ingrandite (lettura "a righello") ----------------
LABBRO_SUP = (180.0, 143.0)        # spigolo del carenaggio superiore
LABBRO_INF = (180.0, 440.0)
BASE = (469.0, 278.0, 297.0)       # x della base troncata, y bordo superiore e inferiore
Y_ASSE = 0.5 * (0.5 * (LABBRO_SUP[1] + LABBRO_INF[1]) + 0.5 * (BASE[1] + BASE[2]))
h = 0.5 * (LABBRO_INF[1] - LABBRO_SUP[1])
# pareti del canale in gola (metà superiore): carenaggio interno e spina
GOLA_CARENA = [(150, 113), (160, 122), (170, 132), (178, 140)]
GOLA_SPINA = [(150, 160), (160, 170), (170, 180), (180, 190), (190, 200)]
# contorno della spina superiore (linea sottile grigia della parete)
SPINA_SUP = [(170, 180), (190, 196), (210, 209), (242, 224), (270, 234), (300, 243),
             (330, 252), (360, 259), (400, 267), (424, 271), (469, 278)]

_, _, th_carena, _ = fit_retta(GOLA_CARENA)
_, _, th_spina_g, _ = fit_retta(GOLA_SPINA)
theta0_mis = 0.5 * (th_carena + th_spina_g)
# larghezza di gola: distanza fra le due pareti (parallele) misurata in verticale x cos45°
a_c, b_c, _, _ = fit_retta(GOLA_CARENA)
a_s, b_s, _, _ = fit_retta(GOLA_SPINA)
xg = 178.0
gap_v = (a_s + b_s * xg) - (a_c + b_c * xg)
h_gola = gap_v * math.cos(rad(abs(theta0_mis)))
print(f"asse y = {Y_ASSE:.1f} px, h (semidistanza fra i labbri) = {h:.1f} px")
print(f"gola: parete del carenaggio {th_carena:.1f}°, parete della spina {th_spina_g:.1f}° "
      f"-> θ0 ≈ {theta0_mis:.1f}° ; larghezza di gola ≈ {h_gola:.0f} px = {h_gola/h:.2f} h")

# pendenza della spina a tratti
print("spina superiore, inclinazione a tratti:")
for (x1, y1), (x2, y2) in zip(SPINA_SUP[:-1], SPINA_SUP[1:]):
    print(f"   x {x1:3.0f}-{x2:3.0f}: {-deg(math.atan2(y2 - y1, x2 - x1)):6.1f}°")
_, _, th_spina_fine, dth_sf = fit_retta(SPINA_SUP[-3:])
print(f"   ultimo tratto (x 400-469): {th_spina_fine:.1f}°")

# --- rilevamenti automatici (gradienti di luminosità lungo le colonne) ---------
B = gaussian_filter(A, 2.0)
GY = np.gradient(B, axis=0)
B15 = gaussian_filter(A, 1.5)
GY15 = np.gradient(B15, axis=0)


def scansione(xs, yc, w, campo, modo):
    out = []
    for x in xs:
        y0 = int(yc(x) - w)
        col = campo[y0:int(yc(x) + w), x - 1:x + 2].mean(1)
        i = np.argmin(col) if modo == "min" else np.argmax(col)
        out.append((x, y0 + i))
    return out


# contorno del getto: transizione chiaro (ambiente) -> scuro (strato di taglio) sopra,
# chiaro (getto) -> scuro (ambiente) sotto
bordo_sup_lab = scansione(range(190, 251, 10), lambda x: 140 + 0.13 * (x - 180), 12, -GY, "max")
bordo_sup = scansione(range(250, 531, 20), lambda x: 140 + 0.13 * (x - 180), 12, -GY, "max")
bordo_inf = scansione(range(250, 531, 20), lambda x: 440 - 0.2 * (x - 180), 25, -GY, "max")
_, _, th_b_lab, dth_b_lab = fit_retta(bordo_sup_lab)
_, _, th_b_sup, dth_b_sup = fit_retta(bordo_sup)
_, _, th_b_inf, dth_b_inf = fit_retta(bordo_inf)
print(f"contorno 1 sup. vicino al labbro (x 190-250): {th_b_lab:.1f}° ± {dth_b_lab:.1f}")
print(f"contorno 1 sup. a valle (x 250-530): {th_b_sup:.1f}° ± {dth_b_sup:.1f}")
print(f"contorno 1 inf. a valle (x 250-530): {th_b_inf:.1f}° ± {dth_b_inf:.1f}")
theta_C_mis = 0.5 * (abs(th_b_sup) + abs(th_b_inf))      # verso l'asse, in modulo

# linea 3: superiore = trattini scuri (minimi di luminosità); inferiore = bordo scuro->chiaro
tre_sup = scansione(range(240, 411, 10), lambda x: 143 + 0.5 * (x - 180), 6, B15, "min")
tre_inf = scansione(range(200, 401, 10), lambda x: 440 - 0.52 * (x - 180), 8, GY15, "max")
a3s, b3s, th3_sup, dth3_sup = fit_retta(tre_sup)
a3i, b3i, th3_inf, dth3_inf = fit_retta(tre_inf)
print(f"linea 3 sup.: {th3_sup:.1f}° ± {dth3_sup:.1f}, passa a x=180 per y={a3s + 180*b3s:.0f} "
      f"(labbro {LABBRO_SUP[1]:.0f})")
print(f"linea 3 inf.: {th3_inf:.1f}° ± {dth3_inf:.1f}, passa a x=180 per y={a3i + 180*b3i:.0f} "
      f"(labbro {LABBRO_INF[1]:.0f})")
phi3_mis = 0.5 * (abs(th3_sup) + abs(th3_inf))

# punto R in cui 3 incontra la spina superiore (ultimo tratto della spina)
a_sp, b_sp, _, _ = fit_retta(SPINA_SUP[-3:])
xR = (a_sp - a3s) / (b3s - b_sp)
yR = a3s + b3s * xR
# metà inferiore: la linea 3 finisce sulla spina (y della spina inferiore ≈ 2 Y_asse - y_sup)
a_spi, b_spi = 2 * Y_ASSE - a_sp, -b_sp
xRi = (a_spi - a3i) / (b3i - b_spi)
print(f"R (3 sulla spina): sup. x={xR:.0f}, y={yR:.0f};  inf. x={xRi:.0f}  "
      f"-> X_R ≈ {(0.5*(xR+xRi)-LABBRO_SUP[0])/h:.2f} h ; base a X = {(BASE[0]-LABBRO_SUP[0])/h:.2f} h, "
      f"altezza base {(BASE[2]-BASE[1])/h:.2f} h")

# riflesso 4 (metà inferiore): linea scura che lascia la spina poco a monte della base
quattro_inf = scansione(range(450, 541, 10), lambda x: 303 + 0.18 * (x - 420), 7, gaussian_filter(A, 1.2), "min")
_, _, th4_inf, dth4_inf = fit_retta(quattro_inf)
print(f"riflesso 4 inf.: {th4_inf:.1f}° ± {dth4_inf:.1f}")

# --- orientazione dominante delle linee in finestre (tensore di struttura) -----
Gs = gaussian_filter(A, 1.2)
gx, gy = sobel(Gs, 1), sobel(Gs, 0)


def orientazione(x0, y0, x1, y1):
    Jxx = (gx[y0:y1, x0:x1] ** 2).sum()
    Jyy = (gy[y0:y1, x0:x1] ** 2).sum()
    Jxy = (gx[y0:y1, x0:x1] * gy[y0:y1, x0:x1]).sum()
    ang = 0.5 * deg(math.atan2(2 * Jxy, Jxx - Jyy))      # direzione del gradiente (y in giù)
    coer = math.hypot(Jxx - Jyy, 2 * Jxy) / (Jxx + Jyy)
    m = -(ang + 90.0)                                     # direzione della linea, y in su
    m = (m + 90.0) % 180.0 - 90.0
    return m, coer


finestre = {
    "C+ zona C sup. (fra 3 e 1)": [(340, 180, 400, 210), (380, 190, 440, 225)],
    "C- zona C inf. (fra 3 e 1)": [(330, 365, 370, 400), (380, 350, 430, 385)],
    "urto di ricompressione sup.": [(500, 250, 560, 270), (540, 240, 600, 262)],
    "urto di ricompressione inf.": [(500, 300, 560, 320), (540, 305, 600, 330)],
    "scia sull'asse": [(500, 282, 600, 292)],
}
ORI = {}
for nome, ww in finestre.items():
    vals = [orientazione(*w) for w in ww]
    ORI[nome] = float(np.mean([v[0] for v in vals]))
    print(f"{nome:30s}: " + ", ".join(f"{v[0]:6.1f}° (coer. {v[1]:.2f})" for v in vals))
phi_Cp = 0.5 * (abs(ORI["C+ zona C sup. (fra 3 e 1)"]) + abs(ORI["C- zona C inf. (fra 3 e 1)"]))
phi_ric = 0.5 * (abs(ORI["urto di ricompressione sup."]) + abs(ORI["urto di ricompressione inf."]))

# raggio labbro -> estremità della spina
phi_raggio = deg(math.atan2(BASE[1] - LABBRO_SUP[1], BASE[0] - LABBRO_SUP[0]))
print(f"raggio dal labbro allo spigolo della base: {phi_raggio:.1f}° sotto l'orizzontale")

# valori arrotondati usati nello svolgimento (misure "a goniometro")
THETA0 = -45.0          # inclinazione della corrente in gola
THETA_C = -8.5          # contorno del getto a valle (metà superiore)
THETA_LAB = -4.0        # contorno del getto vicino al labbro
PHI3 = -28.0            # linea 3 (metà superiore)
PHI_CP = 12.0           # linee C+ nella zona C (metà superiore)
PHI_RIC = 18.0          # urti di ricompressione nella scia
THETA_SPINA = -9.0      # ultimo tratto della spina
PHI_RAGGIO = -25.0      # raggio labbro -> estremità della spina
PHI4 = 12.0             # riflesso 4 (in modulo)
print(f"\nvalori adottati: θ0={THETA0}, θ_C={THETA_C}, θ(labbro)={THETA_LAB}, linea 3={PHI3}, "
      f"C+={PHI_CP}, ricompressione=±{PHI_RIC}, spina={THETA_SPINA}, raggio={PHI_RAGGIO}")
print(f"(medie misurate: θ0={theta0_mis:.1f}, |θ_C|={theta_C_mis:.1f}, |3|={phi3_mis:.1f}, "
      f"|C±|={phi_Cp:.1f}, |ricompr.|={phi_ric:.1f}, spina={th_spina_fine:.1f}, raggio={phi_raggio:.1f})")

# =============================================================================
# 2. CALCOLO
# =============================================================================
print("\n" + "=" * 78)
print("2. CALCOLO (metà superiore; la inferiore è simmetrica)")
print("=" * 78)
NPR = 25.0

# Gola (G): sonica al labbro.  p*/p0 = (2/(γ+1))^(γ/(γ-1))  (F §8)
pstar_p0 = 1.0 / p0_p(1.0)
pG = NPR * pstar_p0
print(f"G: M=1, θ0={THETA0}°, p*/p0={pstar_p0:.4f}, p_G/pa={pG:.2f}  (sottoespansione)")

# Stato adattato sul contorno: p = pa  ->  M_a da p0/pa (F §8), ν, μ (F §10-11)
Ma = gd.mach_da_rapporto(NPR, aria, "p0/p")
nua, mua = nu(Ma), mu(Ma)
print(f"contorno (p=pa): M_a={Ma:.3f}, ν_a={nua:.2f}°, μ_a={mua:.2f}°, T/T0={1/(1+0.2*Ma**2):.3f}")

# A: ventaglio centrato al labbro (onde C-): lungo le C+ θ-ν = cost  ->  θ_A = θ0 + ν_a
thA = THETA0 + nua
coda = thA - mua
print(f"A: θ_A = θ0 + ν_a = {thA:.2f}°  ; coda del ventaglio (C-) a θ_A-μ_a = {coda:.1f}°")
print(f"   ugello di progetto (spina ideale): θ0 = -ν(M_a) = {-nua:.1f}° ≈ misurato {theta0_mis:.1f}°")
L_ideale = 1.0 / math.tan(rad(mua))
print(f"   la coda del ventaglio incontrerebbe l'asse a X = h/tan μ_a = {L_ideale:.2f} h "
      f"(spina reale fino a {(BASE[0]-LABBRO_SUP[0])/h:.2f} h)")

# Contorno che si piega da θ_A a θ_C: espansioni C+ (dalla spina) riflesse come compressioni C-.
# All'arrivo di un'espansione C+ di ampiezza Δν: θ+ν = cost lungo C- -> θ diminuisce di Δν;
# la riflessa C- riporta p = pa: θ-ν = cost lungo C+ -> θ diminuisce di altri Δν  => Δθ = 2Δν
dnu = (thA - THETA_C) / 2.0
nuB = nua + dnu
MB = M_da_nu(nuB)
thB = thA - dnu
pB = NPR / p0_p(MB)
print(f"Δν = (θ_A - θ_C)/2 = {dnu:.2f}°")
print(f"B (prima di 3): ν_B={nuB:.2f}°, M_B={MB:.3f}, θ_B={thB:.2f}°, p_B/pa={pB:.3f}")

# Urto 3 (famiglia C-) da B, deviazione Δν (stessa deviazione delle compressioni che lo formano)
sB = gd.State.da_M(aria, MB, theta=thB, p=pB, a=1.0)
u3 = gd.urto_obliquo(sB, theta=dnu, famiglia="C-")
sC = u3.valle
print(f"urto 3: σ={u3.sigma:.2f}°, fronte a θ_B-σ = {thB - u3.sigma:.1f}° (misurato {PHI3}°), "
      f"Mn1={u3.Mn1:.3f}, p_C/p_B={u3.p2_p1:.3f}, p_C/pa={sC.p:.3f}, M_C={sC.M:.3f}, θ_C={sC.theta:.2f}°, "
      f"p0C/p0B={u3.p02_p01:.4f}")
print(f"   confronto: coda del ventaglio (onda di Mach) a {coda:.1f}°; inclinazione di 3 rispetto "
      f"alla corrente B: {abs(PHI3 - thB):.1f}° (μ_B={mu(MB):.1f}°)")
print(f"   compressione isentropica equivalente da B fino a pa: θ = θ_B - (ν_B - ν_a) = {thB - dnu:.2f}°")

# Mach in C dalle linee C+ (θ+μ) e dal contorno (θ_C): μ = φ_C+ - θ_C
mu_C_mis = PHI_CP - THETA_C
M_C_mis = 1.0 / math.sin(rad(mu_C_mis))
print(f"C: μ = φ(C+) - θ_C = {mu_C_mis:.1f}° -> M ≈ {M_C_mis:.2f} (M_a = {Ma:.2f}); "
      f"C+ previste a θ_C+μ_a = {THETA_C + mua:.1f}°")

# La spina non è la linea di corrente dell'onda semplice: sul raggio dal labbro allo spigolo
# della base, nel ventaglio θ - μ = φ_raggio con θ = θ0 + ν(M)  ->  ν(M) - μ(M) = φ_raggio - θ0
f = lambda M: nu(M) - mu(M) - (PHI_RAGGIO - THETA0)
M_raggio = brentq(f, 1.01, 5.0)
th_ideale = THETA0 + nu(M_raggio)
print(f"raggio labbro-spigolo {PHI_RAGGIO}°: onda semplice -> M={M_raggio:.2f}, θ_ideale={th_ideale:.1f}° "
      f"contro spina {THETA_SPINA}° : differenza {th_ideale - THETA_SPINA:.1f}°")

# Riflessione 4 sulla spina: regolare se la deviazione è < θ_D(M)
thD = gd.theta_distacco(MB, aria)
thmax_B = gd.theta_max(MB, aria)[0]
print(f"riflessione 4: deviazione ≈ {dnu:.1f}° contro θ_D(M_B)={thD:.1f}° (θmax={thmax_B:.1f}°) -> regolare")
# ordine di grandezza: riflessione regolare con lo stato B diretto come la parete
sBw = gd.State.da_M(aria, MB, theta=THETA_SPINA, p=pB, a=1.0)
i4 = gd.urto_obliquo(sBw, theta=dnu, famiglia="C-")
r4 = gd.urto_obliquo(i4.valle, theta=dnu, famiglia="C+")
print(f"   (stima) incidente σ={i4.sigma:.1f}°, riflesso σ={r4.sigma:.1f}°, fronte del riflesso "
      f"{i4.valle.theta + r4.sigma:.1f}° (misurato ±{PHI4}°), p_D/pa={r4.valle.p:.2f}")

# Scia (2): getto C (θ_C) riportato in direzione assiale da un urto C+ (stima a urto unico)
sCa = gd.State.da_M(aria, Ma, theta=THETA_C, p=1.0, a=1.0)
ur = gd.urto_obliquo(sCa, theta=abs(THETA_C), famiglia="C+")
print(f"ricompressione (stima): σ={ur.sigma:.1f}°, fronte {THETA_C + ur.sigma:.1f}° (misurato ±{PHI_RIC}°), "
      f"p_F/pa={ur.p2_p1:.2f}, M_F={ur.valle.M:.2f}")
mu_ric = mu(Ma)
print(f"   un'onda di Mach della stessa corrente sarebbe a {THETA_C + mu_ric:.1f}°")

# Riepilogo pressioni
print("\nRiepilogo: p/pa  G={:.2f}  A=1  B={:.2f}  C={:.2f}  F={:.2f}".format(pG, pB, sC.p, ur.p2_p1))

# =============================================================================
# 3. FIGURE
# =============================================================================
plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans"})

# --- 3a. Immagine annotata e dettagli -----------------------------------------
RGB = np.asarray(Image.open(os.path.join(FIG, "originale.jpg")).convert("RGB"))
PARETE_SUP = np.array(GOLA_SPINA[:-1] + SPINA_SUP[1:], float)


def raggio_fino_alla_spina(ang_deg):
    """lunghezza del raggio dal labbro superiore, con angolo ang_deg, fino alla parete della spina."""
    d = np.array([math.cos(rad(ang_deg)), -math.sin(rad(ang_deg))])
    P0 = np.array(LABBRO_SUP)
    for q1, q2 in zip(PARETE_SUP[:-1], PARETE_SUP[1:]):
        e = q2 - q1
        M_ = np.array([d, -e]).T
        if abs(np.linalg.det(M_)) < 1e-9:
            continue
        t, u = np.linalg.solve(M_, q1 - P0)
        if t > 0 and 0 <= u <= 1:
            return t
    return 60.0


def seg(p, ang_deg, L, **kw):
    """segmento da p con angolo 'matematico' (y in su) e lunghezza L px (sugli assi correnti)."""
    ax = plt.gca()
    x0, y0 = p
    x1 = x0 + L * math.cos(rad(ang_deg))
    y1 = y0 - L * math.sin(rad(ang_deg))
    ax.plot([x0, x1], [y0, y1], **kw)
    return x1, y1




def disegna_annotazioni(ax, fs=1.0):
    """strutture, misure e stati sovrapposti all'immagine (assi in pixel)."""
    plt.sca(ax)
    ax.imshow(RGB)
    ax.axis("off")
    for sgn, lab in ((1, LABBRO_SUP), (-1, LABBRO_INF)):
        # gola
        ag = rad(THETA0 * sgn - 90 * sgn)
        ax.plot([lab[0], lab[0] + h_gola * math.cos(ag)], [lab[1], lab[1] - h_gola * math.sin(ag)],
                color="cyan", lw=2.2)
        # ventaglio al labbro: alcuni raggi (tratteggio giallo) fra la testa e la coda ideale
        for M in (1.15, 1.4, 1.7, 2.0, 2.3):
            th = THETA0 + nu(M)
            seg(lab, sgn * (th - mu(M)), raggio_fino_alla_spina(th - mu(M)), color="yellow", lw=0.9,
                ls="--", alpha=0.95)
        # coda ideale (onda di Mach θ_A - μ_a): non coincide con 3
        seg(lab, sgn * coda, 330, color="yellow", lw=1.2, ls=":", alpha=0.95)
        # contorno del getto 1
        seg(lab, sgn * THETA_LAB, 70, color="lime", lw=2.0)
        x1, y1 = lab[0] + 70 * math.cos(rad(THETA_LAB)), lab[1] - sgn * 70 * math.sin(rad(THETA_LAB))
        seg((x1, y1), sgn * THETA_C, 400, color="lime", lw=2.0)
    # linee 3 (misurate)
    xx = np.array([180, xR])
    ax.plot(xx, a3s + b3s * xx, color="red", lw=2.0)
    xx = np.array([180, xRi])
    ax.plot(xx, a3i + b3i * xx, color="red", lw=2.0)
    # riflessioni 4 (direzione stimata del riflesso)
    for sgn, (x0, y0), p4 in ((1, (xR, yR), (xR, yR)), (-1, (xRi, a3i + b3i * xRi), quattro_inf[0])):
        ax.plot([x0], [y0], "o", ms=5, mfc="none", mec="magenta", mew=1.6)
        seg(p4, sgn * (i4.valle.theta + r4.sigma), 95, color="magenta", lw=2.0)
    # urti di ricompressione della scia e scia
    for sgn in (1, -1):
        seg((BASE[0] + 4, Y_ASSE - sgn * 4), sgn * PHI_RIC, 180, color="orange", lw=2.0)
    ax.plot([BASE[0], W - 5], [Y_ASSE, Y_ASSE], color="white", lw=1.0, ls=":")
    # linee C+ / C- di esempio nella zona C
    for (x0, y0) in ((330, 214), (365, 228)):
        seg((x0, y0), PHI_CP, 90, color="deepskyblue", lw=1.2)
    for (x0, y0) in ((330, 370), (365, 357)):
        seg((x0, y0), -PHI_CP, 90, color="deepskyblue", lw=1.2)

    kws = dict(fontsize=9 * fs, fontweight="bold", color="black", ha="center", va="center", clip_on=True,
               bbox=dict(boxstyle="round,pad=0.15", fc="khaki", ec="none", alpha=0.9))
    for t, (x, y) in {"G": (160, 165), "A": (232, 152), "B": (300, 222), "C": (420, 200), "D": (452, 262),
                      "W": (490, 288), "F": (585, 286), "a": (470, 95)}.items():
        ax.text(x, y, t, **kws)


fig, ax = plt.subplots(figsize=(7.2, 6.84), dpi=200)
disegna_annotazioni(ax)
ax.set_xlim(0, W)
ax.set_ylim(H, 0)
ax.text(8, 598, "giallo: ventaglio al labbro (punteggiato: coda prevista)   verde: contorno 1   rosso: urto 3\n"
        "magenta: riflessione 4   arancio: urti di ricompressione   azzurro: onde di Mach",
        fontsize=6.5, color="white", va="bottom",
        bbox=dict(boxstyle="round,pad=0.25", fc="black", ec="none", alpha=0.65))
fig.savefig(os.path.join(FIG, "annotata.png"), bbox_inches="tight", pad_inches=0.02)
plt.close(fig)

# dettagli (close-up) da affiancare al testo: labbro e inizio di 3; zona C; estremità della spina e scia
DETTAGLI = {"dettaglio_labbro": (135, 335, 110, 260), "dettaglio_zonaC": (300, 470, 150, 290),
            "dettaglio_base": (380, 620, 215, 360)}
for nome, (x0, x1, y0, y1) in DETTAGLI.items():
    fig, ax = plt.subplots(figsize=(4.0, 4.0 * (y1 - y0) / (x1 - x0)), dpi=250)
    disegna_annotazioni(ax, fs=1.25)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    fig.savefig(os.path.join(FIG, nome + ".png"), bbox_inches="tight", pad_inches=0.01)
    plt.close(fig)

# --- 3b. Schema della metà superiore in unità di h ------------------------------
def XY(p):
    return ((p[0] - LABBRO_SUP[0]) / h, (Y_ASSE - p[1]) / h)


fig, ax = plt.subplots(figsize=(7.4, 4.3))
sp = np.array([XY(p) for p in PARETE_SUP])
# corpo della spina e carenaggio (ombreggiati)
ax.fill_between(np.r_[sp[:, 0], sp[-1, 0]], np.r_[sp[:, 1], 0], 0, color="0.75")
car = np.array([XY(p) for p in GOLA_CARENA])
ax.plot(car[:, 0], car[:, 1], "k", lw=1.5)
ax.plot([0, 0], [car[-1, 1], 1.35], "k", lw=1.5)
ax.fill_betweenx([car[-1, 1], 1.35], -0.25, 0, color="0.75")
ax.plot(sp[:, 0], sp[:, 1], "k", lw=1.5)
ax.plot([sp[-1, 0], sp[-1, 0]], [sp[-1, 1], 0], "k", lw=1.5)
L0 = XY(LABBRO_SUP)
# gola
gx_ = L0[0] + (h_gola / h) * math.cos(rad(THETA0 - 90))
gy_ = L0[1] + (h_gola / h) * math.sin(rad(THETA0 - 90))
ax.plot([L0[0], gx_], [L0[1], gy_], color="c", lw=2)
ax.text(gx_ - 0.13, gy_ + 0.02, "G", color="c", fontsize=10, fontweight="bold")
# ventaglio (raggi) e coda prevista
for M in (1.3, 1.7, 2.1, 2.45):
    th = THETA0 + nu(M)
    r = raggio_fino_alla_spina(th - mu(M)) / h
    ax.plot([L0[0], L0[0] + r * math.cos(rad(th - mu(M)))], [L0[1], L0[1] + r * math.sin(rad(th - mu(M)))],
            color="goldenrod", lw=0.9, ls="--")
ax.plot([L0[0], L0[0] + 2.4 * math.cos(rad(coda))], [L0[1], L0[1] + 2.4 * math.sin(rad(coda))],
        color="goldenrod", lw=1.0, ls=":")
ax.text(1.95, 0.32, "coda prevista\ndel ventaglio", color="goldenrod", fontsize=8)
# contorno ideale (tratteggio) e reale
ax.plot([L0[0], 2.9], [L0[1], L0[1] + 2.9 * math.tan(rad(thA))], color="g", ls="--", lw=1)
xb = 0.45
yb = L0[1] + xb * math.tan(rad(THETA_LAB))
ax.plot([L0[0], xb, 2.9], [L0[1], yb, yb + (2.9 - xb) * math.tan(rad(THETA_C))], color="g", lw=2)
ax.text(2.35, 1.06, "contorno ideale ($\\theta_A$)", color="g", fontsize=8)
ax.text(1.0, 1.02, "1: contorno del getto, $p=p_a$", color="g", fontsize=9)
# espansioni C+ dalla spina (blu) e compressioni C- riflesse che convergono in 3
# espansioni C+ dal lato della spina (blu) che raggiungono il contorno e si riflettono come
# compressioni C- (rosso tratteggiato) che convergono nell'urto 3 (schema qualitativo)
X3 = (xR - LABBRO_SUP[0]) / h
Y3 = (Y_ASSE - yR) / h
spp = np.array([XY(p) for p in PARETE_SUP])


def y_contorno(x):
    return L0[1] + (x * math.tan(rad(THETA_LAB)) if x <= xb else
                    xb * math.tan(rad(THETA_LAB)) + (x - xb) * math.tan(rad(THETA_C)))


for xs0, ang, x3 in ((-0.12, 38.0, 0.75), (0.0, 26.0, 1.0), (0.12, 17.0, 1.25)):
    ys0 = np.interp(xs0, spp[:, 0], spp[:, 1])
    xe = brentq(lambda x: ys0 + (x - xs0) * math.tan(rad(ang)) - y_contorno(x), xs0, 3.0)
    ye = y_contorno(xe)
    ax.annotate("", xy=(xe, ye), xytext=(xs0, ys0),
                arrowprops=dict(arrowstyle="->", color="royalblue", lw=1.0))
    ax.plot([xe, x3], [ye, L0[1] + (Y3 - L0[1]) * (x3 - L0[0]) / (X3 - L0[0])],
            color="r", lw=0.7, ls="--")
X3 = (xR - LABBRO_SUP[0]) / h
Y3 = (Y_ASSE - yR) / h
ax.plot([L0[0], X3], [L0[1], Y3], color="r", lw=2.2)
ax.text(0.75, 0.62, "3", color="r", fontsize=12, fontweight="bold")
# riflessione 4
ang4 = rad(i4.valle.theta + r4.sigma)
ax.plot([X3, X3 + 0.6 * math.cos(ang4)], [Y3, Y3 + 0.6 * math.sin(ang4)], color="m", lw=2)
ax.text(X3 + 0.05, Y3 + 0.2, "4", color="m", fontsize=12, fontweight="bold")
# base, scia, ricompressione
xbse = sp[-1, 0]
ax.add_patch(plt.Polygon([[xbse, 0], [xbse, sp[-1, 1]], [xbse + 0.22, 0]], closed=True,
                         fc="lightgray", ec="k", ls=":", lw=0.8))
ax.text(xbse + 0.03, 0.025, "W", fontsize=9)
ax.plot([xbse + 0.22, 3.0], [0, 0], "k:", lw=1)
ax.plot([xbse + 0.1, 3.0], [0.0, (3.0 - xbse - 0.1) * math.tan(rad(PHI_RIC))], color="darkorange", lw=2)
ax.text(2.7, 0.13, "2", color="darkorange", fontsize=12, fontweight="bold")
# stati
for t, (x, y) in {"A": (0.35, 0.93), "B": (0.75, 0.42), "C": (1.75, 0.78), "D": (1.95, 0.28),
                  "F": (2.75, 0.03), "a": (1.3, 1.25)}.items():
    ax.text(x, y, t, fontsize=10, fontweight="bold", bbox=dict(boxstyle="round,pad=0.12", fc="khaki", ec="none"))
ax.axhline(0, color="k", lw=0.6)
ax.set_xlim(-0.3, 3.0)
ax.set_ylim(-0.05, 1.35)
ax.set_aspect("equal")
ax.set_xlabel("$x/h$")
ax.set_ylabel("$y/h$")
ax.set_title("Metà superiore (la inferiore è simmetrica)", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "schema.pdf"))
plt.close(fig)

# --- 3c. Piano p-θ (log e lineare) ------------------------------------------------
Ms = np.linspace(1.0, 3.1, 300)
th_pm = THETA0 + np.array([nu(M) for M in Ms])
p_pm = NPR / np.array([p0_p(M) for M in Ms])
# A -> B: espansione C+ (θ + ν = cost)
nus = np.linspace(nua, nuB, 50)
th_AB = thA - (nus - nua)
p_AB = NPR / np.array([p0_p(M_da_nu(v)) for v in nus])
# polare dello stato B (C-) fino a C; polare di C (C+) verso θ = 0
polB = gd.polare_p_theta(sB, famiglia="C-")
polC = gd.polare_p_theta(sCa, famiglia="C+")

fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.3))
for ax, scala in zip(axs, ("log", "lin")):
    ax.plot(th_pm, p_pm, color="saddlebrown", lw=1.6, label="ventaglio al labbro ($\\theta-\\nu=\\theta_0$)")
    ax.plot(th_AB, p_AB, color="royalblue", lw=1.6, label="espansioni dalla spina ($\\theta+\\nu$ cost.)")
    u = polB["urto"]
    sel = np.asarray(u["ramo"]) == "debole"
    ax.plot(np.asarray(u["theta"])[sel], np.asarray(u["p"])[sel], color="r", lw=1.2,
            label="polare $C^-$ di B (urto 3)")
    u = polC["urto"]
    sel = np.asarray(u["ramo"]) == "debole"
    ax.plot(np.asarray(u["theta"])[sel], np.asarray(u["p"])[sel], color="darkorange", lw=1.2,
            label="polare $C^+$ di C (ricompressione)")
    ax.axhline(1.0, color="g", lw=0.8, ls="--")
    pts = {"G": (THETA0, pG), "A": (thA, 1.0), "B": (thB, pB), "C": (sC.theta, sC.p), "F": (0.0, ur.p2_p1)}
    for t, (x, y) in pts.items():
        ax.plot(x, y, "ko", ms=4)
        ax.annotate(t, (x, y), textcoords="offset points", xytext=(5, 4), fontsize=10, fontweight="bold")
    ax.set_xlabel("$\\theta$ [°]")
    ax.set_ylabel("$p/p_a$")
    ax.grid(alpha=0.3)
    if scala == "log":
        ax.set_yscale("log")
        ax.set_xlim(-50, 12)
        ax.set_ylim(0.5, 20)
        ax.set_title("scala logaritmica", fontsize=10)
        ax.text(-20, 1.06, "$p=p_a$ (contorno del getto)", color="g", fontsize=8)
    else:
        ax.set_xlim(-50, 12)
        ax.set_ylim(0, 14)
        ax.set_title("scala lineare", fontsize=10)
axs[0].legend(fontsize=7.5, loc="upper right")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "p_theta.pdf"))
plt.close(fig)

# dettaglio lineare vicino a p_a
fig, ax = plt.subplots(figsize=(5.2, 3.6))
ax.plot(th_pm, p_pm, color="saddlebrown", lw=1.6)
ax.plot(th_AB, p_AB, color="royalblue", lw=1.6)
u = polB["urto"]
sel = np.asarray(u["ramo"]) == "debole"
ax.plot(np.asarray(u["theta"])[sel], np.asarray(u["p"])[sel], color="r", lw=1.2)
u = polC["urto"]
sel = np.asarray(u["ramo"]) == "debole"
ax.plot(np.asarray(u["theta"])[sel], np.asarray(u["p"])[sel], color="darkorange", lw=1.2)
ax.axhline(1.0, color="g", lw=0.8, ls="--")
for t, (x, y) in {"A": (thA, 1.0), "B": (thB, pB), "C": (sC.theta, sC.p), "F": (0.0, ur.p2_p1)}.items():
    ax.plot(x, y, "ko", ms=4)
    ax.annotate(t, (x, y), textcoords="offset points", xytext=(5, 4), fontsize=10, fontweight="bold")
ax.set_xlim(-14, 4)
ax.set_ylim(0.5, 2.2)
ax.set_xlabel("$\\theta$ [°]")
ax.set_ylabel("$p/p_a$")
ax.grid(alpha=0.3)
ax.set_title("dettaglio vicino a $p_a$ (scala lineare)", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "p_theta_dettaglio.pdf"))
plt.close(fig)
print("\nfigure salvate in", FIG)

# =============================================================================
# 4. VERIFICHE
# =============================================================================
print("\n" + "=" * 78)
print("4. VERIFICHE")
print("=" * 78)
tmB, smB = gd.theta_max(MB, aria)
tma, sma = gd.theta_max(Ma, aria)
print(f"urto 3: deviazione {dnu:.1f}° < θmax(M_B)={tmB:.1f}°, σ={u3.sigma:.1f}° < σ(θmax)={smB:.1f}° -> ramo debole; "
      f"Mn1={u3.Mn1:.2f} (>1, < 1,2: quasi isentropico), Δp0/p0={1-u3.p02_p01:.4f}")
print(f"ricompressione: deviazione {abs(THETA_C):.1f}° < θmax(M_a)={tma:.1f}°, σ={ur.sigma:.1f}° < {sma:.1f}°; "
      f"Mn1={ur.Mn1:.2f}, p02/p01={ur.p02_p01:.3f}, M_F={ur.valle.M:.2f} > 1")
print(f"riflesso 4: Mn1={r4.Mn1:.2f}, M_D={r4.valle.M:.2f} > 1")
print(f"regioni del getto supersoniche: M da 1 (G) a {MB:.2f} (B); M_C={sC.M:.2f}, M_D={r4.valle.M:.2f}, M_F={ur.valle.M:.2f}")
print(f"p_C/pa dopo l'urto 3 = {sC.p:.3f} (condizione di contorno p = pa)")
