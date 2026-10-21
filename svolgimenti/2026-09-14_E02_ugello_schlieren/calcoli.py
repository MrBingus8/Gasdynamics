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

# --- stazioni per la legge delle aree (tratto attaccato, fra gola e W)
X_staz = [0.0, 0.05, 0.10, 0.15, 0.20]
staz = []
for Xs in X_staz:
    x = x_gola + Xs * h_star
    if Xs == 0.0:
        h = h_star
        th = 0.0
    else:
        h = np.polyval(c_inf, x) - np.polyval(c_sup, x)
        th = 0.5 * (theta_sup(x) - theta_inf(x))
    staz.append((Xs, h / h_star, th))

# =============================================================================
# 3. CALCOLO
# =============================================================================
print("\n--- Legge delle aree (valori medi sulla sezione) e Prandtl-Meyer a parete")
tab = []
for Xs, AA, th in staz:
    M1d = 1.0 if AA <= 1.0 + 1e-9 else gd.mach_da_area(AA, aria, "supersonico")
    Mw = 1.0 if th <= 0 else M_da_nu(th)
    tab.append((Xs, AA, M1d, th, Mw))
    print(f"X = {Xs:.2f}: A/A* = {AA:.3f}  M_1D = {M1d:.2f}  theta_w = {th:.1f}  M_parete = {Mw:.2f}")

# stato a parete in W (media dei due lati) e confronto con l'angolo dell'onda
th_W = 0.5 * (th_Ws - th_Wi)
X_W = 0.5 * (X(x_Ws) + X(x_Wi))
AA_W = (np.polyval(c_inf, 0.5 * (x_Ws + x_Wi)) - np.polyval(c_sup, 0.5 * (x_Ws + x_Wi))) / h_star
M1d_W = gd.mach_da_area(AA_W, aria, "supersonico")
sig_Ws = th_Ws - phi_s             # angolo fra onda e corrente a monte (lato superiore)
sig_Wi = phi_i - th_Wi             # lato inferiore
sig_W = 0.5 * (sig_Ws + sig_Wi)
M_Ws, M_Wi = M_da_nu(th_Ws), M_da_nu(-th_Wi)
M_W = M_da_nu(th_W)
print(f"\nIn W: X = {X_W:.2f}, A/A* = {AA_W:.3f} -> M_1D = {M1d_W:.2f}; theta_w = {th_W:.1f} -> nu = theta_w -> "
      f"M_parete = {M_W:.2f} (sup {M_Ws:.2f}, inf {M_Wi:.2f})")
print(f"Angolo onda-parete: sup {sig_Ws:.1f}, inf {sig_Wi:.1f}, medio {sig_W:.1f} gradi; "
      f"angolo di Mach mu(M_parete) = {mu(M_W):.1f}, mu(M_1D) = {mu(M1d_W):.1f}")
print(f"Se fosse un'onda di Mach: M = 1/sin(sigma) = {1 / math.sin(rad(sig_W)):.2f}, "
      f"e la sua inclinazione sull'asse sarebbe theta_w - mu = {th_W - mu(M_W):.1f} gradi (misurata {-0.5 * apertura:.1f})")

# urto in W (stati 1 -> 2), lato superiore e inferiore
sW_s = gd.State.da_M(aria, M_Ws, theta=th_Ws, p=p_p0(M_Ws), a=1.0, nome="1")
uW_s = gd.urto_obliquo(sW_s, sigma=sig_Ws, famiglia="C-", nome="2")
sW_i = gd.State.da_M(aria, M_Wi, theta=th_Wi, p=p_p0(M_Wi), a=1.0, nome="1'")
uW_i = gd.urto_obliquo(sW_i, sigma=sig_Wi, famiglia="C+", nome="2'")
sW = gd.State.da_M(aria, M_W, theta=th_W, p=p_p0(M_W), a=1.0, nome="1")
uW = gd.urto_obliquo(sW, sigma=sig_W, famiglia="C-", nome="2")
for nome, u in (("sup", uW_s), ("inf", uW_i), ("medio", uW)):
    print(f"Urto in W ({nome}): M1 = {u.monte.M:.2f}, sigma = {u.sigma:.1f}, Mn1 = {u.Mn1:.3f}, "
          f"deviazione = {u.deflessione:.1f}, theta2 = {u.valle.theta:+.1f}, p2/p1 = {u.p2_p1:.2f}, "
          f"M2 = {u.valle.M:.2f}, p2/p0 = {u.valle.p:.3f}, p02/p01 = {u.p02_p01:.4f}")
print(f"Direzione dello strato di scorrimento inferiore: prevista {uW_i.valle.theta:+.1f}, misurata {th_strato:+.1f}")
p_B = uW.valle.p          # pressione nella zona separata (p/p0)
print(f"Pressione nelle zone separate p_B/p0 = {p_B:.3f}; p_W/p0 = {sW.p:.3f}; p_W/p_B = {sW.p / p_B:.2f}")

# --- asse: M dall'invariante theta+nu lungo la C- che parte dalla parete in Q
print("\n--- Asse: C- da un punto Q della parete superiore (onda semplice a parete: nu_Q = theta_Q)")
Y_asse_px = y_asse


def da_Q(xQ):
    """C- da Q (parete superiore) fino all'asse: rettilinea con pendenza media."""
    thQ = theta_sup(xQ)
    nQ = thQ
    MQ = M_da_nu(nQ)
    nA = thQ + nQ                 # theta + nu costante lungo la C-, theta = 0 sull'asse
    MA = M_da_nu(nA)
    incl = 0.5 * ((thQ - mu(MQ)) + (0.0 - mu(MA)))
    YQ = (Y_asse_px - np.polyval(c_sup, xQ)) / h_star
    XA = X(xQ) + YQ / math.tan(rad(-incl))
    return dict(xQ=xQ, XQ=X(xQ), YQ=YQ, thQ=thQ, MQ=MQ, nA=nA, MA=MA, muQ=mu(MQ), muA=mu(MA), incl=incl, XA=XA)


XP = X(x_P)
# primo tentativo come a mano: Q sul piede W (la C- coinciderebbe con l'onda)
t0 = da_Q(x_Ws)
print(f"Primo tentativo, Q = W: theta_Q = {t0['thQ']:.1f}, M_Q = {t0['MQ']:.2f}, nu_asse = {t0['nA']:.1f} -> "
      f"M = {t0['MA']:.2f}; la C- raggiunge l'asse a X = {t0['XA']:.2f}, oltre P (X_P = {XP:.2f})")
xQ = brentq(lambda x: da_Q(x)["XA"] - XP, 556, x_Ws)
QP = da_Q(xQ)
print(f"Punto Q la cui C- arriva in P: X_Q = {QP['XQ']:.3f}, Y_Q = {QP['YQ']:.3f}, theta_Q = {QP['thQ']:.1f}, M_Q = {QP['MQ']:.2f}, "
      f"mu_Q = {QP['muQ']:.1f}; in P nu = {QP['nA']:.1f}, M_P = {QP['MA']:.2f}, mu_P = {QP['muA']:.1f}; "
      f"inclinazione media {QP['incl']:.1f}")
asse = [da_Q(x) for x in np.linspace(556, xQ, 7)]
for d in asse:
    print(f"   Q a X = {d['XQ']:.3f} (theta {d['thQ']:.1f}) -> asse X = {d['XA']:.2f}, M = {d['MA']:.2f}")
M_P = QP["MA"]

# --- incrocio degli urti in P (stati 1 -> 2 -> 3), schema simmetrico con l'angolo medio
sig_P = 0.5 * apertura
print(f"\nIn P: sigma = {sig_P:.1f} (metà dell'angolo fra le onde); 1/sin(sigma) = {1 / math.sin(rad(sig_P)):.2f}")
s1P = gd.State.da_M(aria, M_P, theta=0.0, p=p_p0(M_P), a=1.0, nome="1")
iP = gd.urto_obliquo(s1P, sigma=sig_P, famiglia="C-", nome="2")
rr = gd.riflessione_regolare(s1P, iP.deflessione)
s2P, s3P = rr["2"], rr["3"]
ang_rifl = rr["angolo_riflesso_parete"]
print(f"Urto in P: Mn1 = {iP.Mn1:.3f}, deviazione {iP.deflessione:.1f}, p2/p1 = {iP.p2_p1:.2f}, M2 = {s2P.M:.2f}, p2/p0 = {s2P.p:.3f}")
print(f"Dopo l'incrocio: M3 = {s3P.M:.2f}, p3/p0 = {s3P.p:.3f}, theta3 = {s3P.theta:.1f}; "
      f"urto trasmesso inclinato di {ang_rifl:.1f} gradi sull'asse (misurato {-deg(math.atan(b_i2)):.1f} e {deg(math.atan(b_s2)):.1f})")
print(f"Pressione: p3/p0 = {s3P.p:.3f} contro p_B/p0 = {p_B:.3f}")

# =============================================================================
# 4. FUORI MODULO: verifica con il metodo delle caratteristiche (solo per la sintesi)
# =============================================================================
def moc(Yw, thw, X0, M_ax0, N=30, Xmax=1.8):
    """Rete di caratteristiche piana dalla linea iniziale X = X0 (theta lineare da 0 a theta_w,
    nu uniforme) con parete Yw(X), theta_w(X) e asse di simmetria. Coordinate in unità di h*."""
    def Mn(v):
        return M_da_nu(max(v, 1e-6))

    def interno(A, B):
        thP = 0.5 * ((A[2] + A[3]) + (B[2] - B[3]))
        nP = 0.5 * ((A[2] + A[3]) - (B[2] - B[3]))
        P = [0, 0, thP, nP]
        for _ in range(3):
            mA = math.tan(rad(0.5 * ((A[2] - mu(Mn(A[3]))) + (thP - mu(Mn(nP))))))
            mB = math.tan(rad(0.5 * ((B[2] + mu(Mn(B[3]))) + (thP + mu(Mn(nP))))))
            P[0] = (B[1] - A[1] + mA * A[0] - mB * B[0]) / (mA - mB)
            P[1] = A[1] + mA * (P[0] - A[0])
        return P

    def sull_asse(A):
        nP = A[2] + A[3]
        m = math.tan(rad(0.5 * ((A[2] - mu(Mn(A[3]))) + (-mu(Mn(nP))))))
        return [A[0] - A[1] / m, 0.0, 0.0, nP]

    def a_parete(B):
        Xp = B[0]
        for it in range(4):
            th = thw(Xp) if it else B[2]
            nP = th - (B[2] - B[3])
            m = math.tan(rad(0.5 * ((B[2] + mu(Mn(B[3]))) + (th + mu(Mn(nP))))))
            Xp = brentq(lambda Z: B[1] + m * (Z - B[0]) - Yw(Z), B[0] + 1e-9, B[0] + 1.5)
        th = thw(Xp)
        return [Xp, Yw(Xp), th, th - (B[2] - B[3])]

    n0 = nu(M_ax0)
    riga = [[X0, j / N * Yw(X0), j / N * thw(X0), n0] for j in range(N + 1)]
    tutti = [p[:] for p in riga]
    while riga[0][0] < Xmax:
        mezza = [interno(riga[j + 1], riga[j]) for j in range(N)]
        riga = [sull_asse(mezza[0])] + [interno(mezza[j], mezza[j - 1]) for j in range(1, N)] + [a_parete(mezza[-1])]
        tutti += mezza + riga
    return np.array(tutti)


X_vis = X(905)                     # oltre la finestra: parete prolungata in linea retta (serve solo a
                                   # chiudere la rete, non influenza l'asse fino a P)


def Yw_sup(Xs):
    if Xs > X_vis:
        return Yw_sup(X_vis) + (Xs - X_vis) * math.tan(rad(thw_sup(X_vis)))
    return (y_asse - np.polyval(c_sup, x_gola + Xs * h_star)) / h_star


def thw_sup(Xs):
    return theta_sup(x_gola + min(Xs, X_vis) * h_star)


# caso limite: pareti parallele e corrente uniforme a M = 1,5 -> deve restare uniforme
rete0 = moc(lambda Z: 0.5, lambda Z: 0.0, 0.0, 1.5, N=10, Xmax=1.0)
print(f"\n[fuori modulo] Caso limite MOC (canale a pareti parallele, M = 1,5): "
      f"M min/max = {M_da_nu(rete0[:, 3].min()):.4f}/{M_da_nu(rete0[:, 3].max()):.4f}")
X0m = X(556)
rete = moc(Yw_sup, thw_sup, X0m, M_da_nu(1e-3) * 1.02, N=30, Xmax=XP + 0.05)
ax_m = rete[np.abs(rete[:, 1]) < 1e-12]
ax_m = ax_m[np.argsort(ax_m[:, 0])]
par_m = rete[np.abs(rete[:, 1] - np.array([Yw_sup(z) for z in rete[:, 0]])) < 1e-9]
par_m = par_m[np.argsort(par_m[:, 0])]
M_P_moc = M_da_nu(np.interp(XP, ax_m[:, 0], ax_m[:, 3]))
M_W_moc = M_da_nu(np.interp(X(x_Ws), par_m[:, 0], par_m[:, 3]))
print(f"[fuori modulo] MOC con corrente quasi sonica uniforme in X = {X0m:.3f}: M in P (asse) = {M_P_moc:.2f} "
      f"(stima a mano {M_P:.2f}); M a parete in W = {M_W_moc:.2f} (onda semplice {M_Ws:.2f})")

# =============================================================================
# 5. FIGURE
# =============================================================================
plt.rcParams.update({"font.size": 9})


def px(Xs, Ys):
    return x_gola + Xs * h_star, y_asse - Ys * h_star


# --- figura annotata
fig, ax = plt.subplots(figsize=(7.2, 6.0))
ax.imshow(rett, cmap="gray")
ax.set_xlim(0, Wr)
ax.set_ylim(Hr, 0)
xx = np.linspace(x_gola, 1240, 2)
ax.plot(xx, [y_asse] * 2, ls="-.", color="yellow", lw=0.8)
ax.text(30, y_asse - 10, "asse", color="yellow", fontsize=8)
ax.plot([x_gola, x_gola], [np.polyval(c_sup, 556), np.polyval(c_inf, 558)], color="cyan", lw=1.2)
ax.text(x_gola + 12, y_asse + 200, "gola\n(fori)\nM = 1", color="cyan", fontsize=8)
xs_p = np.linspace(555, 905, 60)
ax.plot(xs_p, np.polyval(c_sup, xs_p), color="lime", lw=1.2)
xs_p = np.linspace(558, 835, 60)
ax.plot(xs_p, np.polyval(c_inf, xs_p), color="lime", lw=1.2)
for (a, b, xa, xb) in [(a_s, b_s, x_Ws, x_P), (a_i, b_i, x_Wi, x_P)]:
    ax.plot([xa, xb], [a + b * xa, a + b * xb], color="red", lw=1.4)
for (a, b, xa, xb) in [(a_s2, b_s2, x_P, 1215), (a_i2, b_i2, x_P, 1205)]:
    ax.plot([xa, xb], [a + b * xa, a + b * xb], color="red", lw=1.0, ls="--")
ax.plot([x_Wi, 1120], [a_l + b_l * x_Wi, a_l + b_l * 1120], color="orange", lw=1.2, ls="--")
yS = np.polyval(c_sup, x_Ws)
ax.plot([x_Ws, 1120], [yS, yS - (1120 - x_Ws) * math.tan(rad(uW_s.valle.theta))], color="orange", lw=1.0, ls=":")
ax.plot(*px(QP["XQ"], QP["YQ"]), "o", color="magenta", ms=4)
ax.plot([px(QP["XQ"], QP["YQ"])[0], x_P], [px(QP["XQ"], QP["YQ"])[1], y_asse], color="magenta", lw=0.9, ls=":")
for (x, y, t, dx, dy) in [(x_Ws, np.polyval(c_sup, x_Ws), "W", -10, -22), (x_Wi, np.polyval(c_inf, x_Wi), "W'", -18, 40),
                          (x_P, y_P, "P", 10, -14), (px(QP["XQ"], QP["YQ"])[0], px(QP["XQ"], QP["YQ"])[1], "Q", -30, -18)]:
    ax.text(x + dx, y + dy, t, color="white", fontsize=10, fontweight="bold")
for (x, y, t) in [(800, 560, "1"), (930, 300, "2"), (900, 830, "2'"), (1180, 540, "3"),
                  (960, 115, "B"), (1000, 985, "B'")]:
    ax.text(x, y, t, color="yellow", fontsize=11, fontweight="bold", ha="center", va="center",
            bbox=dict(boxstyle="circle,pad=0.15", fc="black", ec="yellow", lw=0.6))
ax.annotate(f"σ ≈ {sig_P:.1f}°".replace(".", ","), (x_P - 60, y_P + 10), xytext=(860, 640), color="white", fontsize=8,
            arrowprops=dict(arrowstyle="-", color="white", lw=0.6))
ax.axis("off")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "annotata.png"), dpi=220, bbox_inches="tight", pad_inches=0.02)
plt.close(fig)

# --- schema del campo (metà superiore e inferiore), unità h*
fig, ax = plt.subplots(figsize=(7.0, 4.6))
Xg = np.linspace(X(555), X(905), 100)
ax.plot(Xg, [Yw_sup(z) for z in Xg], "k", lw=1.6)
ax.plot(Xg, [-Yw_sup(z) for z in Xg], "k", lw=1.6)
Xc = np.linspace(-0.25, X(555), 20)
ax.plot(Xc, 0.5 + 0.35 * Xc ** 2, "k", lw=1.6, ls=(0, (2, 2)))
ax.plot(Xc, -(0.5 + 0.35 * Xc ** 2), "k", lw=1.6, ls=(0, (2, 2)))
ax.axhline(0, color="gray", ls="-.", lw=0.7)
ax.plot([0, 0], [-0.5, 0.5], color="C0", lw=1.0)
ax.text(0.01, 0.06, "M = 1", color="C0", fontsize=8)
XWm, YWm = X_W, Yw_sup(X_W)
ax.plot([XWm, XP], [YWm, 0], "r", lw=1.6)
ax.plot([XWm, XP], [-YWm, 0], "r", lw=1.6)
Xe = XP + 0.2
ax.plot([XP, Xe], [0, -0.2 * math.tan(rad(ang_rifl))], "r", lw=1.2)
ax.plot([XP, Xe], [0, 0.2 * math.tan(rad(ang_rifl))], "r", lw=1.2)
ax.plot([XWm, Xe], [YWm, YWm + (Xe - XWm) * math.tan(rad(uW.valle.theta))], color="C1", lw=1.2, ls="--")
ax.plot([XWm, Xe], [-YWm, -YWm - (Xe - XWm) * math.tan(rad(uW.valle.theta))], color="C1", lw=1.2, ls="--")
ax.plot([QP["XQ"], XP], [QP["YQ"], 0], color="C4", lw=0.9, ls=":")
ax.text(QP["XQ"] - 0.02, QP["YQ"] + 0.03, "Q", color="C4")
ax.text(XWm - 0.03, YWm + 0.03, "W")
ax.text(XWm - 0.04, -YWm - 0.07, "W'")
ax.text(XP + 0.01, 0.025, "P")
for (x, y, t) in [(0.45, 0.0, "1"), (0.62, 0.38, "2"), (0.62, -0.38, "2'"), (XP + 0.12, 0.0, "3"),
                  (0.78, 0.62, "B"), (0.78, -0.62, "B'")]:
    ax.text(x, y, t, fontsize=11, ha="center", va="center", bbox=dict(boxstyle="circle,pad=0.2", fc="white", ec="k", lw=0.6))
ax.text(0.86, 0.75, "zona separata", fontsize=7, ha="center")
ax.text(0.2, -0.12, "corrente attaccata,\nin espansione", fontsize=7, ha="center")
ax.set_aspect("equal")
ax.set_xlim(-0.3, Xe + 0.05)
ax.set_ylim(-0.95, 0.95)
ax.set_xlabel("x / h*")
ax.set_ylabel("y / h*")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "schema.pdf"))
plt.close(fig)

# --- distribuzione del numero di Mach
fig, ax = plt.subplots(figsize=(6.6, 4.2))
Xp = np.linspace(X(556), X_W, 40)
Mpar = [M_da_nu(max(thw_sup(z), 1e-6)) for z in Xp]
ax.plot(Xp, Mpar, "C3", lw=1.6, label="parete (Prandtl–Meyer, ν = θ$_w$)")
ax.plot([0, Xp[0]], [1, Mpar[0]], "C3", lw=1.0, ls=":")
Xm = np.linspace(X(560), X_W, 40)
M1d = [gd.mach_da_area(max((np.polyval(c_inf, x_gola + z * h_star) - np.polyval(c_sup, x_gola + z * h_star)) / h_star,
                           1.0 + 1e-6), aria, "supersonico") for z in Xm]
ax.plot(np.r_[0, Xm], np.r_[1, M1d], "C0", lw=1.6, label="media sulla sezione (legge delle aree)")
XA = [d["XA"] for d in asse]
MA = [d["MA"] for d in asse]
ax.plot(XA, MA, "C2", lw=1.6, label="asse (invariante θ+ν lungo le C$^-$)")
ax.plot([0, XA[0]], [1, MA[0]], "C2", lw=1.0, ls=":")
ax.text(0.30, 1.13, "regione transonica\n(non risolta)", fontsize=7, color="C2")
ax.plot([XP - 0.02], [1 / math.sin(rad(sig_P))], "v", color="k", ms=5, mfc="none",
        label="1/sin σ in P (vale solo come limite inferiore)")
ax.plot([X_W, X_W], [uW.valle.M, M_W], "C3", lw=0.8, ls="--")
ax.plot([X_W], [uW.valle.M], "s", color="C3", ms=5, label="dietro l'urto in W (2)")
ax.plot([XP, XP], [s2P.M, M_P], "C2", lw=0.8, ls="--")
ax.plot([XP], [s2P.M], "s", color="C2", ms=5, label="dietro l'urto in P (2)")
ax.plot([XP + 0.03], [s3P.M], "D", color="C2", ms=5, label="dopo l'incrocio (3)")
ax.axvline(X_W, color="gray", lw=0.6, ls=":")
ax.axvline(XP, color="gray", lw=0.6, ls=":")
ax.text(X_W + 0.01, 1.02, "W", fontsize=8)
ax.text(XP + 0.01, 1.02, "P", fontsize=8)
ax.set_xlabel("x / h*  (dalla gola)")
ax.set_ylabel("M")
ax.set_xlim(0, XP + 0.12)
ax.set_ylim(0.95, 1.85)
ax.grid(alpha=0.3)
ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.47, 1.0))
fig.tight_layout()
fig.savefig(os.path.join(FIG, "mach.pdf"))
plt.close(fig)

# --- piano p-theta (p/p0, scala lineare: le pressioni stanno in poco più di un fattore 2)
fig, ax = plt.subplots(figsize=(6.6, 4.4))
for st, fam, col, nm in [(sW, "C-", "C3", "W"), (s1P, "C-", "C2", "P")]:
    u = gd.polare_p_theta(st, famiglia=fam, p_min_rel=0.5)["urto"]
    deb = np.array(u["ramo"]) == "debole"
    ax.plot(np.array(u["theta"])[deb], np.array(u["p"])[deb], color=col, lw=1.3, label=f"polare {fam} di 1 in {nm}")
u = gd.polare_p_theta(s2P, famiglia="C+", p_min_rel=0.8)["urto"]
deb = np.array(u["ramo"]) == "debole"
ax.plot(np.array(u["theta"])[deb], np.array(u["p"])[deb], color="C2", lw=1.0, ls="--", label="polare C+ di 2 in P")
th_is = np.linspace(0, 20, 60)
ax.plot(th_is, [p_p0(M_da_nu(max(t, 1e-6))) for t in th_is], color="gray", lw=1.0, ls=":",
        label="corrente attaccata a parete (ν = θ)")
for st, t, off in [(sW, "1 (W)", (5, -10)), (uW.valle, "2 (W)", (5, 4)), (s1P, "1 (P)", (5, -10)), (s2P, "2 (P)", (-28, 4)),
                   (s3P, "3", (5, 3))]:
    ax.plot(st.theta, st.p, "ko", ms=3.5)
    ax.annotate(t, (st.theta, st.p), textcoords="offset points", xytext=off, fontsize=8)
ax.axhline(p_B, color="C1", lw=0.9, ls="--")
ax.text(8.0, p_B + 0.006, "zone separate B, B'", color="C1", fontsize=8)
ax.set_xlabel("θ [°]")
ax.set_ylabel("p / p$_0$")
ax.set_xlim(-12, 22)
ax.grid(alpha=0.3)
ax.legend(fontsize=7, loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "p_theta.pdf"))
plt.close(fig)
print("\nFigure scritte in", FIG)
