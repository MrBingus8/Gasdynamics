"""Prova del 14/09/2021, esercizio 2 (B65-E02): corrente a M = 5 e incidenza α = 11° su un corpo piano
formato da un segmento sottile di lunghezza h davanti a un corpo rettangolare di altezza b.
Forza sul segmento h (spessore nullo).

Riproduce tutti i numeri e le figure dello svolgimento.
Uso:  python calcoli.py   (dalla cartella dello svolgimento)

Convenzioni (gasdyn_guida.md §1): angoli in gradi, θ antiorario dall'asse x, che coincide con il segmento h;
la corrente indisturbata arriva da sotto con θ1 = +α. Il lato superiore del corpo è una parete sotto la
corrente (onde C+), il lato inferiore una parete sopra la corrente (onde C−). Grandezze adimensionali:
p1 = 1, a1 = 1 (ρ1 = γ); lunghezze in unità di h. Forze per unità di apertura (corpo piano).

Modello (interpretazione della traccia, dichiarata nel PDF): il segmento è in mezzeria del corpo e la corrente
separa allo spigolo d'attacco del segmento, formando due zone di fluido morto (bolle) che vanno dallo
spigolo agli spigoli anteriori del corpo rettangolare. Le bolle hanno pressione uniforme, uguale a quella della
corrente esterna lungo il loro contorno (linea di scorrimento). b/2 = h tan α: il contorno superiore è parallelo
alla corrente indisturbata.
"""
import os
import sys
import math

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

QUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(QUI, "..", ".."))
import gasdyn as gd  # noqa: E402

FIG = os.path.join(QUI, "figure")
G = gd.Gas(1.4)
g, d = G.gamma, G.delta
tand = lambda x: math.tan(math.radians(x))
sind = lambda x: math.sin(math.radians(x))
cosd = lambda x: math.cos(math.radians(x))


def nu(M):
    return gd.prandtl_meyer(M, G)


def riga(nome, st):
    print(f"  {nome:>2}: M = {st.M:.4f}  p = {st.p:.4f}  θ = {st.theta:8.3f}°  ν = {nu(st.M):7.3f}°  "
          f"μ = {st.mu:6.3f}°  p0 = {st.p0:.4g}  ρ = {st.rho:.4f}")


# =============================================================================
# 1. DATI E GEOMETRIA
# =============================================================================
M1, ALFA, B, H = 5.0, 11.0, 0.38876, 1.0
s1 = gd.State.da_M(G, M1, theta=ALFA, p=1.0, a=1.0, nome="1")
yc = B / 2                                   # semi-altezza del corpo (segmento in mezzeria)
th_sup = math.degrees(math.atan2(yc, H))     # contorno della bolla superiore: spigolo -> (h, +b/2)
th_inf = -th_sup                             # contorno della bolla inferiore: spigolo -> (h, -b/2)
print("=== Dati e geometria ===")
riga("1", s1)
print(f"  b/(2h) = {yc:.5f},  tan α = {tand(ALFA):.5f}  =>  contorni delle bolle a ±{th_sup:.4f}°")
print(f"  q1 = γ/2 p1 M1² = {g / 2 * M1**2:.3f} p1")

# misura sul disegno della traccia (solo per la sintesi in chat): colonne di pixel nella foto
# segmento: righe 1987–2003 (centro 1995); corpo: bordo superiore 1958–1973, inferiore 2030–2045
y_seg, y_top, y_bot = 1995.0, 1965.5, 2037.5
print(f"  [disegno] posizione del segmento dall'alto: {(y_seg - y_top) / (y_bot - y_top):.2f} b "
      f"(il disegno non è in scala: b/h disegnato ≈ 0,8)")
for frac in (0.41,):
    a_up = math.degrees(math.atan(frac * B / H))
    a_lo = math.degrees(math.atan((1 - frac) * B / H))
    print(f"  [sensibilità] segmento a {frac:.2f} b dall'alto: contorni a +{a_up:.2f}° e −{a_lo:.2f}°")

# =============================================================================
# 2. BOLLA SUPERIORE A: contorno a θ = +α = θ1 => onda nulla, p_A = p1
# =============================================================================
print("\n=== Bolla superiore A ===")
print(f"  deviazione richiesta θ_sup − θ1 = {th_sup - ALFA:+.5f}° => onda nulla, p_A = p1 = 1")
pA = s1.p

# =============================================================================
# 3. URTO DALLO SPIGOLO, LATO INFERIORE (F §12): deviazione θ1 − θ_inf = 2α = 22°, famiglia C−
#    Tentativi su σ con la relazione θ–σ–M, poi R–H sulla componente normale.
# =============================================================================
dth = ALFA - th_inf
print("\n=== Urto inferiore (1 -> 2) ===")
print(f"  deviazione = {dth:.4f}°;  θmax(M=5) = {gd.theta_max(M1, G)[0]:.2f}° (σ = {gd.theta_max(M1, G)[1]:.2f}°)")
print("  tentativi (relazione θ–σ–M):")
tent = [30.0, 32.0]
for sg in tent:
    print(f"    σ = {sg:.2f}° -> θ = {gd.theta_da_sigma(M1, sg, G):.3f}°")
th_a, th_b = gd.theta_da_sigma(M1, tent[0], G), gd.theta_da_sigma(M1, tent[1], G)
sg_int = tent[0] + (dth - th_a) * (tent[1] - tent[0]) / (th_b - th_a)
print(f"    interpolazione lineare: σ = {sg_int:.3f}° -> θ = {gd.theta_da_sigma(M1, sg_int, G):.3f}°")
u12 = gd.urto_obliquo(s1, theta=dth, famiglia="C-", nome="2")
s2 = u12.valle
Mn1 = u12.Mn1
Mn2 = math.sqrt((1 + d * Mn1**2) / (g * Mn1**2 - d))
r21 = (g + 1) * Mn1**2 / ((g - 1) * Mn1**2 + 2)
print(f"  σ = {u12.sigma:.4f}°, Mn1 = M1 sinσ = {Mn1:.4f}")
print(f"  p2/p1 = 1 + 2γ/(γ+1)(Mn1²−1) = {u12.p2_p1:.4f};  ρ2/ρ1 = {r21:.4f} (modulo: {s2.rho / s1.rho:.4f})")
print(f"  Mn2 = {Mn2:.4f};  M2 = Mn2/sin(σ−θ) = {Mn2 / sind(u12.sigma - dth):.4f} (modulo: {s2.M:.4f})")
print(f"  p02/p01 = {u12.p02_p01:.4f};  Δs/R = {gd.salto_entropia(s1, s2)['ds_R']:.4f}")
print(f"  fronte dell'urto (assoluto) = {u12.angolo_fronte:.3f}°;  θmax(M2) = {gd.theta_max(s2.M, G)[0]:.2f}°")
print(f"  ramo forte (scartato): σ = {gd.sigma_da_theta(M1, dth, G, 'forte'):.2f}°")
riga("2", s2)
pB = s2.p

# =============================================================================
# 4. FORZA SUL SEGMENTO (per unità di apertura): F = (p_B − p_A) h, normale al segmento, verso +y
# =============================================================================
F = (pB - pA) * H
q1 = g / 2 * s1.p * M1**2
print("\n=== Forza sul segmento h ===")
print(f"  F/(p1 h) = p_B − p_A = {F:.4f};  F/(q1 h) = {F / q1:.4f}")
print(f"  assi vento: portanza L = F cos α = {F * cosd(ALFA):.4f} p1 h, resistenza D = F sin α = {F * sind(ALFA):.4f} p1 h")
print(f"  coefficienti: c_L = {F * cosd(ALFA) / q1:.4f}, c_D = {F * sind(ALFA) / q1:.4f}")

# =============================================================================
# 5. SPIGOLI DEL CORPO (solo topologia e piano p–θ): ritorno a θ = 0
#    superiore: C+ (θ + ν costante);  inferiore: C− (θ − ν costante)
# =============================================================================
o3 = gd.onda_verso_theta_2d(s1, 0.0, "C+", nome="3")
o4 = gd.onda_verso_theta_2d(s2, 0.0, "C-", nome="4")
s3, s4 = o3.valle, o4.valle
print("\n=== Spigoli del corpo ===")
print(f"  superiore: {o3.tipo}, ν3 = ν1 + (θ1 − θ3) = {nu(M1):.3f} + {ALFA:.0f} = {nu(s3.M):.3f}°")
print(f"     testa θ1+μ1 = {s1.theta + s1.mu:.3f}°, coda θ3+μ3 = {s3.theta + s3.mu:.3f}°")
riga("3", s3)
print(f"  inferiore: {o4.tipo}, ν4 = ν2 + (θ4 − θ2) = {nu(s2.M):.3f} + {ALFA:.0f} = {nu(s4.M):.3f}°")
print(f"     testa θ2−μ2 = {s2.theta - s2.mu:.3f}°, coda θ4−μ4 = {s4.theta - s4.mu:.3f}°")
riga("4", s4)


def incrocio(P, a1, Q, a2):
    """Intersezione delle rette per P con inclinazione a1 e per Q con inclinazione a2."""
    t1, t2 = tand(a1), tand(a2)
    x = (Q[1] - P[1] + t1 * P[0] - t2 * Q[0]) / (t1 - t2)
    return np.array([x, P[1] + t1 * (x - P[0])])


O = np.array([0.0, 0.0])
Csup, Cinf = np.array([H, yc]), np.array([H, -yc])
y_urto_h = H * tand(u12.angolo_fronte)
X = incrocio(O, u12.angolo_fronte, Cinf, s2.theta - s2.mu)
print(f"  urto a x = h: y = {y_urto_h:.4f} h (spigolo inferiore a {-yc:.4f} h)")
print(f"  testa del ventaglio inferiore sull'urto in x = {X[0]:.3f} h, y = {X[1]:.3f} h")

# =============================================================================
# 6. CONFRONTO: SEGMENTO ISOLATO (senza corpo), stessa incidenza
# =============================================================================
lo = gd.onda_verso_theta_2d(s1, 0.0, "C-")     # faccia inferiore: urto da 11°
up = gd.onda_verso_theta_2d(s1, 0.0, "C+")     # faccia superiore: espansione di 11°
Fi = lo.valle.p - up.valle.p
print("\n=== Confronto: lastra isolata ===")
print(f"  inferiore {lo.tipo} σ = {lo.sigma:.2f}°, p = {lo.valle.p:.4f};  superiore {up.tipo}, p = {up.valle.p:.4f}")
print(f"  F/(p1 h) = {Fi:.4f}, F/(q1 h) = {Fi / q1:.4f};  rapporto corpo/isolata = {F / Fi:.3f}")

# =============================================================================
# 7. FIGURE
# =============================================================================
# 7a. traccia originale (ridotta)
im = Image.open(os.path.join(FIG, "originale.png"))
if im.width > 1400:
    im = im.resize((1400, round(im.height * 1400 / im.width)), Image.LANCZOS)
    im.save(os.path.join(FIG, "originale.png"))

# 7b. piano fisico
L_corpo = 3.2
fig, ax = plt.subplots(figsize=(10.5, 6.2))
ax.fill([0, H, H], [0, yc, 0], color="#d9e6f2", zorder=1)
ax.fill([0, H, H], [0, -yc, 0], color="#f6d9d4", zorder=1)
ax.add_patch(plt.Rectangle((H, -yc), L_corpo, 2 * yc, fc="0.75", ec="k", lw=1.5, zorder=3))
ax.plot([0, H], [0, 0], "k", lw=3, zorder=4)
LS = dict(color="0.25", lw=1.1, ls=(0, (5, 2)), zorder=2)
ax.plot([0, H], [0, yc], **LS)
ax.plot([0, H], [0, -yc], **LS)
ax.plot([0, X[0]], [0, X[1]], color="k", lw=2.0)                         # urto rettilineo fino a X
tt = np.linspace(0, 1, 40)                                                # oltre X: urto che si indebolisce
ang = u12.angolo_fronte + 8 * tt
xx = X[0] + np.cumsum(np.cos(np.radians(ang))) * 0.035
yy = X[1] + np.cumsum(np.sin(np.radians(ang))) * 0.035
ax.plot(xx, yy, color="k", lw=1.6, ls="--")
for k, a_ in enumerate(np.linspace(s1.theta + s1.mu, s3.theta + s3.mu, 5)):   # ventaglio C+ superiore
    E = Csup + 1.9 * np.array([cosd(a_), sind(a_)])
    ax.plot([Csup[0], E[0]], [Csup[1], E[1]], color="tab:blue", lw=1.0 if k in (0, 4) else 0.6)
for k, a_ in enumerate(np.linspace(s2.theta - s2.mu, s4.theta - s4.mu, 5)):   # ventaglio C− inferiore
    E = incrocio(Cinf, a_, O, u12.angolo_fronte) if k == 0 else Cinf + 1.15 * np.array([cosd(a_), sind(a_)])
    ax.plot([Cinf[0], E[0]], [Cinf[1], E[1]], color="tab:red", lw=1.0 if k in (0, 4) else 0.6)
for x0 in (-1.25, -0.95):                                                 # corrente indisturbata
    y0 = -0.75 if x0 < -1 else -0.35
    ax.annotate("", xy=(x0 + 0.5 * cosd(ALFA), y0 + 0.5 * sind(ALFA)), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="->", lw=1.2))
ax.text(-1.3, -0.95, "M$_1$ = 5, θ$_1$ = α = 11°", fontsize=9)
ax.annotate("", xy=(0.5, 0.17), xytext=(0.5, 0.02), arrowprops=dict(arrowstyle="-|>", lw=2, color="tab:green"))
ax.text(0.53, 0.115, "F", color="tab:green", fontsize=12, fontweight="bold")
lab = [("1", (-0.6, 0.45)), ("2", (0.97, -0.31)), ("3", (2.6, 0.55)), ("4", (2.55, -0.45)),
       ("A", (0.72, 0.075)), ("B", (0.72, -0.085))]
for n, (x, y) in lab:
    ax.text(x, y, n, fontsize=10, ha="center", va="center",
            bbox=dict(boxstyle="circle,pad=0.18", fc="w", ec="0.4", lw=0.6), zorder=6)
ax.text(0.03, 0.03, "O", fontsize=9)
ax.text(0.75, -0.62, f"urto, σ = {u12.sigma:.1f}°", fontsize=8, rotation=u12.angolo_fronte)
ax.text(1.95, 0.95, "espansione C$^+$", fontsize=8, color="tab:blue")
ax.text(1.85, -1.08, "espansione C$^-$", fontsize=8, color="tab:red")
ax.set_aspect("equal")
ax.set_xlim(-1.4, 4.2)
ax.set_ylim(-1.2, 1.3)
ax.set_xlabel("x / h")
ax.set_ylabel("y / h")
ax.grid(alpha=0.15)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "campo.pdf"))
plt.close(fig)


# 7c. piano p–θ, scala lineare e logaritmica
def piano_pt(ax, log):
    pol = gd.polare_p_theta(s1, "entrambe", p_min_rel=0.15)
    for fam, col in (("C-", "tab:red"), ("C+", "tab:blue")):
        u = pol[fam]["urto"]
        deb = u["ramo"] == "debole"
        ax.plot(u["theta"][deb], u["p"][deb], color=col, lw=1.6, label=f"polare di 1, {fam}")
        ax.plot(u["theta"][~deb], u["p"][~deb], color=col, lw=0.8, ls=":")
    e1 = pol["C+"]["espansione"]
    m = e1["theta"] >= -1.5
    ax.plot(e1["theta"][m], e1["p"][m], color="tab:blue", lw=1.2, ls="--", label="espansione C+ da 1")
    e2 = gd.polare_p_theta(s2, "C-", p_min_rel=0.3)["espansione"]
    m = e2["theta"] <= 1.5
    ax.plot(e2["theta"][m], e2["p"][m], color="tab:red", lw=1.2, ls="--", label="espansione C− da 2")
    for th, txt in ((-ALFA, "contorno B"), (0.0, "segmento e corpo"), (ALFA, "contorno A")):
        ax.axvline(th, color="0.55", lw=0.8, ls="-.")
    for st, n, off in ((s1, "1 = A", (6, 4)), (s2, "2 = B", (6, 2)), (s3, "3", (6, 2)), (s4, "4", (7, -13))):
        ax.plot(st.theta, st.p, "ko", ms=4.5, zorder=5)
        ax.annotate(n, (st.theta, st.p), textcoords="offset points", xytext=off, fontsize=10)
    ax.plot(0, lo.valle.p, marker="s", ms=5, mfc="w", mec="k", zorder=5)
    ax.annotate("lastra isolata", (0, lo.valle.p), textcoords="offset points", xytext=(7, 6), fontsize=8)
    ax.annotate("", xy=(-ALFA, pB), xytext=(-ALFA, pA),
                arrowprops=dict(arrowstyle="<->", color="tab:green", lw=1.4))
    ax.text(-ALFA + 0.6, 0.5 * (pA + pB) if not log else math.sqrt(pA * pB), "F/(p$_1$h)",
            color="tab:green", fontsize=9)
    ax.set_xlabel("θ [°]")
    ax.set_ylabel("p / p$_1$")
    ax.set_xlim(-32, 40)
    ax.grid(alpha=0.25)
    if log:
        ax.set_yscale("log")
        ax.set_ylim(0.15, 20)
    else:
        ax.set_ylim(0, 14)


fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.2))
piano_pt(axs[0], log=False)
piano_pt(axs[1], log=True)
axs[0].set_title("scala lineare")
axs[1].set_title("scala logaritmica")
axs[0].legend(fontsize=7.5, loc="upper right")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "p_theta.pdf"))
plt.close(fig)

# 7d. pressione sulle due facce del segmento (corpo vs lastra isolata)
fig, ax = plt.subplots(figsize=(7.5, 3.6))
xs = np.array([0, H])
ax.plot(xs, [pB, pB], color="tab:red", lw=2.2, label="faccia inferiore (bolla B)")
ax.plot(xs, [pA, pA], color="tab:blue", lw=2.2, label="faccia superiore (bolla A)")
ax.fill_between(xs, pA, pB, color="tab:green", alpha=0.15)
ax.plot(xs, [lo.valle.p] * 2, color="tab:red", lw=1.2, ls="--", label="lastra isolata, inferiore")
ax.plot(xs, [up.valle.p] * 2, color="tab:blue", lw=1.2, ls="--", label="lastra isolata, superiore")
ax.text(0.5, 0.5 * (pA + pB), f"F/(p$_1$h) = {F:.2f}", ha="center", color="tab:green", fontsize=10)
ax.set_xlabel("x / h")
ax.set_ylabel("p / p$_1$")
ax.set_xlim(-0.05, 1.05)
ax.set_ylim(0, 9)
ax.grid(alpha=0.25)
ax.legend(fontsize=7.5, loc="center right", bbox_to_anchor=(1.0, 0.62))
fig.tight_layout()
fig.savefig(os.path.join(FIG, "carico.pdf"))
plt.close(fig)
print("\nFigure salvate in", FIG)
