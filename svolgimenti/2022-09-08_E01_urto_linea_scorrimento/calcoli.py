"""Prova dell'08/09/2022, esercizio 1: urto generato da un cuneo che incide su una linea di scorrimento
fra due correnti supersoniche alla stessa pressione (sopra M = 3, sotto M = 8, cuneo di 20° sotto).

Riproduce tutti i numeri e le figure dello svolgimento.
Uso:  python calcoli.py   (dalla cartella dello svolgimento)

Convenzioni (gasdyn_guida.md §1): angoli in gradi, θ antiorario dall'asse x (verso delle correnti a monte);
la corrente inferiore emette onde C+ (parete sotto), la superiore riceve onde C+ trasmesse; le onde riflesse
dalla linea di scorrimento verso il basso sono C−. Grandezze adimensionali: p2 = 1, ρ2 = 1 (dati); p1 = 1.
ρ1 non è assegnata: non serve, perché p–θ e odografo normalizzato (V/a*) dipendono solo da M, p, θ.
Per costruire lo stato 1 si pone ρ1 = 1 (valore arbitrario, ininfluente su tutti i risultati).
"""
import os
import sys
import math

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

QUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(QUI, "..", ".."))
import gasdyn as gd  # noqa: E402

FIG = os.path.join(QUI, "figure")
G = gd.Gas(1.4)
d = G.delta
tand = lambda x: math.tan(math.radians(x))
sind = lambda x: math.sin(math.radians(x))


def nu(M):
    return gd.prandtl_meyer(M, G)


def p0_p(M):
    return gd.isentropiche(M, G)["p0_p"]


def riga(nome, st):
    print(f"  {nome:>3}: M = {st.M:.4f}  p = {st.p:.4f}  θ = {st.theta:8.4f}°  ν = {nu(st.M):7.3f}°  "
          f"μ = {st.mu:6.3f}°  M* = {st.M_star:.4f}  p0 = {st.p0:.4g}")


# =============================================================================
# 1. DATI
# =============================================================================
TH_W = 20.0
s1 = gd.State.da_M(G, 3.0, p=1.0, rho=1.0, nome="1")      # sopra (ρ1 arbitraria)
s2 = gd.State.da_M(G, 8.0, p=1.0, rho=1.0, nome="2")      # sotto (dati p2 = ρ2 = 1)
print("=== Dati ===")
riga("1", s1)
riga("2", s2)
print(f"  θmax(M=8) = {gd.theta_max(8.0, G)[0]:.2f}°,  θmax(M=3) = {gd.theta_max(3.0, G)[0]:.2f}°")
print(f"  a2 = {s2.a:.4f}, V2 = {s2.V:.4f};  ν_max = {gd.nu_max(G):.2f}°")

# =============================================================================
# 2. URTO DEL CUNEO: STATO 3 (F §12)
# =============================================================================
u3 = gd.urto_obliquo(s2, theta=TH_W, famiglia="C+", nome="3")
s3 = u3.valle
Mn2 = s2.M * sind(u3.sigma)
Mn3 = math.sqrt((1 + d * Mn2**2) / (G.gamma * Mn2**2 - d))
print("\n=== Urto del cuneo (2 -> 3) ===")
print(f"  σ = {u3.sigma:.3f}°, Mn = {Mn2:.4f}, p3/p2 = {u3.p2_p1:.4f}, Mn3 = {Mn3:.4f}, "
      f"M3 = Mn3/sin(σ−θ) = {Mn3 / sind(u3.sigma - TH_W):.4f}")
print(f"  ρ3/ρ2 = {s3.rho / s2.rho:.4f}, p03/p02 = {u3.p02_p01:.4f}, p03/p3 = {p0_p(s3.M):.3f}")
riga("3", s3)

# =============================================================================
# 3. PRIMA INTERAZIONE CON LA LINEA DI SCORRIMENTO: STATI 4 (sopra) E 5 (sotto)
#    p4 = p5, θ4 = θ5 (F §13). Iterazione sull'angolo σ dell'urto trasmesso:
#    σ -> p4, θ4 (urto obliquo, esplicito); p5 = p4 -> M5 (isentropica da p03) -> ν5 -> θ5 = θ3 + ν5 − ν3
#    (onda semplice C−: θ − ν costante). Residuo θ4 − θ5.
# =============================================================================
nu3 = nu(s3.M)


def tentativo(sig):
    Mn = s1.M * sind(sig)
    p4 = s1.p * (1 + 2 * G.gamma / (G.gamma + 1) * (Mn**2 - 1))
    th4 = gd.theta_da_sigma(s1.M, sig, G)
    M5 = gd.mach_da_rapporto(s3.p0 / p4, G, "p0/p")
    th5 = s3.theta + nu(M5) - nu3
    return dict(sigma=sig, Mn=Mn, p4=p4, th4=th4, M5=M5, nu5=nu(M5), th5=th5, res=th4 - th5)


print("\n=== Prima interazione: iterazione sull'angolo dell'urto trasmesso ===")
sig0 = gd.sigma_da_theta(s1.M, TH_W, G)       # primo tentativo: linea di scorrimento parallela al cuneo
print(f"  primo tentativo θ = {TH_W}° => σ(0) = {sig0:.3f}°")
iters = [tentativo(sig0), tentativo(50.0)]
while abs(iters[-1]["res"]) > 1e-3:
    a, b = iters[-2], iters[-1]
    s_new = b["sigma"] - b["res"] * (b["sigma"] - a["sigma"]) / (b["res"] - a["res"])   # secante
    iters.append(tentativo(s_new))
for k, it in enumerate(iters):
    print(f"  k={k}: σ = {it['sigma']:.3f}°, Mn = {it['Mn']:.4f}, p4 = {it['p4']:.4f}, θ4 = {it['th4']:.3f}°, "
          f"M5 = {it['M5']:.4f}, ν5 = {it['nu5']:.3f}°, θ5 = {it['th5']:.3f}°, θ4−θ5 = {it['res']:+.4f}°")

r1 = gd.riemann_stazionario(s1, s3)           # controllo con il solutore del modulo
s4 = r1.onda_sopra.valle.con(nome="4")
s5 = r1.onda_sotto.valle.con(nome="5")
print("  controllo (riemann_stazionario):", f"p = {r1.p:.5f}, θ = {r1.theta:.5f}°")
print("  ", r1.onda_sopra)
print("  ", r1.onda_sotto)
riga("4", s4)
riga("5", s5)
print("  tutte le intersezioni delle polari di 1 (C+) e 3 (C−):")
for x in gd.intersezioni_polari(s1, "C+", s3, "C-"):
    print(f"    p = {x['p']:.4f}, θ = {x['theta']:.3f}°, sopra {x['onda_A'].tipo} {x['onda_A'].ramo}, "
          f"sotto {x['onda_B'].tipo}")
u_up_20 = gd.urto_obliquo(s1, theta=TH_W, famiglia="C+")
print(f"  stima a priori: urto in 1 con θ = 20° => p = {u_up_20.valle.p:.3f} < p3 = {s3.p:.3f} "
      f"=> la corrente 3 deve espandersi")
Mn1t = r1.onda_sopra.Mn1
print(f"  Mn urto trasmesso = {Mn1t:.4f}, p4/p1 = {s4.p:.4f}, p04/p01 = {r1.onda_sopra.p02_p01:.4f}")
print(f"  ventaglio riflesso C− in P: testa θ3−μ3 = {s3.theta - s3.mu:.3f}°, coda θ5−μ5 = {s5.theta - s5.mu:.3f}°")

# limiti: parete superiore rigida (θ = 0) e contorno a pressione costante p = 1
rig = gd.onda_verso_theta_2d(s3, 0.0, "C-")
pc = gd.onda_verso_p_2d(s3, 1.0, "C-")
print(f"  limite parete rigida (θ=0): {rig.tipo}, p = {rig.valle.p:.3f}")
print(f"  limite contorno isobaro p=1: {pc.tipo}, θ = {pc.valle.theta:.3f}°, M = {pc.valle.M:.3f}")

# =============================================================================
# 4. RIFLESSIONE DEL VENTAGLIO SULLA PARETE DEL CUNEO: STATO 6 (θ6 = 20°, onda C+: θ + ν costante)
# =============================================================================
o6 = gd.onda_verso_theta_2d(s5, TH_W, "C+", nome="6")
s6 = o6.valle
print("\n=== Riflessione sulla parete (5 -> 6) ===")
print(f"  ν6 = ν5 + (θ5 − θ6) = {nu(s5.M):.3f} + {s5.theta - TH_W:.3f} = {nu(s6.M):.3f}°")
print(" ", o6)
riga("6", s6)
print(f"  p6/p03 = 1/{p0_p(s6.M):.1f};  testa riflessa θ5+μ5 = {s5.theta + s5.mu:.3f}°, "
      f"coda θ6+μ6 = {s6.theta + s6.mu:.3f}°;  testa da Q1 con θ3+μ3 = {s3.theta + s3.mu:.3f}°")

# =============================================================================
# 5. SECONDA INTERAZIONE CON LA LINEA DI SCORRIMENTO: STATI 7 (sopra) E 8 (sotto)
# =============================================================================
r2 = gd.riemann_stazionario(s4, s6)
s7 = r2.onda_sopra.valle.con(nome="7")
s8 = r2.onda_sotto.valle.con(nome="8")
print("\n=== Seconda interazione (4 | 6) ===")
print(" ", r2)
riga("7", s7)
riga("8", s8)
o9 = gd.onda_verso_theta_2d(s8, TH_W, "C+", nome="9")    # riflessione della compressione sulla parete
s9 = o9.valle
print(" ", o9)
riga("9", s9)
K = lambda st: G.gamma * st.p * st.M**2 / math.sqrt(st.M**2 - 1)    # dp/dθ per onde deboli [1/rad]
print(f"  pendenze dp/dθ: K1 = {K(s1):.2f}, K3 = {K(s3):.2f}, K4 = {K(s4):.2f}, K5 = {K(s5):.2f}; "
      f"(K5−K4)/(K5+K4) = {(K(s5) - K(s4)) / (K(s5) + K(s4)):.3f}")
# controllo isentropico della compressione riflessa (urto debole)
oi8 = gd.onda_verso_p_2d(s6, s8.p, "C-", modello="isentropico")
print(f"  compressione riflessa isentropica: θ = {oi8.valle.theta:.3f}° (urto: {s8.theta:.3f}°)")

# =============================================================================
# 6. SUCCESSIONE DELLE RIFLESSIONI (solo per la sintesi) E STATO ASINTOTICO
# =============================================================================
print("\n=== Successione delle riflessioni (fuori modulo: catena di riemann_stazionario) ===")
up, lo = s4, s6
seq = [(r1.p, r1.theta)]
for k in range(10):
    rr = gd.riemann_stazionario(up, lo)
    up = rr.onda_sopra.valle
    lo = gd.onda_verso_theta_2d(rr.onda_sotto.valle, TH_W, "C+").valle
    seq.append((rr.p, rr.theta))
for k, (p, th) in enumerate(seq):
    print(f"  interazione {k + 1}: p = {p:.4f}, θ linea di scorrimento = {th:.3f}°, θ−20 = {th - 20:+.3f}°")
rap = [(seq[k + 1][1] - 20) / (seq[k][1] - 20) for k in range(len(seq) - 1)]
print("  rapporto fra deviazioni successive (θ−20):", ", ".join(f"{x:.3f}" for x in rap[:6]))

print("\n=== Stato asintotico (θ = 20°) ===")
ninf_up = nu(s4.M) + (s4.theta - TH_W)          # onde C+ nella corrente superiore: θ + ν costante
Minf_up = gd.mach_da_pm(ninf_up, G)
pinf = s4.p0 / p0_p(Minf_up)
Minf_lo = gd.mach_da_rapporto(s3.p0 / pinf, G, "p0/p")
print(f"  sopra: ν∞ = ν4 + (θ4 − 20) = {nu(s4.M):.3f} + {s4.theta - TH_W:.3f} = {ninf_up:.3f}° => M∞ = {Minf_up:.4f}")
print(f"  p04 = p4·(1+δM4²)^3.5 = {s4.p0:.4f};  p∞ = p04/(1+δM∞²)^3.5 = {pinf:.4f}")
print(f"  sotto: p03/p∞ = {s3.p0 / pinf:.2f} => M∞,sotto = {Minf_lo:.4f}, ν = {nu(Minf_lo):.3f}°")
print(f"  confronto: urto unico da 1 con θ = 20°: p = {u_up_20.valle.p:.4f}, M = {u_up_20.valle.M:.4f}, "
      f"p0/p01 = {u_up_20.p02_p01:.4f}  (percorso reale: p04/p01 = {r1.onda_sopra.p02_p01:.4f})")
print(f"  confronto: limite della successione p = {seq[-1][0]:.4f}")
sinf_up = gd.State.da_M(G, Minf_up, theta=TH_W, p=pinf, rho=1.0, nome="∞s")
sinf_lo = gd.State.da_M(G, Minf_lo, theta=TH_W, p=pinf, rho=1.0, nome="∞i")

# =============================================================================
# 7. GEOMETRIA NEL PIANO FISICO (lunghezze in unità h = quota della linea di scorrimento sul vertice)
#    caratteristiche come semirette con le pendenze degli stati uniformi (metodo_docente c11)
# =============================================================================
def incrocio(P, a1, Q, a2):
    """Intersezione delle rette per P con inclinazione a1 e per Q con inclinazione a2 (F §15)."""
    t1, t2 = tand(a1), tand(a2)
    x = (Q[1] - P[1] + t1 * P[0] - t2 * Q[0]) / (t1 - t2)
    return np.array([x, P[1] + t1 * (x - P[0])])


O = np.array([0.0, 0.0])
P = np.array([1.0 / tand(u3.sigma), 1.0])
Q1 = incrocio(P, s3.theta - s3.mu, O, TH_W)
Q2 = incrocio(P, s5.theta - s5.mu, O, TH_W)
S1 = incrocio(Q1, s5.theta + s5.mu, P, s5.theta)
th_corda = 0.5 * (s5.theta + s7.theta)          # corda della linea di scorrimento nel ventaglio riflesso
S2 = incrocio(Q2, s6.theta + s6.mu, S1, th_corda)
T = incrocio(P, r1.onda_sopra.sigma, S1, s4.theta + s4.mu)   # testa del ventaglio trasmesso sull'urto
print("\n=== Geometria (unità h) ===")
for n, X in [("P", P), ("Q1", Q1), ("Q2", Q2), ("S1", S1), ("S2", S2), ("T", T)]:
    print(f"  {n}: x = {X[0]:.3f}, y = {X[1]:.3f}")
print(f"  x_P = h/tan σ3 = {P[0]:.3f} h")

# =============================================================================
# 8. FIGURE
# =============================================================================
COL_SU, COL_GIU = "tab:blue", "tab:red"
STATI = [s1, s2, s3, s4, s5, s6, s7, s8]


def polari_pt(ax, log):
    # corrente superiore: polare di 1 (C+) e ramo d'espansione C+ da 4
    pol1 = gd.polare_p_theta(s1, "C+", p_min_rel=0.5)
    u = pol1["urto"]
    deb = u["ramo"] == "debole"
    ax.plot(u["theta"][deb], u["p"][deb], color=COL_SU, lw=1.6, label="polare di 1 (C+, sopra)")
    ax.plot(u["theta"][~deb], u["p"][~deb], color=COL_SU, lw=0.9, ls=":")
    e4 = gd.polare_p_theta(s4, "C+", p_min_rel=0.45)["espansione"]
    ax.plot(e4["theta"], e4["p"], color=COL_SU, lw=1.2, ls="--", label="espansione C+ da 4")
    # corrente inferiore: polare di 2 (C+) fino a 3, polare di 3 (C−), espansione C+ da 5, polare C− di 6
    pol2 = gd.polare_p_theta(s2, "C+", espansione=False)
    u = pol2["urto"]
    m = (u["ramo"] == "debole") & (u["theta"] <= 32)
    ax.plot(u["theta"][m], u["p"][m], color=COL_GIU, lw=1.6, label="polare di 2 (C+, sotto)")
    pol3 = gd.polare_p_theta(s3, "C-", p_min_rel=0.25)
    u = pol3["urto"]
    m = (u["ramo"] == "debole") & (u["p"] < 60)
    ax.plot(u["theta"][m], u["p"][m], color="tab:orange", lw=1.4, label="polare di 3 (C−)")
    m = (u["ramo"] == "forte") & (u["p"] < 60)
    ax.plot(pol3["espansione"]["theta"], pol3["espansione"]["p"], color="tab:orange", lw=1.2, ls="--")
    e5 = gd.polare_p_theta(s5, "C+", p_min_rel=0.3)["espansione"]
    m = e5["theta"] >= 19.0
    ax.plot(e5["theta"][m], e5["p"][m], color="tab:purple", lw=1.2, ls="--", label="espansione C+ da 5")
    pol6 = gd.polare_p_theta(s6, "C-", p_min_rel=0.8)
    u = pol6["urto"]
    m = (u["ramo"] == "debole") & (u["theta"] > 12)
    ax.plot(u["theta"][m], u["p"][m], color="tab:green", lw=1.4, label="polare di 6 (C−)")
    # rette θ = 20° (parete) e stati
    ax.axvline(TH_W, color="0.5", lw=0.8, ls="-.")
    for st in STATI:
        ax.plot(st.theta, st.p, "ko", ms=4, zorder=5)
    for nome, st, off in [("1, 2", s1, (6, 4)), ("3", s3, (6, 2)), ("4, 5", s4, (6, 3)), ("6", s6, (6, -10)),
                          ("7, 8", s7, (-30, -11))]:
        ax.annotate(nome, (st.theta, st.p), textcoords="offset points", xytext=off, fontsize=10)
    ax.plot(TH_W, pinf, marker="s", color="k", ms=4, mfc="w", zorder=5)
    ax.annotate("∞", (TH_W, pinf), textcoords="offset points", xytext=(-14, 2), fontsize=11)
    xs = gd.intersezioni_polari(s1, "C+", s3, "C-")
    forte = [x for x in xs if x["onda_A"].ramo == "forte"][0]
    ax.plot(forte["theta"], forte["p"], "x", color="0.3", ms=6)
    ax.annotate("ramo forte\n(scartato)", (forte["theta"], forte["p"]), textcoords="offset points",
                xytext=(-62, 6), fontsize=7, color="0.3")
    ax.set_xlabel("θ [°]")
    ax.set_ylabel("p / p$_1$")
    ax.set_xlim(-2, 36)
    ax.grid(alpha=0.25)
    if log:
        ax.set_yscale("log")
        ax.set_ylim(0.8, 25)
    else:
        ax.set_ylim(0, 17)


fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.2))
polari_pt(axs[0], log=False)
polari_pt(axs[1], log=True)
axs[0].set_title("scala lineare")
axs[1].set_title("scala logaritmica")
axs[0].legend(fontsize=7.5, loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "p_theta.pdf"))
plt.close(fig)


def odo(ax, zoom):
    # cerchi sonico e limite
    t = np.linspace(-0.2, 1.2, 300)
    rl = math.sqrt((G.gamma + 1) / (G.gamma - 1))
    ax.plot(np.cos(t), np.sin(t), color="0.5", lw=0.8, ls="--")
    ax.plot(rl * np.cos(t), rl * np.sin(t), color="0.5", lw=0.8, ls="-.")
    if not zoom:
        ax.annotate("M* = 1", (math.cos(1.0), math.sin(1.0)), xytext=(-40, 0), textcoords="offset points",
                    fontsize=8, color="0.4")
        ax.annotate("M* limite", (rl * math.cos(0.62), rl * math.sin(0.62)), xytext=(4, 0),
                    textcoords="offset points", fontsize=8, color="0.4")
    o1 = gd.polare_odografo(s1, "C+", p_min_rel=0.5)
    u = o1["urto"]
    deb = u["ramo"] == "debole"
    ax.plot(u["us"][deb], u["vs"][deb], color=COL_SU, lw=1.6, label="polare di 1 (sopra)")
    ax.plot(u["us"][~deb], u["vs"][~deb], color=COL_SU, lw=0.9, ls=":")
    e4 = gd.polare_odografo(s4, "C+", p_min_rel=0.45)["espansione"]
    ax.plot(e4["us"], e4["vs"], color=COL_SU, lw=1.2, ls="--", label="espansione C+ da 4")
    o2 = gd.polare_odografo(s2, "C+", espansione=False)["urto"]
    deb = o2["ramo"] == "debole"
    ax.plot(o2["us"][deb], o2["vs"][deb], color=COL_GIU, lw=1.6, label="polare di 2 (sotto)")
    ax.plot(o2["us"][~deb], o2["vs"][~deb], color=COL_GIU, lw=0.9, ls=":")
    e3 = gd.polare_odografo(s3, "C-", p_min_rel=0.25)["espansione"]
    ax.plot(e3["us"], e3["vs"], color="tab:orange", lw=1.3, ls="--", label="espansione C− da 3")
    e5 = gd.polare_odografo(s5, "C+", p_min_rel=0.3)["espansione"]
    th5 = np.degrees(np.arctan2(e5["vs"], e5["us"]))
    m = th5 >= 19.0
    ax.plot(e5["us"][m], e5["vs"][m], color="tab:purple", lw=1.3, ls="--", label="espansione C+ da 5")
    o6p = gd.polare_odografo(s6, "C-", p_min_rel=0.9)["urto"]
    th6 = np.degrees(np.arctan2(o6p["vs"], o6p["us"]))
    m = (o6p["ramo"] == "debole") & (th6 > 12)
    ax.plot(o6p["us"][m], o6p["vs"][m], color="tab:green", lw=1.3, label="polare di 6 (C−)")
    for th, lab in [(TH_W, "20°"), (s4.theta, f"{s4.theta:.1f}°"), (s7.theta, f"{s7.theta:.1f}°")]:
        ax.plot([0, 2.6 * math.cos(math.radians(th))], [0, 2.6 * math.sin(math.radians(th))],
                color="0.6", lw=0.7, ls="-.")
        if not zoom:
            ax.annotate(lab, (0.75 * math.cos(math.radians(th)), 0.75 * math.sin(math.radians(th))),
                        textcoords="offset points", xytext=(-4, 4), fontsize=7, color="0.4")
    for st in STATI + [sinf_up, sinf_lo]:
        x, y = st.M_star * math.cos(math.radians(st.theta)), st.M_star * math.sin(math.radians(st.theta))
        ax.plot(x, y, "ko" if "∞" not in st.nome else "ks", ms=4, mfc="k" if "∞" not in st.nome else "w",
                zorder=5)
        lab = {"∞s": "∞ (sopra)", "∞i": "∞ (sotto)"}.get(st.nome, st.nome)
        offs = {"1": (4, 4), "2": (4, 4), "3": (-12, 4), "4": (5, 2), "5": (-12, 3), "6": (5, 2),
                "7": (4, -10), "8": (5, 2), "∞s": (4, 3), "∞i": (-52, 3)}
        if zoom and st.nome in ("1", "4", "7", "∞s"):
            continue
        if (not zoom) and st.nome in ("3", "5", "6", "8", "∞i"):
            continue
        ax.annotate(lab, (x, y), textcoords="offset points", xytext=offs[st.nome], fontsize=9)
    ax.set_aspect("equal")
    ax.grid(alpha=0.25)
    ax.set_xlabel("u* = u/a*")
    ax.set_ylabel("v* = v/a*")
    if zoom:
        ax.set_xlim(1.85, 2.45)
        ax.set_ylim(0.55, 1.25)
        ax.set_title("ingrandimento: corrente inferiore")
    else:
        ax.set_xlim(0, 2.6)
        ax.set_ylim(-0.05, 1.65)
        ax.set_title("odografo (ogni corrente con il proprio a*)")
        ax.legend(fontsize=7, loc="upper left")


fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.0), gridspec_kw={"width_ratios": [1.5, 1]})
odo(axs[0], zoom=False)
odo(axs[1], zoom=True)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "odografo.pdf"))
plt.close(fig)

# --- piano fisico (scala in h) ---
fig, ax = plt.subplots(figsize=(11, 5.6))
xmax = 12.5
ax.fill([0, xmax, xmax], [0, 0, xmax * tand(TH_W)], color="0.85")
ax.plot([0, xmax], [0, xmax * tand(TH_W)], "k", lw=1.5)
ax.plot([0, xmax], [0, 0], "k", lw=0.8)
LS = dict(color="0.2", lw=1.3, ls=(0, (6, 2)))
ax.plot([-3, P[0]], [1, 1], **LS)                                   # linea di scorrimento iniziale
ax.plot([P[0], S1[0]], [1, S1[1]], **LS)                            # tratto a θ5
ax.plot([S1[0], S2[0]], [S1[1], S2[1]], **LS)                       # corda nel ventaglio riflesso
S3 = S2 + 2.4 * np.array([math.cos(math.radians(s7.theta)), math.sin(math.radians(s7.theta))])
ax.plot([S2[0], S3[0]], [S2[1], S3[1]], **LS)                       # tratto a θ7
th_osc = np.linspace(0, 1, 60)                                      # oscillazione residua (qualitativa)
xo = S3[0] + th_osc * (xmax - S3[0])
yo = S3[1] + (xo - S3[0]) * tand(TH_W) + 0.12 * np.sin(np.pi * th_osc * 1.6) * np.exp(-2 * th_osc)
ax.plot(xo, yo, color="0.2", lw=1.0, ls=":")
ax.plot([0, P[0]], [0, 1], color="k", lw=2.0)                       # urto del cuneo
ax.plot([P[0], T[0]], [1, T[1]], color="k", lw=2.0)                 # urto trasmesso, tratto rettilineo
tt = np.linspace(0, 1, 30)                                          # oltre T: urto che si indebolisce
ang = r1.onda_sopra.sigma - 12 * tt
xx = T[0] + np.cumsum(np.cos(np.radians(ang))) * 0.07
yy = T[1] + np.cumsum(np.sin(np.radians(ang))) * 0.07
ax.plot(xx, yy, color="k", lw=1.6, ls="--")
for a_ in np.linspace(s3.theta - s3.mu, s5.theta - s5.mu, 5):       # ventaglio riflesso C−
    Q = incrocio(P, a_, O, TH_W)
    ax.plot([P[0], Q[0]], [1, Q[1]], color="tab:orange", lw=0.9)
for k in range(5):                                                  # ventaglio riflesso dalla parete (C+)
    A = Q1 + k / 4 * (Q2 - Q1)
    aa = (s5.theta + s5.mu) + k / 4 * ((s6.theta + s6.mu) - (s5.theta + s5.mu))
    B = incrocio(A, aa, S1, th_corda) if k else S1
    if k == 4:
        B = S2
    ax.plot([A[0], B[0]], [A[1], B[1]], color="tab:purple", lw=0.9 if k in (0, 4) else 0.6)
for B, a_, Lt in [(S1, s4.theta + s4.mu, None), (S2, s7.theta + s7.mu, 3.0)]:   # ventaglio trasmesso
    E = T if Lt is None else B + Lt * np.array([math.cos(math.radians(a_)), math.sin(math.radians(a_))])
    ax.plot([B[0], E[0]], [B[1], E[1]], color="tab:blue", lw=0.9)
for k in range(1, 4):
    B = S1 + k / 4 * (S2 - S1)
    a_ = (s4.theta + s4.mu) + k / 4 * ((s7.theta + s7.mu) - (s4.theta + s4.mu))
    E = incrocio(B, a_, P, r1.onda_sopra.sigma)
    if E[0] > xx[-1] or E[1] > 6.3:
        E = B + 2.5 * np.array([math.cos(math.radians(a_)), math.sin(math.radians(a_))])
    ax.plot([B[0], E[0]], [B[1], E[1]], color="tab:blue", lw=0.5)
for B in [S1 + k / 3 * (S2 - S1) for k in range(1, 4)]:            # compressioni riflesse (C−)
    a_ = s8.theta - s8.mu
    E = incrocio(B, a_, O, TH_W)
    ax.plot([B[0], E[0]], [B[1], E[1]], color="tab:green", lw=0.8)
for n, X, off in [("O", O, (-12, -12)), ("P", P, (-4, 8)), ("Q$_1$", Q1, (2, -14)), ("Q$_2$", Q2, (2, -14)),
                  ("S$_1$", S1, (-16, 5)), ("S$_2$", S2, (-4, 8)), ("T", T, (6, -4))]:
    ax.plot(*X, "k.", ms=5)
    ax.annotate(n, X, textcoords="offset points", xytext=off, fontsize=10)
txt = [("1", (-1.5, 2.0)), ("2", (-1.5, 0.45)), ("3", (1.9, 0.76)), ("4", (3.4, 2.75)), ("5", (3.3, 1.55)),
       ("6", (7.7, 2.95)), ("7", (11.6, 5.45)), ("8", (11.7, 4.72))]
for n, (x, y) in txt:
    ax.text(x, y, n, fontsize=10, ha="center", va="center",
            bbox=dict(boxstyle="circle,pad=0.15", fc="w", ec="0.4", lw=0.6))
ax.annotate("", xy=(-1.9, 2.6), xytext=(-2.9, 2.6), arrowprops=dict(arrowstyle="->"))
ax.text(-2.4, 2.75, "M$_1$ = 3", ha="center", fontsize=9)
ax.annotate("", xy=(-1.9, 0.75), xytext=(-2.9, 0.75), arrowprops=dict(arrowstyle="->"))
ax.text(-2.4, 0.18, "M$_2$ = 8", ha="center", fontsize=9)
ax.text(1.5, 0.17, "20°", fontsize=9)
ax.text(2.6, 4.3, "urto trasmesso", fontsize=8, rotation=49)
ax.text(0.55, 0.62, "urto del cuneo", fontsize=8, rotation=27)
ax.set_aspect("equal")
ax.set_xlim(-3.2, xmax)
ax.set_ylim(-0.4, 6.4)
ax.set_xlabel("x / h")
ax.set_ylabel("y / h")
ax.grid(alpha=0.15)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "campo.pdf"))
plt.close(fig)
print("\nFigure salvate in", FIG)
