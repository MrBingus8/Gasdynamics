"""Esame 14/09/2026, esercizio 1: storia della pressione sulla parete di testa di un tubo
chiuso a sinistra che si svuota da un orifizio A2/A1 = 0,506 dopo l'apertura del setto.

Riproduce tutti i numeri e le figure dello svolgimento.
Uso:  python calcoli.py   (dalla cartella dello svolgimento)

Parti fuori modulo (scritte qui riusando gasdyn): reticolo delle caratteristiche con le
condizioni di parete e di uscita strozzata; schema ai volumi finiti quasi-1D usato solo come
controllo indipendente del reticolo.
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
ATM = gd.ATM

# =============================================================================
# 1. DATI (aria: gamma = 1,4, R = 287 J/(kg K), ipotesi: la traccia non dà il gas)
# =============================================================================
aria = gd.Gas(1.4, 287.0)
g, d = aria.gamma, aria.delta
L = 10.0                  # m
A2_A1 = 0.506
p0, T0 = 10 * ATM, 293.15
pa = 1 * ATM

s0 = gd.State.da_pT(aria, p0, T0, nome="0")
a0 = s0.a
print("=== Dati ===")
print(f"a0 = {a0:.2f} m/s, rho0 = {s0.rho:.3f} kg/m3, L/a0 = {1e3*L/a0:.2f} ms")

# =============================================================================
# 2. APERTURA DEL SETTO: orifizio strozzato, stato 1 nel tubo (F §16 legge delle aree, F §4)
# =============================================================================
crit = gd.critiche(aria)
pstar_p0 = crit["p*_p0"]
print("\n=== Critiche ===", crit)
# ipotesi: A2 sonica  =>  A1/A* = 1/0,506 e ramo subsonico nel tubo
Me = gd.mach_da_area(1 / A2_A1, aria, "subsonico")
a1 = a0 / (1 + d * Me)                      # R2 = a + delta u uniforme nell'onda u-a
s1 = gd.onda_semplice(s0, "u-a", a=a1, nome="1")
assert abs(s1.u / s1.a - Me) < 1e-10
p01 = s1.p0
pe1 = pstar_p0 * p01                        # pressione nella sezione A2 (sonica)
print("\n=== Stato 1 (dietro l'espansione incidente) ===")
print(f"A1/A* = {1/A2_A1:.4f}, Me = {Me:.4f}")
print(f"a1 = {s1.a:.2f} m/s, u2 = {s1.u:.2f} m/s, p2 = {s1.p/ATM:.4f} atm, T2 = {s1.T:.2f} K")
print(f"p01 = {p01/ATM:.4f} atm, p_A2 = p* = {pe1/ATM:.4f} atm  (> pa = 1 atm: strozzato, getto sottoespanso)")
print(f"fronti: testa u-a = {-a0:.1f} m/s, coda u2-a1 = {s1.u - s1.a:.1f} m/s")

# Controllo: con A2 sonica, la portata in A1 (stato 1) coincide con quella critica in A2
mdotA1 = s1.rho * s1.u                       # per unità di A1
mdotA2 = gd.portata_max(A2_A1, p01, s1.T0, aria)
print(f"portata per unità di A1: rho2 u2 = {mdotA1:.2f}, critica in A2 = {mdotA2:.2f} kg/(s m2)")

# =============================================================================
# 3. RIFLESSIONI: parete (u = 0) e uscita strozzata (u = Me a)  (F §9)
# =============================================================================
r = (1 - d * Me) / (1 + d * Me)             # a_parete,n+1 / a_parete,n
print("\n=== Riflessioni ===")
print(f"r = (1 - dM2)/(1 + dM2) = {r:.5f},  r^7 = {r**7:.5f}")
s2v = gd.contorno_parete(s1, "sinistra")    # parete a sinistra: onda u+a
s2 = s2v.valle
s2 = s2.con(nome="2")
print(f"stato 2 (parete): a2 = {s2.a:.2f} = a1 - d u1 = {s1.a - d*s1.u:.2f}, p2 = {s2.p/ATM:.4f} atm, tipo {type(s2v).__name__}")
# uscita: R2 = a3 dalla caratteristica u+a uscente, condizione u = Me a
a3 = s2.a / (1 + d * Me)
s3 = gd.onda_semplice(s2, "u-a", a=a3, nome="3")
assert abs(s3.u / s3.a - Me) < 1e-10
s4 = gd.contorno_parete(s3, "sinistra").valle.con(nome="4")
a5 = s4.a / (1 + d * Me)
s5 = gd.onda_semplice(s4, "u-a", a=a5, nome="5")
for s in (s0, s1, s2, s3, s4, s5):
    print(f"  stato {s.nome or '3/5'}: p = {s.p/ATM:.4f} atm, u = {s.u:7.2f} m/s, a = {s.a:.2f} m/s, "
          f"T = {s.T:.1f} K, M = {s.u/s.a:.3f}")
# soglia di strozzamento: p*(stato di uscita) = pa con M = Me nel tubo
p0_p_M2 = gd.isentropiche(Me, aria)["p0_p"]
p_soglia = pa / pstar_p0 / p0_p_M2
print(f"strozzamento finché p(uscita, lato tubo) >= {p_soglia/ATM:.4f} atm  "
      f"(p0 uscita >= {pa/pstar_p0/ATM:.4f} atm)")
print(f"p3 = {s3.p/ATM:.3f} (strozzato: {s3.p >= p_soglia}), p5 = {s5.p/ATM:.3f} (strozzato: {s5.p >= p_soglia})")


# =============================================================================
# 4. RETICOLO DELLE CARATTERISTICHE (fuori modulo, riusa gd.punto_caratteristiche come controllo)
#    Flusso omoentropico (solo espansioni): nel nodo (k, m) d'incontro della u-a k con la u+a m
#    lo stato è noto a priori, a = (R1_k + R2_m)/2, u = (R2_m - R1_k)/(2 delta);
#    la posizione si trova con le pendenze medie sui due lati del segmento.
#    u-a k = 0..N: ventaglio centrato in (L, 0), R1 da a0 a a3.
#    u+a m: riflessa a parete della u-a m, R2_m = R1_m.
#    u-a k > N: riflessa all'uscita della u+a m = k - N - 1, R1_k = r R2_m.
# =============================================================================
def reticolo(N, k_max):
    R1 = list(np.linspace(a0, s2.a, N + 1))
    while len(R1) <= k_max:
        R1.append(r * R1[len(R1) - N - 1])
    R1 = np.array(R1)
    R2 = R1.copy()
    nodi = {}

    def st(k, m):
        return 0.5 * (R1[k] + R2[m]), (R2[m] - R1[k]) / (2 * d)

    def incrocio(P1, P2, aP, uP):
        x1, t1, a1, u1 = P1        # sulla u-a
        x2, t2, a_2, u2 = P2       # sulla u+a
        s1 = 0.5 * ((u1 - a1) + (uP - aP))
        s2_ = 0.5 * ((u2 + a_2) + (uP + aP))
        t = (x2 - x1 + s1 * t1 - s2_ * t2) / (s1 - s2_)
        return x1 + s1 * (t - t1), t

    for k in range(k_max + 1):
        if k <= N:
            aC, uC = 0.5 * (R1[k] + a0), (a0 - R1[k]) / (2 * d)
            prec = (L, 0.0, aC, uC)          # centro del ventaglio
            m0 = 0
        else:
            m0 = k - N - 1                   # nodo d'uscita (k, m0)
            aP, uP = st(k, m0)
            x2, t2, a_2, u2 = nodi[(k - 1, m0)]
            s2_ = 0.5 * ((u2 + a_2) + (uP + aP))
            nodi[(k, m0)] = (L, t2 + (L - x2) / s2_, aP, uP)
            prec = nodi[(k, m0)]
            m0 += 1
        for m in range(m0, k):
            aP, uP = st(k, m)
            x, t = incrocio(prec, nodi[(k - 1, m)], aP, uP)
            nodi[(k, m)] = (x, t, aP, uP)
            prec = nodi[(k, m)]
        aP = R1[k]                           # parete: u = 0, a = R1
        x1, t1, a1, u1 = prec
        s1 = 0.5 * ((u1 - a1) + (-aP))
        nodi[(k, k)] = (0.0, t1 + (0.0 - x1) / s1, aP, 0.0)
    return nodi, R1


def p_da_a(a):
    return p0 * (a / a0) ** (2 * g / (g - 1))


def storie(N, ngen=4):
    k_max = ngen * (N + 1)
    nodi, R1 = reticolo(N, k_max)
    tw = np.array([0.0] + [nodi[(k, k)][1] for k in range(k_max + 1)])
    pw = np.array([p0] + [p_da_a(nodi[(k, k)][2]) for k in range(k_max + 1)])
    # uscita (lato tubo): stato 1 da t = 0+, poi nodi d'uscita
    te, pe = [0.0], [s1.p]
    for k in range(N + 1, k_max + 1):
        x, t, a, u = nodi[(k, k - N - 1)]
        te.append(t)
        pe.append(p_da_a(a))
    te, pe = np.array(te), np.array(pe)
    # istante di fine strozzamento: p(uscita) = p_soglia (interpolazione fra nodi)
    i = np.argmax(pe < p_soglia)
    t_u = te[i - 1] + (p_soglia - pe[i - 1]) * (te[i] - te[i - 1]) / (pe[i] - pe[i - 1])
    return dict(nodi=nodi, R1=R1, tw=tw, pw=pw, te=te, pe=pe, t_u=t_u, k_max=k_max)


def istanti(S, N):
    """Arrivi delle teste e delle code dei ventagli a parete e all'uscita."""
    n = S["nodi"]
    out = {}
    out["t1 testa incidente a parete"] = n[(0, 0)][1]
    out["t2 coda incidente a parete"] = n[(N, N)][1]
    out["t3 testa riflessa all'uscita"] = n[(N + 1, 0)][1]
    out["t4 coda riflessa all'uscita"] = n[(2 * N + 1, N)][1]
    out["t5 testa 2a a parete"] = n[(N + 1, N + 1)][1]
    out["t6 coda 2a a parete"] = n[(2 * N + 1, 2 * N + 1)][1]
    out["t7 testa 2a riflessa all'uscita"] = n[(2 * N + 2, N + 1)][1]
    out["t_u fine strozzamento"] = S["t_u"]
    out["t8 coda 2a riflessa all'uscita"] = n[(3 * N + 2, 2 * N + 1)][1]
    return out


print("\n=== Reticolo delle caratteristiche: convergenza degli istanti (ms) ===")
ris = {}
for N in (3, 10, 40, 160):
    S = storie(N)
    ris[N] = (S, istanti(S, N))
chiavi = list(ris[160][1].keys())
print(f"{'':34s}" + "".join(f"{'N=' + str(N):>10s}" for N in ris))
for c in chiavi:
    print(f"{c:34s}" + "".join(f"{1e3*ris[N][1][c]:10.2f}" for N in ris))
NF = 160
S = ris[NF][0]
# verifica: caratteristiche della stessa famiglia non si incrociano (niente urti), arrivi ordinati
assert np.all(np.diff(S["tw"]) >= 0) and np.all(np.diff(S["te"]) >= 0)
for k in range(1, S["k_max"]):
    xs = [S["nodi"][(k, m)][0] for m in range(max(0, k - NF - 1), k + 1)]
    assert np.all(np.diff(xs) < 0)
print("verifica reticolo: arrivi a parete e all'uscita ordinati, nessun incrocio fra caratteristiche della stessa famiglia")

IST = ris[NF][1]
t_u = S["t_u"]
pw_tu = np.interp(t_u, S["tw"], S["pw"])
print(f"p parete a t_u = {pw_tu/ATM:.4f} atm")

# Controllo con gd.punto_caratteristiche su un nodo interno (stato del nodo a priori)
n = S["nodi"]
k, m = NF // 2, 0
A = gd.isentropica(s0, a=S["R1"][m]).con(u=0.0)
B = gd.isentropica(s0, a=0.5 * (S["R1"][k] + a0)).con(u=(a0 - S["R1"][k]) / (2 * d))
P = gd.punto_caratteristiche(A, B)
print(f"controllo nodo ({k},{m}): a = {n[(k, m)][2]:.4f} / {P.a:.4f}, u = {n[(k, m)][3]:.4f} / {P.u:.4f}")

# Il reticolo grossolano N = 3 (a mano): stati e istanti della prima riflessione
S3 = ris[3][0]
print("\n=== Reticolo a mano, N = 3 (ventaglio diviso in 3 onde) ===")
for (k, m), (x, t, a, u) in sorted(S3["nodi"].items()):
    if k <= 4:
        print(f"  nodo ({k},{m}): x = {x:6.3f} m, t = {1e3*t:7.2f} ms, a = {a:7.2f}, u = {u:7.2f}, "
              f"p = {p_da_a(a)/ATM:.3f} atm")


# =============================================================================
# 5. CONTROLLO INDIPENDENTE: volumi finiti quasi-1D (HLL), parete a sinistra,
#    contorno d'uscita con R2 dall'interno e u = Me a (fuori modulo)
# =============================================================================
def volumi_finiti(nc=800, t_fine=0.14, cfl=0.8):
    dx = L / nc
    rho = np.full(nc, s0.rho); u = np.zeros(nc); p = np.full(nc, p0)

    def cons(rho, u, p):
        return np.array([rho, rho * u, p / (g - 1) + 0.5 * rho * u**2])

    def flux(rho, u, p):
        E = p / (g - 1) + 0.5 * rho * u**2
        return np.array([rho * u, rho * u**2 + p, u * (E + p)])

    U = cons(rho, u, p)
    t, storia = 0.0, [(0.0, p0)]
    while t < t_fine:
        rho = U[0]; u = U[1] / rho; p = (g - 1) * (U[2] - 0.5 * rho * u**2)
        a = np.sqrt(g * p / rho)
        dt = min(cfl * dx / np.max(np.abs(u) + a), t_fine - t)
        # ghost: parete (riflessione) a sinistra; a destra stato di contorno caratteristico
        R2b = a[-1] + d * u[-1]
        ab = R2b / (1 + d * Me); ub = Me * ab
        pb = p[-1] * (ab / a[-1]) ** (2 * g / (g - 1)); rhob = g * pb / ab**2
        rL = np.concatenate(([rho[0]], rho)); uL = np.concatenate(([-u[0]], u)); pL = np.concatenate(([p[0]], p))
        rR = np.concatenate((rho, [rhob])); uR = np.concatenate((u, [ub])); pR = np.concatenate((p, [pb]))
        aL = np.sqrt(g * pL / rL); aR = np.sqrt(g * pR / rR)
        SL = np.minimum(uL - aL, uR - aR); SR = np.maximum(uL + aL, uR + aR)
        FL = flux(rL, uL, pL); FR = flux(rR, uR, pR)
        UL = cons(rL, uL, pL); UR = cons(rR, uR, pR)
        Fh = (SR * FL - SL * FR + SL * SR * (UR - UL)) / (SR - SL)
        Fh = np.where(SL >= 0, FL, np.where(SR <= 0, FR, Fh))
        Fh[:, -1] = flux(rhob, ub, pb)        # uscita: flusso dello stato di contorno
        Fh[:, 0] = [0.0, pL[0], 0.0]          # parete
        U = U - dt / dx * (Fh[:, 1:] - Fh[:, :-1])
        t += dt
        pw_ = (g - 1) * (U[2, 0] - 0.5 * U[1, 0]**2 / U[0, 0])
        storia.append((t, pw_))
    return np.array(storia)


VF = volumi_finiti()
print("\n=== Controllo volumi finiti (800 celle) ===")
for tt in (0.020, 0.040, 0.060, 0.080, 0.100, 0.110, 0.120):
    print(f"  t = {1e3*tt:5.0f} ms: p parete reticolo = {np.interp(tt, S['tw'], S['pw'])/ATM:.3f}, "
          f"volumi finiti = {np.interp(tt, VF[:, 0], VF[:, 1])/ATM:.3f} atm")

# =============================================================================
# 6. CONFRONTO: modello a parametri concentrati (F §17), V = A1 L, orifizio A2
# =============================================================================
t_lump = np.linspace(0, 0.16, 400)
LU = gd.serbatoio_svuotamento(t_lump, p0, a0, L * 1.0, A2_A1, aria, p_a=pa)
print("\n=== Parametri concentrati ===")
print({k: v for k, v in LU.items() if not isinstance(v, np.ndarray)})
Cs = (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1)))
K = A2_A1 * a0 * Cs / L
print(f"C* = {Cs:.4f}, K = {K:.3f} 1/s, 1/K = {1e3/K:.1f} ms, 2L/a0 = {2e3*L/a0:.1f} ms")
t_lu = (((pa / pstar_p0) / p0) ** (-(g - 1) / (2 * g)) - 1) / (d * K)
print(f"fine strozzamento (parametri concentrati): p = {pa/pstar_p0/ATM:.3f} atm a t = {1e3*t_lu:.1f} ms")
p_lump_tu = p0 * (1 + d * K * t_u) ** (-2 * g / (g - 1))
print(f"p (parametri concentrati) a t_u = {p_lump_tu/ATM:.3f} atm")

# =============================================================================
# 7. FIGURE
# =============================================================================
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3})
ms = 1e3

# --- 7a. Piano x-t (reticolo con un sottoinsieme di caratteristiche)
Nd = 10
Sd = ris[Nd][0]
nd = Sd["nodi"]
fig, ax = plt.subplots(figsize=(6.4, 7.4))
km = Sd["k_max"]
for k in range(km + 1):
    pts = sorted([(m, v) for (kk, m), v in nd.items() if kk == k])
    xs = [v[0] for _, v in pts]; ts = [v[1] for _, v in pts]
    if k <= Nd:                              # parte dal centro del ventaglio (L, 0)
        xs = [L] + xs; ts = [0.0] + ts
    testa_coda = (k % (Nd + 1) == 0) or (k % (Nd + 1) == Nd)
    ax.plot(xs, np.array(ts) * ms, color="tab:blue", lw=1.4 if testa_coda else 0.5)
for m in range(km + 1):
    pts = sorted([(k, v) for (k, mm), v in nd.items() if mm == m])
    xs = [v[0] for _, v in pts]; ts = [v[1] for _, v in pts]
    testa_coda = (m % (Nd + 1) == 0) or (m % (Nd + 1) == Nd)
    ax.plot(xs, np.array(ts) * ms, color="tab:red", lw=1.4 if testa_coda else 0.5)
ax.axhline(t_u * ms, color="k", ls="--", lw=1)
ax.text(4.2, t_u * ms + 1.5, r"$t_u$: fine dello strozzamento", fontsize=9)
etich = [("0", 3.0, 12), ("1", 8.6, 20), ("2", 1.6, 70), ("3", 9.0, 95), ("4", 1.0, 131)]
for s_, x_, t_ in etich:
    ax.text(x_, t_, s_, fontsize=13, ha="center", va="center",
            bbox=dict(boxstyle="circle,pad=0.2", fc="w", ec="k", lw=0.8))
ax.set_xlim(0, L); ax.set_ylim(0, 145)
ax.set_xlabel("x [m]   (parete di testa a x = 0, orifizio A$_2$ a x = L)")
ax.set_ylabel("t [ms]")
ax.plot([], [], color="tab:blue", label="caratteristiche $u-a$ (onde $u-a$)")
ax.plot([], [], color="tab:red", label="caratteristiche $u+a$ (onde $u+a$)")
ax.legend(loc="lower left", bbox_to_anchor=(0.0, 0.33), fontsize=8, framealpha=0.95)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "xt.pdf"))
plt.close(fig)

# --- 7b. Storia della pressione alla parete (e all'uscita, per confronto)
fig, ax = plt.subplots(figsize=(7.0, 4.4))
tw, pw = S["tw"], S["pw"]
sel = tw <= t_u
ax.plot(np.append(tw[sel], t_u) * ms, np.append(pw[sel], pw_tu) / ATM, color="k", lw=2,
        label="parete di testa (caratteristiche)")
te, pe = S["te"], S["pe"]
sel = te <= t_u
pe_tu = p_soglia
ax.plot(np.append(te[sel], t_u) * ms, np.append(pe[sel], pe_tu) / ATM, color="tab:green", lw=1.2,
        ls="-.", label="uscita del tubo, a monte di A$_2$")
ax.plot(t_lump[t_lump <= t_lu] * ms, LU["p_pi"][t_lump <= t_lu] * p0 / ATM, color="tab:orange", lw=1.2, ls="--",
        label="parametri concentrati")
ax.axhline(p_soglia / ATM, color="gray", lw=0.8, ls=":")
ax.text(1, p_soglia / ATM + 0.15, f"soglia di strozzamento all'uscita, {p_soglia/ATM:.2f} atm", fontsize=8, color="gray")
ax.axvline(t_u * ms, color="k", lw=0.8, ls="--")
ax.text(t_u * ms + 0.8, 8.5, r"$t_u$", fontsize=10)
for nome, (tt, pp) in {"0": (12, p0), "2": (68, s2.p), "4": (129, s4.p)}.items():
    ax.text(tt, pp / ATM + 0.3, nome, fontsize=11, ha="center",
            bbox=dict(boxstyle="circle,pad=0.2", fc="w", ec="k", lw=0.8))
ax.set_xlim(0, 142); ax.set_ylim(0, 10.8)
ax.set_xlabel("t [ms]"); ax.set_ylabel("p [atm]")
ax.legend(loc="upper right", fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "p_t.pdf"))
plt.close(fig)

# --- 7c. Piano p-u: stati di plateau fra le rette di parete e di uscita
fig, ax = plt.subplots(figsize=(6.4, 4.6))
pp = np.linspace(1.3 * ATM, p0, 300)
aa = a0 * (pp / p0) ** ((g - 1) / (2 * g))
ax.plot(np.zeros_like(pp), pp / ATM, color="k", lw=1.2, label="parete: $u=0$")
ax.plot(Me * aa, pp / ATM, color="tab:green", lw=1.2, label=f"uscita strozzata: $u=M_2a$, $M_2={Me:.2f}$")
seq = [s0, s1, s2, s3, s4, s5]
nomi = ["0", "1", "2", "3", "4", "5"]
for i in range(len(seq) - 1):
    A_, B_ = seq[i], seq[i + 1]
    fam = "u-a" if i % 2 == 0 else "u+a"
    pc = np.linspace(A_.p, B_.p, 60)
    uc = [gd.onda_semplice(A_, fam, p=x).u for x in pc]
    ax.plot(uc, pc / ATM, color="tab:blue" if fam == "u-a" else "tab:red", lw=1.0,
            ls="-" if i < 4 else ":")
for s_, n_ in zip(seq, nomi):
    ax.plot(s_.u, s_.p / ATM, "ko", ms=4)
    ax.annotate(n_, (s_.u, s_.p / ATM), textcoords="offset points", xytext=(6, 3), fontsize=10)
ax.axhline(p_soglia / ATM, color="gray", lw=0.8, ls=":")
ax.text(20, p_soglia / ATM + 0.15, "soglia di strozzamento", fontsize=8, color="gray")
ax.plot([], [], color="tab:blue", label="espansioni $u-a$ ($a+\\delta u$ cost.)")
ax.plot([], [], color="tab:red", label="espansioni $u+a$ ($a-\\delta u$ cost.)")
ax.set_xlabel("u [m/s]"); ax.set_ylabel("p [atm]")
ax.set_xlim(-5, 115); ax.set_ylim(0, 10.8)
ax.legend(loc="upper right", fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "p_u.pdf"))
plt.close(fig)

# --- 7d. Schema del condotto (non in scala)
fig, ax = plt.subplots(figsize=(7, 2.2))
h1, h2 = 1.6, 1.6 * 0.506
ax.plot([0, L, L], [h1 / 2, h1 / 2, h2 / 2], "k", lw=2)
ax.plot([0, L, L], [-h1 / 2, -h1 / 2, -h2 / 2], "k", lw=2)
ax.plot([0, 0], [-h1 / 2, h1 / 2], "k", lw=4)
ax.plot([L, L], [-h2 / 2, h2 / 2], color="tab:red", ls="--", lw=1.5)
ax.text(L / 2, 0.15, "aria in quiete", ha="center", va="center")
ax.text(L / 2, -0.25, "$p_0=10$ atm, $T_0=293{,}15$ K", ha="center", va="center")
ax.text(0.2, -0.65, "parete di testa", fontsize=8)
ax.annotate("", (2.8, h1 / 2), (2.8, -h1 / 2), arrowprops=dict(arrowstyle="<->"))
ax.text(2.95, 0.3, "A$_1$", fontsize=9)
ax.text(L + 0.2, 0.0, "A$_2$ = 0,506 A$_1$ (setto)\nambiente $p_a=1$ atm", fontsize=8, va="center")
ax.annotate("", (L, -h1 / 2 - 0.15), (0, -h1 / 2 - 0.15), arrowprops=dict(arrowstyle="<->"))
ax.text(L / 2, -h1 / 2 - 0.38, "L = 10 m", ha="center", fontsize=8)
ax.set_xlim(-0.3, L + 3.8); ax.set_ylim(-1.35, 0.95); ax.axis("off")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "schema.pdf"))
plt.close(fig)
print("\nFigure salvate in", FIG)
