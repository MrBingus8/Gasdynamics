"""Livello 3 di validazione: ricalcolo con gasdyn dei risultati numerici degli svolgimenti del corso.

Ogni :class:`Caso` contiene il valore riportato nello svolgimento (``fonte``), il valore
calcolato con il modulo, la tolleranza relativa adeguata alle cifre della fonte e l'esito
atteso: 'concorda', 'discorda' (sempre con la causa individuata in ``nota``) oppure
'non confrontabile'. I test (tests/test_gasdyn.py) verificano che l'esito calcolato
coincida con quello atteso, così una discrepanza nota resta marcata come tale e un
cambiamento del modulo che la facesse sparire (o ne creasse una nuova) viene segnalato.

Le pagine citate sono pagine fisiche di fonti/esercizi_2003_2023.pdf (prefisso "G p.")
oppure dei PDF in soluzioni_A/ ("A p.").
"""
from __future__ import annotations

import math
import pathlib
import sys
from dataclasses import dataclass

from scipy.optimize import brentq, minimize_scalar

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import gasdyn as gd  # noqa: E402

ARIA = gd.Gas(1.4)
ARIA_SI = gd.Gas(1.4, 287.0)


@dataclass
class Caso:
    id: str
    grandezza: str
    fonte: object
    calcolato: object
    tol: float
    atteso: str
    nota: str = ""
    pagine: str = ""

    @property
    def diff(self):
        if not isinstance(self.fonte, (int, float)) or not isinstance(self.calcolato, (int, float)):
            return None
        if self.fonte == 0:
            return abs(self.calcolato)
        return abs(self.calcolato - self.fonte) / abs(self.fonte)

    @property
    def esito(self):
        if self.atteso == "non confrontabile" or self.diff is None:
            return "non confrontabile"
        return "concorda" if self.diff <= self.tol else "discorda"


def C(id, grandezza, fonte, calcolato, tol, atteso="concorda", nota="", pagine=""):
    return Caso(id, grandezza, fonte, float(calcolato) if calcolato is not None else None, tol, atteso,
                nota, pagine)


# =============================================================================
# Prove 2021–2023 con soluzione impaginata (soluzioni_A)
# =============================================================================
def b61_e01():
    st0 = gd.State.da_pa(ARIA, 1.0, 1.0, u=4.0)
    st1 = gd.onda_semplice(st0, "u+a", u=0.0)
    V = gd.onda_verso(st0, "u+a", u=0.0)
    son = V.stato_xi(0.25 * 1.2 - 0.25 * 0.0)  # non usato: il punto sonico ha u = a
    a_son = (st0.a - 0.2 * st0.u) / (1 - 0.2)
    r = gd.contorno_serbatoio(st1, "destra", 1.0, a0=1.0)
    u, st2, st3 = r["onda"], r["stato_interno"], r["stato"]
    rif = gd.urto_riflesso_parete(st2, "u+a")
    pg = "A p.2–3 (G p.236–237)"
    nota_rif = ("lo svolgimento legge Ms = 2,51 dalle tabelle e lo combina con a2 ≈ 1,16 (w = 0,88); con lo "
                "stato 2 calcolato (u2 = −2,03, a2 = 1,10) le RH danno Ms = 2,60, w = 0,83. Prima "
                "approssimazione dichiarata ('errore inferiore al 2%')")
    return [
        C("B61-E01", "a1 dopo la chiusura", 0.2, st1.a, 1e-9, pagine=pg),
        C("B61-E01", "p1", 1.28e-5, st1.p, 1e-6, pagine=pg),
        C("B61-E01", "a = u nel punto sonico del ventaglio", 0.25, a_son, 1e-9, pagine=pg),
        C("B61-E01", "w urto riflesso dall'ambiente", -2.45, u.w, 0.01, pagine=pg),
        C("B61-E01", "u2", -2.03, st2.u, 0.01, pagine=pg),
        C("B61-E01", "p2", 2.24e-3, st2.p, 0.01, pagine=pg),
        C("B61-E01", "a2", 1.1, st2.a, 0.01, pagine=pg),
        C("B61-E01", "p3 (contorno, pressione totale = p_amb)", 2.27e-3, st3.p, 0.02,
          nota="lo svolgimento accetta p2 ≠ p3 entro il 2% (prima iterazione); il modulo impone p2 = p3",
          pagine=pg),
        C("B61-E01", "a3", 0.42, st3.a, 0.01, pagine=pg),
        C("B61-E01", "w urto riflesso dal setto", 0.88, rif.w, 0.02, "discorda", nota_rif, pg),
        C("B61-E01", "Ms urto riflesso dal setto", 2.51, rif.Ms, 0.02, "discorda", nota_rif, pg),
        C("B61-E01", "p4", 1.65e-2, rif.valle.p, 0.02, "discorda", nota_rif, pg),
        C("B61-E01", "a4", 1.63, rif.valle.a, 0.02, pagine=pg),
    ]


def b63():
    pg = "A p.2–3 (G p.242–244)"
    q = gd.State.da_pa(ARIA, 1.0, 1.0)
    u = gd.urto_mobile(q, "u+a", M_valle=1.0)
    u207 = gd.urto_mobile(q, "u+a", Ms=2.07)
    M_sep = brentq(lambda m: gd.theta_da_sigma(m, 28.0, ARIA) - 18.0, 2.2, 50)
    un10 = gd.urto_normale(10.0, ARIA)
    p_p1 = un10["p02_p1"] / gd.isentropiche(4.8, ARIA)["p0_p"]
    return [
        C("B63-E02", "Ms per M assoluto a valle = 1", 2.07, u.Ms, 0.005,
          nota="valore esatto 2,0681; lo svolgimento si ferma a 2,07 'dopo alcuni tentativi'", pagine=pg),
        C("B63-E02", "ρ2/ρ1 a Ms = 2,07", 2.7689, u207.rho2_rho1, 1e-4, pagine=pg),
        C("B63-E02", "a2/a1 a Ms = 2,07", 1.3211, u207.valle.a, 1e-4, pagine=pg),
        C("B63-E02", "M2 a Ms = 2,07", 1.001, u207.valle.u / u207.valle.a, 1e-3, pagine=pg),
        C("B63-E01", "M da Δθ = 18°, σ = 28°", 4.8, M_sep, 0.005, pagine=pg),
        C("B63-E01", "ν(4,8) [°]", 75.0, gd.prandtl_meyer(4.8, ARIA), 0.005, pagine=pg),
        C("B63-E01", "M da ν = 105°", 11.0, gd.mach_da_pm(105.0, ARIA), 0.01, pagine=pg),
        C("B63-E01", "p/p0 a M = 4,8", 0.24e-2, 1 / gd.isentropiche(4.8, ARIA)["p0_p"], 0.01, pagine=pg),
        C("B63-E01", "p02/p01 urto normale M = 10", 0.3e-2, un10["p02_p01"], 0.02, pagine=pg),
        C("B63-E01", "p01/p1 a M = 10", 4.35e4, gd.isentropiche(10.0, ARIA)["p0_p"], 0.01, "discorda",
          "valore di tabella arrotondato nello svolgimento (esatto 4,244·10⁴, −2,4%)", pg),
        C("B63-E01", "p/p1 nella bolla di ricircolo", 0.312, p_p1, 0.015,
          nota="concorda nonostante l'arrotondamento di p01/p1", pagine=pg),
    ]


def b67_e01():
    pg = "A p.3–4 (G p.251–252)"
    M = gd.mach_da_area(1.11, ARIA, "supersonico")
    Mc = brentq(lambda m: gd.theta_da_sigma(m, 51.0, ARIA) - 22.0, 1.3, 10)
    R = 82.5 / 60.0  # R ≈ 80–85 mm, semialtezza di gola 60 mm
    xf = -gd.sauer_meyer(R, ARIA, "favini")["xi_gola"](1.0) * 60
    xF = -gd.sauer_meyer(R, ARIA, "formulario")["xi_gola"](1.0) * 60
    return [
        C("B67-E01", "M medio con A/A* = 1,11", 1.39, M, 0.005, pagine=pg),
        C("B67-E01", "θmax(1,39) [°] (svolgimento: 9–10°)", 9.5, gd.theta_max(1.39, ARIA)[0], 0.06,
          nota="la fonte dà un intervallo 9–10°", pagine=pg),
        C("B67-E01", "M a monte dell'urto σ = 51°, θ = 22°", 2.2, Mc, 0.01, pagine=pg),
        C("B67-E01", "M da linea di Mach a 35°", 1.75, 1 / math.sin(math.radians(35)), 0.005, pagine=pg),
        C("B67-E01", "x_sonica Sauer–Meyer [mm], versione 'favini'", 8.5, xf, 0.06,
          nota="R = 82,5 mm, semialtezza 60 mm; distanza gola–punto sonico sull'asse = d/6 (scalata). "
               "Torna con la versione senza (γ+1) (S p.65 'Favini')", pagine=pg),
        C("B67-E01", "x_sonica Sauer–Meyer [mm], versione formulario", 8.5, xF, 0.06, "discorda",
          "discrepanza 10 di Fase 0: il formulario (con γ+1) dà 5,5 mm; lo svolgimento del 2022 usa "
          "la versione senza (γ+1)", pg),
    ]


def b68_e01():
    pg = "A p.3–4 (G p.256–257)"
    q = gd.State(ARIA, 1.0, 0.0, 1.0)  # p = 1, a = √γ
    out = []
    i = gd.urto_mobile(q, "u+a", M_valle=1.0)
    r = gd.urto_riflesso_parete(i.valle, "u-a")
    out += [C("B68-E01", "b1) w incidente (M valle = 1)", 2.45, i.w, 0.005, pagine=pg),
            C("B68-E01", "b1) w riflesso", -1.2, r.w, 0.005, pagine=pg),
            C("B68-E01", "b1) p/p1 al ristagno", 17.0, r.valle.p, 0.02, nota="valori arrotondati", pagine=pg)]
    try:
        gd.urto_mobile(q, "u+a", M_valle=1.89)
        lim = "risolto"
    except gd.NessunaSoluzione:
        lim = "nessuna soluzione"
    out.append(Caso("B68-E01", "b2) M valle = 1,89 (oltre il limite 1,8898)", "nessuna soluzione", lim, 0,
                    "non confrontabile",
                    "la fonte riconosce il caso limite e usa 1,87; il modulo solleva NessunaSoluzione", pg))
    i = gd.urto_mobile(q, "u+a", M_valle=1.87)
    r = gd.urto_riflesso_parete(i.valle, "u-a")
    nota = ("vicino al limite M valle è malcondizionato: lo svolgimento usa Ms ≈ 20 dalle tabelle "
            "(M valle = 1,8738 ≈ 1,87); il valore esatto per 1,87 è Ms = 17,99 (w = 21,3). "
            "Con Ms = 20 il modulo riproduce i valori della fonte (vedi riga seguente)")
    i20 = gd.urto_mobile(q, "u+a", Ms=20.0)
    r20 = gd.urto_riflesso_parete(i20.valle, "u-a")
    out += [C("B68-E01", "b2) w incidente (M valle = 1,87)", 23.66, i.w, 0.02, "discorda", nota, pg),
            C("B68-E01", "b2) w riflesso", -8.0, r.w, 0.05, "discorda", nota, pg),
            C("B68-E01", "b2) p/p1", 3700.0, r.valle.p, 0.05, "discorda", nota, pg),
            C("B68-E01", "b2) w incidente con Ms = 20", 23.66, i20.w, 0.002, pagine=pg),
            C("B68-E01", "b2) w riflesso con Ms = 20", -8.0, r20.w, 0.03, pagine=pg),
            C("B68-E01", "b2) p/p1 con Ms = 20", 3700.0, r20.valle.p, 0.03, pagine=pg),
            C("B68-E01", "σ a θmax per M → ∞ [°]", 67.0, gd.theta_max(1e6, ARIA)[1], 0.015,
              nota="'circa 67°' (esatto 67,79°)", pagine=pg)]
    return out


def b69():
    pg1, pg2 = "A p.3 (G p.261)", "A p.2–3 (G p.262–263)"
    R = gd.State.da_pa(ARIA, 1.0, 1.0)
    L = gd.onda_semplice(R, "u+a", u=2.0)
    s = gd.riemann(L, R)
    P = gd.pistone_legge(lambda t: t * t, R, 1.0)
    return [
        C("B69-E01", "θ_D(M = 1,3) [°] (soglia 'circa 3,2')", 3.2, gd.theta_distacco(1.3, ARIA), 0.02,
          pagine=pg1),
        C("B69-E02", "t_BU = 1/(γ+1)", 1 / 2.4, P.t_bu, 1e-9, pagine=pg2),
        C("B69-E02", "x_BU = t_BU", 1 / 2.4, P.x_bu, 1e-9, pagine=pg2),
        C("B69-E02", "p1 stato compresso = γ^(2γ/(γ−1))", 1.4**7, L.p, 1e-12, pagine=pg2),
        C("B69-E02", "w urto (pseudo-Riemann)", 2.87, s.onda_dx.w, 0.002, pagine=pg2),
        C("B69-E02", "u2 = u3", 2.1, s.u_star, 0.005, pagine=pg2),
        C("B69-E02", "p2 = p3", 9.44, s.p_star, 0.005, nota="'roughly'", pagine=pg2),
        C("B69-E02", "w asintotico (pistone a 2)", 2.76, gd.pistone_impulsivo(R, 2.0).w, 0.002, pagine=pg2),
    ]


def b70_e01():
    pg = "A p.3–4 (G p.267–268)"
    G = gd.Gas(1.4, 288.0)
    s1 = gd.State.da_pT(G, 1e5, 293.0)
    i = gd.urto_mobile(s1, "u+a", Ms=4.0)
    s2 = i.valle
    s3 = gd.espansione_pm(s2, 44.0, "C+").valle
    s4 = gd.espansione_pm(s3, 17.0, "C+").valle
    wB = i.w / 2.5
    u5 = gd.urto_mobile(s1, "u+a", w=wB)
    m3 = (s3.V - wB) / s3.a
    m4 = (s4.V - wB) / s4.a
    nR = "lo svolgimento usa a1 = 343,7 m/s, cioè R = 288 J/(kg K) (con R = 287: 343,1 m/s)"
    return [
        C("B70-E01", "a1 [m/s]", 343.7, s1.a, 1e-3, nota=nR, pagine=pg),
        C("B70-E01", "w [m/s]", 1375.0, i.w, 1e-3, pagine=pg),
        C("B70-E01", "ρ2/ρ1", 4.5714, i.rho2_rho1, 1e-4, pagine=pg),
        C("B70-E01", "u2 [m/s]", 1074.0, s2.u, 1e-3, pagine=pg),
        C("B70-E01", "a2 [m/s]", 691.4, s2.a, 1e-3, pagine=pg),
        C("B70-E01", "M2", 1.55, s2.M, 0.005, pagine=pg),
        C("B70-E01", "p2 [Pa]", 18.5e5, s2.p, 1e-3, pagine=pg),
        C("B70-E01", "ν(M2) [°]", 13.58, gd.prandtl_meyer(s2.M, G), 0.005, "discorda",
          "ν(1,55) da tabella; esatto ν(1,5534) = 13,48° (−0,1°), effetto trascurabile sui risultati", pg),
        C("B70-E01", "M3", 3.43, s3.M, 0.005, pagine=pg),
        C("B70-E01", "p3 [Pa]", 1.059e5, s3.p, 0.005, pagine=pg),
        C("B70-E01", "M4", 4.74, s4.M, 0.005, pagine=pg),
        C("B70-E01", "p4 [Pa]", 0.2e5, s4.p, 0.08, nota="valore dato a una cifra (0,2·10⁵)", pagine=pg),
        C("B70-E01", "Ms piede B", 1.6, u5.Ms, 1e-9, pagine=pg),
        C("B70-E01", "p5 [Pa]", 2.82e5, u5.valle.p, 1e-3, pagine=pg),
        C("B70-E01", "M5 = u5/a5", 0.77, u5.valle.u / u5.valle.a, 0.02, "discorda",
          "le RH con Ms = 1,6 danno u5/a5 = 0,690; il '0,77' del manoscritto non è coerente con p5 = 2,82·10⁵ "
          "(che invece torna): probabile errore di trascrizione", pg),
        C("B70-E01", "a3 [m/s]", 458.0, s3.a, 0.005, pagine=pg),
        C("B70-E01", "|V3| [m/s]", 1571.0, s3.V, 0.005, pagine=pg),
        C("B70-E01", "M3 relativo a RC", 2.23, m3, 0.005, pagine=pg),
        C("B70-E01", "p dietro RC lato 3 [Pa]", 6e5, s3.p * gd.urto_normale(m3, G)["p2_p1"], 0.01, pagine=pg),
        C("B70-E01", "a4 [m/s]", 358.0, s4.a, 0.005, pagine=pg),
        C("B70-E01", "|V4| [m/s]", 1697.0, s4.V, 0.005, pagine=pg),
        C("B70-E01", "M4 relativo a RC", 3.2, m4, 0.01, pagine=pg),
        C("B70-E01", "p dietro RC lato 4 [Pa]", 3.2e5, s4.p * gd.urto_normale(m4, G)["p2_p1"], 0.05, "discorda",
          "con M = 3,22 e p4 = 0,187·10⁵ (o 0,2·10⁵) le RH danno 2,2–2,4·10⁵ Pa: il '3,2·10⁵' del manoscritto "
          "sembra ripetere il Mach relativo (3,2)", pg),
    ]


def b71():
    pg1, pg2 = "A p.2 (G p.273)", "A p.2 (G p.274)"
    s1 = gd.State.da_M(ARIA, 2.0, p=1.0, a=1.0)
    s2 = gd.urto_obliquo(s1, theta=11.31).valle
    s3 = gd.onda_verso_theta_2d(s2, 0.0, "C+").valle
    s4 = gd.onda_verso_p_2d(s3, 0.296, "C+").valle
    g = 1.4
    q = gd.State.da_pa(ARIA, 1.0, 1.0)
    xp = lambda t: 1 + (2 * (1 - t) - (g + 1) * (1 - t) ** (2 / (g + 1))) / (g - 1)
    P = gd.pistone_legge(xp, q, 0.5)
    st = P.stato_pistone(0.5)
    r = gd.urto_riflesso_parete(st, "u-a")
    x0 = xp(0.5)
    ts = brentq(lambda t: x0 + st.u * (t - 0.5) - (1 + r.w * (t - 1)), 0.5, 5)
    return [
        C("B71-E01", "M2", 1.59, s2.M, 0.003, pagine=pg1),
        C("B71-E01", "p2/p1", 1.8, s2.p, 0.015, nota="arrotondato a 2 cifre (esatto 1,824)", pagine=pg1),
        C("B71-E01", "M3", 1.98, s3.M, 0.003, pagine=pg1),
        C("B71-E01", "p3/p1", 0.99, s3.p, 0.015, nota="esatto 1,004", pagine=pg1),
        C("B71-E01", "M4", 2.76, s4.M, 0.005, pagine=pg1),
        C("B71-E01", "|θ4| [°]", 19.0, abs(s4.theta), 0.01, pagine=pg1),
        C("B71-E02", "x_p(t_v = 0,5)", 0.13, x0, 0.03, pagine=pg2),
        C("B71-E02", "u dopo la compressione", 0.61, st.u, 0.005, pagine=pg2),
        C("B71-E02", "p dopo la compressione", 2.22, st.p, 0.015, nota="esatto 2,245", pagine=pg2),
        C("B71-E02", "a dopo la compressione", 1.12, st.a, 0.005, pagine=pg2),
        C("B71-E02", "w urto riflesso dalla parete", -0.93, r.w, 0.01, pagine=pg2),
        C("B71-E02", "Mach relativo", 1.38, r.Ms, 0.005, pagine=pg2),
        C("B71-E02", "t_s arresto pistone", 1.36, ts, 0.005, pagine=pg2),
        C("B71-E02", "x_s", 0.66, x0 + st.u * (ts - 0.5), 0.005, pagine=pg2),
    ]


def b72_e02():
    pg = "A p.2 (G p.277)"
    q = gd.State.da_pa(ARIA, 1.0, 1.0)
    i = gd.urto_mobile(q, "u+a", M_valle=1.87)
    s2 = i.valle
    r = gd.riemann_area(s2, q, 1.5, 1.0)
    s5, sol = r["stato_5"], r["riemann_valle"]
    return [
        C("B72-E02", "w", 18.0, i.w, 0.01, pagine=pg),
        C("B72-E02", "u dietro l'urto", 15.0, s2.u, 0.01, pagine=pg),
        C("B72-E02", "p dietro l'urto", 380.0, s2.p, 0.01, pagine=pg),
        C("B72-E02", "ρ dietro l'urto", 8.0, s2.rho, 0.05, pagine=pg),
        C("B72-E02", "u nella sezione stretta", 12.0, s5.u, 0.05, "discorda",
          "la massa impone ρ5u5 = 1,5 ρ2u2 ⇒ u5 = 10,5 con ρ5 = 17,7; il '12' non è coerente con ρ5 ≈ 18 e p5 ≈ 1100 "
          "(che tornano) né con u* ≈ 16 del Riemann successivo (che torna): refuso", pg),
        C("B72-E02", "p nella sezione stretta", 1100.0, s5.p, 0.01, pagine=pg),
        C("B72-E02", "ρ nella sezione stretta", 18.0, s5.rho, 0.02, pagine=pg),
        Caso("B72-E02", "configurazione", "RCS", sol.configurazione, 0, "non confrontabile",
             "coincide (RCS)", pg),
        C("B72-E02", "u* del Riemann", 16.0, sol.u_star, 0.01, pagine=pg),
        C("B72-E02", "p* del Riemann", 440.0, sol.p_star, 0.01, pagine=pg),
        C("B72-E02", "ρ a sinistra del contatto", 9.0, sol.stati["2"].rho, 0.03, pagine=pg),
        C("B72-E02", "ρ a destra del contatto", 8.0, sol.stati["3"].rho, 0.04, pagine=pg),
    ]


def soluzioni_A_qualitative():
    nota = "svolgimento qualitativo o simbolico: nessun numero da ricalcolare"
    out = [Caso(i, "—", None, None, 0, "non confrontabile", nota, "A")
           for i in ("B61-E02", "B67-E02", "B68-E02", "B69-E03", "B72-E01")]
    g = 1.4
    K = lambda M: (g + 1) / 2 * M ** (4 / (g + 1)) + M ** (-2 * (g - 1) / (g + 1))
    r = minimize_scalar(K, bounds=(0.05, 0.99), method="bounded", options={"xatol": 1e-12})
    out.append(C("B70-E02", "M stazionario di K(M) (escluso nel manoscritto)", math.sqrt((g - 1) / (g + 1)),
                 r.x, 1e-6, nota="controllo della derivazione simbolica", pagine="A p.2 (G p.270)"))
    return out


# =============================================================================
# Raccolta 2003–2020 (fonti/esercizi_2003_2023.pdf)
# =============================================================================
def b01_e02():
    pg = "G p.2–3"
    Mu = 2.5
    thD = gd.prandtl_meyer(Mu, ARIA) / 2
    hR = gd.area_su_area_critica(Mu, ARIA)
    L_semplice = ((hR - 1) / 2 + 1) / math.tan(math.radians(gd.angolo_mach(Mu)))
    moc = gd.ugello_minima_lunghezza(Mu, ARIA, n=40)
    return [
        C("B01-E02", "θ_D = ν(2,5)/2 [°]", 19.56, thD, 5e-4, pagine=pg),
        C("B01-E02", "h_R/h_L = A/A*(2,5)", 2.63672, hR, 1e-5, pagine=pg),
        C("B01-E02", "L/h_L (modello dello svolgimento)", 4.16638, L_semplice, 1e-5,
          nota="formula dello svolgimento: caratteristica d'uscita rettilinea inclinata di μ_u", pagine=pg),
        C("B01-E02", "L/h_L con il metodo delle caratteristiche (n = 40)", 4.16638, moc["L_su_h_gola"] / 2, 0.02,
          "discorda",
          "ipotesi diversa: lo svolgimento assume caratteristiche rette anche nella zona a due famiglie; il "
          "MOC del modulo (ugello_minima_lunghezza, semialtezza → altezza) dà L/h ≈ 4,58", pg),
    ]


def b02_e01():
    pg = "G p.5–6"
    G = gd.Gas(1.3, 1200 * 0.3 / 1.3)
    Me = gd.mach_da_area(3.0, G, "supersonico")
    st_e = gd.State.da_totali(G, Me, p0=3e6, T0=1690.0)
    st_t = gd.State.da_totali(G, 1.0, p0=3e6, T0=1690.0)
    mdot = st_t.rho * st_t.V
    Fcd = gd.spinta(mdot, st_e.V, st_e.p, 101.3e3, 3.0)
    Fc = gd.spinta(mdot, st_t.V, st_t.p, 101.3e3, 1.0)
    return [
        C("B02-E01", "Me", 2.5141, Me, 1e-4, pagine=pg),
        C("B02-E01", "p* [Pa]", 16.37183e5, st_t.p, 1e-5, pagine=pg),
        C("B02-E01", "F_C/A_g [N/m²]", 3664221.36332, Fc, 1e-6, pagine=pg),
        C("B02-E01", "ΔF = (F_CD − F_C)/F_CD [%]", 14.93635, (Fcd - Fc) / Fcd * 100, 1e-4,
          nota="la riga dati dello svolgimento scrive p0 = 1,013·10⁵ Pa ma i conti usano 3 MPa; "
               "denominatore F_CD come nello svolgimento", pagine=pg),
    ]


def b03_e01():
    pg = "G p.13"
    V, T1, p1, p2 = 0.25, 293.0, 1e5, 0.5e5
    g, R = 1.4, 287.0
    m1 = p1 * V / (R * T1)
    m2_fonte = m1 * (g - 1) / g + V * p2 / (g * R * T1)
    m2_form = m1 * (p2 / p1) ** (1 / g)
    m2_atm = gd.ATM * V / (R * T1) * 0.5 ** (1 / g)
    nota = ("lo svolgimento usa h_uscita = cp T1 costante (forma del riempimento) e 1 atm = 10⁵ Pa; il "
            "formulario (F §17, svuotamento adiabatico) dà m/m_i = (p/p_i)^(1/γ): 0,1810 kg con 10⁵ Pa, "
            "0,1834 kg con 1 atm = 101325 Pa")
    return [
        C("B03-E01", "m2 [kg] (modello dello svolgimento)", 0.1911, m2_fonte, 1e-3, pagine=pg),
        C("B03-E01", "m2 [kg] (F §17, p1 = 10⁵ Pa)", 0.1911, m2_form, 0.01, "discorda", nota, pg),
        C("B03-E01", "m2 [kg] (F §17, p1 = 1 atm)", 0.1911, m2_atm, 0.01, "discorda", nota, pg),
    ]


def b04_e01():
    pg = "G p.17–18"
    s = gd.riemann(gd.State(ARIA, 1, 0, 1), gd.State(ARIA, 0.125, 0, 0.1))
    k, u0, p2_0, p3_0 = s.iterazioni[0]
    st5 = gd.contorno_pressione(s.stati["3"], "destra", 0.1).valle
    r = gd.riemann(s.stati["2"], st5)
    nota = ("lo svolgimento scrive a5 = a3 (p5/p3)^(δ/γ) = 1,102, ma la stessa formula dà 1,0789; anche "
            "u5 = u3 + (a3−a5)/δ = 1,876 non è coerente con a5 = 1,102. L'errore si propaga a u6 (−0,5%) e p6 (−1,5%)")
    return [
        C("B04-E01", "u⁰ (innesco del formulario)", 0.918948, u0, 2e-6, pagine=pg),
        C("B04-E01", "p2 alla prima iterata", 0.306767, p2_0, 1e-5, pagine=pg),
        C("B04-E01", "p3 alla prima iterata", 0.30041, p3_0, 1e-4, pagine=pg),
        C("B04-E01", "u* a convergenza", 0.927454, s.u_star, 1e-5, pagine=pg),
        C("B04-E01", "p* a convergenza", 0.303131, s.p_star, 1e-5, pagine=pg),
        C("B04-E01", "w urto", 1.75216, s.onda_dx.w, 1e-5, pagine=pg),
        C("B04-E01", "a2", 0.997729, s.stati["2"].a, 1e-5, pagine=pg),
        C("B04-E01", "a3", 1.26411, s.stati["3"].a, 1e-5, pagine=pg),
        C("B04-E01", "a5 (riflessione a p = p4)", 1.102, st5.a, 0.005, "discorda", nota, pg),
        C("B04-E01", "u5", 1.876, st5.u, 0.005, "discorda", nota, pg),
        C("B04-E01", "u6 = u7 (velocità del contatto)", 1.753, r.u_star, 0.003, "discorda", nota, pg),
        C("B04-E01", "a6", 0.833, r.stati["2"].a, 0.003, pagine=pg),
        C("B04-E01", "p6", 0.0854, r.p_star, 0.005, "discorda", nota, pg),
    ]


def b05_e01():
    pg = "G p.21–23"
    s0 = gd.State.da_pa(ARIA, 1e5, 400.0)
    s2 = gd.pistone_impulsivo(s0, -300.0).valle
    S = gd.contorno_parete(s2, "sinistra")
    s3 = S.valle
    r = gd.riemann(s3, s0)
    rr = gd.urto_riflesso_parete(r.stati["2"], "u+a")
    return [
        C("B05-E01", "p2 (espansione) [Pa]", 0.32058e5, s2.p, 1e-4, pagine=pg),
        C("B05-E01", "p3 (urto d'arresto) [Pa]", 0.97831e5, s3.p, 1e-4, pagine=pg),
        Caso("B05-E01", "configurazione interazione", "SCR", r.configurazione, 0, "non confrontabile",
             "coincide (SCR)", pg),
        C("B05-E01", "u4 = u5 [m/s]", -3.15845, r.u_star, 0.002, nota="seconda iterazione dello svolgimento",
          pagine=pg),
        C("B05-E01", "p4 = p5 [Pa]", 0.98899e5, r.p_star, 1e-4, pagine=pg),
        C("B05-E01", "w2 [m/s]", -408.55, r.onda_sx.w, 1e-4, pagine=pg),
        C("B05-E01", "w3 urto riflesso dal pistone [m/s]", 405.73103, rr.w, 0.001,
          nota="la fonte usa a* da relazioni approssimate (diff. 0,07%)", pagine=pg),
        C("B05-E01", "p6 [Pa]", 0.9998e5, rr.valle.p, 1e-4, pagine=pg),
    ]


def b06_e02():
    pg = "G p.25–27"
    r = gd.shock_expansion([(0, 0), (1, 0)], [(0, 0), (1, 0)], 3.0, 10.0, ARIA)
    dorso, ventre = r["pannelli"]
    s = gd.State.da_M(ARIA, 3.0, p=1, a=1)
    o = gd.urto_obliquo(s, theta=10.0)
    return [
        C("B06-E02", "σ [°]", 27.38269, o.sigma, 1e-5, pagine=pg),
        C("B06-E02", "M2 (ventre)", 2.505, ventre["M"], 1e-3, pagine=pg),
        C("B06-E02", "p2/p1", 2.05447, ventre["p_pinf"], 1e-5, pagine=pg),
        C("B06-E02", "Cp ventre", 0.16738, ventre["Cp"], 1e-4, pagine=pg),
        C("B06-E02", "M3 (dorso)", 3.57829, dorso["M"], 1e-5, pagine=pg),
        C("B06-E02", "p3/p1", 0.43114, dorso["p_pinf"], 1e-4, pagine=pg),
        C("B06-E02", "Cp dorso", -0.09029, dorso["Cp"], 1e-3, pagine=pg),
        C("B06-E02", "CL", 0.25376, r["CL"], 1e-4, pagine=pg),
        C("B06-E02", "CD di forma", -0.04474, r["CD"], 1e-3, "discorda",
          "errore di segno nello svolgimento: CD_forma = (Cp_ventre − Cp_dorso) sinα = +0,04474; anche "
          "l'attrito proiettato va sommato con segno + (CD totale ≈ +0,0523, non −0,05232)", pg),
    ]


def b07_e01():
    pg = "G p.28–29"
    aria, elio = gd.Gas(1.4, 287.0), gd.Gas(1.6, 2078.0)
    s1 = gd.State.da_pT(aria, 1e5, 300.0)
    s5 = gd.State.da_pT(elio, 1e5, 300.0)
    u = gd.urto_mobile(s5, "u-a", p_valle=100e5)
    r = gd.riemann(s1, u.valle)
    return [
        C("B07-E01", "a5 elio [m/s]", 998.71918, s5.a, 1e-6, pagine=pg),
        C("B07-E01", "M_R5", 9.02427, u.Ms, 1e-5, pagine=pg),
        C("B07-E01", "M_R4", 0.44230, abs(u.Mr_valle), 1e-4, pagine=pg),
        C("B07-E01", "w [m/s]", -9012.71153, u.w, 1e-6, pagine=pg),
        C("B07-E01", "a4 [m/s]", 4894.89248, u.valle.a, 1e-6, pagine=pg),
        C("B07-E01", "u4 [m/s]", -6847.70058, u.valle.u, 1e-5, pagine=pg),
        C("B07-E01", "T4 [K]", 6126.45, u.valle.T, 0.01, "discorda",
          "lo svolgimento scrive T4/T5 = 20,4215 ma a4/a5 = 4,90117 (che torna) implica T4/T5 = 24,02: "
          "lettura di tabella errata; T4 corretto 7206 K", pg),
        C("B07-E01", "u* (prima approssimazione SCS)", -3849.07, r.u_star, 0.01, "discorda",
          "lo svolgimento si ferma alla prima approssimazione con γ medio = 1,5 e dichiara necessarie una "
          "decina di iterate; soluzione esatta SCS: u* = −4054 m/s, p* = 2,313·10⁷ Pa", pg),
    ]


def b08_e02():
    pg = "G p.33"
    s2 = gd.State.da_pa(ARIA, 1.0, 1.0, u=3.0)
    u = gd.urto_mobile(s2, "u+a", p_valle=4.5)
    return [
        C("B08-E02", "w", 5.0, u.w, 1e-9, pagine=pg),
        C("B08-E02", "M1 dietro l'urto", 3.271, u.valle.Mu, 1e-3, pagine=pg),
        C("B08-E02", "a1 dietro l'urto", 1.299, u.valle.a, 1e-3, pagine=pg),
        C("B08-E02", "A/A*(3,271)", 5.47669, gd.area_su_area_critica(3.271, ARIA), 1e-5, pagine=pg),
        C("B08-E02", "A/A*(3)", 4.23457, gd.area_su_area_critica(3.0, ARIA), 1e-5, pagine=pg),
        C("B08-E02", "A/A*(2,26)", 2.11535, gd.area_su_area_critica(2.26, ARIA), 1e-5, pagine=pg),
    ]


def b09_e01():
    pg = "G p.34"
    st = gd.State.da_M(ARIA, 3.0, p=1, a=1)
    o = gd.urto_obliquo(st, theta=34.0)
    dsC = gd.urto_normale(3.0, ARIA)["ds_R"]
    ts, ss = gd.sonico(3.0, ARIA)
    dsB_ex = gd.urto_normale(3.0 * math.sin(math.radians(ss)), ARIA)["ds_R"]
    return [
        C("B09-E01", "θmax(3) [°]", 34.07344, gd.theta_max(3.0, ARIA)[0], 1e-6, pagine=pg),
        C("B09-E01", "σ per θ = 34°", 63.67317, o.sigma, 1e-6, pagine=pg),
        C("B09-E01", "M_B", 1.00295, o.valle.M, 1e-5, pagine=pg),
        C("B09-E01", "(S_C − S∞)/R", 1.11369, dsC, 1e-5, pagine=pg),
        C("B09-E01", "(S_A − S_B)/R con θ = 34°", 0.264005, dsC + math.log(o.p02_p01), 2e-5,
          nota="la fonte approssima la deviazione sonica (34,008°) con 34°", pagine=pg),
        C("B09-E01", "(S_A − S_B)/R con il punto sonico esatto", 0.264005, dsC - dsB_ex, 2e-3, "discorda",
          "approssimazione dichiarata: con M_B = 1 esatto la differenza è 0,2622 (−0,7%)", pg),
    ]


def b36_e01():
    pg = "G p.89–90"
    g = gd.Gas(1.1)
    q = gd.State.da_pa(g, 1.0, 150.0)
    P = gd.pistone_legge(lambda t: (0.1 * t) ** 3, q, 1000.0, dx_p=lambda t: 0.003 * t * t,
                         ddx_p=lambda t: 0.006 * t)
    return [C("B36-E01", "t_BU [s]", 312.26, P.t_bu, 1e-5,
              nota=f"x_BU (lasciato in formula nella fonte) = {P.x_bu:.1f} m", pagine=pg)]


def b37_e01():
    pg = "G p.93–94"
    sf6, aria = gd.Gas(1.0984, 56.928), gd.Gas(1.4, 287.0)
    s1 = gd.State.da_pT(sf6, 100 * gd.ATM, 373.15)
    s4 = gd.State.da_pT(aria, gd.ATM, 298.15)
    r = gd.riemann(s1, s4)
    nota = ("la relazione d'urto usata per tentativi è errata: (−w)(u3−w) = 2a0²/(γ+1) senza il termine δw² in a*² "
            "(relazione di Prandtl nel riferimento dell'urto); con le RH esatte e gas diversi (SF6 | aria a 25 °C) "
            "u* = 410 m/s, p* = 4,22 atm. Lo svolgimento assume anche T_amb = 25 °C, assente nella traccia")
    return [
        C("B37-E01", "a4 aria a 25 °C [m/s]", 346.0, s4.a, 1e-3, pagine=pg),
        C("B37-E01", "u2 = u3 [m/s] (440–445)", 442.5, r.u_star, 0.01, "discorda", nota, pg),
        C("B37-E01", "p2 = p3 [atm]", 3.2, r.p_star / gd.ATM, 0.02, "discorda", nota, pg),
    ]


def b38_e02():
    pg = "G p.96–97"
    s1 = gd.State(ARIA, 1.4, 0.0, 1.0)
    s2 = gd.onda_semplice(s1, "u-a", a=2 / 3)
    sol = gd.riemann(s1, s2)
    s3 = sol.stato_xi(0.0)
    rho, u, p, dt = gd.godunov_passo([1.4, 1.4, s2.rho], [0, 0, s2.u], [1, 1, s2.p], ARIA, 1.0, cfl=0.8)
    st = gd.State(ARIA, rho[1], u[1], p[1])
    return [
        C("B38-E02", "u2", 1.66667, s2.u, 1e-5, pagine=pg),
        C("B38-E02", "p2", 0.05858, s2.p, 1e-3, pagine=pg),
        C("B38-E02", "a3 = u3 (stato sonico all'interfaccia)", 0.83333, s3.a, 1e-5, pagine=pg),
        C("B38-E02", "p3", 0.2791, s3.p, 1e-4, pagine=pg),
        C("B38-E02", "Δt (CFL 0,8)", 0.343, dt, 1e-3, pagine=pg),
        C("B38-E02", "ρ^(n+1)", 1.23947, rho[1], 5e-4, nota="flusso di massa arrotondato nella fonte", pagine=pg),
        C("B38-E02", "variazione di s = p/ρ^γ [%]", -0.16, (st.s / s1.s - 1) * 100, 0.1, "discorda",
          "nel flusso d'energia lo svolgimento usa ρe = p/(γ−1) + ρu² (manca ½): flusso 1,13365 invece di "
          "0,97678. Con il flusso corretto s cresce del 2,5% (lo schema di Godunov aumenta l'entropia, S p.42): "
          "la conclusione qualitativa 'non conserva l'entropia' resta valida", pg),
    ]


def b40_e03():
    pg = "G p.100–101"
    st = gd.State.da_M(ARIA, 3.0, p=1, a=1)
    f = lambda th: gd.urto_obliquo(st, theta=th).valle.p - gd.onda_isentropica_2d(st, "C+", theta2=th).valle.p
    th = brentq(f, 30, 34.07)
    return [
        C("B40-E03", "ν(3) [°]", 49.76, gd.prandtl_meyer(3.0, ARIA), 1e-3, pagine=pg),
        C("B40-E03", "θ con p_urto = p_compressione [°]", 24.17, th, 0.01, "discorda",
          "lo svolgimento usa ν2 = ν1 + θ (segno d'espansione) per la compressione e inverte il rapporto di "
          "pressione; con ν2 = ν1 − θ la curva di compressione isentropica incontra la polare solo vicino a θmax, "
          "sul ramo debole, a θ = 33,99° (p/p1 = 8,26). La fonte si interrompe ('Mi fermo qui')", pg),
    ]


def b41_e03():
    pg = "G p.104–109"
    s1 = gd.State.da_pa(ARIA, 1.0, 300.0)
    s2 = gd.State.da_pa(ARIA, 35.72, 500.0, u=1000.0)
    r = gd.riemann(s2, s1)
    rs = gd.riemann(s1, s2)
    return [
        C("B41-E03", "u6 = u7 [m/s]", 1118.20, r.u_star, 1e-4, pagine=pg),
        C("B41-E03", "w urto [m/s]", 1405.0, r.onda_dx.w, 0.001, pagine=pg),
        C("B41-E03", "p6 = p7 [atm]", 25.44, r.p_star, 0.001, pagine=pg),
        C("B41-E03", "u* a x = L/3 (onde quasi nulle)", 0.0, rs.u_star, 0.01,
          nota="il primo Riemann (1 | 2) dà onde di intensità ~1e-5: p* = 0,99997 atm", pagine=pg),
    ]


def b42_e02():
    pg = "G p.112"
    g = 1 / (1 - math.log(400 / 300) / math.log(2))
    R = (473 * (1.71 - 1) / 2) ** 2 / (1.71 * (math.sqrt(400) - math.sqrt(300)) ** 2)
    sup = gd.State.da_pT(gd.Gas(1.71, R), 1e5, 300.0)
    valle = gd.onda_semplice(sup, "u+a", u=473.0)
    return [
        C("B42-E02", "γ (ipotesi isentropica)", 1.71, g, 0.001, pagine=pg),
        C("B42-E02", "R [J/(kg K)] con γ = 1,71", 2296.32, R, 1e-3, pagine=pg),
        C("B42-E02", "massa molare [g/mol]", 3.62, 8314.46 / R, 0.002, pagine=pg),
        C("B42-E02", "T a valle dell'onda semplice [K]", 400.0, valle.T, 0.002, pagine=pg),
    ]


def b43_e01():
    pg = "G p.114–115"
    return [C("B43-E01", f"p0/p a M = 0,3 ({m})", v, gd.pressione_totale(0.3, ARIA, m), 6e-4, pagine=pg)
            for m, v in (("incomprimibile", 1.063), ("isentropico", 1.064), ("isotermo", 1.065))]


def b46_e02():
    pg = "G p.121–123, 127"
    s1 = gd.State(ARIA, 0.40, 8.25, 0.125)
    u = gd.urto_mobile(s1, "u-a", w=6.5)
    s2 = u.valle
    s3 = gd.variazione_area(s2, 1.0, 0.5)
    return [
        C("B46-E02", "a1", 0.661438, s1.a, 1e-6, pagine=pg),
        C("B46-E02", "M1", 12.47, s1.M, 0.001, pagine=pg),
        C("B46-E02", "M_r1 = (u1 − w)/a1", 2.645751, u.Mr_monte, 1e-6, pagine=pg),
        C("B46-E02", "p2", 1.0, s2.p, 1e-9, pagine=pg),
        C("B46-E02", "u2 (seconda versione, p.127)", 7.0, s2.u, 1e-9, pagine=pg),
        C("B46-E02", "u2 (prima versione, p.121)", 2.4286, s2.u, 0.01, "discorda",
          "la prima versione applica u2/u1 = ρ1/ρ2 alle velocità assolute invece che a quelle relative; "
          "la seconda (u2 = 7) è corretta", pg),
        C("B46-E02", "a2", 1.0, s2.a, 1e-9, pagine=pg),
        C("B46-E02", "A/A*(7)", 104.14, gd.area_su_area_critica(7.0, ARIA), 1e-4, pagine=pg),
        C("B46-E02", "p02", 4139.83607, s2.p0, 1e-6, pagine=pg),
        C("B46-E02", "M3 dopo la contrazione", 5.95, s3.M, 0.005, nota="Mach letto da tabella", pagine=pg),
        C("B46-E02", "u3", 6.879186, s3.u, 0.001, pagine=pg),
        C("B46-E02", "a3", 1.151166, s3.a, 0.002, pagine=pg),
        C("B46-E02", "p3", 2.764376, s3.p, 0.01, "discorda",
          "con M3 = 5,95 (tabella) invece di 5,970 la pressione cambia del 2,2%; gli altri valori tornano", pg),
    ]


def b48_e02():
    pg = "G p.135–136"
    a = gd.State(ARIA, 3.8125, 3.9161, 10.0)
    b = gd.State(ARIA, 1.0, 6.4928, 1.0)
    r = gd.natura_salto(a, b)
    return [C("B48-E02", "w urto (famiglia u−a)", 2.99, r["w"], 0.005,
              nota=f"natura riconosciuta dal modulo: {r['tipo']}", pagine=pg)]


def b49_e01():
    pg = "G p.138–141"
    out = []
    for g, nome in ((1.4, "aria"), (1.66, "elio"), (1.013, "SF6")):
        gas = gd.Gas(g)
        q = gd.State(gas, 1.0, 0.0, 1.0)
        P = gd.pistone_legge(lambda t: 0.74 * t * t, q, 1.0)
        if nome == "aria":
            out += [C("B49-E01", "t_BU aria", 0.66622, P.t_bu, 1e-5, pagine=pg),
                    C("B49-E01", "x_BU aria", 0.78827, P.x_bu, 1e-4, nota="a0 = 1,1832 arrotondato nella fonte", pagine=pg)]
            stp = gd.onda_semplice(q, "u+a", u=1.48)
            u_p = gd.urto_mobile(q, "u+a", p_valle=stp.p)
            out += [C("B49-E01", "p_compressione/p1 a u = 1,48", 4.77288, stp.p, 1e-4, pagine=pg),
                    C("B49-E01", "w finale (urto con p = p_compressione)", 2.4374, u_p.w, 0.002,
                      nota="ipotesi dello svolgimento p_urto = p_compressione; M letto da tabella 2,06 "
                           "(esatto 2,0577 ⇒ w = 2,4346)", pagine=pg),
                    C("B49-E01", "w finale asintotico (pistone a 1,48)", 2.4374, gd.pistone_impulsivo(q, 1.48).w,
                      0.002, "discorda",
                      "regime asintotico come in B69-E02: urto da pistone impulsivo a v_p = 1,48, w = 2,367. "
                      "L'ipotesi p_urto = p_compressione dello svolgimento ignora l'onda riflessa (pseudo-Riemann)",
                      pg)]
        else:
            out.append(Caso("B49-E01", f"t_BU, x_BU {nome}", None, None, 0, "non confrontabile",
                            f"non calcolati nella fonte; modulo: t_BU = {P.t_bu:.5f}, x_BU = {P.x_bu:.5f}", pg))
    return out


def b51_e01():
    pg = "G p.149–152"
    G = gd.Gas(1.4, 287.0)
    inf = gd.State.da_M(G, 6.0, p=0.01, T=260.0)
    sup = gd.State.da_M(G, 4.0, p=0.50, T=350.0)
    r = gd.riemann_stazionario(sup, inf)
    ri = gd.riemann_stazionario(sup, inf, modello="isentropico")
    m = gd.miscelamento([(inf, 1.0), (sup, 1.0)])
    return [
        C("B51-E01", "θ slip line [°]", -16.5, r.theta, 0.01, pagine=pg),
        C("B51-E01", "p comune [atm]", 0.07, r.p, 0.005, nota="la fonte: 0,0699–0,0703 a θ = −16,5°", pagine=pg),
        C("B51-E01", "M corrente superiore espansa", 5.63, r.onda_sopra.valle.M, 0.005, pagine=pg),
        C("B51-E01", "σ urto nella corrente inferiore [°]", 24.3, r.onda_sotto.sigma, 0.01, pagine=pg),
        C("B51-E01", "θ modello isentropico [°] (fra −16 e −17)", -16.5, ri.theta, 0.04, pagine=pg),
        Caso("B51-E01", "stato miscelato (bilanci globali)", None, None, 0, "non confrontabile",
             f"svolgimento solo impostato (p.152); modulo: ramo supersonico M = {m['supersonico'].M:.4f}, "
             f"p = {m['supersonico'].p:.4f} atm; ramo subsonico M = {m['subsonico'].M:.4f}, "
             f"p = {m['subsonico'].p:.4f} atm", pg),
    ]


def b52_e02():
    pg = "G p.161–163"
    q = gd.State.da_pa(ARIA, 1.0, 1.0)
    r = minimize_scalar(lambda uv: gd.onda_semplice(q, "u+a", u=uv).a0, bounds=(-3, 3), method="bounded",
                        options={"xatol": 1e-12})
    return [C("B52-E02", "u_v del minimo di a0,v (a_m = 1)", -1 / 1.2, r.x, 1e-8,
              nota="prima versione; la seconda versione ha il segno opposto (errore)", pagine=pg),
            C("B52-E02", "a0,v minimo = a_m/√(1+δ)", 1 / math.sqrt(1.2), r.fun, 1e-10, pagine=pg)]


def b53():
    out = []
    pg = "G p.165–167"
    s1 = gd.State.da_M(ARIA, 10.0, p=1.0, a=1.0)
    s2 = gd.urto_obliquo(s1, theta=10.0).valle
    s3 = gd.onda_isentropica_2d(s2, "C+", theta2=30.0).valle
    r = [x for x in gd.intersezioni_polari(s1, "C+", s3, "C-") if x["onda_A"].ramo == "debole"][0]
    out += [C("B53-E01", "M2", 6.65, s2.M, 0.002, pagine=pg),
            C("B53-E01", "p2 [atm]", 7.0868, s2.p, 0.002, pagine=pg),
            C("B53-E01", "M3 dopo la compressione a 30°", 4.32, s3.M, 0.005, "discorda",
              "valori 'provvisori' da tabelle: esatto M3 = 4,262 (−1,3%)", pg),
            C("B53-E01", "θ dell'interazione [°] (34–35)", 34.5, r["theta"], 0.01, pagine=pg),
            C("B53-E01", "p4 = p5 [atm]", 57.0, r["p"], 0.03, nota="valore provvisorio della fonte", pagine=pg),
            C("B53-E01", "M5", 4.63, r["onda_B"].valle.M, 0.015, pagine=pg)]
    pg = "G p.175–176"
    s1 = gd.State.da_pa(ARIA, 1.0, 1.0, u=7.0)
    s2 = gd.variazione_area(s1, 1.0, 2.0)
    u = gd.urto_mobile(s1, "u+a", w=9.0)
    s3 = u.valle
    rp = gd.riemann(s3, s2)
    ra = gd.riemann_area(s3, s2, 1.0, 2.0)
    out += [C("B53-E05", "M2 (sezione larga)", 8.15, s2.M, 0.002, pagine=pg),
            C("B53-E05", "u2", 7.08681539, s2.u, 1e-3, pagine=pg),
            C("B53-E05", "p2", 0.376139295, s2.p, 0.015, nota="M2 da tabella (8,15 invece di 8,161)", pagine=pg),
            C("B53-E05", "p3 dietro l'urto", 4.5, s3.p, 1e-9, pagine=pg),
            C("B53-E05", "u3", 8.25, s3.u, 1e-9, pagine=pg),
            C("B53-E05", "a3", 1.299038106, s3.a, 1e-8, pagine=pg),
            C("B53-E05", "u* Riemann (3 | 2)", 8.68, rp.u_star, 0.002, pagine=pg),
            C("B53-E05", "p* Riemann (3 | 2)", 2.82, rp.p_star, 0.002, pagine=pg),
            C("B53-E05", "p* con la variazione d'area applicata allo stato 3", 2.82, ra["riemann_valle"].p_star, 0.01,
              "discorda",
              "ipotesi diversa: lo svolgimento risolve il Riemann fra lo stato 3 (sezione stretta) e lo stato 2 "
              "senza far passare 3 attraverso l'allargamento; con il collegamento stazionario (come in B72-E02) "
              "la configurazione diventa SCS con p* = 1,88", pg)]
    pg = "G p.181"
    nota_M = ("il catalogo trascrive M1 = 5,45, ma il manoscritto ha u1 = 5,95 con a1 = 1 e tutti i conti "
              "(p01 = 1499,8; ν = 84,61°) corrispondono a M1 = 5,95: la cifra '4' è un '9'. Deviazione usata: 19°")
    s1 = gd.State.da_M(ARIA, 5.95, p=1.0, a=1.0)
    o2 = gd.espansione_pm(s1, 19.0, "C+")
    o3 = gd.onda_verso_theta_2d(o2.valle, 0.0, "C+")
    out += [C("B53-E08", "ν(M1) [°] con M1 = 5,95", 84.61, gd.prandtl_meyer(5.95, ARIA), 2e-4, nota=nota_M,
              pagine=pg),
            C("B53-E08", "p01", 1499.801629, s1.p0, 1e-6, pagine=pg),
            C("B53-E08", "M2 dopo l'espansione di 19°", 10.50, o2.valle.M, 0.002, pagine=pg),
            C("B53-E08", "p2", 0.025508534, o2.valle.p, 0.01, nota="M2 arrotondato a 10,50", pagine=pg),
            C("B53-E08", "a2", 0.540844, o2.valle.a, 0.003, "discorda",
              "refuso nel manoscritto: lo stesso foglio dà u2 = a2·M2 = 6,2117 = 10,50 × 0,5916, coerente con "
              "a2 = 0,5916 calcolato dal modulo", pg),
            C("B53-E08", "u2 = a2 M2", 6.21165, o2.valle.V, 0.003, pagine=pg),
            C("B53-E08", "σ urto di ritorno a θ = 0 [°]", 24.47, o3.sigma, 0.003, pagine=pg),
            C("B53-E08", "p3/p2", 21.910, o3.valle.p / o2.valle.p, 0.005, pagine=pg),
            C("B53-E08", "p3", 0.559, o3.valle.p, 0.01, pagine=pg)]
    pg = "G p.204–207"
    q = gd.State.da_pa(ARIA, 1.0, 1.0)
    s2 = gd.onda_semplice(q, "u+a", u=-2.0)
    S = gd.contorno_parete(s2, "sinistra")
    out += [C("B53-E14", "a2 dietro l'espansione", 0.6, s2.a, 1e-9, pagine=pg),
            C("B53-E14", "w urto d'arresto", 0.541661, S.w, 5e-5, pagine=pg),
            C("B53-E14", "M relativo", -4.2360683, S.Mr_monte, 1e-6, pagine=pg),
            C("B53-E14", "p3", 0.581380063, S.valle.p, 1e-6, pagine=pg),
            C("B53-E14", "a3", 1.2623, S.valle.a, 1e-4, pagine=pg)]
    return out


def b54_e02():
    pg = "G p.210–212"
    B = gd.State.da_M(ARIA, 1.0, p=1.0, rho=1.4)
    Ast = gd.State.da_M(ARIA, 2.0, p=1.0, a=0.5)
    Cst = gd.variazione_area(B, 1.0, 53.18, "supersonico")
    D = gd.urto_normale_stato(Cst)
    Cp = gd.variazione_area(Ast, 1.0, 53.18)
    r = gd.riemann(Cp, D)
    return [
        C("B54-E02", "M_C", 6.0, Cst.M, 1e-4, pagine=pg),
        C("B54-E02", "p_C", 0.0011989, Cst.p, 1e-4, pagine=pg),
        C("B54-E02", "p_D", 0.0501538, D.p, 1e-5, pagine=pg),
        C("B54-E02", "a_D", 1.07795, D.a, 1e-4, pagine=pg),
        C("B54-E02", "u_D", 0.4356810, D.u, 1e-5, pagine=pg),
        C("B54-E02", "M_C'", 6.75, Cp.M, 0.005, nota="Mach da tabella (esatto 6,769)", pagine=pg),
        C("B54-E02", "p_C'", 0.00238223, Cp.p, 0.01, "discorda",
          "conseguenza del Mach di tabella 6,75 invece di 6,769 (−1,8%)", pg),
        Caso("B54-E02", "configurazione", "SCS", r.configurazione, 0, "non confrontabile", "coincide (SCS)", pg),
        C("B54-E02", "p comune", 0.063, r.p_star, 0.002, pagine=pg),
        C("B54-E02", "u comune", 0.6149, r.u_star, 0.002, pagine=pg),
    ]


def tutti():
    gruppi = [b61_e01, b63, b67_e01, b68_e01, b69, b70_e01, b71, b72_e02, soluzioni_A_qualitative,
              b01_e02, b02_e01, b03_e01, b04_e01, b05_e01, b06_e02, b07_e01, b08_e02, b09_e01, b36_e01,
              b37_e01, b38_e02, b40_e03, b41_e03, b42_e02, b43_e01, b46_e02, b48_e02, b49_e01, b51_e01,
              b52_e02, b53, b54_e02]
    out = []
    for f in gruppi:
        out += f()
    return out


if __name__ == "__main__":
    for c in tutti():
        flag = "" if c.esito == c.atteso else "   <-- ESITO INATTESO"
        d = "" if c.diff is None else f"{c.diff:.2e}"
        print(f"{c.id:8s} {c.grandezza[:55]:55s} fonte={c.fonte!s:>14} calc={c.calcolato!s:>22.10} "
              f"diff={d:>9} {c.esito}{flag}")
