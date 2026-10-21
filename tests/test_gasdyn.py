"""Test di gasdyn.py.

Livello 1: coerenza interna (inverse, limiti, conservazione).
Livello 2: valori di riferimento classici.
Livello 3: risultati del corso (catalogo e fonti/soluzioni_A), in fondo al file.
"""
import math
import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import gasdyn as gd  # noqa: E402

ARIA = gd.Gas(1.4, 287.0)
ARIA_ADIM = gd.Gas(1.4)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# =============================================================================
# Area 1 — termodinamica
# =============================================================================
class TestArea1:
    def test_gas_derivate(self):
        g = gd.Gas(1.3, 300.0)
        assert g.delta == pytest.approx(0.15)
        assert g.q == pytest.approx(0.3 / 2.6)
        assert g.cp - g.cv == pytest.approx(300.0)
        assert g.cp / g.cv == pytest.approx(1.3)

    def test_gas_senza_R(self):
        st = gd.State.da_pa(ARIA_ADIM, 1.0, 1.0)
        assert st.rho == pytest.approx(1.4)
        with pytest.raises(gd.GasdynError):
            _ = st.T

    def test_gamma_invalido(self):
        with pytest.raises(gd.GasdynError):
            gd.Gas(1.0)

    def test_stato_proprieta(self):
        st = gd.State.da_pT(ARIA, 101325.0, 288.15, u=100.0)
        assert st.a == pytest.approx(math.sqrt(1.4 * 287 * 288.15))
        assert st.T == pytest.approx(288.15)
        assert st.T0 == pytest.approx(288.15 + 100**2 / (2 * ARIA.cp))
        assert st.a0**2 == pytest.approx(st.a**2 + 0.2 * 100**2)
        assert st.h0 == pytest.approx(ARIA.cp * st.T0)
        assert st.R2 - st.R1 == pytest.approx(2 * 0.2 * 100)

    def test_da_M_e_da_totali(self):
        st = gd.State.da_M(ARIA_ADIM, 2.0, 30.0, p=1.0, a=1.0)
        assert st.M == pytest.approx(2.0)
        assert st.theta == pytest.approx(30.0)
        st2 = gd.State.da_totali(ARIA_ADIM, 2.0, p0=st.p0, a0=st.a0, theta=30.0)
        assert st2.p == pytest.approx(1.0) and st2.a == pytest.approx(1.0)

    def test_isentropiche_inverse(self):
        for M in [0.1, 0.5, 1.0, 2.0, 5.0]:
            r = gd.isentropiche(M, ARIA)
            for tipo, k in [("p0/p", "p0_p"), ("T0/T", "T0_T"), ("rho0/rho", "rho0_rho"),
                            ("a0/a", "a0_a")]:
                assert gd.mach_da_rapporto(r[k], ARIA, tipo) == pytest.approx(M, rel=1e-12)

    def test_isentropiche_valori(self):
        # Valori tabulati γ = 1.4, M = 2: p/p0 = 0.12780, T/T0 = 0.55556
        r = gd.isentropiche(2.0, ARIA)
        assert 1 / r["p0_p"] == pytest.approx(0.127805, rel=1e-5)
        assert 1 / r["T0_T"] == pytest.approx(0.555556, rel=1e-5)
        c = gd.critiche(ARIA)
        assert c["p*_p0"] == pytest.approx(0.528282, rel=1e-5)
        assert c["T*_T0"] == pytest.approx(0.833333, rel=1e-5)

    def test_pressione_totale_modelli(self):
        # B43-E01: M = 0.3 → 1.063 (incompr.), 1.064 (isentr.), 1.065 (isot.)
        M = 0.3
        assert gd.pressione_totale(M, ARIA, "incomprimibile") == pytest.approx(1.063)
        assert gd.pressione_totale(M, ARIA, "isentropico") == pytest.approx(1.0644, abs=5e-5)
        assert gd.pressione_totale(M, ARIA, "isotermo") == pytest.approx(1.0650, abs=5e-5)
        # limite M → 0: i tre modelli coincidono al primo ordine
        for mod in ["incomprimibile", "isentropico", "isotermo"]:
            assert gd.pressione_totale(1e-4, ARIA, mod) - 1 == pytest.approx(0.7e-8, rel=1e-3)

    def test_entropia(self):
        st1 = gd.State(ARIA, 1.2, 0.0, 1e5)
        st2 = gd.isentropica(st1, p=3e5)
        assert gd.salto_entropia(st1, st2)["ds_R"] == pytest.approx(0.0, abs=1e-12)
        st3 = st1.con(p=2e5)  # stessa densità: Δs = cv ln 2
        d = gd.salto_entropia(st1, st3)
        assert d["ds"] == pytest.approx(ARIA.cv * math.log(2))
        # forma cp ln T − R ln p
        alt = ARIA.cp * math.log(st3.T / st1.T) - ARIA.R * math.log(st3.p / st1.p)
        assert d["ds"] == pytest.approx(alt)

    def test_isentropica_T_e_a(self):
        st1 = gd.State.da_pT(ARIA, 1e5, 300.0)
        st2 = gd.isentropica(st1, T=200.0)
        assert st2.p / st1.p == pytest.approx((200 / 300) ** 3.5)
        st3 = gd.isentropica(st1, a=st2.a)
        assert st3.p == pytest.approx(st2.p)

    def test_velocita_max(self):
        st = gd.State.da_M(ARIA_ADIM, 3.0, p=1.0, a=1.0)
        assert gd.velocita_max(st.a0, ARIA_ADIM) == pytest.approx(st.V_max)
        assert st.V_max == pytest.approx(st.a0 * math.sqrt(5))

    def test_U_F(self):
        st = gd.State(ARIA_ADIM, 1.0, 0.5, 1.0)
        st2 = gd.stato_da_U(ARIA_ADIM, st.U())
        assert st2.p == pytest.approx(1.0) and st2.u == pytest.approx(0.5)


# =============================================================================
# Area 2 — isentropico e quasi-1D
# =============================================================================
class TestArea2:
    def test_area_mach_valori(self):
        assert gd.area_su_area_critica(2.0, ARIA) == pytest.approx(1.687500, rel=1e-6)
        assert gd.area_su_area_critica(3.0, ARIA) == pytest.approx(4.234568, rel=1e-6)
        assert gd.area_su_area_critica(0.5, ARIA) == pytest.approx(1.339844, rel=1e-6)
        assert gd.area_su_area_critica(1.0, ARIA) == pytest.approx(1.0)

    def test_area_mach_inversa(self):
        for g in [1.1, 1.3, 1.4, 5 / 3]:
            gas = gd.Gas(g)
            for M in [0.05, 0.3, 0.9, 1.2, 2.5, 6.0]:
                r = gd.area_su_area_critica(M, gas)
                ramo = "subsonico" if M < 1 else "supersonico"
                assert gd.mach_da_area(r, gas, ramo) == pytest.approx(M, rel=1e-10)

    def test_strozzamento(self):
        with pytest.raises(gd.Strozzamento):
            gd.mach_da_area(0.9, ARIA, "subsonico")

    def test_portata(self):
        # ṁ_max: costante 0.040418 √(K)/ (m/s)… per aria: ṁ = 0.0404 A p0/√T0
        m = gd.portata_max(1.0, 1e5, 300.0, ARIA)
        assert m / (1e5 / math.sqrt(300)) == pytest.approx(0.040418, rel=1e-4)
        # coerenza con ρ* a* A*
        st = gd.State.da_totali(ARIA, 1.0, p0=1e5, T0=300.0)
        assert m == pytest.approx(st.rho * st.a)
        # portata conservata lungo l'ugello
        st2 = gd.State.da_totali(ARIA, 2.0, p0=1e5, T0=300.0)
        A2 = gd.area_su_area_critica(2.0, ARIA)
        assert st2.rho * st2.V * A2 == pytest.approx(m)

    def test_convergente(self):
        c = gd.critiche(ARIA)["p*_p0"]
        r = gd.ugello_convergente(1.0, 0.3, ARIA)
        assert r["regime"] == "strozzato" and r["pe"] == pytest.approx(c)
        r = gd.ugello_convergente(1.0, 0.9, ARIA)
        assert r["regime"] == "subsonico"
        assert gd.isentropiche(r["Me"], ARIA)["p0_p"] == pytest.approx(1 / 0.9)

    def test_ugello_cd_regimi(self):
        Ae = 2.0
        r0 = gd.ugello_cd(1.0, 0.5, Ae, ARIA)
        pe1, pe2, pe3 = r0["pe1"], r0["pe2"], r0["pe3"]
        assert pe3 == pytest.approx(0.09393, abs=1e-4)
        assert pe1 == pytest.approx(0.9372, abs=2e-4)
        assert pe3 < pe2 < pe1
        for pb in np.linspace(pe2 * 1.001, pe1 * 0.999, 7):
            r = gd.ugello_cd(1.0, pb, Ae, ARIA)
            assert r["regime"] == "urto nel divergente"
            assert 1.0 < r["As_At"] < Ae
            # ricalcolo indipendente della pressione d'uscita
            assert gd.urto_in_ugello(r["As_At"], Ae, ARIA)["pe_p01"] == pytest.approx(pb, rel=1e-9)
        assert gd.ugello_cd(1.0, 0.95, Ae, ARIA)["regime"] == "subsonico"
        assert gd.ugello_cd(1.0, 0.2, Ae, ARIA)["regime"] == "sovraespanso"
        assert gd.ugello_cd(1.0, pe3, Ae, ARIA)["regime"] == "adattato"
        assert gd.ugello_cd(1.0, 0.05, Ae, ARIA)["regime"] == "sottoespanso"

    def test_variazione_area_conservazione(self):
        st = gd.State.da_M(ARIA_ADIM, 1.87, p=380.0, rho=8.0)
        for A1, A2 in [(1.5, 1.0), (1.0, 2.0)]:
            st5 = gd.variazione_area(st, A1, A2)
            assert st5.rho * st5.u * A2 == pytest.approx(st.rho * st.u * A1, rel=1e-12)
            assert st5.h0 == pytest.approx(st.h0, rel=1e-12)
            assert st5.s == pytest.approx(st.s, rel=1e-12)
            assert (st5.M > 1) == (st.M > 1)
        sub = gd.State.da_M(ARIA_ADIM, 0.5, p=1.0, a=1.0)
        with pytest.raises(gd.Strozzamento):
            gd.variazione_area(sub, 1.0, 0.5)

    def test_svuotamento_vs_ode(self):
        from scipy.integrate import solve_ivp
        gas = ARIA
        V, A, p_i, T_i = 1.0, 1e-4, 1e6, 300.0
        a_i = math.sqrt(1.4 * 287 * T_i)
        Cs = (2 / 2.4) ** (2.4 / 0.8)
        rho_i = p_i / (287 * T_i)

        def ode(t, y):  # y = m
            m = y[0]
            rho = m / V
            a = a_i * (rho / rho_i) ** 0.2
            return [-A * rho * a * Cs]

        sol = solve_ivp(ode, [0, 5], [rho_i * V], rtol=1e-10, atol=1e-12)
        r = gd.serbatoio_svuotamento(5.0, p_i, a_i, V, A, gas)
        assert r["m_mi"] == pytest.approx(sol.y[0, -1] / (rho_i * V), rel=1e-7)
        r = gd.serbatoio_svuotamento(1.0, p_i, a_i, V, A, gas, p_a=1e5)
        pf = 1e5 / gd.critiche(gas)["p*_p0"]
        rt = gd.serbatoio_svuotamento(r["t_fine_strozzamento"], p_i, a_i, V, A, gas)
        assert rt["p_pi"] * p_i == pytest.approx(pf)

    def test_riempimento(self):
        r = gd.serbatoio_riempimento(1.0, 1e5, 300.0, 300.0, ARIA, p=2e5)
        # bilancio: m cv T − m_i cv T_i = (m − m_i) cp T0
        lhs = r["m"] * ARIA.cv * r["T"] - r["m_i"] * ARIA.cv * 300.0
        assert lhs == pytest.approx((r["m"] - r["m_i"]) * ARIA.cp * 300.0)
        assert r["p"] == pytest.approx(2e5)

    def test_spinta(self):
        assert gd.spinta(10.0, 2000.0, 2e5, 1e5, 0.1) == pytest.approx(30000.0)


# =============================================================================
# Area 3 — urto normale
# =============================================================================
class TestArea3:
    def test_valori_classici(self):
        r = gd.urto_normale(2.0, ARIA)
        assert r["p2_p1"] == pytest.approx(4.5)
        assert r["M2"] == pytest.approx(0.57735, abs=1e-5)
        assert r["p02_p01"] == pytest.approx(0.72087, abs=1e-5)
        assert r["rho2_rho1"] == pytest.approx(2.66667, abs=1e-5)
        assert r["T2_T1"] == pytest.approx(1.6875)

    def test_limiti(self):
        r = gd.urto_normale(1.0, ARIA)
        assert r["p2_p1"] == 1 and r["M2"] == pytest.approx(1) and r["p02_p01"] == pytest.approx(1)
        r = gd.urto_normale(1e4, ARIA)
        assert r["rho2_rho1"] == pytest.approx(6.0, rel=1e-6)
        assert r["M2"] == pytest.approx(math.sqrt(0.4 / 2.8), rel=1e-6)
        # urto debole: Δs = O((M−1)³)
        r = gd.urto_normale(1.01, ARIA)
        assert 0 < r["ds_R"] < 1e-5
        with pytest.raises(gd.NessunaSoluzione):
            gd.urto_normale(0.8, ARIA)

    def test_entropia_coerente(self):
        r = gd.urto_normale(3.0, ARIA)
        # Δs/R = −ln(p02/p01)
        assert r["ds_R"] == pytest.approx(-math.log(r["p02_p01"]))

    def test_inverse(self):
        for M in [1.1, 2.0, 4.5]:
            r = gd.urto_normale(M, ARIA)
            for tipo, k in [("p2/p1", "p2_p1"), ("rho2/rho1", "rho2_rho1"), ("T2/T1", "T2_T1"),
                            ("M2", "M2"), ("p02/p01", "p02_p01")]:
                assert gd.mach_da_urto(r[k], ARIA, tipo) == pytest.approx(M, rel=1e-9)

    def test_conservazione_stato(self):
        st1 = gd.State.da_M(ARIA, 2.5, p=1e5, T=250.0)
        st2 = gd.urto_normale_stato(st1)
        assert st2.rho * st2.u == pytest.approx(st1.rho * st1.u)
        assert st2.p + st2.rho * st2.u**2 == pytest.approx(st1.p + st1.rho * st1.u**2)
        assert st2.h0 == pytest.approx(st1.h0)
        # Prandtl v1 v2 = a*²
        assert st1.u * st2.u == pytest.approx(st1.a_star**2)


# =============================================================================
# Area 4 — urti mobili
# =============================================================================
class TestArea4:
    quiete = gd.State.da_pa(ARIA_ADIM, 1.0, 1.0)

    def test_pistone_quadratica(self):
        # F §7: w² − (1+δ) vp w − a0² = 0
        for vp in [0.1, 1.0, 3.0]:
            u = gd.urto_mobile(self.quiete, "u+a", u_valle=vp)
            assert u.w**2 - 1.2 * vp * u.w - 1.0 == pytest.approx(0.0, abs=1e-12)
            assert u.valle.u == pytest.approx(vp)

    def test_rh_residui(self):
        st = gd.State(ARIA_ADIM, 0.8, 0.3, 1.7)
        for fam in ["u+a", "u-a"]:
            for Ms in [1.0, 1.5, 4.0]:
                u = gd.urto_mobile(st, fam, Ms=Ms)
                assert max(u.residui_rh().values()) < 1e-12
                assert u.Ms == pytest.approx(Ms)
                assert u.valle.s >= st.s

    def test_segni(self):
        u = gd.urto_mobile(self.quiete, "u+a", Ms=2.0)
        assert u.Mr_monte == pytest.approx(-2.0) and u.Mr_valle < 0 and u.w > 0
        u = gd.urto_mobile(self.quiete, "u-a", Ms=2.0)
        assert u.Mr_monte == pytest.approx(2.0) and u.Mr_valle > 0 and u.w < 0

    def test_inverse_dati(self):
        st = gd.State(ARIA_ADIM, 1.3, -0.4, 2.0)
        for fam in ["u+a", "u-a"]:
            ref = gd.urto_mobile(st, fam, Ms=2.3)
            for kw in [dict(w=ref.w), dict(p_valle=ref.valle.p), dict(u_valle=ref.valle.u),
                       dict(M_valle=ref.valle.u / ref.valle.a)]:
                u = gd.urto_mobile(st, fam, **kw)
                assert u.Ms == pytest.approx(2.3, rel=1e-9), kw
            assert gd.urto_mobile(st, w=ref.w).famiglia == fam

    def test_galileo(self):
        U = 0.7
        a = gd.urto_mobile(self.quiete, "u+a", Ms=3.0)
        b = gd.urto_mobile(self.quiete.con(u=U), "u+a", Ms=3.0)
        assert b.w == pytest.approx(a.w + U) and b.valle.u == pytest.approx(a.valle.u + U)

    def test_simmetria_famiglie(self):
        a = gd.urto_mobile(self.quiete, "u+a", Ms=2.0)
        b = gd.urto_mobile(self.quiete, "u-a", Ms=2.0)
        assert b.w == pytest.approx(-a.w) and b.valle.u == pytest.approx(-a.valle.u)
        assert b.valle.p == pytest.approx(a.valle.p)

    def test_M_valle_limite(self):
        assert gd.M_valle_limite(ARIA_ADIM) == pytest.approx(1.88982, abs=1e-5)
        with pytest.raises(gd.NessunaSoluzione):
            gd.urto_mobile(self.quiete, "u+a", M_valle=1.89)
        u = gd.urto_mobile(self.quiete, "u+a", M_valle=1.87)
        assert u.valle.u / u.valle.a == pytest.approx(1.87)

    def test_errori(self):
        with pytest.raises(gd.NessunaSoluzione):
            gd.urto_mobile(self.quiete, "u+a", p_valle=0.5)
        with pytest.raises(gd.NessunaSoluzione):
            gd.urto_mobile(self.quiete, "u+a", u_valle=-0.5)
        with pytest.raises(gd.NessunaSoluzione):
            gd.urto_mobile(self.quiete, "u-a", w=0.5)
        with pytest.raises(gd.NessunaSoluzione):
            gd.urto_mobile(self.quiete, w=0.5)  # |M_r1| < 1

    def test_riflesso_parete(self):
        inc = gd.urto_mobile(self.quiete, "u+a", Ms=2.0)
        rif = gd.urto_riflesso_parete(inc.valle, "u-a")
        assert rif.valle.u == pytest.approx(0.0, abs=1e-12)
        assert rif.w < 0 and max(rif.residui_rh().values()) < 1e-12
        # formula classica del riflesso (fuori formulario): p5/p2 = [(3γ−1)Ms² − 2(γ−1)]/[(γ−1)Ms² + 2]
        g, Ms = 1.4, 2.0
        p5_p2 = ((3 * g - 1) * Ms**2 - 2 * (g - 1)) / ((g - 1) * Ms**2 + 2)
        assert rif.valle.p / inc.valle.p == pytest.approx(p5_p2, rel=1e-12)

    def test_natura_salto(self):
        a = gd.State(ARIA_ADIM, 3.8125, 3.9161, 10.0, nome="t<t*")
        b = gd.State(ARIA_ADIM, 1.0, 6.4928, 1.0, nome="t>t*")
        r = gd.natura_salto(a, b)
        assert r["tipo"] == "urto u-a" and r["w"] == pytest.approx(3.0, abs=5e-3)
        st = gd.State(ARIA_ADIM, 1.0, 0.0, 1.0)
        assert gd.natura_salto(st, st.con(rho=2.0))["tipo"] == "contatto"


# =============================================================================
# Area 5 — onde semplici, pistoni, formazione dell'urto
# =============================================================================
class TestArea5:
    st = gd.State.da_pa(ARIA_ADIM, 1.0, 1.0, u=0.3)

    def test_punto_caratteristiche(self):
        P = gd.punto_caratteristiche(self.st, self.st)
        assert P.u == pytest.approx(self.st.u) and P.a == pytest.approx(self.st.a)
        A = gd.onda_semplice(self.st, "u+a", u=0.5)  # stesso R1
        B = gd.onda_semplice(self.st, "u-a", u=0.1)  # stesso R2
        P = gd.punto_caratteristiche(A, B)
        assert P.R2 == pytest.approx(A.R2) and P.R1 == pytest.approx(B.R1)

    def test_onda_semplice_invarianti(self):
        for fam, inv in [("u+a", "R1"), ("u-a", "R2")]:
            for kw in [dict(u=0.6), dict(p=0.5), dict(a=1.1), dict(rho=1.2)]:
                s2 = gd.onda_semplice(self.st, fam, **kw)
                assert getattr(s2, inv) == pytest.approx(getattr(self.st, inv))
                assert s2.s == pytest.approx(self.st.s)
        with pytest.raises(gd.Vuoto):
            gd.onda_semplice(self.st, "u+a", u=self.st.u - 1.0 / 0.2 - 0.01)

    def test_curva_onda(self):
        for fam in ["u+a", "u-a"]:
            # ramo d'espansione = onda semplice; ramo d'urto = Hugoniot
            assert gd.curva_onda(self.st, fam, 0.4) == pytest.approx(
                gd.onda_semplice(self.st, fam, p=0.4).u)
            assert gd.curva_onda(self.st, fam, 3.0) == pytest.approx(
                gd.urto_mobile(self.st, fam, p_valle=3.0).valle.u)
            # contatto del secondo ordine in p_i (urto debole ≈ isentropico)
            h = 1e-3
            up = gd.curva_onda(self.st, fam, 1 + h)
            um = gd.curva_onda(self.st, fam, 1 - h)
            u0 = gd.curva_onda(self.st, fam, 1.0)
            d2 = (up - 2 * u0 + um) / h**2
            iso = (gd.onda_semplice(self.st, fam, p=1 + h).u - 2 * u0
                   + gd.onda_semplice(self.st, fam, p=1 - h).u) / h**2
            assert d2 == pytest.approx(iso, rel=1e-3)

    def test_ventaglio(self):
        for fam in ["u+a", "u-a"]:
            V = gd.onda_verso(self.st, fam, p=0.3)
            assert isinstance(V, gd.Ventaglio) and V.sigla == "R"
            e = 1e-9
            s_testa = V.stato_xi(V.xi_testa - (1 if fam == "u+a" else -1) * e)
            s_coda = V.stato_xi(V.xi_coda + (1 if fam == "u+a" else -1) * e)
            assert s_testa.p == pytest.approx(self.st.p, rel=1e-7)
            assert s_coda.p == pytest.approx(0.3, rel=1e-7)
            xm = 0.5 * (V.xi_testa + V.xi_coda)
            sm = V.stato_xi(xm)
            lam = sm.u + (1 if fam == "u+a" else -1) * sm.a
            assert lam == pytest.approx(xm)  # ξ = u ± a dentro il ventaglio
            assert sm.s == pytest.approx(self.st.s)

    def test_ventaglio_transonico(self):
        quiete = gd.State.da_pa(ARIA_ADIM, 1.0, 1.0)
        # pistone sinistro che arretra con v_p > a0/(1+δ): coda supersonica
        V = gd.pistone_impulsivo(quiete, -1.0, "sinistro")
        assert V.xi_coda == pytest.approx(1.0 - 1.2 * 1.0)
        assert V.transonico
        s0 = V.stato_xi(0.0)
        assert abs(s0.u) == pytest.approx(s0.a)

    def test_pistone_vuoto(self):
        quiete = gd.State.da_pa(ARIA_ADIM, 1.0, 1.0)
        V = gd.pistone_impulsivo(quiete, -6.0, "sinistro")
        assert V.vuoto and V.xi_coda == pytest.approx(-5.0)
        with pytest.raises(gd.Vuoto):
            gd.onda_verso(quiete, "u+a", u=-6.0)

    def test_pistone_compressione(self):
        quiete = gd.State.da_pa(ARIA_ADIM, 1.0, 1.0)
        S = gd.pistone_impulsivo(quiete, 1.0, "sinistro")
        assert isinstance(S, gd.Urto) and S.w == pytest.approx((1.2 + math.sqrt(1.44 + 4)) / 2)
        S = gd.pistone_impulsivo(quiete, -1.0, "destro")
        assert S.w == pytest.approx(-(1.2 + math.sqrt(1.44 + 4)) / 2)

    def test_pistone_legge_parabolica(self):
        # F §6: t_BU = a0/((1+δ)ẍ(0)), x_BU = a0 t_BU (pistone da fermo)
        quiete = gd.State(ARIA_ADIM, 1.0, 0.0, 1.0)  # a0 = √γ
        P = gd.pistone_legge(lambda t: 0.74 * t**2, quiete, 1.0)
        a0 = math.sqrt(1.4)
        assert P.t_bu == pytest.approx(a0 / (1.2 * 1.48), rel=1e-9)
        assert P.x_bu == pytest.approx(a0 * P.t_bu, rel=1e-9)
        # stato sulla faccia del pistone e lungo una caratteristica
        t = 0.5
        sp = P.stato(0.74 * t**2 + 1e-12, t)
        assert sp.u == pytest.approx(1.48 * t, rel=1e-6)
        x0, t0, lam = P.caratteristica(0.2)
        s = P.stato(x0 + lam * (0.6 - t0), 0.6)
        assert s.u == pytest.approx(1.48 * 0.2, rel=1e-6)

    def test_pistone_legge_destro(self):
        quiete = gd.State.da_pa(ARIA_ADIM, 1.0, 1.0)
        P = gd.pistone_legge(lambda t: 2.0 - t**2, quiete, 1.0, lato="destro")
        assert P.t_bu == pytest.approx(1.0 / (1.2 * 2.0), rel=1e-9)
        assert P.x_bu == pytest.approx(2.0 - P.t_bu, rel=1e-9)

    def test_compressione_centrata(self):
        # B71-E02: legge che genera una compressione centrata in (x, t) = (1, 1)
        g = 1.4
        quiete = gd.State.da_pa(ARIA_ADIM, 1.0, 1.0)
        xp = lambda t: 1 + (2 * (1 - t) - (g + 1) * (1 - t) ** (2 / (g + 1))) / (g - 1)
        P = gd.pistone_legge(xp, quiete, 0.5)
        assert P.t_bu == pytest.approx(1.0, rel=1e-6)
        assert P.x_bu == pytest.approx(1.0, rel=1e-6)
        assert P.t_inviluppo(0.1) == pytest.approx(1.0, rel=1e-5)
        assert P.t_inviluppo(0.4) == pytest.approx(1.0, rel=1e-5)

    def test_formazione_ci_onda_semplice(self):
        # B64-E03: R2 = a + δu uniforme ⇒ onda semplice u−a
        u0 = lambda x: 3 - 2.5 * x**2 + 1.667 * x**3
        a0 = lambda x: 1 + 0.5 * x**2 - 0.3334 * x**3
        r = gd.formazione_urto_ci(u0, a0, ARIA_ADIM, (0.0, 1.0))
        assert r["prima"] == "u-a" and r["u-a"]["onda_semplice"]
        v = lambda x: u0(x) - a0(x)
        dvmin = -6 * 0.5 + 6.0012 * 0.25  # circa, controllo grossolano
        assert r["u-a"]["t_bu"] == pytest.approx(-1 / dvmin, rel=2e-3)
        xi = r["u-a"]["xi"]
        assert r["u-a"]["x_bu"] == pytest.approx(xi + v(xi) * r["u-a"]["t_bu"])
        # con array
        xs = np.linspace(0, 1, 20001)
        r2 = gd.formazione_urto_ci(u0(xs), a0(xs), ARIA_ADIM, x=xs)
        assert r2["u-a"]["t_bu"] == pytest.approx(r["u-a"]["t_bu"], rel=1e-6)

    def test_t_bu_lineare(self):
        assert gd.t_bu_lineare(0.0, 2.0, 1.0, 1.0) == pytest.approx(1.0)
        assert gd.t_bu_lineare(0.0, 1.0, 1.0, 2.0) == math.inf

    def test_urto_isotermo(self):
        r = gd.urto_isotermo(3.0)
        assert r["rho2_rho1"] == 9 and r["v2_v1"] * 9 == pytest.approx(1)
        assert r["M2"] * 3 == pytest.approx(1)  # Prandtl v1 v2 = a²


# =============================================================================
# Area 6 — problema di Riemann
# =============================================================================
def _st(rho, u, p, gas=ARIA_ADIM):
    return gd.State(gas, rho, u, p)


class TestArea6:
    def test_sod(self):
        s = gd.riemann(_st(1.0, 0.0, 1.0), _st(0.125, 0.0, 0.1))
        assert s.p_star == pytest.approx(0.30313, abs=1e-5)
        assert s.u_star == pytest.approx(0.92745, abs=1e-5)
        assert s.stati["2"].rho == pytest.approx(0.42632, abs=1e-5)
        assert s.stati["3"].rho == pytest.approx(0.26557, abs=1e-5)
        assert s.onda_dx.w == pytest.approx(1.75216, abs=1e-5)
        assert s.configurazione == "RCS" and s.metodo == "newton"
        assert len(s.iterazioni) < 10

    def test_toro(self):
        # E. F. Toro, Riemann Solvers…, tab. 4.2 (test 2, 3, 5)
        s = gd.riemann(_st(1.0, -2.0, 0.4), _st(1.0, 2.0, 0.4))
        assert s.p_star == pytest.approx(0.00189, abs=1e-5) and s.u_star == pytest.approx(0.0, abs=1e-10)
        assert s.configurazione == "RNR"
        s = gd.riemann(_st(1.0, 0.0, 1000.0), _st(1.0, 0.0, 0.01))
        assert s.p_star == pytest.approx(460.894, abs=1e-3) and s.u_star == pytest.approx(19.5975, abs=1e-4)
        s = gd.riemann(_st(5.99924, 19.5975, 460.894), _st(5.99242, -6.19633, 46.0950))
        # dati d'ingresso arrotondati a 6 cifre: confronto a 1e-5 relativo
        assert s.p_star == pytest.approx(1691.64, rel=1e-5) and s.u_star == pytest.approx(8.68975, rel=1e-5)
        assert s.configurazione == "SCS"

    def test_vuoto(self):
        s = gd.riemann(_st(1.0, -7.0, 1.0), _st(1.0, 7.0, 1.0))
        assert s.vuoto and s.configurazione == "RVR"
        a = math.sqrt(1.4)
        assert s.onda_sx.xi_coda == pytest.approx(-7 + a / 0.2)
        assert s.stato_xi(0.0).p == 0.0

    def test_coerenza_curve(self):
        # p* e u* stanno su entrambe le curve del piano p–u (F §8)
        casi = [(_st(1, 0, 1), _st(0.125, 0, 0.1)), (_st(1, 1, 1), _st(1, -1, 1)),
                (_st(2, 0.3, 5), _st(0.5, -0.2, 0.7)), (_st(1, -1, 1), _st(1, 0.5, 2))]
        for L, R in casi:
            s = gd.riemann(L, R)
            assert gd.curva_onda(L, "u-a", s.p_star) == pytest.approx(s.u_star, abs=1e-10)
            assert gd.curva_onda(R, "u+a", s.p_star) == pytest.approx(s.u_star, abs=1e-10)
            for o in (s.onda_sx, s.onda_dx):
                if isinstance(o, gd.Urto):
                    assert max(o.residui_rh().values()) < 1e-11

    def test_specchio(self):
        L, R = _st(2, 0.3, 5), _st(0.5, -0.2, 0.7)
        a = gd.riemann(L, R)
        b = gd.riemann(R.con(u=-R.u), L.con(u=-L.u))
        assert b.p_star == pytest.approx(a.p_star) and b.u_star == pytest.approx(-a.u_star)
        assert b.configurazione == a.configurazione[::-1]

    def test_gas_diversi(self):
        aria, elio = gd.Gas(1.4), gd.Gas(5 / 3)
        L, R = gd.State(aria, 1.0, 0.0, 10.0), gd.State(elio, 0.2, 0.0, 1.0)
        s = gd.riemann(L, R)
        assert s.stati["2"].gas == aria and s.stati["3"].gas == elio
        assert gd.curva_onda(L, "u-a", s.p_star) == pytest.approx(s.u_star, abs=1e-10)
        assert gd.curva_onda(R, "u+a", s.p_star) == pytest.approx(s.u_star, abs=1e-10)

    def test_isentropico(self):
        L, R = _st(1, 0, 1), _st(1, 0, 0.5)
        s = gd.riemann(L, R, modello="isentropico")
        # forma chiusa F §8: u⁰
        g, d = 1.4, 0.2
        aL, aR = math.sqrt(1.4), math.sqrt(0.7)
        z = aR / aL * (1 / 0.5) ** ((g - 1) / (2 * g))
        u0 = (z * aL - aR) / (d * (1 + z))
        assert s.u_star == pytest.approx(u0)
        e = gd.riemann(L, R)
        assert s.p_star == pytest.approx(e.p_star, rel=2e-3)  # urto debole ≈ isentropico

    def test_conservazione_integrale(self):
        L, R = _st(1, 0.2, 1), _st(0.3, -0.5, 0.2)
        s = gd.riemann(L, R)
        t, X = 0.4, 3.0
        x = np.linspace(-X, X, 400001)
        c = s.campiona(x, t)
        U = np.array([c["rho"], c["rho"] * c["u"], c["p"] / 0.4 + 0.5 * c["rho"] * c["u"] ** 2])
        I = np.trapezoid(U, x, axis=1)
        I0 = X * (L.U() + R.U())
        atteso = I0 + t * (L.F() - R.F())
        assert np.allclose(I, atteso, rtol=2e-5)

    def test_parete_e_pressione(self):
        st = _st(1, 0.5, 1)
        o = gd.contorno_parete(st, "destra")
        assert isinstance(o, gd.Urto) and o.valle.u == pytest.approx(0, abs=1e-12) and o.w < 0
        o = gd.contorno_parete(st, "sinistra")
        assert isinstance(o, gd.Ventaglio) and o.valle.u == pytest.approx(0, abs=1e-12)
        o = gd.contorno_pressione(st, "destra", 0.5)
        assert o.valle.p == pytest.approx(0.5)

    def test_serbatoio(self):
        gas = ARIA_ADIM
        interno = _st(1.0, 0.0, 0.5)
        r = gd.contorno_serbatoio(interno, "sinistra", 1.0, a0=math.sqrt(1.4))
        b = r["stato"]
        assert r["regime"] == "ingresso" and b.u > 0
        assert b.a**2 + 0.2 * b.u**2 == pytest.approx(1.4)
        assert b.p0 == pytest.approx(1.0) and b.s == pytest.approx(1.0)  # s e h0 del serbatoio
        si = r["stato_interno"]
        assert si.p == pytest.approx(b.p) and si.u == pytest.approx(b.u) and r["contatto"]
        assert isinstance(r["onda"], gd.Urto)  # l'ingresso comprime il gas interno
        # S p.45: stessa entropia e modello isentropico ⇒ invariante R1 dall'interno, nessun contatto
        interno = gd.isentropica(gd.State(ARIA_ADIM, 1.0, 0.0, 1.0), p=0.5)
        r = gd.contorno_serbatoio(interno, "sinistra", 1.0, a0=math.sqrt(1.4), modello="isentropico")
        assert r["stato"].R1 == pytest.approx(interno.R1) and not r["contatto"]
        r = gd.contorno_serbatoio(_st(1.0, 0.0, 2.0), "sinistra", 1.0, a0=1.0)
        assert r["regime"] == "uscita" and r["stato"].p == pytest.approx(1.0)

    def test_riemann_area_aree_uguali(self):
        L, R = _st(1.0, 0.3, 1.0), _st(1.0, 0.3, 0.4)
        r = gd.riemann_area(L, R, 1.0, 1.0)
        s = gd.riemann(L, R)
        assert r["stato_5"].p == pytest.approx(s.p_star, rel=1e-9)
        assert r["stato_5"].u == pytest.approx(s.u_star, rel=1e-9)

    # --- riemann_area: regime strozzato (fuori formulario) ---
    @staticmethod
    def _incidente(Mv):
        q = gd.State.da_pa(ARIA_ADIM, 1, 1)
        return gd.urto_mobile(q, "u+a", M_valle=Mv).valle, q

    @staticmethod
    def _p_trasmesso(r):
        return r["onda_trasmessa"].valle.p

    @staticmethod
    def _verifica_collegamento(r, A_sx, A_dx, q):
        """2' → 5 conserva ṁ e h0; 5 (o lo stato dietro l'onda u−a a valle) sta sulla curva u+a di q."""
        s2, s5 = r["stato_2p"], r["stato_5"]
        assert s2.rho * s2.u * A_sx == pytest.approx(s5.rho * s5.u * A_dx, rel=1e-9)
        assert s2.h0 == pytest.approx(s5.h0, rel=1e-10)
        t = r["onda_trasmessa"]
        assert isinstance(t, gd.Urto) and t.monte == q
        sol = r.get("riemann_valle")
        st = sol.stati["2"] if sol is not None else s5
        assert t.valle.p == pytest.approx(st.p, rel=1e-8) and t.valle.u == pytest.approx(st.u, rel=1e-8)

    @pytest.mark.parametrize("Mv,A_sx,A_dx", [(0.5, 1.0, 1.5), (0.7, 1.0, 1.2), (0.7, 1.0, 1.5)])
    def test_riemann_area_strozzato_urto_nel_salto(self, Mv, A_sx, A_dx):
        # A crescente: 2' sonico in A_sx, urto normale stazionario dentro il salto (F §16), 5 subsonico
        L, q = self._incidente(Mv)
        r = gd.riemann_area(L, q, A_sx, A_dx)
        assert r["regime"] == "strozzato" and r["sezione_sonica"] == "A_sx"
        assert r["stato_2p"].Mu == pytest.approx(1.0, abs=1e-9)
        assert isinstance(r["onda_riflessa"], gd.Ventaglio)
        us = r["urto_stazionario"]
        assert A_sx < us["A_s"] < A_dx and us["monte"].Mu > 1 > us["valle"].Mu
        assert r["stato_5"].Mu < 1 and r["stato_5"].p0 < r["stato_2p"].p0
        self._verifica_collegamento(r, A_sx, A_dx, q)
        assert self._p_trasmesso(r) < L.p  # A↑: urto trasmesso più debole

    @pytest.mark.parametrize("A_dx", [1.2, 1.5])
    def test_riemann_area_strozzato_supersonico_a_valle(self, A_dx):
        # A crescente, urto più forte: 5 supersonico, urto u−a che si allontana dalla sezione (S p.49)
        L, q = self._incidente(0.9)
        r = gd.riemann_area(L, q, 1.0, A_dx)
        assert r["regime"] == "strozzato" and r["urto_stazionario"] is None
        assert r["stato_2p"].Mu == pytest.approx(1.0, abs=1e-9) and r["stato_5"].Mu > 1
        o = r["riemann_valle"].onda_sx
        assert isinstance(o, gd.Urto) and o.w >= 0
        self._verifica_collegamento(r, 1.0, A_dx, q)
        assert self._p_trasmesso(r) < L.p

    @pytest.mark.parametrize("Mv,A_sx", [(0.9, 3.0), (1.2, 1.2), (1.2, 1.5)])
    def test_riemann_area_strozzato_area_decrescente(self, Mv, A_sx):
        # A decrescente: urto riflesso, A_dx sonica, espansione u−a con testa ferma in x0
        L, q = self._incidente(Mv)
        r = gd.riemann_area(L, q, A_sx, 1.0)
        assert r["regime"] == "strozzato" and r["sezione_sonica"] == "A_dx"
        assert isinstance(r["onda_riflessa"], gd.Urto) and r["onda_riflessa"].w < 0
        assert r["stato_5"].Mu == pytest.approx(1.0, abs=1e-12)
        M2 = gd.mach_da_area_tra_sezioni(1.0, 1.0, A_sx, ARIA_ADIM, "subsonico")
        assert r["stato_2p"].Mu == pytest.approx(M2, rel=1e-9)
        o = r["riemann_valle"].onda_sx
        assert isinstance(o, gd.Ventaglio) and o.xi_testa == pytest.approx(0.0, abs=1e-12)
        self._verifica_collegamento(r, A_sx, 1.0, q)
        assert self._p_trasmesso(r) > L.p  # A↓: urto trasmesso più forte

    @pytest.mark.parametrize("A_sx", [1.2, 1.5])
    def test_riemann_area_subsonico_vicino_al_limite(self, A_sx):
        # Mv = 0.9, A decrescente: radice subsonica a ridosso dello strozzamento (prima persa)
        L, q = self._incidente(0.9)
        r = gd.riemann_area(L, q, A_sx, 1.0)
        assert r["regime"] == "subsonico" and 0.9 < r["stato_5"].Mu < 1
        self._verifica_collegamento(r, A_sx, 1.0, q)
        assert self._p_trasmesso(r) > L.p

    @pytest.mark.parametrize("A_sx,A_dx", [(1.0, 1.5), (1.5, 1.0)])
    def test_riemann_area_continuita_fra_regimi(self, A_sx, A_dx):
        # al variare di Mv l'urto trasmesso cresce con continuità attraverso i cambi di regime
        Mvs = np.linspace(0.3, 1.5, 25)
        ris = [gd.riemann_area(*self._incidente(Mv), A_sx, A_dx) for Mv in Mvs]
        p = np.array([self._p_trasmesso(r) for r in ris])
        pL = np.array([self._incidente(Mv)[0].p for Mv in Mvs])
        assert np.all(np.diff(p) > 0)
        assert np.all(p < pL) if A_dx > A_sx else np.all(p > pL)
        assert len({r["regime"] for r in ris}) >= 2
        # salto fra campioni adiacenti dello stesso ordine di quello dell'urto incidente
        assert np.max(np.diff(p) / np.diff(pL)) < 2.0


# =============================================================================
# Area 7 — supersonico 2D
# =============================================================================
class TestArea7:
    def test_prandtl_meyer(self):
        assert gd.prandtl_meyer(2.0, ARIA) == pytest.approx(26.3798, abs=1e-4)
        assert gd.nu_max(ARIA) == pytest.approx(130.4541, abs=1e-4)
        for M in [1.0, 1.5, 3.0, 10.0]:
            assert gd.mach_da_pm(gd.prandtl_meyer(M, ARIA), ARIA) == pytest.approx(M, rel=1e-10)
        with pytest.raises(gd.Vuoto):
            gd.mach_da_pm(131.0, ARIA)

    def test_theta_sigma_M(self):
        assert gd.sigma_da_theta(2.0, 10.0, ARIA, "debole") == pytest.approx(39.3139, abs=1e-4)
        assert gd.sigma_da_theta(2.0, 10.0, ARIA, "forte") == pytest.approx(83.70, abs=1e-2)
        tm, sm = gd.theta_max(2.0, ARIA)
        assert tm == pytest.approx(22.9735, abs=1e-4)
        # forma chiusa di σ a θmax (fuori formulario) come controllo
        g, M = 1.4, 2.0
        s2 = ((g + 1) * M**2 - 4 + math.sqrt((g + 1) * ((g + 1) * M**4 + 8 * (g - 1) * M**2 + 16))) / (4 * g * M**2)
        assert sm == pytest.approx(math.degrees(math.asin(math.sqrt(s2))), abs=1e-6)
        with pytest.raises(gd.UrtoStaccato):
            gd.sigma_da_theta(2.0, 23.0, ARIA)
        for th in [0.0, 5.0, 20.0]:
            for ramo in ["debole", "forte"]:
                s = gd.sigma_da_theta(3.0, th, ARIA, ramo)
                assert gd.theta_da_sigma(3.0, s, ARIA) == pytest.approx(th, abs=1e-10)

    def test_mach_da_theta_sigma(self):
        # andata e ritorno con θ(σ; M) su entrambi i rami, anche con γ diverso
        for gas in [ARIA, gd.Gas(1.3)]:
            for M in [1.2, 2.0, 4.0, 10.0]:
                mu = gd.angolo_mach(M)
                for sg in np.linspace(mu + 0.5, 89.5, 7):
                    th = gd.theta_da_sigma(M, sg, gas)
                    assert gd.mach_da_theta_sigma(th, sg, gas) == pytest.approx(M, rel=1e-9)
        # onda di Mach: θ = 0 ⇒ M = 1/sinσ; urto normale
        assert gd.mach_da_theta_sigma(0.0, 30.0, ARIA) == pytest.approx(2.0, rel=1e-12)
        assert gd.mach_da_theta_sigma(0.0, 90.0, ARIA) == pytest.approx(1.0)
        # B63-E01: urto di riattacco, θ = 18°, σ = 30,07°
        assert gd.mach_da_theta_sigma(18.0, 30.07, ARIA) == pytest.approx(4.048, abs=1e-3)
        for th, sg in [(20.0, 15.0), (10.0, 90.0), (-1.0, 30.0), (40.0, 50.0)]:
            with pytest.raises(gd.NessunaSoluzione):
                gd.mach_da_theta_sigma(th, sg, ARIA)

    def test_limiti_polare(self):
        mu = gd.angolo_mach(2.0)
        assert gd.theta_da_sigma(2.0, mu, ARIA) == pytest.approx(0.0, abs=1e-12)
        assert gd.theta_da_sigma(2.0, 90.0, ARIA) == pytest.approx(0.0, abs=1e-12)
        ts, ss = gd.sonico(2.0, ARIA)
        tm, sm = gd.theta_max(2.0, ARIA)
        assert ts < tm and ss < sm  # il punto sonico precede θmax sul ramo debole
        st = gd.State.da_M(ARIA_ADIM, 2.0, p=1.0, a=1.0)
        assert gd.urto_obliquo(st, sigma=ss).valle.M == pytest.approx(1.0, abs=1e-9)

    def test_urto_obliquo_conservazione(self):
        st = gd.State.da_M(ARIA_ADIM, 3.0, 7.0, p=2.0, a=1.3)
        for fam in ["C+", "C-"]:
            o = gd.urto_obliquo(st, theta=12.0, famiglia=fam)
            s = 1 if fam == "C+" else -1
            assert o.deflessione == pytest.approx(s * 12.0)
            n = np.array([-math.sin(math.radians(o.angolo_fronte)), math.cos(math.radians(o.angolo_fronte))])
            t = np.array([math.cos(math.radians(o.angolo_fronte)), math.sin(math.radians(o.angolo_fronte))])
            V1, V2 = np.array([st.u, st.v]), np.array([o.valle.u, o.valle.v])
            assert V1 @ t == pytest.approx(V2 @ t)
            assert st.rho * (V1 @ n) == pytest.approx(o.valle.rho * (V2 @ n))
            assert st.p + st.rho * (V1 @ n) ** 2 == pytest.approx(o.valle.p + o.valle.rho * (V2 @ n) ** 2)
            assert o.valle.h0 == pytest.approx(st.h0)
            assert o.ramo == "debole"

    def test_rotazione(self):
        a = gd.urto_obliquo(gd.State.da_M(ARIA_ADIM, 2.5, 0.0, p=1, a=1), theta=8.0)
        b = gd.urto_obliquo(gd.State.da_M(ARIA_ADIM, 2.5, 25.0, p=1, a=1), theta=8.0)
        assert b.valle.p == pytest.approx(a.valle.p) and b.valle.theta == pytest.approx(a.valle.theta + 25)

    def test_pm_invarianti(self):
        st = gd.State.da_M(ARIA_ADIM, 2.0, 5.0, p=1, a=1)
        o = gd.espansione_pm(st, 10.0, "C+")
        assert o.valle.theta == pytest.approx(-5.0)
        assert o.valle.theta + gd.prandtl_meyer(o.valle.M, ARIA_ADIM) == pytest.approx(
            st.theta + gd.prandtl_meyer(2.0, ARIA_ADIM))
        assert o.valle.p0 == pytest.approx(st.p0)
        o = gd.espansione_pm(st, 10.0, "C-")
        assert o.valle.theta == pytest.approx(15.0)

    def test_onda_verso_p(self):
        st = gd.State.da_M(ARIA_ADIM, 2.5, p=1, a=1)
        for p2 in [0.3, 1.0, 2.0, 6.0]:
            for fam in ["C+", "C-"]:
                assert gd.onda_verso_p_2d(st, p2, fam).valle.p == pytest.approx(p2)
        with pytest.raises(gd.NessunaSoluzione):
            gd.onda_verso_p_2d(st, 8.0, "C+")  # oltre l'urto normale (7.125)

    def test_riemann_stazionario(self):
        sopra = gd.State.da_M(ARIA_ADIM, 4.0, p=0.5, a=1.0)
        sotto = gd.State.da_M(ARIA_ADIM, 6.0, p=0.01, a=0.8)
        r = gd.riemann_stazionario(sopra, sotto)
        assert r.onda_sopra.valle.p == pytest.approx(r.onda_sotto.valle.p)
        assert r.onda_sopra.valle.theta == pytest.approx(r.onda_sotto.valle.theta)
        assert r.onda_sopra.tipo == "espansione" and r.onda_sotto.tipo == "urto"
        assert r.theta < 0
        ri = gd.riemann_stazionario(sopra, sotto, modello="isentropico")
        assert ri.onda_sotto.tipo == "compressione isentropica"
        assert ri.p != pytest.approx(r.p, rel=1e-3)

    def test_riflessione_regolare(self):
        st = gd.State.da_M(ARIA_ADIM, 2.8, p=1, a=1)
        r = gd.riflessione_regolare(st, 10.0)
        assert r["3"].theta == pytest.approx(0.0, abs=1e-10)
        assert gd.theta_da_sigma(2.8, r["incidente"].sigma, ARIA) == pytest.approx(10.0)
        assert gd.theta_da_sigma(r["2"].M, r["riflesso"].sigma, ARIA) == pytest.approx(10.0)
        assert r["angolo_riflesso_parete"] == pytest.approx(r["riflesso"].sigma - 10.0)
        with pytest.raises(gd.UrtoStaccato):
            gd.riflessione_regolare(gd.State.da_M(ARIA_ADIM, 2.0, p=1, a=1), 14.0)

    def test_criteri_rr_mr(self):
        M = 4.0
        tD = gd.theta_distacco(M, ARIA)
        st = gd.State.da_M(ARIA_ADIM, M, p=1, a=1)
        M2 = gd.urto_obliquo(st, theta=tD, famiglia="C-").valle.M
        assert gd.theta_max(M2, ARIA)[0] == pytest.approx(tD, abs=1e-8)
        tN = gd.theta_von_neumann(M, ARIA)
        assert tN < tD
        pn = 1.0 * gd.urto_normale(M, ARIA)["p2_p1"]
        assert gd.riflessione_regolare(st, tN)["3"].p == pytest.approx(pn, rel=1e-9)
        mr = gd.riflessione_mach(st, tN)["principale"]
        assert mr["theta"] == pytest.approx(0.0, abs=1e-5) and mr["p"] == pytest.approx(pn, rel=1e-6)
        mr = gd.riflessione_mach(st, 0.5 * (tN + tD))["principale"]
        assert mr["3"].p == pytest.approx(mr["4"].p) and mr["4"].M < 1
        with pytest.raises(gd.NessunaSoluzione):
            gd.theta_von_neumann(2.0, ARIA)  # esiste solo per M1 > 2.2 circa (γ = 1.4)

    def test_shock_expansion_lastra(self):
        M, al = 3.0, 2.0
        r = gd.shock_expansion([(0, 0), (1, 0)], [(0, 0), (1, 0)], M, al, ARIA)
        lin = 4 * math.radians(al) / math.sqrt(M**2 - 1)
        assert r["CL"] == pytest.approx(lin, rel=0.03)
        assert r["CD"] == pytest.approx(r["CL"] * math.tan(math.radians(al)), rel=1e-10)
        # profilo a diamante simmetrico a incidenza nulla: CL = 0, CD > 0
        r = gd.shock_expansion([(0, 0), (0.5, 0.05), (1, 0)], [(0, 0), (0.5, -0.05), (1, 0)], M, 0.0, ARIA)
        assert r["CL"] == pytest.approx(0.0, abs=1e-12) and r["CD"] > 0

    def test_pitot(self):
        assert gd.pitot_supersonico(2.0, ARIA) == pytest.approx(5.6404, abs=1e-4)
        for M in [0.5, 1.0, 1.5, 4.0]:
            assert gd.mach_da_pitot(gd.pitot_supersonico(M, ARIA), ARIA) == pytest.approx(M, rel=1e-9)
        r = gd.urto_normale(3.0, ARIA)
        assert gd.pitot_supersonico(3.0, ARIA) == pytest.approx(r["p02_p1"])

    def test_caratteristiche_2d(self):
        st = gd.State.da_M(ARIA_ADIM, 2.0, p=1, a=1)
        r = gd.punto_caratteristiche_2d(st, st, (0, 1), (0, -1))
        assert r["M"] == pytest.approx(2.0) and r["theta"] == pytest.approx(0.0)
        assert r["x"] == pytest.approx(1 / math.tan(math.radians(30.0))) and r["y"] == pytest.approx(0, abs=1e-12)

    def test_getto(self):
        r = gd.getto_adattamento(2.0, 2.0, 1.0, ARIA)
        st = gd.State.da_M(ARIA, 2.0, p=2.0, a=1.0)
        assert r["regime"] == "sottoespanso"
        assert r["deflessione"] == pytest.approx(abs(gd.onda_verso_p_2d(st, 1.0, "C+").deflessione))
        r = gd.getto_adattamento(2.0, 0.8, 1.0, ARIA)
        assert r["regime"] == "sovraespanso" and r["onda"].valle.p == pytest.approx(1.0)

    def test_sauer_meyer(self):
        s = gd.sauer_meyer(50.0, ARIA)
        assert s["u"](0.0, 0.0) == 0
        for eta in [0.3, 1.0]:
            assert s["u"](s["xi_sonica"](eta), eta) == pytest.approx(0.0, abs=1e-15)
            assert s["v"](s["xi_gola"](eta), eta) == pytest.approx(0.0, abs=1e-15)

    def test_ugello_minima_lunghezza(self):
        r = gd.ugello_minima_lunghezza(2.5, ARIA, n=40)
        assert r["theta_w_max"] == pytest.approx(gd.prandtl_meyer(2.5, ARIA) / 2)
        assert r["h_uscita_su_h_gola"] == pytest.approx(gd.area_su_area_critica(2.5, ARIA), rel=0.01)

    def test_odografo(self):
        st = gd.State.da_M(ARIA_ADIM, 2.0, 10.0, p=1, a=1)
        o = gd.polare_odografo(st, "C+", n=50)
        us, vs = o["urto"]["us"], o["urto"]["vs"]
        assert us[0] == pytest.approx(o["stato"][0]) and vs[0] == pytest.approx(o["stato"][1])
        # punto d'urto normale: M2* = 1/M1* (Prandtl)
        assert math.hypot(us[-1], vs[-1]) == pytest.approx(1 / st.M_star)
        p = gd.polare_p_theta(st, "C+", n=50)
        assert p["urto"]["p"][-1] == pytest.approx(4.5)
        assert p["punti"]["theta_max"][0] == pytest.approx(10.0 + gd.theta_max(2.0, ARIA)[0])


# =============================================================================
# Aree 8–9 — grafici e Godunov
# =============================================================================
class TestArea89:
    def test_godunov_uniforme(self):
        r, u, p, dt = gd.godunov_passo(np.ones(5), 0.3 * np.ones(5), np.ones(5), ARIA_ADIM, 0.1)
        assert np.allclose(r, 1) and np.allclose(u, 0.3) and np.allclose(p, 1)

    def test_godunov_sod(self):
        N = 400
        x = (np.arange(N) + 0.5) / N - 0.5
        r = np.where(x < 0, 1.0, 0.125)
        u = np.zeros(N)
        p = np.where(x < 0, 1.0, 0.1)
        r, u, p = gd.godunov(r, u, p, ARIA_ADIM, 1.0 / N, 0.2)
        es = gd.riemann(gd.State(ARIA_ADIM, 1, 0, 1), gd.State(ARIA_ADIM, 0.125, 0, 0.1)).campiona(x, 0.2)
        assert np.mean(np.abs(r - es["rho"])) < 0.01
        assert np.mean(np.abs(p - es["p"])) < 0.01

    def test_godunov_parete_conserva_massa(self):
        N = 50
        r = np.linspace(1, 2, N)
        u = np.linspace(-0.3, 0.4, N)
        p = np.linspace(1, 1.5, N)
        r2, u2, p2, _ = gd.godunov_passo(r, u, p, ARIA_ADIM, 0.02, bc=("parete", "parete"))
        assert np.sum(r2) == pytest.approx(np.sum(r), rel=1e-13)
        E = lambda r, u, p: np.sum(p / 0.4 + 0.5 * r * u**2)
        assert E(r2, u2, p2) == pytest.approx(E(r, u, p), rel=1e-13)

    def test_grafici_non_falliscono(self, tmp_path):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        s = gd.riemann(gd.State(ARIA_ADIM, 1, 0, 1), gd.State(ARIA_ADIM, 0.125, 0, 0.1))
        gd.grafico_xt(s, 0.2, (-0.5, 0.5))
        gd.grafico_pu(s)
        gd.grafico_profili(s, 0.2, (-0.5, 0.5))
        st = gd.State.da_M(ARIA_ADIM, 2.0, p=1, a=1, nome="1")
        gd.grafico_p_theta([(st, "entrambe")], [st])
        gd.grafico_odografo([(st, "entrambe")], [st])
        plt.close("all")


# =============================================================================
# Livello 3 — risultati del corso (verifica/casi_corso.py)
# =============================================================================
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "verifica"))
import casi_corso  # noqa: E402

CASI_CORSO = casi_corso.tutti()


class TestLivello3:
    """Ogni caso deve avere l'esito atteso: 'concorda', 'discorda' (causa nota in ``nota``) o
    'non confrontabile'. Una discrepanza nuova o una scomparsa fanno fallire il test."""

    @pytest.mark.parametrize("caso", CASI_CORSO,
                             ids=[f"{c.id}-{i}" for i, c in enumerate(CASI_CORSO)])
    def test_esito(self, caso):
        assert caso.esito == caso.atteso, (
            f"{caso.id} {caso.grandezza}: fonte {caso.fonte}, calcolato {caso.calcolato}, "
            f"diff {caso.diff}, tol {caso.tol}")
        if caso.atteso == "discorda":
            assert caso.nota, "ogni discrepanza deve avere la causa documentata"

    def test_nessuna_discrepanza_senza_causa(self):
        assert all(c.nota for c in CASI_CORSO if c.esito == "discorda")
