"""Test di gasdyn.py.

Livello 1: coerenza interna (inverse, limiti, conservazione).
Livello 2: valori di riferimento classici.
Livello 3: risultati del corso (catalogo e soluzioni_A), in fondo al file.
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
