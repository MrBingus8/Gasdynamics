"""gasdyn — calcoli di Gasdinamica nella lingua del corso (Favini).

Fonte primaria delle formule: ``fonti/formulario.tex`` (citato come "F §n");
fonte secondaria: ``fonti/summary_lectures.tex`` (citato come "S p.n").
Le relazioni che non compaiono in nessuna delle due sono marcate
**fuori formulario** nel docstring.

Convenzioni (vedi anche ``gasdyn_guida.md``):

* δ = (γ−1)/2, q = (γ−1)/(2γ); variabile entropica del corso s = p/ρ^γ.
* 1D: x verso destra, u > 0 verso destra. Famiglie λ1 = u−a (lungo cui è
  costante R1 = a−δu), λ2 = u+a (R2 = a+δu), λ0 = u (contatto).
  Un'onda semplice "u+a" ha R1 uniforme; un'onda "u−a" ha R2 uniforme.
* Urti mobili: v = u−w, M_r = (u−w)/a con segno; 1 = monte (davanti
  all'urto), 2 = valle. Urto u+a: M_r1 < 0; urto u−a: M_r1 > 0.
* Riemann: 1 | u−a | 2 | contatto | 3 | u+a | 4; p2 = p3 = p*, u2 = u3 = u*.
* 2D: σ angolo d'urto, θ direzione del flusso (gradi, antioraria da x),
  μ angolo di Mach, ν Prandtl–Meyer. Onda C+: fronte inclinato θ+μ,
  attraverso essa Δ(θ+ν) = 0 (la compressione aumenta θ). Onda C−: fronte
  θ−μ, Δ(θ−ν) = 0 (la compressione diminuisce θ).
* Tutti gli angoli dell'interfaccia pubblica sono in **gradi**.
* Unità qualsiasi purché coerenti (SI o adimensionali delle tracce).
  T richiede R (``Gas(gamma, R)``); con R = None le temperature non sono
  disponibili.

Dipendenze: numpy, scipy, matplotlib (solo per i grafici).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Callable, Optional, Sequence

import numpy as np
from scipy.optimize import brentq, minimize_scalar

__all__ = []  # popolato in fondo al file

#: 1 atm in Pa (alcuni svolgimenti usano 1e5 Pa: vedere verifica_catalogo.md)
ATM = 101325.0
#: tolleranza relativa di default dei solutori iterativi
TOL = 1e-12
#: tolleranza relativa per dichiarare "nulla" un'onda (sigla N)
TOL_NULLA = 1e-9


# =============================================================================
# Eccezioni
# =============================================================================
class GasdynError(Exception):
    """Errore generico del modulo: dati incoerenti o grandezza non definita."""


class NessunaSoluzione(GasdynError):
    """Il problema posto non ha soluzione fisica (es. M a valle oltre il limite)."""


class UrtoStaccato(GasdynError):
    """Deviazione richiesta maggiore di θmax: urto obliquo attaccato impossibile."""


class Vuoto(GasdynError):
    """L'onda richiesta porterebbe a = 0 (generazione di vuoto)."""


class Strozzamento(GasdynError):
    """Rapporto d'area minore di quello critico: la portata non può passare."""


class RamiMultipli(GasdynError):
    """Più soluzioni ammissibili: sono elencate in ``self.soluzioni`` (etichettate)."""

    def __init__(self, messaggio, soluzioni):
        super().__init__(messaggio)
        self.soluzioni = soluzioni


# =============================================================================
# Utilità
# =============================================================================
def _sind(x):
    return math.sin(math.radians(x))


def _cosd(x):
    return math.cos(math.radians(x))


def _tand(x):
    return math.tan(math.radians(x))


def _fmt(x, n=5):
    try:
        return f"{x:.{n}g}"
    except (TypeError, ValueError):
        return str(x)


def _radici(f, a, b, n=400, scala="lin", xtol=1e-14, rtol=1e-13):
    """Tutte le radici di f in [a, b] cercate per cambio di segno su n punti.

    Ogni radice è raffinata con ``brentq``. Usata dove la soluzione può non
    essere unica, per restituire tutti i rami invece di sceglierne uno.
    """
    if scala == "log":
        xs = np.geomspace(a, b, n)
    else:
        xs = np.linspace(a, b, n)
    fs = []
    for x in xs:
        try:
            fs.append(f(x))
        except (GasdynError, ValueError, ZeroDivisionError):
            fs.append(np.nan)
    out = []
    for i in range(len(xs) - 1):
        f0, f1 = fs[i], fs[i + 1]
        if not (np.isfinite(f0) and np.isfinite(f1)):
            continue
        if f0 == 0.0:
            out.append(xs[i])
        elif f0 * f1 < 0:
            out.append(brentq(f, xs[i], xs[i + 1], xtol=xtol, rtol=rtol))
    if np.isfinite(fs[-1]) and fs[-1] == 0.0:
        out.append(xs[-1])
    return out


# =============================================================================
# 1. Gas, stato, termodinamica e grandezze totali
# =============================================================================
@dataclass(frozen=True)
class Gas:
    """Gas perfetto caloricamente (F §1).

    Parametri
    ---------
    gamma : rapporto dei calori specifici cp/cv (> 1).
    R : costante del gas [J/(kg K)] o None per dati adimensionali
        (in quel caso T, cp, cv non sono disponibili).
    nome : etichetta libera.

    Grandezze derivate: ``delta`` = (γ−1)/2, ``q`` = (γ−1)/(2γ),
    ``cp`` = γR/(γ−1), ``cv`` = R/(γ−1).
    """

    gamma: float
    R: Optional[float] = None
    nome: str = ""

    def __post_init__(self):
        if not self.gamma > 1.0:
            raise GasdynError(f"gamma deve essere > 1 (dato {self.gamma})")
        if self.R is not None and not self.R > 0:
            raise GasdynError(f"R deve essere > 0 (dato {self.R})")

    @classmethod
    def aria(cls, R: float = 287.0, gamma: float = 1.4):
        """Aria come gas perfetto: γ = 1.4, R = 287 J/(kg K) salvo diversa indicazione.

        Alcune tracce usano R = 288 (B73-E01) o 287.05: passare il valore della traccia.
        """
        return cls(gamma, R, "aria")

    @property
    def delta(self):
        """δ = (γ−1)/2 (F §1)."""
        return 0.5 * (self.gamma - 1.0)

    @property
    def q(self):
        """q = (γ−1)/(2γ), esponente delle curve d'espansione (F §8)."""
        return (self.gamma - 1.0) / (2.0 * self.gamma)

    @property
    def cv(self):
        self._serve_R("cv")
        return self.R / (self.gamma - 1.0)

    @property
    def cp(self):
        self._serve_R("cp")
        return self.gamma * self.R / (self.gamma - 1.0)

    def _serve_R(self, cosa):
        if self.R is None:
            raise GasdynError(f"{cosa} richiede R: definire Gas(gamma, R)")

    def __repr__(self):
        r = "" if self.R is None else f", R={_fmt(self.R)}"
        n = f" '{self.nome}'" if self.nome else ""
        return f"Gas{n}(γ={_fmt(self.gamma)}{r})"


#: aria standard del corso (γ = 1.4, R = 287): solo una comodità, mai implicita
ARIA = Gas.aria()


@dataclass(frozen=True)
class State:
    """Stato termodinamico e cinematico uniforme di un gas perfetto.

    Variabili primitive: densità ``rho``, componenti di velocità ``u`` (asse x)
    e ``v`` (asse y, solo 2D), pressione ``p``. ``nome`` è un'etichetta
    (es. "1", "2", "monte").

    Proprietà calcolate (F §1): ``a`` velocità del suono, ``T`` (richiede R),
    ``V`` modulo della velocità, ``M`` = V/a, ``Mu`` = u/a (Mach 1D con segno),
    ``theta`` direzione della velocità in gradi, ``s`` = p/ρ^γ (variabile
    entropica del corso), ``R1`` = a−δu, ``R2`` = a+δu (invarianti 1D, F §3),
    ``h``, ``h0``, ``a0``, ``T0``, ``p0``, ``rho0`` (grandezze totali),
    ``a_star``, ``M_star`` (critiche, F §13), ``V_max``.

    Nota: le grandezze totali dipendono dal riferimento (in 1D non
    stazionario si calcolano con la velocità nel riferimento del laboratorio).
    """

    gas: Gas
    rho: float
    u: float
    p: float
    v: float = 0.0
    nome: str = field(default="", compare=False)

    def __post_init__(self):
        if self.rho < 0 or self.p < 0:
            raise GasdynError(f"ρ e p devono essere non negativi (ρ={self.rho}, p={self.p})")

    # ----- costruttori ------------------------------------------------------
    @classmethod
    def da_pa(cls, gas: Gas, p: float, a: float, u: float = 0.0, v: float = 0.0, nome=""):
        """Stato da pressione e velocità del suono: ρ = γp/a²."""
        return cls(gas, gas.gamma * p / a**2, u, p, v, nome)

    @classmethod
    def da_pT(cls, gas: Gas, p: float, T: float, u: float = 0.0, v: float = 0.0, nome=""):
        """Stato da pressione e temperatura: ρ = p/(RT) (richiede R)."""
        gas._serve_R("da_pT")
        return cls(gas, p / (gas.R * T), u, p, v, nome)

    @classmethod
    def da_M(cls, gas: Gas, M: float, theta: float = 0.0, *, p: float,
             a: float = None, rho: float = None, T: float = None, nome=""):
        """Stato 2D (o 1D con theta = 0 / 180) da Mach, direzione e due grandezze termodinamiche.

        Si passa ``p`` e una fra ``a``, ``rho``, ``T``. ``theta`` in gradi.
        """
        if sum(x is not None for x in (a, rho, T)) != 1:
            raise GasdynError("da_M: specificare esattamente una fra a, rho, T")
        if a is not None:
            rho = gas.gamma * p / a**2
        elif T is not None:
            gas._serve_R("da_M con T")
            rho = p / (gas.R * T)
        a = math.sqrt(gas.gamma * p / rho)
        V = M * a
        return cls(gas, rho, V * _cosd(theta), p, V * _sind(theta), nome)

    @classmethod
    def da_totali(cls, gas: Gas, M: float, *, p0: float, a0: float = None, T0: float = None,
                  theta: float = 0.0, nome=""):
        """Stato da Mach e grandezze totali (p0 e una fra a0, T0), relazioni isentropiche F §1."""
        if (a0 is None) == (T0 is None):
            raise GasdynError("da_totali: specificare esattamente una fra a0, T0")
        if T0 is not None:
            gas._serve_R("da_totali con T0")
            a0 = math.sqrt(gas.gamma * gas.R * T0)
        f = 1.0 + gas.delta * M**2
        a = a0 / math.sqrt(f)
        p = p0 * f ** (-gas.gamma / (gas.gamma - 1.0))
        return cls.da_M(gas, M, theta, p=p, a=a, nome=nome)

    def con(self, **modifiche):
        """Copia dello stato con alcuni campi cambiati (es. ``st.con(u=0, nome='5')``)."""
        return replace(self, **modifiche)

    def ruota(self, dtheta: float):
        """Copia con la velocità ruotata di ``dtheta`` gradi (modulo invariato)."""
        c, s = _cosd(dtheta), _sind(dtheta)
        return replace(self, u=c * self.u - s * self.v, v=s * self.u + c * self.v)

    # ----- proprietà --------------------------------------------------------
    @property
    def gamma(self):
        return self.gas.gamma

    @property
    def a(self):
        if self.rho == 0:
            return 0.0
        return math.sqrt(self.gas.gamma * self.p / self.rho)

    @property
    def T(self):
        self.gas._serve_R("T")
        return self.p / (self.rho * self.gas.R)

    @property
    def V(self):
        return math.hypot(self.u, self.v)

    @property
    def M(self):
        return self.V / self.a

    @property
    def Mu(self):
        """Mach 1D con segno, u/a."""
        return self.u / self.a

    @property
    def theta(self):
        """Direzione della velocità [gradi], atan2(v, u)."""
        return math.degrees(math.atan2(self.v, self.u))

    @property
    def mu(self):
        """Angolo di Mach [gradi] (solo M ≥ 1)."""
        return angolo_mach(self.M)

    @property
    def s(self):
        """Variabile entropica del corso s = p/ρ^γ (F §1)."""
        return self.p / self.rho ** self.gas.gamma

    @property
    def R1(self):
        """Invariante R1 = a − δu, costante lungo λ1 = u − a (F §3)."""
        return self.a - self.gas.delta * self.u

    @property
    def R2(self):
        """Invariante R2 = a + δu, costante lungo λ2 = u + a (F §3)."""
        return self.a + self.gas.delta * self.u

    @property
    def h(self):
        """Entalpia specifica h_i = a²/(γ−1) (F §1)."""
        return self.a**2 / (self.gas.gamma - 1.0)

    @property
    def e(self):
        """Energia interna specifica e_i = p/((γ−1)ρ)."""
        return self.p / ((self.gas.gamma - 1.0) * self.rho)

    @property
    def E(self):
        """Energia totale per unità di volume ρe = ρ(e_i + V²/2)."""
        return self.rho * (self.e + 0.5 * self.V**2)

    @property
    def h0(self):
        return self.h + 0.5 * self.V**2

    @property
    def a0(self):
        """a0² = a² + δV² (F §1)."""
        return math.sqrt(self.a**2 + self.gas.delta * self.V**2)

    @property
    def T0(self):
        return self.T * (1.0 + self.gas.delta * self.M**2)

    @property
    def p0(self):
        return self.p * (1.0 + self.gas.delta * self.M**2) ** (self.gamma / (self.gamma - 1.0))

    @property
    def rho0(self):
        return self.rho * (1.0 + self.gas.delta * self.M**2) ** (1.0 / (self.gamma - 1.0))

    @property
    def a_star(self):
        """a* = a0·√(2/(γ+1)) (F §7, §13)."""
        return self.a0 * math.sqrt(2.0 / (self.gamma + 1.0))

    @property
    def M_star(self):
        """M* = V/a* (F §13)."""
        return self.V / self.a_star

    @property
    def V_max(self):
        """Velocità massima (espansione nel vuoto stazionaria) V_max = a0/√δ."""
        return self.a0 / math.sqrt(self.gas.delta)

    # ----- vettori conservativi --------------------------------------------
    def U(self):
        """Vettore conservativo 1D [ρ, ρu, ρe] (F §2)."""
        return np.array([self.rho, self.rho * self.u, self.E])

    def F(self):
        """Flusso 1D [ρu, ρu²+p, u(ρe+p)] (F §2)."""
        return np.array([self.rho * self.u, self.rho * self.u**2 + self.p,
                         self.u * (self.E + self.p)])

    def __repr__(self):
        n = f"'{self.nome}' " if self.nome else ""
        txt = f"State {n}(ρ={_fmt(self.rho)}, u={_fmt(self.u)}"
        if self.v != 0.0:
            txt += f", v={_fmt(self.v)}"
        txt += f", p={_fmt(self.p)}, a={_fmt(self.a)}"
        if self.rho > 0 and self.a > 0:
            txt += f", M={_fmt(self.M)}"
            if self.v != 0.0:
                txt += f", θ={_fmt(self.theta)}°"
        if self.gas.R is not None and self.rho > 0:
            txt += f", T={_fmt(self.T)}"
        return txt + ")"


def stato_da_U(gas: Gas, U, nome="") -> State:
    """Stato 1D dal vettore conservativo [ρ, ρu, ρe]."""
    rho, m, E = (float(x) for x in U)
    u = m / rho
    p = (gas.gamma - 1.0) * (E - 0.5 * rho * u**2)
    return State(gas, rho, u, p, 0.0, nome)


def isentropica(stato: State, *, p=None, rho=None, a=None, T=None, u=None, v=None, nome="") -> State:
    """Stato sulla stessa isentropica (stessa s = p/ρ^γ) con una grandezza assegnata (F §1).

    La velocità del nuovo stato è quella passata in ``u``/``v`` (default: quella di ``stato``).
    """
    g = stato.gas.gamma
    if sum(x is not None for x in (p, rho, a, T)) != 1:
        raise GasdynError("isentropica: specificare esattamente una fra p, rho, a, T")
    if T is not None:
        stato.gas._serve_R("isentropica con T")
        a = math.sqrt(g * stato.gas.R * T)
    if a is not None:
        if a < 0:
            raise Vuoto("a < 0 richiesto")
        rho = stato.rho * (a / stato.a) ** (2.0 / (g - 1.0))
    elif p is not None:
        rho = stato.rho * (p / stato.p) ** (1.0 / g)
    p_new = stato.s * rho**g
    return State(stato.gas, rho, stato.u if u is None else u, p_new,
                 stato.v if v is None else v, nome)


def salto_entropia(st1: State, st2: State) -> dict:
    """Variazione d'entropia da 1 a 2 (F §1).

    Restituisce ``s2_s1`` = (p2/ρ2^γ)/(p1/ρ1^γ), ``ds_cv`` = Δs_th/cv = ln(s2/s1),
    ``ds_R`` = Δs_th/R = ln(s2/s1)/(γ−1) e, se R è noto, ``ds`` in J/(kg K).
    I due stati devono appartenere allo stesso gas.
    """
    if st1.gas != st2.gas:
        raise GasdynError("salto_entropia: i due stati devono avere lo stesso gas")
    r = st2.s / st1.s
    ds_cv = math.log(r)
    out = {"s2_s1": r, "ds_cv": ds_cv, "ds_R": ds_cv / (st1.gamma - 1.0)}
    if st1.gas.R is not None:
        out["ds"] = st1.gas.cv * ds_cv
    return out


def isentropiche(M, gas: Gas) -> dict:
    """Rapporti totale/statico per Mach M (F §1, §16): T0/T, p0/p, ρ0/ρ, a0/a.

    Accetta scalari o array numpy.
    """
    M = np.asarray(M, dtype=float)
    g = gas.gamma
    f = 1.0 + gas.delta * M**2
    out = {"T0_T": f, "p0_p": f ** (g / (g - 1.0)), "rho0_rho": f ** (1.0 / (g - 1.0)),
           "a0_a": np.sqrt(f)}
    return {k: (float(v) if np.ndim(v) == 0 else v) for k, v in out.items()}


def mach_da_rapporto(valore: float, gas: Gas, tipo: str = "p0/p") -> float:
    """Inversa delle isentropiche: M da p0/p, T0/T, rho0/rho o a0/a (F §1).

    ``valore`` deve essere ≥ 1. Forma chiusa: 1 + δM² = (rapporto)^esponente.
    """
    g = gas.gamma
    esp = {"p0/p": (g - 1.0) / g, "T0/T": 1.0, "rho0/rho": g - 1.0, "a0/a": 2.0}
    if tipo not in esp:
        raise GasdynError(f"tipo deve essere uno fra {list(esp)}")
    if valore < 1.0:
        raise NessunaSoluzione(f"{tipo} = {valore} < 1: nessun Mach reale")
    return math.sqrt((valore ** esp[tipo] - 1.0) / gas.delta)


def pressione_totale(M, gas: Gas, modello: str = "isentropico"):
    """Rapporto p0/p fra pressione totale e statica in tre modelli.

    * ``'isentropico'``: p0/p = (1+δM²)^(γ/(γ−1)) (F §1).
    * ``'incomprimibile'``: p0 = p + ½ρV² ⇒ p0/p = 1 + γM²/2 (**fuori formulario**,
      Bernoulli incomprimibile; B43-E01).
    * ``'isotermo'``: dp/ρ = −V dV con p/ρ = cost ⇒ p0/p = exp(γM²/2)
      (**fuori formulario**, Bernoulli isotermo; B43-E01).
    """
    M = np.asarray(M, dtype=float)
    g = gas.gamma
    if modello == "isentropico":
        r = (1.0 + gas.delta * M**2) ** (g / (g - 1.0))
    elif modello == "incomprimibile":
        r = 1.0 + 0.5 * g * M**2
    elif modello == "isotermo":
        r = np.exp(0.5 * g * M**2)
    else:
        raise GasdynError("modello: 'isentropico', 'incomprimibile' o 'isotermo'")
    return float(r) if np.ndim(r) == 0 else r


def critiche(gas: Gas) -> dict:
    """Rapporti critici (M = 1) rispetto ai totali (F §16): T*/T0, p*/p0, ρ*/ρ0, a*/a0."""
    g = gas.gamma
    t = 2.0 / (g + 1.0)
    return {"T*_T0": t, "p*_p0": t ** (g / (g - 1.0)), "rho*_rho0": t ** (1.0 / (g - 1.0)),
            "a*_a0": math.sqrt(t)}


def velocita_max(a0: float, gas: Gas) -> float:
    """Velocità massima stazionaria V_max = a0·√(2/(γ−1)) = a0/√δ (da a0² = a² + δV² con a = 0)."""
    return a0 / math.sqrt(gas.delta)
