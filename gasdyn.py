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


def _radici(f, a, b, n=400, scala="lin", xtol=1e-14, rtol=1e-13, extra=()):
    """Tutte le radici di f in [a, b] cercate per cambio di segno su n punti.

    Ogni radice è raffinata con ``brentq``. Usata dove la soluzione può non
    essere unica, per restituire tutti i rami invece di sceglierne uno.
    ``extra``: punti aggiunti alla griglia (es. il bordo del dominio di f).
    """
    if scala == "log":
        xs = np.geomspace(a, b, n)
    else:
        xs = np.linspace(a, b, n)
    if len(extra):
        xs = np.unique(np.concatenate([xs, np.asarray(extra, dtype=float)]))
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


# =============================================================================
# 2. Isentropico e quasi-1D (F §16–17)
# =============================================================================
_MACH_SUP_MAX = 1e4  # limite superiore del ramo supersonico nelle inversioni


def area_su_area_critica(M, gas: Gas):
    """Legge area–Mach 𝒜(M) = A/A* = (1/M)[2/(γ+1)(1+δM²)]^((γ+1)/(2(γ−1))) (F §16).

    Accetta scalari o array.
    """
    M = np.asarray(M, dtype=float)
    g = gas.gamma
    r = (1.0 / M) * ((2.0 / (g + 1.0)) * (1.0 + gas.delta * M**2)) ** ((g + 1.0) / (2.0 * (g - 1.0)))
    return float(r) if np.ndim(r) == 0 else r


def mach_da_area(A_su_Astar: float, gas: Gas, ramo: str) -> float:
    """Inversa della legge area–Mach (F §16).

    ``ramo`` = ``'subsonico'`` o ``'supersonico'`` (obbligatorio: per A/A* > 1
    esistono sempre due soluzioni). Solleva :class:`Strozzamento` se A/A* < 1.
    Risolto con ``brentq`` su M ∈ (1e-12, 1) o (1, 1e4), tolleranza 1e-14.
    """
    if A_su_Astar < 1.0 - 1e-14:
        raise Strozzamento(f"A/A* = {A_su_Astar} < 1: la sezione è più piccola di quella critica")
    if A_su_Astar <= 1.0 + 1e-14:
        return 1.0
    f = lambda M: area_su_area_critica(M, gas) - A_su_Astar
    if ramo in ("subsonico", "sub"):
        return brentq(f, 1e-12, 1.0, xtol=1e-15, rtol=1e-14)
    if ramo in ("supersonico", "super"):
        return brentq(f, 1.0, _MACH_SUP_MAX, xtol=1e-14, rtol=1e-14)
    raise GasdynError("ramo deve essere 'subsonico' o 'supersonico'")


def mach_da_area_tra_sezioni(M1: float, A1: float, A2: float, gas: Gas, ramo: str = "stesso") -> float:
    """Mach nella sezione A2 noto M1 in A1, flusso isentropico (S p.48, F §16).

    ``ramo``: ``'stesso'`` (mantiene il regime di M1), ``'subsonico'``, ``'supersonico'``.
    """
    A2s = area_su_area_critica(M1, gas) * A2 / A1
    if ramo == "stesso":
        if abs(M1 - 1.0) < 1e-12:
            raise RamiMultipli("M1 = 1: il ramo a valle è ambiguo, specificare ramo", [])
        ramo = "subsonico" if M1 < 1 else "supersonico"
    return mach_da_area(A2s, gas, ramo)


def portata(A: float, p0: float, T0: float, M: float, gas: Gas) -> float:
    """Portata ṁ = A p0/√T0 √(γ/R) M (1+δM²)^(−(γ+1)/(2(γ−1))) (F §16). Richiede R."""
    gas._serve_R("portata")
    g = gas.gamma
    return A * p0 / math.sqrt(T0) * math.sqrt(g / gas.R) * M * (1.0 + gas.delta * M**2) ** (
        -(g + 1.0) / (2.0 * (g - 1.0)))


def portata_max(A_star: float, p0: float, T0: float, gas: Gas) -> float:
    """Portata critica ṁ_max = A* p0/√T0 √(γ/R) (2/(γ+1))^((γ+1)/(2(γ−1))) (F §16)."""
    return portata(A_star, p0, T0, 1.0, gas)


def portata_per_area(stato: State) -> float:
    """Flusso di massa per unità d'area ρ|V| di uno stato."""
    return stato.rho * stato.V


def ugello_convergente(p0: float, pb: float, gas: Gas) -> dict:
    """Convergente alimentato da serbatoio p0 che scarica in pb (F §16).

    Restituisce ``regime`` ('subsonico' o 'strozzato'), ``Me``, ``pe``.
    Strozzato se pb/p0 ≤ (2/(γ+1))^(γ/(γ−1)); allora Me = 1 e pe = p*.
    """
    pc = critiche(gas)["p*_p0"]
    if pb >= p0:
        return {"regime": "fermo", "Me": 0.0, "pe": p0}
    if pb / p0 <= pc:
        return {"regime": "strozzato", "Me": 1.0, "pe": pc * p0}
    return {"regime": "subsonico", "Me": mach_da_rapporto(p0 / pb, gas, "p0/p"), "pe": pb}


def urto_in_ugello(As_At: float, Ae_At: float, gas: Gas) -> dict:
    """Ugello convergente–divergente strozzato con urto normale in A_s (procedura F §16).

    (i) M1 supersonico da A_s/A_t; (ii) urto normale → M2, p02/p01;
    (iii) A2* = A_s/𝒜(M2); (iv) Me subsonico da A_e/A2*, pe/p01 = (p02/p01)(pe/p02).
    """
    M1 = mach_da_area(As_At, gas, "supersonico")
    un = urto_normale(M1, gas)
    A2s_At = As_At / area_su_area_critica(un["M2"], gas)
    Me = mach_da_area(Ae_At / A2s_At, gas, "subsonico")
    pe_p01 = un["p02_p01"] / isentropiche(Me, gas)["p0_p"]
    return {"As_At": As_At, "M1": M1, "M2": un["M2"], "p02_p01": un["p02_p01"],
            "A2*_At": A2s_At, "Me": Me, "pe_p01": pe_p01}


def ugello_cd(p0: float, pb: float, Ae_At: float, gas: Gas) -> dict:
    """Regime di un ugello convergente–divergente con A_e/A_t assegnato (F §16).

    Pressioni caratteristiche: p_e1 (subsonico con gola sonica), p_e2 (urto
    normale nella sezione d'uscita), p_e3 (progetto supersonico). Regimi:
    'subsonico' (pb > p_e1, gola non sonica), 'urto nel divergente'
    (p_e2 < pb < p_e1, posizione A_s/A_t trovata con brentq), 'sovraespanso'
    (p_e3 < pb < p_e2, urti obliqui esterni), 'adattato' (pb = p_e3),
    'sottoespanso' (pb < p_e3). Tutte le pressioni sono riferite a p0.
    """
    Me_sub = mach_da_area(Ae_At, gas, "subsonico")
    Me_sup = mach_da_area(Ae_At, gas, "supersonico")
    pe1 = p0 / isentropiche(Me_sub, gas)["p0_p"]
    pe3 = p0 / isentropiche(Me_sup, gas)["p0_p"]
    pe2 = pe3 * urto_normale(Me_sup, gas)["p2_p1"]
    out = {"pe1": pe1, "pe2": pe2, "pe3": pe3}
    rtol = 1e-10
    if pb >= p0:
        out.update(regime="fermo", Me=0.0, pe=p0)
    elif pb > pe1 * (1 + rtol):
        Me = mach_da_rapporto(p0 / pb, gas, "p0/p")
        Mt = mach_da_area(area_su_area_critica(Me, gas) / Ae_At, gas, "subsonico")
        out.update(regime="subsonico", Me=Me, pe=pb, Mt=Mt)
    elif pb > pe2 * (1 + rtol):
        f = lambda As: urto_in_ugello(As, Ae_At, gas)["pe_p01"] * p0 - pb
        if abs(pb - pe1) <= rtol * pe1:
            As = 1.0
        else:
            As = brentq(f, 1.0 + 1e-12, Ae_At, xtol=1e-14, rtol=1e-13)
        r = urto_in_ugello(As, Ae_At, gas)
        out.update(regime="urto nel divergente", pe=pb, **r)
    elif abs(pb - pe2) <= rtol * pe2:
        r = urto_in_ugello(Ae_At, Ae_At, gas)
        out.update(regime="urto nella sezione d'uscita", pe=pb, **r)
    elif pb > pe3 * (1 + rtol):
        out.update(regime="sovraespanso", Me=Me_sup, pe=pe3)
    elif pb >= pe3 * (1 - rtol):
        out.update(regime="adattato", Me=Me_sup, pe=pe3)
    else:
        out.update(regime="sottoespanso", Me=Me_sup, pe=pe3)
    return out


def spinta(mdot: float, ue: float, pe: float, pa: float, Ae: float, ingresso: dict = None) -> float:
    """Spinta positiva nel verso del getto (F §17).

    F = ṁ ue + (pe − pa) Ae  [− (ṁ_in u_in + (p_in − pa) A_in) se ``ingresso`` =
    dict(mdot=..., u=..., p=..., A=...) è dato, con velocità positive lungo l'asse].
    """
    F = mdot * ue + (pe - pa) * Ae
    if ingresso is not None:
        F -= ingresso["mdot"] * ingresso["u"] + (ingresso["p"] - pa) * ingresso["A"]
    return F


def variazione_area(stato: State, A1: float, A2: float, ramo: str = "stesso", nome="") -> State:
    """Collegamento stazionario attraverso una variazione d'area (F §9 "Variazione d'area", §16).

    Lo stato ``stato`` attraversa la sezione A1 → A2 conservando ṁ, h0 ed s
    (isentropico, quasi stazionario rispetto al condotto fermo). La direzione
    di u si conserva. ``ramo``: 'stesso' (regime invariato, come negli
    svolgimenti del corso: B72-E02, S p.49), 'subsonico', 'supersonico'.
    Solleva :class:`Strozzamento` se A2 è minore dell'area critica di ``stato``
    (non tutta la portata passa: servono onde riflesse).
    """
    M1 = abs(stato.u) / stato.a
    M2 = mach_da_area_tra_sezioni(M1, A1, A2, stato.gas, ramo)
    return _stato_stazionario_M(stato, M2, nome)


def _stato_stazionario_M(stato: State, M2: float, nome="") -> State:
    """Stato con gli stessi h0 ed s di ``stato``, Mach M2 e stesso verso di u."""
    g = stato.gas.gamma
    f = 1.0 + stato.gas.delta * M2**2
    a2 = stato.a0 / math.sqrt(f)
    p2 = stato.p0 * f ** (-g / (g - 1.0))
    rho2 = (p2 / stato.s) ** (1.0 / g)
    return State(stato.gas, rho2, math.copysign(M2 * a2, stato.u if stato.u != 0 else 1.0), p2,
                 0.0, nome)


def miscelamento(correnti: Sequence, nome: str = "m") -> dict:
    """Stato uniforme a valle del miscelamento di correnti parallele in un condotto a sezione costante (F §17).

    ``correnti``: lista di coppie (stato, area), stesso gas, velocità lungo x.
    Bilanci integrali stazionari con pareti senza attrito e adiabatiche:
    ṁ = Σρu A, I = Σ(p + ρu²)A, ṁ h0 = Σ ρuA h0, sezione finale A = ΣA.
    Con p = (I − ṁu)/A e ρ = ṁ/(uA) l'energia dà una quadratica in u con due
    radici: ramo 'subsonico' e 'supersonico' (entrambe soddisfano i bilanci).
    Restituisce dict con i due stati (None se assenti), Δs/R di ciascuno
    rispetto alla media pesata sulla portata, e 'ammissibili' (Δs ≥ 0).
    Solleva NessunaSoluzione se il discriminante è negativo (portata non smaltibile).
    """
    gas = correnti[0][0].gas
    if any(c[0].gas != gas for c in correnti):
        raise GasdynError("miscelamento: tutte le correnti devono avere lo stesso gas")
    g = gas.gamma
    mdot = sum(st.rho * st.u * A for st, A in correnti)
    I = sum((st.p + st.rho * st.u**2) * A for st, A in correnti)
    H = sum(st.rho * st.u * A * st.h0 for st, A in correnti)
    At = sum(A for _, A in correnti)
    h0, J, k = H / mdot, I / mdot, g / (g - 1.0)
    a2, b2, c2 = 0.5 - k, k * J, -h0
    disc = b2 * b2 - 4 * a2 * c2
    if disc < 0:
        raise NessunaSoluzione("miscelamento: nessuno stato uniforme soddisfa i bilanci (Δ < 0)")
    radici = sorted([(-b2 + s * math.sqrt(disc)) / (2 * a2) for s in (1, -1)])
    s_media = sum(st.rho * st.u * A * math.log(st.s) for st, A in correnti) / mdot
    out = {"mdot": mdot, "I": I, "h0": h0, "A": At, "ammissibili": []}
    for u in radici:
        p = (I - mdot * u) / At
        st = State(gas, mdot / (u * At), u, p, 0.0, nome)
        ramo = "subsonico" if st.M < 1 else "supersonico"
        ds_R = (math.log(st.s) - s_media) / (g - 1.0)
        out[ramo] = st
        out[f"ds_R_{ramo}"] = ds_R
        if ds_R >= -1e-12:
            out["ammissibili"].append(ramo)
    return out


def serbatoio_svuotamento(t, p_i: float, a_i: float, V: float, A: float, gas: Gas,
                          p_a: float = None) -> dict:
    """Svuotamento adiabatico di un serbatoio rigido da un orifizio critico (F §17).

    Con C* = (2/(γ+1))^((γ+1)/(2(γ−1))) e K = A a_i C*/V, finché l'orifizio è
    strozzato: m/m_i = [1+δKt]^(−2/(γ−1)), T/T_i = [1+δKt]^(−2),
    p/p_i = [1+δKt]^(−2γ/(γ−1)). Se ``p_a`` è dato restituisce anche
    ``t_fine_strozzamento``, istante in cui p_a/p = (2/(γ+1))^(γ/(γ−1)).
    Le formule non valgono oltre quell'istante (il dict lo segnala).
    """
    g = gas.gamma
    Cs = (2.0 / (g + 1.0)) ** ((g + 1.0) / (2.0 * (g - 1.0)))
    K = A * a_i * Cs / V
    t = np.asarray(t, dtype=float)
    b = 1.0 + gas.delta * K * t
    out = {"K": K, "m_mi": b ** (-2.0 / (g - 1.0)), "T_Ti": b ** (-2.0),
           "p_pi": b ** (-2.0 * g / (g - 1.0))}
    if p_a is not None:
        p_fine = p_a / critiche(gas)["p*_p0"]
        if p_fine >= p_i:
            out["t_fine_strozzamento"] = 0.0
        else:
            out["t_fine_strozzamento"] = ((p_fine / p_i) ** (-(g - 1.0) / (2.0 * g)) - 1.0) / (
                gas.delta * K)
        out["valido"] = t <= out["t_fine_strozzamento"]
    return {k: (float(v) if np.ndim(v) == 0 and not isinstance(v, bool) else v)
            for k, v in out.items()}


def serbatoio_riempimento(V: float, p_i: float, T_i: float, T0_in: float, gas: Gas, *,
                          m: float = None, p: float = None) -> dict:
    """Riempimento adiabatico da una linea a entalpia totale costante h0,in = cp T0,in (F §17).

    Bilancio m cv T − m_i cv T_i = (m − m_i) h0,in ⇒ m T = m_i T_i + (m − m_i) γ T0,in.
    Dato lo stato finale per massa ``m`` oppure pressione ``p``, restituisce m, T, p.
    """
    gas._serve_R("serbatoio_riempimento")
    g, R = gas.gamma, gas.R
    m_i = p_i * V / (R * T_i)
    if (m is None) == (p is None):
        raise GasdynError("specificare esattamente uno fra m e p")
    if p is not None:
        m = (p * V / R - m_i * T_i) / (g * T0_in) + m_i
    T = (m_i * T_i + (m - m_i) * g * T0_in) / m
    return {"m_i": m_i, "m": m, "T": T, "p": m * R * T / V}


# =============================================================================
# 3. Urto normale stazionario (F §7, §12)
# =============================================================================
def urto_normale(M1: float, gas: Gas) -> dict:
    """Rapporti attraverso un urto normale con Mach normale a monte M1 ≥ 1 (F §7).

    Chiavi: ``p2_p1``, ``rho2_rho1`` (= u1/u2), ``T2_T1``, ``a2_a1``, ``M2``,
    ``p02_p01`` (perdita di pressione totale), ``s2_s1`` (p/ρ^γ), ``ds_R``
    (Δs_th/R), ``p02_p1`` (pressione di ristagno a valle su statica a monte,
    utile per il Pitot), ``u2_u1``.
    Solleva :class:`NessunaSoluzione` per M1 < 1 (urto d'espansione, S p.19).
    """
    if M1 < 1.0 - 1e-12:
        raise NessunaSoluzione(f"M1 = {M1} < 1: l'urto d'espansione viola il secondo principio (S p.19)")
    M1 = max(M1, 1.0)
    g, d = gas.gamma, gas.delta
    m2 = M1 * M1
    p21 = 1.0 + 2.0 * g / (g + 1.0) * (m2 - 1.0)
    r21 = (g + 1.0) * m2 / ((g - 1.0) * m2 + 2.0)
    T21 = p21 / r21
    M2 = math.sqrt((1.0 + d * m2) / (g * m2 - d))
    p0r = r21 ** (g / (g - 1.0)) * ((g + 1.0) / (2.0 * g * m2 - (g - 1.0))) ** (1.0 / (g - 1.0))
    s21 = p21 / r21**g
    return {"M1": M1, "p2_p1": p21, "rho2_rho1": r21, "T2_T1": T21, "a2_a1": math.sqrt(T21),
            "M2": M2, "p02_p01": p0r, "s2_s1": s21, "ds_R": math.log(s21) / (g - 1.0),
            "p02_p1": p21 * (1.0 + d * M2**2) ** (g / (g - 1.0)), "u2_u1": 1.0 / r21}


def urto_normale_stato(stato: State, nome="") -> State:
    """Stato a valle di un urto normale stazionario (urto perpendicolare alla velocità di ``stato``)."""
    r = urto_normale(stato.M, stato.gas)
    k = 1.0 / r["rho2_rho1"]
    return State(stato.gas, stato.rho * r["rho2_rho1"], stato.u * k, stato.p * r["p2_p1"],
                 stato.v * k, nome)


def mach_da_urto(valore: float, gas: Gas, tipo: str = "p2/p1") -> float:
    """M1 dell'urto normale da un rapporto noto (inversa di :func:`urto_normale`).

    ``tipo``: 'p2/p1', 'rho2/rho1', 'T2/T1', 'M2', 'p02/p01'. Le prime tre e M2
    hanno forma chiusa o monotona; p02/p01 è risolto con brentq.
    """
    g, d = gas.gamma, gas.delta
    if tipo == "p2/p1":
        if valore < 1:
            raise NessunaSoluzione("p2/p1 < 1")
        return math.sqrt(1.0 + (g + 1.0) / (2.0 * g) * (valore - 1.0))
    if tipo == "rho2/rho1":
        lim = (g + 1.0) / (g - 1.0)
        if not 1 <= valore < lim:
            raise NessunaSoluzione(f"ρ2/ρ1 deve stare in [1, {lim})")
        return math.sqrt(2.0 * valore / ((g + 1.0) - (g - 1.0) * valore))
    if tipo == "M2":
        lim = math.sqrt((g - 1.0) / (2.0 * g))
        if not lim < valore <= 1:
            raise NessunaSoluzione(f"M2 deve stare in ({lim}, 1]")
        return math.sqrt((1.0 + d * valore**2) / (g * valore**2 - d))
    chiave = {"T2/T1": "T2_T1", "p02/p01": "p02_p01"}.get(tipo)
    if chiave is None:
        raise GasdynError("tipo non riconosciuto")
    f = lambda M: urto_normale(M, gas)[chiave] - valore
    return brentq(f, 1.0, _MACH_SUP_MAX, xtol=1e-14, rtol=1e-14)


# =============================================================================
# 4. Urti mobili 1D (F §7, S p.17–23)
# =============================================================================
def _segno_famiglia(famiglia: str) -> int:
    f = str(famiglia).replace("−", "-").replace(" ", "").lower()
    if f in ("u+a", "+", "2", "destra", "dx"):
        return +1
    if f in ("u-a", "-", "1", "sinistra", "sx"):
        return -1
    raise GasdynError(f"famiglia '{famiglia}' non riconosciuta: usare 'u+a' o 'u-a'")


def _nome_famiglia(sgn: int) -> str:
    return "u+a" if sgn > 0 else "u-a"


@dataclass(frozen=True)
class Urto:
    """Urto mobile 1D (F §7).

    ``monte`` è lo stato davanti all'urto (in cui l'urto entra), ``valle`` quello
    dietro. ``w`` velocità dell'urto nel riferimento del laboratorio.
    Urto della famiglia u+a: si muove verso destra rispetto al gas, monte a
    destra, M_r1 = (u1−w)/a1 < 0. Famiglia u−a: speculare, M_r1 > 0.
    Proprietà: ``Ms`` = |M_r1| (Mach dell'urto), ``Mr_monte``, ``Mr_valle``
    (con segno), ``p2_p1``, ``rho2_rho1``. Posizione: ``x(t) = x0 + w (t − t0)``.
    """

    famiglia: str
    monte: State
    valle: State
    w: float
    tipo: str = "urto"

    @property
    def Mr_monte(self):
        return (self.monte.u - self.w) / self.monte.a

    @property
    def Mr_valle(self):
        return (self.valle.u - self.w) / self.valle.a

    @property
    def Ms(self):
        return abs(self.Mr_monte)

    @property
    def p2_p1(self):
        return self.valle.p / self.monte.p

    @property
    def rho2_rho1(self):
        return self.valle.rho / self.monte.rho

    @property
    def sigla(self):
        """'S' (urto) o 'N' se l'urto ha intensità nulla (Ms = 1)."""
        return "N" if abs(self.Ms - 1.0) < TOL_NULLA else "S"

    @property
    def xi(self):
        """Velocità della discontinuità nel piano x–t (pendenza dx/dt)."""
        return self.w

    def residui_rh(self) -> dict:
        """Residui relativi di massa, quantità di moto ed energia (F §2): [ρv], [ρuv+p], [ρev+pu]."""
        a, b, w = self.monte, self.valle, self.w
        va, vb = a.u - w, b.u - w
        m = (a.rho * va, b.rho * vb)
        q = (a.rho * a.u * va + a.p, b.rho * b.u * vb + b.p)
        e = (a.E * va + a.p * a.u, b.E * vb + b.p * b.u)
        sc = lambda x: abs(x[0] - x[1]) / max(abs(x[0]), abs(x[1]), 1e-300)
        return {"massa": sc(m), "qdm": sc(q), "energia": sc(e)}

    def __repr__(self):
        return (f"Urto {self.famiglia}: w={_fmt(self.w)}, Ms={_fmt(self.Ms)}, "
                f"Mr_monte={_fmt(self.Mr_monte)}, Mr_valle={_fmt(self.Mr_valle)}, "
                f"p2/p1={_fmt(self.p2_p1)}, ρ2/ρ1={_fmt(self.rho2_rho1)}\n"
                f"  monte: {self.monte}\n  valle: {self.valle}")


def _urto_da_Ms(monte: State, sgn: int, Ms: float, nome_valle="") -> Urto:
    """Urto della famiglia di segno ``sgn`` con Mach |M_r1| = Ms (F §7)."""
    if Ms < 1.0 - 1e-12:
        raise NessunaSoluzione(f"|M_r1| = {Ms} < 1: urto d'espansione non ammissibile (S p.19)")
    Ms = max(Ms, 1.0)
    g = monte.gas.gamma
    a1 = monte.a
    w = monte.u + sgn * a1 * Ms
    m2 = Ms * Ms
    r21 = (g + 1.0) * m2 / ((g - 1.0) * m2 + 2.0)
    p21 = 1.0 + 2.0 * g / (g + 1.0) * (m2 - 1.0)
    v1 = monte.u - w
    v2 = v1 / r21
    valle = State(monte.gas, monte.rho * r21, w + v2, monte.p * p21, monte.v, nome_valle)
    return Urto(_nome_famiglia(sgn), monte, valle, w)


def M_valle_limite(gas: Gas) -> float:
    """Limite del Mach assoluto a valle di un urto in gas in quiete per Ms → ∞.

    u2/a2 → √(2/(γ(γ−1))) (= 1.8898 per γ = 1.4). Ricavato dalle relazioni di
    F §7 nel limite d'urto forte (**fuori formulario** come formula esplicita).
    """
    g = gas.gamma
    return math.sqrt(2.0 / (g * (g - 1.0)))


def urto_mobile(monte: State, famiglia: str = None, *, w: float = None, Ms: float = None,
                p_valle: float = None, u_valle: float = None, M_valle: float = None,
                nome_valle: str = "") -> Urto:
    """Urto mobile che entra nello stato ``monte`` (F §7), dato UNO fra:

    * ``w``: velocità dell'urto (la famiglia, se omessa, è dedotta dal segno di w−u1);
    * ``Ms``: Mach dell'urto |M_r1| = |u1−w|/a1;
    * ``p_valle``: pressione a valle, Ms² = 1 + (γ+1)/(2γ)(p2/p1 − 1);
    * ``u_valle``: velocità a valle, Ms = k + √(1+k²), k = (γ+1)/4·|u2−u1|/a1
      (formula di M_rel del formulario, F §8);
    * ``M_valle``: Mach assoluto a valle u2/a2 **con segno**; risolto per
      ricerca di tutte le radici in Ms ∈ [1, 1e4] (brentq). Se le radici sono
      più d'una solleva :class:`RamiMultipli`; se nessuna,
      :class:`NessunaSoluzione` (in gas in quiete |u2/a2| < :func:`M_valle_limite`).

    Lo stato a monte può essere in moto (u1 ≠ 0).
    """
    dati = {"w": w, "Ms": Ms, "p_valle": p_valle, "u_valle": u_valle, "M_valle": M_valle}
    dati = {k: v for k, v in dati.items() if v is not None}
    if len(dati) != 1:
        raise GasdynError("urto_mobile: specificare esattamente uno fra w, Ms, p_valle, u_valle, M_valle")
    g = monte.gas.gamma
    a1 = monte.a
    if w is not None:
        sgn = 1 if w > monte.u else -1
        if famiglia is not None and _segno_famiglia(famiglia) != sgn:
            raise NessunaSoluzione(
                f"w = {w} incompatibile con la famiglia {famiglia} (u1 = {monte.u}): "
                "un urto u+a ha w > u1, un urto u−a ha w < u1")
        return _urto_da_Ms(monte, sgn, abs(monte.u - w) / a1, nome_valle)
    if famiglia is None:
        raise GasdynError("urto_mobile: specificare la famiglia ('u+a' o 'u-a')")
    sgn = _segno_famiglia(famiglia)
    if Ms is not None:
        return _urto_da_Ms(monte, sgn, Ms, nome_valle)
    if p_valle is not None:
        if p_valle < monte.p * (1 - 1e-14):
            raise NessunaSoluzione("p_valle < p_monte: un urto deve comprimere")
        return _urto_da_Ms(monte, sgn, math.sqrt(1.0 + (g + 1.0) / (2.0 * g) * (p_valle / monte.p - 1.0)),
                           nome_valle)
    if u_valle is not None:
        du = sgn * (u_valle - monte.u)
        if du < -1e-14 * max(1.0, abs(monte.u)):
            raise NessunaSoluzione(
                f"u_valle = {u_valle}: un urto {famiglia} richiede "
                f"{'u2 > u1' if sgn > 0 else 'u2 < u1'} (u1 = {monte.u})")
        k = 0.25 * (g + 1.0) * max(du, 0.0) / a1
        return _urto_da_Ms(monte, sgn, k + math.sqrt(1.0 + k * k), nome_valle)
    # M_valle
    f = lambda m: (lambda st: st.u / st.a)(_urto_da_Ms(monte, sgn, m).valle) - M_valle
    radici = _radici(f, 1.0, _MACH_SUP_MAX, n=2000, scala="log")
    if not radici:
        lim = M_valle_limite(monte.gas)
        raise NessunaSoluzione(
            f"nessun urto {famiglia} dà u2/a2 = {M_valle} a valle di {monte}. "
            f"Per Ms→∞ u2/a2 → {sgn * lim:.5f} (limite d'urto forte); con gas in quiete "
            f"|u2/a2| < {lim:.5f}")
    sol = [_urto_da_Ms(monte, sgn, m, nome_valle) for m in radici]
    if len(sol) > 1:
        raise RamiMultipli(f"{len(sol)} urti danno u2/a2 = {M_valle}: Ms = "
                           + ", ".join(_fmt(s.Ms) for s in sol), sol)
    return sol[0]


def urto_riflesso_parete(stato: State, famiglia_riflesso: str, u_parete: float = 0.0,
                         nome_valle: str = "") -> Urto:
    """Urto riflesso da una parete (F §9 "Urto–parete").

    ``stato`` è lo stato dietro l'urto incidente (gas che si muove verso la parete),
    che diventa il monte dell'urto riflesso; a valle u = ``u_parete``.
    Parete a destra ⇒ urto riflesso della famiglia 'u-a'; a sinistra ⇒ 'u+a'.
    Se la condizione di parete richiede un'espansione solleva :class:`GasdynError`
    (usare :func:`onda_verso`).
    """
    sgn = _segno_famiglia(famiglia_riflesso)
    if sgn * (u_parete - stato.u) < 0:
        raise GasdynError("la parete genera un'espansione, non un urto: usare onda_verso")
    return urto_mobile(stato, famiglia_riflesso, u_valle=u_parete, nome_valle=nome_valle)


def natura_salto(A: State, B: State, tol: float = 1e-4) -> dict:
    """Identifica quale onda singola può collegare due stati uniformi (B48-E02).

    Prova in ordine: assenza di salto, contatto (u e p uguali), urto
    (Rankine–Hugoniot con w = [ρu]/[ρ], F §7; monte = stato a entropia minore),
    onda semplice isentropica (stessa s e stesso R1 → onda u+a, stesso R2 → u−a).
    ``tol`` è la tolleranza relativa sui residui (i dati delle tracce sono arrotondati).
    Restituisce un dict con ``tipo`` e, per un urto, l'oggetto :class:`Urto`
    ricostruito e i residui.
    """
    if A.gas != B.gas:
        uguali_pu = abs(A.u - B.u) <= tol * max(1.0, abs(A.u)) and abs(A.p - B.p) <= tol * A.p
        return {"tipo": "contatto (gas diversi)" if uguali_pu else "incompatibile"}
    sc = lambda x, y: abs(x - y) / max(abs(x), abs(y), 1e-300)
    ref_u = max(A.a, B.a)
    du = abs(A.u - B.u) / ref_u
    if du <= tol and sc(A.p, B.p) <= tol:
        if sc(A.rho, B.rho) <= tol:
            return {"tipo": "nessun salto"}
        return {"tipo": "contatto", "w": 0.5 * (A.u + B.u)}
    out = {}
    if sc(A.rho, B.rho) > tol:
        w = (B.rho * B.u - A.rho * A.u) / (B.rho - A.rho)
        monte, valle = (A, B) if A.s < B.s else (B, A)
        Ms = abs(monte.u - w) / monte.a
        if Ms >= 1.0 - tol:
            urto = _urto_da_Ms(monte, 1 if w > monte.u else -1, max(Ms, 1.0))
            res = {"p": sc(urto.valle.p, valle.p), "u": abs(urto.valle.u - valle.u) / ref_u,
                   "rho": sc(urto.valle.rho, valle.rho)}
            out = {"urto_candidato": urto, "residui_urto": res}
            if max(res.values()) <= tol:
                return {"tipo": f"urto {urto.famiglia}", "w": w, "urto": urto,
                        "monte": monte.nome or ("A" if monte is A else "B"), "residui": res}
    if sc(A.s, B.s) <= tol:
        if abs(A.R1 - B.R1) / ref_u <= tol:
            return {"tipo": "onda semplice u+a (isentropica)", **out}
        if abs(A.R2 - B.R2) / ref_u <= tol:
            return {"tipo": "onda semplice u-a (isentropica)", **out}
    return {"tipo": "incompatibile: nessuna onda singola collega i due stati", **out}


def urto_isotermo(M1: float) -> dict:
    """Salti attraverso un urto in flusso **omotermo** (p = a²ρ, a costante): **fuori formulario** (B62-E02).

    Con M1 = |v1|/a (Mach relativo a monte): massa ρ1v1 = ρ2v2 e quantità di moto
    p + ρv² ⇒ ρ2/ρ1 = p2/p1 = M1², v2/v1 = 1/M1², relazione di Prandtl v1v2 = a²,
    M2 = 1/M1. L'energia non è conservata (serve un flusso di calore). Per
    M1 → ∞ il rapporto di densità è illimitato, a differenza del caso adiabatico
    ((γ+1)/(γ−1)).
    """
    if M1 < 1:
        raise NessunaSoluzione("M1 < 1")
    return {"rho2_rho1": M1**2, "p2_p1": M1**2, "v2_v1": 1.0 / M1**2, "M2": 1.0 / M1}


# =============================================================================
# 5. Onde semplici 1D, pistoni, formazione dell'urto (F §3–6, S p.6–16)
# =============================================================================
def _stato_da_a_s(gas: Gas, a: float, s: float, u: float, nome="") -> State:
    """Stato da velocità del suono, variabile entropica s = p/ρ^γ e velocità."""
    g = gas.gamma
    if a <= 0:
        return State(gas, 0.0, u, 0.0, 0.0, nome)
    rho = (a * a / (g * s)) ** (1.0 / (g - 1.0))
    return State(gas, rho, u, s * rho**g, 0.0, nome)


def punto_caratteristiche(A: State, B: State, nome: str = "P") -> State:
    """Stato nel punto P in cui la caratteristica u+a uscente da A incontra la u−a uscente da B (F §3).

    a_P = (a_A + δu_A + a_B − δu_B)/2, u_P = (a_A − a_B + δ(u_A + u_B))/(2δ).
    Ipotesi omoentropica: p e ρ si ricavano con la s di A (se s_A ≠ s_B il
    risultato termodinamico vale solo nell'ipotesi del corso e va controllato).
    """
    d = A.gas.delta
    aP = 0.5 * (A.R2 + B.R1)
    uP = (A.R2 - B.R1) / (2.0 * d)
    return _stato_da_a_s(A.gas, aP, A.s, uP, nome)


def curva_onda(stato: State, famiglia: str, p):
    """Curva u(p) degli stati raggiungibili con un'onda della famiglia data (F §8).

    ``famiglia='u-a'``: curva sinistra u_L(p) generata dallo stato 1 a sinistra;
    ``famiglia='u+a'``: curva destra u_R(p) generata dallo stato 4 a destra.
    Per p ≤ p_i espansione isentropica, per p > p_i urto (Hugoniot), con
    A_i = 2/((γ+1)ρ_i), B_i = (γ−1)/(γ+1) p_i. Accetta array.
    """
    sgn = _segno_famiglia(famiglia)
    g = stato.gas.gamma
    p = np.asarray(p, dtype=float)
    Ai = 2.0 / ((g + 1.0) * stato.rho)
    Bi = (g - 1.0) / (g + 1.0) * stato.p
    with np.errstate(invalid="ignore", divide="ignore"):
        R = stato.a / stato.gas.delta * (1.0 - (p / stato.p) ** stato.gas.q)
        S = -(p - stato.p) * np.sqrt(Ai / (p + Bi))
    du = np.where(p <= stato.p, R, S)
    u = stato.u - sgn * du
    return float(u) if np.ndim(u) == 0 else u


def onda_semplice(stato: State, famiglia: str, *, u: float = None, p: float = None,
                  a: float = None, rho: float = None, nome: str = "") -> State:
    """Stato collegato a ``stato`` da un'onda semplice isentropica della famiglia data (F §4).

    Onda u+a: R1 = a − δu uniforme (da = δ du); onda u−a: R2 = a + δu uniforme
    (da = −δ du). Si assegna una fra ``u``, ``p``, ``a``, ``rho``; il resto segue
    dal completamento isentropico p/p_i = (a/a_i)^(2γ/(γ−1)). Vale sia per
    espansioni sia per compressioni (finché non si forma l'urto).
    Solleva :class:`Vuoto` se a ≤ 0.
    """
    sgn = _segno_famiglia(famiglia)
    gas = stato.gas
    d = gas.delta
    if sum(x is not None for x in (u, p, a, rho)) != 1:
        raise GasdynError("onda_semplice: specificare esattamente una fra u, p, a, rho")
    if p is not None:
        a = stato.a * (p / stato.p) ** gas.q
    elif rho is not None:
        a = stato.a * (rho / stato.rho) ** d
    if a is not None:
        u = stato.u + sgn * (a - stato.a) / d
    else:
        a = stato.a + sgn * d * (u - stato.u)
    if a <= 0:
        raise Vuoto(f"onda {famiglia}: a = {a:.5g} ≤ 0, il gas non può seguire (vuoto per "
                    f"u {'≤' if sgn > 0 else '≥'} {stato.u - sgn * stato.a / d:.6g})")
    return _stato_da_a_s(gas, a, stato.s, u, nome)


@dataclass(frozen=True)
class Ventaglio:
    """Onda semplice centrata (ventaglio) nel piano x–t (F §5).

    ``monte``: stato indisturbato in cui avanza la testa; ``valle``: stato dietro
    la coda. Famiglia u+a: monte a destra; u−a: monte a sinistra.
    Pendenze: testa λ(monte), coda λ(valle) con λ = u ± a. Dentro il ventaglio,
    con ξ = x/t (F §5): u+a → u = (ξ − R1)/(1+δ), a = (R1 + δξ)/(1+δ);
    u−a → u = (ξ + R2)/(1+δ), a = (R2 − δξ)/(1+δ).
    ``vuoto=True`` se la coda raggiunge a = 0 (stato valle di vuoto).
    ``compressione=True`` per una compressione isentropica (solo nel modello
    isentropico: le caratteristiche convergono e il ventaglio non è campionabile).
    """

    famiglia: str
    monte: State
    valle: State
    vuoto: bool = False
    compressione: bool = False

    @property
    def tipo(self):
        return "compressione isentropica" if self.compressione else "espansione"

    @property
    def sigla(self):
        if abs(self.valle.p / self.monte.p - 1.0) < TOL_NULLA:
            return "N"
        return "R" if not self.compressione else "S"

    @property
    def _sgn(self):
        return _segno_famiglia(self.famiglia)

    @property
    def xi_testa(self):
        return self.monte.u + self._sgn * self.monte.a

    @property
    def xi_coda(self):
        if self.vuoto:
            return self.valle.u
        return self.valle.u + self._sgn * self.valle.a

    @property
    def transonico(self):
        """True se il ventaglio contiene ξ = 0 (punto sonico fermo in x0)."""
        return self.xi_testa * self.xi_coda < 0

    def stato_xi(self, xi: float, nome="") -> State:
        """Stato in ξ = (x−x0)/(t−t0) (fuori dal ventaglio restituisce monte o valle)."""
        if self.compressione:
            raise GasdynError("compressione isentropica: il ventaglio non è a un sol valore")
        sgn = self._sgn
        if sgn * (xi - self.xi_testa) >= 0:
            return self.monte
        if sgn * (xi - self.xi_coda) <= 0:
            return self.valle
        d = self.monte.gas.delta
        if sgn > 0:
            K = self.monte.R1
            u, a = (xi - K) / (1.0 + d), (K + d * xi) / (1.0 + d)
        else:
            K = self.monte.R2
            u, a = (xi + K) / (1.0 + d), (K - d * xi) / (1.0 + d)
        return _stato_da_a_s(self.monte.gas, a, self.monte.s, u, nome)

    def __repr__(self):
        v = " (VUOTO in coda)" if self.vuoto else ""
        t = " transonico" if self.transonico else ""
        return (f"Ventaglio {self.famiglia} [{self.tipo}{t}]{v}: testa ξ={_fmt(self.xi_testa)}, "
                f"coda ξ={_fmt(self.xi_coda)}, p2/p1={_fmt(self.valle.p / self.monte.p)}\n"
                f"  monte: {self.monte}\n  valle: {self.valle}")


def onda_verso(stato: State, famiglia: str, *, u: float = None, p: float = None,
               consenti_vuoto: bool = False, modello: str = "esatto", nome: str = ""):
    """Onda della famiglia data che porta ``stato`` (monte, indisturbato) a una velocità o pressione assegnata.

    È la condizione al contorno del corso (F §9): parete (u = 0), pistone
    (u = u_p), pressione imposta (p = p_b). Restituisce un :class:`Urto` se
    l'onda comprime (modello 'esatto'), un :class:`Ventaglio` se espande.
    Con ``modello='isentropico'`` anche la compressione è isentropica
    (Ventaglio con compressione=True).
    Famiglia u+a: ``stato`` sta a destra dell'onda (es. gas a destra di un
    pistone o di un contorno sinistro); u−a: a sinistra.
    Se l'espansione richiesta supera il vuoto solleva :class:`Vuoto`, salvo
    ``consenti_vuoto=True`` (ventaglio fino ad a = 0, il gas si stacca).
    """
    sgn = _segno_famiglia(famiglia)
    if (u is None) == (p is None):
        raise GasdynError("onda_verso: specificare esattamente uno fra u e p")
    if u is not None:
        compr = sgn * (u - stato.u) > 0
    else:
        compr = p > stato.p
    if compr and modello == "esatto":
        if u is not None:
            return urto_mobile(stato, famiglia, u_valle=u, nome_valle=nome)
        return urto_mobile(stato, famiglia, p_valle=p, nome_valle=nome)
    try:
        valle = onda_semplice(stato, famiglia, u=u, p=p, nome=nome)
        return Ventaglio(_nome_famiglia(sgn), stato, valle, compressione=compr)
    except Vuoto:
        if not consenti_vuoto:
            raise
        u_vuoto = stato.u - sgn * stato.a / stato.gas.delta
        return Ventaglio(_nome_famiglia(sgn), stato, State(stato.gas, 0.0, u_vuoto, 0.0, 0.0, nome),
                         vuoto=True)


def pistone_impulsivo(stato0: State, u_p: float, lato: str = "sinistro", nome: str = "1"):
    """Pistone messo in moto impulsivamente con velocità ``u_p`` (F §5, §7; S p.9, p.19).

    ``lato='sinistro'``: pistone a sinistra, gas a destra ⇒ onda u+a;
    ``'destro'``: gas a sinistra ⇒ onda u−a. Compressione ⇒ :class:`Urto`
    (w da w² − (1+δ)v_p w − a0² = 0 per gas in quiete); espansione ⇒
    :class:`Ventaglio` centrato. Se il pistone è più veloce della velocità di
    fuga (|u_p − u0| ≥ a0/δ) il ventaglio arriva al vuoto (``vuoto=True``) e
    il pistone si stacca dal gas.
    """
    fam = "u+a" if lato.startswith("sin") else "u-a"
    return onda_verso(stato0, fam, u=u_p, consenti_vuoto=True, nome=nome)


def _derivata(f, t, h=None):
    h = h or 1e-5 * max(1.0, abs(t))
    return (f(t + h) - f(t - h)) / (2.0 * h)


def _derivata2(f, t, h=None):
    h = h or 1e-4 * max(1.0, abs(t))
    return (f(t + h) - 2.0 * f(t) + f(t - h)) / (h * h)


@dataclass
class PistoneLegge:
    """Campo d'onda semplice generato da un pistone con legge x_p(t) (S p.15–16, F §6).

    Ogni caratteristica della famiglia dell'onda parte dalla faccia del pistone
    all'istante t_p con lo stato del pistone (u = ẋ_p, a dall'invariante) e
    pendenza λ(t_p) costante. Il blow-up è il primo punto dell'inviluppo:
    t̂(t_p) = t_p + (λ − ẋ_p)/λ', λ' = (1+δ)ẍ_p (S p.15), minimizzato su t_p.
    Attributi: ``t_bu``, ``x_bu``, ``tp_bu`` (istante di emissione della
    caratteristica che coalesce per prima; inf se non c'è compressione).
    """

    x_p: Callable
    dx_p: Callable
    ddx_p: Callable
    stato0: State
    famiglia: str
    t_max: float
    t_bu: float = math.inf
    x_bu: float = math.nan
    tp_bu: float = math.nan

    @property
    def _sgn(self):
        return _segno_famiglia(self.famiglia)

    def stato_pistone(self, tp: float, nome="") -> State:
        """Stato del gas sulla faccia del pistone all'istante tp (onda semplice dallo stato0)."""
        return onda_semplice(self.stato0, self.famiglia, u=float(self.dx_p(tp)), nome=nome)

    def lam(self, tp: float) -> float:
        """Pendenza dx/dt della caratteristica emessa in tp."""
        st = self.stato_pistone(tp)
        return st.u + self._sgn * st.a

    def t_inviluppo(self, tp: float) -> float:
        """t̂(tp) in cui la caratteristica emessa in tp tocca l'inviluppo (inf se non converge)."""
        st = self.stato_pistone(tp)
        dl = (1.0 + st.gas.delta) * float(self.ddx_p(tp))
        num = self._sgn * st.a  # λ − ẋ
        if dl == 0 or num / dl <= 0:
            return math.inf
        return tp + num / dl

    def caratteristica(self, tp: float):
        """(x0, t0, λ) della caratteristica emessa dal pistone in tp: x = x0 + λ(t − t0)."""
        return float(self.x_p(tp)), tp, self.lam(tp)

    def x_testa(self, t: float) -> float:
        x0, t0, l = self.caratteristica(0.0)
        return x0 + l * (t - t0)

    def stato(self, x: float, t: float, nome="") -> State:
        """Stato in (x, t) per t < t_bu: indisturbato oltre la testa, altrimenti
        lungo la caratteristica che passa per (x, t) (tp trovato con brentq)."""
        if t > self.t_bu * (1 + 1e-12):
            raise GasdynError(f"t = {t} > t_BU = {self.t_bu}: oltre il blow-up le relazioni isentropiche non valgono")
        sgn = self._sgn
        if sgn * (x - self.x_testa(t)) >= 0:
            return self.stato0
        if sgn * (x - float(self.x_p(t))) < 0:
            raise GasdynError(f"x = {x} è dietro il pistone (x_p = {float(self.x_p(t))})")
        tmax = min(t, self.t_max)
        g = lambda tp: (lambda c: c[0] + c[2] * (t - c[1]))(self.caratteristica(tp)) - x
        if t > self.t_max and sgn * g(tmax) < 0:
            raise GasdynError("(x, t) raggiunto da caratteristiche emesse dopo t_max")
        tp = brentq(g, 0.0, tmax, xtol=1e-14, rtol=1e-13)
        return self.stato_pistone(tp, nome)


def pistone_legge(x_p: Callable, stato0: State, t_max: float, lato: str = "sinistro",
                  dx_p: Callable = None, ddx_p: Callable = None, n: int = 4001) -> PistoneLegge:
    """Onda semplice generata da un pistone con legge del moto qualsiasi x_p(t), t ∈ [0, t_max] (S p.15–16).

    ``lato='sinistro'``: pistone a sinistra del gas (onda u+a); ``'destro'``: u−a.
    ``stato0``: gas uniforme iniziale (anche in moto). Le derivate, se non date,
    si calcolano per differenze centrate. Il blow-up si cerca minimizzando
    t̂(t_p) su n punti e raffinando con ``minimize_scalar`` (tolleranza 1e-12).
    Per pistone partito da fermo con ẍ(0) > 0 e gas in quiete si ritrova
    t_BU = a0/((1+δ)ẍ(0)), x_BU = a0 t_BU (F §6).
    """
    fam = "u+a" if lato.startswith("sin") else "u-a"
    dx = dx_p or (lambda t: _derivata(x_p, t))
    ddx = ddx_p or (lambda t: _derivata2(x_p, t))
    P = PistoneLegge(x_p, dx, ddx, stato0, fam, t_max)
    tps = np.linspace(0.0, t_max, n)
    th = np.array([P.t_inviluppo(t) for t in tps])
    if np.all(~np.isfinite(th)):
        return P
    i = int(np.nanargmin(np.where(np.isfinite(th), th, np.nan)))
    lo, hi = tps[max(i - 1, 0)], tps[min(i + 1, n - 1)]
    if hi > lo:
        r = minimize_scalar(P.t_inviluppo, bounds=(lo, hi), method="bounded",
                            options={"xatol": 1e-13 * max(1.0, t_max)})
        tp, tb = (r.x, r.fun) if r.fun <= th[i] else (tps[i], th[i])
    else:
        tp, tb = tps[i], th[i]
    x0, t0, l = P.caratteristica(tp)
    P.t_bu, P.tp_bu, P.x_bu = float(tb), float(tp), float(x0 + l * (tb - t0))
    return P


def t_bu_lineare(x1: float, v1: float, x2: float, v2: float) -> float:
    """Blow-up per profilo lineare di v = u ± a su [x1, x2]: t_BU = |x2−x1|/|v(x2)−v(x1)| (F §6).

    Restituisce inf se v non decresce verso valle (nessuna coalescenza).
    """
    dv = (v2 - v1) / (x2 - x1)
    return math.inf if dv >= 0 else -1.0 / dv


def formazione_urto_ci(u0, a0, gas: Gas, intervallo=None, x=None, n: int = 20001) -> dict:
    """Tempo e posizione di formazione dell'urto da condizioni iniziali (F §6, S p.13–15).

    Dati: ``u0``, ``a0`` come funzioni di x su ``intervallo`` = (x1, x2) oppure
    come array campionati su ``x``. Per ciascuna famiglia (v = u+a e v = u−a)
    calcola t_BU = −1/min v'(ξ), ξ*, x_BU = ξ* + v(ξ*) t_BU.
    Le formule sono esatte per un'onda semplice (l'invariante dell'altra
    famiglia uniforme); per ciascuna famiglia il risultato riporta
    ``onda_semplice`` (True/False) e, se False, va considerato una stima.
    Restituisce {'u+a': {...}, 'u-a': {...}, 'prima': famiglia con t_BU minore}.
    """
    d = gas.delta
    if x is None:
        if intervallo is None:
            raise GasdynError("dare intervallo=(x1, x2) con funzioni, oppure x con array")
        xs = np.linspace(intervallo[0], intervallo[1], n)
        U = np.array([u0(t) for t in xs], dtype=float)
        A = np.array([a0(t) for t in xs], dtype=float)
    else:
        xs = np.asarray(x, float)
        U, A = np.asarray(u0, float), np.asarray(a0, float)
    out = {}
    for fam, sgn in (("u+a", 1), ("u-a", -1)):
        v = U + sgn * A
        dv = np.gradient(v, xs)
        altro = A - sgn * d * U  # R1 per u+a, R2 per u−a
        scala = max(np.max(np.abs(A)), 1e-300)
        semplice = bool((altro.max() - altro.min()) / scala < 1e-6)
        i = int(np.argmin(dv))
        if dv[i] >= 0:
            out[fam] = {"t_bu": math.inf, "x_bu": math.nan, "xi": math.nan, "onda_semplice": semplice}
            continue
        xi, dvm = xs[i], dv[i]
        if x is None and 0 < i < len(xs) - 1:
            vf = lambda t: u0(t) + sgn * a0(t)
            r = minimize_scalar(lambda t: _derivata(vf, t, 1e-6 * max(1.0, abs(t))),
                                bounds=(xs[i - 1], xs[i + 1]), method="bounded",
                                options={"xatol": 1e-12})
            if r.fun < dvm:
                xi, dvm = r.x, r.fun
        vxi = float(np.interp(xi, xs, v)) if x is not None else u0(xi) + sgn * a0(xi)
        tb = -1.0 / dvm
        out[fam] = {"t_bu": tb, "x_bu": xi + vxi * tb, "xi": xi, "v_xi": vxi,
                    "dv_min": dvm, "onda_semplice": semplice}
    out["prima"] = min(("u+a", "u-a"), key=lambda f: out[f]["t_bu"])
    return out


# =============================================================================
# 6. Problema di Riemann 1D (F §8–9)
# =============================================================================
def _lato_sinistro(uk: float, st1: State):
    """Onda u−a in funzione di u^k (F §8): (p2, a2, dp2/du, M_rel,1 o None)."""
    g, d = st1.gas.gamma, st1.gas.delta
    a1, p1 = st1.a, st1.p
    if uk < st1.u:  # urto
        k = 0.25 * (g + 1.0) * (st1.u - uk) / a1
        M = k + math.sqrt(1.0 + k * k)
        p2 = p1 * (1.0 + 2.0 * g / (g + 1.0) * (M * M - 1.0))
        a2 = a1 * math.sqrt((g + 1.0 + (g - 1.0) * p2 / p1) / (g + 1.0 + (g - 1.0) * p1 / p2))
        dp = -2.0 * g * p1 / a1 * M**3 / (1.0 + M * M)
        return p2, a2, dp, M
    a2 = a1 - d * (uk - st1.u)  # espansione
    if a2 <= 0:
        return 0.0, 0.0, 0.0, None
    p2 = p1 * (a2 / a1) ** (2.0 * g / (g - 1.0))
    return p2, a2, -g * p2 / a2, None


def _lato_destro(uk: float, st4: State):
    """Onda u+a in funzione di u^k (F §8): (p3, a3, dp3/du, M_rel,4 o None)."""
    g, d = st4.gas.gamma, st4.gas.delta
    a4, p4 = st4.a, st4.p
    if uk > st4.u:  # urto
        k = 0.25 * (g + 1.0) * (uk - st4.u) / a4
        M = k + math.sqrt(1.0 + k * k)
        p3 = p4 * (1.0 + 2.0 * g / (g + 1.0) * (M * M - 1.0))
        a3 = a4 * math.sqrt((g + 1.0 + (g - 1.0) * p3 / p4) / (g + 1.0 + (g - 1.0) * p4 / p3))
        dp = 2.0 * g * p4 / a4 * M**3 / (1.0 + M * M)
        return p3, a3, dp, M
    a3 = a4 + d * (uk - st4.u)  # espansione
    if a3 <= 0:
        return 0.0, 0.0, 0.0, None
    p3 = p4 * (a3 / a4) ** (2.0 * g / (g - 1.0))
    return p3, a3, g * p3 / a3, None


@dataclass
class SoluzioneRiemann:
    """Soluzione del problema di Riemann 1D (F §8).

    ``stati``: dict {'1','2','3','4'} (in caso di vuoto 2 e 3 sono gli stati di
    vuoto al bordo dei ventagli). ``onda_sx`` (famiglia u−a) e ``onda_dx`` (u+a):
    :class:`Urto` o :class:`Ventaglio`. ``configurazione``: sigla del corso
    (es. 'RCS', 'SCR', 'NCS', 'RNS', 'RVR'). ``iterazioni``: lista di tuple
    (k, u^k, p2^k, p3^k) di Newton (F §8). ``metodo``: 'newton', 'brentq'
    (rete di sicurezza), 'chiuso' (modello isentropico, gas uguali) o 'vuoto'.
    ``x0``, ``t0``: origine della discontinuità iniziale.
    """

    stati: dict
    p_star: float
    u_star: float
    onda_sx: object
    onda_dx: object
    configurazione: str
    vuoto: bool = False
    modello: str = "esatto"
    metodo: str = "newton"
    iterazioni: list = field(default_factory=list)
    x0: float = 0.0
    t0: float = 0.0

    @property
    def contatto(self):
        return self.configurazione[1] == "C"

    def velocita(self) -> dict:
        """Pendenze dx/dt di tutte le linee della soluzione (F §8, tabella testa/coda)."""
        v = {}
        for lato, onda in (("sx", self.onda_sx), ("dx", self.onda_dx)):
            if isinstance(onda, Urto):
                v[f"urto_{lato}"] = onda.w
            else:
                v[f"testa_{lato}"] = onda.xi_testa
                v[f"coda_{lato}"] = onda.xi_coda
        if not self.vuoto:
            v["contatto"] = self.u_star
        return v

    def stato_xi(self, xi: float) -> State:
        """Stato in ξ = (x − x0)/(t − t0)."""
        if self.vuoto:
            if xi <= self.onda_sx.xi_coda:
                return self.onda_sx.stato_xi(xi)
            if xi >= self.onda_dx.xi_coda:
                return self.onda_dx.stato_xi(xi)
            return State(self.stati["1"].gas, 0.0, xi, 0.0, 0.0, "vuoto")
        onda = self.onda_sx if xi < self.u_star else self.onda_dx
        if isinstance(onda, Urto):
            if onda is self.onda_sx:
                return onda.monte if xi < onda.w else onda.valle
            return onda.monte if xi > onda.w else onda.valle
        return onda.stato_xi(xi)

    def campiona(self, x, t) -> dict:
        """Campiona ρ, u, p, a, M, s, p0 in x (array) all'istante t > t0."""
        if t <= self.t0:
            raise GasdynError("campiona: serve t > t0")
        x = np.atleast_1d(np.asarray(x, float))
        sts = [self.stato_xi((xx - self.x0) / (t - self.t0)) for xx in x]
        out = {k: np.array([getattr(s, k) for s in sts]) for k in ("rho", "u", "p")}
        out["a"] = np.array([s.a for s in sts])
        with np.errstate(invalid="ignore", divide="ignore"):
            out["M"] = np.where(out["a"] > 0, out["u"] / np.where(out["a"] > 0, out["a"], 1), np.nan)
            out["s"] = np.array([s.s if s.rho > 0 else np.nan for s in sts])
            out["p0"] = np.array([s.p0 if s.rho > 0 else np.nan for s in sts])
        out["x"] = x
        return out

    def __repr__(self):
        r = [f"Riemann {self.configurazione} (modello {self.modello}, metodo {self.metodo}, "
             f"{len(self.iterazioni)} iterazioni): p* = {_fmt(self.p_star, 7)}, u* = {_fmt(self.u_star, 7)}"]
        for k in ("1", "2", "3", "4"):
            r.append(f"  {k}: {self.stati[k]}")
        for nome, o in (("sinistra", self.onda_sx), ("destra", self.onda_dx)):
            if isinstance(o, Urto):
                r.append(f"  onda {nome} ({o.famiglia}): urto w = {_fmt(o.w, 7)}, Ms = {_fmt(o.Ms, 7)}")
            else:
                r.append(f"  onda {nome} ({o.famiglia}): {o.tipo}, testa {_fmt(o.xi_testa, 7)}, "
                         f"coda {_fmt(o.xi_coda, 7)}")
        if not self.vuoto:
            r.append(f"  contatto: w = {_fmt(self.u_star, 7)}")
        return "\n".join(r)


def _innesco_riemann(st1: State, st4: State) -> float:
    """u⁰ della doppia espansione (F §8). Gas uguali: forma chiusa con z; gas diversi:
    doppia espansione risolta con brentq (estensione: serve solo come punto di partenza)."""
    g1, g4 = st1.gas, st4.gas
    if g1.gamma == g4.gamma:
        d = g1.delta
        z = (st4.a / st1.a) * (st1.p / st4.p) ** ((g1.gamma - 1.0) / (2.0 * g1.gamma))
        return (z * (st1.a + d * st1.u) - (st4.a - d * st4.u)) / (d * (1.0 + z))
    lo = st4.u - st4.a / g4.delta
    hi = st1.u + st1.a / g1.delta
    f = lambda u: (st1.p * max(1.0 - g1.delta * (u - st1.u) / st1.a, 0.0) ** (1.0 / g1.q)
                   - st4.p * max(1.0 + g4.delta * (u - st4.u) / st4.a, 0.0) ** (1.0 / g4.q))
    return brentq(f, lo, hi, xtol=1e-14)


def riemann(sinistra: State, destra: State, *, modello: str = "esatto", tol: float = TOL,
            maxit: int = 50, x0: float = 0.0, t0: float = 0.0) -> SoluzioneRiemann:
    """Problema di Riemann 1D esatto con l'algoritmo del formulario (F §8).

    Stato 1 = ``sinistra``, stato 4 = ``destra`` (anche con gas diversi).
    Algoritmo: innesco u⁰ (doppia espansione, formula con z); per ogni u^k si
    calcolano p2^k, ṗ2^k (onda u−a: urto se u^k < u1, espansione altrimenti) e
    p3^k, ṗ3^k (onda u+a: urto se u^k > u4); Newton
    u^(k+1) = u^k − (p2−p3)/(ṗ2−ṗ3) fino a |1 − p2/p3| < ``tol``.
    Smorzamento: se u^(k+1) esce dall'intervallo fisico (a2, a3 > 0) il passo
    viene dimezzato. Rete di sicurezza (non nel formulario): se Newton non
    converge in ``maxit`` iterazioni si usa ``brentq`` su p2(u) − p3(u), che è
    strettamente decrescente; la soluzione è marcata ``metodo='brentq'``.
    Vuoto se u4 − u1 ≥ a1/δ1 + a4/δ4 (F §8 con gas uguali: (a1+a4)/δ) ⇒ 'RVR'.
    ``modello='isentropico'``: entrambe le onde isentropiche (compressioni
    comprese); con gas uguali è la soluzione chiusa u⁰ del primo passo.
    """
    st1, st4 = sinistra.con(nome="1"), destra.con(nome="4")
    g1, g4 = st1.gas, st4.gas
    lo = st4.u - st4.a / g4.delta
    hi = st1.u + st1.a / g1.delta
    if lo >= hi:  # vuoto
        V1 = onda_verso(st1, "u-a", u=hi, consenti_vuoto=True, nome="2")
        V4 = onda_verso(st4, "u+a", u=lo, consenti_vuoto=True, nome="3")
        V1 = replace(V1, vuoto=True)
        V4 = replace(V4, vuoto=True)
        return SoluzioneRiemann({"1": st1, "2": V1.valle, "3": V4.valle, "4": st4}, 0.0, math.nan,
                                V1, V4, "RVR", vuoto=True, modello=modello, metodo="vuoto",
                                x0=x0, t0=t0)
    iterazioni = []
    if modello == "isentropico":
        u = _innesco_riemann(st1, st4)
        p = st1.p * (1.0 - g1.delta * (u - st1.u) / st1.a) ** (1.0 / g1.q)
        metodo = "chiuso" if g1.gamma == g4.gamma else "brentq"
    elif modello == "esatto":
        u = _innesco_riemann(st1, st4)
        metodo = "newton"
        conv = False
        for k in range(maxit):
            p2, a2, dp2, _ = _lato_sinistro(u, st1)
            p3, a3, dp3, _ = _lato_destro(u, st4)
            iterazioni.append((k, u, p2, p3))
            if p3 > 0 and abs(1.0 - p2 / p3) < tol:
                conv = True
                break
            du = -(p2 - p3) / (dp2 - dp3)
            un = u + du
            while not (lo < un < hi):
                du *= 0.5
                un = u + du
            u = un
        if not conv:
            f = lambda uu: _lato_sinistro(uu, st1)[0] - _lato_destro(uu, st4)[0]
            span = hi - lo
            u = brentq(f, lo + 1e-15 * span, hi - 1e-15 * span, xtol=1e-15, rtol=1e-15)
            metodo = "brentq"
        p = _lato_sinistro(u, st1)[0]
    else:
        raise GasdynError("modello: 'esatto' o 'isentropico'")
    onda_sx = onda_verso(st1, "u-a", u=u, modello=modello, nome="2")
    onda_dx = onda_verso(st4, "u+a", u=u, modello=modello, nome="3")
    st2, st3 = onda_sx.valle, onda_dx.valle
    sigla_c = "C"
    if g1 == g4 and abs(st2.rho / st3.rho - 1.0) < TOL_NULLA * 10:
        sigla_c = "N"
    conf = onda_sx.sigla + sigla_c + onda_dx.sigla
    return SoluzioneRiemann({"1": st1, "2": st2, "3": st3, "4": st4}, p, u, onda_sx, onda_dx, conf,
                            modello=modello, metodo=metodo, iterazioni=iterazioni, x0=x0, t0=t0)


def interazione(stato_sx: State, stato_dx: State, **kw) -> SoluzioneRiemann:
    """Interazione d'onde come pseudo-problema di Riemann fra i due stati esterni (F §9).

    Al primo incontro (urto–contatto, urto–urto, onda–contatto, ...) si
    trascura il supporto della zona d'interazione e si risolve il Riemann fra
    gli stati esterni noti. Per interazioni distribuite fra ventagli il
    risultato dà solo gli stati asintotici (F §9). Argomenti come :func:`riemann`.
    """
    return riemann(stato_sx, stato_dx, **kw)


def contorno_parete(stato: State, lato_parete: str, u_parete: float = 0.0, nome: str = ""):
    """Onda generata da una parete (ferma o in moto con ``u_parete``) sullo stato adiacente (F §9).

    ``lato_parete``: 'destra' (gas a sinistra, onda u−a) o 'sinistra' (onda u+a).
    Restituisce :class:`Urto` (compressione) o :class:`Ventaglio` (espansione).
    """
    fam = "u-a" if lato_parete.startswith("d") else "u+a"
    return onda_verso(stato, fam, u=u_parete, nome=nome)


def contorno_pressione(stato: State, lato: str, p_b: float, modello: str = "esatto", nome: str = ""):
    """Estremo aperto con pressione imposta p_b (F §9 "Pressione imposta").

    ``lato``: estremo del condotto ('destra' ⇒ onda riflessa u−a, 'sinistra' ⇒ u+a).
    Una compressione incidente su p costante riflette un'espansione e viceversa.
    """
    fam = "u-a" if lato.startswith("d") else "u+a"
    return onda_verso(stato, fam, p=p_b, modello=modello, nome=nome)


def contorno_pistone(stato: State, lato: str, u_p: float, nome: str = ""):
    """Pistone all'estremo ``lato`` ('sinistra' o 'destra') con velocità u_p (F §9)."""
    return contorno_parete(stato, lato, u_parete=u_p, nome=nome)


def contorno_serbatoio(stato: State, lato: str, p0: float, *, a0: float = None, T0: float = None,
                       modello: str = "esatto", nome: str = "") -> dict:
    """Estremo del condotto collegato a un grande serbatoio/ambiente (p0, T0, u = 0) (S p.44–45, F §9).

    ``lato``: estremo del condotto a cui sta il serbatoio ('sinistra' o 'destra').
    * Ingresso (gas che entra dal serbatoio): lo stato al contorno ha s e h0 del
      serbatoio, a_b² + δu_b² = a0², p_b = p0 (a_b/a0)^(2γ/(γ−1)); verso
      l'interno parte un'onda (famiglia u+a se il serbatoio è a sinistra) che
      porta ``stato`` a (p_b, u_b); fra i due gas c'è un contatto se le entropie
      differiscono. Con ``modello='esatto'`` l'onda interna è un urto se comprime
      (estensione di S p.45, che usa solo relazioni isentropiche);
      con ``modello='isentropico'`` è sempre un'onda semplice, come in S p.45.
    * Uscita (gas che entra nel serbatoio): p_b = p0 (F §9 "Pressione imposta").
    La scelta fra ingresso e uscita è automatica (segno di u al contorno).
    Restituisce dict: 'regime', 'stato' (gas al contorno), 'stato_interno'
    (gas interno dietro l'onda, stessi p e u), 'onda', 'contatto'.
    """
    gas = stato.gas
    if (a0 is None) == (T0 is None):
        raise GasdynError("specificare esattamente uno fra a0 e T0 del serbatoio")
    if T0 is not None:
        gas._serve_R("contorno_serbatoio con T0")
        a0 = math.sqrt(gas.gamma * gas.R * T0)
    g, d = gas.gamma, gas.delta
    sin = lato.startswith("s")
    sgn_in = 1.0 if sin else -1.0
    fam = "u+a" if sin else "u-a"
    s0 = p0 / (g * p0 / a0**2) ** g

    def p_ser(u):
        ab = math.sqrt(max(a0 * a0 - d * u * u, 0.0))
        return p0 * (ab / a0) ** (2.0 * g / (g - 1.0))

    def p_int(u):
        return onda_verso(stato, fam, u=u, modello=modello, consenti_vuoto=True).valle.p

    if p0 > p_int(0.0):  # il gas entra dal serbatoio
        umax = a0 / math.sqrt(d)
        f = lambda m: p_ser(sgn_in * m) - p_int(sgn_in * m)
        m = brentq(f, 0.0, umax * (1 - 1e-12), xtol=1e-14, rtol=1e-14)
        u_b = sgn_in * m
        ab = math.sqrt(a0 * a0 - d * u_b * u_b)
        stb = _stato_da_a_s(gas, ab, s0, u_b, nome)
        onda = onda_verso(stato, fam, u=u_b, modello=modello, nome=nome)
        nota = "ingresso subsonico" if m < ab else "ingresso supersonico: il serbatoio impone tutto lo stato"
        return {"regime": "ingresso", "stato": stb, "stato_interno": onda.valle, "onda": onda,
                "contatto": abs(stb.s / onda.valle.s - 1.0) > 1e-9, "nota": nota}
    onda = onda_verso(stato, fam, p=p0, modello=modello, nome=nome)
    return {"regime": "uscita", "stato": onda.valle, "stato_interno": onda.valle, "onda": onda,
            "contatto": False, "nota": "p = p0 al contorno"}

def riemann_area(sinistra: State, destra: State, A_sx: float, A_dx: float, tol: float = 1e-12) -> dict:
    """Interazione con una brusca variazione d'area in x0 (urto/flusso che attraversa la sezione).

    Modello del corso (F §9 "Variazione d'area", S p.49, B72-E02): il gas che
    attraversa la sezione segue il collegamento stazionario isentropico
    (:func:`variazione_area`); il resto è un problema di Riemann.
    Il flusso deve andare da sinistra a destra (u_sx > 0).

    * Flusso a monte supersonico (M_sx ≥ 1): nessuna onda può risalire;
      lo stato 5 = variazione_area(sinistra) e si risolve riemann(5, destra).
      Si verifica che l'onda sinistra di questo Riemann si muova verso destra
      (altrimenti l'ipotesi non regge e si solleva GasdynError).
    * Flusso subsonico: si cerca p' dello stato 2' dietro un'onda u−a riflessa
      a sinistra della sezione, tale che lo stato 5 (stazionario a destra)
      sia collegato a ``destra`` da contatto + onda u+a (brentq su p').
      **Estensione** del procedimento del corso (dichiarata).
      Vale anche per monte supersonico quando il caso precedente non regge
      (Strozzamento in A_dx, o onda sinistra che risalirebbe): l'onda riflessa
      è allora un urto u−a che risale la corrente (w ≤ 0) e lascia 2' subsonico.
    * Regime strozzato (``regime='strozzato'``), quando il caso subsonico non
      ha soluzione: la sezione minima diventa sonica. Il corso lo descrive solo
      a parole (S p.48 punto 4, F §16 "avviamento": gola sonica, a valle
      supersonico, l'informazione non risale); la costruzione seguente è
      **fuori formulario**:
        - A_dx < A_sx (o aree uguali): M5 = 1 in A_dx; 2' è subsonico con
          M2' dalla legge area–Mach fra A_dx (sonica) e A_sx; a valle
          riemann(5, destra) deve avere a sinistra un'espansione con testa
          ferma in x0 (u5 − a5 = 0).
        - A_dx > A_sx: M2' = 1 in A_sx; 5 sul ramo supersonico della legge
          area–Mach; riemann(5, destra) deve avere l'onda sinistra che si
          allontana dalla sezione (urto con w ≥ 0, come in S p.49, o espansione).
          Se invece quell'urto risalirebbe la sezione (w < 0), l'urto resta
          stazionario dentro il salto d'area, nella sezione A_s ∈ [A_sx, A_dx]:
          è la procedura "urto normale nel divergente" di F §16 (ramo
          supersonico fino ad A_s, urto normale, nuova A*, ramo subsonico fino
          ad A_dx) applicata a un divergente di lunghezza nulla; A_s si trova
          con brentq imponendo che 5 sia collegato a ``destra`` da contatto +
          onda u+a ('urto_stazionario' nel risultato).

    Restituisce dict con 'regime', gli stati ('stato_2p', 'stato_5'), le onde
    ('onda_riflessa' u−a, 'onda_trasmessa' u+a) e, quando a valle della sezione
    c'è anche un'onda u−a, 'riemann_valle' (Riemann fra 5 e ``destra`` nella
    sezione A_dx). Nel regime strozzato: 'sezione_sonica' ('A_sx' o 'A_dx') e
    'urto_stazionario' (None o dict con 'A_s', 'monte', 'valle').
    """
    if sinistra.u <= 0:
        raise GasdynError("riemann_area: richiede flusso verso destra a sinistra della sezione (u_sx > 0)")

    def v_min(o):
        return o.w if isinstance(o, Urto) else min(o.xi_testa, o.xi_coda)

    if sinistra.M >= 1.0:
        try:
            st5 = variazione_area(sinistra, A_sx, A_dx, "stesso", nome="5")
            sol = riemann(st5, destra)
            if v_min(sol.onda_sx) >= 0:
                return {"regime": "supersonico", "stato_5": st5, "riemann_valle": sol,
                        "onda_riflessa": None, "onda_trasmessa": sol.onda_dx}
        except Strozzamento:
            pass  # la portata non passa in A_dx: serve un urto riflesso che risale la corrente

    def onda_2p(pp):
        o = onda_verso(sinistra, "u-a", p=pp)  # onda riflessa che risale verso sinistra
        if isinstance(o, Urto) and o.w > 0:  # (per le espansioni basta 2' non supersonico)
            raise GasdynError("l'urto riflesso non risale la corrente")
        return o, o.valle

    # Limite di strozzamento: stato 2'c con la sezione minima sonica.
    # A_dx ≤ A_sx: M5 = 1 in A_dx, M2' dalla legge area–Mach; A_dx > A_sx: M2' = 1 in A_sx.
    gas = sinistra.gas
    if A_dx <= A_sx:
        sezione = "A_dx"
        M2c = 1.0 if A_dx == A_sx else mach_da_area_tra_sezioni(1.0, A_dx, A_sx, gas, "subsonico")
        M5c = 1.0
    else:
        sezione = "A_sx"
        M2c = 1.0
        M5c = mach_da_area_tra_sezioni(1.0, A_sx, A_dx, gas, "subsonico")
    lo, hi = sinistra.p * 1e-6, sinistra.p * 1e3
    rad_c = _radici(lambda pp: onda_2p(pp)[1].Mu - M2c, lo, hi, n=600, scala="log")
    p_c = rad_c[0] if len(rad_c) == 1 else None

    def stato5(pp):
        o, st2p = onda_2p(pp)
        if pp == p_c:  # bordo del dominio subsonico: 5 al limite di strozzamento
            return o, st2p, _stato_stazionario_M(st2p, M5c, nome="5")
        if not 0 < st2p.u < st2p.a:  # 2' deve essere subsonico ed entrante nella sezione
            raise GasdynError("stato 2' non subsonico")
        return o, st2p, variazione_area(st2p, A_sx, A_dx, "stesso", nome="5")

    def f(pp):
        _, _, st5 = stato5(pp)
        return st5.u - curva_onda(destra, "u+a", st5.p)

    rad = _radici(f, lo, hi, n=600, scala="log", extra=[p_c] if p_c is not None else ())
    if len(rad) > 1:
        raise GasdynError(f"riemann_area subsonico: {len(rad)} soluzioni trovate")
    if len(rad) == 1:
        o, st2p, st5 = stato5(rad[0])
        trasmessa = onda_verso(destra, "u+a", p=st5.p, nome="6")
        return {"regime": "subsonico", "onda_riflessa": o, "stato_2p": st2p, "stato_5": st5,
                "onda_trasmessa": trasmessa}

    # --- regime strozzato (fuori formulario) ---
    if p_c is None:
        raise NessunaSoluzione(f"riemann_area: nessuna soluzione subsonica e {len(rad_c)} stati 2' "
                               f"con M = {M2c:.6g} raggiungibili da un'onda riflessa che risale la corrente")
    o, st2p = onda_2p(p_c)
    st2p = st2p.con(nome="2'")
    risultato = {"regime": "strozzato", "sezione_sonica": sezione, "onda_riflessa": o,
                 "stato_2p": st2p, "urto_stazionario": None}
    if sezione == "A_dx":
        st5 = _stato_stazionario_M(st2p, 1.0, nome="5")
        sol = riemann(st5, destra)
        if not (isinstance(sol.onda_sx, Ventaglio) and not sol.onda_sx.compressione):
            raise NessunaSoluzione("riemann_area strozzato: a valle della sezione sonica l'onda u−a "
                                   "non è un'espansione (incoerente con il limite di strozzamento)")
        risultato.update(stato_5=st5, riemann_valle=sol, onda_trasmessa=sol.onda_dx)
        return risultato
    st5 = variazione_area(st2p, A_sx, A_dx, "supersonico", nome="5")
    sol = riemann(st5, destra)
    if v_min(sol.onda_sx) >= -1e-9 * st5.a:  # l'onda u−a si allontana dalla sezione
        risultato.update(stato_5=st5, riemann_valle=sol, onda_trasmessa=sol.onda_dx)
        return risultato

    # urto normale stazionario dentro il salto d'area (F §16 con divergente di lunghezza nulla)
    def urto_in(As):
        s1 = variazione_area(st2p, A_sx, As, "supersonico", nome="monte urto")
        s2 = urto_normale_stato(s1, nome="valle urto")
        return s1, s2, variazione_area(s2, As, A_dx, "subsonico", nome="5")

    def h(As):
        st5 = urto_in(As)[2]
        return st5.u - curva_onda(destra, "u+a", st5.p)

    h0, h1 = h(A_sx), h(A_dx)
    if h0 * h1 > 0:
        raise NessunaSoluzione("riemann_area strozzato: nessuna posizione dell'urto stazionario nel "
                               "salto d'area soddisfa la condizione a valle")
    As = brentq(h, A_sx, A_dx, xtol=1e-14 * A_dx, rtol=1e-14)
    s1, s2, st5 = urto_in(As)
    risultato.update(stato_5=st5, onda_trasmessa=onda_verso(destra, "u+a", p=st5.p, nome="6"),
                     urto_stazionario={"A_s": As, "monte": s1, "valle": s2})
    return risultato


# =============================================================================
# 7. Supersonico stazionario 2D (F §10–15, S p.50–61)
# =============================================================================
def _segno_C(famiglia: str) -> int:
    f = str(famiglia).replace("−", "-").replace(" ", "").upper()
    if f in ("C+", "+"):
        return +1
    if f in ("C-", "-"):
        return -1
    raise GasdynError(f"famiglia 2D '{famiglia}' non riconosciuta: usare 'C+' o 'C-'")


def angolo_mach(M: float) -> float:
    """Angolo di Mach μ = asin(1/M) [gradi] (F §10). Solleva NessunaSoluzione per M < 1."""
    if M < 1.0 - 1e-12:
        raise NessunaSoluzione(f"angolo di Mach non definito per M = {M} < 1")
    return math.degrees(math.asin(min(1.0, 1.0 / M)))


def prandtl_meyer(M, gas: Gas):
    """Funzione di Prandtl–Meyer ν(M) [gradi] (F §11), ν(1) = 0. Accetta array."""
    M = np.asarray(M, dtype=float)
    if np.any(M < 1.0 - 1e-12):
        raise NessunaSoluzione("ν(M) definita solo per M ≥ 1")
    g = gas.gamma
    k = math.sqrt((g + 1.0) / (g - 1.0))
    x = np.sqrt(np.maximum(M * M - 1.0, 0.0))
    nu = np.degrees(k * np.arctan(x / k) - np.arctan(x))
    return float(nu) if np.ndim(nu) == 0 else nu


def nu_max(gas: Gas) -> float:
    """ν_max = (π/2)(√((γ+1)/(γ−1)) − 1) [gradi] (F §11)."""
    return math.degrees(0.5 * math.pi * (math.sqrt((gas.gamma + 1.0) / (gas.gamma - 1.0)) - 1.0))


def mach_da_pm(nu: float, gas: Gas) -> float:
    """Inversa di Prandtl–Meyer: M da ν [gradi] (brentq su M ∈ [1, 1e8]).

    Solleva NessunaSoluzione per ν < 0 e :class:`Vuoto` per ν ≥ ν_max.
    """
    if nu < -1e-12:
        raise NessunaSoluzione(f"ν = {nu} < 0: la compressione porterebbe a M < 1")
    if nu <= 0:
        return 1.0
    if nu >= nu_max(gas):
        raise Vuoto(f"ν = {nu}° ≥ ν_max = {nu_max(gas):.4f}°: espansione nel vuoto")
    return brentq(lambda M: prandtl_meyer(M, gas) - nu, 1.0, 1e8, xtol=1e-14, rtol=1e-14)


def theta_da_sigma(M: float, sigma, gas: Gas):
    """Relazione θ–σ–M (F §12): tanθ = 2cotσ (M²sin²σ − 1)/(M²(γ + cos2σ) + 2). Angoli in gradi."""
    s = np.radians(np.asarray(sigma, dtype=float))
    num = 2.0 / np.tan(s) * (M * M * np.sin(s) ** 2 - 1.0)
    den = M * M * (gas.gamma + np.cos(2.0 * s)) + 2.0
    th = np.degrees(np.arctan(num / den))
    return float(th) if np.ndim(th) == 0 else th


def _M2_obliquo(M: float, sigma: float, gas: Gas) -> float:
    """Mach a valle M2 = M_n2/sin(σ − θ) (F §12)."""
    Mn1 = M * _sind(sigma)
    Mn2 = urto_normale(Mn1, gas)["M2"]
    return Mn2 / _sind(sigma - theta_da_sigma(M, sigma, gas))


def theta_max(M: float, gas: Gas) -> tuple:
    """Deviazione massima (θmax, σ a θmax) [gradi], da dθ/dσ = 0 (F §12), con minimize_scalar (xatol 1e-12°)."""
    mu = angolo_mach(M)
    r = minimize_scalar(lambda s: -theta_da_sigma(M, s, gas), bounds=(mu, 90.0), method="bounded",
                        options={"xatol": 1e-12})
    return -r.fun, r.x


def sonico(M: float, gas: Gas) -> tuple:
    """Punto sonico della polare (θ_s, σ_s): M2 = 1 (F §12), brentq su σ."""
    mu = angolo_mach(M)
    if M <= 1.0:
        return 0.0, 90.0
    s = brentq(lambda s: _M2_obliquo(M, s, gas) - 1.0, mu + 1e-12, 90.0, xtol=1e-13)
    return theta_da_sigma(M, s, gas), s


def sigma_da_theta(M: float, theta: float, gas: Gas, ramo: str = "debole") -> float:
    """Angolo d'urto σ [gradi] per deviazione |θ| [gradi] (F §12).

    ``ramo``: 'debole' (μ < σ < σ_max) o 'forte' (σ_max < σ ≤ 90°).
    Solleva :class:`UrtoStaccato` se |θ| > θmax(M).
    """
    th = abs(theta)
    tm, sm = theta_max(M, gas)
    if th > tm * (1 + 1e-12) + 1e-13:
        raise UrtoStaccato(f"θ = {th:.6g}° > θmax(M = {M:.6g}) = {tm:.6g}°: urto staccato")
    mu = angolo_mach(M)
    if th >= tm:
        return sm
    f = lambda s: theta_da_sigma(M, s, gas) - th
    if ramo == "debole":
        return mu if th == 0 else brentq(f, mu, sm, xtol=1e-13, rtol=1e-15)
    if ramo == "forte":
        return 90.0 if th == 0 else brentq(f, sm, 90.0, xtol=1e-13, rtol=1e-15)
    raise GasdynError("ramo: 'debole' o 'forte'")


@dataclass(frozen=True)
class Onda2D:
    """Onda stazionaria 2D: urto obliquo o ventaglio di Prandtl–Meyer (F §11–12).

    ``famiglia`` 'C+' (fronte θ+μ / θ+σ) o 'C-'. ``deflessione`` = θ2 − θ1
    [gradi, con segno]. Per gli urti: ``sigma`` (angolo rispetto alla velocità
    a monte) e ``ramo`` ('debole'/'forte'); ``angolo_fronte`` = θ1 ± σ
    (direzione assoluta della linea d'urto). Per i ventagli: ``angolo_testa``,
    ``angolo_coda`` = θ ± μ di monte e di valle.
    """

    tipo: str
    famiglia: str
    monte: State
    valle: State
    sigma: float = math.nan
    ramo: str = ""

    @property
    def deflessione(self):
        return self.valle.theta - self.monte.theta

    @property
    def _s(self):
        return _segno_C(self.famiglia)

    @property
    def angolo_fronte(self):
        return self.monte.theta + self._s * self.sigma

    @property
    def angolo_testa(self):
        return self.monte.theta + self._s * angolo_mach(self.monte.M)

    @property
    def angolo_coda(self):
        return self.valle.theta + self._s * angolo_mach(self.valle.M)

    @property
    def Mn1(self):
        return self.monte.M * _sind(self.sigma)

    @property
    def p2_p1(self):
        return self.valle.p / self.monte.p

    @property
    def p02_p01(self):
        return self.valle.p0 / self.monte.p0

    def __repr__(self):
        if self.tipo == "urto":
            info = f"σ={_fmt(self.sigma, 7)}° ({self.ramo}), Mn1={_fmt(self.Mn1)}, fronte {_fmt(self.angolo_fronte, 7)}°"
        else:
            info = f"testa {_fmt(self.angolo_testa, 7)}°, coda {_fmt(self.angolo_coda, 7)}°"
        return (f"Onda2D {self.tipo} {self.famiglia}: Δθ={_fmt(self.deflessione, 7)}°, {info}, "
                f"p2/p1={_fmt(self.p2_p1, 7)}, M {_fmt(self.monte.M)}→{_fmt(self.valle.M)}\n"
                f"  monte: {self.monte}\n  valle: {self.valle}")


def urto_obliquo(stato: State, *, theta: float = None, sigma: float = None, famiglia: str = "C+",
                 ramo: str = "debole", nome: str = "") -> Onda2D:
    """Urto obliquo su uno stato a direzione qualsiasi (F §12).

    Dati: ``theta`` = modulo della deviazione [gradi] (con ``ramo``) oppure
    ``sigma`` = angolo d'urto rispetto alla velocità a monte [gradi].
    Famiglia C+: θ2 = θ1 + θ; C−: θ2 = θ1 − θ (convenzione confermata in Fase 0).
    Stato a valle: componente tangenziale conservata, normale da Rankine–Hugoniot
    con M_n1 = M1 sinσ. Solleva :class:`UrtoStaccato` se θ > θmax.
    """
    s = _segno_C(famiglia)
    M1 = stato.M
    if M1 < 1.0:
        raise NessunaSoluzione(f"urto obliquo in corrente subsonica (M = {M1})")
    if (theta is None) == (sigma is None):
        raise GasdynError("urto_obliquo: specificare esattamente uno fra theta e sigma")
    if sigma is None:
        sigma = sigma_da_theta(M1, theta, stato.gas, ramo)
    mu = angolo_mach(M1)
    if not mu - 1e-9 <= sigma <= 90.0 + 1e-9:
        raise NessunaSoluzione(f"σ = {sigma}° fuori da [μ, 90°] = [{mu:.6g}, 90]")
    sigma = min(max(sigma, mu), 90.0)
    r = urto_normale(M1 * _sind(sigma), stato.gas)
    V1 = stato.V
    Vn1, Vt = V1 * _sind(sigma), V1 * _cosd(sigma)
    Vn2 = Vn1 / r["rho2_rho1"]
    dth = sigma - math.degrees(math.atan2(Vn2, Vt))
    th2 = stato.theta + s * dth
    V2 = math.hypot(Vn2, Vt)
    valle = State(stato.gas, stato.rho * r["rho2_rho1"], V2 * _cosd(th2), stato.p * r["p2_p1"],
                  V2 * _sind(th2), nome)
    sm = theta_max(M1, stato.gas)[1]
    return Onda2D("urto", "C+" if s > 0 else "C-", stato, valle, sigma,
                  "debole" if sigma <= sm else "forte")


def onda_isentropica_2d(stato: State, famiglia: str, *, theta2: float = None, M2: float = None,
                        p2: float = None, nome: str = "") -> Onda2D:
    """Onda semplice di Prandtl–Meyer (espansione o compressione isentropica) (F §11, §13).

    Attraverso un'onda C+ vale Δ(θ+ν) = 0, attraverso un'onda C− Δ(θ−ν) = 0:
    θ2 = θ1 ∓ (ν2 − ν1). Si assegna la direzione finale ``theta2`` [gradi]
    oppure ``M2`` o ``p2``. Grandezze totali conservate.
    """
    s = _segno_C(famiglia)
    gas = stato.gas
    M1 = stato.M
    nu1 = prandtl_meyer(M1, gas)
    if sum(x is not None for x in (theta2, M2, p2)) != 1:
        raise GasdynError("specificare esattamente uno fra theta2, M2, p2")
    if theta2 is not None:
        nu2 = nu1 - s * (theta2 - stato.theta)
        M2 = mach_da_pm(nu2, gas)
    else:
        if p2 is not None:
            M2 = mach_da_rapporto(stato.p0 / p2, gas, "p0/p")
        if M2 < 1:
            raise NessunaSoluzione("compressione isentropica oltre M = 1")
        nu2 = prandtl_meyer(M2, gas)
        theta2 = stato.theta - s * (nu2 - nu1)
    valle = State.da_totali(gas, M2, p0=stato.p0, a0=stato.a0, theta=theta2, nome=nome)
    tipo = "espansione" if nu2 >= nu1 else "compressione isentropica"
    return Onda2D(tipo, "C+" if s > 0 else "C-", stato, valle)


def espansione_pm(stato: State, dtheta: float, famiglia: str, nome: str = "") -> Onda2D:
    """Espansione di Prandtl–Meyer di modulo |Δθ| [gradi]: ν(M2) − ν(M1) = |Δθ| (F §11).

    C+: la corrente ruota in senso orario (θ2 = θ1 − |Δθ|); C−: antiorario.
    """
    s = _segno_C(famiglia)
    return onda_isentropica_2d(stato, famiglia, theta2=stato.theta - s * abs(dtheta), nome=nome)


def onda_verso_theta_2d(stato: State, theta2: float, famiglia: str, ramo: str = "debole",
                        modello: str = "esatto", nome: str = "") -> Onda2D:
    """Onda della famiglia data che porta la corrente alla direzione assoluta ``theta2`` [gradi].

    Condizione di parete (θ = θ_w, F §10). Compressione (s(θ2−θ1) > 0) ⇒ urto
    obliquo (ramo scelto) nel modello 'esatto', compressione isentropica nel
    modello 'isentropico'; espansione ⇒ Prandtl–Meyer.
    """
    s = _segno_C(famiglia)
    dth = theta2 - stato.theta
    if s * dth > 0 and modello == "esatto":
        return urto_obliquo(stato, theta=abs(dth), famiglia=famiglia, ramo=ramo, nome=nome)
    return onda_isentropica_2d(stato, famiglia, theta2=theta2, nome=nome)


def onda_verso_p_2d(stato: State, p2: float, famiglia: str, modello: str = "esatto",
                    nome: str = "") -> Onda2D:
    """Onda della famiglia data che porta la corrente alla pressione ``p2`` (F §13, §15).

    È la condizione di contorno a pressione assegnata (getto, slip line).
    p2 > p1 nel modello esatto: urto con M_n1² = 1 + (γ+1)/(2γ)(p2/p1 − 1)
    (σ unico; il ramo debole/forte segue). Altrimenti onda isentropica.
    Solleva NessunaSoluzione se p2 supera la pressione dell'urto normale.
    """
    if p2 > stato.p and modello == "esatto":
        g = stato.gamma
        Mn1 = math.sqrt(1.0 + (g + 1.0) / (2.0 * g) * (p2 / stato.p - 1.0))
        if Mn1 > stato.M * (1 + 1e-14):
            raise NessunaSoluzione(f"p2/p1 = {p2 / stato.p:.6g} oltre l'urto normale (M = {stato.M:.6g})")
        sigma = math.degrees(math.asin(min(1.0, Mn1 / stato.M)))
        return urto_obliquo(stato, sigma=sigma, famiglia=famiglia, nome=nome)
    return onda_isentropica_2d(stato, famiglia, p2=p2, nome=nome)


def _p_normale(stato: State) -> float:
    return stato.p * urto_normale(stato.M, stato.gas)["p2_p1"]


def _p_sonica_isentropica(stato: State) -> float:
    return stato.p0 * critiche(stato.gas)["p*_p0"]


def polare_p_theta(stato: State, famiglia: str = "entrambe", n: int = 400, espansione: bool = True,
                   p_min_rel: float = 0.05, modello: str = "esatto") -> dict:
    """Polare d'urto di Busemann nel piano p–θ con il ramo d'espansione (F §13).

    σ ∈ [μ1, 90°]: θ2 = θ1 ± θ(σ) (C+: +, C−: −), p2 = p1[1 + 2γ/(γ+1)(M1²sin²σ − 1)].
    Espansione parametrizzata da M2 ≥ M1 fino a p2 = ``p_min_rel``·p1.
    Restituisce per ogni famiglia un dict con 'urto' (theta, p, sigma, M2,
    ramo), 'espansione' (theta, p, M) e 'punti' notevoli (θmax, sonico,
    normale) come (θ, p). ``modello='isentropico'``: ramo di compressione
    isentropica (fino a M = 1) al posto dell'urto.
    """
    fams = ["C+", "C-"] if famiglia == "entrambe" else [famiglia]
    out = {}
    M1, g = stato.M, stato.gamma
    for fam in fams:
        s = _segno_C(fam)
        d = {}
        if modello == "esatto":
            mu = angolo_mach(M1)
            sig = np.linspace(mu, 90.0, n)
            th = theta_da_sigma(M1, sig, stato.gas)
            p = stato.p * (1.0 + 2.0 * g / (g + 1.0) * (M1**2 * np.sin(np.radians(sig)) ** 2 - 1.0))
            M2 = np.array([_M2_obliquo(M1, x, stato.gas) for x in sig])
            tm, sm = theta_max(M1, stato.gas)
            d["urto"] = {"theta": stato.theta + s * th, "p": p, "sigma": sig, "M2": M2,
                         "ramo": np.where(sig <= sm, "debole", "forte")}
            ts, ss = sonico(M1, stato.gas)
            pf = lambda sg: stato.p * (1.0 + 2.0 * g / (g + 1.0) * (M1**2 * _sind(sg) ** 2 - 1.0))
            d["punti"] = {"theta_max": (stato.theta + s * tm, pf(sm)), "sonico": (stato.theta + s * ts, pf(ss)),
                          "normale": (stato.theta, _p_normale(stato)), "origine": (stato.theta, stato.p)}
        else:
            Ms = np.linspace(M1, 1.0, n)
            nu1 = prandtl_meyer(M1, stato.gas)
            nus = prandtl_meyer(Ms, stato.gas)
            d["compressione"] = {"theta": stato.theta - s * (nus - nu1),
                                 "p": stato.p0 / isentropiche(Ms, stato.gas)["p0_p"], "M": Ms}
            d["punti"] = {"origine": (stato.theta, stato.p)}
        if espansione:
            p_lo = stato.p * p_min_rel
            M_hi = mach_da_rapporto(stato.p0 / p_lo, stato.gas, "p0/p")
            Ms = np.linspace(M1, M_hi, n)
            nu1 = prandtl_meyer(M1, stato.gas)
            d["espansione"] = {"theta": stato.theta - s * (prandtl_meyer(Ms, stato.gas) - nu1),
                               "p": stato.p0 / isentropiche(Ms, stato.gas)["p0_p"], "M": Ms}
        out[fam] = d
    return out if famiglia == "entrambe" else out[famiglia]


def polare_odografo(stato: State, famiglia: str = "entrambe", n: int = 400, espansione: bool = True,
                    riferimento: str = "assoluto", p_min_rel: float = 0.05) -> dict:
    """Polare d'urto nel piano odografo (u*, v*) = (V/a*)(cosθ, sinθ) (F §13, S p.56).

    ``riferimento='assoluto'``: componenti nel riferimento x–y (come chiesto
    nelle tracce, es. B51-E01); ``'monte'``: assi ruotati sulla velocità a
    monte (piano (V_n/a*, V_t/a*) di F §13 e S p.57, stato a monte sull'asse).
    a* è quello della corrente (si conserva attraverso urti ed espansioni).
    Restituisce per famiglia i rami 'urto' ed 'espansione' con chiavi
    'us', 'vs', più 'stato' = (u*, v*) del monte e i raggi dei cerchi
    'cerchio_sonico' = 1 e 'cerchio_limite' = √((γ+1)/(γ−1)).
    """
    pol = polare_p_theta(stato, famiglia, n, espansione, p_min_rel)
    if famiglia != "entrambe":
        pol = {famiglia: pol}
    ast = stato.a_star
    rot = stato.theta if riferimento == "monte" else 0.0
    out = {}
    for fam, d in pol.items():
        o = {}
        if "urto" in d:
            u = d["urto"]
            V2 = []
            for sg in u["sigma"]:
                Vn1, Vt = stato.V * _sind(sg), stato.V * _cosd(sg)
                Vn2 = Vn1 / urto_normale(stato.M * _sind(sg), stato.gas)["rho2_rho1"]
                V2.append(math.hypot(Vn2, Vt))
            V2 = np.array(V2) / ast
            th = np.radians(u["theta"] - rot)
            o["urto"] = {"us": V2 * np.cos(th), "vs": V2 * np.sin(th), "ramo": u["ramo"]}
        if "espansione" in d:
            e = d["espansione"]
            Ms = e["M"]
            Mst = np.sqrt((stato.gamma + 1.0) * Ms**2 / (2.0 + (stato.gamma - 1.0) * Ms**2))
            th = np.radians(e["theta"] - rot)
            o["espansione"] = {"us": Mst * np.cos(th), "vs": Mst * np.sin(th)}
        out[fam] = o
    th0 = math.radians(stato.theta - rot)
    res = out if famiglia == "entrambe" else out[famiglia]
    meta = {"stato": (stato.M_star * math.cos(th0), stato.M_star * math.sin(th0)), "cerchio_sonico": 1.0,
            "cerchio_limite": math.sqrt((stato.gamma + 1.0) / (stato.gamma - 1.0))}
    res = dict(res)
    res.update(meta)
    return res


def _theta_di_p(stato: State, fam: str, p: float, modello: str) -> float:
    return onda_verso_p_2d(stato, p, fam, modello).valle.theta


def intersezioni_polari(A: State, famA: str, B: State, famB: str, modello: str = "esatto",
                        n: int = 800) -> list:
    """Tutte le intersezioni delle polari p–θ di A (famiglia famA) e B (famB) (F §13).

    Su ogni polare di una famiglia θ è funzione a un sol valore di p (espansione
    per p < p_i, urto per p_i < p ≤ p_normale), quindi le intersezioni sono le
    radici di θ_A(p) − θ_B(p), cercate per cambio di segno su n punti in scala
    logaritmica e raffinate con brentq. ``famA``/``famB`` possono essere
    'entrambe'. Restituisce una lista di dict con p, theta, onda_A, onda_B
    (oggetti :class:`Onda2D`, con ramo dell'urto) ordinata per p.
    """
    fA = ["C+", "C-"] if famA == "entrambe" else [famA]
    fB = ["C+", "C-"] if famB == "entrambe" else [famB]
    lim = _p_normale if modello == "esatto" else _p_sonica_isentropica
    p_lim = min(lim(A), lim(B))
    p_hi = p_lim * (1 - 1e-12)
    p_lo = 1e-6 * min(A.p, B.p)
    out = []
    for fa in fA:
        for fb in fB:
            f = lambda p: _theta_di_p(A, fa, p, modello) - _theta_di_p(B, fb, p, modello)
            radici = _radici(f, p_lo, p_hi, n=n, scala="log")
            # estremo superiore (urto normale di una delle due polari): intersezione nel punto
            # d'urto normale (es. criterio di von Neumann), tolleranza 1e-6°
            if abs(f(p_lim)) < 1e-6 and not any(abs(r / p_lim - 1) < 1e-6 for r in radici):
                radici.append(p_lim)
            for p in radici:
                oA = onda_verso_p_2d(A, p, fa, modello)
                oB = onda_verso_p_2d(B, p, fb, modello)
                out.append({"p": p, "theta": oA.valle.theta, "fam_A": fa, "fam_B": fb,
                            "onda_A": oA, "onda_B": oB})
    out.sort(key=lambda r: r["p"])
    return out


@dataclass
class SoluzioneRiemann2D:
    """Incontro di due correnti supersoniche separate da una slip line (F §13, B51-E01).

    ``onda_sopra`` (famiglia C+, nella corrente superiore) e ``onda_sotto`` (C−).
    A valle: p e θ comuni (``p``, ``theta``), slip line inclinata di θ.
    """

    sopra: State
    sotto: State
    onda_sopra: Onda2D
    onda_sotto: Onda2D
    p: float
    theta: float
    modello: str

    def __repr__(self):
        return (f"Riemann stazionario ({self.modello}): p = {_fmt(self.p, 7)}, θ slip = {_fmt(self.theta, 7)}°\n"
                f"  sopra: {self.onda_sopra}\n  sotto: {self.onda_sotto}")


def riemann_stazionario(sopra: State, sotto: State, modello: str = "esatto",
                        ramo: str = "debole") -> SoluzioneRiemann2D:
    """Problema di Riemann stazionario 2D: due correnti supersoniche che si incontrano (F §13).

    La corrente superiore emette un'onda C+, l'inferiore un'onda C−; si impone
    p e θ uguali (slip line). Sui rami deboli (ed espansioni) θ_sopra(p) cresce
    e θ_sotto(p) decresce, quindi la soluzione debole è unica; altre
    intersezioni possono esistere sui rami forti. ``ramo='debole'`` (default,
    come in campo aperto, F §14) scarta le intersezioni con urti forti;
    ``ramo='tutti'`` solleva :class:`RamiMultipli` se ce n'è più d'una (tutte in
    ``.soluzioni``; vedi anche :func:`intersezioni_polari`).
    ``modello='isentropico'``: compressioni isentropiche invece degli urti
    (confronto richiesto in B51-E01).
    """
    r = intersezioni_polari(sopra, "C+", sotto, "C-", modello)
    if ramo == "debole":
        r = [x for x in r if x["onda_A"].ramo != "forte" and x["onda_B"].ramo != "forte"]
    if not r:
        raise NessunaSoluzione("le polari non si intersecano: nessuna configurazione regolare "
                               "(possibile disco di Mach / urto staccato)")
    if len(r) > 1:
        raise RamiMultipli(f"{len(r)} intersezioni", r)
    x = r[0]
    return SoluzioneRiemann2D(sopra, sotto, x["onda_A"], x["onda_B"], x["p"], x["theta"], modello)


def riflessione_regolare(stato1: State, theta_w: float, ramo_riflesso: str = "debole") -> dict:
    """Riflessione regolare di un urto incidente su una parete orizzontale (F §14).

    Corrente 1 parallela alla parete (θ1 = 0), parete sotto; urto incidente C−
    (ramo debole) che devia di −θ_w; urto riflesso C+ che riporta θ3 = 0.
    Restituisce gli stati 1, 2, 3, i due urti e l'angolo fra urto riflesso e
    parete. Solleva :class:`UrtoStaccato` se θ_w > θmax(M2) (niente RR: MR).
    """
    st1 = stato1.con(nome="1")
    inc = urto_obliquo(st1, theta=theta_w, famiglia="C-", ramo="debole", nome="2")
    try:
        rif = onda_verso_theta_2d(inc.valle, st1.theta, "C+", ramo=ramo_riflesso, nome="3")
    except UrtoStaccato as e:
        raise UrtoStaccato(f"riflessione regolare impossibile (M2 = {inc.valle.M:.5g}): {e}") from None
    return {"1": st1, "2": inc.valle, "3": rif.valle, "incidente": inc, "riflesso": rif,
            "angolo_riflesso_parete": rif.angolo_fronte - st1.theta}


def riflessione_mach(stato1: State, theta_w: float) -> dict:
    """Teoria dei tre urti al punto triplo (F §14).

    Stato 2 dietro l'urto incidente (C−, debole); stato 3 sulla polare riflessa
    di 2 (C+); stato 4 dietro lo stelo, sul ramo forte della polare di 1;
    slip line: p3 = p4, θ3 = θ4. Si cercano tutte le intersezioni della polare
    riflessa con la polare di 1 (entrambe le famiglie) escludendo l'origine
    della polare riflessa; restano quelle sul ramo forte della polare di 1.
    'soluzioni' contiene tutte quelle trovate, 'principale' quella con |θ| minima.
    """
    st1 = stato1.con(nome="1")
    inc = urto_obliquo(st1, theta=theta_w, famiglia="C-", ramo="debole", nome="2")
    st2 = inc.valle
    sols = []
    for r in intersezioni_polari(st2, "C+", st1, "entrambe"):
        if r["onda_A"].tipo != "urto" or abs(r["p"] / st2.p - 1) < 1e-8:
            continue
        stem = r["onda_B"]
        if stem.tipo == "urto" and stem.ramo == "forte":
            sols.append({"p": r["p"], "theta": r["theta"], "riflesso": r["onda_A"], "stelo": stem,
                         "3": r["onda_A"].valle.con(nome="3"), "4": stem.valle.con(nome="4")})
    if not sols:
        raise NessunaSoluzione("nessuna soluzione a tre urti con stelo sul ramo forte")
    princ = min(sols, key=lambda s: abs(s["theta"]))
    return {"1": st1, "2": st2, "incidente": inc, "soluzioni": sols, "principale": princ}


def theta_von_neumann(M1: float, gas: Gas) -> float:
    """Angolo θ_N del criterio di von Neumann [gradi] (F §14).

    θ_w per cui la riflessione regolare (riflesso debole) dà p3 uguale alla
    pressione dell'urto normale di 1 (le polari si incontrano sull'asse θ = 0).
    Solleva NessunaSoluzione se non esiste (M1 bassi, polare riflessa che non
    raggiunge l'asse). Risolto con brentq su θ_w ∈ (0, θ_D).
    """
    st1 = State.da_M(gas, M1, p=1.0, a=1.0)
    pn = _p_normale(st1)
    tD = theta_distacco(M1, gas)
    f = lambda tw: riflessione_regolare(st1, tw)["3"].p - pn
    a, b = 1e-6, tD * (1 - 1e-9)
    if f(a) * f(b) > 0:
        raise NessunaSoluzione(f"criterio di von Neumann assente per M1 = {M1} (θ_N > θ_D)")
    return brentq(f, a, b, xtol=1e-12)


def theta_distacco(M1: float, gas: Gas) -> float:
    """Angolo θ_D del criterio di distacco [gradi] (F §14): la polare riflessa è
    tangente alla retta di parete, θmax(M2(θ_w)) = θ_w. brentq su θ_w."""
    st1 = State.da_M(gas, M1, p=1.0, a=1.0)
    tm1 = theta_max(M1, gas)[0]

    def f(tw):
        M2 = urto_obliquo(st1, theta=tw, famiglia="C-").valle.M
        if M2 <= 1.0:
            return -tw
        return theta_max(M2, gas)[0] - tw

    xs = np.linspace(1e-6, tm1 * (1 - 1e-9), 400)
    fs = [f(x) for x in xs]
    for i in range(len(xs) - 1):
        if fs[i] > 0 >= fs[i + 1]:
            return brentq(f, xs[i], xs[i + 1], xtol=1e-12)
    raise NessunaSoluzione(f"criterio di distacco non trovato per M1 = {M1}")


def diagramma_rr_mr(M_valori, gas: Gas) -> dict:
    """Curve θ_N(M1) e θ_D(M1) del diagramma θ–M1 (F §14); NaN dove il criterio non esiste."""
    tn, td = [], []
    for M in M_valori:
        try:
            td.append(theta_distacco(M, gas))
        except NessunaSoluzione:
            td.append(math.nan)
        try:
            tn.append(theta_von_neumann(M, gas))
        except (NessunaSoluzione, GasdynError):
            tn.append(math.nan)
    return {"M": np.asarray(M_valori, float), "theta_N": np.array(tn), "theta_D": np.array(td)}


def shock_expansion(superiore: Sequence, inferiore: Sequence, M: float, alpha: float, gas: Gas,
                    p_inf: float = 1.0) -> dict:
    """Teoria urto–espansione per lastre e profili a spigoli (F §11–12; teoria standard).

    ``superiore``, ``inferiore``: liste di punti (x, y) del profilo nel suo
    riferimento, dal bordo d'attacco al bordo d'uscita (stessi estremi).
    ``alpha``: incidenza [gradi] (positiva a cabrare: nel riferimento del corpo
    la corrente ha direzione θ∞ = +α). Dorso: onde C+; ventre: onde C−.
    Ogni spigolo: urto obliquo se la parete comprime, Prandtl–Meyer se espande.
    Coefficienti riferiti a q∞ = ½γp∞M² e alla corda c = |TE − LE|;
    portanza ⟂ e resistenza ∥ alla corrente asintotica.
    """
    a_inf = 1.0
    st_inf = State.da_M(gas, M, alpha, p=p_inf, a=a_inf, nome="∞")
    q = 0.5 * gas.gamma * p_inf * M * M
    pannelli = []
    F = np.zeros(2)
    for nome, pts, fam, sgn_n in (("dorso", superiore, "C+", 1.0), ("ventre", inferiore, "C-", -1.0)):
        st = st_inf
        for i in range(len(pts) - 1):
            (x0, y0), (x1, y1) = pts[i], pts[i + 1]
            L = math.hypot(x1 - x0, y1 - y0)
            phi = math.degrees(math.atan2(y1 - y0, x1 - x0))
            onda = onda_verso_theta_2d(st, phi, fam)
            st = onda.valle
            n = sgn_n * np.array([-(y1 - y0), x1 - x0]) / L  # normale uscente
            F += -(st.p - p_inf) * n * L
            pannelli.append({"superficie": nome, "i": i, "phi": phi, "onda": onda.tipo, "M": st.M,
                             "p_pinf": st.p / p_inf, "Cp": (st.p - p_inf) / q})
    c = math.hypot(superiore[-1][0] - superiore[0][0], superiore[-1][1] - superiore[0][1])
    ca, sa = _cosd(alpha), _sind(alpha)
    D = F[0] * ca + F[1] * sa
    Lf = -F[0] * sa + F[1] * ca
    return {"pannelli": pannelli, "CL": Lf / (q * c), "CD": D / (q * c), "CN": F[1] / (q * c),
            "CA": F[0] / (q * c), "F": F}


def pitot_supersonico(M1: float, gas: Gas) -> float:
    """Rapporto p02/p1 letto da un Pitot (formula di Rayleigh, **fuori formulario**).

    M1 > 1: urto normale davanti alla presa e compressione isentropica,
    p02/p1 = [(γ+1)M1²/2]^(γ/(γ−1)) / [(2γM1² − (γ−1))/(γ+1)]^(1/(γ−1)).
    M1 ≤ 1: p0/p isentropico.
    """
    g = gas.gamma
    if M1 <= 1.0:
        return isentropiche(M1, gas)["p0_p"]
    return (0.5 * (g + 1.0) * M1**2) ** (g / (g - 1.0)) / (
        (2.0 * g * M1**2 - (g - 1.0)) / (g + 1.0)) ** (1.0 / (g - 1.0))


def mach_da_pitot(p02_p1: float, gas: Gas) -> float:
    """Inversa di :func:`pitot_supersonico` (ramo subsonico in forma chiusa, supersonico con brentq)."""
    soglia = isentropiche(1.0, gas)["p0_p"]
    if p02_p1 <= soglia:
        return mach_da_rapporto(p02_p1, gas, "p0/p")
    return brentq(lambda M: pitot_supersonico(M, gas) - p02_p1, 1.0, _MACH_SUP_MAX, xtol=1e-14)


def punto_caratteristiche_2d(A: State, B: State, posA: tuple = None, posB: tuple = None,
                             nome: str = "P") -> dict:
    """Punto P all'incontro della C− uscente da A e della C+ uscente da B (F §10).

    θ_P = [(θ_A+ν_A) + (θ_B−ν_B)]/2, ν_P = [(θ_A+ν_A) − (θ_B−ν_B)]/2; M_P da ν_P,
    stato dai totali di A (flusso omoentropico e omoentalpico). Se si danno le
    posizioni, P è l'intersezione delle rette con pendenze medie
    tan((θ−μ)_A,P) e tan((θ+μ)_B,P) (metodo delle caratteristiche standard).
    """
    gas = A.gas
    Kp = A.theta + prandtl_meyer(A.M, gas)
    Km = B.theta - prandtl_meyer(B.M, gas)
    thP, nuP = 0.5 * (Kp + Km), 0.5 * (Kp - Km)
    MP = mach_da_pm(nuP, gas)
    P = State.da_totali(gas, MP, p0=A.p0, a0=A.a0, theta=thP, nome=nome)
    out = {"theta": thP, "nu": nuP, "M": MP, "mu": angolo_mach(MP), "stato": P}
    if posA is not None and posB is not None:
        mA = _tand(0.5 * ((A.theta - angolo_mach(A.M)) + (thP - angolo_mach(MP))))
        mB = _tand(0.5 * ((B.theta + angolo_mach(B.M)) + (thP + angolo_mach(MP))))
        (xa, ya), (xb, yb) = posA, posB
        x = (yb - ya + mA * xa - mB * xb) / (mA - mB)
        out["x"], out["y"] = x, ya + mA * (x - xa)
    return out


def ugello_minima_lunghezza(Me: float, gas: Gas, n: int = 60, h_gola: float = 1.0) -> dict:
    """Ugello piano di lunghezza minima con il metodo delle caratteristiche (F §10; schema standard).

    Semi-ugello sopra l'asse di simmetria (y = 0), semialtezza di gola
    ``h_gola``, espansione centrata nello spigolo di gola con n caratteristiche
    C− di direzione θ_j ∈ (0, θ_w,max], θ_w,max = ν(Me)/2. Nel ventaglio θ = ν
    (K+ = θ+ν = 2θ_j lungo la C− j); la C+ k nasce sull'asse con K− = −2θ_k,
    quindi nel punto (k, j): θ = θ_j − θ_k, ν = θ_j + θ_k. Rette fra i nodi
    con pendenze medie; parete non riflettente: ogni tratto di parete ha la
    pendenza media fra le direzioni del flusso ai suoi estremi.
    Restituisce θ_w,max, il profilo (x, y), L/h_gola e h_uscita/h_gola (→ A/A*(Me)
    per n → ∞). **Discretizzazione**: L converge come O(1/n).
    """
    twmax = 0.5 * prandtl_meyer(Me, gas)
    th = np.linspace(twmax / n, twmax, n)
    mu = lambda nu: angolo_mach(mach_da_pm(nu, gas)) if nu > 0 else 90.0
    xc, yc = 0.0, h_gola
    P = {}
    wx, wy, wth = [xc], [yc], [twmax]
    for k in range(n):
        for j in range(k, n):
            t_, nu_ = th[j] - th[k], th[j] + th[k]
            if k == 0:
                xa, ya, ta, na = xc, yc, th[j], th[j]
            else:
                xa, ya, ta, na = P[(k - 1, j)]
            mA = _tand(0.5 * ((ta - mu(na)) + (t_ - mu(nu_))))
            if j == k:  # punto sull'asse
                x = xa - ya / mA
                P[(k, j)] = (x, 0.0, t_, nu_)
            else:
                xb, yb, tb, nb = P[(k, j - 1)]
                mB = _tand(0.5 * ((tb + mu(nb)) + (t_ + mu(nu_))))
                x = (yb - ya + mA * xa - mB * xb) / (mA - mB)
                P[(k, j)] = (x, ya + mA * (x - xa), t_, nu_)
        xb, yb, tb, nb = P[(k, n - 1)]
        mC = _tand(tb + mu(nb))
        mW = _tand(0.5 * (wth[-1] + tb))
        x = (yb - wy[-1] + mW * wx[-1] - mC * xb) / (mW - mC)
        y = wy[-1] + mW * (x - wx[-1])
        wx.append(x)
        wy.append(y)
        wth.append(tb)
    return {"theta_w_max": twmax, "x_parete": np.array(wx), "y_parete": np.array(wy),
            "L_su_h_gola": wx[-1] / h_gola, "h_uscita_su_h_gola": wy[-1] / h_gola,
            "A_su_Astar_teorico": area_su_area_critica(Me, gas), "nodi": P}

def getto_adattamento(Me: float, pe: float, pa: float, gas: Gas) -> dict:
    """Primo fenomeno al labbro di un getto piano che esce a (Me, pe) in ambiente pa (F §15).

    Sottoespanso (pe > pa): espansione di Prandtl–Meyer fino a pa,
    M_a = {[(p0/pa)^((γ−1)/γ) − 1]/δ}^½, |Δθ| = ν(M_a) − ν(Me); apertura
    iniziale totale 2|Δθ| per un getto piano simmetrico.
    Sovraespanso (pe < pa): urto obliquo al labbro con p2 = pa (se compatibile).
    """
    st = State.da_M(gas, Me, p=pe, a=1.0)
    if abs(pe - pa) <= 1e-12 * pa:
        return {"regime": "adattato", "Ma": Me, "deflessione": 0.0}
    if pe > pa:
        Ma = mach_da_rapporto(st.p0 / pa, gas, "p0/p")
        dth = prandtl_meyer(Ma, gas) - prandtl_meyer(Me, gas)
        return {"regime": "sottoespanso", "Ma": Ma, "deflessione": dth, "apertura_totale": 2 * dth}
    o = onda_verso_p_2d(st, pa, "C-")
    return {"regime": "sovraespanso", "Ma": o.valle.M, "deflessione": abs(o.deflessione),
            "sigma": o.sigma, "onda": o}


def sauer_meyer(R: float, gas: Gas, versione: str = "formulario") -> dict:
    """Soluzione di Sauer–Meyer nella gola di raggio adimensionale R ≫ 1 (F §19).

    ``versione='formulario'`` (default, F §19 e S p.65 "esercizio 2024"):
    d = 1/√((γ+1)R); u = dξ + (γ+1)d²η²/2; v = (γ+1)d²ξη + (γ+1)d³η³/6;
    linea di gola ξ_g = −dη²/6; linea sonica ξ_s = −(γ+1)dη²/2.
    ``versione='favini'`` (S p.65, appunti senza (γ+1)): d = 1/√R,
    u = dξ + d²η²/2, v = d²ξη + d³η³/6, ξ_g = −dη²/6, ξ_s = −dη²/2.
    È la discrepanza 10 del report di Fase 0; lo svolgimento B67-E01 (2022)
    torna con la versione 'favini' (vedi verifica_catalogo.md).
    In entrambe M ≈ √((1+u)² + v²); lunghezze scalate con la semialtezza di gola.
    Restituisce d e funzioni u(ξ,η), v(ξ,η), M(ξ,η), xi_gola(η), xi_sonica(η).
    """
    g = gas.gamma
    if versione == "formulario":
        k = g + 1.0
    elif versione == "favini":
        k = 1.0
    else:
        raise GasdynError("versione: 'formulario' o 'favini'")
    d = 1.0 / math.sqrt(k * R)
    u = lambda xi, eta: d * xi + 0.5 * k * d * d * eta**2
    v = lambda xi, eta: k * d * d * xi * eta + k * d**3 * eta**3 / 6.0
    return {"d": d, "u": u, "v": v, "M": lambda xi, eta: np.sqrt((1 + u(xi, eta)) ** 2 + v(xi, eta) ** 2),
            "xi_gola": lambda eta: -d * eta**2 / 6.0, "xi_sonica": lambda eta: -0.5 * k * d * eta**2}


# =============================================================================
# 8. Grafici (matplotlib, importato solo quando serve)
# =============================================================================
def _plt():
    import matplotlib.pyplot as plt
    return plt


_STILE = {"urto": dict(color="k", lw=2.0), "contatto": dict(color="k", lw=1.2, ls="--"),
          "ventaglio": dict(color="0.35", lw=0.8), "pistone": dict(color="tab:blue", lw=2.0),
          "parete": dict(color="0.2", lw=3.0)}


def grafico_xt(soluzioni, t_max: float, x_range: tuple = None, ax=None, n_ventaglio: int = 7,
               etichette: bool = True, linee: Sequence = (), pistone: Callable = None,
               titolo: str = ""):
    """Diagramma x–t di una o più soluzioni di Riemann (F §8: rette x = x0 + λt).

    ``soluzioni``: una :class:`SoluzioneRiemann` o una lista di elementi
    ``SoluzioneRiemann`` oppure ``(SoluzioneRiemann, t_fine)`` (ognuna parte
    dai propri x0, t0 e si disegna fino a t_fine o t_max). Urti: linea spessa;
    contatto: tratteggio; ventagli: ``n_ventaglio`` caratteristiche.
    ``linee``: segmenti aggiuntivi ((x0, t0), (x1, t1), stile, etichetta) con
    stile in {'urto', 'contatto', 'ventaglio', 'pistone', 'parete'} o dict.
    ``pistone``: funzione x_p(t) da tracciare. Etichette degli stati nelle regioni.
    """
    plt = _plt()
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4.5))
    if isinstance(soluzioni, SoluzioneRiemann):
        soluzioni = [soluzioni]
    xs_all = []
    for el in soluzioni:
        sol, tf = (el if isinstance(el, tuple) else (el, t_max))
        x0, t0 = sol.x0, sol.t0
        dt = tf - t0
        if dt <= 0:
            continue
        v = sol.velocita()
        for nome, lam in v.items():
            stile = "urto" if nome.startswith("urto") else "contatto" if nome == "contatto" else "ventaglio"
            ax.plot([x0, x0 + lam * dt], [t0, tf], **_STILE[stile])
            xs_all += [x0, x0 + lam * dt]
        for lato in ("sx", "dx"):
            if f"testa_{lato}" in v:
                for lam in np.linspace(v[f"testa_{lato}"], v[f"coda_{lato}"], n_ventaglio)[1:-1]:
                    ax.plot([x0, x0 + lam * dt], [t0, tf], **_STILE["ventaglio"])
        if etichette:
            tl = t0 + 0.8 * dt
            confini = sorted(v.values())
            xi_reg = ([confini[0] - 0.3 * max(1.0, abs(confini[0]))]
                      + [0.5 * (a + b) for a, b in zip(confini[:-1], confini[1:])]
                      + [confini[-1] + 0.3 * max(1.0, abs(confini[-1]))])
            nomi = []
            for xi in xi_reg:
                st = sol.stato_xi(xi)
                nomi.append(st.nome)
            for xi, nm in zip(xi_reg, nomi):
                if nm and not nm.startswith("vuoto"):
                    xl = x0 + xi * (tl - t0)
                    if x_range is not None:  # etichetta dentro la finestra
                        m = 0.04 * (x_range[1] - x_range[0])
                        xl = min(max(xl, x_range[0] + m), x_range[1] - m)
                    ax.text(xl, tl, nm, ha="center", va="center", fontsize=10,
                            bbox=dict(fc="w", ec="none", alpha=0.8, pad=1))
    for seg in linee:
        (xa, ta), (xb, tb), st = seg[0], seg[1], seg[2]
        st = _STILE.get(st, st) if isinstance(st, str) else st
        ax.plot([xa, xb], [ta, tb], label=seg[3] if len(seg) > 3 else None, **st)
        xs_all += [xa, xb]
    if pistone is not None:
        tt = np.linspace(0, t_max, 200)
        ax.plot([pistone(t) for t in tt], tt, **_STILE["pistone"], label="pistone")
    if x_range is None and xs_all:
        x_range = (min(xs_all), max(xs_all))
    if x_range:
        ax.set_xlim(*x_range)
    ax.set_ylim(0, t_max)
    ax.set_xlabel("x")
    ax.set_ylabel("t")
    if titolo:
        ax.set_title(titolo)
    ax.grid(alpha=0.25)
    return ax


def grafico_pu(sol: SoluzioneRiemann = None, stati: Sequence = (), curve: Sequence = (), ax=None,
               p_max: float = None, titolo: str = ""):
    """Piano p–u (u in ascissa, p in ordinata) con curve d'onda e stati etichettati (F §8).

    ``sol``: disegna la curva u−a dallo stato 1 e la u+a dallo stato 4 e gli
    stati 1, 2(=3), 4. ``curve``: coppie (stato, famiglia) aggiuntive.
    ``stati``: stati da segnare (etichetta = ``stato.nome``). Ramo d'urto
    continuo, ramo d'espansione tratteggiato.
    """
    plt = _plt()
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4.5))
    cc, ss = list(curve), list(stati)
    if sol is not None:
        cc = [(sol.stati["1"], "u-a"), (sol.stati["4"], "u+a")] + cc
        ss = [sol.stati["1"], sol.stati["4"]] + ([] if sol.vuoto else [sol.stati["2"].con(nome="2=3")]) + ss
    pm = p_max or 1.6 * max([s.p for s in ss] + [c[0].p for c in cc])
    for st, fam in cc:
        pe = np.linspace(1e-6 * st.p, st.p, 200)
        ps = np.linspace(st.p, pm, 200)
        col = "tab:red" if _segno_famiglia(fam) < 0 else "tab:blue"
        ax.plot(curva_onda(st, fam, pe), pe, color=col, ls="--", lw=1.2)
        ax.plot(curva_onda(st, fam, ps), ps, color=col, lw=1.6,
                label=f"{fam} da {st.nome}" if st.nome else fam)
    for st in ss:
        ax.plot(st.u, st.p, "ko", ms=4)
        ax.annotate(st.nome, (st.u, st.p), textcoords="offset points", xytext=(5, 5))
    ax.set_xlabel("u")
    ax.set_ylabel("p")
    ax.set_ylim(0, pm)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    if titolo:
        ax.set_title(titolo)
    return ax


def grafico_p_theta(polari: Sequence = (), stati: Sequence = (), ax=None, p_min_rel: float = 0.05,
                    titolo: str = "", log: bool = False):
    """Piano p–θ (θ in ascissa, gradi) con polari d'urto e d'espansione (F §13).

    ``polari``: coppie (stato, famiglia) con famiglia 'C+', 'C-' o 'entrambe'.
    Ramo debole continuo, forte a tratti fini, espansione tratteggiata;
    il punto di θmax è segnato. ``stati``: stati da segnare in (θ, p).
    """
    plt = _plt()
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4.5))
    colori = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for k, (st, fam) in enumerate(polari):
        pol = polare_p_theta(st, fam, p_min_rel=p_min_rel)
        pol = pol if fam == "entrambe" else {fam: pol}
        c = colori[k % len(colori)]
        for f, d in pol.items():
            u = d["urto"]
            deb = u["ramo"] == "debole"
            ax.plot(u["theta"][deb], u["p"][deb], color=c, lw=1.6,
                    label=f"polare di {st.nome}" if f == list(pol)[0] else None)
            ax.plot(u["theta"][~deb], u["p"][~deb], color=c, lw=1.0, ls=":")
            if "espansione" in d:
                ax.plot(d["espansione"]["theta"], d["espansione"]["p"], color=c, lw=1.2, ls="--")
            tm = d["punti"]["theta_max"]
            ax.plot(*tm, marker="|", color=c, ms=8)
    for st in stati:
        ax.plot(st.theta, st.p, "ko", ms=4)
        ax.annotate(st.nome, (st.theta, st.p), textcoords="offset points", xytext=(5, 5))
    ax.set_xlabel("θ [°]")
    ax.set_ylabel("p")
    if log:
        ax.set_yscale("log")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    if titolo:
        ax.set_title(titolo)
    return ax


def grafico_odografo(polari: Sequence = (), stati: Sequence = (), ax=None, riferimento: str = "assoluto",
                     p_min_rel: float = 0.05, titolo: str = ""):
    """Piano odografo (u*, v*) con polari, cerchio sonico M* = 1 e stati numerati (F §13).

    ``polari``: coppie (stato, famiglia). ``stati``: stati da segnare
    (coordinate V/a* del proprio a*).
    """
    plt = _plt()
    if ax is None:
        _, ax = plt.subplots(figsize=(5.5, 5))
    colori = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for k, (st, fam) in enumerate(polari):
        o = polare_odografo(st, fam, riferimento=riferimento, p_min_rel=p_min_rel)
        fams = ["C+", "C-"] if fam == "entrambe" else [None]
        c = colori[k % len(colori)]
        for f in fams:
            d = o[f] if f else o
            ax.plot(d["urto"]["us"], d["urto"]["vs"], color=c, lw=1.4,
                    label=f"polare di {st.nome}" if f in (None, "C+") else None)
            if "espansione" in d:
                ax.plot(d["espansione"]["us"], d["espansione"]["vs"], color=c, lw=1.0, ls="--")
    t = np.linspace(0, 2 * np.pi, 300)
    ax.plot(np.cos(t), np.sin(t), color="0.5", lw=0.8, ls="-.", label="M* = 1")
    for st in stati:
        th = math.radians(st.theta - (0 if riferimento == "assoluto" else stati[0].theta))
        x, y = st.M_star * math.cos(th), st.M_star * math.sin(th)
        ax.plot(x, y, "ko", ms=4)
        ax.annotate(st.nome, (x, y), textcoords="offset points", xytext=(5, 5))
    ax.set_xlabel("u*")
    ax.set_ylabel("v*")
    ax.set_aspect("equal", adjustable="datalim")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    if titolo:
        ax.set_title(titolo)
    return ax


def grafico_profili(sol: SoluzioneRiemann, t: float, x_range: tuple, grandezze=("p", "rho", "u", "M", "p0"),
                    n: int = 1000, fig=None, titolo: str = ""):
    """Profili spaziali a un istante t delle grandezze scelte ('p', 'rho', 'u', 'a', 'M', 's', 'p0')."""
    plt = _plt()
    x = np.linspace(x_range[0], x_range[1], n)
    c = sol.campiona(x, t)
    if fig is None:
        fig, axs = plt.subplots(len(grandezze), 1, figsize=(6, 1.6 * len(grandezze)), sharex=True)
    else:
        axs = fig.axes
    axs = np.atleast_1d(axs)
    nomi = {"rho": "ρ", "p0": "p0", "M": "M = u/a"}
    for ax, g in zip(axs, grandezze):
        ax.plot(x, c[g], "k", lw=1.3)
        ax.set_ylabel(nomi.get(g, g))
        ax.grid(alpha=0.25)
    axs[-1].set_xlabel("x")
    axs[0].set_title(titolo or f"t = {t:g}")
    fig.tight_layout()
    return fig


# =============================================================================
# 9. Numerica: schema di Godunov con il solutore esatto (S p.41–47)
# =============================================================================
def flusso_godunov(sinistra: State, destra: State) -> np.ndarray:
    """Flusso numerico di Godunov: F(U) dello stato del problema di Riemann esatto in ξ = 0 (S p.42)."""
    st = riemann(sinistra, destra).stato_xi(0.0)
    if st.rho == 0:
        return np.zeros(3)
    return st.F()


def godunov_passo(rho, u, p, gas: Gas, dx: float, cfl: float = 0.9, bc=("trasmissivo", "trasmissivo"),
                  dt: float = None):
    """Un passo dello schema di Godunov 1D (S p.42): Ū_j^(n+1) = Ū_j^n − Δt/Δx (F_(j+1/2) − F_(j−1/2)).

    Δt = cfl·Δx/max(|u|+a) se ``dt`` non è dato (S p.37, p.42). Condizioni al
    contorno con celle fantasma (S p.46–47): 'parete' (ρ, p uguali, u
    ribaltata: u_g = −u_0) o 'trasmissivo' (copia). Restituisce
    (rho, u, p, dt) aggiornati.
    """
    rho, u, p = (np.asarray(x, float) for x in (rho, u, p))
    sts = [State(gas, r, v, q) for r, v, q in zip(rho, u, p)]

    def fantasma(st, tipo):
        if tipo == "parete":
            return st.con(u=-st.u)
        if tipo == "trasmissivo":
            return st
        raise GasdynError("bc: 'parete' o 'trasmissivo'")

    tutti = [fantasma(sts[0], bc[0])] + sts + [fantasma(sts[-1], bc[1])]
    if dt is None:
        dt = cfl * dx / max(abs(s.u) + s.a for s in sts)
    Fi = [flusso_godunov(tutti[i], tutti[i + 1]) for i in range(len(tutti) - 1)]
    U = np.array([s.U() for s in sts])
    Un = U - dt / dx * (np.array(Fi[1:]) - np.array(Fi[:-1]))
    new = [stato_da_U(gas, x) for x in Un]
    return (np.array([s.rho for s in new]), np.array([s.u for s in new]),
            np.array([s.p for s in new]), dt)


def godunov(rho, u, p, gas: Gas, dx: float, t_fine: float, cfl: float = 0.9,
            bc=("trasmissivo", "trasmissivo")):
    """Integra con :func:`godunov_passo` fino a ``t_fine`` (ultimo passo accorciato)."""
    t = 0.0
    while t < t_fine * (1 - 1e-14):
        sts = [State(gas, r, v, q) for r, v, q in zip(rho, u, p)]
        dtc = cfl * dx / max(abs(s.u) + s.a for s in sts)
        rho, u, p, dt = godunov_passo(rho, u, p, gas, dx, bc=bc, dt=min(dtc, t_fine - t))
        t += dt
    return rho, u, p


__all__ = [n for n in list(globals()) if not n.startswith("_") and n not in (
    "annotations", "math", "np", "brentq", "minimize_scalar", "dataclass", "field", "replace",
    "Callable", "Optional", "Sequence")]
