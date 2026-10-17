# Guida a `gasdyn.py`

Modulo di calcolo per il corso di Gasdinamica (Favini). Serve a ottenere numeri corretti e a
risolvere in modo esatto i calcoli iterativi (Riemann, polari, formazione d'urto). Il
ragionamento fisico (quale onda, quale modello, quale stato) resta a chi risolve l'esercizio.
Notazione e algoritmi sono quelli del formulario del corso (citato nei docstring come "F §n") e
del riassunto delle lezioni ("S p.n"). Ogni funzione ha un docstring in italiano: in caso di
dubbio usare `help(gd.nome_funzione)`.

```python
import gasdyn as gd
aria = gd.Gas(1.4)            # adimensionale (T non disponibile)
aria_SI = gd.Gas(1.4, 287.0)  # oppure gd.Gas.aria(R=287.0)
```

Dipendenze: numpy, scipy, matplotlib (solo per i grafici).

## 1. Convenzioni (da rispettare sempre)

| Cosa | Convenzione |
|---|---|
| Gas | `Gas(gamma, R=None)`: **non esiste un γ implicito**. δ = (γ−1)/2 (`gas.delta`), q = (γ−1)/(2γ) (`gas.q`). |
| Unità | Qualsiasi sistema coerente. Dati tipici delle tracce: p = 1, a = 1 ⇒ ρ = γ (`State.da_pa(gas, 1, 1)`); p = 1, ρ = 1 ⇒ a = √γ. T richiede R. `gd.ATM = 101325` Pa (alcuni svolgimenti usano 10⁵). |
| Entropia | `stato.s` = p/ρ^γ (variabile entropica del corso); `salto_entropia` dà Δs/R e Δs/c_v. |
| 1D | x verso destra, u > 0 verso destra. Famiglie: **u−a** (lungo cui è costante R1 = a−δu, `stato.R1`) e **u+a** (R2 = a+δu, `stato.R2`). Un'onda semplice u+a ha R1 uniforme; un'onda u−a ha R2 uniforme. |
| Urti mobili | v = u − w, **M_r = (u−w)/a con segno**; monte = stato davanti all'urto, valle = dietro. Urto u+a (va a destra rispetto al gas): M_r monte < 0. `Ms` = \|M_r monte\|. |
| Riemann | stato **1** a sinistra, **4** a destra; 1 \| onda u−a \| 2 \| contatto \| 3 \| onda u+a \| 4; p2 = p3 = p*, u2 = u3 = u*. Sigle: R (espansione), S (urto), N (onda nulla), C (contatto), vuoto `RVR`. |
| 2D | Angoli in **gradi**; θ = direzione del flusso (antiorario dall'asse x); σ = angolo d'urto rispetto alla velocità a monte (non β); μ = angolo di Mach; ν = Prandtl–Meyer. |
| Famiglie 2D | Onda **C+** = fronte inclinato θ+μ (nasce da una parete *sotto* il flusso e sale): attraverso di essa Δ(θ+ν) = 0, la compressione **aumenta** θ, polare θ2 = θ1 + θ(σ). Onda **C−** = fronte θ−μ (nasce da una parete sopra): la compressione **diminuisce** θ. Due correnti affiancate: la superiore emette C+, l'inferiore C−. |
| Rami | Urto obliquo: `ramo='debole'` (default, campo aperto) o `'forte'`. |

**Errori invece di numeri sbagliati.** Il modulo solleva eccezioni esplicite:
`NessunaSoluzione` (es. urto con M a valle oltre il limite), `UrtoStaccato` (θ > θmax),
`Vuoto` (a ≤ 0), `Strozzamento` (A/A* < 1), `RamiMultipli` (più soluzioni, elencate in
`e.soluzioni`). Non catturarle per "forzare" un risultato: indicano un problema fisico.

## 2. Stati e oggetti restituiti

- `State(gas, rho, u, p, v=0, nome='')` — costruttori `State.da_pa(gas, p, a, u=0)`,
  `State.da_pT(gas, p, T, u=0)`, `State.da_M(gas, M, theta=0, p=..., a=|rho=|T=)`,
  `State.da_totali(gas, M, p0=..., a0=|T0=, theta=0)`; `st.con(u=0, nome='5')` copia con modifiche.
  Proprietà: `a, T, V, M, Mu (=u/a con segno), theta, mu, s, R1, R2, h, h0, a0, T0, p0, rho0,
  a_star, M_star, V_max, U(), F()`.
- `Urto` (1D): `famiglia, monte, valle, w, Ms, Mr_monte, Mr_valle, p2_p1, rho2_rho1, sigla, residui_rh()`.
- `Ventaglio` (1D): `famiglia, monte, valle, xi_testa, xi_coda, transonico, vuoto, stato_xi(xi)`.
- `SoluzioneRiemann`: `stati['1'..'4'], p_star, u_star, onda_sx, onda_dx, configurazione, metodo,
  iterazioni (k, u^k, p2^k, p3^k), velocita(), stato_xi(xi), campiona(x, t)`.
- `Onda2D`: `tipo ('urto'|'espansione'|'compressione isentropica'), famiglia, monte, valle, sigma,
  ramo, deflessione, angolo_fronte, angolo_testa, angolo_coda, Mn1, p2_p1, p02_p01`.

Tutti gli oggetti hanno un `print()` leggibile con gli stati.

## 3. Funzioni ed esempi per area

### 3.1 Termodinamica e grandezze totali

`isentropiche(M, gas)` → T0/T, p0/p, ρ0/ρ, a0/a · `mach_da_rapporto(valore, gas, tipo='p0/p'|'T0/T'|'rho0/rho'|'a0/a')` ·
`critiche(gas)` · `velocita_max(a0, gas)` · `pressione_totale(M, gas, modello='isentropico'|'incomprimibile'|'isotermo')` ·
`isentropica(stato, p=|rho=|a=|T=)` (stessa s) · `salto_entropia(st1, st2)` · `stato_da_U(gas, U)`.

```python
r = gd.isentropiche(2.0, aria)                       # r['p0_p'] = 7.824
M = gd.mach_da_rapporto(7.824, aria, "p0/p")         # 2.0
gd.pressione_totale(0.3, aria, "isotermo")           # 1.0650 (B43-E01)
st = gd.State.da_pT(aria_SI, 1e5, 300.0, u=200.0)
print(st.T0, st.p0, st.a_star)
gd.salto_entropia(st, st.con(p=2e5))["ds_R"]          # (s2−s1)/R
```

### 3.2 Isentropico, quasi-1D, ugelli, serbatoi

`area_su_area_critica(M, gas)` · `mach_da_area(A_su_Astar, gas, ramo='subsonico'|'supersonico')` ·
`mach_da_area_tra_sezioni(M1, A1, A2, gas, ramo='stesso')` · `portata(A, p0, T0, M, gas)` ·
`portata_max(A_star, p0, T0, gas)` · `ugello_convergente(p0, pb, gas)` ·
`ugello_cd(p0, pb, Ae_At, gas)` (regime e posizione dell'urto nel divergente) ·
`urto_in_ugello(As_At, Ae_At, gas)` · `spinta(mdot, ue, pe, pa, Ae, ingresso=None)` ·
`variazione_area(stato, A1, A2, ramo='stesso')` (collegamento stazionario ṁ, h0, s) ·
`miscelamento([(stato, A), ...])` · `serbatoio_svuotamento(t, p_i, a_i, V, A, gas, p_a=None)` ·
`serbatoio_riempimento(V, p_i, T_i, T0_in, gas, m=|p=)`.

```python
gd.mach_da_area(3.0, gd.Gas(1.3), "supersonico")     # 2.5141 (B02-E01)
r = gd.ugello_cd(1.0, 0.7, 2.0, aria)                # pb/p0 = 0.7, Ae/At = 2
print(r["regime"], r["As_At"], r["M1"])              # urto nel divergente
st5 = gd.variazione_area(gd.State.da_M(aria, 1.87, p=380.0, rho=8.0), 1.5, 1.0)  # contrazione, ramo supersonico
```

### 3.3 Urto normale stazionario

`urto_normale(M1, gas)` → dict p2_p1, rho2_rho1, T2_T1, a2_a1, M2, p02_p01, s2_s1, ds_R, p02_p1 ·
`urto_normale_stato(stato)` · `mach_da_urto(valore, gas, tipo='p2/p1'|'rho2/rho1'|'T2/T1'|'M2'|'p02/p01')`.

```python
r = gd.urto_normale(2.0, aria)    # p2/p1 = 4.5, M2 = 0.5774, p02/p01 = 0.7209
gd.mach_da_urto(0.7209, aria, "p02/p01")
```

### 3.4 Urti mobili 1D

`urto_mobile(monte, famiglia, w=|Ms=|p_valle=|u_valle=|M_valle=)` (un solo dato; `M_valle` = u2/a2
con segno, risolto cercando tutte le radici) · `urto_riflesso_parete(stato_dietro_incidente,
famiglia_riflesso, u_parete=0)` · `M_valle_limite(gas)` · `natura_salto(A, B)` (quale onda collega due
stati) · `urto_isotermo(M1)` (fuori formulario).

```python
q = gd.State.da_pa(aria, 1.0, 1.0)                    # gas fermo, p = 1, a = 1
inc = gd.urto_mobile(q, "u+a", M_valle=1.87)          # B72-E02: w = 17.99
rif = gd.urto_riflesso_parete(inc.valle, "u-a")      # parete a destra
pist = gd.urto_mobile(q, "u+a", u_valle=2.0)          # pistone a 2: w = 2.762
print(inc, rif.valle.p)
gd.M_valle_limite(aria)                               # 1.8898: oltre ⇒ NessunaSoluzione
```

### 3.5 Onde semplici, pistoni, formazione dell'urto

`onda_semplice(stato, famiglia, u=|p=|a=|rho=)` · `curva_onda(stato, famiglia, p)` (curve p–u del formulario) ·
`punto_caratteristiche(A, B)` (u+a da A, u−a da B) · `onda_verso(stato, famiglia, u=|p=, consenti_vuoto=False,
modello='esatto')` → `Urto` o `Ventaglio` · `pistone_impulsivo(stato0, u_p, lato='sinistro'|'destro')` ·
`pistone_legge(x_p, stato0, t_max, lato='sinistro')` → `PistoneLegge` (t_bu, x_bu, stato(x, t),
caratteristica(tp)) · `formazione_urto_ci(u0, a0, gas, intervallo=(x1, x2))` (o `x=` con array) ·
`t_bu_lineare(x1, v1, x2, v2)`.

Attenzione alla famiglia: un pistone **a sinistra** del gas (o un contorno sinistro) genera onde
**u+a**; a destra onde **u−a**.

```python
q = gd.State(aria, 1.0, 0.0, 1.0)                     # p = 1, ρ = 1
P = gd.pistone_legge(lambda t: 0.74 * t**2, q, 1.0)   # B49-E01
print(P.t_bu, P.x_bu)                                 # 0.66623, 0.78829
V = gd.pistone_impulsivo(gd.State.da_pa(aria, 1, 1), -1.0)   # espansione, coda transonica
print(V.xi_testa, V.xi_coda, V.transonico)
r = gd.formazione_urto_ci(lambda x: 3 - 2.5*x**2 + 1.667*x**3,
                          lambda x: 1 + 0.5*x**2 - 0.3334*x**3, aria, (0, 1))   # B64-E03
print(r["prima"], r["u-a"]["t_bu"], r["u-a"]["x_bu"], r["u-a"]["onda_semplice"])
```

### 3.6 Problema di Riemann 1D, contorni, interazioni

`riemann(sinistra, destra, modello='esatto'|'isentropico', x0=0, t0=0)` — algoritmo del formulario
(innesco u⁰ con z, Newton su u con ṗ analitiche; `brentq` solo come rete di sicurezza, segnalato in
`metodo`); gas diversi ai due lati ammessi · `interazione(sx, dx)` (pseudo-Riemann fra gli stati
esterni, F §9) · `contorno_parete(stato, lato_parete, u_parete=0)` · `contorno_pistone(stato, lato, u_p)` ·
`contorno_pressione(stato, lato, p_b)` · `contorno_serbatoio(stato, lato, p0, a0=|T0=)` (ingresso con
urto + contatto oppure uscita a p0) · `riemann_area(sx, dx, A_sx, A_dx)` (urto/flusso attraverso una
brusca variazione d'area, flusso verso destra).

```python
sod = gd.riemann(gd.State(aria, 1, 0, 1), gd.State(aria, 0.125, 0, 0.1))
print(sod)                         # RCS, p* = 0.30313, u* = 0.92745
sod.iterazioni[0]                  # (0, u⁰ = 0.918947, p2⁰, p3⁰) come nello svolgimento B04-E01
c = sod.campiona([-0.3, 0.0, 0.3], t=0.2)       # dict di array rho, u, p, a, M, s, p0
o = gd.contorno_parete(sod.stati["3"], "destra")   # urto riflesso da una parete a destra
q = gd.State.da_pa(aria, 1, 1)
r = gd.riemann_area(gd.urto_mobile(q, "u+a", M_valle=1.87).valle, q, 1.5, 1.0)   # B72-E02
print(r["riemann_valle"].configurazione, r["riemann_valle"].p_star)
```

Per una catena di interazioni: risolvere un Riemann, prendere gli stati che si incontrano
(`sol.stati[...]` o `onda.valle`) e risolvere il successivo con `interazione`. Per il diagramma
x–t usare `x0`, `t0` del punto d'incontro.

### 3.7 Supersonico stazionario 2D

`angolo_mach(M)` · `prandtl_meyer(M, gas)` · `mach_da_pm(nu, gas)` · `nu_max(gas)` ·
`theta_da_sigma(M, sigma, gas)` · `sigma_da_theta(M, theta, gas, ramo)` · `theta_max(M, gas)` → (θmax, σ) ·
`sonico(M, gas)` → (θ_s, σ_s) · `urto_obliquo(stato, theta=|sigma=, famiglia='C+'|'C-', ramo)` ·
`espansione_pm(stato, dtheta, famiglia)` · `onda_isentropica_2d(stato, famiglia, theta2=|M2=|p2=)` ·
`onda_verso_theta_2d(stato, theta2, famiglia, ramo, modello)` (condizione di parete) ·
`onda_verso_p_2d(stato, p2, famiglia, modello)` (contorno a pressione assegnata, getti) ·
`polare_p_theta(stato, famiglia='entrambe')` · `polare_odografo(stato, famiglia, riferimento='assoluto'|'monte')` ·
`intersezioni_polari(A, famA, B, famB)` (tutte, etichettate) · `riemann_stazionario(sopra, sotto,
modello, ramo='debole')` · `riflessione_regolare(stato1, theta_w)` · `riflessione_mach(stato1, theta_w)` ·
`theta_von_neumann(M1, gas)` · `theta_distacco(M1, gas)` · `diagramma_rr_mr(M_valori, gas)` ·
`shock_expansion(superiore, inferiore, M, alpha, gas)` · `pitot_supersonico(M1, gas)` · `mach_da_pitot(r, gas)` ·
`punto_caratteristiche_2d(A, B, posA, posB)` · `ugello_minima_lunghezza(Me, gas, n)` ·
`getto_adattamento(Me, pe, pa, gas)` · `sauer_meyer(R, gas, versione='formulario'|'favini')`.

```python
gd.sigma_da_theta(2.0, 10.0, aria)              # 39.31° (ramo debole)
gd.prandtl_meyer(2.0, aria)                     # 26.38°
s1 = gd.State.da_M(aria, 2.0, p=1.0, a=1.0)
o2 = gd.urto_obliquo(s1, theta=11.31)           # B71-E01: rampa sotto il flusso ⇒ C+
o3 = gd.onda_verso_theta_2d(o2.valle, 0.0, "C+")  # espansione allo spigolo, ritorno a θ = 0
G = gd.Gas(1.4, 287.0)
sotto = gd.State.da_M(G, 6.0, p=0.01, T=260.0)
sopra = gd.State.da_M(G, 4.0, p=0.50, T=350.0)
r = gd.riemann_stazionario(sopra, sotto)        # B51-E01: θ = −16.58°, p = 0.0699
rr = gd.riflessione_regolare(gd.State.da_M(aria, 2.8, p=1, a=1), 10.0)
print(rr["3"].p, rr["angolo_riflesso_parete"])
gd.theta_distacco(1.3, aria)                    # 3.25° (B69-E01)
```

### 3.8 Grafici

`grafico_xt(sol | [sol, (sol, t_fine), ...], t_max, x_range, linee=[((x0,t0),(x1,t1),'urto'|...)], pistone=f)` ·
`grafico_pu(sol, stati=[...], curve=[(stato, famiglia)])` · `grafico_p_theta([(stato, famiglia)], [stati])` ·
`grafico_odografo([(stato, famiglia)], [stati])` · `grafico_profili(sol, t, x_range, grandezze=('p','rho','u','M','p0'))`.
Le etichette degli stati sono i loro `nome`. Restituiscono `ax` (o `fig`): salvare con
`ax.figure.savefig(...)`.

```python
import matplotlib; matplotlib.use("Agg")
ax = gd.grafico_xt(sod, 0.25, (-0.5, 0.5), titolo="Sod")
ax2 = gd.grafico_pu(sod)
ax3 = gd.grafico_p_theta([(sopra.con(nome="sup"), "C+"), (sotto.con(nome="inf"), "C-")], log=True)
```

### 3.9 Numerica

`flusso_godunov(sx, dx)` · `godunov_passo(rho, u, p, gas, dx, cfl=0.9, bc=('trasmissivo'|'parete', ...), dt=None)` →
(rho, u, p, dt) · `godunov(rho, u, p, gas, dx, t_fine)`.

```python
s1 = gd.State(aria, 1.4, 0.0, 1.0)
s2 = gd.onda_semplice(s1, "u-a", a=2/3)          # espansione transonica simmetrica (B38-E02)
rho, u, p, dt = gd.godunov_passo([1.4, 1.4, s2.rho], [0, 0, s2.u], [1, 1, s2.p], aria, 1.0, cfl=0.8)
```

## 4. Funzioni fuori formulario

Dichiarate nei docstring; usarle sapendo che non sono nelle fonti del corso:
`pressione_totale` (modelli incomprimibile e isotermo), `M_valle_limite`, `urto_isotermo`,
`pitot_supersonico`/`mach_da_pitot` (Rayleigh), `shock_expansion`, `ugello_minima_lunghezza` (MOC completo),
`riemann_area` nel caso subsonico (estensione del procedimento del corso), `contorno_serbatoio` con
urto in ingresso (S p.45 lo tratta solo isentropico), innesco del Riemann con γ diversi.

## 5. Discrepanze note fra le fonti

Il modulo segue il formulario; le varianti sono in `fase0_report.md`. Le due che si vedono negli
svolgimenti:
- **Sauer–Meyer**: default con (γ+1) come nel formulario; lo svolgimento del 2022 (B67-E01) usa la
  versione senza (γ+1): `sauer_meyer(R, gas, versione='favini')`.
- **Iterazione di Riemann**: Newton del formulario (il summary propone una secante).

Il confronto con gli svolgimenti d'esame (affidabili, con riserva, errati) è in
`verifica/verifica_catalogo.md`: prima di usare uno svolgimento come riferimento numerico,
controllare lì il suo giudizio.

## 6. Limiti noti

- Gas perfetto caloricamente, flusso non viscoso. Strato limite, attrito e dissociazione non sono trattati.
- 1D piano. La sezione variabile si tratta solo con il collegamento stazionario fra due sezioni
  (`variazione_area`, `riemann_area`). In `riemann_area` i casi misti (sub/supersonico) sollevano un errore.
- `riemann(..., modello='isentropico')`: le compressioni sono isentropiche, quindi `campiona` non è definito.
- `formazione_urto_ci` è esatta solo per un'onda semplice (`onda_semplice=True` nel risultato); altrimenti è una stima.
- `pistone_legge` vale fino a t_BU. Gas iniziale uniforme, pistone senza salto impulsivo iniziale di velocità.
  Dopo il blow-up servono gli strumenti di Riemann (pseudo-Riemann, regime asintotico).
- `riflessione_mach`: teoria dei tre urti a stelo rettilineo; `principale` è la soluzione con |θ| minima,
  le altre sono in `soluzioni`.
- `shock_expansion`: trascura le onde riflesse e la scia, nessuna pressione di base; profili a spigoli vivi.
- `ugello_minima_lunghezza`: piano, discretizzazione con O(1/n) di errore; con n grandi è lento (n ≈ 40 basta).
- Polari e intersezioni: θ è funzione a un sol valore di p per ciascuna famiglia; i punti di contatto
  tangenti (es. criterio di von Neumann) sono trattati con una tolleranza di 10⁻⁶°.
- Schema di Godunov al primo ordine. Contorni solo 'parete' e 'trasmissivo'.
- Grafici essenziali: per figure d'esame complesse comporre più chiamate sullo stesso `ax`.
