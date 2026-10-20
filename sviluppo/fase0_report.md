# Fase 0 — Convenzioni, algoritmo di Riemann, discrepanze, struttura di `gasdyn.py`

Fonti lette per intero: `fonti/formulario.tex` (F), `fonti/summary_lectures.tex` (S).
Consultati: `fonti/catalogo.md` (tutti i campi "svolgimento" e le schede 2021–2023 più quelle citate nel prompt),
`fonti/rapporto.md`, e le soluzioni scansionate B63-E02, B71-E02, B72-E02 (per la notazione degli svolgimenti).

## (a) Convenzioni

| Aspetto | Convenzione del corso | Fonte |
|---|---|---|
| Gas perfetto | p = ρRT, cv = R/(γ−1), cp = γR/(γ−1), a² = γp/ρ = γRT, **δ = (γ−1)/2** | F §1 |
| Entropia | variabile entropica **s := p/ρ^γ**; entropia termodinamica Δs_th = cv ln(s2/s1) = cp ln(T2/T1) − R ln(p2/p1) | F §1 |
| Totali/critiche | T0/T = 1+δM², a0² = a²+δV², a*² = 2a0²/(γ+1), M*² = (γ+1)M²/(2+(γ−1)M²) | F §1, §7, §13 |
| Asse 1D | x verso destra, u > 0 verso destra | F §3 |
| Famiglie 1D | λ1 = u−a con **R1 = a−δu** costante lungo λ1; λ2 = u+a con **R2 = a+δu**; λ0 = u (contatto, Ds/Dt = 0). Attenzione: l'indice 1 va con u−a (non come J± di molti libri) | F §3, S p.6 |
| Onde semplici | onda u+a: R1 uniforme, da = δ du; onda u−a: R2 uniforme, da = −δ du; completamento isentropico p/pi = (a/ai)^(2γ/(γ−1)) | F §4 |
| Pistone | pistone a sinistra, gas a destra ⇒ onde u+a; u_p velocità del pistone (espansione u_p = −v_p) | F §5, S p.9 |
| Urti mobili | riferimento solidale v = u−w, **M_r = (u−w)/a con segno**; 1 = monte (davanti), 2 = valle; [q] = q_valle − q_monte; nelle formule M_r1² ("usare il modulo"). Urto u+a: w > u1, M_r1 < 0; urto u−a: M_r1 > 0. "Mach dell'urto" = \|M_r1\| | F §7, S p.18–20; soluzioni B63-E02 ("M = −w/a1"), B72-E02 ("M2r − M2a + w/a2 = 0") |
| Riemann | **1** (sinistra) \| onda u−a \| **2** \| contatto \| **3** \| onda u+a \| **4** (destra); p2 = p3 = p*, u2 = u3 = u*, w_C = u* | F §8 |
| Configurazioni | tre lettere sinistra–contatto–destra: RCR, RCS, SCR, SCS (F); con onda nulla RCN, SCN, NCR, NCS (S p.26–28); "RNS" senza contatto (S p.34). Il vuoto non ha sigla nelle fonti | F §8, S |
| Altri nomi di stati | pistone: 0 indisturbato, 1 dietro l'onda (S p.9, p.20); interazioni successive numerate 5, 6, 7… (S p.34; B04-E01 "u6=u7"; B41-E03 "u6=u7") | S, catalogo |
| Adimensionalizzazioni delle tracce | p=1, a=1 (⇒ ρ=γ): B61, B71, B72; p=1, ρ=γ: B69-E02; p=1, ρ=1 (⇒ a=√γ): B49-E01, B04-E01; p=1, a=√γ: B68-E01. T disponibile solo se R è dato | catalogo |
| Unità SI nelle tracce | Pa, K, m/s; spesso atm (B41-E03, B51-E01) o bar (B44-E03, B64-E01) | catalogo |
| 2D: angoli | σ = angolo d'urto rispetto alla velocità a monte (**σ, non β**), θ direzione/deviazione del flusso (positiva antioraria dall'asse x), μ = asin(1/M), ν(M) di Prandtl–Meyer | F §10–12 |
| 2D: famiglie | C+: dy/dx = tan(θ+μ), θ−ν costante lungo C+; C−: dy/dx = tan(θ−μ), θ+ν costante lungo C−. Attraverso un'onda C+ (fronte θ+μ, nasce da una parete inferiore) vale Δ(θ+ν) = 0 | F §10–11, S p.51–53 |
| Polari | p–θ di Busemann parametrizzata in σ ∈ [μ1, π/2]; rami: onda nulla, debole, θmax, punto sonico, forte, urto normale; slip line: p e θ continui | F §13 |
| Odografo | F e S: piano (V_n/a*, V_t/a*) solidale all'urto; le tracce (B51-E01) chiedono (u*, v*) = (u, v)/a* assoluti | F §13, S p.56 |
| Riflessioni | 1 monte, 2 dietro l'urto incidente, 3 dietro il riflesso; in MR 4 dietro lo stelo (S p.59: "M4<1 e M3>1") | F §14, S p.58–59 |
| Spinta | positiva nel verso del getto, F = ṁ ue + (pe − pa)Ae | F §17 |
| Contorni | tante condizioni quante caratteristiche entranti | F §9 |

## (b) Algoritmo di Riemann come lo leggo nel formulario (§8)

0. Dati: stati 1 (sinistra) e 4 (destra) con (ρ, u, p), γ.
1. **Innesco** (due onde isentropiche): z = (a4/a1)(p1/p4)^((γ−1)/(2γ)), u⁰ = [z(a1+δu1) − (a4−δu4)] / [δ(1+z)].
2. **Onda sinistra u−a** in funzione di u^k:
   - u^k < u1 (urto): k = (γ+1)/4·(u1−u^k)/a1, M_rel,1 = k + √(1+k²); p2 = p1[1 + 2γ/(γ+1)(M²−1)];
     a2 = a1·√{[γ+1+(γ−1)p2/p1]/[γ+1+(γ−1)p1/p2]}; ṗ2 = −2γ(p1/a1)·M³/(1+M²).
   - u^k ≥ u1 (espansione): a2 = a1 − δ(u^k−u1); p2 = p1(a2/a1)^(2γ/(γ−1)); ṗ2 = −γ p2/a2.
3. **Onda destra u+a**, simmetrica: urto se u^k > u4 (k con u^k−u4, ṗ3 = +2γ(p4/a4)M³/(1+M²));
   espansione se u^k ≤ u4 (a3 = a4 + δ(u^k−u4), ṗ3 = +γ p3/a3).
4. **Newton–Raphson**: u^(k+1) = u^k − (p2−p3)/(ṗ2−ṗ3); arresto |1 − p2/p3| < ε.
5. **Controlli**: a2, a3, p2, p3 > 0; se Newton esce dall'intervallo fisico, smorzare il passo.
6. **Ricostruzione**: p* = p2 = p3, u* = u2 = u3; ρ2, ρ3 isentropici (R) o di Hugoniot (S); T, a.
   Velocità: urti w_L = u1 − a1 M_rel,1, w_R = u4 + a4 M_rel,4; ventagli testa/coda u1−a1 / u*−a2 e u4+a4 / u*+a3; contatto u*.
7. **Vuoto** (doppia espansione): u4 − u1 ≥ (a1+a4)/δ.

Dettagli di implementazione che propongo (da confermare):
- Newton esattamente come sopra, con la storia delle iterazioni restituita (utile per mostrare i passaggi).
- Smorzamento: se u^(k+1) esce da (u4 − a4/δ, u1 + a1/δ) (dove a2, a3 > 0) il passo viene dimezzato.
  Rete di sicurezza: f(u) = p2(u) − p3(u) è strettamente decrescente, quindi se Newton non converge entro N iterazioni
  si usa `brentq` sull'intervallo fisico e il risultato viene marcato come tale.
- Gas diversi: le formule usano γ1 a sinistra e γ4 a destra. La formula di u⁰ vale solo per γ uguali: con γ diversi
  l'innesco è la doppia espansione risolta numericamente (influisce solo sul punto di partenza).
  Condizione di vuoto estesa: u4 − u1 ≥ a1/δ1 + a4/δ4.
- Onda nulla "N" se |p*/p_k − 1| < tolleranza dichiarata; contatto "N" se ρ2 = ρ3 con lo stesso gas.
- Opzione `modello='isentropico'`: entrambe le onde isentropiche (è la soluzione chiusa del passo 1). Serve per i confronti
  isentropico/esatto (es. B51-E01 in 2D) e coincide con il primo passo del summary.
- Campionamento in ξ = x/t incluso il ventaglio transonico: per u+a, u(ξ) = [ξ − (a4 − δu4)]/(1+δ),
  a(ξ) = [a4 − δu4 + δξ]/(1+δ) (F §5); per u−a, u(ξ) = [ξ + (a1 + δu1)]/(1+δ), a(ξ) = [a1 + δu1 − δξ]/(1+δ).

## (c) Discrepanze tra le fonti (non risolte da me: per ciascuna una proposta)

| # | Punto | Formulario | Summary | Proposta |
|---|---|---|---|---|
| 1 | Iterazione di Riemann | Newton su u^k con ṗ analitiche | Secante su v con e_k = \|1 − p2/p3\| (p.25): col modulo l'errore non cambia segno e la secante non è garantita | Newton del formulario |
| 2 | Natura delle onde | u^k < u1 ⇒ urto a sinistra; u^k > u4 ⇒ urto a destra | p.24: "se u2 < v0 espansione", "se u4 > v0 espansione": simboli confusi, il senso coincide leggendo u1 al posto di u2 | Formulario |
| 3 | Urto da pistone impulsivo | w = [(1+δ)vp + √((1+δ)²vp² + 4a0²)]/2 | p.23: w = [−(δ+1)vp + √(…)]/2, che non è radice dell'equazione scritta a p.20 e p.23 | Formulario |
| 4 | Espansione riflessa da parete | Solo procedura (invariante entrante) | p.30: a3 = a2 − δ\|vp\| (corretto) e poi "a3 = a0 + δ\|vp\|" (errato: vale a0 − 2δ\|vp\|) | Invarianti |
| 5 | Serbatoio a p costante dietro un urto | Assente | p.31: a3 = a0·e^(s2−s1), dimensionalmente incoerente con s = p/ρ^γ | Fuori formulario: a3/a0 = (s3/s0)^(1/(2γ)) = exp[Δs_th/(2γ cv)] |
| 6 | Tubo d'urto chiuso con a1 = a4 | Il contatto c'è salvo ρ2 = ρ3 | p.34: "nessun contatto perché a1 = a4 ⇒ s1 = s4": falso se p1 ≠ p4, e l'urto genera comunque entropia | Contatto deciso dal calcolo, mai assunto |
| 7 | Ventaglio da pistone | a − δu = a0 | p.10: "u − δa = a0" (refuso; la riga dopo è corretta) | Formulario |
| 8 | σ a θmax per M → ∞ | Calcolo con dθ/dσ = 0 | p.55: σ ≈ 61.8°; per γ = 1.4 risulta 67.79° (θmax = 45.58° coincide) | Calcolo numerico |
| 9 | Limite di M* | M*² = (γ+1)/(γ−1+2/M²) | p.56: "M → ∞ ⇒ M* → (γ+1)/(γ−1)" manca la radice (p.57 usa correttamente 2.45) | Formulario |
| 10 | Sauer–Meyer | d = 1/√((γ+1)R), u = dξ + (γ+1)d²η²/2, ξ_sonica = −(γ+1)dη²/2 | p.65: due versioni, una senza (γ+1) (d = 1/√R, ξ_s = −dη²/2, "Favini/Franco") e una con (γ+1) da un esercizio 2024 | Formulario, con nota nel docstring |
| 11 | Segno nelle polari | θ2 = θ1 ± θ(σ), "± deviazioni opposte", senza associarlo alla famiglia | p.53: espansione su parete che ruota verso il basso con fronti θ+α | C+ ⇒ θ2 = θ1 + θ(σ) (la compressione ruota in senso antiorario); C− ⇒ θ1 − θ(σ) |
| 12 | Odografo | (V_n/a*, V_t/a*) nel riferimento dell'urto | Idem | Le tracce chiedono (u*, v*): implemento entrambi, ciascuna corrente con il proprio a* |
| 13 | Refusi minori | — | p.53 "dy/dx = θ ± α" senza tan; p.59 valori indicativi θN ≈ 20°, θD ≈ 30° per M ≥ 5 | Nessuna conseguenza |

Limiti interni del formulario (estensioni che propongo, dichiarate nel docstring):
- §8 condizione di vuoto valida per γ uguale; con gas diversi a1/δ1 + a4/δ4.
- §6 t_BU = a0/((1+δ)\|ẍp(0)\|) vale per pistone partito da fermo con u0 = 0. Per una legge generica uso l'inviluppo delle
  caratteristiche (S p.15–16: (1+δ)ẍ(t̂ − tp) − δẋ − a0 = 0), minimizzato numericamente.
- §6 profilo lineare: il blow-up esiste solo se v decresce con x; il modulo lo controlla.

Anticipazioni utili per la Fase 2 (trovate controllando i numeri):
- **B68-E01** (caso Mach assoluto a valle 1.89): con γ = 1.4 e gas in quiete il Mach assoluto a valle è limitato da
  √(2/(γ(γ−1))) = 1.8898, quindi 1.89 non ha soluzione. La soluzione usa 1.87, e questo è coerente.
  Il modulo solleverà un'eccezione esplicita.
- **B65-E01**: t = 0.4246 = 1/2.3555 è il blow-up della famiglia u−a dei profili lineari dati, anche se i dati non sono
  un'onda semplice. La funzione di formazione d'urto calcolerà t_BU per famiglia e segnalerà se l'altro invariante non è uniforme.
- **B49-E01**: t_BU = 0.66622 e x_BU = 0.78827 tornano subito; la "w finale ≈ 2.4375" non torna con il calcolo veloce
  (≈ 2.367 con vp = 1.48). Da indagare in Fase 2.

## (d) Struttura proposta di `gasdyn.py`

Convenzioni d'interfaccia proposte: **angoli in gradi** in tutte le funzioni pubbliche (come nelle tracce); grandezze in
unità coerenti qualsiasi (SI o adimensionali); costante `ATM = 101325.0` Pa; nessun γ di default.

```python
# --- eccezioni ---
class GasdynError(Exception)            # base
class NessunaSoluzione(GasdynError)     # es. M assoluto a valle oltre il limite
class UrtoStaccato(GasdynError)         # θ > θmax
class Vuoto(GasdynError)                # generazione di vuoto dove non ammessa
class Strozzamento(GasdynError)         # A/A* < 1 richiesto

# --- 1. termodinamica ---
@dataclass(frozen=True) class Gas(gamma, R=None, nome='')
    .delta .q .cp .cv ;  Gas.aria(R=287.0)
@dataclass(frozen=True) class State(gas, rho, u, p, v=0.0)
    .a .T .V .M .theta .s(=p/ρ^γ) .R1 .R2 .h .h0 .T0 .p0 .rho0 .a0 .a_star .M_star .V_max
    State.da_pa(gas, p, a, u=0, v=0) / da_pT / da_rhoT / da_M(gas, M, p, *, a|rho|T, theta=0)
    .con(**modifiche)
salto_entropia(s1, s2) -> dict(s2_su_s1, ds_su_R, ds_su_cv, ds)      # F §1
isentropiche(M, gas) -> dict(T0_T, p0_T, rho0_rho, a0_a)             # F §1
mach_da_rapporto(valore, gas, tipo='p0/p'|'T0/T'|'rho0/rho')
pressione_totale(M, gas, modello='isentropico'|'isotermo'|'incomprimibile')   # isot./incompr. fuori formulario
critiche(gas) -> dict(T*/T0, p*/p0, rho*/rho0, a*/a0) ; velocita_max(a0, gas)

# --- 2. quasi-1D ---
area_su_area_critica(M, gas)                                  # F §16
mach_da_area(A_su_Astar, gas, ramo='subsonico'|'supersonico')
portata(A, p0, T0, M, gas) ; portata_max(A_star, p0, T0, gas)
ugello_convergente(p0, T0, pb, Ae, gas) -> regime e stato d'uscita
ugello_cd(p0, pb, At, Ae, gas) -> regime: subsonico | urto nel divergente (As, M1, M2, p02/p01) | sovra/sottoespanso | adattato
spinta(mdot, ue, pe, pa, Ae, ingresso=None)                   # F §17
variazione_area(stato, A1, A2, ramo='stesso')                 # collegamento stazionario ṁ, h0, s
serbatoio_svuotamento(t, p_i, T_i, V, A, gas) ; serbatoio_riempimento(...)   # F §17 (proposta)

# --- 3. urto normale stazionario ---
urto_normale(M1, gas) -> RisultatoUrto(p2_p1, rho2_rho1, T2_T1, M2, p02_p01, ds_su_R)   # F §7, §12
urto_normale_stato(stato) -> State

# --- 4. urti mobili ---
@dataclass class Urto(famiglia, monte, valle, w, Mr_monte, Mr_valle)
urto_mobile(monte, famiglia, *, w=None, Ms=None, p_valle=None, u_valle=None, M_valle=None) -> Urto
urto_riflesso_parete(valle_incidente, famiglia_riflesso, u_parete=0.0) -> Urto
natura_salto(stato_a, stato_b) -> 'urto u+a' | 'urto u−a' | 'contatto' | 'onda semplice …' | 'incompatibile'   # B48-E02
verifica_rh(urto) -> residui di massa, quantità di moto, energia

# --- 5. onde semplici 1D ---
punto_caratteristiche(A, B, gas)                               # F §3
curva_onda(stato, famiglia, p) -> u(p)  (ramo R e S, F §8)
onda_semplice(stato, famiglia, *, u=None|p=None|a=None) -> State
@dataclass class Ventaglio(famiglia, monte, valle) .xi_testa .xi_coda .stato(xi) .transonico
pistone_impulsivo(stato0, u_p, lato='sinistro') -> Urto | Ventaglio (con controllo del vuoto)
pistone_legge(x_p, stato0, t_max, lato='sinistro', dx_p=None, ddx_p=None) -> caratteristiche, t_BU, x_BU, tp*
formazione_urto_ci(x, u0, a0, gas, famiglia='auto') -> t_BU, x_BU, xi*, avvisi   # F §6
t_bu_lineare(x1, v1, x2, v2)

# --- 6. Riemann 1D ---
@dataclass class SoluzioneRiemann(stati{1,2,3,4}, p_star, u_star, onde[L, C, R], configurazione, vuoto, iterazioni)
    .campiona(x, t) / .stato_xi(xi) / .velocita_onde()
riemann(sinistra, destra, *, tol=1e-12, maxit=50, modello='esatto'|'isentropico') -> SoluzioneRiemann   # F §8
onda_verso(stato, famiglia, *, u=None|p=None) -> Urto | Ventaglio     # parete (u=0), pistone (u=up), pressione imposta
bc_serbatoio(stato_interno, lato, p0, T0, gas)                        # S p.44–45 (proposta)
interazione(stato_a, stato_b) -> SoluzioneRiemann                     # pseudo-Riemann, F §9 (alias documentato)

# --- 7. supersonico 2D ---
angolo_mach(M) ; prandtl_meyer(M, gas) ; mach_da_pm(nu, gas)                 # F §11
theta_da_sigma(M, sigma, gas) ; sigma_da_theta(M, theta, gas, ramo='debole'|'forte')   # UrtoStaccato se θ>θmax
theta_max(M, gas) -> (theta_max, sigma) ; sonico(M, gas) -> (theta_s, sigma_s)          # F §12
urto_obliquo(stato, *, theta=None|sigma=None, famiglia='C+'|'C-', ramo='debole') -> (State, info)
espansione_pm(stato, dtheta, famiglia) ; compressione_isentropica(stato, dtheta, famiglia)
polare_p_theta(stato, famiglia='entrambe', n=400, espansione=True) -> array etichettati per ramo
polare_odografo(stato, piano='uv'|'nt')
onda_verso_p_2d(stato, p, famiglia, modello='esatto'|'isentropico') ; onda_verso_theta_2d(...)
riemann_stazionario(sopra, sotto, modello='esatto'|'isentropico') -> p, θ, onde, stati   # B51-E01, B74-E01
intersezioni_polari(A, famA, B, famB) -> tutte le intersezioni etichettate per ramo
riflessione_regolare(M1, theta_w, gas) ; riflessione_mach(M1, theta_w, gas)  # tre urti
theta_von_neumann(M1, gas) ; theta_distacco(M1, gas) ; diagramma_rr_mr(M_range, gas)
shock_expansion(vertici, M, alpha, gas) -> Cp per pannello, CL, CD           # lastre e profili
pitot_supersonico(M1, gas) ; mach_da_pitot(p02_su_p1, gas)                   # Rayleigh, fuori formulario
punto_caratteristiche_2d(A, B, gas)                                          # F §10
getto_adattamento(Me, pe, pa, gas)                                           # F §15
sauer_meyer(R, gas)                                                          # F §19

# --- 8. grafici (matplotlib, stile sobrio) ---
grafico_xt(sol, t_max, ...) ; grafico_pu(stati, curve=..., etichette=...)
grafico_p_theta(polari, stati, ...) ; grafico_odografo(...) ; grafico_profili(sol, t, grandezze=('p','rho','u','M','p0'))

# --- 9. numerica ---
flusso_godunov(sinistra, destra) ; godunov_passo(U, dx, cfl, bc=('parete'|'trasmissivo', ...)) -> U_nuovo, dt
```

Ogni oggetto restituito porta le grandezze con i nomi del corso (p, ρ, u, a, M, M_r, w, σ, θ, ν) e un `__repr__` leggibile.

## Proposte aggiuntive (calcoli ricorrenti nelle tracce, non nell'elenco del prompt)

- Serbatoi: svuotamento strozzato e riempimento adiabatico (F §17; B03-E01, B37-E01, B40-E02).
- Contorno a serbatoio in ingresso/uscita (S p.44–45; B04-E01, B37-E01).
- Punto di intersezione delle caratteristiche 2D e ugello di lunghezza minima planare (F §10; B01-E02, B15-E03).
- Sauer–Meyer in gola (F §19; B67-E01).
- Relazioni di salto omoterme (p = a²ρ con a costante; B62-E02), fuori formulario.
- Identificazione della natura di un salto tra due stati dati (B48-E02).
- Diagramma θ–M1 delle zone RR/MR (F §14).

## Punti su cui serve la tua decisione

1. **git non è installato.** Posso installarlo con `winget install --id Git.Git` (modifica di sistema: serve il tuo ok), oppure lo installi tu.
2. Discrepanze 1–12: confermi le proposte? In particolare la 1 (Newton), la 10 (Sauer–Meyer) e la 11 (segno delle polari).
3. Angoli in gradi nell'interfaccia pubblica.
4. `Gas.aria()` con R = 287.0 J/(kg K) di default (B64-E01 usa 287, B73-E01 usa 288; R è sempre sovrascrivibile).
5. Sigla del vuoto: propongo "RVR" (R–vuoto–R), non presente nelle fonti.
6. Quali proposte aggiuntive includere.
