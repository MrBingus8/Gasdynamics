# 10 — Lezione 10, 12/10/2021 (file 32): Cauchy–K. 1D, autovettori destri, flux splitting

## Detto esplicitamente
- Il fine della teoria delle caratteristiche è passare da IVP a **IBVP** (condizioni iniziali + al contorno):
  "Le sole equazioni differenziali non bastano"; il flusso reale è sempre un problema misto.
- Distingue sempre discontinuità **deboli** (sulle derivate) da **forti** (urti, contatti: forma integrale,
  Rankine–Hugoniot, aumento di entropia).
- Riduzione 1D: v è disaccoppiata e trasportata lungo la traiettoria.
- Autovalori "ordinati rigorosamente in ordine crescente": λ1=u−a, λ2=u (entropia), λ3=u (scorrimento), λ4=u+a.
  Compatibilità: dp − ρa du = 0 lungo u−a; ds=0 e dv=0 lungo u; dp + ρa du = 0 lungo u+a.
- Lettura di R: p e u vedono solo le due acustiche; ρ vede acustiche + entropia; v solo se stessa.
- Problema a 4 regioni: 8 incognite (derivate temporali in 2 e 3) = 3 (attraverso u−a) + 2 (traiettoria:
  p_t e u_t continue) + 3 (attraverso u+a) ⇒ ben posto. Poi si guarda **dove cade l'asse t** (in 1,2,3 o 4) per
  scegliere ∂U/∂t: è il germe dello schema upwind/Godunov.
- Steger–Warming: Λ± = (Λ ± |Λ|)/2, derivata a monte per le onde con λ>0, a valle per λ<0.

## Convenzioni — confronto con il formulario
- Numerazione delle **onde** per autovalore crescente 1..4 (u−a, u, u, u+a) — nel formulario §3 invece
  λ1=u−a, λ0=u, λ2=u+a. Differenza solo di etichetta.
- Regioni: 1 | u−a | 2 | u | 3 | u+a | 4 (coerente col formulario: 1 a sinistra). La frase "Stato 1 a destra di x̄,
  Stato 4 a sinistra" ripetuta anche qui è **incoerente** con la sequenza 1–(u−a)–2 descritta subito dopo:
  quasi certamente un errore della trascrizione (o un lapsus); la struttura vera è 1 a sinistra, 4 a destra.

## Esempi d'esame
- Collegabile alle domande teoriche/numeriche "schema upwind / splitting" (B05-E02, B07-E02, B11-E05, B14-E02).
