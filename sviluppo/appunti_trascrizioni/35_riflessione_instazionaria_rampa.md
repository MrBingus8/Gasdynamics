# 35 — Lezione 36, 29/11/2021 (file 16): urto piano su rampa (Schardin), diffrazione, Richtmyer–Meshkov

## Detto esplicitamente — metodo
- Lettura delle immagini (Album di Van Dyke n. 235–236, Schardin): "raccogliere immagini senza una chiave
  interpretativa è come **collezionare santini** nella speranza che risolvano i problemi"; serve la teoria per
  interrogare l'immagine, e l'immagine "risponde alla teoria".
- **Tolleranza all'esame**: nel problema di Riemann un residuo sulla velocità intermedia ~10⁻³ "è fisicamente
  equivalente a zero".
- In un serbatoio infinito p0, T0 costanti; all'imbocco condizioni quasi stazionarie in prima approssimazione.
- a = √(γp/ρ) è puntuale; le variazioni attraverso un'onda isentropica si calcolano con Poisson
  p2/p1 = (a2/a1)^(2γ/(γ−1)).
- Successione di infiniti urti deboli → limite isentropico (massimo recupero di pressione).

## Esempio fenomenologico — urto piano (W) su rampa inclinata di φ (rispetto alla verticale)
- Serie a φ = 60°: **RR**; φ = 25°: **MR** (Mach stem quasi dritto, punto triplo).
- **Riduzione a stazionario** (procedura esplicita): il piede risale con V_T = W/cos φ; struttura **autosimile**;
  riferimento solidale al piede (velocità −V_T lungo la parete) ⇒ corrente relativa M_T = (W/a0)/cos φ > 1, urto
  incidente obliquo con σ = 90° − φ ⇒ teoria della riflessione stazionaria (polari).
  φ grande ⇒ M_T alto e σ piccolo ⇒ deflessione modesta ⇒ RR; φ piccolo (σ = 65°) ⇒ deflessione > θ_max ⇒ MR.
- Le relazioni di salto sono invarianti per cambio di riferimento galileiano: "equivalenza tra urto normale
  instazionario e urto obliquo stazionario".
- **Diffrazione** su spigolo divergente (muzzle flow): urto che si indebolisce e si curva restando normale alla parete;
  allo spigolo espansione centrata, **stazionaria (Prandtl–Meyer) se la valle dell'urto è supersonica nel riferimento
  assoluto, instazionaria se subsonica**.
- Urto su contatto inclinato ⇒ **Richtmyer–Meshkov** (supernove, ICF): nel riferimento del punto d'intersezione
  (V_T = W/cos φ) si ha un **problema di Riemann 2D stazionario supersonico**: urto trasmesso, contatto deviato,
  onda riflessa (urto o ventaglio).

## Differenze di programma
- Riflessione instazionaria su rampa, diffrazione, Richtmyer–Meshkov, Schardin: **assenti** da formulario e
  summary (grep negativo). Presenti nel catalogo: **B30-E01** (diffrazione su scalino), **B70-E01** (diffrazione
  attorno a uno spigolo di 105°, M_s = 4), B59-E03 (urto su cavità).
