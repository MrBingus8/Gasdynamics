# 43 — Lezioni 46–47, 20/12/2021 (file 01): viscosità turbolenta/LES; esame 14/09/2021; getto con ago

## Turbolenza e numerica (fuori dal programma del summary)
- ν_T: trasferimento di q.d.m. dalle scale sub-grid; senza ν_T l'energia della cascata si accumula a k = π/Δx ⇒
  esplosione numerica ("buco di scarico"). Δx → scala di Kolmogorov ⇒ ν_T → 0, LES → DNS.
- Modelli: lunghezza di mescolamento di Prandtl, Smagorinsky, **ILES** (schema di alto ordine, errore di troncamento
  diffusivo come ν_T implicita; es. studio ONERA C1xb, motori a solido con soffiaggio).
- Deluge systems sulle rampe di lancio (Vega C, Ariane 5): abbattimento acustico; getti caldi multifase.
→ **Solo nelle trascrizioni** (LES/ILES/turbolenza non sono nel summary né nel formulario).

## Esempio d'esame — 14/09/2021 = **B65-E01** (catalogo: problema 1D con dati iniziali isentropici, grafici a
t = 0,4246+ε, soluzione asintotica nei piani x–t e p–u)
Ragionamento del docente:
1. **Discontinuità in x = 0**: calcola a + δu a sinistra (1,522) e a destra (1 + 0,2·2,61 = 1,522): stesso invariante,
   a scende e u sale ⇒ **espansione centrata transonica** (riconoscimento del tipo d'onda dai dati tramite invarianti).
2. **Tratto lineare 0<x<1**: u decrescente, a crescente, a + δu costante ⇒ **compressione centrata** sulla famiglia
   u − a (profilo lineare ⇒ tutte le caratteristiche convergono in un punto).
3. **t_bu = 0,4246**; a t_bu + ε dal punto di coalescenza nasce urto + contatto + espansione debole di raccordo
   (problema di Riemann al centro della compressione, cfr. lez. 20).
4. (lez. 48) Asintoto: l'interazione espansione–compressione "annulla la struttura d'urto locale, lasciando
   un'espansione residua". [Formulazione della trascrizione poco chiara: da verificare con il calcolo.]
- **Discordanza dei dati (non risolta)**: la trascrizione dà a = 1 + 0,3711x (a(1) = 1,3711), coerente con a + δu = 1,522
  costante ma che darebbe t_bu = 1/(1,8557+0,3711) = 0,449; il catalogo dà a = 1 + 0,4995x (a(1) = 1,4995), che
  riproduce t_bu = 0,4246 ma **non** mantiene a + δu costante (1,650 a x=1). Uno dei due testi è trascritto male, o il
  tratto non è un'onda semplice esatta. Il catalogo data la prova 14/09/2021 come B65 (coerente con la trascrizione).

## Getto sottoespanso con ago sull'asse (M_e = 1, p_e/p_a = 9) — non trovato nel catalogo
- Senza ago: espansione all'uscita, riflessione come compressione sul contorno isobaro ⇒ **barrel shock** + **disco di
  Mach**. Con ago: separazione dello strato limite sull'ago (gradiente avverso) ⇒ **bolla di ricircolo** subsonica ⇒
  il disco di Mach sparisce, sostituito da urti obliqui attorno alla bolla; cambia la pressione sull'ago e l'acustica.
