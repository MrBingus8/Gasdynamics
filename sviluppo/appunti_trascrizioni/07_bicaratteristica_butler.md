# 07 — Lezione 8 parte 1, 07/10/2021 (file 34): compatibilità acustica, teorema di Butler

## Detto esplicitamente
- AFS "L'arte della modellistica numerica in gasdinamica" (1 CFU): onda lineare scalare → Burgers/Riemann
  elementare → Eulero 1D instazionario; "la simulazione numerica è una vera e propria arte".
- **Distinzione da non sbagliare**: "il regime supersonico instazionario è profondamente diverso dal regime
  supersonico stazionario". Stazionario: caratteristiche fisse nello spazio, h0 e p0 costanti lungo linea di
  corrente (isentropico). Instazionario: bicaratteristiche nello spazio-tempo; **h0 e p0 non si conservano**
  (∂p/∂t compie lavoro sulle particelle). → Errore tipico: usare grandezze totali costanti in problemi instazionari.
- Segno nell'autovettore acustico dipende dall'orientamento di ν (convenzione).

## Teoria
- Derivata lungo la bicaratteristica D_a/Dt = ∂/∂t + (V + aν)·∇; resta il termine non direzionale
  ρa²{∇·V − ν·[(ν·∇)V]} (divergenza laterale).
- **Teorema di Butler**: sul cono al più 3 direzioni bicaratteristiche indipendenti; 4 compatibilità acustiche
  sono dipendenti ⇒ esiste combinazione che elimina i termini non direzionali. Non presente nel summary
  (che passa direttamente al caso 1D) — **solo nelle trascrizioni**.

## Esempi d'esame
- Nessuno.
