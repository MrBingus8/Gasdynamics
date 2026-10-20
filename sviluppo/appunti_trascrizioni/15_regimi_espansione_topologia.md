# 15 — Lezione 15, 20/10/2021 (file 46): limiti dell'espansione, regimi, topologia a zone, riflessione su parete

## Detto esplicitamente
- Pistone che arretra: a_p = a0 − δ|V_p| ≥ 0 ⇒ |V_p,max| = a0/δ (vuoto, "cavitazione del gas"); a quel limite le tre
  caratteristiche collassano sulla traiettoria del pistone.
- **Confronto stazionario/instazionario** (che lui ama): u_max staz = a0/√δ < u_max instaz = a0/δ.
- **Regimi del ventaglio** (coda dx/dt = a0 − (1+δ)|V_p|): subsonico (|V_p| < a0/(1+δ)), **limite sonico** (coda
  verticale, zona uniforme a M=1), **transonico** (coda inclinata a sinistra, il ventaglio scavalca l'asse t, valle
  supersonica). "Nei flussi non stazionari le condizioni soniche possono essere raggiunte in un condotto a sezione
  costante grazie alla sola propagazione ondosa" (contro il riflesso "serve una gola").
- **Teorema di Courant–Friedrichs**: confine esteso tra onda semplice e zona uniforme = caratteristica rettilinea.
- Contorno aperto: uscita ⇒ p statica ambiente; **ingresso ⇒ p0, a0** con a² + δu² = a0² al contorno.
  Parete chiusa: u=0.
- Riflessione di un'espansione su parete: espansione riflessa (stessa natura).
- **Topologia del piano x–t** (regole che usa sempre):
  - zone a 0 famiglie (uniformi), a 1 famiglia (onde semplici), a 2 famiglie (interazione, regione "romboidale");
  - confine 0↔1 = **semiretta**; 1↔2 = **curva**; 0↔2 = **un punto**.
- Calcolo dello stato finale: si attraversa ciascun ventaglio **lungo una caratteristica dell'altra famiglia**
  (A→B lungo u−a con R1; B→C lungo u+a con R2) + condizione al contorno u_C = 0 ⇒ a_C = a0 − 2δ|V_p|
  ("riduzione doppia").

## Catalogo
- Pistone + parete di fondo / riflessione dell'espansione: B10-E01, B15-E02, B20-E01, B33-E01, B60-E01, B69-E02.
