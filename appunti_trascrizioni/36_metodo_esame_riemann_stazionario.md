# 36 — Lezione 37, 30/11/2021 (file 17): metodo vs "formulette"; casi d'esame; Riemann stazionario supersonico

## Detto esplicitamente — atteggiamento richiesto all'esame
- Superare "la logica scolastica dell'applicazione meccanica di **formulette**": nessun esercizio è standard, si
  costruisce una strategia (teoria + fenomenologia + calcolo).
- "Il tempo di risoluzione dipende essenzialmente dalla **comprensione qualitativa del campo d'onda nello
  spazio-tempo**"; leggere bene il testo fa trovare **scorciatoie esatte** (es. condizione critica).
- Dati scarsi/incerti: l'ingegnere formula **ipotesi fisiche coerenti** (e le dichiara).

## Casi d'esame discussi
1. **Tracciato di un trasduttore nel tubo d'urto**: salto seguito subito da calo senza plateau ⇒ l'espansione riflessa
   ha già raggiunto l'urto e ne ha consumato la zona uniforme (forma tipo Friedlander); la velocità dell'urto si
   ricava dai salti di pressione, considerando oscillazioni e incertezze.
   (Catalogo: tubo d'urto/tracciati, B34-E01, B36-E03, B39-E02, B41-E03, B44-E01; B67-E02/B73 grafici.)
2. **Pistone con legge cubica x_p = c t³**: t_p(x,t) non invertibile in forma chiusa ⇒ costruire la soluzione
   **tracciando a ritroso le semirette u+a** fino al pistone, punto per punto (soluzione esatta ma per punti).
3. **Imbocco da serbatoio infinito in condotto a sezione costante**: l'ingresso è l'area minima; una corrente non può
   diventare supersonica senza divergente ⇒ se il rapporto di pressione è basso, **imbocco critico** (M = 1, p*, T*)
   ⇒ disaccoppia il calcolo dell'urto (ripresa della lez. 19, esame 24/01/2019).

## Urto piano su contatto inclinato → Riemann stazionario supersonico
- Riferimento del punto d'interazione (V_T = W/cos φ). Stati: C (monte, M_C = V_T/a_C > 1), A (valle dell'urto
  incidente obliquo, σ = 90° − φ, soluzione debole, deflesso θ_A), B (oltre il contatto: p_B = p_C, ρ diversa,
  M_B = V_T/a_B).
- **Condizioni di raccordo sulla linea di scorrimento**: **stessa p e stessa θ** (analogo di p, u nel caso 1D).
- Da A (alta p) parte un'**espansione di Prandtl–Meyer**; da B (bassa p) un **urto obliquo**; in mezzo linea di
  scorrimento. Algoritmo: **iterare su θ** fino a uguagliare le pressioni (PM sopra, polare sotto).
- **Condizione di validità esplicita**: tutte le regioni devono restare **supersoniche**; se una diventa subsonica
  il modello stazionario/autosimile cade ⇒ Eulero 2D instazionario numerico. → verifica obbligatoria.

## Catalogo
- Riemann stazionario / interazione di correnti: **B33-E02**, **B74-E01**, B56-E01 (domanda 2), getti a contatto.
