# 32 — Lezione 33, 23/11/2021 (file 20 = 21): piano p–θ ("curva del cuore"), σ–θ, isentropica vs urto obliquo

## Detto esplicitamente
- Urto debole ⇒ valle iperbolica, "l'analisi locale dello spigolo è autoconsistente"; forte ⇒ ellittico, serve tutta la
  geometria a valle.
- Polare simmetrica; diedro simmetrico a incidenza nulla: urto superiore (famiglia θ+α) e inferiore (θ−α),
  **regioni disaccoppiate** (una lamina sulla corda non cambia nulla). A incidenza: estradosso urto più forte (se
  θ < θ_max), intradosso **ventaglio di Prandtl–Meyer** dallo spigolo di prua.
- Il piano odografo non distingue i contatti stazionari (salto di v_t a parità di p) ⇒ si usa il **piano p–θ**
  ("polare della pressione", "**curva del cuore**"), analogo del p–u.
  - parte da (0, p1) (urto nullo) e torna sull'asse in (0, p_normale) (urto normale); debole = p minore, forte = p
    maggiore; tangenza verticale a θ_max, θ_sonico a 1°–2°;
  - ramo debole quasi rettilineo e inclinato; ramo forte molto piatto (p vicino all'urto normale);
  - **Avvertenze grafiche esplicite**: con scala logaritmica la concavità del ramo inferiore si deforma; **la curva non
    deve mai partire dall'origine** (p1 > 0, al più p/p1 = 1).
- Diagramma σ–θ al variare di M1: campane tra α1 e 90°, vertice (θ_max, σ); M1→1 collasso in (0, 90°); M1→∞:
  45,58° / 67,79°. Debole = ramo con σ minore.
- Chiarimento: variazione di sola T a pari p ⇒ nessuna onda acustica, è un **contatto**; le onde di pressione nascono
  solo da gradienti/salti di p.

## Esempio numerico svolto (confronto modelli) — M1 = 2, θ = 20°, aria
- **Isentropica** (parete concava): ν(2) = 26,38° ⇒ ν2 = 6,38° ⇒ M2 ≈ 1,31; p0 costante ⇒ p2/p1 ≈ 0,356/0,1278 ≈ 2,79.
- **Urto obliquo debole** (spigolo): θ_max(2) = 22,97° > 20° ⇒ attaccato; σ = 53,42°, M1n ≈ 1,61 ⇒ M2n ≈ 0,66,
  p2/p1 ≈ 2,80; M2 = M2n / sin(σ − θ) ≈ 1,32.
- Conclusione: valori quasi uguali perché Δs ∝ (M1n² − 1)³; "finché il Mach a monte rimane inferiore a circa 3" le due
  previsioni sono sovrapponibili. [Affermazione del docente: soglia orientativa.]
- Uso delle **tabelle** (ν, isentropiche, urto obliquo/normale) con interpolazione: è il modo di lavorare in aula.

## Catalogo
- Confronto isentropico/urto, rampe e cunei a incidenza: B06-E02, B10-E03, B23-E01, B27-E03, B28-E01, B34-E03,
  B45-E03, B46-E03 (aerodinamica supersonica).
