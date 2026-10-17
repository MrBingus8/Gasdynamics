# 16 — Esercitazione 21/10/2021 (file 42): pistone con legge parabolica (espansione)

## Atteggiamento didattico (esplicito)
- "Lo scopo di fare esercitazione è far uscire soprattutto le cose sbagliate in modo tale che poi le buttiamo via";
  "non sbagliare non serve".

## Esempio svolto — pistone x_p = b t², b < 0 (espansione), gas in quiete (u=0, a=a0)
Corrisponde alle tracce con legge parabolica: **B53-E14** (x_P = −t² poi fermo), B20-E01, B33-E01, B49-E01,
B69-E02 (casi di compressione o con tratto successivo a velocità costante).
Ragionamento, nell'ordine del docente:
1. Classificare: **problema misto** (IC + BC). BC sul pistone: u(x_p,t) = ẋ_p (unico dato cinematico; il pistone
   "subisce la pressione"). IC: quiete.
2. **Verificare il limite di cavitazione**: a_p ≥ 0 ⇒ |ẋ_p| ≤ a0/δ ⇒ limitare i tempi.
3. Osservare che ẋ_p(0) = 0: **nessuna singolarità** in (0,0) (≠ caso impulsivo).
4. Verificare l'omoentropia (R1 = R2 = a0 in quiete) — digressione su dati sperimentali: calcolare p/ρ^γ, mediana e
   dispersione; scarti di qualche % ⇒ omoentropico accettabile; 10–15% ⇒ no.
5. Topologia x–t: fronte anteriore x = a0 t (testa u+a); a destra zona a 0 famiglie. **Nessun fronte posteriore**
   (velocità del pistone sempre diversa ⇒ nessuna zona uniforme a ridosso): tutto il campo tra pistone e testa è
   un'onda semplice di espansione.
6. R1 = a0 uniforme ⇒ a = a0 + δu; caratteristiche u+a semirette dal pistone.
7. Per un punto (x, t̄): caratteristica a ritroso fino al pistone al **tempo di emissione** t_p:
   x − b t_p² = (a0 + (γ+1) b t_p)(t̄ − t_p) ⇒ γ b t_p² − (a0 + (γ+1) b t̄) t_p + (a0 t̄ − x) = 0; si sceglie la radice
   con 0 ≤ t_p ≤ t̄. Poi u = 2b t_p, a = a0 + 2δ b t_p.
8. Domanda tipica: "c'è discontinuità di gradiente sui fronti?" Sul fronte anteriore (∂u/∂x)_FA =
   2b/((γ+1)b t̄ − a0) ≠ 0 mentre a monte è 0 ⇒ **discontinuità debole** sul fronte anteriore.
9. Estensione: b > 0 ⇒ compressione, caratteristiche convergenti ⇒ tempo e ascissa di formazione dell'urto.

## Inferenza
- Schema ricorrente: (classificazione IC/BC) → (limiti fisici: cavitazione/vuoto) → (singolarità sì/no) →
  (omoentropia) → (topologia x–t: zone 0/1/2 famiglie, fronti) → (invariante uniforme) → (calcolo lungo le
  semirette) → (discussione delle discontinuità deboli).
