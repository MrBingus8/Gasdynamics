# 40 — Lezione 41, 07/12/2021 (file 05 = 09): retroazione aeroacustica/termoacustica; equazione TSD

## Detto esplicitamente
- L'oscillazione dell'urto transonico richiede la **retroazione** (scia o sorgente di calore → acustica → urto →
  sorgente): con corpo affusolato senza scia o con lamina divisoria che stabilizza la scia, l'urto resterebbe
  stazionario. Termoacustica nei motori: risonanza con le frequenze della camera; rimedio = rompere le simmetrie
  (iniettori, piastre).
- **Fallimento della linearizzazione** per M∞ → 1: Prandtl–Glauert (subsonico) e Ackeret (supersonico) hanno
  coefficienti (1 − M²) → 0 ⇒ C_p → ∞; il termine non lineare in u_x è dello stesso ordine.
- Derivazione TSD: Φ = V0 x + εLV0 φ(ξ,η), ξ = x/L, η = λy/L con λ ≪ 1 (vicino a M=1 le caratteristiche sono quasi
  normali alla corrente: perturbazioni a grande distanza trasversale); BC al corpo y = S f(x) ⇒ ελ/S = 1;
  bilanciamento λ² = M0²(γ+1)ε ⇒ ε = S^{2/3}/[M0²(γ+1)]^{1/3}; **parametro di similitudine transonica**
  K = (1 − M0²)/[S(γ+1)]^{2/3}; forma ((1−M0²)/K' − φ_ξ)φ_ξξ + φ_ηη = 0; a M0 = 1: −φ_ξ φ_ξξ + φ_ηη = 0.
- Forma conservativa in (u, v) e urti deboli transonici con giacitura W = dη/dξ = −[v]/[u] = −((γ+1)/2)(u1 + u2).

## Differenze con il formulario (§19)
- Il formulario scrive [k − (γ+1)φ_ξ]φ_ξξ + φ_ηη = 0 (γ+1 esplicito, k = 0 al sonico; sistema (γ+1)u u_ξ − v_η = 0).
  La trascrizione assorbe (γ+1) nella scalatura (K, ε) ma poi lo fa ricomparire nella forma conservativa
  (−(γ+1)/2 u²): la scalatura della trascrizione non è internamente coerente (probabile confusione di trascrizione).
  **Seguire il formulario**; la sostanza (struttura non lineare, cambio di tipo con il segno di u) è la stessa.
- Il parametro di similitudine K e la formula della giacitura degli urti transonici: solo nelle trascrizioni.

## Catalogo
- B08-E01, B61-E02, B67-E01 (flussi transonici).
