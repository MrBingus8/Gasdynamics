# 44 — Lezioni 48–49, 22/12/2021 (file 02): modalità d'esame; svuotamento di un serbatoio; chiusura

## Modalità d'esame (detto esplicitamente)
- Ammessi: **tabelle** della gasdinamica, righello, **goniometro**, calcolatrice (⇒ misure di angoli su figure,
  interpolazione da tabelle).
- Orale di norma entro 15 giorni dallo scritto.
- **Consiglio di redazione**: scrivere prima il corpo e i risultati, **l'introduzione per ultima** (deve sintetizzare ciò
  che è stato effettivamente dimostrato); vale per relazioni e prove d'esame.

## Esempio d'esame — svuotamento di un serbatoio = **B66-E01** (grafico dello stato medio e della portata nel tempo)
1. Adimensionalizzazione rispetto ai valori iniziali (tutte le curve partono da 1).
2. **Portata adimensionale costante ≈ 0,578**: efflusso critico (Saint-Venant–Wantzel / Fliegner); adimensionalizzando
   con A* p0(t)/√(R T0(t)) la dipendenza dal tempo si elide: ṁ* = √γ (2/(γ+1))^((γ+1)/(2(γ−1))) = 0,578 (γ = 1,4).
   Riconoscere il numero nel grafico ⇒ dedurre il regime sonico.
3. Vale finché p0(t) ≥ p_a((γ+1)/2)^(γ/(γ−1)); poi "ginocchio" e portata → 0 all'equilibrio.
4. **Scelta del modello** (esplicita): per lo svuotamento **non** si usa la teoria delle caratteristiche (transitori
   d'onda) ma il **bilancio integrale a parametri concentrati**: d(ρV)/dt = −ṁ(t) con ṁ critica e trasformazione
   isentropica nel serbatoio ⇒ ODE integrabile (Runge–Kutta).
   [Differenza di scala dei tempi: onde (L/a) vs svuotamento (V/(A* a)) — inferenza mia, coerente con F §17.]

## Chiusura
- Rifrazione di un urto obliquo su un contatto tra due correnti supersoniche a impedenza diversa: continuità di p e θ
  decide la natura dell'onda riflessa e trasmessa.
