# 41 — Lezione 42, 09/12/2021 (file 04 = 08): tipo dell'equazione TSD; soluzione in gola ("Meyer-Sänger" = Sauer–Meyer)

## Detto esplicitamente
- TSD: u > 0 iperbolico (caratteristiche dη/dξ = ±1/√((γ+1)u)), u < 0 ellittico, u = 0 parabolico.
- Urti transonici deboli: W = −(γ+1)ū; Δs ~ [u]³ ⇒ la TSD (irrotazionale) descrive bene anche gli urti deboli.
- **Gola di un ugello simmetrico** (altezza 2H, parete approssimata col cerchio osculatore R ≫ H):
  sviluppo φ = f0 + f2 η² + f4 η⁴ (simmetria: u pari, v dispari, v(ξ,0) = 0); bilanci per potenze di η;
  ipotesi u(ξ,0) = αξ lineare con origine nel punto sonico dell'asse; α = 1/√((γ+1)RH).
  u = αξ + ((γ+1)/2)α²η², v = (γ+1)α²ξη (+ termine cubico, vedi sotto).
- Linea sonica ξ = −((γ+1)/2)αη²: **parabola**, a parete più a monte (nel convergente), sull'asse in ξ = 0.
- Conclusioni "fisiche" esplicite: la transizione sonica **non avviene su una sezione piana**; **la gola geometrica
  non coincide con la linea sonica**; isobare/isomach curve.

## Differenze / errori rispetto al formulario (§19, Sauer–Meyer)
- Nome: la trascrizione dice "Meyer-Sänger": quasi certamente storpiatura di **Sauer–Meyer** (nome del formulario e
  del summary).
- La trascrizione conclude f4 = 0 derivando come costante il prodotto f0'f2' = αξ·(γ+1)α²/2: è un errore algebrico
  (della lavagna o della trascrizione): d/dξ[f0'f2'] = (γ+1)α³/2 ≠ 0 ⇒ compare il termine cubico in v. Il formulario
  ha v = (γ+1)d²ξη + ((γ+1)d³/6)η³ e la linea di gola ξ_gola = −(d/6)η² (la linea v = 0 **non** è ξ = 0, è curva).
  Con la scalatura della trascrizione (γ+1 non assorbito) il cubico sarebbe (γ+1)²α³η³/6: differenza di
  normalizzazione. **Seguire il formulario.** (Nota: `gasdyn.sauer_meyer` ha già le due versioni 'formulario' e
  'favini'.)
- d = 1/√((γ+1)R) nel formulario (R adimensionale) ⇔ α = 1/√((γ+1)RH) dimensionale: coerenti.

## Catalogo
- **B67-E01** (flussi transonici, Sauer–Meyer: variante usata in soluzione), B08-E01, B61-E02.
