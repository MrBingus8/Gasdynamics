# 24 — Lezione 25, 09/11/2021 (file 11): shock fitting vs shock capturing, schema di Godunov

## Detto esplicitamente
- Due filosofie: **shock fitting** (caratteristiche nel continuo + R–H su discontinuità tracciate: precisissimo, ma
  topologia ingestibile in 2D/3D con punti tripli, interazioni) e **shock capturing** (forma integrale conservativa
  su griglia fissa; urti "catturati" su poche celle). Oggi domina il capturing.
- **Godunov (1959)**: volumi finiti + Riemann esatto locale.
  - U = (ρ, ρu, ρE), F = (ρu, ρu²+p, ρuH); media di cella = valore al centro a meno di O(Δx²).
  - Avanzamento esatto: Ū^{n+1} = Ū^n − (1/Δx)∫[F_{j+1/2} − F_{j−1/2}]dt (nessuna approssimazione finché…).
  - Interfaccia j+1/2: RP locale con stato 1 = cella j (sinistra), 4 = cella j+1 (destra); stato all'interfaccia =
    soluzione RP **lungo x/t = 0**; casi: in 1, in 2 o 3, in 4, **dentro un ventaglio transonico (u = a)**, su urto/contatto
    stazionario. F^Godunov = F(U^RP).
  - **CFL**: le onde da interfacce adiacenti non devono raggiungere l'interfaccia in Δt ⇒ Δt = α Δx/λ_max,
    λ_max = max(|u|+a), 0 < α ≤ 1.
- **Homework** (tipico da esame "calcolare a mano il primo passo"): condotto [0,1], u=1, p=1 uniformi, contatto puro
  in mezzeria ρ_L = 1,4 (a_L = 1), ρ_R = 0,35 (a_R = 2), primo passo di Godunov con α = 1.
  → il RP all'interfaccia è un contatto puro che si muove a u=1: l'asse t cade nello stato di sinistra ⇒ flusso
  F(U_L) ⇒ trasporto upwind. [Inferenza mia.]

## Differenze con le fonti
- Il summary ha Godunov (stessa derivazione, Δt = αΔx/λ_max, senza specificare α ≤ 1; aggiunge che Godunov produce
  entropia e le ghost cells per i contorni). **Il formulario non ha una sezione di numerica.**
- Shock fitting vs capturing come confronto esplicito: solo nelle trascrizioni.

## Catalogo
- Godunov / metodi numerici: B02-E03, B05-E02, B07-E02, B11-E02, B11-E05, B14-E02, B38-E02, B44-E03, B57-E01.
