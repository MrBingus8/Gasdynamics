# 25 — Lezione 26, 10/11/2021 (file 15): Godunov a mano su un contatto, upwind, dissipazione numerica

## Esempio svolto (homework della lez. 25) — tipologia "metodi numerici" del catalogo (B44-E03, B38-E02, B14-E02…)
Procedura completa del docente (da replicare negli esercizi "primo passo di Godunov"):
1. a_L = 1, a_R = 2; λ_max = max(|u|+a) = 3; con α = 1 ⇒ Δt/Δx = 1/3.
2. RP all'interfaccia con il solutore: Z = 1/2, u⁽⁰⁾ = 1 = u1 = u4 ⇒ onde acustiche di intensità nulla, solo il
   contatto, W_C = 1 > 0 ⇒ l'asse t cade nello **stato 1** ⇒ U^RP = Ū_j.
3. Le altre interfacce: RP banali.
4. Aggiornamento: cella j invariata; cella j+1: ρ = 0,35 − (1/3)(0,35 − 1,4) = 0,70 ⇒ profilo 1,4 | 1,4 | 0,70 | 0,35:
   **diffusione numerica** (il gradino si sdoppia).
5. Verifica: q.d.m. ed energia danno u = 1, p = 1 esatti (nessuna oscillazione spuria).

## Detto esplicitamente
- Con u>0 Godunov = **upwind al 1° ordine**; "Godunov non è altro che la generalizzazione non lineare della teoria
  delle caratteristiche upwind" (legame con Cauchy–Kowalevski).
- Ū_j^{n+1} = (1 − c)Ū_j + c Ū_{j−1}, c = uΔt/Δx: media convessa; c = 1 trasla esattamente; c < 1 equazione modificata
  con ν_num = (uΔx/2)(1 − c).
- **Contatto = linearmente degenere**: nessuna caratteristica vi entra ⇒ si allarga indefinitamente.
  **Urto = genuinamente non lineare**: le caratteristiche della sua famiglia entrano da entrambi i lati ⇒ resta
  compatto su **2–3 celle**. (Domanda d'esame tipica nell'interpretazione di risultati numerici.)
- Secondo homework: tubo d'urto con Godunov, p_L=1, ρ_L=1,4, p_R=0,1, ρ_R=0,35, quiete.

## Ambiguità
- Condizione di ingresso delle caratteristiche scritta "u1 + a1 > W > u2 + a2": con 1 = monte di un urto u+a la
  condizione corretta è u1 + a1 < W < u2 + a2. Probabilmente qui 1/2 indicano sinistra/destra (valle/monte) o è un
  errore di trascrizione. Da non riportare come formula.
