# 08 — Lezione 8 parte 2, 07/10/2021 (file 30): Butler con 4 direzioni, diagonalizzazione per assi

Teoria pura, non presente nel summary in questa forma.

- Termine **direzionale** = termine spaziale associato a una velocità di propagazione (modello: Ds/Dt);
  **non direzionale** = es. ρa²∇·V (isotropo).
- Butler "discretizzato" su 4 versori ±x, ±y: sommando le coppie di compatibilità acustiche si eliminano i termini
  incrociati e ∂p/∂t risulta somma di sole onde piane lungo x (u±a) e lungo y (v±a).
- Forma canonica: U_t + R_x Λ_x L_x U_x + R_y Λ_y L_y U_y = 0, Λ_x = diag(u−a,u,u,u+a); A_x e A_y non commutano ma
  sono singolarmente diagonalizzabili ⇒ "sovrapposizione di onde piane dirette lungo x e lungo y".
- Uso dichiarato: base per numerica a cattura d'onda e per **condizioni al contorno** (scrivere ∂U/∂t come
  somma di contributi d'onda → scegliere quali onde entrano dal contorno). → collegare a lezioni su BC/Godunov.
- Nessun esempio d'esame.
