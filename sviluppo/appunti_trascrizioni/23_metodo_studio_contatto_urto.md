# 23 — Lezione 24, 08/11/2021 (file 14): come studiare; esercizio d'esame 16/03/2018 (contatto contro urto stazionario)

## Detto esplicitamente (metodo e aspettative)
- "In gasdinamica **tutto ciò che viene trattato è importante**": nessuna graduatoria. Bisogna governare tutte le
  variabili insieme: ρ, u, v, w, p e le derivate M, s, T, h0, p0 — "se trascuriamo una sola di queste variabili, la
  comprensione del fenomeno collassa".
- La soluzione generale è numerica, ma per non farne "una scatola nera" bisogna scomporre in onde semplici, urti,
  contatti, problemi di Riemann.
- "Le **prove scritte** come filo di Arianna": esercizi = casi-studio per verificare la comprensione della struttura
  fisica, non formule meccaniche.
- **Errore tipico denunciato**: "Mach in gola = 1 ⇒ portata costante?" — No: Fliegner ṁ/A* ∝ p0/√T0; se il serbatoio
  si svuota la gola resta sonica ma la portata varia.

## Esempio svolto — esame 16/03/2018 = **B56-E01**
Urto normale stazionario in condotto; A: p=1, ρ=1,4, u=2; nella metà inferiore arriva un contatto con B: p=1, u=2,
ρ=0,082 (bolla calda).
1. Stato C a valle dell'urto indisturbato: a_A = 1, M_A = 2 ⇒ p_C = 4,5, ρ_C ≈ 3,733, u_C ≈ 0,75, M_C ≈ 0,577.
2. All'istante dell'impatto: **problema di Riemann** 1 ≡ B (sinistra), 4 ≡ C (destra).
3. Osservazione fisica chiave: a_B ≈ 4,13 ⇒ M_B ≈ 0,48 — la corrente calda è **subsonica** con la stessa u: l'urto
   stazionario non può più sussistere.
4. **Validità del modello 1D** discussa esplicitamente: alla parete inferiore il moto è 1D; sulla mezzeria c'è un
   contatto di scorrimento ⇒ per t>0 lo squilibrio di pressione tra metà superiore e inferiore genera onde oblique e
   **punti tripli** (2D). 1D vale per t = 0⁺.
5. **Limite del piano p–u**: A e B sono lo stesso punto in p–u ma hanno impedenze diverse (ρa = 1,4 vs 0,339) ⇒ curve
   uscenti con pendenze diverse; non si legge a vista, serve il solutore iterativo. (Il catalogo: soluzione SCR.)
6. Un urto può **risalire una corrente supersonica** nel riferimento assoluto (solo il Mach relativo conta);
   un'espansione **non può mai** risalire una corrente supersonica.
- Seconda parte (B56-E01 domanda 2): urto obliquo e punto triplo sul piano p–θ — non svolta in questa lezione.
