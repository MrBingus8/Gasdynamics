# 11 — Lezione 11, 13/10/2021 (file 44): upwind, MOC al 1° ordine, IBVP e regole per le condizioni al contorno

## Detto esplicitamente
- Superfici caratteristiche: per i matematici "un blocco analitico; per noi la via" per risolvere i flussi comprimibili.
- Schema: per λ>0 derivata a sinistra (monte), per λ<0 a destra (valle). Motivo dello splitting anche senza
  discontinuità: **stabilità** (CFL).
- **MOC al 1° ordine**: Δp ± (ρa)^n Δu = 0, Δx = (u±a)^n Δt con impedenza e pendenze congelate a t_n;
  Δt troppo grande ⇒ "partire per la tangente". Oggi il MOC geometrico "non viene quasi più impiegato" nei codici.
- Tubo d'urto (Vieille): le equazioni valgono "incluso rigorosamente il bordo".
- **Condizioni al contorno (regole del docente)**:
  - n. di condizioni = n. di caratteristiche **entranti** dall'esterno; varia nel tempo. In 1D-con-v (4 onde):
    ingresso supersonico 4, ingresso subsonico 3, uscita supersonica 0 (uscita subsonica: 1).
  - Divieti: (1) niente velocità infinite (caratteristiche mai orizzontali); (2) niente sovraspecificazione;
    (3) "**Divieto assoluto di prescrivere variabili d'onda uscenti**" (quelle che vengono dall'interno);
    (4) efflusso verso ambiente: si impone p statica ambiente; **immissione dall'ambiente: si impone la pressione
    totale p0 (o h0)**, non la statica.
- Equazione al contorno: ∂U/∂t + RΛ⁺L(U_x)⁺ = −RΛ⁻L(U_x)⁻ (noto dall'interno a destra, incognito a sinistra).

## Collegamento al catalogo
- Tipologia "Condizioni al contorno" (B47-E01, B53-E02/E04/E09, B67-E02, B68-E02, B70-E02) e "Metodo delle
  caratteristiche e condizioni al contorno" (B34-E02): contare le caratteristiche entranti, imporre solo quelle.
