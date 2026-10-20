# 42 — Lezione 43, 13/12/2021 (file 03 = 07): esame — getto sottoespanso che impatta su parete normale

## Esercizio = **B64-E04** (15/06/2021): getto da ugello conico, M_e = 3, p_e/p_amb = 1,36, p0/p_amb = 50; schlieren con
richiami 1–10 e mappe p_w/p_e. Affine: B42-E03 (sequenza schlieren getto impattante), B63-E02.

## Detto esplicitamente
- **Analogia 2D/3D**: la topologia delle discontinuità in un campo assialsimmetrico è la stessa del 2D; differisce
  per il termine 1/r (attenuazione più rapida, come 1/A dA/dx nel quasi-1D). Le tabelle 2D si usano localmente;
  il disco di Mach (urto normale) segue R–H esattamente anche in 3D.
- **Descrizione per numeri, nell'ordine dell'immagine** (lessico del docente):
  1 bocca d'uscita; 2 regione del getto sottoespanso (espansione centrata allo spigolo) a monte dell'urto;
  3 onde deboli dalla parete interna dell'ugello; 4 **disco di Mach**; 5 urto riflesso (conico) dal **punto triplo**;
  6 **shear layer** dal punto triplo (separa corrente subsonica dietro il disco da quella supersonica dietro i due
  obliqui); 7 regione subsonica dietro il disco fino al ristagno a parete; 8–10 **lama radiale** lungo la parete
  (deviazione di 90°, accelerazione transonica attraverso una gola efficace, riespansione, urto di ricompressione
  radiale = "ginocchio" al n. 10).
- **Profilo di pressione a parete**:
  - picco sull'asse (p_w/p_e ≈ 6): urto normale (disco di Mach) M2 → M7 + compressione isentropica subsonica fino al
    ristagno: p0/p_e = (p7/p2)(1 + δM7²)^(γ/(γ−1)) [con p2 riferita a p_e tramite l'espansione 1→2];
  - allontanandosi, la lama accelera, la pressione scende **sotto** l'ambiente (regime supersonico locale);
  - urto di ricompressione radiale ⇒ andamento **non monotono** prima di tornare all'ambiente.
- Conclusione metodologica: teorie 2D locali (R–H, Prandtl–Meyer, caratteristiche) combinate prevedono bene i
  carichi su strutture 3D, con margine adeguato al progetto.

## Ambiguità
- Trascrizione: "supera la pressione ambiente (P_W < P_E)": nel testo d'esame p_e è la pressione di **uscita
  dell'ugello** (p_e/p_amb = 1,36), non l'ambiente. Scendere sotto l'ambiente significa p_w/p_e < 1/1,36 ≈ 0,74.
  Attenzione a non confondere p_e (exit) con p_amb nelle normalizzazioni.
