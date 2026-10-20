# 26 — Lezione 27, 11/11/2021 (file 26): urto vs scorrimento nella realtà; interferogramma; quasi-1D instazionario

## Detto esplicitamente
- "Bisogna avere una **grande fiducia sulle equazioni di bilancio**" ("un atto di fede ampiamente ripagato").
- In Eulero urto e contatto sono entrambi discontinuità; nella realtà l'urto ha spessore di pochi cammini liberi
  medi, lo strato di scorrimento si allarga (instabilità, vortici, von Kármán; il **potential core** di un getto
  scompare entro ~10 D).
- Spiegazione con la teoria delle onde: nell'urto le caratteristiche convergono ed entrano ⇒ meccanismo di
  **chiusura** che bilancia la diffusione; nel contatto (linearmente degenere) solo apertura. Stesso fenomeno visto in
  Godunov (contatto che si allarga).

## Interpretazione di immagini — interferogramma dell'avviamento di un ugello (Saitō & Takayama)
Ordine e terminologia del docente:
1. **Primary shock wave**: urto d'avviamento della famiglia u+a, "veloce", nell'aria in quiete a valle; separa 4
   (quiete) da 3 (gas accelerato); linea netta e sottile.
2. **Interface** (contatto): tra gas inizialmente a valle (urtato) e gas proveniente da monte della gola espanso nel
   divergente; spessore visibile, sfumature di grigio = diffusione.
3. **Secondary shock wave**: famiglia opposta u−a, urto "lento" di **ricompressione** nel divergente: il flusso
   supersonico sovraespanso non sente la contropressione finché non vi "sbatte contro"; serve per adattare la
   pressione a quella imposta dal contatto.
- "Questo fenomeno non può essere spiegato con un modello 1D a sezione costante: richiede la variazione d'area".
Collegamento al catalogo: avviamento ugello (S "Nozzle start up"), immagini B55-E03, B57-E02, B58-E01, B59-E03,
B35-E02.

## Quasi-1D instazionario
- Sorgente B = (ρuA'/A, 0, ρa²uA'/A); autovalori e autovettori **identici** al caso a sezione costante (B non salta).
- Compatibilità: dp ± ρa du + ρa²(u/A)(dA/dx)dt = 0 lungo u±a; omoentropico: d(a ± δu) + δ(ua/A)(dA/dx)dt = 0;
  entropia non vede la geometria.
- Conseguenze (esplicite): (1) intensità non costante: in un **divergente** con u>0 l'onda si **attenua**, in un
  **convergente** si **intensifica**; (2) **non esistono onde semplici** in condotti a sezione variabile: entrambe le
  famiglie sono sempre attive (riflessioni continue verso monte se subsonico); (3) **una bolla termica che attraversa
  un ugello genera onde acustiche** — "i motori suonano" (instabilità termoacustiche, feedback con la combustione).

## Catalogo
- Onde in condotti a sezione variabile: B72-E02 (urto mobile e variazione di sezione), B08-E02, B48-E03, B51-E01.
