# 28 — Lezione 29, 16/11/2021 (file 27): lettura dettagliata dell'interferogramma d'avviamento (Saitō–Takayama)

## Detto esplicitamente — concetti
- In un **divergente supersonico** nulla torna indietro (nessun choking: A > A*); in un **convergente** chiedersi sempre
  se A_min < A* ⇒ se sì, urto che risale per riorganizzare il campo in subsonico.
- **Urti instazionari**: deve essere supersonico il Mach **relativo** (u1 − W)/a1; un urto può risalire una corrente
  supersonica (es. urto riflesso da parete). "Attenzione a non confondere il regime della corrente con la cinematica
  degli urti instazionari."
- Urto obliquo: v_t conservata, v_n diminuisce ⇒ la velocità ruota verso il fronte; a valle di un obliquo il flusso
  può restare supersonico.

## Interpretazione di immagini — come la fa (ordine e lessico)
1. **Primary shock wave**: linea scura, netta, sottilissima (caratteristiche convergenti la mantengono sottile).
2. **Interface**: transizione sfumata, spessore visibile (diffusione: caratteristiche non convergenti).
3. **Frange bianche/nere**: l'interferometro vede variazioni di indice di rifrazione ∝ gradienti di ρ; le frange sono
   **curve di livello quantizzate dalla soglia dello strumento**, non bande fisiche. Analogia con le medie di cella
   di Godunov; interpolare (spline/Fourier) non aggiunge informazione (N punti ⇒ N/2 frequenze).
4. **Secondary shock**: non è un urto normale ma una struttura 2D: due obliqui dalla gola che convergono
   ("clessidra"), raccordati da un tratto normale ⇒ **riflessione di Mach** con **punto triplo** (incidente,
   riflesso, disco di Mach, linea di scorrimento); qui un doppio punto triplo simmetrico.
5. **Grafici di pressione assiale** (esperimento vs quasi-1D, tre istanti):
   - salto dell'urto primario che **diminuisce** nel tempo (urto in un divergente si indebolisce);
   - tra secondario e primario p cresce verso valle in un divergente ⇒ dalla legge delle aree
     dp/(ρu²) = dA/A /(1−M²) ⇒ **regione subsonica** (deduzione del regime dal segno del gradiente);
   - il **contatto non si vede nella pressione** ([p] = 0), solo nella densità;
   - distanza primario–secondario cresce: primario "veloce" (u+a), secondario "lento" (u−a, "retrogrado" nella
     letteratura russa) che contrasta la corrente supersonica;
   - terzo salto vicino alla gola nella simulazione quasi-1D: **non è un errore**, è la proiezione mediata sull'asse
     della struttura 2D di urti obliqui a parete.

## Inferenza (metodo per le domande "commentare un'immagine/grafico")
- Identificare ogni struttura per tipo d'onda e famiglia; giustificare l'aspetto visivo con la fisica (spessore
  ⇔ convergenza delle caratteristiche); dedurre il regime (sub/super) dai gradienti con la legge delle aree;
  distinguere ciò che è fisico da ciò che è artefatto dello strumento o del modello (frange, salto mediato).
- Relazione di uno studente sull'avviamento con quasi-1D: stato a monte fortemente supersonico ⇒ tutte le onde
  verso valle ⇒ trattazione a mano semplificata.
