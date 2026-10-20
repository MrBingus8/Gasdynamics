# Analisi statistica delle prove di Gas Dinamica

**Data dell’analisi:** 9 settembre 2026. **Unità:** esercizio, non pagina, richiesta o svolgimento. **Scopo:** orientare lo studio verso competenze trasferibili; non ricostruire un programma ufficiale o predire gli appelli futuri.

## 1. Risultato dell’audit e piano eseguito

La sequenza è stata: lettura delle due tabelle complete → controllo delle anomalie → riconciliazione per scenario, dati e richieste → definizione del sottoinsieme principale → codifica multilabel con provenienza esplicita/inferita → conteggi e co-occorrenze → confronto dei periodi → priorità ed esercizi rappresentativi. Nessuno svolgimento è stato validato fisicamente.

**Correzione documentale importante:** il catalogo cumulativo 2017–2020 dichiara 30 esercizi, ma la tabella effettiva contiene 29 righe (B01-E01–B11-E02): 27 tracce identificabili e due sole intestazioni vuote. Il PDF originale conferma la sequenza di esercizi 2.1–5.11. Il totale dichiarato 30 e il totale dichiarato 20 svolgimenti non sono quindi utilizzati come conteggi verificati. Le righe con materiale di svolgimento sono 19, incluso il materiale parziale di 5.10.

| Passaggio | Conteggio | Criterio |
| --- | --- | --- |
| Catalogo 2003–2023 | 197 | Schede incluse ricostruzioni e frammenti. |
| Catalogo 2017–2020 | 29 righe | 27 contenuti + 2 intestazioni vuote. |
| Sovrapposizioni confermate | 13 | Identità di scenario, dati e richieste; non semplice analogia. |
| Nuovi esercizi identificabili | 14 | 27 − 13. |
| Totale unico esteso | 211 | 197 + 27 − 13. |
| Analisi principale | 201 | Originale disponibile e identificabile; esclude 8 ricostruzioni residue e 2 tracce tagliate. |
| Analisi di sensibilità | 211 | Include i 10 casi incompleti/ricostruiti, senza inventare parti mancanti. |

Le tre tracce B53-E04, B53-E05 e B53-E06 della raccolta ampia sono state completate documentalmente attraverso la raccolta Franco, rispettivamente con F/B03-E01, F/B03-E02 e F/B03-E03. La corrispondenza di B53-E05 è sostenuta anche dai dati M=7, p=1, ρ=1,4, w=9 e A₂/A₁=2. Per B53-E06, il frammento M=6 e i due angoli limite coincide con la traccia di marzo. La data precisa non è stata inventata.

Le lettere **G/** e **F/** nel rapporto distinguono, rispettivamente, gli ID originali del catalogo 2003–2023 e di quello 2017–2020. Non sono nuovi numeri d’esame. Il crosswalk completo conserva entrambi gli ID e le pagine fisiche.

### Crosswalk delle 13 sovrapposizioni

| Catalogo F | Catalogo G | Evidenza essenziale |
| --- | --- | --- |
| F/B01-E01 | G/B52-E01 | Stesso getto e condotto; M=5,5, p1=0,1149 atm, H=10 cm e L=8,692 cm. |
| F/B01-E02 | G/B52-E02 | Stessa onda u+a e quattro richieste sulla velocità del suono totale. |
| F/B01-E03 | G/B52-E03 | Stesso urto Mi=1,55 e sequenza di diffrazione sullo spigolo. |
| F/B02-E01 | G/B53-E01 | Stessa parete concava e richieste sui profili M, p, p0 e sull’urto esterno. |
| F/B02-E02 | G/B53-E02 | Stessi vincoli rho u costante e T costante; stesse tre richieste. |
| F/B02-E03 | G/B53-E03 | Stessa retropropulsione e richiesta su urti e contatti nel campo di Mach. |
| F/B03-E01 | G/B53-E04 | Stesso problema di portata transitoria al contorno con pressione costante; il catalogo ampio conserva solo il manoscritto. |
| F/B03-E02 | G/B53-E05 | Corrispondenza specifica: M1=7, p1=1, rho1=1,4, w=9 e A2/A1=2; la traccia integra il manoscritto. |
| F/B03-E03 | G/B53-E06 | Stesso M=6 e richiesta degli angoli di von Neumann e di distacco; il frammento è completato dalla traccia. |
| F/B04-E01 | G/B55-E01 | Stessa discontinuità di gradiente con primitive continue e richiesta di derivata temporale di T. |
| F/B04-E02 | G/B55-E02 | Stesso contatto inclinato di 42,83 gradi nel tubo d’urto e stessa rappresentazione p-theta. |
| F/B04-E03 | G/B55-E03 | Stessa richiesta di descrivere la visualizzazione e stimare il Mach del proiettile. |
| F/B05-E02 | G/B59-E02 | Stesso ugello a spina piano, M=3 e quattro rapporti/angoli richiesti. |

I due segnaposto vuoti F/B05-E01 e F/B05-E03 restano nel registro di audit, ma non sono esercizi classificabili. Restano inoltre otto ricostruzioni B53-E07–E14 senza originale integrale e due tracce B28-E01–E02 con margine tagliato. I casi con data incerta sono utilizzabili per le statistiche tematiche quando il contenuto è leggibile, ma non vengono attribuiti ad anni per supposizione. Le varianti, le traduzioni e le copie di pagina già individuate nei registri non generano nuove unità.

## 2. Tassonomia e regole di codifica

Le categorie sono non esclusive. **E** indica un fenomeno o uno strumento nominato nell’enunciato o nelle richieste; **S** un’etichetta supportata dai sottoargomenti del catalogo, che sono già inferenze documentali; **I** un modello considerato applicabile per ragionamento fisico. I conteggi degli strumenti sono quindi conteggi di applicabilità operativa, non di formule esplicitamente richieste. Ogni riga del dataset contiene i tag e la rispettiva provenienza. Le etichette non certificano che un modello sia sempre l’unico, o il migliore, per risolvere il problema.

Le famiglie sono: **fenomeni** (urti, onde continue, contatti, interazioni/riflessioni, formazione di urti, transonico, choking, acustica, viscosità e miscelamento); **contesti** (pistoni, serbatoi, condotti/ugelli, corpi, getti, misure); **modelli** (Eulero, caratteristiche, Riemann, RH, polari, PM, quasi-1D, bilanci integrali, potenziale, numerica, termodinamica e sorgenti); **azioni richieste** (derivare, impostare, calcolare, rappresentare, interpretare, verificare, progettare, costruire algoritmi e confrontare modelli). Un modello comune, come Eulero, non è automaticamente attribuito a tutti i problemi di gas ideale.

Le categorie “teorico”, “quantitativo”, “qualitativo”, “interpretativo” e “numerico” possono sovrapporsi. “Numerico” significa richiesta di metodo di calcolo, non semplice sostituzione numerica. “Interpretativo” richiede una visualizzazione o misura fornita. Un quesito su un urto e un’espansione è contato una volta in ciascuna categoria, ma una sola volta nel denominatore.

La codifica è riproducibile e controllata sulle eccezioni documentali; non equivale a una seconda correzione indipendente di tutti gli enunciati. Perciò le percentuali descrivono **questo protocollo di catalogazione**, non un censimento fisico infallibile. I campi originali e tutte le regole sono inclusi nel pacchetto.

## 3. Frequenze osservate

Tutte le percentuali seguenti usano **N=201 esercizi**. La somma può superare il 100% perché le categorie sono multiple. La colonna E distingue le etichette direttamente supportate dal testo da quelle del catalogo o inferite.

### Fenomeni e oggetti fisici

| Fenomeno | n | % | E | S/I |
| --- | --- | --- | --- | --- |
| Urti | 87 | 43.3% | 61 | 26 |
| Interazioni e riflessioni | 72 | 35.8% | 40 | 32 |
| Onde continue | 52 | 25.9% | 27 | 25 |
| Contatti | 23 | 11.4% | 16 | 7 |
| Getti e adattamento | 21 | 10.4% | 17 | 4 |
| Formazione di urti | 11 | 5.5% | 6 | 5 |
| Transonico e linea sonica | 11 | 5.5% | 8 | 3 |
| Strato limite e viscosita | 8 | 4.0% | 6 | 2 |
| Acustica | 5 | 2.5% | 3 | 2 |
| Choking e condizioni critiche | 5 | 2.5% | 4 | 1 |
| Miscelamento | 3 | 1.5% | 3 | 0 |

![Frequenze dei fenomeni](fenomeni.png)

L’evidenza più robusta è la ricorrenza di **urti e interazioni**, spesso insieme a onde continue e contatti. Le categorie “formazione di urti” e “choking” sono volutamente più strette: non includono ogni problema in cui un urto o un Mach sonico potrebbe comparire incidentalmente.

### Modelli e strumenti: applicabilità operativa

| Strumento | n | % | E | S/I |
| --- | --- | --- | --- | --- |
| Termodinamica e grandezze totali | 95 | 47.3% | 71 | 24 |
| Rankine-Hugoniot | 93 | 46.3% | 5 | 88 |
| Eulero e bilanci locali | 80 | 39.8% | 11 | 69 |
| Caratteristiche 1D | 76 | 37.8% | 7 | 69 |
| Prandtl-Meyer e caratteristiche 2D | 63 | 31.3% | 2 | 61 |
| Polari e odografa | 41 | 20.4% | 23 | 18 |
| Problema di Riemann | 40 | 19.9% | 10 | 30 |
| Bilanci integrali | 28 | 13.9% | 2 | 26 |
| Quasi-1D e area-Mach | 23 | 11.4% | 1 | 22 |
| Metodi numerici | 13 | 6.5% | 13 | 0 |
| Sorgenti e modelli non standard | 11 | 5.5% | 10 | 1 |
| Potenziale e modelli ridotti | 6 | 3.0% | 5 | 1 |
| Adimensionalizzazione | 4 | 2.0% | 3 | 1 |

Il dato sui modelli non va letto come “93 tracce chiedono esplicitamente Rankine–Hugoniot”: solo 5 sono classificate come menzione diretta, mentre il resto è applicabilità inferita. Analogamente, le caratteristiche 1D e il metodo di Prandtl–Meyer sono strumenti che collegano molte situazioni diverse. La centralità di bilanci ed equazioni fondamentali è superiore a quanto suggerirebbe il numero delle sole domande teoriche.

### Che cosa viene chiesto di fare

| Azione richiesta | n | % | E |
| --- | --- | --- | --- |
| Determinare stati e grandezze | 139 | 69.2% | 136 |
| Costruire diagrammi e profili | 56 | 27.9% | 56 |
| Interpretare osservazioni | 36 | 17.9% | 36 |
| Ammissibilita e verifica | 31 | 15.4% | 31 |
| Impostare modello e condizioni | 30 | 14.9% | 20 |
| Confrontare modelli e parametri | 25 | 12.4% | 22 |
| Derivare e giustificare | 25 | 12.4% | 19 |
| Costruire algoritmo numerico | 17 | 8.5% | 16 |
| Progettare e ottimizzare | 7 | 3.5% | 7 |

![Frequenze delle competenze](competenze.png)

**Nota metodologica:** “determinare stati e grandezze” include grandezze simboliche, andamenti e stati fisici; non significa che tutti questi esercizi richiedano un numero finale. “Ammissibilità e verifica” include limiti, regimi e confronti. Le richieste abbreviate nel catalogo sono integrate dal testo originale quando disponibile; le azioni non leggibili non vengono inventate.

### Tipologia dei problemi

| Tipo | n | % |
| --- | --- | --- |
| Quantitativo | 141 | 70.1% |
| Qualitativo | 120 | 59.7% |
| Teorico | 41 | 20.4% |
| Interpretativo | 36 | 17.9% |
| Numerico | 17 | 8.5% |

Le categorie teorica, quantitativa, qualitativa, interpretativa e numerica non sono alternative: un esercizio può chiedere una derivazione, un calcolo e una discussione. Il conteggio per categoria non è una distribuzione percentuale di domande d’esame.

## 4. Combinazioni e trasversalità

Per due categorie A e B, il conteggio è il numero di esercizi che possiedono entrambe. La percentuale è n(A∩B)/201. È riportata anche la quota condizionata n(A∩B)/n(A), che misura quanto spesso B accompagna A nel campione, non la probabilità di un futuro appello.

| Combinazione | n | % | B tra gli A (%) | ID di esempio |
| --- | --- | --- | --- | --- |
| Interazioni e riflessioni + Urti | 49 | 24.4% | 68.06 | G/B03-E02, G/B05-E01, G/B06-E01, G/B07-E01, G/B11-E05 … |
| Onde continue + Urti | 23 | 11.4% | 44.23 | G/B05-E01, G/B06-E01, G/B06-E02, G/B23-E02, G/B26-E01 … |
| Contatti + Urti | 18 | 9.0% | 78.26 | G/B07-E01, G/B26-E01, G/B27-E02, G/B30-E01, G/B30-E02 … |
| Contatti + Interazioni e riflessioni | 16 | 8.0% | 69.57 | G/B04-E01, G/B07-E01, G/B13-E03, G/B16-E03, G/B27-E02 … |
| Formazione di urti + Onde continue | 9 | 4.5% | 81.82 | G/B06-E01, G/B35-E03, G/B36-E01, G/B45-E02, G/B49-E01 … |
| Caratteristiche 1D + Rankine-Hugoniot | 40 | 19.9% | 52.63 | G/B03-E03, G/B05-E01, G/B06-E01, G/B07-E01, G/B13-E01 … |
| Prandtl-Meyer e caratteristiche 2D + Rankine-Hugoniot | 34 | 16.9% | 53.97 | G/B03-E02, G/B06-E02, G/B11-E04, G/B14-E03, G/B19-E01 … |
| Polari e odografa + Rankine-Hugoniot | 33 | 16.4% | 80.49 | G/B03-E02, G/B19-E01, G/B20-E02, G/B21-E01, G/B27-E03 … |
| Costruire diagrammi e profili + Determinare stati e grandezze | 42 | 20.9% | 75.0 | G/B01-E02, G/B03-E02, G/B07-E01, G/B09-E02, G/B17-E01 … |
| Derivare e giustificare + Impostare modello e condizioni | 14 | 7.0% | 56.0 | G/B02-E03, G/B10-E04, G/B13-E01, G/B16-E02, G/B31-E02 … |
| Costruire algoritmo numerico + Impostare modello e condizioni | 11 | 5.5% | 64.71 | G/B02-E03, G/B04-E02, G/B05-E02, G/B07-E02, G/B11-E02 … |
| Costruire diagrammi e profili + Interpretare osservazioni | 14 | 7.0% | 25.0 | G/B29-E01, G/B29-E02, G/B34-E03, G/B39-E02, G/B54-E03 … |

La matrice completa di tutte le coppie, con i rispettivi ID, è in `cooccurrences.csv`. Le coppie tra modelli sono inferenze di applicabilità, mentre le coppie tra richieste sono basate sulle azioni catalogate. I valori alti di co-occorrenza non dimostrano causalità né indipendenza statistica.

### Interpretazione fisica dei collegamenti

**Onde 1D → urti → interazioni.** Le caratteristiche descrivono propagazione e onde continue; quando la soluzione regolare non è più valida, servono salti e condizioni di ammissibilità. Le interazioni introducono nuovi stati e compatibilità, senza cambiare i principi di base. È il collegamento più utile per pistoni, tubi d’urto e condizioni al contorno.

**Urti obliqui → polari → compatibilità → espansioni.** Le polari e le relazioni di salto determinano gli stati possibili; la compatibilità di pressione e direzione seleziona le interazioni. Prandtl–Meyer e caratteristiche 2D completano il trattamento delle regioni continue. La geometria può cambiare senza cambiare il nucleo di ragionamento.

**Bilanci → condizioni critiche → modello ridotto.** Condotti, ugelli, serbatoi e spinta richiedono bilanci con ipotesi esplicite, relazioni termodinamiche e verifica del regime. Il passaggio da modello locale a integrale o quasi-1D va giustificato.

**Equazioni → discretizzazione → controlli.** Godunov, stabilità e conservazione non sono un blocco separato dalla fisica: dipendono dal problema di Riemann, dai salti e dalla scelta delle condizioni al contorno.

## 5. Copertura cronologica e confronto descrittivo

Solo **191/201** esercizi principali hanno un anno attribuibile. Di questi, **174** hanno una data giornaliera attestata (inclusi due con asterisco non spiegato) e **17** hanno soltanto una collocazione parziale, soprattutto i mesi 2017–2020. I 174 esercizi con giorno attestato appartengono a **75 date distinte**: questa è una verifica di copertura, non il numero totale degli appelli tenuti nel periodo. Dieci esercizi non sono assegnati a un anno certo.

| Anno | Esercizi principali | Data giornaliera | Data parziale |
| --- | --- | --- | --- |
| 2002.0 | 12 | 12 | 0 |
| 2003.0 | 9 | 9 | 0 |
| 2009.0 | 3 | 3 | 0 |
| 2010.0 | 16 | 16 | 0 |
| 2011.0 | 11 | 11 | 0 |
| 2012.0 | 5 | 5 | 0 |
| 2013.0 | 12 | 12 | 0 |
| 2014.0 | 14 | 14 | 0 |
| 2015.0 | 19 | 19 | 0 |
| 2016.0 | 17 | 17 | 0 |
| 2017.0 | 9 | 6 | 3 |
| 2018.0 | 19 | 19 | 0 |
| 2019 | 3 | 0 | 3 |
| 2020 | 11 | 0 | 11 |
| 2021.0 | 15 | 15 | 0 |
| 2022.0 | 4 | 4 | 0 |
| 2023.0 | 12 | 12 | 0 |
| unknown | 10 | 0 | 0 |

Il nome della raccolta ampia non delimita cronologicamente il contenuto: include anche il 2002. Mancano interi anni nella copertura di quella raccolta; il blocco Franco aggiunge il 2019 e 2020. La diversa quantità di materiale non permette di trattare gli anni come un campione uniforme di tutte le prove effettivamente svolte.

### Periodi confrontabili soltanto in senso descrittivo

| Periodo | N | Urti | Interazioni | Onde continue | Getti |
| --- | --- | --- | --- | --- | --- |
| 2002–2003 | 21 | 7 (33.3%) | 6 (28.6%) | 7 (33.3%) | 1 (4.8%) |
| 2009–2012 | 35 | 10 (28.6%) | 9 (25.7%) | 7 (20.0%) | 2 (5.7%) |
| 2013–2016 | 62 | 27 (43.5%) | 22 (35.5%) | 14 (22.6%) | 6 (9.7%) |
| 2017–2020 | 42 | 24 (57.1%) | 20 (47.6%) | 14 (33.3%) | 7 (16.7%) |
| 2021–2023 | 31 | 17 (54.8%) | 11 (35.5%) | 8 (25.8%) | 3 (9.7%) |

![Confronto fra periodi](periodi.png)

Nel materiale 2017–2020 la quota codificata di urti è 24/42, contro 27/62 nel 2013–2016; nel 2021–2023 è 17/31. Si osserva quindi una maggiore presenza relativa di problemi con urti nelle raccolte più recenti considerate, ma il confronto è **sensibile alla selezione delle prove, alle tracce mancanti e al tipo di esercizi raccolti**. Non si conclude che il docente abbia modificato il programma o che un urto sia più probabile al prossimo appello. Non vengono applicati test di significatività che presupporrebbero un campionamento non dimostrato.

## 6. Sensibilità ai frammenti

L’analisi estesa usa N=211. I dieci record aggiuntivi non sono imputazioni di dati mancanti: sono casi parzialmente leggibili, mantenuti con la loro incertezza. Le categorie non identificabili restano assenti, non negative.

| Categoria | Principale n/201 | Estesa n/211 |
| --- | --- | --- |
| Urti | 87 (43.3%) | 90 (42.7%) |
| Interazioni e riflessioni | 72 (35.8%) | 76 (36.0%) |
| Onde continue | 52 (25.9%) | 56 (26.5%) |
| Contatti | 23 (11.4%) | 23 (10.9%) |
| Getti e adattamento | 21 (10.4%) | 21 (10.0%) |
| Formazione di urti | 11 (5.5%) | 12 (5.7%) |
| Transonico e linea sonica | 11 (5.5%) | 11 (5.2%) |
| Strato limite e viscosita | 8 (4.0%) | 8 (3.8%) |
| Acustica | 5 (2.5%) | 5 (2.4%) |
| Choking e condizioni critiche | 5 (2.5%) | 5 (2.4%) |
| Miscelamento | 3 (1.5%) | 3 (1.4%) |

L’ordine dei principali fenomeni non cambia. Questo rende la priorità generale relativamente stabile rispetto ai frammenti, ma non elimina l’incertezza documentale né quella dovuta alla selezione degli appelli.

## 7. Mappa delle priorità di studio

Le tre dimensioni sono distinte: **F** = frequenza osservata nella codifica, **C** = centralità teorica, **T** = trasferibilità a scenari nuovi. I giudizi C e T sono valutazioni didattiche, non statistiche. Non vengono trasformati in una fittizia probabilità di voto.

| Priorità | Blocco | F osservata | C | T | Obiettivo di padronanza |
| --- | --- | --- | --- | --- | --- |
| 0 — Fondamenti | Bilanci, termodinamica, Eulero e ipotesi | Trasversale: equazioni locali 80 e termodinamica 95 come applicabilità | Molto alta | Molto alta | Ricavare le equazioni, distinguere forma integrale/locale, riferimento e grandezze totali. |
| 1 — Onde 1D | Caratteristiche, Riemann, urti mobili, contatti, pistoni e contorni | Urti 87; interazioni 72; caratteristiche 1D applicabili 76 | Molto alta | Molto alta | Saper costruire gli stati e il piano x-t prima di usare formule tabulate. |
| 1 — Supersonico 2D | Urti obliqui, polari, PM, compatibilità, riflessioni e getti | Polari 41 e PM/char2D 63 come applicabilità | Molto alta | Molto alta | Separare regioni uniformi/continue, selezionare rami e controllare la geometria. |
| 2 — Condotti e sistemi | Quasi-1D, choking, serbatoi, spinta e bilanci integrali | Quasi-1D/area 23; bilanci integrali 28 | Alta | Alta | Riconoscere regime, vincoli e validità del modello ridotto. |
| 2 — Numerica essenziale | Godunov, conservazione, stabilità, condizioni al contorno | 17 richieste di algoritmo/metodo | Alta | Alta | Capire il solutore locale e l’aggiornamento, non memorizzare uno schema isolato. |
| 3 — Temi mirati | Transonico/potenziale, viscosità, acustica, modelli non standard | Categorie meno rappresentate ma presenti | Variabile | Alta per alcuni | Un esercizio rappresentativo per ciascun modello; non eliminarli dal programma personale. |

### Strategia quando il tempo è limitato

Il nucleo iniziale non deve essere una lista di formule: per ogni scenario occorre riconoscere regime e fenomeno, formulare le ipotesi, scegliere le leggi, identificare le incognite e i vincoli, costruire una soluzione e controllarne la coerenza fisica. L’allenamento va fatto cambiando una sola ipotesi alla volta: parete aperta invece che chiusa, compressione invece che espansione, stato a monte non uniforme, gas diverso, riferimento mobile, oppure onda incidente su un contatto.

È più efficiente studiare insieme caratteristiche e Rankine–Hugoniot, e insieme polari e Prandtl–Meyer, anziché trattarli come capitoli indipendenti. I modelli meno frequenti richiedono almeno una conoscenza di base sufficiente a riconoscerli: in particolare transonico/potenziale, strato limite, acustica e sistemi con sorgenti. La selezione degli esercizi sotto è una proposta di allenamento, non una garanzia di copertura dell’esame.

## 8. Esercizi rappresentativi

Per gli ID G usare il catalogo 2003–2023; per gli ID F usare il catalogo 2017–2020. I riferimenti indicano le pagine del PDF originale, non la numerazione stampata. I problemi sono selezionati per la combinazione di competenze e per la possibilità di variare lo scenario. Nessuna soluzione viene qui svolta.

| ID | Blocco | Pagine del testo | Perché allenarlo |
| --- | --- | --- | --- |
| G/B01-E01 | Fondamenti | 2003_2023: ['1'] | Bilanci locali, quasi-1D, gradienti e confronto comprimibile/incomprimibile. |
| G/B01-E03 | Fondamenti | 2003_2023: ['4'] | Riemann scalare, caratteristiche, salto e ammissibilità. |
| G/B03-E02 | Supersonico 2D | 2003_2023: ['13'] | Urto obliquo, riflessione di Mach, polari e odografa. |
| G/B04-E01 | Onde 1D | 2003_2023: ['17'] | Apertura di setto, condizioni al contorno e interazione con contatto. |
| G/B16-E03 | Onde 1D | 2003_2023: ['50'] | Trasmissione/riflessione onda-contatto e condizione di onda nulla. |
| F/B08-E01 | Onde 1D | 2017_2020: ['52', '56'] | Pistone con leggi alternative, coalescenza e diagrammi; cambiare legge di moto come test di trasferibilità. |
| G/B51-E01 | Supersonico 2D | 2003_2023: ['148'] | Due correnti, urti/espansioni, confronto isentropico/esatto e bilancio di miscelamento. |
| G/B67-E01 | Trasversale | 2003_2023: ['249', '250'] | Transonico, gola, shock, separazione e confronto di modelli: esercizio ponte. |
| F/B11-E02 | Supersonico 2D | 2017_2020: ['69', '71'] | Getto sottoespanso con e senza corpo centrale; stessa fisica, geometria diversa. |
| G/B44-E03 | Numerico | 2003_2023: ['116'] | Godunov e confronto con soluzione esatta; svolgimento associato parziale, da non assumere corretto. |
| G/B33-E03 | Fondamenti | 2003_2023: ['82'] | Adimensionalizzazione di Eulero e scelta delle scale. |
| G/B62-E02 | Modelli alternativi | 2003_2023: ['240'] | Flusso omotermo, salti e confronto con il modello adiabatico. |
| G/B37-E02 | Misure | 2003_2023: ['91'] | Pitot supersonico, relazione di misura e limiti del Mach. |
| G/B74-E02 | Acustica | 2003_2023: ['280'] | Caratteristiche, condizioni al contorno e frequenza fondamentale. |
| F/B11-E01 | Problema avanzato | 2017_2020: ['68'] | Pistone-gas accoppiato, onde e metodo approssimato; usare la traccia, non il solo svolgimento parziale. |

**Percorso minimo consigliato:** iniziare con G/B01-E01, G/B01-E03, G/B04-E01, G/B03-E02, G/B16-E03 e F/B08-E01. Passare poi a G/B51-E01, G/B67-E01 e F/B11-E02. Le prove di numerica e i temi mirati completano il percorso. F/B11-E01 è un problema avanzato da usare dopo aver consolidato gli strumenti, perché lo svolgimento disponibile è soltanto un’impostazione iniziale.

## 9. Limiti e verifiche future

La fonte è una raccolta di materiali, non un registro ufficiale completo degli appelli. Non sono noti i criteri di selezione dei problemi, la distribuzione dei punteggi, il tempo assegnato o eventuali modifiche del corso. Le 75 date attestano soltanto le sessioni riconoscibili nel sottoinsieme con data completa. Le categorie sono multilabel e alcune sono inferite; i conteggi di applicabilità dei modelli non sono conteggi di richieste esplicite. Le soluzioni manoscritte non sono state corrette, quindi nessun valore numerico in esse è assunto come verità fisica.

I controlli documentali prioritari restano B10–B11 per le date, B25 per la data oscurata, B28 per il taglio, le varianti B53 e le associazioni di svolgimento fuori posizione. L’inserto p.133 rimane una copia delle tre tracce B41. Le due righe vuote del catalogo Franco non possono essere attribuite a temi. Un’eventuale futura verifica fisica dovrà essere separata da questa statistica.

## 10. Tracciabilità e file

`dataset.csv` e `dataset.json` conservano ID, date, pagine, scenario, richieste, principi sorgente, svolgimenti, stato documentale, tag ed evidenze. `crosswalk.csv` documenta tutti i 29 record della tabella Franco. `frequencies.csv`, `cooccurrences.csv`, `coverage.csv` e `period_frequencies.csv` contengono i denominatori e **tutti gli ID** che contribuiscono a ogni conteggio. Il codice permette di riprodurre la codifica e modificare le regole senza alterare i testi sorgente.

### Appendice A — ID completi per ciascuna frequenza principale

#### Fenomeni

**Urti — 87/201.** G/B03-E02, G/B05-E01, G/B06-E01, G/B06-E02, G/B07-E01, G/B08-E02, G/B09-E01, G/B11-E04, G/B11-E05, G/B14-E03, G/B15-E02, G/B17-E02, G/B19-E01, G/B20-E02, G/B21-E01, G/B22-E01, G/B23-E02, G/B26-E01, G/B27-E02, G/B29-E01, G/B30-E01, G/B30-E02, G/B31-E01, G/B32-E02, G/B33-E01, G/B33-E02, G/B34-E01, G/B34-E03, G/B35-E03, G/B36-E01, G/B36-E03, G/B38-E01, G/B39-E02, G/B40-E03, G/B41-E02, G/B41-E03, G/B44-E01, G/B45-E02, G/B46-E02, G/B46-E03, G/B47-E02, G/B47-E03, G/B49-E01, G/B50-E02, G/B50-E03, G/B51-E02, G/B52-E03, G/B53-E01, G/B53-E03, G/B53-E05, G/B53-E06, G/B54-E02, G/B54-E03, G/B55-E02, G/B56-E01, G/B56-E02, G/B59-E03, G/B60-E02, G/B62-E01, G/B62-E02, G/B63-E01, G/B63-E02, G/B64-E01, G/B64-E03, G/B65-E02, G/B65-E03, G/B66-E02, G/B67-E01, G/B68-E01, G/B69-E01, G/B69-E02, G/B70-E01, G/B72-E01, G/B72-E02, G/B74-E01, F/B06-E01, F/B06-E02, F/B06-E03, F/B07-E01, F/B07-E02, F/B07-E03, F/B08-E01, F/B08-E02, F/B09-E01, F/B10-E01, F/B11-E01, F/B11-E02.

**Interazioni e riflessioni — 72/201.** G/B01-E02, G/B03-E02, G/B04-E01, G/B05-E01, G/B06-E01, G/B07-E01, G/B10-E05, G/B11-E01, G/B11-E03, G/B11-E05, G/B13-E03, G/B16-E03, G/B17-E02, G/B18-E01, G/B19-E01, G/B19-E02, G/B21-E01, G/B22-E01, G/B27-E02, G/B29-E01, G/B30-E01, G/B30-E02, G/B32-E02, G/B34-E01, G/B34-E02, G/B37-E01, G/B38-E01, G/B39-E01, G/B39-E02, G/B40-E02, G/B40-E03, G/B41-E02, G/B41-E03, G/B46-E02, G/B46-E03, G/B47-E01, G/B47-E03, G/B48-E02, G/B49-E03, G/B50-E03, G/B51-E02, G/B52-E01, G/B52-E03, G/B53-E01, G/B53-E02, G/B53-E05, G/B53-E06, G/B54-E02, G/B55-E02, G/B56-E01, G/B60-E02, G/B61-E01, G/B62-E01, G/B67-E02, G/B68-E01, G/B68-E02, G/B69-E01, G/B70-E01, G/B70-E02, G/B72-E01, G/B72-E02, G/B74-E01, F/B06-E01, F/B06-E03, F/B07-E01, F/B07-E02, F/B07-E03, F/B08-E02, F/B10-E01, F/B10-E02, F/B11-E01, F/B11-E02.

**Onde continue — 52/201.** G/B01-E02, G/B01-E03, G/B02-E01, G/B04-E01, G/B05-E01, G/B06-E01, G/B06-E02, G/B10-E01, G/B10-E02, G/B15-E03, G/B16-E03, G/B18-E02, G/B23-E02, G/B26-E01, G/B27-E03, G/B27-E04, G/B33-E01, G/B33-E02, G/B34-E01, G/B35-E03, G/B36-E01, G/B38-E02, G/B40-E03, G/B42-E01, G/B43-E02, G/B45-E02, G/B47-E01, G/B49-E01, G/B49-E02, G/B50-E02, G/B52-E01, G/B52-E02, G/B53-E02, G/B54-E02, G/B57-E01, G/B58-E02, G/B59-E02, G/B60-E01, G/B61-E02, G/B64-E01, G/B64-E03, G/B68-E02, G/B70-E01, G/B70-E02, G/B71-E02, G/B72-E01, F/B07-E01, F/B07-E02, F/B08-E01, F/B09-E02, F/B10-E02, F/B11-E02.

**Contatti — 23/201.** G/B04-E01, G/B07-E01, G/B09-E02, G/B13-E03, G/B16-E03, G/B26-E01, G/B27-E02, G/B28-E03, G/B30-E01, G/B30-E02, G/B31-E01, G/B33-E02, G/B50-E02, G/B53-E03, G/B54-E02, G/B55-E02, G/B56-E01, G/B62-E01, G/B68-E01, G/B70-E01, F/B06-E03, F/B07-E01, F/B10-E01.

**Getti e adattamento — 21/201.** G/B01-E02, G/B10-E02, G/B12-E02, G/B13-E02, G/B25-E01, G/B28-E03, G/B31-E01, G/B35-E02, G/B41-E01, G/B42-E03, G/B49-E02, G/B52-E01, G/B53-E03, G/B54-E03, G/B58-E02, G/B60-E03, G/B64-E04, G/B69-E03, G/B72-E01, F/B10-E01, F/B11-E02.

**Formazione di urti — 11/201.** G/B06-E01, G/B15-E02, G/B35-E03, G/B36-E01, G/B45-E02, G/B49-E01, G/B64-E03, G/B69-E02, F/B07-E02, F/B08-E01, F/B09-E02.

**Transonico e linea sonica — 11/201.** G/B07-E02, G/B08-E01, G/B09-E01, G/B26-E01, G/B38-E02, G/B60-E01, G/B61-E02, G/B63-E01, G/B67-E01, F/B09-E01, F/B09-E02.

**Strato limite e viscosita — 8/201.** G/B05-E02, G/B06-E02, G/B10-E03, G/B10-E04, G/B12-E02, G/B15-E01, G/B29-E01, G/B50-E01.

**Acustica — 5/201.** G/B58-E03, G/B67-E02, G/B70-E02, G/B74-E02, F/B09-E02.

**Choking e condizioni critiche — 5/201.** G/B02-E01, G/B08-E01, G/B35-E01, G/B59-E02, F/B06-E01.

**Miscelamento — 3/201.** G/B12-E02, G/B15-E01, G/B51-E01.

#### Modelli e strumenti

**Termodinamica e grandezze totali — 95/201.** G/B01-E01, G/B01-E02, G/B02-E01, G/B02-E02, G/B03-E01, G/B04-E01, G/B06-E01, G/B07-E01, G/B07-E02, G/B08-E01, G/B08-E02, G/B09-E01, G/B10-E04, G/B11-E01, G/B11-E03, G/B11-E05, G/B12-E01, G/B13-E03, G/B15-E03, G/B16-E01, G/B16-E02, G/B17-E01, G/B18-E01, G/B19-E02, G/B21-E02, G/B22-E01, G/B24-E02, G/B26-E01, G/B27-E01, G/B27-E02, G/B27-E03, G/B27-E04, G/B28-E03, G/B29-E02, G/B30-E02, G/B31-E01, G/B31-E02, G/B32-E01, G/B34-E01, G/B34-E02, G/B35-E01, G/B35-E02, G/B36-E01, G/B36-E02, G/B37-E01, G/B38-E02, G/B40-E02, G/B43-E01, G/B43-E02, G/B44-E01, G/B44-E02, G/B45-E01, G/B46-E01, G/B47-E01, G/B49-E01, G/B51-E01, G/B52-E01, G/B52-E02, G/B53-E01, G/B53-E02, G/B53-E04, G/B54-E02, G/B55-E02, G/B57-E01, G/B58-E01, G/B58-E03, G/B59-E02, G/B60-E01, G/B62-E02, G/B63-E02, G/B64-E01, G/B64-E03, G/B64-E04, G/B65-E01, G/B65-E03, G/B66-E01, G/B67-E01, G/B68-E01, G/B68-E02, G/B69-E02, G/B69-E03, G/B70-E02, G/B71-E02, G/B72-E01, G/B73-E01, G/B74-E02, F/B06-E01, F/B06-E02, F/B07-E01, F/B07-E02, F/B08-E01, F/B08-E02, F/B10-E01, F/B10-E02, F/B11-E02.

**Rankine-Hugoniot — 93/201.** G/B03-E02, G/B03-E03, G/B05-E01, G/B06-E01, G/B06-E02, G/B07-E01, G/B08-E02, G/B09-E01, G/B11-E04, G/B11-E05, G/B13-E01, G/B14-E03, G/B15-E02, G/B17-E02, G/B19-E01, G/B20-E02, G/B21-E01, G/B22-E01, G/B23-E02, G/B26-E01, G/B27-E02, G/B27-E03, G/B27-E04, G/B29-E01, G/B30-E01, G/B30-E02, G/B31-E01, G/B32-E02, G/B33-E01, G/B33-E02, G/B34-E01, G/B34-E03, G/B35-E03, G/B36-E01, G/B36-E03, G/B38-E01, G/B39-E02, G/B40-E03, G/B41-E02, G/B41-E03, G/B44-E01, G/B45-E02, G/B46-E02, G/B46-E03, G/B47-E02, G/B47-E03, G/B49-E01, G/B50-E02, G/B50-E03, G/B51-E02, G/B52-E03, G/B53-E01, G/B53-E03, G/B53-E05, G/B53-E06, G/B54-E01, G/B54-E02, G/B54-E03, G/B55-E02, G/B56-E01, G/B56-E02, G/B59-E03, G/B60-E02, G/B62-E01, G/B62-E02, G/B63-E01, G/B63-E02, G/B64-E01, G/B64-E02, G/B64-E03, G/B65-E02, G/B65-E03, G/B66-E02, G/B67-E01, G/B68-E01, G/B69-E01, G/B69-E02, G/B70-E01, G/B72-E01, G/B72-E02, G/B74-E01, F/B06-E01, F/B06-E02, F/B06-E03, F/B07-E01, F/B07-E02, F/B07-E03, F/B08-E01, F/B08-E02, F/B09-E01, F/B10-E01, F/B11-E01, F/B11-E02.

**Eulero e bilanci locali — 80/201.** G/B01-E03, G/B02-E02, G/B03-E01, G/B03-E03, G/B04-E01, G/B05-E01, G/B06-E01, G/B07-E01, G/B09-E02, G/B10-E01, G/B10-E04, G/B10-E05, G/B12-E02, G/B13-E01, G/B13-E03, G/B14-E01, G/B14-E02, G/B15-E02, G/B16-E03, G/B17-E02, G/B20-E01, G/B21-E02, G/B22-E01, G/B23-E02, G/B24-E01, G/B26-E01, G/B27-E02, G/B29-E02, G/B30-E02, G/B31-E02, G/B32-E01, G/B33-E01, G/B33-E02, G/B33-E03, G/B35-E03, G/B36-E01, G/B41-E03, G/B42-E01, G/B42-E02, G/B43-E01, G/B43-E02, G/B44-E02, G/B45-E01, G/B45-E02, G/B46-E01, G/B47-E01, G/B47-E03, G/B48-E01, G/B49-E01, G/B50-E01, G/B50-E02, G/B52-E02, G/B53-E02, G/B53-E04, G/B53-E05, G/B54-E01, G/B55-E01, G/B57-E01, G/B58-E02, G/B58-E03, G/B59-E01, G/B59-E02, G/B60-E01, G/B61-E02, G/B62-E02, G/B64-E02, G/B64-E03, G/B68-E02, G/B69-E02, G/B70-E02, G/B71-E02, G/B72-E02, G/B73-E01, F/B06-E01, F/B07-E01, F/B07-E02, F/B08-E01, F/B09-E02, F/B10-E02, F/B11-E01.

**Caratteristiche 1D — 76/201.** G/B01-E03, G/B02-E02, G/B03-E03, G/B04-E01, G/B05-E01, G/B06-E01, G/B07-E01, G/B09-E02, G/B10-E01, G/B10-E05, G/B13-E01, G/B13-E03, G/B14-E02, G/B15-E02, G/B16-E03, G/B17-E02, G/B19-E02, G/B20-E01, G/B21-E02, G/B22-E01, G/B23-E02, G/B24-E01, G/B26-E01, G/B27-E02, G/B29-E02, G/B30-E02, G/B31-E01, G/B33-E01, G/B34-E01, G/B35-E03, G/B36-E01, G/B36-E03, G/B39-E02, G/B41-E03, G/B42-E01, G/B43-E02, G/B44-E01, G/B44-E02, G/B45-E02, G/B47-E01, G/B47-E02, G/B48-E01, G/B48-E03, G/B49-E01, G/B50-E02, G/B51-E02, G/B52-E02, G/B52-E03, G/B53-E02, G/B53-E04, G/B53-E05, G/B55-E01, G/B56-E02, G/B57-E01, G/B58-E03, G/B60-E01, G/B61-E01, G/B62-E02, G/B63-E02, G/B64-E03, G/B65-E01, G/B66-E02, G/B67-E02, G/B68-E02, G/B69-E02, G/B70-E02, G/B71-E02, G/B72-E02, G/B74-E02, F/B06-E01, F/B07-E01, F/B07-E02, F/B08-E01, F/B09-E02, F/B10-E02, F/B11-E01.

**Prandtl-Meyer e caratteristiche 2D — 63/201.** G/B01-E02, G/B03-E02, G/B06-E02, G/B10-E02, G/B11-E03, G/B11-E04, G/B13-E02, G/B14-E03, G/B15-E03, G/B16-E02, G/B17-E01, G/B18-E02, G/B19-E01, G/B20-E02, G/B23-E01, G/B24-E02, G/B25-E01, G/B27-E03, G/B27-E04, G/B29-E01, G/B30-E01, G/B31-E01, G/B33-E02, G/B34-E02, G/B34-E03, G/B35-E02, G/B39-E01, G/B40-E01, G/B41-E01, G/B41-E02, G/B42-E03, G/B45-E03, G/B46-E03, G/B48-E01, G/B49-E02, G/B51-E01, G/B52-E01, G/B53-E01, G/B53-E03, G/B53-E06, G/B54-E03, G/B55-E02, G/B56-E01, G/B58-E02, G/B59-E02, G/B60-E03, G/B63-E01, G/B64-E01, G/B64-E04, G/B65-E02, G/B67-E01, G/B69-E01, G/B69-E03, G/B70-E01, G/B71-E01, G/B72-E01, G/B73-E01, G/B74-E01, F/B07-E03, F/B08-E02, F/B09-E01, F/B10-E01, F/B11-E02.

**Polari e odografa — 41/201.** G/B03-E02, G/B17-E01, G/B19-E01, G/B20-E02, G/B21-E01, G/B25-E01, G/B27-E03, G/B27-E04, G/B29-E01, G/B32-E02, G/B33-E02, G/B34-E03, G/B39-E01, G/B40-E03, G/B41-E02, G/B46-E03, G/B48-E01, G/B51-E01, G/B53-E06, G/B55-E02, G/B56-E01, G/B58-E02, G/B60-E02, G/B62-E01, G/B63-E01, G/B64-E01, G/B65-E02, G/B67-E01, G/B68-E01, G/B69-E01, G/B70-E01, G/B71-E01, G/B73-E01, G/B74-E01, F/B06-E02, F/B06-E03, F/B07-E02, F/B07-E03, F/B09-E01, F/B10-E01, F/B11-E02.

**Problema di Riemann — 40/201.** G/B01-E03, G/B04-E01, G/B07-E01, G/B09-E02, G/B10-E05, G/B11-E02, G/B11-E05, G/B13-E03, G/B14-E02, G/B16-E03, G/B17-E02, G/B22-E01, G/B26-E01, G/B27-E02, G/B30-E02, G/B31-E01, G/B33-E02, G/B34-E01, G/B36-E03, G/B38-E02, G/B39-E02, G/B41-E03, G/B44-E01, G/B44-E02, G/B44-E03, G/B46-E02, G/B47-E02, G/B48-E01, G/B48-E03, G/B50-E02, G/B51-E02, G/B52-E03, G/B53-E05, G/B54-E02, G/B56-E01, G/B57-E01, G/B62-E02, G/B72-E02, F/B06-E01, F/B07-E01.

**Bilanci integrali — 28/201.** G/B02-E01, G/B03-E01, G/B06-E02, G/B10-E03, G/B10-E04, G/B11-E01, G/B12-E01, G/B12-E02, G/B13-E01, G/B15-E01, G/B17-E02, G/B18-E01, G/B23-E01, G/B27-E01, G/B27-E02, G/B35-E01, G/B35-E02, G/B36-E02, G/B37-E01, G/B40-E02, G/B47-E02, G/B48-E03, G/B51-E01, G/B59-E01, G/B65-E02, G/B66-E01, G/B73-E01, F/B06-E01.

**Quasi-1D e area-Mach — 23/201.** G/B01-E01, G/B01-E02, G/B02-E01, G/B08-E02, G/B10-E04, G/B11-E03, G/B15-E03, G/B16-E01, G/B17-E01, G/B24-E02, G/B27-E04, G/B46-E02, G/B48-E03, G/B49-E02, G/B53-E05, G/B54-E02, G/B59-E02, G/B64-E01, G/B66-E01, G/B67-E01, G/B72-E02, F/B06-E01, F/B08-E02.

**Metodi numerici — 13/201.** G/B02-E03, G/B04-E02, G/B05-E02, G/B07-E02, G/B11-E02, G/B11-E05, G/B14-E02, G/B27-E03, G/B33-E01, G/B38-E02, G/B39-E01, G/B44-E03, G/B57-E01.

**Sorgenti e modelli non standard — 11/201.** G/B03-E03, G/B09-E02, G/B10-E04, G/B11-E02, G/B13-E01, G/B27-E01, G/B32-E01, G/B43-E01, G/B45-E01, G/B46-E01, G/B62-E02.

**Potenziale e modelli ridotti — 6/201.** G/B04-E02, G/B07-E02, G/B08-E01, G/B42-E01, G/B61-E02, G/B67-E01.

**Adimensionalizzazione — 4/201.** G/B22-E01, G/B33-E03, G/B43-E02, F/B08-E01.

#### Azioni

**Determinare stati e grandezze — 139/201.** G/B01-E01, G/B01-E02, G/B01-E03, G/B02-E01, G/B02-E02, G/B03-E01, G/B03-E02, G/B04-E01, G/B04-E02, G/B05-E01, G/B06-E01, G/B06-E02, G/B07-E01, G/B08-E01, G/B08-E02, G/B09-E01, G/B09-E02, G/B10-E01, G/B10-E02, G/B10-E03, G/B10-E04, G/B11-E01, G/B11-E03, G/B11-E04, G/B11-E05, G/B12-E01, G/B12-E02, G/B13-E02, G/B13-E03, G/B14-E02, G/B14-E03, G/B15-E01, G/B15-E02, G/B16-E01, G/B16-E03, G/B17-E01, G/B17-E02, G/B18-E01, G/B18-E02, G/B19-E01, G/B19-E02, G/B20-E01, G/B21-E01, G/B21-E02, G/B22-E01, G/B23-E01, G/B23-E02, G/B24-E02, G/B25-E01, G/B26-E01, G/B27-E01, G/B27-E02, G/B27-E03, G/B27-E04, G/B28-E03, G/B29-E01, G/B29-E02, G/B30-E02, G/B31-E01, G/B32-E02, G/B33-E01, G/B34-E01, G/B35-E01, G/B35-E03, G/B36-E01, G/B36-E02, G/B37-E01, G/B37-E02, G/B39-E01, G/B39-E02, G/B40-E02, G/B40-E03, G/B42-E01, G/B42-E03, G/B43-E01, G/B43-E02, G/B44-E01, G/B44-E03, G/B45-E01, G/B45-E02, G/B46-E02, G/B46-E03, G/B47-E01, G/B47-E02, G/B48-E02, G/B48-E03, G/B49-E01, G/B49-E02, G/B49-E03, G/B51-E01, G/B52-E01, G/B52-E02, G/B53-E01, G/B53-E02, G/B53-E04, G/B53-E05, G/B53-E06, G/B54-E02, G/B55-E01, G/B55-E03, G/B56-E01, G/B56-E02, G/B57-E02, G/B58-E01, G/B58-E02, G/B58-E03, G/B59-E01, G/B59-E02, G/B60-E01, G/B60-E02, G/B62-E01, G/B63-E02, G/B64-E01, G/B64-E03, G/B65-E01, G/B65-E02, G/B67-E01, G/B67-E02, G/B68-E01, G/B68-E02, G/B71-E01, G/B71-E02, G/B72-E02, G/B73-E01, G/B74-E01, G/B74-E02, F/B06-E01, F/B06-E02, F/B07-E01, F/B07-E02, F/B07-E03, F/B08-E01, F/B08-E02, F/B09-E01, F/B09-E02, F/B10-E01, F/B10-E02, F/B11-E01, F/B11-E02.

**Costruire diagrammi e profili — 56/201.** G/B01-E02, G/B03-E02, G/B07-E01, G/B09-E02, G/B10-E05, G/B15-E03, G/B17-E01, G/B19-E01, G/B19-E02, G/B20-E01, G/B20-E02, G/B21-E01, G/B22-E01, G/B24-E01, G/B25-E01, G/B26-E01, G/B27-E04, G/B29-E01, G/B29-E02, G/B30-E02, G/B31-E01, G/B33-E01, G/B34-E01, G/B34-E03, G/B39-E01, G/B39-E02, G/B41-E03, G/B46-E02, G/B46-E03, G/B47-E02, G/B48-E02, G/B50-E02, G/B53-E01, G/B53-E05, G/B54-E02, G/B54-E03, G/B55-E02, G/B56-E01, G/B60-E02, G/B61-E01, G/B62-E01, G/B63-E01, G/B65-E01, G/B69-E01, G/B70-E01, G/B71-E01, G/B72-E02, G/B73-E01, G/B74-E01, F/B06-E03, F/B07-E02, F/B07-E03, F/B08-E02, F/B09-E01, F/B10-E01, F/B11-E02.

**Interpretare osservazioni — 36/201.** G/B04-E02, G/B29-E01, G/B29-E02, G/B30-E01, G/B34-E03, G/B35-E02, G/B36-E03, G/B37-E02, G/B39-E02, G/B40-E01, G/B42-E02, G/B42-E03, G/B47-E03, G/B50-E03, G/B52-E03, G/B53-E03, G/B54-E03, G/B55-E03, G/B57-E02, G/B58-E01, G/B59-E03, G/B60-E03, G/B62-E01, G/B63-E01, G/B64-E04, G/B67-E01, G/B68-E01, G/B70-E01, G/B71-E01, G/B72-E01, F/B07-E01, F/B07-E03, F/B08-E02, F/B09-E01, F/B10-E01, F/B11-E02.

**Ammissibilita e verifica — 31/201.** G/B01-E02, G/B01-E03, G/B03-E03, G/B11-E03, G/B14-E01, G/B16-E03, G/B19-E01, G/B24-E02, G/B37-E02, G/B38-E02, G/B43-E01, G/B44-E03, G/B47-E03, G/B49-E02, G/B50-E03, G/B51-E01, G/B52-E02, G/B52-E03, G/B53-E03, G/B53-E06, G/B57-E01, G/B60-E01, G/B62-E01, G/B62-E02, G/B64-E01, G/B67-E01, G/B74-E01, F/B06-E02, F/B07-E02, F/B08-E01, F/B09-E02.

**Impostare modello e condizioni — 30/201.** G/B02-E03, G/B04-E02, G/B05-E02, G/B07-E02, G/B10-E04, G/B11-E02, G/B13-E01, G/B16-E02, G/B20-E02, G/B23-E01, G/B24-E01, G/B27-E03, G/B31-E02, G/B32-E01, G/B33-E03, G/B37-E01, G/B39-E01, G/B44-E02, G/B46-E01, G/B48-E01, G/B50-E01, G/B54-E01, G/B59-E01, G/B61-E02, G/B62-E02, G/B64-E02, G/B66-E01, G/B66-E02, G/B68-E02, F/B11-E01.

**Confrontare modelli e parametri — 25/201.** G/B01-E01, G/B02-E01, G/B14-E01, G/B19-E01, G/B27-E01, G/B37-E02, G/B38-E02, G/B42-E02, G/B43-E01, G/B44-E03, G/B47-E03, G/B49-E02, G/B50-E03, G/B51-E01, G/B53-E06, G/B57-E01, G/B60-E01, G/B62-E01, G/B62-E02, G/B67-E01, G/B68-E01, G/B69-E02, F/B06-E02, F/B07-E02, F/B08-E01.

**Derivare e giustificare — 25/201.** G/B02-E03, G/B03-E03, G/B10-E04, G/B13-E01, G/B14-E01, G/B16-E02, G/B31-E02, G/B32-E01, G/B33-E03, G/B34-E01, G/B38-E02, G/B43-E01, G/B45-E01, G/B46-E01, G/B50-E01, G/B54-E01, G/B56-E02, G/B61-E02, G/B62-E02, G/B63-E02, G/B64-E02, G/B65-E03, G/B66-E02, G/B74-E02, F/B07-E02.

**Costruire algoritmo numerico — 17/201.** G/B02-E03, G/B04-E02, G/B05-E02, G/B07-E02, G/B11-E02, G/B11-E05, G/B14-E02, G/B20-E02, G/B27-E03, G/B33-E01, G/B38-E02, G/B39-E01, G/B44-E03, G/B48-E01, G/B57-E01, G/B66-E01, F/B11-E01.

**Progettare e ottimizzare — 7/201.** G/B01-E02, G/B11-E03, G/B15-E03, G/B16-E02, G/B24-E02, G/B58-E02, G/B64-E01.

#### Tipologie

**Quantitativo — 141/201.** G/B01-E01, G/B01-E02, G/B01-E03, G/B02-E01, G/B02-E02, G/B03-E01, G/B03-E02, G/B04-E01, G/B04-E02, G/B05-E01, G/B06-E01, G/B06-E02, G/B07-E01, G/B08-E01, G/B08-E02, G/B09-E01, G/B09-E02, G/B10-E01, G/B10-E02, G/B10-E03, G/B10-E04, G/B11-E01, G/B11-E03, G/B11-E04, G/B11-E05, G/B12-E01, G/B12-E02, G/B13-E02, G/B13-E03, G/B14-E02, G/B14-E03, G/B15-E01, G/B15-E02, G/B15-E03, G/B16-E01, G/B16-E02, G/B16-E03, G/B17-E01, G/B17-E02, G/B18-E01, G/B18-E02, G/B19-E01, G/B19-E02, G/B20-E01, G/B21-E01, G/B21-E02, G/B22-E01, G/B23-E01, G/B23-E02, G/B24-E02, G/B25-E01, G/B26-E01, G/B27-E01, G/B27-E02, G/B27-E03, G/B27-E04, G/B28-E03, G/B29-E01, G/B29-E02, G/B30-E02, G/B31-E01, G/B32-E02, G/B33-E01, G/B34-E01, G/B35-E01, G/B35-E03, G/B36-E01, G/B36-E02, G/B37-E01, G/B37-E02, G/B39-E01, G/B39-E02, G/B40-E02, G/B40-E03, G/B42-E01, G/B42-E03, G/B43-E01, G/B43-E02, G/B44-E01, G/B44-E03, G/B45-E01, G/B45-E02, G/B46-E02, G/B46-E03, G/B47-E01, G/B47-E02, G/B48-E02, G/B48-E03, G/B49-E01, G/B49-E02, G/B49-E03, G/B51-E01, G/B52-E01, G/B52-E02, G/B53-E01, G/B53-E02, G/B53-E04, G/B53-E05, G/B53-E06, G/B54-E02, G/B55-E01, G/B55-E03, G/B56-E01, G/B56-E02, G/B57-E02, G/B58-E01, G/B58-E02, G/B58-E03, G/B59-E01, G/B59-E02, G/B60-E01, G/B60-E02, G/B62-E01, G/B63-E02, G/B64-E01, G/B64-E03, G/B65-E01, G/B65-E02, G/B67-E01, G/B67-E02, G/B68-E01, G/B68-E02, G/B71-E01, G/B71-E02, G/B72-E02, G/B73-E01, G/B74-E01, G/B74-E02, F/B06-E01, F/B06-E02, F/B07-E01, F/B07-E02, F/B07-E03, F/B08-E01, F/B08-E02, F/B09-E01, F/B09-E02, F/B10-E01, F/B10-E02, F/B11-E01, F/B11-E02.

**Qualitativo — 120/201.** G/B01-E03, G/B02-E03, G/B03-E02, G/B03-E03, G/B04-E02, G/B05-E02, G/B07-E01, G/B09-E02, G/B10-E05, G/B11-E01, G/B11-E02, G/B11-E04, G/B14-E03, G/B16-E03, G/B17-E01, G/B17-E02, G/B18-E02, G/B19-E01, G/B19-E02, G/B20-E01, G/B20-E02, G/B21-E01, G/B22-E01, G/B23-E01, G/B24-E01, G/B25-E01, G/B26-E01, G/B27-E02, G/B27-E04, G/B28-E03, G/B29-E01, G/B29-E02, G/B30-E01, G/B30-E02, G/B31-E01, G/B31-E02, G/B32-E01, G/B33-E01, G/B33-E02, G/B34-E01, G/B34-E02, G/B34-E03, G/B35-E02, G/B35-E03, G/B36-E03, G/B37-E01, G/B37-E02, G/B37-E03, G/B38-E01, G/B39-E01, G/B39-E02, G/B40-E01, G/B40-E02, G/B41-E01, G/B41-E02, G/B41-E03, G/B42-E02, G/B42-E03, G/B44-E01, G/B44-E02, G/B45-E01, G/B45-E02, G/B45-E03, G/B46-E02, G/B46-E03, G/B47-E01, G/B47-E02, G/B47-E03, G/B48-E02, G/B48-E03, G/B49-E03, G/B50-E02, G/B50-E03, G/B51-E01, G/B51-E02, G/B52-E02, G/B52-E03, G/B53-E02, G/B53-E03, G/B54-E01, G/B54-E02, G/B54-E03, G/B55-E01, G/B55-E02, G/B55-E03, G/B56-E01, G/B57-E02, G/B58-E01, G/B59-E03, G/B60-E02, G/B60-E03, G/B61-E01, G/B61-E02, G/B62-E02, G/B63-E01, G/B64-E04, G/B66-E01, G/B67-E01, G/B68-E01, G/B68-E02, G/B69-E01, G/B69-E02, G/B69-E03, G/B70-E02, G/B71-E01, G/B71-E02, G/B72-E01, G/B73-E01, G/B74-E02, F/B06-E02, F/B06-E03, F/B07-E02, F/B07-E03, F/B08-E01, F/B08-E02, F/B09-E01, F/B09-E02, F/B10-E01, F/B10-E02, F/B11-E01.

**Teorico — 41/201.** G/B02-E03, G/B03-E03, G/B04-E02, G/B05-E02, G/B07-E02, G/B10-E04, G/B11-E02, G/B13-E01, G/B14-E01, G/B16-E02, G/B20-E02, G/B23-E01, G/B24-E01, G/B27-E03, G/B31-E02, G/B32-E01, G/B33-E03, G/B34-E01, G/B37-E01, G/B38-E02, G/B39-E01, G/B43-E01, G/B44-E02, G/B45-E01, G/B46-E01, G/B48-E01, G/B50-E01, G/B54-E01, G/B56-E02, G/B59-E01, G/B61-E02, G/B62-E02, G/B63-E02, G/B64-E02, G/B65-E03, G/B66-E01, G/B66-E02, G/B68-E02, G/B74-E02, F/B07-E02, F/B11-E01.

**Interpretativo — 36/201.** G/B04-E02, G/B29-E01, G/B29-E02, G/B30-E01, G/B34-E03, G/B35-E02, G/B36-E03, G/B37-E02, G/B39-E02, G/B40-E01, G/B42-E02, G/B42-E03, G/B47-E03, G/B50-E03, G/B52-E03, G/B53-E03, G/B54-E03, G/B55-E03, G/B57-E02, G/B58-E01, G/B59-E03, G/B60-E03, G/B62-E01, G/B63-E01, G/B64-E04, G/B67-E01, G/B68-E01, G/B70-E01, G/B71-E01, G/B72-E01, F/B07-E01, F/B07-E03, F/B08-E02, F/B09-E01, F/B10-E01, F/B11-E02.

**Numerico — 17/201.** G/B02-E03, G/B04-E02, G/B05-E02, G/B07-E02, G/B11-E02, G/B11-E05, G/B14-E02, G/B20-E02, G/B27-E03, G/B33-E01, G/B38-E02, G/B39-E01, G/B44-E03, G/B48-E01, G/B57-E01, G/B66-E01, F/B11-E01.

### Appendice B — Fonte e significato delle evidenze

La colonna `explicit_n` è calcolata esclusivamente dai tag classificati E. `inferred_or_source_index_n` comprende S e I. Per i tipi di problema la distinzione E/S/I non è applicata: sono etichette derivate dalle azioni. Per le coppie il denominatore è 201; `share_of_a` e `share_of_b` usano le dimensioni delle rispettive categorie; `jaccard` è una misura descrittiva di sovrapposizione, non una significatività statistica. Gli identificativi canonici G:/F: sono introdotti soltanto per evitare collisioni fra gli ID Bxx-Exx delle due raccolte.

**Fonti prioritarie:** catalogo.md, Analisi Struttura PDF.txt, riepilogo_anomalie.md, parti_da_verificare.md, anomalie_e_verifiche.md, registro_pagine.md, esercizi_2003_2023.pdf ed esercizi_2017_2020.pdf. Per la mappa teorica di riferimento: appunti.pdf. Non sono state utilizzate fonti web per alterare o completare le tracce mancanti.
