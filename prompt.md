# Obiettivo

In questa cartella costruirai e validerai `gasdyn.py`, un modulo Python di calcolo per il corso universitario di Gasdinamica. Il modulo sarà poi usato in un progetto claude.ai per risolvere esercizi d'esame: il ragionamento fisico resterà a carico del modello, il modulo deve garantire che i numeri siano corretti e che i calcoli iterativi (intersezioni di polari, problemi di Riemann, formazione di urti) siano risolti esattamente invece che per tentativi.

Il modulo deve parlare la lingua del corso: notazione, convenzioni di segno, nomi degli stati e algoritmi devono coincidere con quelli delle fonti in `fonti/`, non con quelli di un libro generico.

# Fonti

- `fonti/formulario.tex`: formulario basato sugli appunti del corso, include l'algoritmo per il problema di Riemann. **Fonte primaria** per formule e convenzioni.
- `fonti/summary_lectures.tex`: riassunto delle lezioni. Seconda fonte per teoria e notazione.
- `fonti/catalogo.md`: catalogo di tutte le tracce d'esame 2003–2023, con dati, richieste e risultati numerici riportati negli svolgimenti. Le pagine citate sono pagine fisiche di `fonti/esercizi_2003_2023.pdf`.
- `soluzioni_A/*.pdf`: 17 tracce con soluzione impaginata delle prove 2021–2023, probabilmente del docente ma senza certezza. Sono scansioni: aprile (rasterizzando le pagine) solo quando serve verificare un numero o una convenzione specifica.
- `fonti/esercizi_*.pdf` (se presenti): raccolte complete, stessa regola d'uso.

Le formule vengono dalle fonti. Se una funzione richiede una relazione che non trovi nel formulario o nel summary, puoi usare la teoria standard, ma dichiaralo nel docstring ("fuori formulario") e segnalamelo. Se formulario e summary sono in disaccordo, non scegliere da solo: elenca la discrepanza.

# Organizzazione del lavoro

- Inizializza un repository git e fai un commit alla fine di ogni passo significativo, con messaggi chiari.
- Crea un ambiente virtuale (`.venv`) con `numpy`, `scipy`, `matplotlib`, `pytest`, `pypdf` e quello che ti serve per rasterizzare i PDF. Registra le dipendenze in `requirements.txt`.
- Crea `CLAUDE.md` con: obiettivo, struttura della cartella, convenzioni confermate, stato di avanzamento delle fasi, decisioni prese e problemi aperti. Aggiornalo alla fine di ogni fase e ogni volta che cambia qualcosa di rilevante: se la sessione si interrompe, la successiva deve poter ripartire leggendo solo quel file.
- Struttura di output:
  - `gasdyn.py`
  - `tests/test_gasdyn.py`
  - `verifica/verifica_catalogo.md` (+ eventuali script di verifica e figure)
  - `gasdyn_guida.md`

# Requisiti tecnici del modulo

- Un unico file `gasdyn.py`, solo `numpy`, `scipy`, `matplotlib` come dipendenze (deve funzionare anche nell'ambiente di esecuzione di claude.ai).
- Classe `Gas` con γ e R e grandezze derivate (δ = (γ−1)/2 o la notazione del formulario, cp, cv). Nessun γ = 1.4 implicito. Supporto a gas diversi ai due lati di un contatto.
- Funzionamento sia con grandezze SI sia con dati adimensionali tipici delle tracce (p = 1, ρ = γ, a = 1, ecc.).
- Stati rappresentati in modo uniforme (per esempio una dataclass `State` con ρ, u, p, eventualmente v, e proprietà calcolate a, T, M, s, p0, T0).
- Docstring in italiano: cosa calcola, ipotesi, riferimento all'equazione o sezione del formulario, convenzione di segno.
- Quando una soluzione non esiste o è ambigua (θ > θmax, urto staccato, generazione di vuoto, due rami), sollevare un'eccezione chiara o restituire entrambi i rami etichettati. Mai un numero silenziosamente sbagliato.
- Solutori numerici robusti (`brentq` con bracketing esplicito dove possibile, tolleranze dichiarate).

# Contenuto del modulo

Ordine di priorità indicativo. Ogni area deve coprire i casi che compaiono in `fonti/catalogo.md`.

1. **Termodinamica e grandezze totali**: velocità del suono, entropia (differenze s/R e funzione p/ρ^γ), grandezze totali e critiche, velocità massima. Pressione totale nei modelli incomprimibile, isotermo e isentropico.
2. **Isentropico e quasi-1D**: rapporti isentropici diretti e inversi, A/A* e inversa (rami subsonico e supersonico), portata per unità di area, condizioni critiche e choking, ugello con urto nel divergente, spinta. Collegamento stazionario attraverso una brusca variazione d'area.
3. **Urto normale stazionario**: stato a valle, rapporti, perdita di pressione totale, salto di entropia.
4. **Urti non stazionari 1D**: urto mobile per entrambe le famiglie, a partire da uno qualsiasi dei dati che compaiono nelle tracce (velocità w, Mach dell'urto, pressione a valle, velocità a valle, Mach assoluto a valle); urto riflesso da parete; stato in moto a monte.
5. **Onde semplici 1D**: invarianti di Riemann, curve di espansione e compressione nel piano p–u, ventaglio centrato campionabile in x/t (incluso il caso transonico), problemi del pistone (velocità impulsiva, arresto, legge di moto qualsiasi xp(t)), tempo e posizione di formazione dell'urto sia da legge del pistone sia da condizioni iniziali di onda semplice.
6. **Problema di Riemann 1D esatto**: implementare **esattamente l'algoritmo del formulario**. Gas diversi ai due lati, classificazione della soluzione con la nomenclatura del corso (es. RCS, SCR), generazione di vuoto, campionamento in x/t. Varianti per le condizioni al contorno tipiche: parete chiusa, pressione imposta (estremo aperto), pistone. Problemi di interazione (urto–contatto, urto–urto, onda–parete, urto attraverso variazione d'area) risolti come pseudo-Riemann, come negli svolgimenti del corso. Se il formulario prevede modelli isotermo o isentropico, aggiungerli come opzione.
7. **Supersonico stazionario 2D**: relazione θ–β–M (rami debole e forte), θmax, deviazione sonica, angolo di Mach, Prandtl–Meyer diretta e inversa, compressione centrata. Polari p–θ e nel piano odografo (u*, v*), con stato di partenza a direzione qualsiasi. Intersezione di due polari (problema di Riemann stazionario 2D, incontro di due correnti, getti con contorno a pressione assegnata). Riflessione regolare e di Mach (teoria dei tre urti), criteri di von Neumann θN e di distacco θD. Teoria shock-expansion per lastre e profili (coefficienti di portanza e resistenza). Pitot supersonico (Rayleigh).
8. **Grafici**: diagramma x–t di una soluzione (urti, contatti, ventagli), piano p–u con curve d'onda e stati etichettati, polari p–θ e odografo con stati numerati, profili spaziali di p, ρ, u, M, p0 a un istante dato. Stile sobrio, adatto a una soluzione d'esame. Salva figure di esempio in `verifica/figure/` e controllale visivamente.
9. **Numerica (priorità bassa)**: un passo dello schema di Godunov 1D con il solutore esatto, per i confronti con la soluzione esatta richiesti da alcune tracce.

Se dalle fonti emergono altri calcoli ricorrenti non elencati qui, proponili.

# Validazione

Tre livelli, tutti in `tests/` ed eseguibili con `pytest`:

1. **Coerenza interna**: funzioni inverse (A/A* e inversa, ν e inversa, θ–β–M), limiti noti (urto debole → isentropico, M = 1, γ → 1 dove sensato), conservazione di massa, quantità di moto ed energia attraverso ogni salto calcolato.
2. **Valori di riferimento classici**, per esempio: urto normale a M = 2 con γ = 1.4 (p2/p1 = 4.5, M2 = 0.5774, p02/p01 = 0.7209); urto obliquo M = 2, θ = 10° (β debole ≈ 39.3°); ν(2) ≈ 26.38°; problema di Sod (p* ≈ 0.3031, u* ≈ 0.9275).
3. **Risultati del corso**: estrai da `fonti/catalogo.md` tutti i risultati numerici riportati negli svolgimenti (campo "Associazione e contenuto dello svolgimento") e ricalcolali con il modulo. Priorità ai 17 esercizi di `soluzioni_A/`, poi agli altri con svolgimento. Esempi da includere: B01-E02 (θD = 19.56°, hR/hL = 2.63672, L/hL = 4.16638), B02-E01 (Me = 2.5141, ΔF ≈ 14.94%), B09-E01 ((SA−SB)/R = 0.264), B41-E03 (u ≈ 1118.2 m/s, p ≈ 25.44 atm), B49-E01 (tBU ≈ 0.666, xBU ≈ 0.788), B51-E01 (θ ≈ −16.5°, p ≈ 0.07 atm), B63-E02, B69-E02, B71-E01, B71-E02, B72-E02.

Per il livello 3, **non dare per scontato che sbagli il modulo o la fonte**. Molti svolgimenti contengono arrotondamenti, ipotesi diverse dalla traccia (per esempio 1 atm = 10⁵ Pa) o errori dichiarati. Per ogni discrepanza cerca la causa, eventualmente controllando la pagina originale nel PDF. I test di livello 3 devono distinguere tra "concorda" e "discorda per causa nota" (marcati come tali, non come fallimenti silenziosi).

# Consegne

1. `gasdyn.py`
2. `tests/test_gasdyn.py`, con tutti i test eseguiti e l'esito riportato.
3. `verifica/verifica_catalogo.md`: tabella con ID esercizio, grandezza, valore dello svolgimento, valore calcolato, differenza relativa, esito (concorda / discorda / non confrontabile) e nota sulla causa delle discrepanze. Alla fine, l'elenco degli svolgimenti che si possono considerare numericamente affidabili.
4. `gasdyn_guida.md`: guida breve per le chat future su claude.ai, con le convenzioni, l'elenco delle funzioni con firma e un esempio d'uso per ciascuna area, e i limiti noti del modulo. Deve bastare a un modello che non ha visto questa sessione per usare il modulo correttamente.

# Fasi

- **Fase 0 (poi fermati e attendi la mia conferma):** leggi formulario e summary e presentami (a) una tabella delle convenzioni (notazione, segni delle velocità e delle famiglie d'onda, orientamento degli assi, riferimenti, nomi degli stati, adimensionalizzazioni, nomenclatura delle configurazioni), (b) i passi dell'algoritmo di Riemann come li hai capiti dal formulario, (c) eventuali discrepanze tra le fonti, (d) la struttura del modulo con le firme delle funzioni principali. Non scrivere codice del modulo prima della mia conferma.
- **Fase 1:** costruisci il modulo per aree, nell'ordine sopra. Per ogni area esegui subito i test di livello 1 e 2 prima di passare alla successiva. Commit e aggiornamento di `CLAUDE.md` a fine area.
- **Fase 2:** validazione di livello 3 e report.
- **Fase 3:** guida d'uso, pulizia finale, riepilogo di cosa è stato fatto e dei limiti noti.

Lavora in italiano (commenti, docstring, report). Sii esplicito su ogni ipotesi che introduci.