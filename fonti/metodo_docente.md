# Metodo del docente (Favini, Gasdinamica): come impostare e svolgere gli esercizi d'esame

## 0. Scopo e uso

Questo documento descrive **come ragiona il docente**: in che ordine affronta un problema, cosa si aspetta in uno
svolgimento e a cosa dà peso. È ricavato da 44 lezioni trascritte (a.a. 2021/22, programma leggermente diverso).
**Non contiene teoria**: formule e derivazioni stanno nel formulario (fonte primaria) e nel riassunto delle
lezioni (seconda fonte). In caso di conflitto su una formula o un segno vale il formulario; questo documento vale per
il metodo.

Legenda: **[D]** = il docente lo dice esplicitamente; **[I]** = inferito dal suo modo di procedere.
Le citazioni tra virgolette sono sue (dalle trascrizioni automatiche, quindi parafrasi fedeli più che verbatim).
I riferimenti `Bxx-Eyy` rimandano al catalogo delle tracce 2003–2023.

---

## (a) Checklist di impostazione di un esercizio

Il docente segue sempre questa sequenza [I, confermata da molte affermazioni D]:

1. **Leggere bene il testo e i dati.** "Il tempo di risoluzione dipende dalla comprensione qualitativa del campo
   d'onda nello spazio-tempo"; leggendo bene si trovano "scorciatoie esatte" [D]. Annotare il gas **di ogni
   regione** (γ, R possono differire ai lati di un contatto o tra serbatoio e condotto).
2. **Scegliere e dichiarare il modello**, con la sua giustificazione e i suoi limiti [D]:
   - Eulero 1D instazionario (condotto a sezione costante, onde piane);
   - omoentropico (serve per gli invarianti di Riemann) vs solo isentropico lungo la traiettoria;
   - quasi-1D (variazione d'area; dà il comportamento **medio** sulla sezione, non separazioni né urti obliqui);
   - 2D supersonico stazionario (valido solo se **tutte** le regioni restano supersoniche);
   - bilancio a parametri concentrati (svuotamento di serbatoi: non le caratteristiche);
   - approssimazioni: urto debole ≈ onda isentropica, onde a spessore nullo, pseudo-Riemann, brusca variazione d'area.
   Se i dati sono incerti, formulare "ipotesi fisiche coerenti" e dichiararle [D].
3. **Analisi qualitativa prima dei numeri**: ai contorni e in ogni discontinuità iniziale chiedersi quali onde
   nascono, di quale famiglia (u−a, u, u+a; C⁺, C⁻), di quale natura (urto, espansione, contatto), in che verso
   spingono il fluido (§d1–d3).
4. **Topologia del piano x–t** (o x–y): zone uniformi, onde semplici, zone d'interazione; fronti anteriori e
   posteriori; quali caratteristiche sono semirette (§d2).
5. **Contare equazioni e incognite** prima di risolvere (es. Riemann SCS: 8 = 3+3 R–H + 2 di contatto) [D].
6. **Cercare la scorciatoia fisica**: invariante uniforme (onda semplice), condizione critica a un imbocco o a
   un'uscita, uscita supersonica che non riceve informazione, caso limite 1D come stima [D].
7. **Calcolare** con formulario, tabelle (interpolate) e solutori.
8. **Verificare a posteriori** ogni ipotesi (§c, punto 6): natura delle onde, Mach relativo, entropia, a > 0, BC
   compatibili, ramo debole/forte, regioni supersoniche.
9. **Rappresentare**: x–t sempre; p–u per Riemann e interazioni 1D; p–θ (e odografo) per urti obliqui e
   interazioni stazionarie.
10. **Discutere la fisica e confrontare i modelli** (isentropico vs urto, 1D vs 2D, stazionario vs instazionario).

---

## (b) Convenzioni

### Coincidenti con il formulario [D]
- δ = (γ−1)/2; s ~ p/ρ^γ (variabile entropica), ds = dp − a² dρ.
- **R1 = a − δu costante lungo u−a; R2 = a + δu costante lungo u+a.**
- Riemann: **1 sinistra, 4 destra**; 2 a valle dell'onda u−a, 3 a valle dell'onda u+a; p2 = p3, u2 = u3; ρ2 ≠ ρ3.
  Configurazioni RCR, RCS, SCR, SCS (onda sinistra–contatto–onda destra).
- Urti mobili: 1 = monte, 2 = valle; Mach **relativo** M_s = |u1 − W|/a1.
- Solutore di Riemann: primo tentativo isentropico su **u** e iterazione su u (secante o Newton).
- 2D: θ positivo antiorario, θ + ν costante lungo C⁻ (dy/dx = tan(θ−μ)), θ − ν lungo C⁺; angolo d'urto σ;
  stati delle riflessioni 1 (monte), 2 (dopo l'incidente), 3 (dopo il riflesso), 4 (dietro il Mach stem).
- Getti: adattato p_e = p_a, sovraespanso p_e < p_a, sottoespanso p_e > p_a (parametro n = p_e/p_a).

### Differenze solo di etichetta
| Docente | Formulario | Nota |
|---|---|---|
| angolo di Mach **α** | **μ** | stesso significato |
| Z = (p4/p1)^((γ−1)/2γ) a1/a4, u⁽⁰⁾ = [(a1+δu1) − Z(a4−δu4)]/[δ(1+Z)] | z = 1/Z | formula equivalente |
| onde numerate 1–4 per autovalore crescente (u−a, u, u, u+a) | λ1 = u−a, λ0 = u, λ2 = u+a | in 1D con velocità trasversale ci sono 4 onde |
| curve p–u **RCN, SCN** (uscenti da 1), **NCR, NCS** (entranti in 4) | curve u_L(p), u_R(p) | N = onda nulla |
| salto [φ] = φ⁺ − φ⁻, media φ̄ = (φ_L+φ_R)/2 | — | |
| brusca variazione d'area: 4 quiete, 2 dietro l'urto, 1 dopo il salto d'area | summary: 1, 3, 5 | etichette libere |

### Definizioni a cui tiene [D]
- **Monte/valle**: non il tempo, non il verso di x, non il segno dell'autovalore, ma l'attraversamento della
  caratteristica da parte della particella (monte = prima di essere investita). Per **onde di entropia e di
  scorrimento non esistono monte e valle** (solo destra/sinistra).
- **Isentropico** (Ds/Dt = 0 lungo la traiettoria) ≠ **omoentropico** (s uniforme ovunque).
- **Transonico** = corrente lontana subsonica con bolle supersoniche; per un corpo tozzo in corrente supersonica si
  parla di **bolla subsonica**.
- Glossario: onda di entropia ("**bolla termica**"), onda di scorrimento, impedenza acustica ρa, urto "veloce"
  (u+a) e "lento" (u−a), primary/secondary shock, **gambetto di Mach** (Mach stem), punto triplo, **slipstream** o
  linea di scorrimento, barrel shock, disco di Mach, shock diamonds, potential core, linea sonica, stand-off
  distance, **poppa fluida**, lip shock, "**curva del cuore**" (piano p–θ).

---

## (c) Cosa deve contenere uno svolgimento completo

1. **Ipotesi e modello dichiarati** con la giustificazione (punto a.2) e il campo di validità.
2. **Descrizione qualitativa** delle onde generate, con famiglia e natura, prima del calcolo.
3. **Piano x–t** in scala qualitativa con stati numerati, fronti di testa e coda, contatti, urti; per interazioni
   successive, la sequenza degli eventi con tempi e posizioni d'incontro.
4. **Piano p–u** (1D) o **p–θ** (2D) con le curve e le intersezioni; ricordare che il contatto **non si vede** in p–u
   e che due stati coincidenti in p–u possono avere impedenze diverse (non si legge "a vista") [D].
5. **Calcoli ordinati** stato per stato, con il valore di γ, R di ciascuna regione.
6. **Verifiche** [D]:
   - natura delle onde coerente con il risultato (sinistra R se u* > u1, S altrimenti; destra R se u* < u4);
   - M relativo a monte > 1 e Δs > 0 per ogni urto ("gli urti di espansione sono fisicamente impossibili");
   - a > 0, p > 0 (vuoto/cavitazione: |V_p| ≤ a0/δ);
   - pendenza della coda di un'espansione al contorno (se la coda "esce" dal dominio il contorno è **sonico**);
   - numero di condizioni al contorno = caratteristiche entranti;
   - ramo debole/forte giustificato; θ ≤ θ_max; regioni supersoniche nei modelli stazionari;
   - all'esame un residuo ~10⁻³ nel solutore di Riemann "è fisicamente zero".
7. **Discussione fisica** dei risultati e, quando possibile, **confronto tra modelli**.
8. Scrivere **tutte le deduzioni parziali**: "le prove sono calibrate per testare fin dove ciascuno riesce a
   spingersi"; non rinunciare a scrivere perché la soluzione non è conclusa [D].
9. Redazione: "l'introduzione è l'ultima cosa che si scrive", deve riassumere ciò che è stato effettivamente
   dimostrato [D].
10. Avvertenze grafiche [D]: la polare p–θ **non parte dall'origine** (p1 > 0, al più p/p1 = 1); una scala
    logaritmica in p deforma la concavità del ramo debole.

11. **Livello di dettaglio** [D, dalla soluzione ufficiale del 14/09/2026, es. 1]: gli stati uniformi si calcolano
    stato per stato con le due equazioni disponibili (invariante dell'onda + condizione al contorno) e si
    confrontano con il criterio fisico richiesto (es. p_tot di ogni stato con p_tot,lim). **Le zone d'interazione a
    due famiglie non si risolvono**: basta dire che si potrebbero risolvere (compatibilità + caratteristiche +
    condizione al contorno, soluzione non analitica) e qual è la conseguenza qualitativa (es. "la pressione a
    parete è decrescente"). I tempi di arrivo si stimano trattando le caratteristiche come **semirette** anche
    nella zona a due famiglie. Non introdurre parametri ausiliari (es. r = (1−δM)/(1+δM)) se non strettamente
    necessari.

All'esame sono ammessi tabelle, righello, **goniometro** e calcolatrice [D]: le misure di angoli sulle figure
(onde di Mach, urti) fanno parte del lavoro atteso.

---

## (d) Indicazioni per tipologia

### d1. Onde 1D, problema di Riemann, tubo d'urto
- Da un punto singolare dello spazio-tempo possono nascere solo urti, contatti, espansioni centrate: ogni
  singolarità (setto, pistone impulsivo, centro di compressione, incontro di onde) è un problema di Riemann [D].
- **Regole di verso** [D]: un'onda di compressione spinge il fluido nel verso in cui si propaga; un'espansione lo
  richiama nel verso opposto. Nel tubo d'urto (p1 > p4): espansione u−a nel driver, contatto, urto u+a nel driven (RCS).
- Solutore (sua procedura): primo tentativo u⁽⁰⁾ isentropico ("perché la velocità e non la pressione? a pari Δu le
  curve isentropica e R–H differiscono poco in u") → natura delle onde → pressioni → |p2 − p3| → aggiornamento.
- Approssimazione isentropica degli urti: lecita per urti deboli (**M_s ≲ 1,2**) perché Δs ∝ (M_s² − 1)³ e le curve
  R–H e isentropica hanno stessa tangente e curvatura nello stato iniziale [D].
- Un urto **può risalire una corrente supersonica** (conta il Mach relativo); un'espansione **mai** [D].
- **Esempio svolto, B56-E01** (contatto che raggiunge un urto normale stazionario, A: p=1, ρ=1,4, u=2; B: ρ=0,082):
  stato C dietro l'urto (M=2 → p=4,5, u≈0,75); all'impatto Riemann tra B (sinistra) e C (destra); osservazione
  chiave: la bolla calda ha a≈4,13 ⇒ M_B≈0,48, l'urto stazionario non può più esistere; 1D valido solo per t=0⁺
  (lungo la mezzeria nascono onde oblique e punti tripli); p–u non leggibile a vista (ρa = 1,4 vs 0,34).

### d2. Pistoni, onde semplici, formazione dell'urto, interazioni
- **Ordine di analisi di un campo assegnato** [D]: (1) M = u/a **con segno**; (2) entropia p/ρ^γ (omoentropico?
  con dati misurati: mediana e dispersione, pochi % ok, 10–15% no); (3) quali combinazioni caratteristiche sono
  attive (riconoscere le onde dagli invarianti: stesso a+δu ai due lati ⇒ onda u−a).
- **Topologia** [D]: zone a 0 famiglie (uniformi), 1 famiglia (onde semplici, caratteristiche della famiglia attiva
  **semirette**), 2 famiglie (interazione, regione "romboidale"); confine 0↔1 semiretta, 1↔2 curva, 0↔2 un punto.
- Attraversare un ventaglio **lungo una caratteristica dell'altra famiglia** (es. pistone che arretra + parete di
  fondo: A→B con R1, B→C con R2 e u_C = 0 ⇒ a_C = a0 − 2δ|V_p|).
- Espansione da pistone impulsivo: regimi subsonico / limite sonico (coda verticale) / transonico (coda verso
  sinistra); "in instazionario si arriva a M=1 in un condotto a sezione costante" [D]. Limite |V_p| = a0/δ (vuoto),
  maggiore della velocità limite stazionaria a0/√δ.
- **Pistone con legge non lineare** (esercitazione, tipo **B53-E14**, B20-E01, B33-E01, B49-E01, B69-E02):
  classificare il problema misto, controllare la cavitazione, notare se c'è singolarità in (0,0) (ẋ_p(0)=0 ⇒ no
  fronte posteriore), R1 uniforme ⇒ a = a0 + δu, tracciare a ritroso la semiretta u+a fino al pistone e risolvere per il
  **tempo di emissione** t_p (scegliere la radice in [0, t̄]); discutere le discontinuità di gradiente sui fronti.
  Se t_p non è invertibile (legge cubica) si costruisce la soluzione per punti [D].
- **Formazione dell'urto**: compressione centrata (urto finito subito) vs coalescenza (zona vietata con cuspide);
  "la compressione degenera, l'espansione vive per l'eternità" [D]. t_bu = 1/max|V0'|. Al centro di una compressione
  centrata nascono urto + contatto + espansione (la curva R–H sta sotto l'isentropica) [D].
  Esempio **B65-E01** (14/09/2021): espansione transonica in x=0 riconosciuta da a+δu uguale ai due lati, compressione
  centrata u−a nel tratto lineare, analisi a t_bu+ε.
- **Problema di Friedrichs e pseudo-Riemann** [D]: pistone avviato (urto u+a) e arrestato ⇒ espansione **della stessa
  famiglia** che insegue e indebolisce l'urto. Stato finale con pseudo-RP tra C (gas fermo sul pistone, p_C > p_A
  perché l'espansione è più "efficiente" dell'urto) e A (quiete): RCS. Il pseudo-RP dà gli stati asintotici, non la
  struttura interna della zona d'interazione (formulario §9). Tipi: B05-E01, B52-E03, B57-E01, B71-E02.
- Riflessioni [D]: su **parete** un'onda si riflette con la **stessa natura**; su un contorno a **pressione costante**
  (estremo aperto, confine di un getto) con **natura opposta**. Più in generale ragionare con l'impedenza ρa.

### d3. Condizioni al contorno e serbatoi
- Regole [D]: si impongono tante condizioni quante sono le caratteristiche **entranti**; il numero cambia nel
  tempo; vietato sovraspecificare; "**divieto assoluto di prescrivere variabili d'onda uscenti**"; niente
  velocità infinite.
- Estremo aperto: in **uscita** si impone la pressione statica ambiente; in **ingresso** dall'ambiente o da un serbatoio
  si impongono **p0 e T0** (a² + δu² = a0²), non la statica [D]. Parete: u = 0.
- **Verificare sempre il contorno a posteriori**: se l'espansione verso l'ambiente ha la coda con pendenza uscente,
  l'uscita è **strozzata** (u = a); se la pressione critica dell'imbocco p0(2/(γ+1))^(γ/(γ−1)) supera la pressione
  nel condotto, l'imbocco è **sonico** e lo stato è noto senza iterare [D].
- **Esempio svolto, esame 24/01/2019** (non nel catalogo; affini B37-E01, B11-E01, B27-E02): condotto tra serbatoio
  (γ=1,148) e ambiente. A destra: espansione u−a, coda con u−a = +0,7 > 0 ⇒ uscita sonica. A sinistra: il gas del
  serbatoio fa da "pistone aerodinamico" (urto + contatto); 7 incognite = 3 R–H + 2 contatto + 2 BC; scorciatoia:
  imbocco critico (p* ≈ 14,9 bar ≫ 3,6 bar) ⇒ p_B, u_B noti ⇒ M_s ≈ 1,91.
- **Svuotamento di un serbatoio, B66-E01**: la portata adimensionale costante 0,578 = √γ(2/(γ+1))^((γ+1)/(2(γ−1)))
  rivela l'efflusso critico; poi "ginocchio" sotto la pressione critica. Modello: **bilancio integrale a parametri
  concentrati** d(ρV)/dt = −ṁ, non caratteristiche [D].
- **Esempio svolto, esame 14/09/2026 es. 1** (non nel catalogo; soluzione ufficiale): tubo L = 10 m chiuso a
  sinistra, a destra brusco restringimento A2/A1 = 0,506 con setto; p0 = 10 atm, T0 = 293,15 K, pa = 1 atm;
  pressione sulla parete di testa fino alla fine dell'efflusso critico. Impostazione: ipotesi diametro ≪ L ⇒ onde
  piane, modello 1D nel solo tubo; contorno destro M_e = 0,31 (legge delle aree con A2 critica) finché l'uscita è
  strozzata, sinistro u = 0. Stati: 1 (u−a dall'apertura: a1+δu1 = a0, u1 = M_e a1 ⇒ p1 = 6,56 atm, p1,tot = 7,01),
  2 (u+a riflessa a parete, u2 = 0 ⇒ 4,2 atm), 3 (u−a dall'uscita ⇒ 2,76 atm, p3,tot = 2,95), 4 (seconda
  riflessione a parete ⇒ 1,76 atm < p_tot,lim = 1/0,528 = 1,892 atm) ⇒ l'uscita diventa adattata durante la seconda
  riflessione all'uscita e lo studio termina. Tempi a parete con caratteristiche semirette: 0,029–0,045 s e
  0,085–0,127 s; nelle zone di riflessione solo "pressione decrescente". Grafico p–t a gradini e piano x–t.
- "Gola sonica ⇒ portata costante?" No: ṁ ∝ p0/√T0 varia se il serbatoio si svuota [D].

### d4. Quasi-1D e condotti
- Il quasi-1D è "il modello più semplice per il comportamento medio"; cattura posizione e velocità degli urti
  d'avviamento e di ricompressione, non i dettagli 2D [D].
- Onde in sezione variabile [D]: in un **divergente** (u>0) l'onda si **attenua**, in un **convergente** si
  **intensifica**; **non esistono onde semplici** (entrambe le famiglie sempre attive); una **bolla termica che
  attraversa un ugello genera onde acustiche** ("i motori suonano": instabilità termoacustiche).
- **Brusca variazione d'area** (stima dell'avviamento senza CFD) [D]: R–H per lo stato dietro l'urto → se è
  supersonico anche in assoluto, l'attraversamento del salto d'area è istantaneo, stazionario e isentropico (legge
  delle aree) → problema di Riemann tra lo stato espanso e la quiete. Affini (contrazione invece di allargamento): B72-E02, B08-E02.
- Divergente supersonico: niente blocco; convergente: chiedersi sempre se A_min < A* (urto che risale) [D].
- Coerenza tra modelli: quasi-1D e 2D devono dare lo stesso risultato quando le pareti non generano riflessioni
  (vedi aerospike, d6) [D].

### d5. Urti obliqui, polari, riflessioni
- Analogia completa x–t ↔ x–y: "ogni sistema iperbolico 2×2 ammette invarianti di Riemann" [D]; la compressione su
  parete concava degenera in urto obliquo come quella del pistone in x–t.
- **Debole o forte** [D]: la debole lascia valle supersonica (iperbolico, l'analisi dello spigolo è autoconsistente);
  la forte richiede una contropressione imposta a valle. **Flussi esterni: sempre la debole**; flussi interni (prese,
  condotti, ugelli) possono stare sulla forte. θ > θ_max ⇒ urto staccato.
- Diedro: superiore e inferiore **disaccoppiati**; a incidenza, sul lato che diverge un ventaglio di Prandtl–Meyer.
- **Confronto quantitativo** che ama (M1 = 2, θ = 20°): compressione isentropica M2 ≈ 1,31, p2/p1 ≈ 2,79; urto debole
  M2 ≈ 1,32, p2/p1 ≈ 2,80 ⇒ equivalenti finché M1 ≲ 3 [D, soglia orientativa]; oltre, l'isentropica sta sopra la polare.
- **Riflessione su parete o asse** (procedura p–θ) [D]: polare di M1 → stato 2 sul ramo debole → polare dello stato 2
  (più piccola) percorsa verso θ = 0 → intersezione = **RR**. θ_N (von Neumann): la seconda polare incontra la prima
  su θ = 0. Oltre θ_N: punto triplo staccato e Mach stem **curvo**. θ_D = θ_max(M2): oltre, RR impossibile.
  Tra θ_N e θ_D **isteresi** (crescendo resta RR fino a θ_D, decrescendo resta MR fino a θ_N). L'angolo del riflesso non
  è uguale a quello dell'incidente ("non vale la legge ottica").
- **Urto piano mobile su rampa** (Schardin, Album n. 235–236; catalogo B30-E01, B70-E01): struttura autosimile;
  riferimento solidale al piede, V_T = W/cos φ, M_T = (W/a0)/cos φ > 1, σ = 90° − φ ⇒ teoria stazionaria. Rampa ripida ⇒
  RR, poco inclinata ⇒ MR [D].
- Urto curvo ⇒ salti diversi su ogni linea di corrente ⇒ campo rotazionale (Crocco).

### d6. Getti e interazioni stazionarie
- Confine del getto = **strato di scorrimento**, condizione **p = p_a** (non geometrica): "sul contorno del getto non
  esiste alcuno strato limite" [D].
- Sottoespanso: ventagli allo spigolo → riflessione sul contorno isobaro **come compressioni** → coalescenza in urto
  (barrel shock) → riflessione sull'asse RR o MR (disco di Mach, punto triplo, slipstream) [D].
- **Riemann stazionario supersonico** (urto su contatto inclinato, getti affiancati) [D]: sulla linea di scorrimento
  **stessa p e stessa θ**; dallo stato ad alta p parte un'espansione, da quello a bassa p un urto; iterare su θ.
  Condizione: tutte le regioni supersoniche, altrimenti il modello cade. Tipi: B33-E02, B74-E01, B56-E01 (domanda 2).
- **Rifrazione dell'urto** su un'interfaccia: urto trasmesso + onda riflessa + nuovo contatto; partire dal caso piano
  inclinato, trascurare le pareti nella finestra utile; il caso **1D a incidenza normale** si risolve in pochi minuti
  ed è spesso il **caso peggiore** per i carichi [D]. Natura della riflessa: impedenza crescente ⇒ urto, decrescente ⇒
  espansione.
- **Esempio svolto, B59-E02 (aerospike)**: onda semplice centrata nello spigolo ⇒ θ0 = ν(3) = 49,76°; testa del
  ventaglio ortogonale alla gola, coda a α = 19,47°; parete **tangente punto per punto alla velocità** dell'onda semplice
  (non riflettente); h_e/h0 = A/A*(3) = 4,234; L/h0 = 4,234/tan 19,47° ≈ 12; p_a/p_c = 0,0272. Conclusione: coerenza tra
  modello 2D e quasi-1D [D].
- **Esempio svolto, B64-E04 (getto impattante su parete)**: descrizione per numeri (vedi d8); picco di pressione al
  ristagno = urto normale (disco di Mach) + compressione isentropica subsonica; lungo la parete la lama accelera, scende
  sotto l'ambiente e si ricomprime con un urto radiale (andamento non monotono). Attenzione: p_e è la pressione di
  uscita dell'ugello, non l'ambiente. Topologia 3D assialsimmetrica = 2D con il termine 1/r [D].

### d7. Corpi tozzi e transonico
- Corpo tozzo = impone θ > θ_max; inoltre l'urto da solo non azzera la velocità (u2/u1 ≥ (γ−1)/(γ+1)) ⇒ urto staccato
  a distanza finita e decelerazione subsonica fino al ristagno [D] (è il "paradosso del ristagno" risolto dal quasi-1D).
- Lungo il bow shock il campo a valle percorre **tutta la polare**; la linea sonica **non è un'isobara** (p0 diversa su
  ogni linea di corrente); campo rotazionale (Crocco) [D].
- "Efficienza vs efficacia" (H. J. Allen): l'urto dissipa ma riduce i carichi termici ⇒ capsule tozze [D].
- Transonico: M∞ critico (cilindro ≈ 0,45, lamina → 1); bolla supersonica chiusa da urto di ricompressione (Morawetz:
  il flusso senza urti è eccezionale); buffeting per retroazione scia–acustica–urto [D].
- Linearizzazione (Prandtl–Glauert, Ackeret) fallisce per M → 1 ⇒ TSD; in gola (Sauer–Meyer) la linea sonica è una
  parabola, **la gola geometrica non coincide con la linea sonica** [D]. Formule: formulario §19 (vedi f3).

### d8. Interpretazione di immagini (schlieren, shadowgraph, interferogrammi, mappe CFD)
Atteggiamento [D]: le immagini vanno studiate insieme alle equazioni (Van Dyke, *An Album of Fluid Motion*), ma
"raccogliere immagini senza una chiave interpretativa è come collezionare santini"; l'immagine si interroga con la
teoria.

**Protocollo** [I, dagli esempi]:
1. Stazionario o instazionario? (onde sferiche concentriche = sorgenti impulsive; strutture che cambiano tra
   fotogrammi; un'immagine **mediata nel tempo** cancella le discontinuità mobili).
2. Riferimento e verso della corrente; regime (dove ci sono onde di Mach rettilinee il campo è supersonico).
3. Elencare le strutture **nell'ordine dell'immagine o dei numeri del testo**, nominando per ciascuna tipo e
   famiglia: urto (linea netta e sottilissima), contatto/strato di scorrimento (sfumato, che si allarga, vortici),
   ventagli, onde di Mach, punti tripli, dischi di Mach, scie, zone di separazione.
4. **Giustificare l'aspetto con la fisica**: l'urto resta sottile perché le caratteristiche vi convergono e lo
   "richiudono"; il contatto si allarga perché è linearmente degenere [D].
5. **Quantificare**: Mach dall'angolo delle onde, sin α = 1/M, misurando con il goniometro (proiettile: α ≈ 27°–30° ⇒
   M ≈ 2–2,2); velocità degli urti da posizioni a tempi noti; regime dal segno dei gradienti con la legge delle aree
   (in un divergente p che cresce verso valle ⇒ subsonico) [D].
6. Separare il **fisico dall'artefatto**: le frange di un interferogramma sono curve di livello quantizzate dalla
   soglia dello strumento (come le medie di cella di Godunov); un salto di pressione vicino alla gola in una
   simulazione quasi-1D è la media sull'asse di urti obliqui a parete, non un errore [D].

Esempi commentati a lezione:
- **Avviamento di un ugello (interferogramma Saitō–Takayama)**: primary shock (u+a, veloce, sottile), interface
  (contatto, sfumato), secondary shock (u−a, lento, di ricompressione del flusso sovraespanso nel divergente; in realtà
  una riflessione di Mach con punti tripli). Nei grafici di pressione: il salto del primario cala nel tempo (divergente),
  il contatto non appare nella pressione, la distanza tra i due urti cresce.
- **Proiettile attraverso un cilindro forato (schlieren)**: scia turbolenta dietro la poppa tronca; onde sferiche dalle
  asole; cono di Mach come inviluppo ⇒ stima di M (cfr. B55-E03, B57-E02).
- **Buffeting** (Album): urti la cui posizione dipende dall'istante; immagine mediata senza discontinuità.
- **Getto impattante** (B64-E04): bocca, regione sottoespansa, onde deboli dall'ugello, disco di Mach, urto riflesso
  dal punto triplo, shear layer, regione subsonica fino al ristagno, lama radiale con urto di ricompressione.
- Tracciato di un trasduttore in un tubo d'urto: salto seguito subito da un calo senza plateau ⇒ l'espansione riflessa
  ha già raggiunto l'urto [D].

### d9. Numerica
- Shock fitting (preciso, topologicamente ingestibile in 2D/3D) vs shock capturing (forma integrale) [D].
- **Godunov a mano** (procedura del suo esempio, contatto puro ρ_L = 1,4, ρ_R = 0,35, u = p = 1):
  a_L, a_R → λ_max = max(|u|+a) → Δt/Δx = α/λ_max (α ≤ 1) → RP a ogni interfaccia (qui u⁽⁰⁾ = u1 = u4: solo il contatto,
  W = 1 > 0, l'asse t cade nello stato 1) → flussi F(U^RP) → aggiornamento delle celle (1,4 | 0,70 | 0,35) →
  verificare che u e p restino uniformi.
- Interpretazione [D]: con λ > 0 Godunov = upwind al 1° ordine ("generalizzazione non lineare delle caratteristiche
  upwind"); media convessa ⇒ diffusione numerica ν_num = (uΔx/2)(1 − uΔt/Δx); il **contatto** (linearmente degenere) si
  allarga indefinitamente, l'**urto** (genuinamente non lineare) resta su 2–3 celle.
- Quasi-1D con Godunov: flussi dal RP a sezione costante + sorgente di cella −p̄ dA/dx [D].
- Condizioni al contorno nei metodi numerici: stesse regole delle caratteristiche entranti (d3).

---

## (e) Errori da evitare

Dichiarati dal docente [D]:
1. "Divergenza nulla ⇒ incomprimibile ⇒ densità costante": no, la densità resta costante solo lungo la traiettoria.
2. "Gola sonica ⇒ portata costante": no, dipende da p0/√T0.
3. Usare p0 e h0 costanti in problemi **instazionari** ("il supersonico instazionario è profondamente diverso dallo
   stazionario").
4. Confondere il regime della corrente con quello dell'urto: conta il Mach **relativo**; un gas in quiete può essere
   urtato, un urto può risalire una corrente supersonica.
5. Ammettere urti di espansione.
6. Prescrivere al contorno una variabile d'onda uscente, o sovraspecificare; imporre la pressione statica in ingresso
   da un serbatoio (va p0, T0).
7. Parlare di strato limite al bordo di un getto.
8. Applicare "formulette" senza aver capito il campo d'onda.
9. Prolungare le caratteristiche oltre il blow-up (soluzione a più valori).
10. Leggere a vista il piano p–u quando gli stati hanno impedenze diverse.
11. Disegnare la polare p–θ partendo dall'origine.
12. Usare il modello stazionario 2D quando una regione diventa subsonica.

Inferiti dai suoi controlli [I]:
13. Non verificare la natura delle onde ipotizzate (R/S) dopo il calcolo.
14. Usare un solo γ o R quando ai lati di un contatto ci sono gas diversi.
15. Dare per scontata la condizione al contorno senza controllarla (uscita o imbocco sonico).
16. Confondere le normalizzazioni (p_e di uscita vs p_a ambiente; grandezze iniziali vs di ristagno).
17. Dimenticare il limite di vuoto (a ≥ 0) in espansioni forti e nei pistoni che arretrano.

---

## (f) Differenze rispetto al tuo corso

### f1. Presenti solo nelle trascrizioni
- Teoria delle caratteristiche **2D instazionaria** completa: cono di Mach nello spazio-tempo, bicaratteristiche,
  **teorema di Butler**, diagonalizzazione separata per assi; salti delle derivate (autovettori destri) e
  Cauchy–Kowalevski in forma estesa; splitting di Steger–Warming esplicito.
- **Formula di Kentzer** e salti in variabili primitive con coefficienti mediati.
- Shock fitting vs shock capturing; equazione modificata e viscosità numerica; LES/ILES e viscosità turbolenta.
- Onde in condotti a sezione variabile come fenomeno (attenuazione/intensificazione, rumore da entropia).
- Interpretazione dettagliata di interferogrammi e schlieren; buffeting, M∞ critico, Morawetz, condensazione.
- Riflessione **instazionaria** di un urto su rampa (Schardin), **diffrazione**, **Richtmyer–Meshkov**, **rifrazione**
  dell'urto su interfacce e bolle (queste ultime compaiono nel catalogo: B30-E01, B70-E01, B59-E03).
- Parametro di similitudine transonica K; getto con ago sull'asse; getto impattante (B64-E04).

### f2. Presenti solo nel tuo corso (summary)
- Problema di Friedrichs **inverso** e l'esempio del tubo d'urto chiuso; distribuzione **cubica simmetrica** per la
  formazione dell'urto.
- Condizione di serbatoio con il criterio AL/V ≪ 1 (il docente qui usa solo p0, T0 costanti all'imbocco).
- Ghost cells per le condizioni al contorno in Godunov; il summary osserva che Godunov produce entropia.
- Il formulario non ha una sezione di numerica; le trascrizioni sì.
- Conteggio delle condizioni al contorno: il formulario usa Eulero 1D a 3 equazioni (3/2/1/0 condizioni); il docente,
  includendo la velocità trasversale, conta 4/3/…/0. Stessa regola, modello diverso.

### f3. Discrepanze ed errori probabili (non seguire le trascrizioni in questi punti)
- **Sauer–Meyer** (lezione 42, chiamata "Meyer-Sänger"): la trascrizione azzera f4 per un errore di derivazione e
  conclude che la linea v = 0 è ξ = 0; il formulario ha il termine cubico in v e la gola ξ = −(d/6)η². **Seguire il
  formulario** (il modulo ha anche la variante 'favini').
- **TSD**: la scalatura della trascrizione assorbe (γ+1) e poi lo fa ricomparire; usare la forma del formulario
  [k − (γ+1)φ_ξ]φ_ξξ + φ_ηη = 0.
- **Rifrazione dell'urto**: le etichette slow–fast/fast–slow risultano invertite rispetto alla letteratura; usare la
  regola delle impedenze.
- **B65-E01**: la trascrizione dà a = 1 + 0,3711x (coerente con a+δu costante ma con t_bu ≈ 0,449), il catalogo
  a = 1 + 0,4995x (coerente con t_bu = 0,4246 ma non con un'onda semplice). Discrepanza non risolta: controllare il testo.
- In due lezioni (Cauchy–K.) la trascrizione scrive "stato 1 a destra": incoerente con la sequenza descritta subito
  dopo; vale 1 a sinistra.
- "Singolarità di Richardson" al piede dell'urto transonico (lez. 40): altrove il docente la chiama di **Zierep**.
- Condizione di vuoto citata a voce come "u_lim,1 + u_lim,4": vale quella del formulario u4 − u1 ≥ (a1+a4)/δ.

---

## Appendice — esempi svolti a lezione e corrispondenza con il catalogo

| Esempio a lezione | Catalogo | Sezione |
|---|---|---|
| Contatto che raggiunge un urto stazionario (esame 16/03/2018) | B56-E01 | d1 |
| Pistone con legge parabolica (esercitazione) | B53-E14; affini B20-E01, B33-E01, B49-E01, B69-E02 | d2 |
| Dati iniziali isentropici, espansione + compressione centrata (14/09/2021) | B65-E01 | d2 |
| Pistone avviato e arrestato, pseudo-Riemann (Friedrichs) | B05-E01, B52-E03, B57-E01, B71-E02 | d2 |
| Condotto tra serbatoio e ambiente (24/01/2019) | non in catalogo; affini B37-E01, B11-E01 | d3 |
| Svuotamento di un serbatoio dal grafico | B66-E01 | d3 |
| Svuotamento di un tubo da un orifizio, pressione in testa (14/09/2026, soluzione ufficiale) | non in catalogo | d3, c11 |
| Brusca variazione d'area (avviamento) | affini (contrazione): B72-E02, B08-E02 | d4 |
| Rampa: isentropica vs urto obliquo (M=2, θ=20°) | affini: B06-E02, B23-E01, B27-E03 | d5 |
| Urto piano su rampa, RR/MR | B30-E01, B70-E01 | d5 |
| Ugello aerospike con onda semplice centrata | B59-E02; immagini B72-E01, B69-E03 | d6 |
| Getto sottoespanso impattante | B64-E04; affini B42-E03 | d6, d8 |
| Proiettile, stima del Mach dalla schlieren | B55-E03, B57-E02 | d8 |
| Godunov a mano su un contatto | B44-E03, B38-E02, B14-E02 | d9 |
