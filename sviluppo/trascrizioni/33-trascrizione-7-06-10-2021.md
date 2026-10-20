# Trascrizione integrale della lezione del 6 ottobre 2021
## Aerotermodinamica e Propagazione Ondosa: Cono di Mach, Bicaratteristiche, Piani Acustici e Impedenza Acustica

---

### 1. Considerazioni preliminari sul versore normale $\mathbf{\nu}$ e questioni di notazione

Forse è ora di cominciare. C'era una domanda al volo:
«Non è molto volante...»
Vabbè, la domanda può essere volante, la risposta di solito non lo è; se è molto corta è un problema.
«Può darsi che sia una questione di notazione: qui il versore $\mathbf{\nu}$ l'avevamo preso sullo spazio $(x, y)$, sostanzialmente sul piano...»

Questo versore $\mathbf{\nu}$ è una cosa che ci crea qualche difficoltà: una volta lo prendiamo direttamente sul piano bidimensionale, adesso invece lo collochiamo nello spazio $(x, y, t)$. Diciamo che uno può giostrarsela usando due simboli, aggiungendo un segno distintivo in modo tale da memorizzare se stiamo parlando di un oggetto che vive solo in $(x, y)$ oppure di un oggetto che vive anche nel tempo. Però lui è planare: quando parliamo di $\mathbf{\nu}$ parliamo sempre di una grandezza planare. Solo che ci può far comodo, come vedrete adesso, una volta guardarlo nello spazio tridimensionale esteso e una volta guardarlo nel sottospazio fisico. Per cui aggiungeremo un cappelletto, $\hat{\mathbf{\nu}}$, per distinguere chiaramente la notazione, anche se dal contesto si capisce. Ma non è un vero problema concettuale: non serve attaccarsi alla notazione per risolvere i problemi di comprensione della fisica.

---

### 2. Riflessione epistemologica: la gasdinamica, i modelli continui e le metafore

Conoscete la gasdinamica? La risposta più onesta, senza voler togliere nulla a nessuno, è no. Non escludo che conosciate già dei pezzi, non escludo che qualcuno già la conosca abbastanza bene ma debba ancora superare l'esame; tuttavia, è ragionevole dire che normalmente non conosciamo la gasdinamica. 

Abbiamo dei riferimenti di vita quotidiana: mettetevi nei panni di una persona normale, che non sia nata con la fissa per le bombe, le esplosioni, i meteoriti, il vento solare, i lanciatori o i missili. Uno studente che decide di iscriversi a ingegneria aerospaziale conosce la gasdinamica? No, non è ragionevole pensarlo. Qualche fenomeno lo si è visto: il getto d'acqua che colpisce un piatto, il suono che si propaga, una sirena, le onde del mare. Sono tutti appigli su cui costruire analogie e metafore. 

È giusto farlo? **Assolutamente sì**. Il nostro modo di ragionare procede in primo luogo per analogie e metafore. Ma sono tutte buone? **No, non sono tutte buone**. E come facciamo a distinguere quelle buone da quelle cattive? C'è un unico modo, per fortuna o per sfortuna: **equazioni alla mano**. O quelle analogie riusciamo a tradurle rigorosamente nei bilanci e nelle equazioni differenziali derivate, con le opportune ipotesi e semplificazioni, oppure non ce ne facciamo nulla. Una metafora può sembrare bellissima, ma se non si traduce in formule operative resta priva di utilità scientifica.

Dobbiamo armarci di santa pazienza per imparare le equazioni. Non per impararle in sé come sterile esercizio matematico, ma perché descrivono fedelmente i fenomeni fisici. Manipolare le equazioni, «farle parlare», significa smontarle e rimontarle. La teoria delle caratteristiche per un sistema multidimensionale instazionario la stiamo affrontando per la prima volta. Dobbiamo familiarizzare con concetti precisi: variabile d'onda, direzione di propagazione, condizioni a monte e a valle. Ieri qualcuno chiedeva: «Ma non stiamo ricavando l'equazione di d'Alembert?». L'equazione delle onde di d'Alembert è solo un sottocaso particolare, acustico e lineare, di una formulazione molto più generale. La propagazione di un'onda d'urto non è la propagazione del suono: un'onda d'urto smuove l'intero campo termofluidodinamico, non induce una semplice fluttuazione infinitesima di pressione.

Noi ci muoviamo interamente all'interno del **modello del continuo**. Esiste il continuo nella realtà microscopica della materia? No, a livello microscopico ci sono particelle discrete e vuoto; il continuo è una costruzione ideale per descrivere la realtà macroscopica. Non discutiamo le equazioni cinetiche di Boltzmann o la dinamica molecolare: assumiamo il modello del continuo e pretendiamo che esso risponda alle nostre domande finché è applicabile. Quando arriviamo ai suoi limiti intrinseci, saremo autorizzati a cambiare modello, ma non prima di averne compreso appieno la portata.

Il nostro obiettivo ultimo non è la teoria delle caratteristiche in sé, bensì **dimensionare i sistemi fisici e ingegneristici**. Per progettare la camera di combustione o l'ugello del motore di un lanciatore, non basta dire che in regime nominale la pressione scende lentamente per cento secondi. Al momento dell'accensione, la pressione passa da 1 bar a 100 bar, e la temperatura da 300 K a 3000 K, in appena un decimo di secondo. Questo transitorio violentissimo genera onde di pressione che sollecitano la struttura: se non si dimensiona tenendo conto di questi fenomeni dinamici, il motore collassa.

---

### 3. Ripasso del sistema bidimensionale instazionario: onda di entropia e onda di scorrimento

Riprendiamo l'analisi del sistema bidimensionale instazionario per un gas ideale non viscoso. Dallo studio dell'equazione caratteristica per la matrice omogenea:
$$\det(\mathbf{Q} \cdot \mathbf{N} \, \mathbf{I} - \mathbf{M}) = 0$$
abbiamo individuato le radici reali. 

La prima radice, avente **molteplicità algebrica 2**, è data da:
$$\mathbf{Q} \cdot \mathbf{N} = 0$$
dove $\mathbf{Q} = (\mathbf{V}, 1)^T = (u, v, 1)^T$ è il vettore velocità esteso nello spazio-tempo e $\mathbf{N} = (\nu_x, \nu_y, -\lambda)^T$. Questa condizione geometrica individua una **stella di piani** avente come asse la **traiettoria della particella fluida** $\frac{d\mathbf{x}}{dt} = \mathbf{V}$.

A questa radice doppia corrispondono due autovettori sinistri linearmente indipendenti, che generano due distinte famiglie di onde:

1. **Onda di entropia**:
   $$ds = dp - a^2 d\rho$$
   L'equazione di compatibilità associata esprime l'evoluzione dell'entropia lungo la traiettoria:
   $$\frac{Ds}{Dt} = \frac{\partial s}{\partial t} + \mathbf{V} \cdot \nabla s = 0$$
   Questa è un'equazione di compatibilità «pura»: la proprietà si conserva identicamente lungo la direzione della traiettoria fluidodinamica:
   $$ds = 0 \quad \text{lungo} \quad \frac{d\mathbf{x}}{dt} = \mathbf{V}$$
   Se in un punto la derivata temporale $\frac{\partial s}{\partial t}$ è diversa da zero, ciò accade a causa del trasporto convettivo associato a un gradiente spaziale preesistente:
   $$\frac{\partial s}{\partial t} = -\mathbf{V} \cdot \nabla s$$
   Ad esempio, in un processo di combustione reale non uniforme nello spazio, si crea un gradiente di entropia $\nabla s \neq 0$, che evolve nel tempo per pura convezione lungo la traiettoria della particella.

2. **Onda di scorrimento (shear wave / vorticità)**:
   L'autovettore sinistro corrispondente è:
   $$\mathbf{L} = (0, -\nu_y, \nu_x, 0)$$
   che definisce la variabile d'onda:
   $$d\sigma = -\nu_y du + \nu_x dv = \mathbf{\tau} \cdot d\mathbf{V}$$
   dove $\mathbf{\tau} = (-\nu_y, \nu_x)$ rappresenta il versore tangenziale al fronte d'onda, perpendicolare al versore normale $\mathbf{\nu} = (\nu_x, \nu_y)$.
   
   Moltiplicando il sistema originario per questo autovettore, si ottiene l'equazione di compatibilità per lo scorrimento:
   $$\mathbf{\tau} \cdot \left( \frac{\partial \mathbf{V}}{\partial t} + (\mathbf{V} \cdot \nabla)\mathbf{V} \right) + \frac{1}{\rho} \mathbf{\tau} \cdot \nabla p = 0$$
   ovvero:
   $$\mathbf{\tau} \cdot \frac{D\mathbf{V}}{Dt} + \frac{1}{\rho} \mathbf{\tau} \cdot \nabla p = 0$$
   
   A differenza dell'onda di entropia, questa equazione **non** asserisce che la variabile d'onda $\mathbf{\tau} \cdot \mathbf{V}$ si mantenga costante lungo la traiettoria. Compare infatti il termine residuo $\frac{1}{\rho} \mathbf{\tau} \cdot \nabla p$. Se il gradiente di pressione possiede una componente tangenziale lungo $\mathbf{\tau}$, essa agisce come termine sorgente che modifica la quantità di moto tangenziale. Tale termine si annullerebbe unicamente scegliendo $\mathbf{\tau}$ ortogonale al gradiente di pressione (ossia allineando $\mathbf{\nu}$ con $\nabla p$), ma poiché la direzione di $\nabla p$ varia punto per punto e istante per istante nel campo di moto, tale scelta non è praticabile in termini generali. 

Le equazioni di compatibilità ci forniscono le **direzioni di propagazione** e le **variabili caratteristiche associate**, ma non implicano necessariamente che la variabile d'onda sia un invariante differenziale privo di termini di sorgente o di accoppiamento multidimensionale.

---

### 4. La radice acustica, il cono di Mach e le bicaratteristiche

Prendiamo ora l'altra coppia di radici dell'equazione caratteristica, corrispondenti ai fenomeni acustici:
$$\mathbf{Q} \cdot \mathbf{N} = \pm a$$
dove $a = \sqrt{\frac{\gamma p}{\rho}}$ è la velocità locale del suono.

Per analizzare il significato geometrico nello spazio-tempo $(x, y, t)$, definiamo:
* $\mathbf{Q} = (\mathbf{V}, 1)^T = (u, v, 1)^T$;
* $\mathbf{N} = (\nu_x, \nu_y, -\lambda)^T$;
* il versore planare $\hat{\mathbf{\nu}} = (\nu_x, \nu_y, 0)^T$, per cui $\hat{\mathbf{\nu}} \cdot \mathbf{N} = \nu_x^2 + \nu_y^2 = 1$.

La relazione caratteristica può essere riscritta nella forma:
$$(\mathbf{Q} \mp a \hat{\mathbf{\nu}}) \cdot \mathbf{N} = 0$$
Poiché il versore normale $\mathbf{\nu}$ sul piano fisico può ruotare liberamente su un angolo giro ($0 \le \theta \le 2\pi$), l'insieme di tutti i possibili vettori direzionali:
$$\frac{d\mathbf{x}}{dt} = \mathbf{V} \pm a \mathbf{\nu}$$
costituisce le **generatrici di una superficie conica** nello spazio-tempo tridimensionale: il **cono di Mach** (o conoide caratteristico). 
* L'asse del cono è orientato lungo la traiettoria estesa $\mathbf{Q} = (\mathbf{V}, 1)^T$;
* La circonferenza di raggio $a$ rappresenta la direttrice;
* Le singole direzioni generatrici tangenti al cono prendono il nome di **bicaratteristiche**;
* I piani tangenti al cono inviluppati dalle bicaratteristiche sono i **piani acustici**.

L'orientamento del cono rispetto all'asse temporale $t$ dipende direttamente dal **numero di Mach** locale:
$$M = \frac{|\mathbf{V}|}{a}$$
1. **Regime subsonico ($M < 1$)**: $|\mathbf{V}| < a$. L'apertura del cono è tale che l'asse temporale $t$ è contenuto **all'interno** del cono di Mach. Le onde possono propagarsi in tutte le direzioni spaziali (a monte e a valle), riflettendo la natura ellittica della porzione spaziale del problema.
2. **Regime sonico ($M = 1$)**: $|\mathbf{V}| = a$. L'asse temporale giace esattamente sulla superficie del cono.
3. **Regime supersonico ($M > 1$)**: $|\mathbf{V}| > a$. Il cono è inclinato in avanti a tal punto che l'asse temporale cade **all'esterno** del cono. Tutti i fronti d'onda sono confinati all'interno del cono rivolto a valle: nessun disturbo acustico può risalire la corrente.
4. **Regime ipersonico ($M \gtrsim 5 - 10$)**: all'aumentare estremo del numero di Mach, le variazioni di temperatura ed energia cinetica attraverso i fronti d'onda divengono tali da attivare fenomeni di dissociazione molecolare, ionizzazione e non-equilibrio chimico-vibrazionale, richiedendo estensioni del modello termodinamico.

---

### 5. Calcolo dell'autovettore sinistro e definizione della variabile d'onda acustica

Sostituiamo la radice $\mathbf{Q} \cdot \mathbf{N} = a$ nella matrice del problema omogeneo per determinare l'autovettore sinistro riga $\mathbf{L} = (L_1, L_2, L_3, L_4)$:
$$\mathbf{L} \begin{pmatrix}
a & \rho \nu_x & \rho \nu_y & 0 \\
0 & a & 0 & \frac{\nu_x}{\rho} \\
0 & 0 & a & \frac{\nu_y}{\rho} \\
0 & \rho a^2 \nu_x & \rho a^2 \nu_y & a
\end{pmatrix} = \mathbf{0}^T$$

Sviluppando il prodotto per le quattro colonne:
1. Prima colonna:
   $$a L_1 = 0 \implies L_1 = 0$$
   L'autovettore acustico non contiene componenti dirette sulla prima variabile primitiva (la densità $\rho$). La variabile d'onda acustica non vede direttamente la densità.
2. Seconda colonna:
   $$\rho \nu_x L_1 + a L_2 + \rho a^2 \nu_x L_4 = 0 \implies a L_2 + \rho a^2 \nu_x L_4 = 0$$
3. Terza colonna:
   $$\rho \nu_y L_1 + a L_3 + \rho a^2 \nu_y L_4 = 0 \implies a L_3 + \rho a^2 \nu_y L_4 = 0$$
4. Quarta colonna:
   $$\frac{\nu_x}{\rho} L_2 + \frac{\nu_y}{\rho} L_3 + a L_4 = 0$$

Il sistema è omogeneo ed è determinato a meno di una costante arbitraria moltiplicativa. Poniamo convenzionalmente:
$$L_4 = 1$$
Semplificando per $a \neq 0$, otteniamo immediatamente:
$$L_2 = -\rho a \nu_x, \qquad L_3 = -\rho a \nu_y$$
Verifichiamo la quarta colonna:
$$\frac{\nu_x}{\rho}(-\rho a \nu_x) + \frac{\nu_y}{\rho}(-\rho a \nu_y) + a(1) = -a(\nu_x^2 + \nu_y^2) + a = -a(1) + a = 0$$
L'autovettore sinistro associato alla propagazione acustica è dunque:
$$\mathbf{L} = (0, -\rho a \nu_x, -\rho a \nu_y, 1)$$

La **variabile d'onda acustica** $d\sigma$ si costruisce proiettando il vettore delle variazioni delle variabili primitive $d\mathbf{U} = (d\rho, du, dv, dp)^T$ sull'autovettore sinistro $\mathbf{L}$:
$$d\sigma = \mathbf{L} \, d\mathbf{U} = dp - \rho a (\nu_x du + \nu_y dv) = dp - \rho a (\mathbf{\nu} \cdot d\mathbf{V})$$

---

### 6. Il concetto fondamentale di impedenza acustica $\rho a$ ed esempi applicativi

La quantità fondamentale emersa dal calcolo è:
$$Z = \rho a$$
denominata **impedenza acustica** del mezzo.

Essa stabilisce il fattore di proporzionalità dimensionale tra le variazioni di pressione e le variazioni di velocità attraverso il fronte d'onda:
$$dp \sim \rho a \, |\mathbf{\nu} \cdot d\mathbf{V}|$$
A parità di disturbo di velocità imposto, maggiore è l'impedenza acustica $\rho a$, maggiore sarà il salto di pressione generato.

Un'applicazione pratica cruciale nella propulsione spaziale riguarda il lanciatore europeo Vega. Nei motori a propellente solido (come lo Zefiro o il P80), per pressurizzare il sistema o gestire l'accensione è necessario impiegare un gas pressurizzante. Se si utilizzasse l'azoto ($N_2$, massa molecolare 28), a fronte di un gas di combustione ad alta temperatura (circa 3000 K), l'impedenza acustica dell'azoto freddo risulterebbe fortemente disaccoppiata rispetto a quella dei gas caldi, fungendo da barriera acustica riflettente che innesca violente risonanze e instabilità termoacustiche distruttive. 

Utilizzando invece l'**elio** ($He$, massa molecolare 4), il rapporto tra le masse molecolari è:
$$\frac{\mathcal{M}_{N_2}}{\mathcal{M}_{He}} = \frac{28}{4} = 7$$
Poiché la velocità del suono scala inversamente con la radice della massa molecolare ($a \propto \sqrt{T/\mathcal{M}}$) e la densità scala direttamente con la massa molecolare ($\rho \propto p \mathcal{M}/T$), l'impedenza acustica $\rho a \propto p \sqrt{\mathcal{M}/T}$ dell'elio a temperatura ambiente risulta straordinariamente vicina a quella dei gas caldi a 3000 K. In questo modo il disadattamento di impedenza viene quasi completamente annullato, il fronte d'onda non subisce riflessioni interne coerenti e le oscillazioni acustiche vengono soppresse, salvaguardando l'integrità strutturale del lanciatore.

---

### 7. Derivazione dell'equazione di compatibilità acustica e insorgenza di termini non direzionali

Per ricavare l'equazione di compatibilità che governa la dinamica della variabile d'onda acustica, moltiplichiamo le equazioni di bilancio di partenza per i coefficienti dell'autovettore sinistro:
1. Moltiplichiamo l'equazione della pressione (energia interna):
   $$\frac{\partial p}{\partial t} + \mathbf{V} \cdot \nabla p + \rho a^2 \nabla \cdot \mathbf{V} = 0$$
   per il coefficiente $L_4 = 1$.
2. Moltiplichiamo l'equazione vettoriale della quantità di moto:
   $$\frac{\partial \mathbf{V}}{\partial t} + (\mathbf{V} \cdot \nabla)\mathbf{V} + \frac{1}{\rho}\nabla p = \mathbf{0}$$
   scalarmente per $-\rho a \mathbf{\nu}$.

Sommando le due espressioni:
$$\left( \frac{\partial p}{\partial t} + \mathbf{V} \cdot \nabla p \right) - \rho a \mathbf{\nu} \cdot \left( \frac{\partial \mathbf{V}}{\partial t} + (\mathbf{V} \cdot \nabla)\mathbf{V} \right) + \rho a^2 \nabla \cdot \mathbf{V} - a \mathbf{\nu} \cdot \nabla p = 0$$

Raggruppando i termini temporali nella variazione della variabile d'onda:
$$\left( \frac{\partial p}{\partial t} - \rho a \mathbf{\nu} \cdot \frac{\partial \mathbf{V}}{\partial t} \right) + \mathbf{V} \cdot \nabla p - \rho a \mathbf{\nu} \cdot [(\mathbf{V} \cdot \nabla)\mathbf{V}] + \rho a^2 \nabla \cdot \mathbf{V} - a \mathbf{\nu} \cdot \nabla p = 0$$

In questa forma emergono chiaramente due difficoltà:
1. **Mancanza della direzione bicaratteristica**: la derivata convettiva spaziale è calcolata lungo la velocità di convezione del fluido $\mathbf{V}$, anziché lungo la velocità di propagazione del segnale acustico $\mathbf{V} + a\mathbf{\nu}$.
2. **Termini non direzionali residui**: sono presenti il termine $\rho a^2 \nabla \cdot \mathbf{V}$ e $-a \mathbf{\nu} \cdot \nabla p$, che non sono associabili a una derivata direzionale lungo una singola linea caratteristica. La divergenza del campo di velocità $\nabla \cdot \mathbf{V}$ è un operatore isotropo che riceve contributi da tutte le direzioni dello spazio.

L'operazione successiva consisterà nell'introdurre formalmente la direzione della bicaratteristica $\mathbf{V} + a\mathbf{\nu}$ aggiungendo e sottraendo i termini mancanti, e nell'utilizzare il teorema di dipendenza lineare delle direzioni sul cono (Teorema di Butler) per eliminare rigorosamente tutti i termini non direzionali.
