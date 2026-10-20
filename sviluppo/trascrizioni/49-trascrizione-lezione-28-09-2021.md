# Trascrizione integrale della lezione del 28 settembre 2021
## Aerotermodinamica e Propagazione Ondosa nei Fluidi Comprimibili

---

### 1. Fenomeni aeroacustici, instabilità di combustione e limiti dei modelli teorici

Tutti quelli che vanno sotto la categoria di fenomeni di tipo aeroacustico o termoaeroacustico, insomma, si possono declinare in vari modi. Intendiamo fenomeni nei quali sono presenti dinamiche di strati di scorrimento, di mulinelli di vorticità turbolenta, che generano rumore. 

Tipicamente dietro un motore c'è rumore. Tipicamente un motore a razzo fa rumore, stando all'esterno. Un po' meno tipicamente — nel senso che magari è una cosa meno nota — i motori suonano: suonano all'interno. In particolare i motori a razzo suonano se vengono fatti male, se non vengono ben disegnati, ben controllati; suonano. Tutta la vicenda dello sviluppo dei motori a combustione è stata segnata proprio da questo problema, cioè dal problema dell'instabilità di combustione, che non è un problema semplicemente che brucia male, no: induce vibrazioni, induce fenomeni di tipo ondoso.

Bisogna staccare le due cose: bisogna staccare le onde dalla combustione, bisogna staccare le onde dalla dinamica della vorticità, altrimenti funzionano male, o al limite usiamo le onde per comandare la combustione. Si cerca di fare queste cose, non so adesso con quanto successo siano state fatte, però nel campo della produzione dell'energia, per esempio, si cerca di fare anche questo: migliorare la combustione comandandola attraverso fenomeni di propagazione ondosa, attraverso onde che viaggiano a una certa velocità, a una certa frequenza, in modo tale da parlarsi con la combustione. 

Quindi tutto questo, ovviamente, se uno dice: «Ma sto facendo un modello ideale, tutto questo non si può vedere». Non è poi così vero che non si può vedere, ma qui entriamo in un mondo che è quello poi delle applicazioni dei modelli teorici. I modelli teorici vengono stiracchiati, vengono usati lì dove non si possono usare.

Per dire una cosa che è molto vicina a voi, il modello quasi-unidimensionale viene stiracchiato per farlo funzionare il più possibile perché è un modello molto semplice ma molto ricco, perché l'effetto di variazione d'area è una delle cose fondamentali di cui dobbiamo tenere conto. Quindi, se riusciamo ad applicarlo anche laddove nominalmente non si dovrebbe applicare, beh, siamo contenti, lo facciamo, ci proviamo. Però, è vero, qui entriamo in quel mondo che è l'applicazione del modello, che è un'arte proprio delle arti che più arti non potrebbero esistere, quindi lo mettiamo da parte. 

Veniamo a noi: facciamo questa ipotesi consapevoli del fatto che quello che studiamo è l'aspetto della propagazione ondosa in sé. Ok? Quindi arriveremo poi a capire questo, e a non sapere tutto quello che ho nominato; non lo sappiamo e quindi, vabbè, quello poi lo metteremo a posto in un altro momento. Da dove partiamo?

---

### 2. Bilancio della quantità di moto in forma integrale

Partiamo dall'equazione di conservazione della quantità di moto, che è l'equazione fondamentale. Ovviamente l'equazione del moto è l'equazione di Newton del moto riscritta per un continuo. 

Di cosa abbiamo bisogno? Abbiamo bisogno di identificare il nostro continuo, di dargli l'immagine giusta. Quindi prendiamo una regione del continuo dove vive il fluido, dove vive il nostro gas, che viene chiamata per comodità, per semplicità, particella, no? Si dice la particella fluidodinamica, gasdinamica... Che cos'è? Boh, a seconda delle situazioni è una cosa piuttosto che un'altra. 

Ecco, allora noi prendiamo un volume $V$, c'abbiamo un volume. Dentro questo volume è contenuto del gas, gas ideale. Che proprietà ha questo gas? Beh, le proprietà alcune le abbiamo già dette: ovviamente sta a una certa pressione, sta a una certa temperatura, quindi sta a una certa densità. Abbiamo una massa contenuta in questo volume, che è la massa della particella: un po' estesa, la particella è grossa, ma non si dice da nessuna parte che debba essere piccola. Quindi abbiamo una massa. Oltre a queste proprietà — cioè $\rho \cdot V$ è la massa — la particella si muove; questo volume contiene le particelle, le prendiamo in media. Questa particella si muove con una velocità $\mathbf{V}$. Sono tutte quante funzioni di spazio e tempo. Lo spazio in generale è $\mathbb{R}^3$, c'abbiamo tre componenti della velocità, quindi è un vettore ($\underline{V}$ con la barra sotto è un vettore) e quindi abbiamo una quantità di moto.

Quantità di moto che è un concetto abbastanza elaborato rispetto a dire semplicemente $m \mathbf{a}$. È un po' più sottile il fatto di parlare di quantità di moto, ma il concetto sostanzialmente è quello, non è che stiamo parlando di un'altra cosa, no? Quindi abbiamo la quantità di moto $\rho \mathbf{V}$. Che gli succede a questa quantità di moto? Cioè, vado rapido: sto riprendendo ovviamente cose che voi conoscete in teoria, ma insomma riguardiamole insieme. Quindi un po' di pazienza per chi se ne ricorda molto bene, e spero che per quelli che non se lo ricordino sia sufficiente per fissare le idee. 

Stiamo guardando, quindi dicevo, un $m \mathbf{a}$, e nel discorso elementare: dato un corpo con una certa massa, c'ha una certa velocità... Parlo di accelerazione — scusate, poi smetterò di parlare di $a$ come accelerazione, perché noi usiamo $a$ per indicare il suono, quindi poi dobbiamo smettere di parlare di $a$ come accelerazione — però parlando di Newton si dice: «Vabbè, c'ho una massa, c'ho una velocità, c'è un'accelerazione». Quindi quello che stiamo facendo è dire: «Vabbè, prendi questo oggetto che è $m \mathbf{V}$, la quantità di moto, e vediamo che cosa gli succede». 

E che cosa gli succede alla quantità di moto? Alla quantità di moto succede che in generale può variare. Perché può variare? Perché ci sono delle forze che sono imposte in quel punto, su quella regione, sul fluido. 

E quali sono le forze che noi dobbiamo considerare, visto che abbiamo fatto un modello ideale? Sono le forze di pressione. Le uniche che consideriamo sono forze di pressione, le quali agiscono sulla frontiera del volume. Quindi prendiamo questo volume che è la particella: questo volume $V$ ha una frontiera, quindi abbiamo il volume $V$ e abbiamo la frontiera $\partial V$. Sulla frontiera di $V$ fissiamo la normale in un punto, uscente. La regola è che deve essere uscente, $\mathbf{n}$. Lungo questa normale agisce la pressione. Qualunque pezzo della frontiera prendo, comunque orientata questa frontiera, agisce la pressione, quindi è diretta normalmente alla superficie. 

Come lo prendiamo questo volume? Le forze le abbiamo dette. Come lo prendiamo questo volume? Lo prendiamo fisso. Prendiamo un'immagine semplice: prendiamo il volume fisso, non si muove nel tempo. E quindi attraverso la frontiera abbiamo flussi. Siccome abbiamo una velocità diversa da zero, vuol dire che abbiamo attraverso la frontiera di questo volume un flusso, entrante o uscente, quello che è, dipende dal segno della velocità. 

Velocità che è diretta come le pare, ovviamente: è una proprietà del campo, sarà fatta come è fatta. Nel momento in cui definisco il volume e definisco la frontiera del volume, ovviamente il flusso è dato dalla proiezione di $\mathbf{V}$ lungo la normale $\mathbf{n}$. Quello che scorre tangenzialmente sulla superficie non è flusso, né entra né esce; devo misurare quello che entra ed esce. Quindi il flusso è comandato da $\mathbf{V} \cdot \mathbf{n}$. Flusso di chi? Beh, flusso della quantità di moto, ok?

Quello che entra ed esce, la proprietà che stiamo guardando, è la quantità di moto. Quindi devo controllare quello che entra ed esce, che fluisce attraverso la frontiera: perché se lo svuoto tutto non vedrò più la stessa quantità di moto, se lo riempio tutto non vedrò più la stessa quantità di moto. Quindi abbiamo questi due effetti. È chiaro che se il volume fosse mobile e si muovesse come il fluido questo pezzo non ci sarebbe, ma non ci piace: la cosiddetta visione lagrangiana ne facciamo volentieri a meno, ok? Quindi stiamo dentro quella che viene chiamata la visione euleriana, con un tipico abuso di etichette che non corrispondono quasi mai a quello che è avvenuto nella storia della scienza, ma non importa. Non è vero: uno ha fatto Lagrange, l'altro Euler... Come per il numero di Mach: non è vero che si chiama numero di Mach perché Mach ha detto «si chiama numero di Mach», non è vero che l'ha inventato lui; l'ha inventato Mach e poi gli altri hanno detto «vabbè, diamo il numero di Mach a Mach». Vabbè, tutte cose di poca importanza.

Quindi bilancio: a questo punto possiamo fare il bilancio della quantità di moto. La variazione nel tempo — questo è quello che vogliamo misurare — sul volume $V$ della quantità di moto: quindi questo è l'integrale su $V$. Questa cosa varia. Qui c'è la convenzione dei segni, ok? Più cosa? Più l'integrale sulla superficie di $V$ del flusso di quantità di moto: quindi abbiamo un $\rho \mathbf{V} (\mathbf{V} \cdot \mathbf{n}) \, dS$. $dV$ è l'elementino di volume, $dS$ è l'elementino della superficie della frontiera del volume. Più le forze applicate sulla frontiera, quindi più l'integrale sulla frontiera di $p \mathbf{n} \, dS$:

$$\frac{d}{dt} \int_V \rho \mathbf{V} \, dV + \int_{\partial V} \rho \mathbf{V} (\mathbf{V} \cdot \mathbf{n}) \, dS + \int_{\partial V} p \mathbf{n} \, dS = 0$$

E questo fa zero. Fa zero, ovviamente, nel senso che la variazione nel tempo della quantità di moto è data dalla somma di questi due; questa è la scrittura classica, tutto questo fa zero. Dal punto di vista pratico operativo non è che uno lo dovrebbe scrivere così, ma lo dovrebbe scrivere nel senso di dire: «Questo varia perché ci stanno questi». Ok, già questo va bene, è più operativo. È nel senso di $\mathbf{F} = m \mathbf{a}$: che significa $\mathbf{F} = m \mathbf{a}$ se non scrivere $m \mathbf{a} = \mathbf{F}$? Questa $\mathbf{F}$ è un po' più ricca perché non ci sono solo le forze, ma abbiamo anche il flusso della quantità di moto, perché abbiamo preso questo volume, siamo fermi e quindi vediamo quello che entra ed esce. $m \mathbf{a} = \mathbf{F}$ è la visione lagrangiana: prendo la massa, la seguo con le forze applicate, vedo come cambia la velocità tenendo ferma la massa. Ok? Quindi questa è l'equazione di conservazione della quantità di moto scritta in forma integrale, cioè è il bilancio.

Diciamo che quando siamo fortunati e riusciamo a scrivere dei bilanci in forma integrale, siamo su un punto molto molto saldo. Da qua si parte per fare tutto il resto. Non è detto che si riesca a fare sempre. Esistono problemi fisici nei quali fare il bilancio dal punto di vista integrale non è una cosa che riesce, o per meglio dire non riesce tanto facilmente; la prima volta non ce la facciamo. Però in questo caso invece è semplice. Nel costruirlo ci dice come varia la quantità di moto nel tempo: devo considerare due elementi, due fattori, ovvero il flusso di questa proprietà (la quantità di moto) e le forze applicate, che siccome abbiamo preso il modello ideale sono solo la pressione. E se ci mettevavamo pure gli sforzi dissipativi? Vabbè, ci mettevamo gli sforzi dissipativi: invece della pressione, che è un oggetto scalare moltiplicato per $\mathbf{n}$, avremmo avuto un tensore e avremmo fatto una cosa un pochettino più complicata, ma insomma il concetto è quello. Non cambia dal punto di vista concettuale.

---

### 3. Ipotesi di regolarità e passaggio alla forma a divergenza

Che succede se queste proprietà — densità, pressione, temperatura, velocità — sono proprietà derivabili? Quando chiediamo che siano derivabili almeno, come dire, non dappertutto ma in alcune zone, quello che stiamo chiedendo è che abbiano dei comportamenti molto regolari. La derivata c'è, riesco a calcolarla. Ripeto: posso accettare che in un punto ce ne siano due di derivate. L'importante è che sappiamo dove sta questo salto: vengono fuori due valori, poi metteremo insieme i due valori e sistemiamo le cose, ma insomma non vogliamo trovarci in una situazione nella quale la grandezza non sia derivabile. 

È una definizione di tipo matematico che chiaramente corrisponde a un qualche cosa di fisico: noi ci aspettiamo che in ogni punto siamo in grado di misurarla, questa proprietà. In ogni punto sappiamo determinare quant'è la temperatura, quant'è la pressione, quant'è la velocità. Non vorremmo trovarci in condizioni nelle quali in quel punto non siamo capaci di misurarla. E cerchiamo che in quel punto ci sia un gradiente, quindi non siano proprietà costanti: hanno delle variazioni. Possiamo anche accettare, come vedremo, che queste variazioni siano diverse a seconda di dove guardo: mentre in un punto guardo da una parte e vedo che la derivata ha un certo valore, mi metto a guardare dall'altra parte e vedo che ne ha un altro. Chiaramente abbiamo due valori della derivata, però la funzione è continua ed è derivabile a tratti. L'unico problema che abbiamo è che quando abbiamo due valori diciamo: «Poi che facciamo? Lo sistemiamo».

Supponiamo di stare in condizioni derivabili: allora si può applicare il teorema di Gauss-Green, il teorema della divergenza, e ci ricordiamo che questo volume è fisso, quindi questi due operatori possono essere scambiati. Invece di fare la variazione nel tempo dell'integrale, faccio l'integrale della variazione nel tempo. Porto tutto sotto lo stesso simbolo di integrazione di volume. 

Qui abbiamo usato il $d$ piccolo perché è una derivata di un oggetto che non ha più dimensioni spaziali, ma ha una dipendenza solo dal tempo. Nel momento in cui lo portiamo dentro dobbiamo passare dal $d$ piccolo al $d$ distorto ($\partial$) o, equivalentemente, all'indice in basso, che è quello che faccio: quindi abbiamo $(\rho \mathbf{V})_t$. La divergenza applicata a questo integrale di superficie lo trasforma in un integrale di volume, quindi abbiamo la divergenza applicata a $\rho \mathbf{V} \mathbf{V}$. Questo è un tensore, quindi quando facciamo la divergenza viene fuori un oggetto che è un vettore; deve essere tale perché deve bilanciare questo termine temporale che è un vettore.

Più applico anche qua il teorema della divergenza, quello scalare, quindi la pressione diventa un gradiente: $+\nabla p$. Tutto questo adesso è andato a finire in un integrale in $dV$ uguale a zero:

$$\int_V \left[ \frac{\partial (\rho \mathbf{V})}{\partial t} + \nabla \cdot (\rho \mathbf{V} \mathbf{V}) + \nabla p \right] dV = 0$$

Inversione dei due operatori nel primo termine e applicazione del teorema della divergenza di Gauss-Green sugli integrali di superficie: gli integrali di superficie li portiamo a volume. 

Questo sotto l'ipotesi che sia possibile fare le derivate. Quindi tutte queste derivate qua sono evidentemente fattibili, il campo è regolare, non ci sono punti dove non sappiamo qual è la derivata, cioè un salto. Se noi abbiamo un salto netto, una discontinuità, la derivata qual è? Non la sappiamo. È vero che abbiamo un valore a destra e uno a sinistra, ma c'è il problema che la funzione è discontinua; in realtà uno potrebbe definire una derivata in un senso particolare del termine, alla Heaviside o con le distribuzioni, ma stiamo escludendo quello. 

Questa scrittura qua in forma integrale non esclude la presenza delle discontinuità; da qui noi ricaviamo poi la teoria delle discontinuità (le onde d'urto). Questa scrittura differenziale qua invece non vale nel caso in cui si presentino discontinuità. 

Assumendo un campo regolare, ovviamente l'integrando lo possiamo liberare dall'integrale, perché deve essere vero per ogni volume arbitrario. Quindi scriviamo direttamente, facciamo vivere separatamente questa espressione qua:

$$\frac{\partial (\rho \mathbf{V})}{\partial t} + \nabla \cdot (\rho \mathbf{V} \mathbf{V}) + \nabla p = 0$$

Questa è quella che si chiama la **forma di divergenza**, che viene chiamata anche forma conservativa, anche se ha un po' poco significato fisico chiamarla così; sicuramente «forma di divergenza» ha significato perché stiamo operando con delle divergenze, quindi letteralmente è vero. 

---

### 4. Accoppiamento con la conservazione della massa e forme quasi-lineari

Questa espressione la possiamo sviluppare: possiamo sviluppare le derivate. Sviluppando abbiamo:

$$\rho \frac{\partial \mathbf{V}}{\partial t} + \mathbf{V} \frac{\partial \rho}{\partial t} + \rho (\mathbf{V} \cdot \nabla)\mathbf{V} + \mathbf{V} \nabla \cdot (\rho \mathbf{V}) + \nabla p = 0$$

Abbiamo spacchettato sia la derivata temporale che il termine convettivo, quello associato per l'appunto al flusso di convezione. E riconosciamo che possiamo mettere a fattor comune $\mathbf{V}$:

$$\rho \frac{\partial \mathbf{V}}{\partial t} + \rho (\mathbf{V} \cdot \nabla)\mathbf{V} + \mathbf{V} \left[ \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{V}) \right] + \nabla p = 0$$

Per chi è già imparato, ovviamente riconosce che dentro quella parentesi c'è un'altra equazione: l'equazione di conservazione della massa in forma di divergenza! 

Che dobbiamo fare per ricavarla? Facciamo la stessa cosa che abbiamo fatto prima. Presto detto: la proprietà che consideriamo non è la quantità di moto, ma è la massa della particella. La massa della particella può variare? Certo che può variare per effetto del campo di moto: non è detto che troviamo sempre lo stesso valore della densità di massa nel punto. Facciamo lo stesso ragionamento: prendo il volume, prendo la massa contenuta nel volume. Come può variare questa massa? Escludo la presenza di sorgenti e pozzi, quindi l'unica possibilità che abbiamo è il flusso attraverso la frontiera del volume. Quindi:

$$\frac{d}{dt} \int_V \rho \, dV + \int_{\partial V} \rho (\mathbf{V} \cdot \mathbf{n}) \, dS = 0$$

Se non c'è flusso, se il volume si muove con la velocità che hanno le particelle (se prendo un volume la cui frontiera è sempre la stessa, sono sempre le stesse particelle), ovviamente la massa contenuta là dentro non può variare. Siccome abbiamo preso un volume fisso, la massa può variare perché abbiamo un flusso. Se questo flusso dovesse essere zero nel netto, non vuol dire che non entra e non esce, ma che tanto entra quanto esce, cioè la massa contenuta nel volumetto non varia. 

È lo stesso ragionamento fatto prima applicato alla proprietà densità invece che alla quantità di moto; gli integrali si costruiscono sempre nello stesso modo. Certo, cambiando la proprietà dobbiamo considerare effetti diversi: sulla quantità di moto abbiamo detto il flusso più le forze; quando parleremo di energia col primo principio della termodinamica parleremo di lavoro delle forze. A seconda della quantità considerata si mettono dentro pezzi diversi, ma la strada che si segue è sempre la stessa.

Applichiamo il teorema di Gauss anche su questo, invertiamo gli operatori e arriviamo alla conservazione della massa in forma di divergenza:

$$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{V}) = 0$$

Ma allora quella parentesi moltiplicata per $\mathbf{V}$ nell'equazione della quantità di moto è uguale a zero! Vuol dire che posso togliere quei termini e portare l'equazione in una forma ancora più trasparente:

$$\rho \frac{\partial \mathbf{V}}{\partial t} + \rho (\mathbf{V} \cdot \nabla)\mathbf{V} + \nabla p = 0$$

Questa è la **forma quasi-lineare** pulita. Abbiamo tolto quello zero che era presente nell'equazione di divergenza e siamo arrivati a questa scrittura per la velocità. Palesemente la variabile che vive dentro questa equazione è la velocità, mentre la densità vive dentro l'equazione di conservazione della massa. Certo, non sono separabili tra loro: è un sistema non lineare accoppiato.

Anche l'equazione di conservazione della massa la possiamo sviluppare:

$$\frac{\partial \rho}{\partial t} + \mathbf{V} \cdot \nabla \rho + \rho \nabla \cdot \mathbf{V} = 0$$

Queste si chiamano forme quasi-lineari, dove il termine «quasi» sta ad indicare che non sono lineari, ma i coefficienti delle derivate non dipendono dalle derivate stesse: dipendono solo dalle variabili dipendenti incognite ($\mathbf{V}$, $\rho$). I coefficienti dipendono da queste variabili, ma non dalle loro derivate. «Quasi-lineare» è un termine matematico preciso.

*(Uno studente chiede chiarimenti sui passaggi tra un'equazione e l'altra).*

Il professore risponde:  
«Abbiamo sviluppato la forma di divergenza, abbiamo visto che c'erano termini proporzionali a $\mathbf{V}$, abbiamo messo a fattor comune la velocità $\mathbf{V}$ e abbiamo riconosciuto che la somma dentro la parentesi non è altro che uno zero, perché è l'equazione di continuità della massa in forma di divergenza. Come vedrete, le equazioni di bilancio contengono degli zeri. Vedremo poi quando scriveremo la conservazione dell'energia totale che in realtà sono due cose combinate insieme ciascuna uguale a zero: l'energia cinetica e l'energia interna. È una particolarità della struttura delle equazioni di bilancio. Sviluppando $\nabla \cdot (\rho \mathbf{V})$ abbiamo il termine convettivo $\mathbf{V} \cdot \nabla \rho$ e l'altro termine $\rho \nabla \cdot \mathbf{V}$».

---

### 5. Operatore derivata materiale e variazione lungo la traiettoria

Notiamo che entrambe le equazioni hanno la stessa struttura convettiva: $\frac{\partial}{\partial t} + \mathbf{V} \cdot \nabla$. Questa struttura ci permette di introdurre un operatore che chiamiamo $D$ grande:

$$\frac{D}{Dt} = \frac{\partial}{\partial t} + \mathbf{V} \cdot \nabla$$

Che cosa significa introdurre questo operatore? Significa introdurre il concetto di una **derivata direzionale**: derivata lagrangiana, derivata convettiva, derivata sostanziale. Quello che stiamo facendo è misurare le variazioni lungo una particolare direzione. 

Stiamo dicendo: che succede al differenziale totale di $\rho$? Mettiamoci in un punto, dimentichiamoci un attimo la fisica. Dato un sistema di coordinate nello spazio-tempo $(x, y, z, t)$, il differenziale totale vale:

$$d\rho = \frac{\partial \rho}{\partial t} dt + \nabla \rho \cdot d\mathbf{x}$$

Note le variazioni nel tempo e nello spazio, sappiamo quanto vale il differenziale totale, cioè la variazione totale della grandezza. Posso misurare questa cosa lungo una particolare direzione? Certo che sì: chi me lo impedisce? Scelgo come sistema di coordinate uno particolare, allineato lungo la traiettoria. 

Che cos'è la traiettoria? Si definisce traiettoria della particella l'equazione:

$$\frac{d\mathbf{x}}{dt} = \mathbf{V}$$

In ogni punto conosco $\mathbf{V}$: se prendo una curva che in ogni punto è tangente a $\mathbf{V}$, definisco una particolare direzione nello spazio-tempo che si chiama traiettoria (se il campo di moto fosse stazionario coinciderebbe con la linea di corrente). La particella fluidodinamica nello spazio-tempo viene trasportata con la velocità $\mathbf{V}$: la sua densità, la sua temperatura e la sua velocità stessa vengono trasportate con $\mathbf{V}$.

Scegliamo un sistema di coordinate ruotato, in cui una coordinata è allineata lungo la traiettoria e le altre giacciono sul piano normale. Ci chiediamo: che cosa succede muovendoci lungo questa direzione? Dividendo per $dt$ il differenziale totale lungo la direzione $d\mathbf{x}/dt = \mathbf{V}$, otteniamo:

$$\frac{d\rho}{dt}\Bigg|_{\text{traiettoria}} = \frac{\partial \rho}{\partial t} + \mathbf{V} \cdot \nabla \rho = \frac{D\rho}{Dt}$$

Guardando l'equazione di conservazione della massa, essa si riscrive immediatamente:

$$\frac{D\rho}{Dt} = - \rho \nabla \cdot \mathbf{V}$$

Come varia la densità lungo la traiettoria? Varia per effetto della divergenza della velocità, con coefficiente di proporzionalità pari a $\rho$. Questo concetto di allinearsi con le direzioni di variazione ritornerà continuamente, perché la teoria delle caratteristiche — che è la teoria della propagazione ondosa — si appoggia esattamente su questo ragionamento.

Perché la densità varia in questo modo lungo la traiettoria? È elementare se pensiamo al volumetto chiuso: immaginiamo un palloncino. Se la frontiera si allarga, la velocità sulla superficie ha una componente positiva rispetto alla normale uscente $\mathbf{n}$. La divergenza della velocità misura la variazione istantanea di volume nell'unità di tempo: se il volume si espande, la divergenza è positiva ($\nabla \cdot \mathbf{V} > 0$), quindi col segno meno davanti $\frac{D\rho}{Dt} < 0$, ovvero la densità diminuisce. Se invece il palloncino si comprime, la velocità è opposta alla normale, la divergenza è negativa, e con il segno meno la densità aumenta. Questo è il bilancio di massa: dobbiamo imparare a leggere la fisica dentro i simboli delle equazioni.

---

### 6. Campi a divergenza nulla e densità non costante

E se il campo di moto è solenoidale, cioè a divergenza nulla ($\nabla \cdot \mathbf{V} = 0$)?  
Spesso si dice: «Il campo di moto è a divergenza nulla, quindi il fluido è incompressibile, quindi la densità è costante». No! Non è vero, non funziona così. 

Se la divergenza è nulla, l'equazione ci dice unicamente che la densità non varia lungo la traiettoria:

$$\frac{D\rho}{Dt} = 0$$

Non dice che la densità è costante ovunque nello spazio! Dice che la densità si mantiene costante lungo ciascuna traiettoria. Se il campo di moto parte da una distribuzione iniziale di densità che non è uniforme nello spazio, nel corso del tempo avremo un campo a densità non costante. 

È di un qualche interesse questo discorso? Altroché: pensiamo alla combustione a bassissima velocità (a basso numero di Mach). Ci sono zone calde e zone fredde: dove è caldo la densità è bassa, dove è freddo la densità è alta. La densità non è costante nello spazio, quindi non è corretto usare un modello incompressibile a densità costante, anche se il campo di velocità può essere trattato a divergenza nulla. Discutere bene le equazioni aiuta moltissimo a evitare errori concettuali.

Guardando invece l'equazione della quantità di moto con la derivata materiale, abbiamo:

$$\rho \frac{D\mathbf{V}}{Dt} = - \nabla p$$

La variazione della velocità della particella lungo la sua traiettoria è direttamente proporzionale al gradiente delle forze di pressione. Se il gradiente di pressione lungo quella traiettoria è nullo, la velocità della particella non può cambiare.

---

### 7. Interruzione tecnica in aula

*(A questo punto della lezione si verifica lo spegnimento improvviso del computer del docente e la caduta del collegamento Zoom)*

**Docente:** «Non so che è successo... È morto tutto. Non ci sono... Non c'è... Ci sta il computer questa volta... Dove sono? Se ne sono andati tutti a spasso. Pensavo di aver infilato bene la spina... Non lo so. D'accordo, bene... Da dove...?»  
**Studente:** «Forse è un'altra aula...»  
**Intervento dall'aula:** «Domande? Domanda un c\*\*\*\*! Non ci stiamo a capire un c\*\*\*\* manco in aula!»  
**Docente:** «Mannaggia... Se non ve ne siete accorti era Luca Padrignani, presidente rappresentante... È ripartita la connessione Zoom, vediamo se ci sono ancora sopravvissuti... Se volete potete entrare... Oddio, sto senza mascherina... Vabbè, mancano 5 minuti».

---

### 8. Chiusura termodinamica: trasformazione isentropica e velocità del suono

Riprendiamo. Questo sistema costituito dalla conservazione della massa e della quantità di moto, scritto in forma quasi-lineare, è aperto. È aperto perché abbiamo due equazioni (una scalare per la massa e una vettoriale per la quantità di moto) ma abbiamo più incognite: densità $\rho$, velocità $\mathbf{V}$ e pressione $p$. Ci manca un'equazione per chiuderlo.

Possiamo operare in due modi:
1. Assumere o imporre una trasformazione termodinamica che leghi pressione e densità;
2. Introdurre il primo principio della termodinamica facendo il bilancio dell'energia.

Consideriamo la prima via: possiamo assumere una legge politropica generica:

$$\frac{p}{\rho^k} = \text{costante}$$

A seconda dell'esponente $k$ abbiamo diverse trasformazioni:
* Se $k = 1$, abbiamo la trasformazione isoterma (che per noi sulla Terra ha poco interesse nella dinamica dei gas, mentre ne ha per chi studia lo spazio intergalattico).
* Molto più comune e fondamentale è il caso $k = \gamma = \frac{c_p}{c_v}$ (rapporto tra i calori specifici a pressione e volume costante), che rappresenta la **relazione isentropica**, detta anche legge di Poisson.

In forma elementare possiamo definire l'entropia legata alla grandezza $s \sim \frac{p}{\rho^\gamma}$. Con questa relazione tra $p$ e $\rho$ il sistema è chiuso.

In gasdinamica si fa una distinzione fondamentale tra:
* **Moto isentropico**: l'entropia si mantiene costante lungo la traiettoria di ciascuna particella:
  $$\frac{Ds}{Dt} = \frac{\partial s}{\partial t} + \mathbf{V} \cdot \nabla s = 0$$
* **Moto omoentropico**: l'entropia è identicamente uniforme e costante ovunque nello spazio e nel tempo:
  $$ds = 0 \quad \forall \, \mathbf{x}, t \implies s(\mathbf{x}, t) = \text{costante}$$

Chiediamo che il moto sia isentropico lungo la traiettoria: $\frac{Ds}{Dt} = 0$. Differenziando la legge isentropica otteniamo la derivata logaritmica:

$$\frac{dp}{p} - \gamma \frac{d\rho}{\rho} = 0 \implies dp = \left(\frac{\gamma p}{\rho}\right) d\rho$$

Definiamo la **velocità del suono** secondo la formula standard dell'acustica:

$$a^2 = \frac{\gamma p}{\rho}$$

Quindi lungo la traiettoria:

$$dp = a^2 d\rho \quad \Longleftrightarrow \quad \frac{Dp}{Dt} = a^2 \frac{D\rho}{Dt}$$

Moltiplicando l'equazione di conservazione della massa per $a^2$:

$$a^2 \frac{\partial \rho}{\partial t} + a^2 \mathbf{V} \cdot \nabla \rho + \rho a^2 \nabla \cdot \mathbf{V} = 0$$

Poiché $a^2 \frac{D\rho}{Dt} = \frac{Dp}{Dt}$ e notando che $\rho a^2 = \gamma p$, l'equazione si riscrive in termini di pressione:

$$\frac{\partial p}{\partial t} + \mathbf{V} \cdot \nabla p + \gamma p \nabla \cdot \mathbf{V} = 0$$

Ovvero, in forma compatta con la derivata materiale:

$$\frac{Dp}{Dt} + \gamma p \nabla \cdot \mathbf{V} = 0$$

Lavorando con un modello isentropico chiuso, scrivere l'equazione per la densità o scriverla per la pressione è esattamente la stessa cosa: descrivono la medesima fisica.

---

### 9. Bilancio dell'energia interna e lavoro delle forze di pressione

Ma che cos'è la pressione? La pressione è legata direttamente all'energia interna del gas. L'energia interna per unità di massa per un gas ideale è $e_i = c_v T$. Quando parliamo di un fluido in movimento dobbiamo considerare l'energia interna per unità di volume: $\rho e_i = \rho c_v T$. Tramite l'equazione di stato dei gas perfetti e la definizione di $\gamma$, sappiamo che:

$$\rho e_i = \frac{p}{\gamma - 1}$$

Possiamo allora riscrivere l'equazione della pressione dividendo tutti i termini per $(\gamma - 1)$ (assumendo $\gamma$ costante):

$$\frac{\partial}{\partial t}\left(\frac{p}{\gamma - 1}\right) + \mathbf{V} \cdot \nabla \left(\frac{p}{\gamma - 1}\right) + \frac{\gamma p}{\gamma - 1} \nabla \cdot \mathbf{V} = 0$$

Vogliamo ora risalire concettualmente alla forma di divergenza. Riscriviamo il termine convettivo combinando pressione e velocità nel flusso $\nabla \cdot \left(\frac{p}{\gamma - 1} \mathbf{V}\right)$:

$$\frac{\partial}{\partial t}\left(\frac{p}{\gamma - 1}\right) + \nabla \cdot \left(\frac{p}{\gamma - 1} \mathbf{V}\right) + p \nabla \cdot \mathbf{V} = 0$$

Notate la struttura:
* Il primo termine è la variazione temporale dell'energia interna;
* Il secondo termine è la divergenza del flusso di energia interna;
* Ma compare un termine residuo: $+ p \nabla \cdot \mathbf{V}$.

Questo termine non è esprimibile come divergenza pura, quindi non possiamo applicare il teorema di Gauss-Green all'indietro per ricondurlo a un integrale di superficie sulla frontiera: rimane un termine di volume. 

Che cos'è questo termine $p \nabla \cdot \mathbf{V}$? È un pezzo del **lavoro delle forze di pressione** ($p \mathbf{V}$). Quando parliamo di lavoro meccanico nel continuo, consideriamo il prodotto tra forza e velocità. Il lavoro delle forze di pressione si scinde in due contributi:
1. Una parte legata alla variazione di volume (la divergenza della velocità), che va a modificare l'energia interna del fluido comprimendolo o espandendolo;
2. L'altra parte andrà ad alimentare la variazione di energia cinetica.

---

### 10. Conclusioni e cenni storici

In conclusione, se scriviamo il sistema formato da:
1. Conservazione della massa;
2. Conservazione della quantità di moto;
3. Conservazione dell'energia interna con vincolo isentropico;

abbiamo scritto esattamente lo stesso sistema di partenza chiuso dalla legge politropica. Tutto questo ci prepara ad affrontare il primo principio della termodinamica nella sua formulazione generale, dove considereremo sia l'energia interna sia l'energia cinetica. 

E chiudo veramente — giuro che finisco — ricordando che storicamente, agli albori della scienza moderna, ci fu una disputa accesissima tra i cartesiani, che sostenevano che tutto andasse descritto in termini di quantità di moto ($m \mathbf{V}$), e i leibniziani, che volevano descrivere il moto attraverso la «forza viva» (l'energia cinetica, $m V^2$). Alla fine non ha vinto nessuno dei due, perché avevano entrambi ragione: si trattava solo di capire di cosa stavamo parlando.

*(Fine lezione: ore 14:05).*
