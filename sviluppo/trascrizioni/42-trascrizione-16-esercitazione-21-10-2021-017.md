# Trascrizione dell'Esercitazione del 21 ottobre 2021
## Esercitazione: Moto del Pistone con Legge Parabolica, Espansione ed Equazioni di Compatibilità

---

### 1. Introduzione all'esercitazione e modalità didattica

Avrete ricevuto la comunicazione che da adesso in avanti bisogna fare la cosiddetta somministrazione dei questionari di valutazione della didattica, degli OPIS, di quegli affari là, e penso che ve la possiate cavare da soli. In teoria dovrebbe essere dedicato un pezzo di una lezione per questi questionari. Se c'è bisogno, se ritenete sia l'occasione per poter fare anche delle osservazioni in diretta, si può fare, però diciamo non voglio certo io deviarvi per questa cosa la cui utilità c'è, per carità, non voglio dire che non c'è, però insomma non è tale da pensare di dedicare poi un tempo a lezione per questo. 

Proviamo oggi ad applicare quello che abbiamo fatto. Cioè voi applicate quello che abbiamo fatto, il che è un po' difficile ovviamente nel sistema misto, quindi mi rendo conto che per quelli che sono in remoto può essere un po' più faticoso bucare la distanza, arrivare a dire, a fare, a commentare. Ci si mette anche il fatto che faccio fatica a sentire, non che voi fate rumore, solo io faccio fatica a sentire. Però siete caldamente invitati, anche quelli che sono in remoto, ma insomma tutti voi, a sforzarvi di partecipare, cioè di dire le cose che saranno sbagliate ovviamente, ma non è questo il punto. 

Lo scopo di fare esercitazione è quello, per l'appunto, di far uscire soprattutto le cose sbagliate in modo tale che poi le buttiamo via, ok? Ma solo se ce le mettiamo davanti, solo se le guardiamo negli occhi direttamente, riusciamo poi a farne a meno. Di questo ne possiamo essere certi. Peraltro quelli che sono in remoto ovviamente c'hanno anche il vantaggio che non ci mettono la faccia: è solo una voce, al più posso controllare il nome, ma insomma se perdo pure il tempo a controllare il nome non si finisce più. Mentre voi in aula ci mettete la faccia e a un certo punto ovviamente io associo faccia e nome. Dovrebbe essere chiaro che non sbagliare non serve, serve esattamente che adesso voi sbagliate. Quindi dite le cose. 

---

### 2. Definizione del problema: moto del pistone con legge parabolica

In che senso applicare quello che abbiamo fatto? Nel senso di passare al caso nel quale il pistone, invece di muoversi con moto impulsivo, si muove con **legge parabolica**. Abbiamo il condotto a sezione costante. Di nuovo non ci interessa adesso specificare il contorno: l'abbiamo discusso e tutto sommato la discussione che abbiamo fatto sul contorno poi non cambia molto; viene di conseguenza da quello che vedremo insieme che qualcosa cambia ovviamente nel cambiare la legge del pistone, però non l'essenza del fenomeno (pressione assegnata, velocità assegnata, pressione totale). 

Cosa succede? Quello succede chiaramente se arriva un'onda che è fatta in un certo modo. Se l'onda è fatta in una maniera differente, sarà qualcosa di differente, però l'essenza è quella, quindi di questo non ci preoccupiamo. 

Qui c'abbiamo il pistone. La legge del moto del pistone è del tipo:
$$x_p(t) = b t^2$$
dove consideriamo $b < 0$ (o se volete possiamo metterci il valore assoluto $-|b|t^2$ e vedere il segno). Questo è il tempo $t$, questo è lo spazio $x$, e quindi la traccia del moto del pistone si muove nel piano $(x, t)$ svergando verso sinistra perché stiamo guardando un'**espansione**.

La condizione iniziale all'interno del condotto è quella di **quiete**:
$$u(x, 0) = 0, \qquad a(x, 0) = a_0$$
Per la velocità del pistone, indicando con il punto la derivata temporale $\dot{x}_p = rac{dx_p}{dt}$ e con due punti la derivata seconda $\ddot{x}_p = rac{d^2 x_p}{dt^2} = 2b$:
$$\dot{x}_p(t) = 2b t$$
Notate che al tempo $t = 0$ abbiamo $\dot{x}_p(0) = 0$. Quindi **non abbiamo una discontinuità sulla velocità**: la velocità non passa istantaneamente da valore nullo a un valore finito come nel caso del moto impulsivo. La parabola nel piano $(x, t)$ parte dall'origine ed è tangente all'asse temporale $t$.

Quali sono le domande che ci facciamo su questo problema? Le stesse che ci siamo fatti ieri. Abbiamo integrato le equazioni di compatibilità, abbiamo gli **invarianti di Riemann**:
$$R_1 = a - \delta u = 	ext{costante lungo } rac{dx}{dt} = u - a$$
$$R_2 = a + \delta u = 	ext{costante lungo } rac{dx}{dt} = u + a$$
con $\delta = rac{\gamma - 1}{2}$. Il campo rimane sempre **omoentropico**. Vogliamo applicare questo modello senza essere costretti a prendere equazioni differenziali e integrarle numericamente, ma costruendo carta e matita analiticamente la soluzione per $t > 0$, cioè ricavare come dipendono da $x$ la velocità $u$ e la velocità del suono $a$ dal pistone in avanti.

Scegliamo una scala temporale tale per cui guardiamo il fenomeno prima che arrivi all'altra estremità del condotto. Queste soluzioni sono soluzioni simili, nel senso che le grandezze si mantengono costanti lungo ciascuna linea caratteristica. 

Una delle domande fondamentali che dobbiamo porci è: *sui fronti anteriore e posteriore, o da qualche altra parte, ci sarà una discontinuità di gradiente?* Nel caso impulsivo avevamo una discontinuità di gradiente sul fronte anteriore e posteriore. Qui siamo passati da una velocità a gradino a una legge parabolica continua: ci saranno ancora discontinuità sulle derivate?

---

### 3. Analisi del problema misto: condizioni al contorno e condizioni iniziali

Questo è un **problema misto**, ovvero con condizioni iniziali e condizioni al contorno. 

#### Condizione al contorno sul pistone
Il pistone costituisce il contorno mobile $x = x_p(t)$. Finché il gas rimane a contatto con il pistone e non va in cavitazione, il gas deve avere istantaneamente la stessa velocità del pistone:
$$u(x_p(t), t) = \dot{x}_p(t) = 2b t$$
Il pistone è infinitamente rigido, non ha una pressione propria ma subisce la pressione del gas. L'unico parametro cinematico imposto è la velocità. 

Dobbiamo però verificare il limite di cavitazione: poiché $b < 0$, al crescere del tempo la velocità $\dot{x}_p$ diventa sempre più negativa verso $-\infty$. Ma l'espansione massima possibile (espansione nel vuoto) si ha quando la velocità del suono sul pistone si annulla:
$$a_p \ge 0 \implies |\dot{x}_{p,	ext{max}}| = rac{a_0}{\delta} = rac{2}{\gamma - 1} a_0$$
Quindi dobbiamo considerare tempi sufficientemente piccoli tali che $a_p > 0$.

#### Condizioni iniziali e assenza di singolarità
Al tempo $t = 0$, sia il gas nel condotto sia il pistone hanno velocità rigorosamente nulla:
$$u(x, 0) = 0, \qquad \dot{x}_p(0) = 0$$
A differenza del moto impulsivo, non vi è alcuna singolarità o salto istantaneo di velocità all'istante iniziale nell'origine $(0, 0)$.

---

### 4. Digressione metodologica: verifica dell'omoentropia e trattamento dei dati sperimentali

Come si analizza la condizione iniziale se ci trovassimo di fronte a dati sperimentali reali? 
Se vi fornisco una distribuzione misurata di pressione $p$ e temperatura $T$ (o densità $ho$), la prima cosa da verificare per poter utilizzare gli invarianti di Riemann è calcolare l'entropia:
$$rac{p}{ho^\gamma} = 	ext{costante}$$
Se il campo è omoentropico, questa quantità è rigorosamente uniforme nello spazio.

Nella realtà sperimentale, i dati non sono espressi da una funzione continua analitica, ma da misure per punti tramite trasduttori e sonde termiche affette da rumore e incertezza. Lo sperimentatore dovrebbe fornire le barre d'errore (le forchette d'incertezza). In assenza di specifiche strumentali, si calcola la **mediana** (o il valore medio) della grandezza $rac{p}{ho^\gamma}$ e si valuta la dispersione massima (*scattering*). 
Se la fluttuazione rispetto alla mediana è contenuta entro una soglia modesta (qualche percento), si assume lecitamente l'ipotesi di campo omoentropico; se la variazione supera una soglia significativa (ad esempio il 10-15%), l'ipotesi omoentropica cade e non è più lecito applicare gli invarianti globali.

Nel nostro problema analitico ideale, la condizione iniziale assegna $u = 0$ e $a = a_0$ costanti su tutto il semiasse $x \ge 0$: il campo è perfettamente omoentropico con:
$$R_1 = a_0, \qquad R_2 = a_0$$

---

### 5. Topologia delle caratteristiche nel piano $(x, t)$ e assenza del fronte posteriore

All'istante iniziale, in ogni punto del condotto le velocità caratteristiche sono:
$$\left(rac{dx}{dt}ight)_{u-a} = -a_0, \qquad \left(rac{dx}{dt}ight)_u = 0, \qquad \left(rac{dx}{dt}ight)_{u+a} = +a_0$$
Dal punto $(0, 0)$ parte la linea caratteristica di testa della famiglia $u+a$, che si propaga con pendenza:
$$rac{dx}{dt} = a_0 \implies x_{FA}(t) = a_0 t$$
Questa retta costituisce il **fronte anteriore** dell'onda. A destra di essa ($x > a_0 t$), tracciando a ritroso le caratteristiche si intercetta unicamente lo stato imperturbato di quiete: è una **zona a zero famiglie**, dove $u = 0$ e $a = a_0$.

**Esiste un fronte posteriore delimitante una zona uniforme a ridosso del pistone?**
La risposta è **no**. Nel moto impulsivo la velocità del pistone era costante ($V_p$), per cui si formava una regione uniforme a ridosso del pistone. Qui la velocità del pistone varia continuamente nel tempo ($\dot{x}_p = 2b t$): ad ogni istante di tempo il pistone assume una velocità diversa, generando continuamente nuove caratteristiche con inclinazioni sempre diverse. Non può quindi esistere una zona a zero famiglie a ridosso del pistone: il campo è costituito da un'**onda semplice di espansione** su tutta la regione compresa tra la traiettoria del pistone $x_p(t)$ e il fronte anteriore $x_{FA}(t) = a_0 t$.

Poiché le caratteristiche della famiglia opposta $u-a$ provengono tutte dalla regione imperturbata di monte dove $a = a_0$ e $u = 0$, l'invariante $R_1$ è costante su tutto il campo:
$$R_1 = a - \delta u = a_0 \implies a = a_0 + \delta u$$
Trattandosi di un'onda semplice, le caratteristiche della famiglia attiva $u+a$ sono **semirette** a pendenza costante: ciascuna semiretta trasporta valori costanti di $u$ e $a$. Poiché si tratta di un'espansione, man mano che il tempo cresce la velocità del gas diventa sempre più negativa e la velocità del suono decresce, per cui le caratteristiche divergono.

---

### 6. Calcolo analitico: determinazione del tempo di emissione $t_p$

Fissiamo un generico istante di osservazione $ar{t} > 0$ e consideriamo una posizione spaziale $x$ compresa nell'intervallo:
$$x_p(ar{t}) \le x \le a_0 ar{t}$$
Per determinare lo stato $u(x, ar{t})$ e $a(x, ar{t})$, tracciamo a ritroso la semiretta caratteristica della famiglia $u+a$ passante per il punto $(x, ar{t})$. Questa caratteristica incontra la traccia del pistone in un punto di emissione $(x_p(t_p), t_p)$ all'istante incognito $t_p \in [0, ar{t}]$.

L'equazione della semiretta caratteristica che congiunge $(x_p(t_p), t_p)$ con $(x, ar{t})$ è:
$$rac{x - x_p(t_p)}{ar{t} - t_p} = u_p + a_p$$
Dalla cinematica del pistone e dall'invariante di Riemann:
$$x_p(t_p) = b t_p^2$$
$$u_p(t_p) = \dot{x}_p(t_p) = 2b t_p$$
$$a_p(t_p) = a_0 + \delta u_p = a_0 + 2\delta b t_p$$
Sommando la velocità del fluido e la velocità del suono sul pistone:
$$u_p + a_p = a_0 + 2(1 + \delta) b t_p$$
Sostituendo nell'equazione della caratteristica:
$$x - b t_p^2 = \left[ a_0 + 2(1 + \delta) b t_p ight] (ar{t} - t_p)$$

Sviluppando il prodotto algebrico:
$$x - b t_p^2 = a_0 ar{t} - a_0 t_p + 2(1 + \delta) b ar{t} t_p - 2(1 + \delta) b t_p^2$$
Riorganizzando i termini rispetto alle potenze di $t_p$, si ottiene un'**equazione algebrica di secondo grado** in $t_p$:
$$\left[ 2(1 + \delta) - 1 ight] b \, t_p^2 - \left[ a_0 + 2(1 + \delta) b ar{t} ight] t_p + (a_0 ar{t} - x) = 0$$
Ricordando che $2(1 + \delta) - 1 = 2\left(1 + rac{\gamma - 1}{2}ight) - 1 = \gamma$:
$$\gamma b \, t_p^2 - \left( a_0 + (\gamma + 1) b ar{t} ight) t_p + (a_0 ar{t} - x) = 0$$

Risolvendo l'equazione di secondo grado si ottengono due radici algebriche. Il criterio fisico di selezione impone che il tempo di emissione debba appartenere rigorosamente all'intervallo temporale:
$$0 \le t_p \le ar{t}$$
* Quando $x = a_0 ar{t}$ (sul fronte anteriore), il termine noto si annulla ($a_0 ar{t} - x = 0$) e si ottiene $t_p = 0$.
* Quando $x = x_p(ar{t})$ (sul pistone), l'equazione fornisce $t_p = ar{t}$.

Una volta calcolato univocamente $t_p(x, ar{t})$, lo stato del fluido nel punto $(x, ar{t})$ è immediatamente noto:
$$u(x, ar{t}) = 2b t_p(x, ar{t})$$
$$a(x, ar{t}) = a_0 + 2\delta b t_p(x, ar{t})$$

---

### 7. Gradiente di velocità sul fronte anteriore

Rispondiamo ora alla domanda sulla presenza di una discontinuità di gradiente sul fronte anteriore ($x = a_0 ar{t}$).
Poiché $u$ dipende da $x$ attraverso la funzione implicita $t_p(x)$:
$$rac{\partial u}{\partial x} = rac{du_p}{dt_p} rac{\partial t_p}{\partial x} = \ddot{x}_p(t_p) rac{\partial t_p}{\partial x} = 2b rac{\partial t_p}{\partial x}$$
Differenziando rispetto a $x$ l'equazione della caratteristica:
$$rac{\partial}{\partial x} \left[ x - x_p(t_p) ight] = rac{\partial}{\partial x} \left[ (u_p + a_p)(ar{t} - t_p) ight]$$
$$1 - \dot{x}_p rac{\partial t_p}{\partial x} = rac{d(u_p + a_p)}{dt_p} rac{\partial t_p}{\partial x} (ar{t} - t_p) - (u_p + a_p) rac{\partial t_p}{\partial x}$$
Sul fronte anteriore si ha $t_p = 0$, dove $\dot{x}_p(0) = 0$ e $u_p + a_p = a_0$. Ricavando $rac{\partial t_p}{\partial x}$ in $t_p = 0$:
$$1 = \left[ rac{d(u_p + a_p)}{dt_p} ar{t} - a_0 ight] \left( rac{\partial t_p}{\partial x} ight)_{x = a_0 ar{t}}$$
Poiché $rac{d(u_p + a_p)}{dt_p} = 2(1+\delta)b = (\gamma + 1)b$:
$$\left( rac{\partial t_p}{\partial x} ight)_{FA} = rac{1}{(\gamma + 1)b ar{t} - a_0}$$
Di conseguenza, la derivata spaziale della velocità sul fronte anteriore vale:
$$\left( rac{\partial u}{\partial x} ight)_{FA} = rac{2b}{(\gamma + 1)b ar{t} - a_0}$$
Poiché nello stato imperturbato imperturbato a monte ($x > a_0 ar{t}$) la velocità è identicamente nulla ($\partial u / \partial x = 0$), mentre immediatamente a valle del fronte anteriore la derivata spaziale assume il valore finito calcolato sopra, **sul fronte anteriore permane una discontinuità di gradiente (discontinuità debole)**.

---

### 8. Conclusioni didattiche ed estensione alla compressione

Questo procedimento analitico mostra come il modello degli invarianti consenta di risolvere in forma chiusa problemi con contorni mobili non lineari. 
Lo stesso apparato matematico è direttamente applicabile al caso opposto in cui $b > 0$, corrispondente a una **compressione parabolica**: in quel caso le caratteristiche convergono, consentendo di determinare analiticamente il tempo e lo spazio di coalescenza dove si formerà l'onda d'urto.
