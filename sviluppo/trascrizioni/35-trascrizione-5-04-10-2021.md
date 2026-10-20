# Trascrizione integrale della lezione del 4 ottobre 2021
## Propagazione Ondosa, Variabili d'Onda e Teoria delle Caratteristiche: Formulazione Matriciale

---

### 1. Fenomeni transitori nello spazio: Coronal Mass Ejection (CME) e modelli ondosi

Perché ci sono persone che se ne occupano? Perché non possono fare altro, sono pagati per quello! Si studia cosa arriva verso la Terra... se arriva qualche cosa. Come ha detto lei? «Ondata». Va bene, andato giù pesante! «Friggere» va benissimo, sull'«ondata» sono rimasto un po' così, è un difetto tecnico: «ondata» non fa parte delle parole che si usano nell'ambiente, ma rende benissimo l'idea. 

Associati a questi fenomeni di **Coronal Mass Ejection (CME)** che inducono una perturbazione sul vento solare (che normalmente spira ed è plasma, da cui siamo protetti grazie al campo magnetico terrestre), partono particelle ad altissima energia. Queste particelle bucano il nostro scudo magnetico e possono letteralmente friggere sistemi elettrici e satelliti. E se un astronauta si avventura nello spazio interplanetario, viene fritto rapidamente. Quindi è fondamentale studiarlo.

Come si studia? È un fenomeno transitorio, una **propagazione ondosa**: bisogna capire come va modellato. Non abbiamo misure dirette facili; abbiamo misure indirette e visualizzazioni della superficie solare. Ma la visualizzazione da sola richiede un'interpretazione: cos'è quella macchia? Cos'è quel colore? Senza un modello matematico non riusciamo a interpretarlo. E il modello sono le equazioni, che parlano come l'oracolo di Delfi: bisogna farle parlare in modo tale da capirle.

---

### 2. Chiarimenti sull'entropia, equazione di stato e moto omoentropico

Ho ricevuto una domanda su quanto fatto giovedì scorso. Riguardava la relazione tra pressione, densità ed entropia e le relative equazioni. 

Definiamo l'entropia attraverso la legge di Poisson:

$$s \sim \frac{p}{\rho^\gamma}$$

Questa definizione è comoda e nasconde le costanti del gas ($\gamma$, $R$), pur essendo del tutto equivalente alla formulazione termodinamica classica. Differenziando questa espressione, si ottiene la variazione di entropia in forma differenziale:

$$ds = dp - a^2 d\rho$$

dove la velocità del suono è definita da $a^2 = \frac{\gamma p}{\rho}$. 

Quando scriviamo un differenziale $d$, intendiamo una variazione generica, che può essere una variazione nel tempo, nello spazio, o totale:

$$\frac{\partial s}{\partial t} = \frac{\partial p}{\partial t} - a^2 \frac{\partial \rho}{\partial t}, \quad \nabla s = \nabla p - a^2 \nabla \rho$$

Richiamiamo l'equazione di conservazione della massa in termini di derivata sostanziale:

$$\frac{D\rho}{Dt} = -\rho \nabla \cdot \mathbf{V}$$

Se consideriamo un **processo isentropico lungo la traiettoria della particella**, abbiamo:

$$\frac{Ds}{Dt} = 0 \implies \frac{Dp}{Dt} - a^2 \frac{D\rho}{Dt} = 0$$

Sostituendo l'equazione della massa:

$$\frac{Dp}{Dt} + \rho a^2 \nabla \cdot \mathbf{V} = 0 \quad \Longleftrightarrow \quad \frac{\partial p}{\partial t} + \mathbf{V} \cdot \nabla p + \gamma p \nabla \cdot \mathbf{V} = 0$$

Ribadisco la distinzione terminologica fondamentale:
* **Moto isentropico**: l'entropia si conserva unicamente lungo la traiettoria (se non stazionario) o lungo la linea di corrente (se stazionario): $\frac{Ds}{Dt} = 0$.
* **Moto omoentropico**: se la condizione iniziale presenta un'entropia uniforme in ogni punto dello spazio, l'entropia è costante in tutto il campo ($s = \text{costante}$ ovunque, $\nabla s = 0$).

*(Parentesi didattica del docente sugli appunti, sui refusi nei libri di testo — citando il classico Courant & Friedrichs — e sull'importanza dell'impedenza acustica $\rho a$).*

---

### 3. Definizione generale di onda e concetto di monte e valle

Che cosa significa $ds = 0$ lungo la traiettoria $\frac{d\mathbf{x}}{dt} = \mathbf{V}$?
Significa che:

$$dp - a^2 d\rho = 0$$

Abbiamo trovato che **una combinazione delle variabili primitive si mantiene costante lungo una specifica direzione nello spazio-tempo**. Questa combinazione $dp - a^2 d\rho$ prende il nome di **variabile d'onda**, e la traiettoria è la sua **direzione di propagazione**.

Perché usiamo il termine "onda"? Perché il fenomeno fisico che stiamo osservando è una propagazione. 
Che cos'è un'onda nel senso più generale possibile? È un fenomeno che **si propaga con velocità finita**. 
Se facciamo cadere un sasso in uno specchio d'acqua calma, vediamo i fronti d'onda circolari che avanzano allontanandosi dal punto dell'impatto. In ogni punto, l'arrivo dell'onda altera la quota del pelo libero.

Per registrare e misurare l'effetto di un'onda, dobbiamo passare **da monte a valle**:
* **Monte**: la regione non ancora raggiunta dal fronte d'onda;
* **Valle**: la regione già investita dall'onda, dove lo stato del sistema è cambiato.

Se ci muoviamo esattamente lungo la direzione dell'onda (o cavalchiamo la cresta dell'onda con una tavola da surf), non vediamo alcuna variazione: veniamo semplicemente trasportati, e la variabile d'onda rimane costante! 
Per misurare l'**intensità dell'onda**, dobbiamo **attraversare il fronte d'onda**, muovendoci in direzione trasversale rispetto alla propagazione.

Vogliamo ora cercare sistematicamente: **esistono altre direzioni caratteristiche e altre variabili d'onda oltre alla traiettoria e all'entropia?**

---

### 4. Struttura del problema bidimensionale e coordinate rotanti

Consideriamo un flusso **bidimensionale non stazionario** nel piano $(x, y)$ al variare del tempo $t$. 
*(Nota: la differenza fondamentale tra 2D e 3D sta nella dispersione geometrica dell'energia: nel cerchio $2\pi r$ l'energia decade come $1/r$, nella sfera $4\pi r^2$ decade come $1/r^2$).*

Nello spazio-tempo abbiamo:
* 3 variabili indipendenti: $x, y, t$;
* 4 variabili dipendenti (primitive): $\rho, u, v, p$.

Un fronte d'onda nello spazio bidimensionale possiede una direzione normale individuata da un versore $\mathbf{n} = (\nu_x, \nu_y)$, con $\|\mathbf{n}\| = 1$. Poiché l'onda può propagarsi in qualunque direzione del piano, $\mathbf{n}$ è un parametro continuo che può ruotare di $2\pi$.

Associata a ciascuna possibile onda con velocità di propagazione $\lambda_k$ (con $k = 1, 2, 3, 4$), definiamo la direzione normale nello spazio-tempo:

$$\xi_k = \mathbf{x} \cdot \mathbf{n}_k - \lambda_k t = x \nu_{kx} + y \nu_{ky} - \lambda_k t$$

La coordinata $\xi_k$ misura lo spostamento nella direzione normale che **attraversa il fronte d'onda**. 
Se $\xi_k = \text{costante}$ ($d\xi_k = 0$), ci muoviamo lungo la direzione di propagazione dell'onda:

$$d\mathbf{x} \cdot \mathbf{n}_k - \lambda_k dt = 0 \implies \frac{d\mathbf{x} \cdot \mathbf{n}_k}{dt} = \lambda_k$$

Ipotizziamo che le perturbazioni ondogene siano funzioni di una **variabile d'onda scalare** $\sigma_k$, la quale varia unicamente attraverso il fronte d'onda:

$$\sigma_k = \sigma_k(\xi_k)$$

Lo stato del sistema (il vettore delle variabili primitive $\mathbf{U}$) dipende quindi da $\xi_k$ tramite la variabile d'onda:

$$\mathbf{U} = \mathbf{U}(\sigma_k(\xi_k))$$

L'**intensità dell'onda** è data dalla derivata della variabile d'onda rispetto alla direzione trasversale:

$$\frac{d\sigma_k}{d\xi_k} \ne 0$$

Se questa derivata fosse nulla, non vi sarebbe alcuna onda.

---

### 5. Equazioni di Eulero in variabili primitive e formulazione quasi-lineare

Scriviamo il sistema completo delle equazioni dei gas perfetti bidimensionali non stazionari in forma quasi-lineare per le variabili primitive $\mathbf{U} = (\rho, u, v, p)^T$:

1. **Conservazione della massa**:
   $$\frac{\partial \rho}{\partial t} + u \frac{\partial \rho}{\partial x} + v \frac{\partial \rho}{\partial y} + \rho \frac{\partial u}{\partial x} + \rho \frac{\partial v}{\partial y} = 0$$

2. **Quantità di moto lungo $x$**:
   $$\frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} + \frac{1}{\rho}\frac{\partial p}{\partial x} = 0$$

3. **Quantità di moto lungo $y$**:
   $$\frac{\partial v}{\partial t} + u \frac{\partial v}{\partial x} + v \frac{\partial v}{\partial y} + \frac{1}{\rho}\frac{\partial p}{\partial y} = 0$$

4. **Energia interna / Pressione**:
   $$\frac{\partial p}{\partial t} + u \frac{\partial p}{\partial x} + v \frac{\partial p}{\partial y} + \gamma p \frac{\partial u}{\partial x} + \gamma p \frac{\partial v}{\partial y} = 0$$

In forma compatta matriciale:

$$\frac{\partial \mathbf{U}}{\partial t} + \mathbf{A}_x \frac{\partial \mathbf{U}}{\partial x} + \mathbf{A}_y \frac{\partial \mathbf{U}}{\partial y} = 0$$

dove $\mathbf{U} = \begin{pmatrix} \rho \\ u \\ v \\ p \end{pmatrix}$ e le matrici giacobiane $4 \times 4$ sono:

$$\mathbf{A}_x = \begin{pmatrix}
u & \rho & 0 & 0 \\
0 & u & 0 & 1/\rho \\
0 & 0 & u & 0 \\
0 & \gamma p & 0 & u
\end{pmatrix}, \quad
\mathbf{A}_y = \begin{pmatrix}
v & 0 & \rho & 0 \\
0 & v & 0 & 0 \\
0 & 0 & v & 1/\rho \\
0 & 0 & \gamma p & v
\end{pmatrix}$$

---

### 6. Problema algebrico agli autovalori: autovettori destri e sinistri

Calcoliamo le derivate parziali di $\mathbf{U}$ sfruttando la regola della catena:

$$\frac{\partial \mathbf{U}}{\partial t} = \frac{d\mathbf{U}}{d\sigma_k} \frac{d\sigma_k}{d\xi_k} \frac{\partial \xi_k}{\partial t} = -\lambda_k \frac{d\mathbf{U}}{d\sigma_k} \frac{d\sigma_k}{d\xi_k}$$

$$\frac{\partial \mathbf{U}}{\partial x} = \frac{d\mathbf{U}}{d\sigma_k} \frac{d\sigma_k}{d\xi_k} \frac{\partial \xi_k}{\partial x} = \nu_{kx} \frac{d\mathbf{U}}{d\sigma_k} \frac{d\sigma_k}{d\xi_k}$$

$$\frac{\partial \mathbf{U}}{\partial y} = \frac{d\mathbf{U}}{d\sigma_k} \frac{d\sigma_k}{d\xi_k} \frac{\partial \xi_k}{\partial y} = \nu_{ky} \frac{d\mathbf{U}}{d\sigma_k} \frac{d\sigma_k}{d\xi_k}$$

Sostituendo queste espressioni nel sistema matriciale delle equazioni di Eulero:

$$\left[ -\lambda_k \mathbf{I} + \mathbf{A}_x \nu_{kx} + \mathbf{A}_y \nu_{ky} \right] \frac{d\mathbf{U}}{d\sigma_k} \left( \frac{d\sigma_k}{d\xi_k} \right) = 0$$

Poiché l'intensità dell'onda è non nulla ($\frac{d\sigma_k}{d\xi_k} \ne 0$), possiamo semplificarla, ottenendo un **problema algebrico omogeneo agli autovalori**:

$$\left[ -\lambda_k \mathbf{I} + \mathbf{A}_x \nu_{kx} + \mathbf{A}_y \nu_{ky} \right] \mathbf{R}_k = 0$$

dove $\mathbf{R}_k = \frac{d\mathbf{U}}{d\sigma_k}$ è l'**autovettore destro** (vettore colonna).
L'autovettore destro stabilisce come l'onda agisce sullo stato del sistema:

$$d\mathbf{U} = \mathbf{R}_k \, d\sigma_k$$

Definiamo inoltre l'**autovettore sinistro** $\mathbf{L}_k$ (vettore riga) associato, normalizzato in modo che:

$$\mathbf{L}_j \mathbf{R}_k = \delta_{jk} \quad (\text{delta di Kronecker})$$

Moltiplicando a sinistra la relazione differenziale, l'autovettore sinistro proietta le variazioni delle grandezze di stato sulla variabile d'onda:

$$d\sigma_k = \mathbf{L}_k \, d\mathbf{U}$$

---

### 7. L'equazione caratteristica

Il sistema algebrico omogeneo ammette soluzioni non banali ($\mathbf{R}_k \ne 0$) se e solo se il determinante della matrice dei coefficienti è nullo:

$$\det\left( -\lambda_k \mathbf{I} + \mathbf{A}_x \nu_{kx} + \mathbf{A}_y \nu_{ky} \right) = 0$$

Questa è l'**equazione caratteristica** del sistema. Essendo un sistema di ordine 4, essa fornisce un polinomio di quarto grado nell'incognita $\lambda_k$. 

Le radici $\lambda_k$ rappresentano le **velocità di propagazione delle onde**:
* Esse **devono essere reali** affinché il sistema sia iperbolico e le onde esistano fisicamente nello spazio-tempo reale;
* Dobbiamo disporre di un **insieme completo di 4 autovettori linearmente indipendenti**, per poter descrivere qualunque perturbazione dello stato del gas come combinazione lineare di onde elementari.

Domani calcoleremo esplicitamente il determinante, determineremo le velocità caratteristiche $\lambda_k$ e identificheremo gli autovettori e le corrispondenti variabili d'onda.

*(Fine lezione: ore 12:54).*
