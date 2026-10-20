# Trascrizione della Lezione del 28 ottobre 2021
## Risoluzione di un Esercizio d'Esame: Tubo d'Urto Asimmetrico con Serbatoio, Espansione Centrata e Condizioni Critiche

---

### 1. Introduzione e presentazione del problema d'esame (24 gennaio 2019)

Oggi tiriamo le somme su quanto sviluppato nelle lezioni precedenti affrontando un problema reale tratto dalla prova d'esame del **24 gennaio 2019**. 
*(Viene menzionato il testo di riferimento per i flussi comprimibili del "Potto Project" di Genick Bar-Meir e il celebre saggio di Carlo M. Cipolla sulle leggi fondamentali della stupidità umana, analizzato tramite piano cartesiano).*

#### Testo del problema
Consideriamo un condotto a sezione costante di area $A = 1 \text{ m}^2$ e lunghezza $L = 3 \text{ m}$, inizialmente chiuso ad entrambe le estremità da due setti o diaframmi che vengono rimossi istantaneamente e contemporaneamente all'istante $t = 0$.

* **Condizione iniziale all'interno del condotto ($0 \le x \le L$)**:
  Il condotto è riempito di **aria** ($\gamma = 1.4$, $R = 287.05 \text{ J/(kg K)}$) in condizioni di quiete:
  $$u(x, 0) = 0$$
  $$p(x, 0) = 3.6306 \text{ bar} = 363060 \text{ Pa}$$
  $$T(x, 0) = 293 \text{ K} \implies a_0 = \sqrt{\gamma R T} = \sqrt{1.4 \times 287.05 \times 293} \approx 343.1 \text{ m/s}$$

* **Condizione al contorno a destra ($x = L$)**:
  La sezione di destra, una volta rimossa la membrana, è a contatto con l'**ambiente esterno** alla pressione atmosferica standard:
  $$p_{\text{amb}} = 1.01325 \text{ bar} = 101325 \text{ Pa}$$
  (aria in condizioni standard, $T_{\text{amb}} = 288.15 \text{ K}$).

* **Condizione al contorno a sinistra ($x = 0$)**:
  La sezione di sinistra, una volta rimossa la membrana, è collegata a un **serbatoio di capacità molto elevata** contenente un gas di combustione/propellente con proprietà termodinamiche:
  $$\gamma_s = 1.148, \qquad R_s = 299 \text{ J/(kg K)}$$
  mantenuto a pressione totale e temperatura totale costanti:
  $$p_{0s} = 26.0703 \text{ bar} = 2607030 \text{ Pa}$$
  $$T_{0s} = 500.6593 \text{ K}$$

L'obiettivo è determinare l'evoluzione ondosa completa nel piano spazio-temporale $(x, t)$ dallo sblocco iniziale fino al raggiungimento delle condizioni asintotiche stazionarie. Poiché il serbatoio ha capacità molto elevata rispetto al volume del tubo, durante l'intero transitorio iniziale $p_0$ e $T_0$ si mantengono rigorosamente costanti.

---

### 2. Analisi qualitativa al contorno di destra ($x = L$): espansione centrata

All'istante iniziale $t = 0^+$, sul contorno destro $x = L$ il gas all'interno del condotto si trova alla pressione di $3.6306 \text{ bar}$, mentre l'ambiente esterno è a $1.01325 \text{ bar}$. Poiché $p_{\text{interno}} > p_{\text{amb}}$, il gas del condotto tende a espandersi verso l'ambiente esterno.

Poiché la rimozione del setto è istantanea, dal punto $(L, 0)$ si genera un'**onda di espansione centrata** sulla famiglia di caratteristiche:
$$\frac{dx}{dt} = u - a$$
che entra nel condotto propagandosi da destra verso sinistra.

* **Fronte anteriore dell'espansione**:
  Si propaga verso sinistra nello stato imperturbato di quiete ($u = 0$, $a = a_0$):
  $$\left(\frac{dx}{dt}\right)_{FA} = -a_0 \implies x_{FA}(t) = L - a_0 t$$

* **Attraversamento dell'espansione e invariante di Riemann**:
  Attraversando l'onda sulla famiglia $u-a$, l'invariante associato alla famiglia opposta $u+a$ rimane costante:
  $$R_2 = a + \delta u = \text{costante} = a_0$$
  con $\delta = \frac{\gamma - 1}{2} = 0.2$.
  Se l'espansione si espandesse liberamente fino alla pressione ambiente $p_{\text{amb}}$, la velocità del suono all'uscita per trasformazione isentropica varrebbe:
  $$a_{\text{amb}} = a_0 \left(\frac{p_{\text{amb}}}{p_0}\right)^{\frac{\gamma - 1}{2\gamma}} = 343.1 \left(\frac{1.01325}{3.6306}\right)^{\frac{0.4}{2.8}} \approx 285.8 \text{ m/s}$$
  La velocità indotta verso l'esterno ($u > 0$) risulterebbe:
  $$u_{\text{valle}} = \frac{a_0 - a_{\text{amb}}}{\delta} = \frac{343.1 - 285.8}{0.2} \approx 286.5 \text{ m/s}$$

* **Verifica della pendenza del fronte posteriore e strozzamento sonico**:
  Calcoliamo la pendenza della caratteristica del fronte posteriore dell'espansione:
  $$\left(\frac{dx}{dt}\right)_{FP} = u_{\text{valle}} - a_{\text{amb}} = 286.5 - 285.8 = +0.7 \text{ m/s} > 0$$
  Poiché la pendenza è positiva, **il fronte posteriore dell'espansione punterebbe verso l'esterno del condotto ($x > L$)!**
  Ciò significa che il flusso raggiungerebbe la condizione transonica/supersonica prima di scaricare l'intera caduta di pressione all'interno del tubo. 
  Di conseguenza, la sezione terminale $x = L$ agisce come una gola sonica strozzata: **la condizione al contorno su $x = L$ non può essere la pressione ambiente, ma diventa la condizione critica sonica**:
  $$u(L, t) = a(L, t) \implies u - a = 0$$
  La caratteristica sonica $u - a = 0$ rimane ancorata verticalmente alla sezione terminale $x = L$. L'espansione residua fino alla pressione atmosferica si sviluppa come getto supersonico all'esterno del condotto. 
  Sulla sezione di uscita $x = L$, sostituendo $u = a$ nell'invariante $a + \delta u = a_0$:
  $$(1 + \delta) a^* = a_0 \implies a^* = \frac{a_0}{1 + \delta} = \frac{2}{\gamma + 1} a_0 = \frac{343.1}{1.2} \approx 285.9 \text{ m/s}$$
  $$u^* = a^* \approx 285.9 \text{ m/s}$$
  $$p^* = p_0 \left(\frac{a^*}{a_0}\right)^{\frac{2\gamma}{\gamma - 1}} = 3.6306 \times (0.8333)^{7} \approx 1.01 \text{ bar}$$

---

### 3. Analisi qualitativa al contorno di sinistra ($x = 0$): propagazione dell'urto e discontinuità di contatto

All'istante $t = 0^+$, sul contorno sinistro $x = 0$ il serbatoio si trova a una pressione totale $p_{0s} \approx 26.07 \text{ bar}$, enormemente superiore alla pressione iniziale dell'aria nel tubo ($3.63 \text{ bar}$).
Il gas del serbatoio irrompe istantaneamente nel condotto agendo come un pistone aerodinamico. 

Poiché la compressione è istantanea, dal punto $(0, 0)$ si originano due fronti d'onda distinti che avanzano verso destra:
1. Un'**onda d'urto** sulla famiglia $u+a$, che viaggia con velocità $W > a_0$ nello stato imperturbato dell'aria.
2. Una **discontinuità di contatto** che separa il gas originariamente contenuto nel tubo (aria compressa dall'urto) dal gas fresco immesso dal serbatoio ($\gamma_s = 1.148$).

Nello spazio-tempo compreso tra $x = 0$ e il fronte dell'urto si distinguono quindi tre regioni:
* **Regione A (stato imperturbato)**: aria iniziale in quiete ($u_A = 0$, $p_A = 3.6306 \text{ bar}$, $a_A = 343.1 \text{ m/s}$).
* **Regione B (stato a valle dell'urto)**: aria iniziale compressa e accelerata dal passaggio dell'urto ($p_B$, $u_B$, $\rho_B$).
* **Regione C (gas del serbatoio)**: gas espanso dal serbatoio all'interno del condotto ($p_C$, $u_C$, $\rho_C$).

---

### 4. Bilancio del sistema: equazioni e incognite

Per determinare completamente lo stato nei domini B e C occorre determinare **7 incognite**:
* 3 variabili di stato nella regione B: $(p_B, u_B, \rho_B)$
* 3 variabili di stato nella regione C: $(p_C, u_C, \rho_C)$
* La velocità di propagazione dell'urto: $W$

A fronte di queste 7 incognite, disponiamo di **7 equazioni indipendenti**:

1. **Relazioni di Rankine-Hugoniot attraverso l'urto (3 equazioni)**:
   Legano lo stato A con lo stato B attraverso la velocità dell'urto $W$:
   $$M_{sA} = \frac{W - u_A}{a_A} = \frac{W}{a_A}$$
   $$\frac{p_B}{p_A} = 1 + \frac{2\gamma_A}{\gamma_A + 1} (M_{sA}^2 - 1)$$
   $$\frac{\rho_B}{\rho_A} = \frac{(\gamma_A + 1) M_{sA}^2}{(\gamma_A - 1) M_{sA}^2 + 2}$$
   $$u_B = u_A + \frac{a_A}{\gamma_A} \left( \frac{p_B}{p_A} - 1 \right) \left[ \frac{\frac{2\gamma_A}{\gamma_A + 1}}{\frac{p_B}{p_A} + \frac{\gamma_A - 1}{\gamma_A + 1}} \right]^{1/2}$$

2. **Relazioni di salto sulla discontinuità di contatto (2 equazioni)**:
   Attraverso la superficie di contatto tra il gas A e il gas C, la pressione e la velocità sono continue:
   $$p_C = p_B$$
   $$u_C = u_B$$
   mentre la densità salta ($\rho_C \ne \rho_B$) e l'entropia salta.

3. **Condizioni al contorno all'imbocco dal serbatoio ($x = 0$) (2 equazioni)**:
   Il gas C effluisce dal serbatoio a $p_{0s}$ e $T_{0s}$ nel condotto. Poiché il processo d'imbocco è un'espansione quasi-stazionaria isentropica:
   $$T_{0s} = T_C + \frac{u_C^2}{2 c_{ps}} = \text{costante}$$
   $$p_{0s} = p_C \left(1 + \frac{\gamma_s - 1}{2} M_C^2\right)^{\frac{\gamma_s}{\gamma_s - 1}} = \text{costante}$$

---

### 5. Individuazione della condizione critica all'imbocco

In linea generale, un sistema non lineare di 7 equazioni in 7 incognite richiede una procedura iterativa (ad esempio ipotizzando $W$, ricavando lo stato B da Rankine-Hugoniot, uguagliando $(p_C, u_C)$ e verificando se soddisfano le condizioni di ristagno del serbatoio).

Tuttavia, esaminando i valori numerici del problema, emerge una scorciatoia analitica fondamentale.
Il rapporto di pressione tra serbatoio e condotto è enorme:
$$\frac{p_{0s}}{p_A} = \frac{26.07 \text{ bar}}{3.63 \text{ bar}} \approx 7.18$$
La massima pressione che un flusso espanso da un serbatoio a condotto a sezione costante può sostenere corrisponde alle **condizioni critiche di strozzamento sonico ($M_C = 1$)**:
$$p_s^* = p_{0s} \left( \frac{2}{\gamma_s + 1} \right)^{\frac{\gamma_s}{\gamma_s - 1}}$$
Poiché $\gamma_s = 1.148$:
$$\frac{2}{\gamma_s + 1} = \frac{2}{2.148} \approx 0.9311$$
$$\frac{\gamma_s}{\gamma_s - 1} = \frac{1.148}{0.148} \approx 7.757$$
$$p_s^* = 26.0703 \times (0.9311)^{7.757} \approx 26.0703 \times 0.5702 \approx 14.86 \text{ bar}$$

Poiché la pressione critica del serbatoio ($p_s^* \approx 14.86 \text{ bar}$) è di gran lunga superiore alla pressione nel condotto (sia iniziale $3.63 \text{ bar}$, sia a valle di un urto moderato), **l'imbocco del condotto a $x = 0$ lavora rigorosamente in condizioni critiche soniche ($M_C = 1$)**:
$$p_C = p_s^* \approx 14.86 \text{ bar}$$
$$T_C = T_s^* = T_{0s} \left( \frac{2}{\gamma_s + 1} \right) = 500.6593 \times 0.9311 \approx 466.16 \text{ K}$$
$$a_C = \sqrt{\gamma_s R_s T_s^*} = \sqrt{1.148 \times 299 \times 466.16} \approx 400.0 \text{ m/s}$$
$$u_C = a_C \approx 400.0 \text{ m/s}$$

Essendo lo stato C completamente determinato dalle condizioni critiche dell'imbocco:
$$p_B = p_C = p_s^* \approx 14.86 \text{ bar}$$
$$u_B = u_C \approx 400.0 \text{ m/s}$$

A questo punto il problema dell'urto si disaccoppia completamente: conoscendo $p_A = 3.6306 \text{ bar}$ e $p_B = 14.86 \text{ bar}$, il rapporto di pressione attraverso l'urto è noto:
$$\frac{p_B}{p_A} = \frac{14.86}{3.6306} \approx 4.093$$
Dalla relazione di Rankine-Hugoniot della pressione:
$$\frac{p_B}{p_A} = 1 + \frac{2\gamma_A}{\gamma_A + 1}(M_{sA}^2 - 1) \implies 4.093 = 1 + \frac{2.8}{2.4}(M_{sA}^2 - 1)$$
$$M_{sA}^2 - 1 = 3.093 \times \frac{2.4}{2.8} \approx 2.651 \implies M_{sA} = \sqrt{3.651} \approx 1.911$$
La velocità di propagazione dell'urto è univocamente determinata:
$$W = M_{sA} a_A = 1.911 \times 343.1 \approx 655.6 \text{ m/s}$$

---

### 6. Sintesi dello scenario asintotico

L'evoluzione completa si riassume come segue:
1. L'onda d'urto viaggia a circa $655.6 \text{ m/s}$ verso destra comprimendo l'aria fino a circa $14.86 \text{ bar}$.
2. La superficie di contatto avanza verso destra alla velocità del fluido $u_B = u_C \approx 400 \text{ m/s}$.
3. Da destra, l'espansione centrata sonicamente strozzata avanza verso sinistra con fronte anteriore a $-343.1 \text{ m/s}$.
4. A un istante calcolabile dall'intersezione tra le traiettorie dell'urto e del fronte anteriore dell'espansione, le due onde si scontreranno dando origine a un'ulteriore interazione di rifrazione ondosa non lineare, fino al raggiungimento dello stato asintotico stazionario di svuotamento del serbatoio attraverso il condotto.
