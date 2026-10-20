# 19 — Lezione 19, 28/10/2021 (file 43): esercizio d'esame 24/01/2019 — tubo tra serbatoio e ambiente

Prova del 24/01/2019: **non presente nel catalogo** (il catalogo salta il 2019–2020). Tipologie affini:
B37-E01 (serbatoio di SF6 + condotto + ambiente, condizioni all'apertura e asintotiche), B11-E01, B27-E02,
B35-E01/B36-E02 (serbatoi e transitori), B47-E01/B67-E02 (condizioni al contorno).

## Testo (dalla trascrizione)
Condotto A = 1 m², L = 3 m, chiuso da due setti rimossi a t=0. Dentro: aria in quiete, p = 3,6306 bar, T = 293 K.
A destra: ambiente 1,01325 bar. A sinistra: serbatoio grande, gas γ_s = 1,148, R_s = 299, p0s = 26,0703 bar,
T0s = 500,6593 K (costanti). Chiede l'evoluzione nel piano x–t fino all'asintoto.

## Ragionamento del docente (ordine)
1. **Contorno destro, analisi qualitativa**: p_int > p_amb ⇒ espansione centrata u−a da (L,0); testa
   dx/dt = −a0.
2. Ipotesi di espansione fino a p_amb lungo R2 = a + δu = a0 ⇒ a_amb ≈ 285,8, u ≈ 286,5 ⇒ **verifica della pendenza
   della coda**: u − a = +0,7 > 0 ⇒ la coda uscirebbe dal tubo ⇒ **la condizione al contorno non è p_amb ma la
   condizione sonica u = a** in x = L (uscita strozzata, getto sottoespanso all'esterno).
   Allora a* = a0/(1+δ) = 2a0/(γ+1), p* ≈ 1,01 bar. [Nota: qui il caso è al limite; il numero 3,6306 è scelto
   apposta per stare vicino alla soglia sonica.]
3. **Contorno sinistro**: p0s ≫ p_A ⇒ il gas del serbatoio agisce "come un pistone aerodinamico": urto u+a +
   contatto. Regioni A (aria in quiete), B (aria urtata), C (gas del serbatoio).
4. **Conteggio equazioni/incognite**: 7 incognite (p,u,ρ in B e C, W) = 3 R–H + 2 contatto (p, u continui) + 2 BC
   d'imbocco isentropico quasi-stazionario (T0s, p0s). Procedura generale iterativa (ipotizzare W, …).
5. **Scorciatoia da riconoscere guardando i numeri**: pressione critica dell'imbocco p_s* = p0s(2/(γs+1))^(γs/(γs−1))
   ≈ 14,86 bar ≫ pressioni nel tubo ⇒ **imbocco sonico** (M_C = 1) ⇒ stato C noto (T* ≈ 466 K, u_C = a_C ≈ 400 m/s)
   ⇒ p_B, u_B noti ⇒ urto disaccoppiato: p_B/p_A ≈ 4,09 ⇒ M_s ≈ 1,91, W ≈ 656 m/s.
   [Inferenza da controllare: con imbocco sonico u_C = a_C deve essere coerente anche con u_B da R–H; il docente
   non lo verifica esplicitamente nella trascrizione.]
6. Scenario: urto → destra, contatto a u ≈ 400 m/s, espansione strozzata da destra; poi interazione urto–espansione e
   asintoto stazionario.

## Metodo che emerge (esplicito + inferito)
- Prima analisi **qualitativa ai contorni** (che onde nascono e di che famiglia), poi stati, poi numeri.
- Ogni condizione al contorno va **verificata a posteriori** (pendenza della coda ⇒ sonico; p* ⇒ imbocco critico).
- Contare equazioni e incognite prima di risolvere.
- Cercare nei dati la scorciatoia fisica (condizione critica) invece di iterare alla cieca.
- Gas diversi ai due lati del contatto (γ, R diversi) — attenzione a usare le proprietà giuste per ciascuna regione.
- Cita il testo "Potto Project" (Bar-Meir).
