# 29 — Lezione 30, 17/11/2021 (file 28): metodo all'esame; schlieren di un proiettile; supersonico 2D, Prandtl–Meyer

## Detto esplicitamente — esame
- Le prove "sono calibrate per testare fin dove ciascuno studente riesce a spingersi"; procedere "dal semplice al
  complesso ed **esporre chiaramente ogni deduzione fisica parziale**, senza rinunciare a scrivere solo perché non si
  è ancora raggiunta la soluzione conclusiva". → Uno svolgimento deve mostrare i passaggi di ragionamento anche
  incompleti.

## Interpretazione di immagini — schlieren di un proiettile che attraversa un cilindro con 4 asole
Ordine di descrizione del docente:
1. **Scia turbolenta (wake)** dietro la poppa tronca: separazione allo spigolo a 90°, vortici tipo von Kármán.
2. **Onde sferiche concentriche non stazionarie** dalle asole: ogni foro è una sorgente acustica impulsiva (sasso
   nello stagno).
3. **Fasce rettilinee oblique stazionarie dal vertice** = **cono di Mach**, inviluppo dei fronti emessi.
- Proprietà: niente segnale a monte; perturbazioni confinate nel cono.
- **Stima del Mach**: triangolo V t / a t ⇒ sin α = 1/M; angolo misurato col goniometro ≈ 27°–30° ⇒ M ≈ 2,0–2,2,
  V ≈ 700 m/s con a ≈ 340 m/s. (Esattamente le richieste di **B55-E03**, B57-E02: "stimare il numero di Mach".)

## Supersonico stazionario 2D
- Nel riferimento del corpo: Eulero stazionario, A_x invertibile se M>1 ⇒ **x fa da "tempo"**, sistema iperbolico in x.
- Caratteristiche: dy/dx = tan θ (doppia, linea di corrente: h_t e s costanti), tan(θ+α) = C⁺, tan(θ−α) = C⁻;
  compatibilità ±(√(M²−1)/(ρV²)) dp − dθ = 0.
- Omoentropico: θ + ν costante lungo C⁻ (tan(θ−α)), θ − ν lungo C⁺ (tan(θ+α)) — **= formulario §10**.
- ν_max = (π/2)(√((γ+1)/(γ−1)) − 1) = 130,45° (aria): getto sonico nel vuoto si apre al massimo di 130,45°;
  con M1 > 1 resta ν_max − ν(M1).

## Convenzioni — differenze con il formulario
- Angolo di Mach: docente **α**, formulario **μ**. Stessa associazione famiglie/invarianti.
