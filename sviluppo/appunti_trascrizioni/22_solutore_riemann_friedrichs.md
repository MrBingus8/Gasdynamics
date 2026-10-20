# 22 — Lezioni 22–23, 04/11/2021 (file 13): solutore esatto di Riemann, problema di Friedrichs, pseudo-Riemann

## Detto esplicitamente
- Storia: Mach, shadowgraph e schlieren nati per la balistica (Melsens); xilografie delle prime foto.
- Riemann: soluzione autosimile in ξ = x/t, nessuna scala intrinseca; al più 3 onde; stati 2 (sinistra del contatto)
  e 3 (destra); p2 = p3, u2 = u3, ρ2 ≠ ρ3.
- Conteggio: SCS 8 incognite (u,p,ρ in 2 e 3, W_L, W_R) = 3+3 R–H + 2 contatto; RCR 6 = 2+2 invarianti + 2.
- **Algoritmo del docente** (= formulario: Newton su u):
  1. primo tentativo isentropico ("due espansioni"): Z = (p4/p1)^((γ−1)/(2γ)) a1/a4 (indice di non-omoentropia,
     Z = 1 se s1 = s4); u⁽⁰⁾ = [(a1 + δu1) − Z(a4 − δu4)] / [δ(1+Z)];
  2. "**Perché stimare la velocità e non la pressione?**" — a parità di Δu le curve isentropica e R–H differiscono poco in
     u, molto in p ⇒ u⁽⁰⁾ è vicinissimo alla soluzione;
  3. natura delle onde: sinistra R se u⁽⁰⁾ > u1, S se ≤; destra R se u⁽⁰⁾ < u4, S se ≥;
  4. Δp = |p2(u) − p3(u)| < ε; aggiornamento con secante (robusta) o Newton (2–4 iterazioni, "non diverge mai" con
     questo u⁽⁰⁾).
- **Problema di Friedrichs** (blast wave 1D, curva di Friedlander): pistone parte con V_p > 0 (urto u+a), si arresta a
  t_f ⇒ espansione centrata **sulla stessa famiglia u+a** (deve inseguire l'urto e annullare la velocità); le sue
  caratteristiche, più veloci, entrano nell'urto e lo indeboliscono e curvano.
- **Pseudo-problema di Riemann** per l'esito globale senza integrare la zona di interazione: stati A (quiete a monte
  dell'urto) e C (gas fermo sul pistone dopo l'arresto). A→B urto (NCS) con u_B = V_p; B→C espansione (NCR) con u_C = 0;
  poiché l'espansione è "più efficiente" dell'urto, **p_C > p_A**. Pseudo-RP tra C (sinistra) e A (destra) ⇒ **RCS**:
  espansione riflessa u−a verso il pistone, contatto (separa gas urtato forte/debole: salto di entropia), urto residuo.
- Riferimento dell'urto: x' = x − Wt; monte supersonico relativo, valle subsonico relativo; **le caratteristiche
  della famiglia dell'urto entrano nell'urto da entrambi i lati** e determinano W; le altre lo attraversano subendo
  i salti.

## Catalogo
- Pistone avviato e arrestato / urto–espansione stessa famiglia / pseudo-Riemann: B04-E02, B05-E01, B03-E03,
  B52-E03, B53-E14, B56-E02, B57-E01, B69-E02 (la fonte distingue pseudo-Riemann da stato asintotico), B71-E02.
