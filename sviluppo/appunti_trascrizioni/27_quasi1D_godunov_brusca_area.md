# 27 — Lezione 28, 15/11/2021 (file 29): quasi-1D per Godunov, brusca variazione d'area (avviamento)

## Detto esplicitamente
- Studiare le equazioni **insieme alle immagini** (Van Dyke, *An Album of Fluid Motion*): "senza le immagini la teoria
  rischia di apparire astratta ed enigmatica".
- Perché il quasi-1D per un campo multidimensionale (avviamento di un ugello): è il modello più semplice per il
  **comportamento medio** lungo l'asse; cattura posizione/velocità dell'urto d'avviamento e di ricompressione; **non**
  risolve separazione, scia, urti obliqui, disco di Mach. → giustificazione del modello da scrivere sempre.
- Ripete: intensità variabile (divergente attenua, convergente intensifica), niente onde semplici, l'entropia genera
  acustica (instabilità nei motori; nei liquidi si rompono le simmetrie della piastra d'iniezione).
- CK esteso con splitting anche del termine sorgente (I⁺B + I⁻B).
- Forma a divergenza: (ρA)_t + (ρuA)_x = 0; (ρuA)_t + (ρu²A + pA)_x − p dA/dx = 0; (ρEA)_t + (ρuH A)_x = 0.
  Godunov: flussi dal **RP classico a sezione costante**, sorgente media di cella S_j = (0, −p̄_j (dA/dx)_j, 0).

## Esempio svolto — approssimazione a **brusca variazione d'area** (avviamento di un divergente)
Dati: A1 = 1 → A2 = 2; quiete (stato 4) ρ=1,4, p=1, a=1; urto u+a con W = 3.
1. R–H: M_s = 3 ⇒ p2 = 10,33, ρ2 ≈ 5,40, u2 ≈ 2,22, a2 ≈ 1,637 ⇒ **M2 ≈ 1,36 > 1 anche nel riferimento assoluto**.
2. Quando l'urto supera il salto d'area, il gas a valle dell'urto vi passa: transizione istantanea ⇒ **stazionaria**,
   isentropica, legge delle aree (lecita perché supersonica: nulla risale) ⇒ M ≈ 2,2, stato
   (ρ, u, p) = (1,952; 3,0665; 2,4918).
3. Problema di Riemann tra questo stato (sinistra) e la quiete (destra): tutte le onde verso destra; nuovo urto
   trasmesso e onda di ricompressione (SCS o RCS — nel summary risulta un secondo urto).
- Logica: "concentrare il divergente in un'ascissa" + legge delle aree stazionaria + Riemann = stima rapida senza CFD.

## Differenze con il summary
- Il summary ha "Nozzle Start Up (simplified)" con lo **stesso metodo** ma numerazione diversa degli stati
  (1 = ambiente, 3 = dietro l'urto, 5 = dopo il salto d'area) e conclude che la terza onda è un urto.
  Nella trascrizione: 4 = quiete, 2 = dietro l'urto, 1 = dopo il salto d'area. → Etichette libere, metodo identico.
- Il summary (punto 2 del caso non isentropico) scrive "if A↑, it strengthens up": refuso; la trascrizione dice
  chiaramente che in un divergente l'onda si **attenua**.

## Catalogo
- B72-E02 (urto mobile e variazione di sezione), B08-E02 (urti in condotti a sezione variabile), avviamento ugelli.
