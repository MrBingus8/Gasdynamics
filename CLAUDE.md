# CLAUDE.md — progetto `gasdyn.py`

File di ripartenza: una nuova sessione deve poter riprendere leggendo solo questo file.

## Obiettivo

Costruire e validare `gasdyn.py`, modulo Python di calcolo per il corso di Gasdinamica
(Favini). Sarà usato in un progetto claude.ai per risolvere esercizi d'esame: il
ragionamento fisico resta al modello, il modulo garantisce numeri corretti e risolve
esattamente i calcoli iterativi (polari, Riemann, formazione d'urto).
Notazione, segni, nomi degli stati e algoritmi devono essere quelli delle fonti del corso.

Regole fissate dall'utente:
- Fonte primaria: `fonti/formulario.tex`; seconda: `fonti/summary_lectures.tex`.
- Relazioni non presenti nelle fonti: ammesse, ma dichiarate "fuori formulario" nel docstring e segnalate.
- Disaccordo formulario/summary: non scegliere da soli, elencare la discrepanza.
- Un unico file `gasdyn.py`, dipendenze solo numpy, scipy, matplotlib.
- Docstring e report in italiano. Mai numeri silenziosamente sbagliati: eccezioni chiare o rami etichettati.
- Livello 3 dei test: distinguere "concorda" da "discorda per causa nota".

## Struttura della cartella

- `fonti/formulario.tex`, `fonti/summary_lectures.tex` — teoria (sola lettura).
- `fonti/catalogo.md` — catalogo tracce 2003–2023 (pagine = pagine fisiche di `fonti/esercizi_2003_2023.pdf`).
- `fonti/rapporto.md` — analisi statistica delle tracce (non usata per le formule).
- `fonti/esercizi_2003_2023.pdf` — raccolta completa (113 MB, scansioni).
- `soluzioni_A/A_G_Bxx-Eyy.pdf` — 17 soluzioni impaginate 2021–2023 (B61-E01 … B72-E02).
- `verifica/rasterizza.py` — rasterizza pagine PDF in PNG: `python verifica/rasterizza.py <pdf> <dir> [pagine] --scala 2`.
  Le immagini vanno nello scratchpad o in `verifica/raster/` (ignorata da git).
- Da creare: `gasdyn.py`, `tests/test_gasdyn.py`, `verifica/verifica_catalogo.md`, `verifica/figure/`, `gasdyn_guida.md`.

## Ambiente

- Python 3.14 di sistema; venv in `.venv` (`.\.venv\Scripts\python.exe`). Dipendenze in `requirements.txt`.
- Rasterizzazione PDF con `pypdfium2` (solo per la verifica, non serve al modulo).
- **git non è installato** sul sistema (al 2026-10-16): repository non ancora inizializzato. In attesa di decisione dell'utente.

## Stato delle fasi

- [x] Fase 0 — lettura fonti, convenzioni, algoritmo di Riemann, discrepanze, struttura proposta.
      **In attesa di conferma dell'utente** (non scrivere codice del modulo prima).
- [ ] Fase 1 — modulo per aree (1 termodinamica, 2 quasi-1D, 3 urto normale, 4 urti mobili,
      5 onde semplici, 6 Riemann, 7 supersonico 2D, 8 grafici, 9 Godunov), test liv. 1–2 per area.
- [ ] Fase 2 — validazione livello 3 su catalogo e soluzioni_A, `verifica/verifica_catalogo.md`.
- [ ] Fase 3 — `gasdyn_guida.md`, pulizia, riepilogo.

## Convenzioni (proposte in Fase 0, DA CONFERMARE)

- δ = (γ−1)/2; q = (γ−1)/(2γ); A_i = 2/((γ+1)ρ_i), B_i = (γ−1)/(γ+1) p_i (formulario §1, §8).
- Variabile entropica del corso s = p/ρ^γ; entropia termodinamica Δs_th = c_v ln(s2/s1).
- 1D: asse x verso destra; u > 0 verso destra. Famiglie λ1 = u−a (invariante R1 = a−δu costante lungo essa),
  λ2 = u+a (R2 = a+δu), contatto λ0 = u. Onda u+a: R1 uniforme; onda u−a: R2 uniforme.
- Riemann: stato 1 sinistra, 4 destra; onda u−a tra 1 e 2, contatto tra 2 e 3, onda u+a tra 3 e 4;
  p2 = p3 = p*, u2 = u3 = u*. Configurazioni a tre lettere (sinistra–contatto–destra): R, S, N (onda nulla), C.
- Urti mobili: v = u − w, M_r = (u − w)/a con segno; 1 = monte (davanti all'urto), 2 = valle;
  urto u+a: M_r1 < 0 (w > u1); urto u−a: M_r1 > 0. "Mach dell'urto" = |M_r1|. [q] = q_valle − q_monte.
- 2D: σ = angolo d'urto (non β), θ = deviazione/direzione del flusso, μ = asin(1/M), ν(M) Prandtl–Meyer.
  C+: dy/dx = tan(θ+μ), lungo C+ si conserva θ−ν; C−: tan(θ−μ), lungo C− si conserva θ+ν (formulario §10).
  Attraverso un'onda C+ (fronte inclinato θ+μ) vale Δ(θ+ν)=0: compressione ⇒ θ aumenta.
- Report completo di Fase 0: `fase0_report.md`.
- Pistone: u_p = velocità del pistone; legge x_p(t).

## Decisioni prese (confermate dall'utente il 2026-10-16)

- Tutte le proposte sulle discrepanze di `fase0_report.md` (c) sono confermate: si segue il formulario
  (Newton su u per Riemann, Sauer–Meyer con (γ+1), polari C+ ⇒ θ1+θ(σ), C− ⇒ θ1−θ(σ), ecc.).
- Angoli in **gradi** in tutta l'interfaccia pubblica.
- `Gas.aria(R=287.0, gamma=1.4)` disponibile come comodità, sempre sovrascrivibile; le funzioni richiedono
  comunque il gas esplicitamente (nessun γ implicito).
- Sigla del vuoto nel problema di Riemann: **"RVR"**.
- Includere tutte le aggiunte proposte: serbatoi, contorno a serbatoio, caratteristiche 2D e ugello di
  lunghezza minima, Sauer–Meyer, salti omotermi, natura di un salto tra due stati, diagramma θ–M1 RR/MR.
- git installato con winget (`C:\Program Files\Git\cmd\git.exe`; se non è nel PATH della shell:
  `$env:Path = "C:\Program Files\Git\cmd;" + $env:Path`). `fonti/*.pdf` escluso dal repository (113 MB).
- Livelli di test 1–3 tutti in `tests/test_gasdyn.py`.

## Problemi aperti

- (nessuno per ora)
