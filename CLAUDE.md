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
- Commit a ogni passo significativo; aggiornare questo file a fine fase.

## Struttura della cartella

- `gasdyn.py` — il modulo (aree 1–9: termodinamica, quasi-1D/serbatoi/miscelamento, urto normale,
  urti mobili, onde semplici/pistoni/formazione urto, Riemann 1D e contorni, supersonico 2D, grafici, Godunov).
- `tests/test_gasdyn.py` — 82 test di livello 1–2 (classi `TestArea1` … `TestArea89`), tutti verdi.
- `fase0_report.md` — report di Fase 0 (convenzioni, algoritmo di Riemann, 13 discrepanze, struttura).
- `verifica/rasterizza.py` — rasterizza pagine PDF: `python verifica/rasterizza.py <pdf> <dir> [pagine] --scala 2`.
  Immagini nello scratchpad o in `verifica/raster/` (ignorata da git).
- `verifica/figure_esempio.py` → `verifica/figure/*.png` (controllate a vista).
- `verifica/casi_corso.py` — **livello 3**: classe `Caso` (id, grandezza, fonte, calcolato, tol, atteso, nota,
  pagine), una funzione per esercizio, `tutti()` restituisce 257 casi. Eseguibile da solo: stampa la tabella
  e marca `<-- ESITO INATTESO` se esito calcolato ≠ atteso (attualmente nessuno).
- `verifica/verifica_catalogo.py` → genera `verifica/verifica_catalogo.md` (tabella completa, cause, giudizio per
  esercizio, elenco svolgimenti affidabili). Ultimo esito: 211 concorda, 34 discorda per causa nota, 12 non confrontabili.
- `fonti/` — formulario, summary, catalogo, rapporto, `esercizi_2003_2023.pdf` (113 MB, escluso da git).
- `soluzioni_A/A_G_Bxx-Eyy.pdf` — 17 soluzioni 2021–2023 (tutte lette in Fase 2).

## Ambiente

- Python 3.14; venv `.venv` (`.\.venv\Scripts\python.exe`); dipendenze in `requirements.txt`.
- git: `C:\Program Files\Git\cmd\git.exe` (se manca nel PATH: `$env:Path = "C:\Program Files\Git\cmd;" + $env:Path`).
  `core.autocrlf=false`.
- Attenzione: in PowerShell 5.1 `Set-Content -Encoding utf8` aggiunge il BOM e i here-string perdono il
  newline finale: preferire gli strumenti Write/Edit per i file Python.
- Test: `.\.venv\Scripts\python.exe -m pytest -q tests` (~8 s).

## Stato delle fasi

- [x] Fase 0 — confermata dall'utente il 2026-10-16.
- [x] Fase 1 — modulo aree 1–9 completo, test liv. 1–2 verdi, figure di esempio.
- [~] Fase 2 — validazione livello 3. **Fatto**: lette tutte le soluzioni_A e le pagine della raccolta con
  numeri; `casi_corso.py` e `verifica_catalogo.md` generati. **Da fare per chiudere la Fase 2**:
  1. aggiungere a `tests/test_gasdyn.py` una classe `TestLivello3` parametrizzata su `casi_corso.tutti()`
     che asserisce `caso.esito == caso.atteso` (importare con `sys.path` su `verifica/`);
  2. eseguire tutti i test, rigenerare `verifica_catalogo.md`, rileggerlo una volta, commit.
- [ ] Fase 3 — `gasdyn_guida.md` (convenzioni, firme di tutte le funzioni pubbliche con un esempio per area,
  limiti noti), pulizia finale, riepilogo all'utente (cosa fatto, limiti, segnalazioni "fuori formulario").

## Convenzioni (confermate)

- δ = (γ−1)/2; q = (γ−1)/(2γ); A_i = 2/((γ+1)ρ_i), B_i = (γ−1)/(γ+1) p_i (formulario §1, §8).
- Variabile entropica del corso s = p/ρ^γ; entropia termodinamica Δs_th = c_v ln(s2/s1).
- 1D: x verso destra, u > 0 verso destra. λ1 = u−a (R1 = a−δu costante lungo essa), λ2 = u+a (R2 = a+δu),
  λ0 = u. Onda u+a: R1 uniforme; onda u−a: R2 uniforme.
- Riemann: 1 | u−a | 2 | contatto | 3 | u+a | 4; sigle R, S, N (onda nulla), C; vuoto "RVR".
- Urti mobili: M_r = (u−w)/a con segno; 1 = monte, 2 = valle; urto u+a: M_r1 < 0; u−a: M_r1 > 0; Ms = |M_r1|.
- 2D: σ angolo d'urto, θ in gradi antiorario; onda C+ (fronte θ+μ): Δ(θ+ν)=0, compressione ⇒ θ aumenta,
  polare θ2 = θ1+θ(σ); onda C−: θ2 = θ1−θ(σ). Corrente superiore emette C+, inferiore C−.
- Angoli in gradi nell'interfaccia pubblica. `Gas.aria(R=287.0)` è solo una comodità; nessun γ implicito.

## Decisioni prese

- Discrepanze formulario/summary: si segue il formulario (Newton su u per Riemann, Sauer–Meyer con (γ+1)
  come default, ecc.). `sauer_meyer(..., versione='favini')` riproduce la variante senza (γ+1) usata in B67-E01.
- Riemann: Newton del formulario con smorzamento; `brentq` come rete di sicurezza (metodo marcato);
  gas diversi supportati; `modello='isentropico'` disponibile.
- `contorno_serbatoio` esatto (urto + contatto) con opzione `modello='isentropico'` (S p.45).
- `riemann_stazionario` restituisce la soluzione sui rami deboli (default); `intersezioni_polari` dà tutte.
- Fuori formulario (dichiarati nei docstring): pressione totale incomprimibile/isoterma, M_valle_limite,
  urto isotermo, Pitot di Rayleigh, shock-expansion, MOC ugello di minima lunghezza, miscelamento è F §17.

## Risultati notevoli della Fase 2 (dettagli in verifica_catalogo.md)

- Nessuna discrepanza attribuibile al modulo. Cause trovate: arrotondamenti/tabelle, refusi, errori che si
  propagano (B04, B06, B37, B38, B40, B46 v1), ipotesi diverse (B01 L, B03, B49 w finale, B53-E05, B67, B70 R=288).
- B53-E08: errore del catalogo (M1 = 5,95, non 5,45).
- B68-E01 b2: M valle = 1,89 oltre il limite 1,8898 (γ=1,4) → NessunaSoluzione; la fonte usa Ms ≈ 20.

## Problemi aperti

- Nessuno bloccante. Da riportare nel riepilogo finale: le 13 discrepanze di `fase0_report.md` e le funzioni
  fuori formulario elencate sopra.
