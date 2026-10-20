# Verifica di livello 3: risultati del corso ricalcolati con `gasdyn`

Generato da `verifica/verifica_catalogo.py` a partire da `verifica/casi_corso.py` (gli stessi casi sono eseguiti da `tests/test_gasdyn.py::TestLivello3`).

**Metodo.** Per ogni svolgimento con valori numerici (campo *Associazione e contenuto dello svolgimento* di `fonti/catalogo.md`, più i numeri letti direttamente sulle pagine) si ricalcola la grandezza con il modulo usando i dati della traccia. La tolleranza relativa è scelta in base alle cifre riportate dalla fonte (valori a 2 cifre: 1–2%; valori a 5–6 cifre: 10⁻⁵–10⁻⁴). Ogni discrepanza è stata ricondotta a una causa controllando la pagina originale (rasterizzata da `fonti/esercizi_2003_2023.pdf` o da `fonti/soluzioni_A/`). Pagine: 'G p.' = pagina fisica della raccolta, 'A p.' = pagina del PDF in `fonti/soluzioni_A/`.

**Esito complessivo:** 257 confronti — concorda 211, discorda per causa nota 34, non confrontabile 12. Nessuna discrepanza è attribuibile al modulo.

## Tabella dei confronti

| ID | Grandezza | Svolgimento | Calcolato | Diff. rel. | Esito | Nota / causa | Pagine |
|---|---|---|---|---|---|---|---|
| B61-E01 | a1 dopo la chiusura | 0.2 | 0.2 | 9.7e-16 | concorda |  | A p.2–3 (G p.236–237) |
| B61-E01 | p1 | 1.28000e-5 | 1.28000e-5 | 4.9e-15 | concorda |  | A p.2–3 (G p.236–237) |
| B61-E01 | a = u nel punto sonico del ventaglio | 0.25 | 0.25 | 2.2e-16 | concorda |  | A p.2–3 (G p.236–237) |
| B61-E01 | w urto riflesso dall'ambiente | -2.45 | -2.45341 | 1.4e-03 | concorda |  | A p.2–3 (G p.236–237) |
| B61-E01 | u2 | -2.03 | -2.03093 | 4.6e-04 | concorda |  | A p.2–3 (G p.236–237) |
| B61-E01 | p2 | 0.00224 | 0.00224505 | 2.3e-03 | concorda |  | A p.2–3 (G p.236–237) |
| B61-E01 | a2 | 1.1 | 1.09916 | 7.7e-04 | concorda |  | A p.2–3 (G p.236–237) |
| B61-E01 | p3 (contorno, pressione totale = p_amb) | 0.00227 | 0.00224505 | 1.1e-02 | concorda | lo svolgimento accetta p2 ≠ p3 entro il 2% (prima iterazione); il modulo impone p2 = p3 | A p.2–3 (G p.236–237) |
| B61-E01 | a3 | 0.42 | 0.418412 | 3.8e-03 | concorda |  | A p.2–3 (G p.236–237) |
| B61-E01 | w urto riflesso dal setto | 0.88 | 0.828674 | 5.8e-02 | **discorda** | lo svolgimento legge Ms = 2,51 dalle tabelle e lo combina con a2 ≈ 1,16 (w = 0,88); con lo stato 2 calcolato (u2 = −2,03, a2 = 1,10) le RH danno Ms = 2,60, w = 0,83. Prima approssimazione dichiarata ('errore inferiore al 2%') | A p.2–3 (G p.236–237) |
| B61-E01 | Ms urto riflesso dal setto | 2.51 | 2.60163 | 3.7e-02 | **discorda** | lo svolgimento legge Ms = 2,51 dalle tabelle e lo combina con a2 ≈ 1,16 (w = 0,88); con lo stato 2 calcolato (u2 = −2,03, a2 = 1,10) le RH danno Ms = 2,60, w = 0,83. Prima approssimazione dichiarata ('errore inferiore al 2%') | A p.2–3 (G p.236–237) |
| B61-E01 | p4 | 0.0165 | 0.017354 | 5.2e-02 | **discorda** | lo svolgimento legge Ms = 2,51 dalle tabelle e lo combina con a2 ≈ 1,16 (w = 0,88); con lo stato 2 calcolato (u2 = −2,03, a2 = 1,10) le RH danno Ms = 2,60, w = 0,83. Prima approssimazione dichiarata ('errore inferiore al 2%') | A p.2–3 (G p.236–237) |
| B61-E01 | a4 | 1.63 | 1.64507 | 9.2e-03 | concorda |  | A p.2–3 (G p.236–237) |
| B63-E02 | Ms per M assoluto a valle = 1 | 2.07 | 2.06809 | 9.2e-04 | concorda | valore esatto 2,0681; lo svolgimento si ferma a 2,07 'dopo alcuni tentativi' | A p.2–3 (G p.242–244) |
| B63-E02 | ρ2/ρ1 a Ms = 2,07 | 2.7689 | 2.76895 | 1.7e-05 | concorda |  | A p.2–3 (G p.242–244) |
| B63-E02 | a2/a1 a Ms = 2,07 | 1.3211 | 1.32106 | 2.8e-05 | concorda |  | A p.2–3 (G p.242–244) |
| B63-E02 | M2 a Ms = 2,07 | 1.001 | 1.00103 | 3.0e-05 | concorda |  | A p.2–3 (G p.242–244) |
| B63-E01 | M da Δθ = 18°, σ = 28° | 4.8 | 4.78758 | 2.6e-03 | concorda |  | A p.2–3 (G p.242–244) |
| B63-E01 | ν(4,8) [°] | 75 | 74.9863 | 1.8e-04 | concorda |  | A p.2–3 (G p.242–244) |
| B63-E01 | M da ν = 105° | 11 | 11.091 | 8.3e-03 | concorda |  | A p.2–3 (G p.242–244) |
| B63-E01 | p/p0 a M = 4,8 | 0.0024 | 0.00239426 | 2.4e-03 | concorda |  | A p.2–3 (G p.242–244) |
| B63-E01 | p02/p01 urto normale M = 10 | 0.003 | 0.00304475 | 1.5e-02 | concorda |  | A p.2–3 (G p.242–244) |
| B63-E01 | p01/p1 a M = 10 | 4.35000e4 | 4.24392e4 | 2.4e-02 | **discorda** | valore di tabella arrotondato nello svolgimento (esatto 4,244·10⁴, −2,4%) | A p.2–3 (G p.242–244) |
| B63-E01 | p/p1 nella bolla di ricircolo | 0.312 | 0.309379 | 8.4e-03 | concorda | concorda nonostante l'arrotondamento di p01/p1 | A p.2–3 (G p.242–244) |
| B67-E01 | M medio con A/A* = 1,11 | 1.39 | 1.39094 | 6.7e-04 | concorda |  | A p.3–4 (G p.251–252) |
| B67-E01 | θmax(1,39) [°] (svolgimento: 9–10°) | 9.5 | 9.15273 | 3.7e-02 | concorda | la fonte dà un intervallo 9–10° | A p.3–4 (G p.251–252) |
| B67-E01 | M a monte dell'urto σ = 51°, θ = 22° | 2.2 | 2.21119 | 5.1e-03 | concorda |  | A p.3–4 (G p.251–252) |
| B67-E01 | M da linea di Mach a 35° | 1.75 | 1.74345 | 3.7e-03 | concorda |  | A p.3–4 (G p.251–252) |
| B67-E01 | x_sonica Sauer–Meyer [mm], versione 'favini' | 8.5 | 8.52803 | 3.3e-03 | concorda | R = 82,5 mm, semialtezza 60 mm; distanza gola–punto sonico sull'asse = d/6 (scalata). Torna con la versione senza (γ+1) (S p.65 'Favini') | A p.3–4 (G p.251–252) |
| B67-E01 | x_sonica Sauer–Meyer [mm], versione formulario | 8.5 | 5.50482 | 3.5e-01 | **discorda** | discrepanza 10 di Fase 0: il formulario (con γ+1) dà 5,5 mm; lo svolgimento del 2022 usa la versione senza (γ+1) | A p.3–4 (G p.251–252) |
| B68-E01 | b1) w incidente (M valle = 1) | 2.45 | 2.44699 | 1.2e-03 | concorda |  | A p.3–4 (G p.256–257) |
| B68-E01 | b1) w riflesso | -1.2 | -1.19708 | 2.4e-03 | concorda |  | A p.3–4 (G p.256–257) |
| B68-E01 | b1) p/p1 al ristagno | 17 | 16.7492 | 1.5e-02 | concorda | valori arrotondati | A p.3–4 (G p.256–257) |
| B68-E01 | b2) M valle = 1,89 (oltre il limite 1,8898) | nessuna soluzione | nessuna soluzione | — | non confr. | la fonte riconosce il caso limite e usa 1,87; il modulo solleva NessunaSoluzione | A p.3–4 (G p.256–257) |
| B68-E01 | b2) w incidente (M valle = 1,87) | 23.66 | 21.281 | 1.0e-01 | **discorda** | vicino al limite M valle è malcondizionato: lo svolgimento usa Ms ≈ 20 dalle tabelle (M valle = 1,8738 ≈ 1,87); il valore esatto per 1,87 è Ms = 17,99 (w = 21,3). Con Ms = 20 il modulo riproduce i valori della fonte (vedi riga seguente) | A p.3–4 (G p.256–257) |
| B68-E01 | b2) w riflesso | -8 | -7.13753 | 1.1e-01 | **discorda** | vicino al limite M valle è malcondizionato: lo svolgimento usa Ms ≈ 20 dalle tabelle (M valle = 1,8738 ≈ 1,87); il valore esatto per 1,87 è Ms = 17,99 (w = 21,3). Con Ms = 20 il modulo riproduce i valori della fonte (vedi riga seguente) | A p.3–4 (G p.256–257) |
| B68-E01 | b2) p/p1 | 3700 | 2969.65 | 2.0e-01 | **discorda** | vicino al limite M valle è malcondizionato: lo svolgimento usa Ms ≈ 20 dalle tabelle (M valle = 1,8738 ≈ 1,87); il valore esatto per 1,87 è Ms = 17,99 (w = 21,3). Con Ms = 20 il modulo riproduce i valori della fonte (vedi riga seguente) | A p.3–4 (G p.256–257) |
| B68-E01 | b2) w incidente con Ms = 20 | 23.66 | 23.6643 | 1.8e-04 | concorda |  | A p.3–4 (G p.256–257) |
| B68-E01 | b2) w riflesso con Ms = 20 | -8 | -7.92755 | 9.1e-03 | concorda |  | A p.3–4 (G p.256–257) |
| B68-E01 | b2) p/p1 con Ms = 20 | 3700 | 3683.62 | 4.4e-03 | concorda |  | A p.3–4 (G p.256–257) |
| B68-E01 | σ a θmax per M → ∞ [°] | 67 | 67.7923 | 1.2e-02 | concorda | 'circa 67°' (esatto 67,79°) | A p.3–4 (G p.256–257) |
| B69-E01 | θ_D(M = 1,3) [°] (soglia 'circa 3,2') | 3.2 | 3.24922 | 1.5e-02 | concorda |  | A p.3 (G p.261) |
| B69-E02 | t_BU = 1/(γ+1) | 0.416667 | 0.416667 | 0.0e+00 | concorda |  | A p.2–3 (G p.262–263) |
| B69-E02 | x_BU = t_BU | 0.416667 | 0.416667 | 0.0e+00 | concorda |  | A p.2–3 (G p.262–263) |
| B69-E02 | p1 stato compresso = γ^(2γ/(γ−1)) | 10.5414 | 10.5414 | 3.4e-16 | concorda |  | A p.2–3 (G p.262–263) |
| B69-E02 | w urto (pseudo-Riemann) | 2.87 | 2.87482 | 1.7e-03 | concorda |  | A p.2–3 (G p.262–263) |
| B69-E02 | u2 = u3 | 2.1 | 2.10581 | 2.8e-03 | concorda |  | A p.2–3 (G p.262–263) |
| B69-E02 | p2 = p3 | 9.44 | 9.47533 | 3.7e-03 | concorda | 'roughly' | A p.2–3 (G p.262–263) |
| B69-E02 | w asintotico (pistone a 2) | 2.76 | 2.76205 | 7.4e-04 | concorda |  | A p.2–3 (G p.262–263) |
| B70-E01 | a1 [m/s] | 343.7 | 343.712 | 3.3e-05 | concorda | lo svolgimento usa a1 = 343,7 m/s, cioè R = 288 J/(kg K) (con R = 287: 343,1 m/s) | A p.3–4 (G p.267–268) |
| B70-E01 | w [m/s] | 1375 | 1374.85 | 1.1e-04 | concorda |  | A p.3–4 (G p.267–268) |
| B70-E01 | ρ2/ρ1 | 4.5714 | 4.57143 | 6.3e-06 | concorda |  | A p.3–4 (G p.267–268) |
| B70-E01 | u2 [m/s] | 1074 | 1074.1 | 9.2e-05 | concorda |  | A p.3–4 (G p.267–268) |
| B70-E01 | a2 [m/s] | 691.4 | 691.439 | 5.7e-05 | concorda |  | A p.3–4 (G p.267–268) |
| B70-E01 | M2 | 1.55 | 1.55342 | 2.2e-03 | concorda |  | A p.3–4 (G p.267–268) |
| B70-E01 | p2 [Pa] | 1.85000e6 | 1.85000e6 | 0.0e+00 | concorda |  | A p.3–4 (G p.267–268) |
| B70-E01 | ν(M2) [°] | 13.58 | 13.4825 | 7.2e-03 | **discorda** | ν(1,55) da tabella; esatto ν(1,5534) = 13,48° (−0,1°), effetto trascurabile sui risultati | A p.3–4 (G p.267–268) |
| B70-E01 | M3 | 3.43 | 3.43501 | 1.5e-03 | concorda |  | A p.3–4 (G p.267–268) |
| B70-E01 | p3 [Pa] | 1.05900e5 | 1.05599e5 | 2.8e-03 | concorda |  | A p.3–4 (G p.267–268) |
| B70-E01 | M4 | 4.74 | 4.75 | 2.1e-03 | concorda |  | A p.3–4 (G p.267–268) |
| B70-E01 | p4 [Pa] | 2.00000e4 | 1.86664e4 | 6.7e-02 | concorda | valore dato a una cifra (0,2·10⁵) | A p.3–4 (G p.267–268) |
| B70-E01 | Ms piede B | 1.6 | 1.6 | 0.0e+00 | concorda |  | A p.3–4 (G p.267–268) |
| B70-E01 | p5 [Pa] | 2.82000e5 | 2.82000e5 | 2.1e-16 | concorda |  | A p.3–4 (G p.267–268) |
| B70-E01 | M5 = u5/a5 | 0.77 | 0.689658 | 1.0e-01 | **discorda** | le RH con Ms = 1,6 danno u5/a5 = 0,690; il '0,77' del manoscritto non è coerente con p5 = 2,82·10⁵ (che invece torna): probabile errore di trascrizione | A p.3–4 (G p.267–268) |
| B70-E01 | a3 [m/s] | 458 | 459.314 | 2.9e-03 | concorda |  | A p.3–4 (G p.267–268) |
| B70-E01 | |V3| [m/s] | 1571 | 1577.75 | 4.3e-03 | concorda |  | A p.3–4 (G p.267–268) |
| B70-E01 | M3 relativo a RC | 2.23 | 2.2377 | 3.5e-03 | concorda |  | A p.3–4 (G p.267–268) |
| B70-E01 | p dietro RC lato 3 [Pa] | 6.00000e5 | 5.99296e5 | 1.2e-03 | concorda |  | A p.3–4 (G p.267–268) |
| B70-E01 | a4 [m/s] | 358 | 358.588 | 1.6e-03 | concorda |  | A p.3–4 (G p.267–268) |
| B70-E01 | |V4| [m/s] | 1697 | 1703.29 | 3.7e-03 | concorda |  | A p.3–4 (G p.267–268) |
| B70-E01 | M4 relativo a RC | 3.2 | 3.21638 | 5.1e-03 | concorda |  | A p.3–4 (G p.267–268) |
| B70-E01 | p dietro RC lato 4 [Pa] | 3.20000e5 | 2.22179e5 | 3.1e-01 | **discorda** | con M = 3,22 e p4 = 0,187·10⁵ (o 0,2·10⁵) le RH danno 2,2–2,4·10⁵ Pa: il '3,2·10⁵' del manoscritto sembra ripetere il Mach relativo (3,2) | A p.3–4 (G p.267–268) |
| B71-E01 | M2 | 1.59 | 1.59145 | 9.1e-04 | concorda |  | A p.2 (G p.273) |
| B71-E01 | p2/p1 | 1.8 | 1.82383 | 1.3e-02 | concorda | arrotondato a 2 cifre (esatto 1,824) | A p.2 (G p.273) |
| B71-E01 | M3 | 1.98 | 1.98327 | 1.7e-03 | concorda |  | A p.2 (G p.273) |
| B71-E01 | p3/p1 | 0.99 | 1.00407 | 1.4e-02 | concorda | esatto 1,004 | A p.2 (G p.273) |
| B71-E01 | M4 | 2.76 | 2.76843 | 3.1e-03 | concorda |  | A p.2 (G p.273) |
| B71-E01 | |θ4| [°] | 19 | 19.1665 | 8.8e-03 | concorda |  | A p.2 (G p.273) |
| B71-E02 | x_p(t_v = 0,5) | 0.13 | 0.132614 | 2.0e-02 | concorda |  | A p.2 (G p.274) |
| B71-E02 | u dopo la compressione | 0.61 | 0.61231 | 3.8e-03 | concorda |  | A p.2 (G p.274) |
| B71-E02 | p dopo la compressione | 2.22 | 2.24492 | 1.1e-02 | concorda | esatto 2,245 | A p.2 (G p.274) |
| B71-E02 | a dopo la compressione | 1.12 | 1.12246 | 2.2e-03 | concorda |  | A p.2 (G p.274) |
| B71-E02 | w urto riflesso dalla parete | -0.93 | -0.936132 | 6.6e-03 | concorda |  | A p.2 (G p.274) |
| B71-E02 | Mach relativo | 1.38 | 1.37951 | 3.6e-04 | concorda |  | A p.2 (G p.274) |
| B71-E02 | t_s arresto pistone | 1.36 | 1.36245 | 1.8e-03 | concorda |  | A p.2 (G p.274) |
| B71-E02 | x_s | 0.66 | 0.6607 | 1.1e-03 | concorda |  | A p.2 (G p.274) |
| B72-E02 | w | 18 | 17.9857 | 7.9e-04 | concorda |  | A p.2 (G p.277) |
| B72-E02 | u dietro l'urto | 15 | 14.9418 | 3.9e-03 | concorda |  | A p.2 (G p.277) |
| B72-E02 | p dietro l'urto | 380 | 377.235 | 7.3e-03 | concorda |  | A p.2 (G p.277) |
| B72-E02 | ρ dietro l'urto | 8 | 8.27214 | 3.4e-02 | concorda |  | A p.2 (G p.277) |
| B72-E02 | u nella sezione stretta | 12 | 10.4777 | 1.3e-01 | **discorda** | la massa impone ρ5u5 = 1,5 ρ2u2 ⇒ u5 = 10,5 con ρ5 = 17,7; il '12' non è coerente con ρ5 ≈ 18 e p5 ≈ 1100 (che tornano) né con u* ≈ 16 del Riemann successivo (che torna): refuso | A p.2 (G p.277) |
| B72-E02 | p nella sezione stretta | 1100 | 1093.78 | 5.7e-03 | concorda |  | A p.2 (G p.277) |
| B72-E02 | ρ nella sezione stretta | 18 | 17.6948 | 1.7e-02 | concorda |  | A p.2 (G p.277) |
| B72-E02 | configurazione | RCS | RCS | — | non confr. | coincide (RCS) | A p.2 (G p.277) |
| B72-E02 | u* del Riemann | 16 | 16.1482 | 9.3e-03 | concorda |  | A p.2 (G p.277) |
| B72-E02 | p* del Riemann | 440 | 440.249 | 5.7e-04 | concorda |  | A p.2 (G p.277) |
| B72-E02 | ρ a sinistra del contatto | 9 | 9.23712 | 2.6e-02 | concorda |  | A p.2 (G p.277) |
| B72-E02 | ρ a destra del contatto | 8 | 8.2902 | 3.6e-02 | concorda |  | A p.2 (G p.277) |
| B61-E02 | — | — | — | — | non confr. | svolgimento qualitativo o simbolico: nessun numero da ricalcolare | A |
| B67-E02 | — | — | — | — | non confr. | svolgimento qualitativo o simbolico: nessun numero da ricalcolare | A |
| B68-E02 | — | — | — | — | non confr. | svolgimento qualitativo o simbolico: nessun numero da ricalcolare | A |
| B69-E03 | — | — | — | — | non confr. | svolgimento qualitativo o simbolico: nessun numero da ricalcolare | A |
| B72-E01 | — | — | — | — | non confr. | svolgimento qualitativo o simbolico: nessun numero da ricalcolare | A |
| B70-E02 | M stazionario di K(M) (escluso nel manoscritto) | 0.408248 | 0.408248 | 1.5e-08 | concorda | controllo della derivazione simbolica | A p.2 (G p.270) |
| B01-E02 | θ_D = ν(2,5)/2 [°] | 19.56 | 19.5618 | 9.1e-05 | concorda |  | G p.2–3 |
| B01-E02 | h_R/h_L = A/A*(2,5) | 2.63672 | 2.63672 | 4.7e-07 | concorda |  | G p.2–3 |
| B01-E02 | L/h_L (modello dello svolgimento) | 4.16638 | 4.16638 | 1.1e-06 | concorda | formula dello svolgimento: caratteristica d'uscita rettilinea inclinata di μ_u | G p.2–3 |
| B01-E02 | L/h_L con il metodo delle caratteristiche (n = 40) | 4.16638 | 4.57757 | 9.9e-02 | **discorda** | ipotesi diversa: lo svolgimento assume caratteristiche rette anche nella zona a due famiglie; il MOC del modulo (ugello_minima_lunghezza, semialtezza → altezza) dà L/h ≈ 4,58 | G p.2–3 |
| B02-E01 | Me | 2.5141 | 2.5141 | 1.4e-06 | concorda |  | G p.5–6 |
| B02-E01 | p* [Pa] | 1.63718e6 | 1.63718e6 | 1.2e-07 | concorda |  | G p.5–6 |
| B02-E01 | F_C/A_g [N/m²] | 3.66422e6 | 3.66422e6 | 8.1e-13 | concorda |  | G p.5–6 |
| B02-E01 | ΔF = (F_CD − F_C)/F_CD [%] | 14.9363 | 14.936 | 2.5e-05 | concorda | la riga dati dello svolgimento scrive p0 = 1,013·10⁵ Pa ma i conti usano 3 MPa; denominatore F_CD come nello svolgimento | G p.5–6 |
| B03-E01 | m2 [kg] (modello dello svolgimento) | 0.1911 | 0.191119 | 1.0e-04 | concorda |  | G p.13 |
| B03-E01 | m2 [kg] (F §17, p1 = 10⁵ Pa) | 0.1911 | 0.181205 | 5.2e-02 | **discorda** | lo svolgimento usa h_uscita = cp T1 costante (forma del riempimento) e 1 atm = 10⁵ Pa; il formulario (F §17, svuotamento adiabatico) dà m/m_i = (p/p_i)^(1/γ): 0,1810 kg con 10⁵ Pa, 0,1834 kg con 1 atm = 101325 Pa | G p.13 |
| B03-E01 | m2 [kg] (F §17, p1 = 1 atm) | 0.1911 | 0.183605 | 3.9e-02 | **discorda** | lo svolgimento usa h_uscita = cp T1 costante (forma del riempimento) e 1 atm = 10⁵ Pa; il formulario (F §17, svuotamento adiabatico) dà m/m_i = (p/p_i)^(1/γ): 0,1810 kg con 10⁵ Pa, 0,1834 kg con 1 atm = 101325 Pa | G p.13 |
| B04-E01 | u⁰ (innesco del formulario) | 0.918948 | 0.918947 | 1.2e-06 | concorda |  | G p.17–18 |
| B04-E01 | p2 alla prima iterata | 0.306767 | 0.306767 | 1.2e-06 | concorda |  | G p.17–18 |
| B04-E01 | p3 alla prima iterata | 0.30041 | 0.300409 | 2.7e-06 | concorda |  | G p.17–18 |
| B04-E01 | u* a convergenza | 0.927454 | 0.927453 | 1.5e-06 | concorda |  | G p.17–18 |
| B04-E01 | p* a convergenza | 0.303131 | 0.30313 | 2.7e-06 | concorda |  | G p.17–18 |
| B04-E01 | w urto | 1.75216 | 1.75216 | 2.4e-06 | concorda |  | G p.17–18 |
| B04-E01 | a2 | 0.997729 | 0.997725 | 3.6e-06 | concorda |  | G p.17–18 |
| B04-E01 | a3 | 1.26411 | 1.26411 | 2.8e-06 | concorda |  | G p.17–18 |
| B04-E01 | a5 (riflessione a p = p4) | 1.102 | 1.0789 | 2.1e-02 | **discorda** | lo svolgimento scrive a5 = a3 (p5/p3)^(δ/γ) = 1,102, ma la stessa formula dà 1,0789; anche u5 = u3 + (a3−a5)/δ = 1,876 non è coerente con a5 = 1,102. L'errore si propaga a u6 (−0,5%) e p6 (−1,5%) | G p.17–18 |
| B04-E01 | u5 | 1.876 | 1.85351 | 1.2e-02 | **discorda** | lo svolgimento scrive a5 = a3 (p5/p3)^(δ/γ) = 1,102, ma la stessa formula dà 1,0789; anche u5 = u3 + (a3−a5)/δ = 1,876 non è coerente con a5 = 1,102. L'errore si propaga a u6 (−0,5%) e p6 (−1,5%) | G p.17–18 |
| B04-E01 | u6 = u7 (velocità del contatto) | 1.753 | 1.74444 | 4.9e-03 | **discorda** | lo svolgimento scrive a5 = a3 (p5/p3)^(δ/γ) = 1,102, ma la stessa formula dà 1,0789; anche u5 = u3 + (a3−a5)/δ = 1,876 non è coerente con a5 = 1,102. L'errore si propaga a u6 (−0,5%) e p6 (−1,5%) | G p.17–18 |
| B04-E01 | a6 | 0.833 | 0.834327 | 1.6e-03 | concorda |  | G p.17–18 |
| B04-E01 | p6 | 0.0854 | 0.0866774 | 1.5e-02 | **discorda** | lo svolgimento scrive a5 = a3 (p5/p3)^(δ/γ) = 1,102, ma la stessa formula dà 1,0789; anche u5 = u3 + (a3−a5)/δ = 1,876 non è coerente con a5 = 1,102. L'errore si propaga a u6 (−0,5%) e p6 (−1,5%) | G p.17–18 |
| B05-E01 | p2 (espansione) [Pa] | 3.20580e4 | 3.20577e4 | 9.1e-06 | concorda |  | G p.21–23 |
| B05-E01 | p3 (urto d'arresto) [Pa] | 9.78310e4 | 9.78307e4 | 3.0e-06 | concorda |  | G p.21–23 |
| B05-E01 | configurazione interazione | SCR | SCR | — | non confr. | coincide (SCR) | G p.21–23 |
| B05-E01 | u4 = u5 [m/s] | -3.15845 | -3.15889 | 1.4e-04 | concorda | seconda iterazione dello svolgimento | G p.21–23 |
| B05-E01 | p4 = p5 [Pa] | 9.88990e4 | 9.88996e4 | 6.2e-06 | concorda |  | G p.21–23 |
| B05-E01 | w2 [m/s] | -408.55 | -408.551 | 1.5e-06 | concorda |  | G p.21–23 |
| B05-E01 | w3 urto riflesso dal pistone [m/s] | 405.731 | 406.023 | 7.2e-04 | concorda | la fonte usa a* da relazioni approssimate (diff. 0,07%) | G p.21–23 |
| B05-E01 | p6 [Pa] | 9.99800e4 | 9.99785e4 | 1.5e-05 | concorda |  | G p.21–23 |
| B06-E02 | σ [°] | 27.3827 | 27.3827 | 2.3e-08 | concorda |  | G p.25–27 |
| B06-E02 | M2 (ventre) | 2.505 | 2.505 | 2.7e-07 | concorda |  | G p.25–27 |
| B06-E02 | p2/p1 | 2.05447 | 2.05447 | 1.0e-06 | concorda |  | G p.25–27 |
| B06-E02 | Cp ventre | 0.16738 | 0.167377 | 2.1e-05 | concorda |  | G p.25–27 |
| B06-E02 | M3 (dorso) | 3.57829 | 3.57829 | 1.3e-06 | concorda |  | G p.25–27 |
| B06-E02 | p3/p1 | 0.43114 | 0.431148 | 1.7e-05 | concorda |  | G p.25–27 |
| B06-E02 | Cp dorso | -0.09029 | -0.090294 | 4.5e-05 | concorda |  | G p.25–27 |
| B06-E02 | CL | 0.25376 | 0.253756 | 1.6e-05 | concorda |  | G p.25–27 |
| B06-E02 | CD di forma | -0.04474 | 0.044744 | 2.0e+00 | **discorda** | errore di segno nello svolgimento: CD_forma = (Cp_ventre − Cp_dorso) sinα = +0,04474; anche l'attrito proiettato va sommato con segno + (CD totale ≈ +0,0523, non −0,05232) | G p.25–27 |
| B07-E01 | a5 elio [m/s] | 998.719 | 998.719 | 2.5e-10 | concorda |  | G p.28–29 |
| B07-E01 | M_R5 | 9.02427 | 9.02427 | 3.1e-07 | concorda |  | G p.28–29 |
| B07-E01 | M_R4 | 0.4423 | 0.442295 | 1.1e-05 | concorda |  | G p.28–29 |
| B07-E01 | w [m/s] | -9012.71 | -9012.71 | 3.1e-07 | concorda |  | G p.28–29 |
| B07-E01 | a4 [m/s] | 4894.89 | 4894.89 | 4.5e-07 | concorda |  | G p.28–29 |
| B07-E01 | u4 [m/s] | -6847.7 | -6847.73 | 3.7e-06 | concorda |  | G p.28–29 |
| B07-E01 | T4 [K] | 6126.45 | 7206.45 | 1.8e-01 | **discorda** | lo svolgimento scrive T4/T5 = 20,4215 ma a4/a5 = 4,90117 (che torna) implica T4/T5 = 24,02: lettura di tabella errata; T4 corretto 7206 K | G p.28–29 |
| B07-E01 | u* (prima approssimazione SCS) | -3849.07 | -4054.34 | 5.3e-02 | **discorda** | lo svolgimento si ferma alla prima approssimazione con γ medio = 1,5 e dichiara necessarie una decina di iterate; soluzione esatta SCS: u* = −4054 m/s, p* = 2,313·10⁷ Pa | G p.28–29 |
| B08-E02 | w | 5 | 5 | 0.0e+00 | concorda |  | G p.33 |
| B08-E02 | M1 dietro l'urto | 3.271 | 3.27165 | 2.0e-04 | concorda |  | G p.33 |
| B08-E02 | a1 dietro l'urto | 1.299 | 1.29904 | 2.9e-05 | concorda |  | G p.33 |
| B08-E02 | A/A*(3,271) | 5.47669 | 5.47669 | 2.6e-07 | concorda |  | G p.33 |
| B08-E02 | A/A*(3) | 4.23457 | 4.23457 | 5.0e-07 | concorda |  | G p.33 |
| B08-E02 | A/A*(2,26) | 2.11535 | 2.11535 | 6.6e-07 | concorda |  | G p.33 |
| B09-E01 | θmax(3) [°] | 34.0734 | 34.0734 | 6.6e-09 | concorda |  | G p.34 |
| B09-E01 | σ per θ = 34° | 63.6732 | 63.6732 | 6.9e-08 | concorda |  | G p.34 |
| B09-E01 | M_B | 1.00295 | 1.00295 | 1.4e-06 | concorda |  | G p.34 |
| B09-E01 | (S_C − S∞)/R | 1.11369 | 1.11369 | 3.4e-06 | concorda |  | G p.34 |
| B09-E01 | (S_A − S_B)/R con θ = 34° | 0.264005 | 0.264009 | 1.4e-05 | concorda | la fonte approssima la deviazione sonica (34,008°) con 34° | G p.34 |
| B09-E01 | (S_A − S_B)/R con il punto sonico esatto | 0.264005 | 0.262205 | 6.8e-03 | **discorda** | approssimazione dichiarata: con M_B = 1 esatto la differenza è 0,2622 (−0,7%) | G p.34 |
| B36-E01 | t_BU [s] | 312.26 | 312.259 | 3.3e-06 | concorda | x_BU (lasciato in formula nella fonte) = 39213.9 m | G p.89–90 |
| B37-E01 | a4 aria a 25 °C [m/s] | 346 | 346.117 | 3.4e-04 | concorda |  | G p.93–94 |
| B37-E01 | u2 = u3 [m/s] (440–445) | 442.5 | 410.446 | 7.2e-02 | **discorda** | la relazione d'urto usata per tentativi è errata: (−w)(u3−w) = 2a0²/(γ+1) senza il termine δw² in a*² (relazione di Prandtl nel riferimento dell'urto); con le RH esatte e gas diversi (SF6 / aria a 25 °C) u* = 410 m/s, p* = 4,22 atm. Lo svolgimento assume anche T_amb = 25 °C, assente nella traccia | G p.93–94 |
| B37-E01 | p2 = p3 [atm] | 3.2 | 4.21882 | 3.2e-01 | **discorda** | la relazione d'urto usata per tentativi è errata: (−w)(u3−w) = 2a0²/(γ+1) senza il termine δw² in a*² (relazione di Prandtl nel riferimento dell'urto); con le RH esatte e gas diversi (SF6 / aria a 25 °C) u* = 410 m/s, p* = 4,22 atm. Lo svolgimento assume anche T_amb = 25 °C, assente nella traccia | G p.93–94 |
| B38-E02 | u2 | 1.66667 | 1.66667 | 2.0e-06 | concorda |  | G p.96–97 |
| B38-E02 | p2 | 0.05858 | 0.0585277 | 8.9e-04 | concorda |  | G p.96–97 |
| B38-E02 | a3 = u3 (stato sonico all'interfaccia) | 0.83333 | 0.833333 | 4.0e-06 | concorda |  | G p.96–97 |
| B38-E02 | p3 | 0.2791 | 0.279082 | 6.6e-05 | concorda |  | G p.96–97 |
| B38-E02 | Δt (CFL 0,8) | 0.343 | 0.342857 | 4.2e-04 | concorda |  | G p.96–97 |
| B38-E02 | ρ^(n+1) | 1.23947 | 1.23925 | 1.8e-04 | concorda | flusso di massa arrotondato nella fonte | G p.96–97 |
| B38-E02 | variazione di s = p/ρ^γ [%] | -0.16 | 2.48423 | 1.7e+01 | **discorda** | nel flusso d'energia lo svolgimento usa ρe = p/(γ−1) + ρu² (manca ½): flusso 1,13365 invece di 0,97678. Con il flusso corretto s cresce del 2,5% (lo schema di Godunov aumenta l'entropia, S p.42): la conclusione qualitativa 'non conserva l'entropia' resta valida | G p.96–97 |
| B40-E03 | ν(3) [°] | 49.76 | 49.7573 | 5.3e-05 | concorda |  | G p.100–101 |
| B40-E03 | θ con p_urto = p_compressione [°] | 24.17 | 33.9939 | 4.1e-01 | **discorda** | lo svolgimento usa ν2 = ν1 + θ (segno d'espansione) per la compressione e inverte il rapporto di pressione; con ν2 = ν1 − θ la curva di compressione isentropica incontra la polare solo vicino a θmax, sul ramo debole, a θ = 33,99° (p/p1 = 8,26). La fonte si interrompe ('Mi fermo qui') | G p.100–101 |
| B41-E03 | u6 = u7 [m/s] | 1118.2 | 1118.16 | 3.8e-05 | concorda |  | G p.104–109 |
| B41-E03 | w urto [m/s] | 1405 | 1405.81 | 5.8e-04 | concorda |  | G p.104–109 |
| B41-E03 | p6 = p7 [atm] | 25.44 | 25.452 | 4.7e-04 | concorda |  | G p.104–109 |
| B41-E03 | u* a x = L/3 (onde quasi nulle) | 0 | 0.00735111 | 7.4e-03 | concorda | il primo Riemann (1 / 2) dà onde di intensità ~1e-5: p* = 0,99997 atm | G p.104–109 |
| B42-E02 | γ (ipotesi isentropica) | 1.71 | 1.70951 | 2.9e-04 | concorda |  | G p.112 |
| B42-E02 | R [J/(kg K)] con γ = 1,71 | 2296.32 | 2296.56 | 1.0e-04 | concorda |  | G p.112 |
| B42-E02 | massa molare [g/mol] | 3.62 | 3.6204 | 1.1e-04 | concorda |  | G p.112 |
| B42-E02 | T a valle dell'onda semplice [K] | 400 | 400 | 2.8e-16 | concorda |  | G p.112 |
| B43-E01 | p0/p a M = 0,3 (incomprimibile) | 1.063 | 1.063 | 0.0e+00 | concorda |  | G p.114–115 |
| B43-E01 | p0/p a M = 0,3 (isentropico) | 1.064 | 1.06443 | 4.0e-04 | concorda |  | G p.114–115 |
| B43-E01 | p0/p a M = 0,3 (isotermo) | 1.065 | 1.06503 | 2.5e-05 | concorda |  | G p.114–115 |
| B46-E02 | a1 | 0.661438 | 0.661438 | 2.6e-07 | concorda |  | G p.121–123, 127 |
| B46-E02 | M1 | 12.47 | 12.4728 | 2.3e-04 | concorda |  | G p.121–123, 127 |
| B46-E02 | M_r1 = (u1 − w)/a1 | 2.64575 | 2.64575 | 1.2e-07 | concorda |  | G p.121–123, 127 |
| B46-E02 | p2 | 1 | 1 | 2.2e-16 | concorda |  | G p.121–123, 127 |
| B46-E02 | u2 (seconda versione, p.127) | 7 | 7 | 0.0e+00 | concorda |  | G p.121–123, 127 |
| B46-E02 | u2 (prima versione, p.121) | 2.4286 | 7 | 1.9e+00 | **discorda** | la prima versione applica u2/u1 = ρ1/ρ2 alle velocità assolute invece che a quelle relative; la seconda (u2 = 7) è corretta | G p.121–123, 127 |
| B46-E02 | a2 | 1 | 1 | 1.1e-16 | concorda |  | G p.121–123, 127 |
| B46-E02 | A/A*(7) | 104.14 | 104.143 | 2.7e-05 | concorda |  | G p.121–123, 127 |
| B46-E02 | p02 | 4139.84 | 4139.84 | 3.9e-11 | concorda |  | G p.121–123, 127 |
| B46-E02 | M3 dopo la contrazione | 5.95 | 5.97044 | 3.4e-03 | concorda | Mach letto da tabella | G p.121–123, 127 |
| B46-E02 | u3 | 6.87919 | 6.88166 | 3.6e-04 | concorda |  | G p.121–123, 127 |
| B46-E02 | a3 | 1.15117 | 1.15262 | 1.3e-03 | concorda |  | G p.121–123, 127 |
| B46-E02 | p3 | 2.76438 | 2.70277 | 2.2e-02 | **discorda** | con M3 = 5,95 (tabella) invece di 5,970 la pressione cambia del 2,2%; gli altri valori tornano | G p.121–123, 127 |
| B48-E02 | w urto (famiglia u−a) | 2.99 | 2.99994 | 3.3e-03 | concorda | natura riconosciuta dal modulo: urto u-a | G p.135–136 |
| B49-E01 | t_BU aria | 0.66622 | 0.666225 | 7.8e-06 | concorda |  | G p.138–141 |
| B49-E01 | x_BU aria | 0.78827 | 0.788288 | 2.3e-05 | concorda | a0 = 1,1832 arrotondato nella fonte | G p.138–141 |
| B49-E01 | p_compressione/p1 a u = 1,48 | 4.77288 | 4.7728 | 1.7e-05 | concorda |  | G p.138–141 |
| B49-E01 | w finale (urto con p = p_compressione) | 2.4374 | 2.43462 | 1.1e-03 | concorda | ipotesi dello svolgimento p_urto = p_compressione; M letto da tabella 2,06 (esatto 2,0577 ⇒ w = 2,4346) | G p.138–141 |
| B49-E01 | w finale asintotico (pistone a 1,48) | 2.4374 | 2.36737 | 2.9e-02 | **discorda** | regime asintotico come in B69-E02: urto da pistone impulsivo a v_p = 1,48, w = 2,367. L'ipotesi p_urto = p_compressione dello svolgimento ignora l'onda riflessa (pseudo-Riemann) | G p.138–141 |
| B49-E01 | t_BU, x_BU elio | — | — | — | non confr. | non calcolati nella fonte; modulo: t_BU = 0.65455, x_BU = 0.84332 | G p.138–141 |
| B49-E01 | t_BU, x_BU SF6 | — | — | — | non confr. | non calcolati nella fonte; modulo: t_BU = 0.67566, x_BU = 0.68004 | G p.138–141 |
| B51-E01 | θ slip line [°] | -16.5 | -16.5765 | 4.6e-03 | concorda |  | G p.149–152 |
| B51-E01 | p comune [atm] | 0.07 | 0.0699173 | 1.2e-03 | concorda | la fonte: 0,0699–0,0703 a θ = −16,5° | G p.149–152 |
| B51-E01 | M corrente superiore espansa | 5.63 | 5.64277 | 2.3e-03 | concorda |  | G p.149–152 |
| B51-E01 | σ urto nella corrente inferiore [°] | 24.3 | 24.3835 | 3.4e-03 | concorda |  | G p.149–152 |
| B51-E01 | θ modello isentropico [°] (fra −16 e −17) | -16.5 | -16.0514 | 2.7e-02 | concorda |  | G p.149–152 |
| B51-E01 | stato miscelato (bilanci globali) | — | — | — | non confr. | svolgimento solo impostato (p.152); modulo: ramo supersonico M = 4.0334, p = 0.2569 atm; ramo subsonico M = 0.4341, p = 4.8323 atm | G p.149–152 |
| B52-E02 | u_v del minimo di a0,v (a_m = 1) | -0.833333 | -0.833333 | 4.9e-11 | concorda | prima versione; la seconda versione ha il segno opposto (errore) | G p.161–163 |
| B52-E02 | a0,v minimo = a_m/√(1+δ) | 0.912871 | 0.912871 | 1.2e-16 | concorda |  | G p.161–163 |
| B53-E01 | M2 | 6.65 | 6.65726 | 1.1e-03 | concorda |  | G p.165–167 |
| B53-E01 | p2 [atm] | 7.0868 | 7.07489 | 1.7e-03 | concorda |  | G p.165–167 |
| B53-E01 | M3 dopo la compressione a 30° | 4.32 | 4.26231 | 1.3e-02 | **discorda** | valori 'provvisori' da tabelle: esatto M3 = 4,262 (−1,3%) | G p.165–167 |
| B53-E01 | θ dell'interazione [°] (34–35) | 34.5 | 34.7399 | 7.0e-03 | concorda |  | G p.165–167 |
| B53-E01 | p4 = p5 [atm] | 57 | 58.5543 | 2.7e-02 | concorda | valore provvisorio della fonte | G p.165–167 |
| B53-E01 | M5 | 4.63 | 4.68652 | 1.2e-02 | concorda |  | G p.165–167 |
| B53-E05 | M2 (sezione larga) | 8.15 | 8.16134 | 1.4e-03 | concorda |  | G p.175–176 |
| B53-E05 | u2 | 7.08682 | 7.08727 | 6.5e-05 | concorda |  | G p.175–176 |
| B53-E05 | p2 | 0.376139 | 0.372413 | 9.9e-03 | concorda | M2 da tabella (8,15 invece di 8,161) | G p.175–176 |
| B53-E05 | p3 dietro l'urto | 4.5 | 4.5 | 0.0e+00 | concorda |  | G p.175–176 |
| B53-E05 | u3 | 8.25 | 8.25 | 0.0e+00 | concorda |  | G p.175–176 |
| B53-E05 | a3 | 1.29904 | 1.29904 | 2.5e-10 | concorda |  | G p.175–176 |
| B53-E05 | u* Riemann (3 | 2) | 8.68 | 8.66977 | 1.2e-03 | concorda |  | G p.175–176 |
| B53-E05 | p* Riemann (3 | 2) | 2.82 | 2.81907 | 3.3e-04 | concorda |  | G p.175–176 |
| B53-E05 | p* con la variazione d'area applicata allo stato 3 | 2.82 | 1.88417 | 3.3e-01 | **discorda** | ipotesi diversa: lo svolgimento risolve il Riemann fra lo stato 3 (sezione stretta) e lo stato 2 senza far passare 3 attraverso l'allargamento; con il collegamento stazionario (come in B72-E02) la configurazione diventa SCS con p* = 1,88 | G p.175–176 |
| B53-E08 | ν(M1) [°] con M1 = 5,95 | 84.61 | 84.6085 | 1.7e-05 | concorda | il catalogo trascrive M1 = 5,45, ma il manoscritto ha u1 = 5,95 con a1 = 1 e tutti i conti (p01 = 1499,8; ν = 84,61°) corrispondono a M1 = 5,95: la cifra '4' è un '9'. Deviazione usata: 19° | G p.181 |
| B53-E08 | p01 | 1499.8 | 1499.8 | 3.1e-10 | concorda |  | G p.181 |
| B53-E08 | M2 dopo l'espansione di 19° | 10.5 | 10.4985 | 1.4e-04 | concorda |  | G p.181 |
| B53-E08 | p2 | 0.0255085 | 0.0255325 | 9.4e-04 | concorda | M2 arrotondato a 10,50 | G p.181 |
| B53-E08 | a2 | 0.540844 | 0.592164 | 9.5e-02 | **discorda** | refuso nel manoscritto: lo stesso foglio dà u2 = a2·M2 = 6,2117 = 10,50 × 0,5916, coerente con a2 = 0,5916 calcolato dal modulo | G p.181 |
| B53-E08 | u2 = a2 M2 | 6.21165 | 6.21685 | 8.4e-04 | concorda |  | G p.181 |
| B53-E08 | σ urto di ritorno a θ = 0 [°] | 24.47 | 24.4797 | 4.0e-04 | concorda |  | G p.181 |
| B53-E08 | p3/p2 | 21.91 | 21.9125 | 1.1e-04 | concorda |  | G p.181 |
| B53-E08 | p3 | 0.559 | 0.55948 | 8.6e-04 | concorda |  | G p.181 |
| B53-E14 | a2 dietro l'espansione | 0.6 | 0.6 | 1.9e-16 | concorda |  | G p.204–207 |
| B53-E14 | w urto d'arresto | 0.541661 | 0.541641 | 3.7e-05 | concorda |  | G p.204–207 |
| B53-E14 | M relativo | -4.23607 | -4.23607 | 7.6e-08 | concorda |  | G p.204–207 |
| B53-E14 | p3 | 0.58138 | 0.58138 | 1.7e-07 | concorda |  | G p.204–207 |
| B53-E14 | a3 | 1.2623 | 1.26226 | 2.8e-05 | concorda |  | G p.204–207 |
| B54-E02 | M_C | 6 | 6.00001 | 9.5e-07 | concorda |  | G p.210–212 |
| B54-E02 | p_C | 0.0011989 | 0.0011989 | 8.7e-07 | concorda |  | G p.210–212 |
| B54-E02 | p_D | 0.0501538 | 0.0501541 | 6.4e-06 | concorda |  | G p.210–212 |
| B54-E02 | a_D | 1.07795 | 1.07798 | 2.6e-05 | concorda |  | G p.210–212 |
| B54-E02 | u_D | 0.435681 | 0.435677 | 8.3e-06 | concorda |  | G p.210–212 |
| B54-E02 | M_C' | 6.75 | 6.76859 | 2.8e-03 | concorda | Mach da tabella (esatto 6,769) | G p.210–212 |
| B54-E02 | p_C' | 0.00238223 | 0.00233838 | 1.8e-02 | **discorda** | conseguenza del Mach di tabella 6,75 invece di 6,769 (−1,8%) | G p.210–212 |
| B54-E02 | configurazione | SCS | SCS | — | non confr. | coincide (SCS) | G p.210–212 |
| B54-E02 | p comune | 0.063 | 0.0630578 | 9.2e-04 | concorda |  | G p.210–212 |
| B54-E02 | u comune | 0.6149 | 0.614991 | 1.5e-04 | concorda |  | G p.210–212 |

## Cause delle discrepanze (riepilogo)

- **Arrotondamenti e letture di tabella**: B46-E02 (M3), B54-E02 (M_C'), B63-E01 (p01/p1), B70-E01 (ν), B53-E01 (valori provvisori), B68-E01 b2 (Ms ≈ 20 vicino al limite, malcondizionato).
- **Refusi senza conseguenze sui risultati finali**: B53-E08 (a2), B70-E01 (M5, p dietro RC lato 4), B72-E02 (u5), B07-E01 (T4/T5), B02-E01 (p0 nella riga dati).
- **Errori che si propagano**: B04-E01 (a5 ⇒ u6 −0,5%, p6 −1,5%), B06-E02 (segno di CD), B37-E01 (relazione d'urto senza δw²), B38-E02 (energia cinetica senza ½), B40-E03 (segno di ν nella compressione), B46-E02 prima versione (rapporto di densità sulle velocità assolute).
- **Ipotesi o modelli diversi da quelli del formulario**: B01-E02 (caratteristiche rette per L), B03-E01 (svuotamento con h = cp T1 e 1 atm = 10⁵ Pa), B49-E01 ('w finale' con p_urto = p_compressione, invece del regime asintotico usato in B69-E02), B53-E05 (Riemann senza collegamento stazionario attraverso l'allargamento), B67-E01 (Sauer–Meyer senza (γ+1)), B70-E01 (R = 288).
- **Soluzioni dichiaratamente parziali**: B07-E01 (prima approssimazione), B61-E01 (prima approssimazione con errore 2%), B09-E01 (deviazione sonica approssimata a 34°: concorda con l'ipotesi).
- **Errore del catalogo, non della fonte**: B53-E08, M1 = 5,95 (il catalogo legge 5,45).

## Svolgimenti senza numeri da confrontare (soluzioni_A)

- B61-E02: integrali dell'onda semplice nel transonico di piccole perturbazioni (simbolico).
- B67-E02: condizione M = cost all'uscita, relazione p_x = −c ρ_x (simbolico).
- B68-E02: condizioni di portata ed entalpia totale in ingresso (simbolico).
- B69-E03: identificazione delle strutture in una fotografia.
- B72-E01: identificazione delle strutture in una fotografia.

## Giudizio per esercizio

| ID | Giudizio | Motivo |
|---|---|---|
| B61-E01 | con riserva | prima approssimazione dichiarata: l'urto riflesso dal setto ha Ms 2,51 invece di 2,60 (w −6%) |
| B63-E01 | affidabile | stime dalla fotografia; un valore di tabella arrotondato (−2,4%) senza effetto sul risultato |
| B63-E02 | affidabile | Ms = 2,07 è l'arrotondamento di 2,0681 |
| B67-E01 | con riserva | Sauer–Meyer nella versione senza (γ+1) (discrepanza 10 di Fase 0) |
| B68-E01 | con riserva | caso b2: valori ottenuti con Ms ≈ 20 invece dell'esatto 17,99 per M valle = 1,87 |
| B69-E01 | affidabile |  |
| B69-E02 | affidabile |  |
| B70-E01 | con riserva | due valori incoerenti (M5 = 0,77; p dietro RC lato 4 = 3,2·10⁵); R = 288 |
| B70-E02 | affidabile | derivazione simbolica, punto stazionario verificato |
| B71-E01 | affidabile |  |
| B71-E02 | affidabile |  |
| B72-E02 | con riserva | refuso u5 = 12 (corretto 10,5); i risultati finali tornano |
| B01-E02 | con riserva | L/h_L con caratteristiche rette (modello semplificato dichiarato) |
| B02-E01 | affidabile | riga dati con p0 errato, conti corretti |
| B03-E01 | non affidabile | modello di svuotamento con h_uscita = cp T1 e 1 atm = 10⁵ Pa |
| B04-E01 | con riserva | errore in a5 che sposta u6 dello 0,5% e p6 dell'1,5% |
| B05-E01 | affidabile |  |
| B06-E02 | con riserva | CL corretto; CD con segno sbagliato |
| B07-E01 | non affidabile | soluzione ferma alla prima approssimazione; T4 letto male (stato 4 ok) |
| B08-E02 | affidabile |  |
| B09-E01 | affidabile | approssimazione dichiarata della deviazione sonica (34°) |
| B36-E01 | affidabile |  |
| B37-E01 | non affidabile | relazione d'urto errata nei tentativi |
| B38-E02 | non affidabile | flusso d'energia senza il fattore ½ (conclusione qualitativa valida) |
| B40-E03 | non affidabile | segno di ν errato nella compressione; svolgimento interrotto |
| B41-E03 | affidabile | verificata la prima fase (u6, w, p6); fasi successive dichiarate provvisorie dalla fonte |
| B42-E02 | affidabile |  |
| B43-E01 | affidabile |  |
| B46-E02 | con riserva | seconda versione corretta (p3 −2,2% per M di tabella); prima versione errata |
| B48-E02 | affidabile |  |
| B49-E01 | con riserva | t_BU e x_BU corretti; 'w finale' con ipotesi p_urto = p_compressione |
| B51-E01 | affidabile | punti a) e b); c) solo impostato |
| B52-E02 | affidabile | solo la prima versione (la seconda ha il segno opposto) |
| B53-E01 | con riserva | valori provvisori da tabelle (1–3%) |
| B53-E05 | con riserva | Riemann senza far passare lo stato 3 per l'allargamento |
| B53-E08 | affidabile | con M1 = 5,95 (trascrizione del catalogo errata); refuso su a2 |
| B53-E14 | affidabile |  |
| B54-E02 | affidabile | Mach di tabella (p_C' −1,8%) senza effetto sul risultato finale |

## Svolgimenti numericamente affidabili

Si possono usare come riferimento numerico (concordano col modulo entro le cifre riportate):

B63-E01, B63-E02, B69-E01, B69-E02, B70-E02, B71-E01, B71-E02, B02-E01, B05-E01, B08-E02, B09-E01, B36-E01, B41-E03, B42-E02, B43-E01, B48-E02, B51-E01, B52-E02, B53-E08, B53-E14, B54-E02.

**Affidabili con riserva** (usabili, ma con la correzione indicata): B61-E01, B67-E01, B68-E01, B70-E01, B72-E02, B01-E02, B04-E01, B06-E02, B46-E02, B49-E01, B53-E01, B53-E05.

**Non affidabili numericamente** (impostazione utile, numeri da non usare): B03-E01, B07-E01, B37-E01, B38-E02, B40-E03.

Gli esercizi del catalogo senza svolgimento non sono verificabili e non compaiono qui.

## Interventi sul modulo nati dalla verifica

- `contorno_serbatoio` reso esatto (urto interno + contatto fra gas interno e gas del serbatoio): necessario per riprodurre B61-E01.
- `sauer_meyer(..., versione='favini')` per riprodurre lo svolgimento B67-E01 senza cambiare il default del formulario.
- `miscelamento` (bilanci integrali, F §17) per il punto c) di B51-E01.
- `riemann_area` per l'urto che attraversa una brusca variazione d'area (B72-E02, B53-E05, B46-E02).
