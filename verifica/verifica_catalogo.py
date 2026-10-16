"""Genera verifica/verifica_catalogo.md dai casi di verifica/casi_corso.py."""
import collections
import pathlib
import sys

QUI = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(QUI))
import casi_corso as cc  # noqa: E402

# Giudizio complessivo per esercizio (motivato nelle note della tabella)
GIUDIZIO = {
    "B61-E01": ("con riserva", "prima approssimazione dichiarata: l'urto riflesso dal setto ha Ms 2,51 invece di 2,60 (w −6%)"),
    "B63-E01": ("affidabile", "stime dalla fotografia; un valore di tabella arrotondato (−2,4%) senza effetto sul risultato"),
    "B63-E02": ("affidabile", "Ms = 2,07 è l'arrotondamento di 2,0681"),
    "B67-E01": ("con riserva", "Sauer–Meyer nella versione senza (γ+1) (discrepanza 10 di Fase 0)"),
    "B68-E01": ("con riserva", "caso b2: valori ottenuti con Ms ≈ 20 invece dell'esatto 17,99 per M valle = 1,87"),
    "B69-E01": ("affidabile", ""),
    "B69-E02": ("affidabile", ""),
    "B70-E01": ("con riserva", "due valori incoerenti (M5 = 0,77; p dietro RC lato 4 = 3,2·10⁵); R = 288"),
    "B70-E02": ("affidabile", "derivazione simbolica, punto stazionario verificato"),
    "B71-E01": ("affidabile", ""),
    "B71-E02": ("affidabile", ""),
    "B72-E02": ("con riserva", "refuso u5 = 12 (corretto 10,5); i risultati finali tornano"),
    "B01-E02": ("con riserva", "L/h_L con caratteristiche rette (modello semplificato dichiarato)"),
    "B02-E01": ("affidabile", "riga dati con p0 errato, conti corretti"),
    "B03-E01": ("non affidabile", "modello di svuotamento con h_uscita = cp T1 e 1 atm = 10⁵ Pa"),
    "B04-E01": ("con riserva", "errore in a5 che sposta u6 dello 0,5% e p6 dell'1,5%"),
    "B05-E01": ("affidabile", ""),
    "B06-E02": ("con riserva", "CL corretto; CD con segno sbagliato"),
    "B07-E01": ("non affidabile", "soluzione ferma alla prima approssimazione; T4 letto male (stato 4 ok)"),
    "B08-E02": ("affidabile", ""),
    "B09-E01": ("affidabile", "approssimazione dichiarata della deviazione sonica (34°)"),
    "B36-E01": ("affidabile", ""),
    "B37-E01": ("non affidabile", "relazione d'urto errata nei tentativi"),
    "B38-E02": ("non affidabile", "flusso d'energia senza il fattore ½ (conclusione qualitativa valida)"),
    "B40-E03": ("non affidabile", "segno di ν errato nella compressione; svolgimento interrotto"),
    "B41-E03": ("affidabile", "verificata la prima fase (u6, w, p6); fasi successive dichiarate provvisorie dalla fonte"),
    "B42-E02": ("affidabile", ""),
    "B43-E01": ("affidabile", ""),
    "B46-E02": ("con riserva", "seconda versione corretta (p3 −2,2% per M di tabella); prima versione errata"),
    "B48-E02": ("affidabile", ""),
    "B49-E01": ("con riserva", "t_BU e x_BU corretti; 'w finale' con ipotesi p_urto = p_compressione"),
    "B51-E01": ("affidabile", "punti a) e b); c) solo impostato"),
    "B52-E02": ("affidabile", "solo la prima versione (la seconda ha il segno opposto)"),
    "B53-E01": ("con riserva", "valori provvisori da tabelle (1–3%)"),
    "B53-E05": ("con riserva", "Riemann senza far passare lo stato 3 per l'allargamento"),
    "B53-E08": ("affidabile", "con M1 = 5,95 (trascrizione del catalogo errata); refuso su a2"),
    "B53-E14": ("affidabile", ""),
    "B54-E02": ("affidabile", "Mach di tabella (p_C' −1,8%) senza effetto sul risultato finale"),
}

NON_CONFRONTABILI_A = {
    "B61-E02": "integrali dell'onda semplice nel transonico di piccole perturbazioni (simbolico)",
    "B67-E02": "condizione M = cost all'uscita, relazione p_x = −c ρ_x (simbolico)",
    "B68-E02": "condizioni di portata ed entalpia totale in ingresso (simbolico)",
    "B69-E03": "identificazione delle strutture in una fotografia",
    "B72-E01": "identificazione delle strutture in una fotografia",
}


def fmt(x):
    if x is None:
        return "—"
    if isinstance(x, str):
        return x
    if x == 0:
        return "0"
    a = abs(x)
    if a >= 1e4 or a < 1e-3:
        return f"{x:.5e}".replace("e+0", "e").replace("e-0", "e-")
    return f"{x:.6g}"


def main():
    casi = cc.tutti()
    per_id = collections.OrderedDict()
    for c in casi:
        per_id.setdefault(c.id, []).append(c)
    tot = collections.Counter(c.esito for c in casi)
    righe = []
    w = righe.append
    w("# Verifica di livello 3: risultati del corso ricalcolati con `gasdyn`\n")
    w("Generato da `verifica/verifica_catalogo.py` a partire da `verifica/casi_corso.py` "
      "(gli stessi casi sono eseguiti da `tests/test_gasdyn.py::TestLivello3`).\n")
    w("**Metodo.** Per ogni svolgimento con valori numerici (campo *Associazione e contenuto dello "
      "svolgimento* di `fonti/catalogo.md`, più i numeri letti direttamente sulle pagine) si ricalcola la "
      "grandezza con il modulo usando i dati della traccia. La tolleranza relativa è scelta in base alle cifre "
      "riportate dalla fonte (valori a 2 cifre: 1–2%; valori a 5–6 cifre: 10⁻⁵–10⁻⁴). Ogni discrepanza è "
      "stata ricondotta a una causa controllando la pagina originale (rasterizzata da "
      "`fonti/esercizi_2003_2023.pdf` o da `soluzioni_A/`). Pagine: 'G p.' = pagina fisica della raccolta, "
      "'A p.' = pagina del PDF in `soluzioni_A/`.\n")
    w(f"**Esito complessivo:** {len(casi)} confronti — concorda {tot['concorda']}, "
      f"discorda per causa nota {tot['discorda']}, non confrontabile {tot['non confrontabile']}. "
      "Nessuna discrepanza è attribuibile al modulo.\n")
    w("## Tabella dei confronti\n")
    w("| ID | Grandezza | Svolgimento | Calcolato | Diff. rel. | Esito | Nota / causa | Pagine |")
    w("|---|---|---|---|---|---|---|---|")
    for id_, cs in per_id.items():
        for c in cs:
            d = "—" if c.diff is None else f"{c.diff:.1e}"
            esito = {"concorda": "concorda", "discorda": "**discorda**",
                     "non confrontabile": "non confr."}[c.esito]
            nota = c.nota.replace("|", "/")
            w(f"| {c.id} | {c.grandezza} | {fmt(c.fonte)} | {fmt(c.calcolato)} | {d} | {esito} | {nota} | {c.pagine} |")
    w("")
    w("## Cause delle discrepanze (riepilogo)\n")
    w("- **Arrotondamenti e letture di tabella**: B46-E02 (M3), B54-E02 (M_C'), B63-E01 (p01/p1), "
      "B70-E01 (ν), B53-E01 (valori provvisori), B68-E01 b2 (Ms ≈ 20 vicino al limite, malcondizionato).")
    w("- **Refusi senza conseguenze sui risultati finali**: B53-E08 (a2), B70-E01 (M5, p dietro RC lato 4), "
      "B72-E02 (u5), B07-E01 (T4/T5), B02-E01 (p0 nella riga dati).")
    w("- **Errori che si propagano**: B04-E01 (a5 ⇒ u6 −0,5%, p6 −1,5%), B06-E02 (segno di CD), "
      "B37-E01 (relazione d'urto senza δw²), B38-E02 (energia cinetica senza ½), B40-E03 (segno di ν "
      "nella compressione), B46-E02 prima versione (rapporto di densità sulle velocità assolute).")
    w("- **Ipotesi o modelli diversi da quelli del formulario**: B01-E02 (caratteristiche rette per L), "
      "B03-E01 (svuotamento con h = cp T1 e 1 atm = 10⁵ Pa), B49-E01 ('w finale' con p_urto = p_compressione, "
      "invece del regime asintotico usato in B69-E02), B53-E05 (Riemann senza collegamento stazionario "
      "attraverso l'allargamento), B67-E01 (Sauer–Meyer senza (γ+1)), B70-E01 (R = 288).")
    w("- **Soluzioni dichiaratamente parziali**: B07-E01 (prima approssimazione), B61-E01 (prima "
      "approssimazione con errore 2%), B09-E01 (deviazione sonica approssimata a 34°: concorda con l'ipotesi).")
    w("- **Errore del catalogo, non della fonte**: B53-E08, M1 = 5,95 (il catalogo legge 5,45).\n")
    w("## Svolgimenti senza numeri da confrontare (soluzioni_A)\n")
    for k, v in NON_CONFRONTABILI_A.items():
        w(f"- {k}: {v}.")
    w("")
    w("## Giudizio per esercizio\n")
    w("| ID | Giudizio | Motivo |")
    w("|---|---|---|")
    for id_, (g, m) in GIUDIZIO.items():
        w(f"| {id_} | {g} | {m} |")
    w("")
    aff = [k for k, (g, _) in GIUDIZIO.items() if g == "affidabile"]
    ris = [k for k, (g, _) in GIUDIZIO.items() if g == "con riserva"]
    non = [k for k, (g, _) in GIUDIZIO.items() if g == "non affidabile"]
    w("## Svolgimenti numericamente affidabili\n")
    w("Si possono usare come riferimento numerico (concordano col modulo entro le cifre riportate):\n")
    w(", ".join(aff) + ".\n")
    w("**Affidabili con riserva** (usabili, ma con la correzione indicata): " + ", ".join(ris) + ".\n")
    w("**Non affidabili numericamente** (impostazione utile, numeri da non usare): " + ", ".join(non) + ".\n")
    w("Gli esercizi del catalogo senza svolgimento non sono verificabili e non compaiono qui.\n")
    w("## Interventi sul modulo nati dalla verifica\n")
    w("- `contorno_serbatoio` reso esatto (urto interno + contatto fra gas interno e gas del serbatoio): "
      "necessario per riprodurre B61-E01.")
    w("- `sauer_meyer(..., versione='favini')` per riprodurre lo svolgimento B67-E01 senza cambiare il default "
      "del formulario.")
    w("- `miscelamento` (bilanci integrali, F §17) per il punto c) di B51-E01.")
    w("- `riemann_area` per l'urto che attraversa una brusca variazione d'area (B72-E02, B53-E05, B46-E02).")
    (QUI / "verifica_catalogo.md").write_text("\n".join(righe) + "\n", encoding="utf-8")
    print("scritto", QUI / "verifica_catalogo.md", dict(tot))


if __name__ == "__main__":
    main()
