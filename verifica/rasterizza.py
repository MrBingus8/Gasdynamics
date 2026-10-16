"""Rasterizza pagine di un PDF in PNG per il controllo visivo degli svolgimenti.

Uso:
    python verifica/rasterizza.py <file.pdf> <cartella_uscita> [pagine] [--scala S]

`pagine` è un elenco 1-based tipo "1,3-5" (default: tutte). Le pagine sono
numerate come pagine fisiche del PDF, coerentemente con fonti/catalogo.md.
"""
import argparse
import pathlib

import pypdfium2 as pdfium


def _pagine(spec, n):
    if not spec:
        return list(range(1, n + 1))
    out = []
    for parte in spec.split(","):
        if "-" in parte:
            a, b = parte.split("-")
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(parte))
    return [p for p in out if 1 <= p <= n]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdf")
    ap.add_argument("uscita")
    ap.add_argument("pagine", nargs="?", default="")
    ap.add_argument("--scala", type=float, default=2.0)
    args = ap.parse_args()

    pdf = pdfium.PdfDocument(args.pdf)
    uscita = pathlib.Path(args.uscita)
    uscita.mkdir(parents=True, exist_ok=True)
    stem = pathlib.Path(args.pdf).stem
    for p in _pagine(args.pagine, len(pdf)):
        img = pdf[p - 1].render(scale=args.scala).to_pil()
        f = uscita / f"{stem}_p{p:03d}.png"
        img.save(f)
        print(f)


if __name__ == "__main__":
    main()
