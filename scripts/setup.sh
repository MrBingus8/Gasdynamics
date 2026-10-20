#!/usr/bin/env bash
# Installa l'ambiente per risolvere gli esercizi (Ubuntu/Debian, non interattivo, idempotente):
#   - dipendenze Python di requirements.txt (numpy, scipy, matplotlib, pytest, ...);
#   - LaTeX minimo per template/svolgimento.tex, latexmk;
#   - poppler-utils (pdftoppm) per convertire pagine PDF in immagini.
# Salta ciò che è già installato; usa sudo solo se serve ed è disponibile.
# Uso: bash scripts/setup.sh
set -euo pipefail

RADICE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$RADICE"

info() { printf '[setup] %s\n' "$*"; }

# --- Privilegi per apt -------------------------------------------------------
if [ "$(id -u)" -eq 0 ]; then
    SUDO=""
elif command -v sudo >/dev/null 2>&1 && sudo -n true 2>/dev/null; then
    SUDO="sudo -n"
else
    SUDO="__nessuno__"
fi

# --- Pacchetti di sistema ----------------------------------------------------
PACCHETTI_APT=(
    texlive-latex-recommended
    texlive-latex-extra
    texlive-lang-italian
    texlive-science
    latexmk
    poppler-utils
)

mancanti=()
if command -v dpkg-query >/dev/null 2>&1; then
    for p in "${PACCHETTI_APT[@]}"; do
        if ! dpkg-query -W -f='${Status}' "$p" 2>/dev/null | grep -q "install ok installed"; then
            mancanti+=("$p")
        fi
    done
else
    info "dpkg non disponibile: salto i pacchetti di sistema (non è Debian/Ubuntu?)."
fi

if [ "${#mancanti[@]}" -gt 0 ]; then
    if [ "$SUDO" = "__nessuno__" ]; then
        info "ERRORE: servono privilegi di root per installare: ${mancanti[*]}" >&2
        exit 1
    fi
    info "Installo con apt: ${mancanti[*]}"
    export DEBIAN_FRONTEND=noninteractive
    $SUDO apt-get update -qq
    $SUDO apt-get install -y -qq --no-install-recommends "${mancanti[@]}"
else
    info "Pacchetti di sistema già presenti."
fi

# --- Python ------------------------------------------------------------------
PY="$(command -v python3 || command -v python || true)"
if [ -z "$PY" ]; then
    info "Python assente: installo python3 e pip."
    [ "$SUDO" = "__nessuno__" ] && { info "ERRORE: servono privilegi di root." >&2; exit 1; }
    export DEBIAN_FRONTEND=noninteractive
    $SUDO apt-get update -qq
    $SUDO apt-get install -y -qq --no-install-recommends python3 python3-pip
    PY="$(command -v python3)"
fi

if "$PY" -c "import numpy, scipy, matplotlib, pytest, pypdf, pypdfium2, PIL" 2>/dev/null; then
    info "Dipendenze Python già presenti ($PY)."
else
    if ! "$PY" -m pip --version >/dev/null 2>&1; then
        info "pip assente: provo ensurepip, poi apt."
        "$PY" -m ensurepip --upgrade >/dev/null 2>&1 || {
            [ "$SUDO" = "__nessuno__" ] && { info "ERRORE: pip non installabile senza root." >&2; exit 1; }
            export DEBIAN_FRONTEND=noninteractive
            $SUDO apt-get update -qq
            $SUDO apt-get install -y -qq --no-install-recommends python3-pip
        }
    fi
    info "Installo requirements.txt con pip."
    # Su Python di sistema "externally managed" (PEP 668) serve --break-system-packages.
    "$PY" -m pip install -q -r requirements.txt 2>/dev/null \
        || "$PY" -m pip install -q --break-system-packages -r requirements.txt \
        || "$PY" -m pip install -q --user --break-system-packages -r requirements.txt
fi

# --- Controllo finale --------------------------------------------------------
"$PY" -c "import numpy, scipy, matplotlib; print('[setup] numpy', numpy.__version__, '| scipy', scipy.__version__, '| matplotlib', matplotlib.__version__)"
latexmk -v | head -n 1
command -v pdftoppm >/dev/null && info "pdftoppm: $(command -v pdftoppm)"
kpsewhich italian.ldf siunitx.sty booktabs.sty float.sty >/dev/null && info "Pacchetti LaTeX del template presenti."
info "Ambiente pronto."
