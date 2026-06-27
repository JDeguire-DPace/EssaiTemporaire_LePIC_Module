import os
import re
from pathlib import Path

# ====== CONFIG ======
ROOT = Path("./")           # Dossier racine à analyser
OUT_MD = Path("routine_headers.md")
EXTS = {".f90", ".f95", ".f03", ".f08", ".f", ".for", ".ftn"}

# Préfixes/attributs possibles avant le mot-clé
PREFIXES = r"(?:pure|impure|elemental|recursive|module)\s+"
START_RE = re.compile(
    rf"^\s*(?:{PREFIXES})*(subroutine|function)\b", re.IGNORECASE
)

# --- utilitaires ---

def strip_inline_comment_outside_quotes(s: str) -> str:
    """
    Retire le commentaire commençant par '!' s'il est hors guillemets.
    Gère '...' et "..." (sans échapper complexe).
    """
    out = []
    in_sq = False
    in_dq = False
    for i, ch in enumerate(s):
        if ch == "'" and not in_dq:
            in_sq = not in_sq
        elif ch == '"' and not in_sq:
            in_dq = not in_dq
        if ch == "!" and not in_sq and not in_dq:
            break
        out.append(ch)
    return "".join(out)

def paren_balance(s: str) -> int:
    """Balance () hors guillemets, pour savoir si on a fermé tous les ()"""
    bal = 0
    in_sq = False
    in_dq = False
    for ch in s:
        if ch == "'" and not in_dq:
            in_sq = not in_sq
        elif ch == '"' and not in_sq:
            in_dq = not in_dq
        elif not in_sq and not in_dq:
            if ch == "(":
                bal += 1
            elif ch == ")":
                bal -= 1
    return bal

def is_fixed_form(fn: Path) -> bool:
    return fn.suffix.lower() in {".f", ".for", ".ftn"}

def logical_lines_free_form(lines):
    """
    Recolle les lignes avec & (fin de ligne) et & (début de la suivante).
    Conserve aussi le numéro de ligne d'origine (1-based) de la première ligne.
    Retourne [(lineno_start, text), ...]
    """
    out = []
    buffer = ""
    start_no = None
    for idx, raw in enumerate(lines, 1):
        line = raw.rstrip("\n")
        # On ne retire pas encore les commentaires : on les gérera par la suite
        if buffer == "":
            start_no = idx

        # fin avec &
        if line.rstrip().endswith("&"):
            # retire le & final
            buffer += line.rstrip()[:-1] + " "
            continue

        # début avec &
        if line.lstrip().startswith("&"):
            buffer += line.lstrip()[1:].lstrip()
            continue

        # ligne normale sans &
        buffer += line
        out.append((start_no, buffer))
        buffer = ""
        start_no = None

    if buffer:
        out.append((start_no, buffer))
    return out

def logical_lines_fixed_form(lines):
    """
    Recolle les continuations fixed-form :
    - Colonne 1-5 : label/espaces
    - Colonne 6 : caractère non blanc => continuation
    On retourne [(lineno_start, text), ...]
    """
    out = []
    buffer = ""
    start_no = None
    for idx, raw in enumerate(lines, 1):
        line = raw.rstrip("\n")
        # Assure une longueur minimale pour indexer
        padded = line + " " * max(0, 6 - len(line))
        cont = len(padded) >= 6 and padded[5].strip() != ""

        # corps utile (colonnes 7+), mais on garde prudent si la ligne est trop courte
        body = line[6:] if len(line) > 6 else ""

        if not cont:
            # Nouvelle carte
            if buffer:
                out.append((start_no, buffer))
                buffer = ""
            start_no = idx
            # Corps de la ligne "maîtresse" (colonnes 7+)
            buffer = body
        else:
            # continuation
            buffer += " " + body.strip()

    if buffer:
        out.append((start_no, buffer))
    return out

def iter_logical_lines(path: Path):
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        raw = f.readlines()
    if is_fixed_form(path):
        return logical_lines_fixed_form(raw)
    else:
        return logical_lines_free_form(raw)

def normalize_spaces(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip()

# --- extraction principale ---

def extract_headers_from_file(path: Path):
    """
    Retourne une liste de dicts avec:
      - 'file': chemin
      - 'lineno': numéro de la 1re ligne de l'en-tête
      - 'header': en-tête (sur une seule ligne, reconstituée)
      - 'kind': 'subroutine' ou 'function'
    """
    results = []
    for lineno, logical in iter_logical_lines(path):
        # Retire les commentaires hors guillemets
        no_comment = strip_inline_comment_outside_quotes(logical)

        # Si on repère le début d'une routine (même si multi-parens à suivre),
        # on devra peut-être continuer à lire d’autres logical-lines.
        if not START_RE.match(no_comment):
            continue

        # Accumule jusqu’à solder les parenthèses
        acc = no_comment
        bal = paren_balance(no_comment)

        # Si le header contient des parenthèses ouvertes (args, result, bind, etc.)
        # on va chercher les logical-lines suivantes jusqu’à fermer.
        # Ici, on n’a pas d’itérateur « lookahead », donc on refait un second passage.
        # Solution simple : on stocke toutes les logical-lines d'abord? Pour rester simple,
        # on relance un mini-parse à partir de la ligne actuelle.
        # Pour éviter de complexifier, on suppose que tout le header est sur CETTE logical-line
        # OU que les continuations ont déjà été recollées (ce qui est le cas).
        # => Donc bal doit déjà être 0 ici la plupart du temps.
        # Mais il peut rester >0 si le code a des lignes sans & en fixed-form ; déjà géré.

        # Normalisation légère
        header = normalize_spaces(acc)
        kind = START_RE.match(no_comment).group(1).lower()
        results.append({
            "file": str(path),
            "lineno": lineno,
            "header": header,
            "kind": kind
        })
    return results

def walk_and_extract(root: Path):
    all_results = []
    for dirpath, _, filenames in os.walk(root):
        for fn in filenames:
            p = Path(dirpath) / fn
            if p.suffix.lower() in EXTS:
                try:
                    all_results.extend(extract_headers_from_file(p))
                except Exception as e:
                    # On n’arrête pas tout si un fichier pose problème
                    all_results.append({
                        "file": str(p),
                        "lineno": 0,
                        "header": f"# ERROR while parsing: {e}",
                        "kind": "error"
                    })
    return all_results

def write_markdown(results, out_path: Path):
    results = sorted(results, key=lambda r: (r["file"].lower(), r["lineno"]))
    with out_path.open("w", encoding="utf-8") as f:
        f.write("# Fortran routine headers\n\n")
        for r in results:
            f.write(f"- **{r['kind']}** — `{r['header']}`  \n")
            f.write(f"  *File:* `{r['file']}`  —  *Line:* {r['lineno']}\n\n")
    return out_path

if __name__ == "__main__":
    results = walk_and_extract(ROOT)
    out = write_markdown(results, OUT_MD)
    print(f"✅ {len(results)} routines trouvées.")
    print(f"📝 Document écrit: {out}")
