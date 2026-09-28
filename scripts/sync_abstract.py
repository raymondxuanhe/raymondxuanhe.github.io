#!/usr/bin/env python3
"""Sync paper abstracts from their LaTeX sources of truth into the website.

Each paper keeps its abstract in an abstract.tex that is \\input directly into the
paper and the CV (so those update on recompile). The website is HTML and cannot
read .tex, so this script copies the text into the matching marked block in
_pages/research.md.

To add another paper: create its abstract.tex, \\input it in the paper + CV, wrap
the website block in START/END markers with a new id, and add a row to ABSTRACTS.

Usage:
    python3 scripts/sync_abstract.py          # write any changes
    python3 scripts/sync_abstract.py --check  # exit 1 if anything is out of sync (CI)
"""
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]        # raymondxuanhe.github.io
DOCS = REPO.parent                                # ~/Documents (all repos side by side)
RESEARCH_MD = REPO / "_pages" / "research.md"

# marker id  ->  abstract.tex source of truth
ABSTRACTS = {
    "currency-maturity": DOCS / "Currency-Maturity-Composition" / "7-Writing" / "Main" / "abstract.tex",
    "pay-transparency": DOCS / "Pay-Transparency" / "Writing" / "Chapters" / "abstract.tex",
}


def latex_to_text(tex: str) -> str:
    """Turn the plain-prose abstract body into one clean HTML-safe paragraph."""
    lines = [ln for ln in tex.splitlines() if not ln.lstrip().startswith("%")]
    text = " ".join(lines)
    # minimal LaTeX -> text/HTML fixups (extend if an abstract gains new markup)
    text = text.replace("\\%", "%").replace("\\_", "_").replace("\\$", "$")
    text = text.replace("\\&", "&amp;").replace("~", " ")
    text = re.sub(r"\s+", " ", text).strip()
    return text


def build_block(marker_id: str, indent: str, body: str) -> str:
    return (
        f"{indent}<!-- ABSTRACT:{marker_id} START (auto-synced from abstract.tex "
        f"by scripts/sync_abstract.py — do not edit by hand) -->\n"
        f"{indent}{body}\n"
        f"{indent}<!-- ABSTRACT:{marker_id} END -->"
    )


def main() -> int:
    check_only = "--check" in sys.argv

    md = RESEARCH_MD.read_text(encoding="utf-8")
    problems: list[str] = []
    changed = False

    for marker_id, src in ABSTRACTS.items():
        if not src.exists():
            problems.append(f"source not found for '{marker_id}': {src}")
            continue

        body = latex_to_text(src.read_text(encoding="utf-8"))
        start = re.escape(f"<!-- ABSTRACT:{marker_id} START")
        end = re.escape(f"<!-- ABSTRACT:{marker_id} END -->")
        m = re.search(rf"(?m)^([ \t]*){start}.*?{end}", md, re.DOTALL)
        if not m:
            problems.append(f"markers for '{marker_id}' not found in {RESEARCH_MD.name}")
            continue

        new_block = build_block(marker_id, m.group(1), body)
        if md[m.start():m.end()] == new_block:
            print(f"  {marker_id}: already in sync")
            continue

        if check_only:
            problems.append(f"'{marker_id}' is OUT OF SYNC")
            continue

        md = md[: m.start()] + new_block + md[m.end():]
        changed = True
        print(f"  {marker_id}: updated")

    if problems:
        for p in problems:
            print(f"ERROR: {p}", file=sys.stderr)
        return 1 if check_only else 2

    if changed:
        RESEARCH_MD.write_text(md, encoding="utf-8")
        print(f"Wrote {RESEARCH_MD.relative_to(REPO)}.")
    else:
        print("Nothing to write.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
