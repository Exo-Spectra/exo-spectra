"""Plain-English pass for project documents using a local LLM (Ollama).

Usage:
    python tools/plain_english.py INPUT.md [-o OUTPUT.md] [--model gemma4:26b]
    python tools/plain_english.py --check DRAFT.md SIMPLIFIED.md

The rewrite asks the model to simplify language (CEFR B2: short sentences,
common words) while keeping every number, unit, object name and technical
term unchanged.

Hard gate: before accepting the model output, every number token
(integers, decimals, scientific notation, percentages) found in the input
must appear in the output and vice versa. Any mismatch -> nonzero exit code
and the output is NOT written. LLMs paraphrase well but occasionally drop
or mangle a figure; a publication pipeline must fail loudly when that happens.
"""
import argparse
import json
import re
import sys
import urllib.request
from collections import Counter
from pathlib import Path

OLLAMA_URL = "http://localhost:11434/api/generate"

PROMPT = """Rewrite the following Markdown document in plain, simple English \
(CEFR B2 level): short sentences, common words, active voice, no hype and no \
filler. STRICT rules:
- Keep ALL numbers, units, wavelengths, statistical values and percentages \
exactly as they are. Do not round, drop, or add any number.
- Keep all object names (planets, instruments, telescopes), author names, \
technical terms, links, and Markdown structure (headings, tables, lists) unchanged.
- Do not add new claims or explanations that are not in the text.
Return only the rewritten Markdown, no commentary.

Document:

{text}"""

# number tokens: 1.4e-5 | 12.5% | 3,991 | 0.656 | 42
NUM_RE = re.compile(r"\d+(?:[.,]\d+)*(?:[eE][+-]?\d+)?%?")


def numbers(text: str) -> Counter:
    """Multiset of number tokens, normalized (thousands separators stripped)."""
    toks = []
    for t in NUM_RE.findall(text):
        t = t.rstrip("%")
        # "3,991" -> "3991" but keep "1.4" as-is
        if "," in t and "." not in t:
            t = t.replace(",", "")
        toks.append(t)
    return Counter(toks)


def check(src: str, out: str) -> list[str]:
    a, b = numbers(src), numbers(out)
    problems = []
    for tok, n in (a - b).items():
        problems.append(f"missing in output: {tok!r} (x{n})")
    for tok, n in (b - a).items():
        problems.append(f"invented in output: {tok!r} (x{n})")
    return problems


def rewrite(text: str, model: str) -> str:
    req = urllib.request.Request(
        OLLAMA_URL,
        # think=False: gemma4 is a reasoning model; without this the whole
        # generation goes into hidden thinking until the context runs out
        data=json.dumps({"model": model, "prompt": PROMPT.format(text=text),
                         "stream": False, "think": False,
                         "options": {"temperature": 0.3}}).encode(),
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=1200) as resp:
        return json.loads(resp.read())["response"].strip()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("output_for_check", nargs="?",
                    help="with --check: the simplified file to verify")
    ap.add_argument("-o", "--output", help="output path (default: INPUT.plain.md)")
    ap.add_argument("--model", default="gemma4:26b")
    ap.add_argument("--check", action="store_true",
                    help="only verify numbers between two existing files")
    ap.add_argument("--retries", type=int, default=2)
    args = ap.parse_args()

    src = Path(args.input).read_text(encoding="utf-8")

    if args.check:
        out = Path(args.output_for_check).read_text(encoding="utf-8")
        problems = check(src, out)
        for p in problems:
            print("GATE FAIL:", p)
        print("numbers gate:", "FAIL" if problems else "OK")
        return 1 if problems else 0

    for attempt in range(1, args.retries + 2):
        out = rewrite(src, args.model)
        problems = check(src, out)
        if not problems:
            dest = Path(args.output or Path(args.input).with_suffix(".plain.md"))
            dest.write_text(out + "\n", encoding="utf-8")
            print(f"OK -> {dest}")
            return 0
        print(f"attempt {attempt}: numbers gate failed ({len(problems)} problems)")
        for p in problems[:10]:
            print("  ", p)
    print("FAIL: model output kept violating the numbers gate; nothing written")
    return 1


if __name__ == "__main__":
    sys.exit(main())
