#!/usr/bin/env python3
"""审计论文引用、占位符和关键物理表述。"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANUSCRIPT = ROOT / "paper" / "manuscript_v1.md"
DEFAULT_TEX = ROOT / "paper" / "manuscript_v1.tex"
DEFAULT_BIBLIOGRAPHY = ROOT / "results" / "literature_search_v1" / "references.bib"

CITATION_RE = re.compile(r"(?<![\w.])@([A-Za-z0-9_:.+\-]+)")
ENTRY_RE = re.compile(r"(?mi)^\s*@([A-Za-z]+)\s*\{\s*([^,\s]+)\s*,")
DOI_RE = re.compile(r'(?mi)^\s*doi\s*=\s*(?:\{([^}]*)\}|"([^"]*)")\s*,?')
INTERNAL_SOURCE_RE = re.compile(
    r"\b(?:source|source asset|existing panel)\s*:\s*`?(?:\.\.?[/\\])?"
    r"(?:figures|results|scripts|paper|docs)[/\\]",
    re.IGNORECASE,
)
PRODUCTION_TODO_RE = re.compile(
    r"TODO.*(?:production|thermodynamic|finite-size|pooled observables|final SSE)"
    r"|(?:production|thermodynamic|finite-size|pooled observables|final SSE).*TODO",
    re.IGNORECASE,
)
STATIC_METHOD_RE = re.compile(r"\b(?:static\s+)?(?:SSE|QMC|quantum Monte Carlo)\b", re.IGNORECASE)
DYNAMIC_TARGET_RE = re.compile(
    r"\b(?:magnon\s+(?:flat\s+)?band|flat\s+(?:magnon|magnetic|spin-wave)\s+(?:band|branch)"
    r"|magnon\s+dispersion|excitation\s+spectrum)\b",
    re.IGNORECASE,
)
PROOF_VERB_RE = re.compile(
    r"\b(?:prove[sd]?|establish(?:es|ed)?|verif(?:y|ies|ied)|confirm(?:s|ed)?"
    r"|demonstrat(?:e|es|ed)|determin(?:e|es|ed)|resolv(?:e|es|ed))\b",
    re.IGNORECASE,
)
NEGATION_RE = re.compile(r"\b(?:cannot|does\s+not|do\s+not|did\s+not|neither|not)\b", re.IGNORECASE)
LIEB_FLAT_RE = re.compile(
    r"Lieb\s*--?\s*Mattis.*\b(?:protect(?:s|ed)?|prove[sd]?|establish(?:es|ed)?)\b"
    r".*\b(?:flat|dispersion|excitation)\b",
    re.IGNORECASE,
)
RHO_S_RE = re.compile(r"(?:\\rho_(?:s\b|\{s\})|\brho_s\b)")
BETA_NUMERIC_RE = re.compile(
    r"\\beta(?!\s*J)\s*(?:=\s*(?:[0-9L]|\\)|/\s*(?!\\nu\b)(?:[0-9L]|\\))"
)
AMBIGUOUS_ENERGY_RE = re.compile(r"(?:\$|\\\()[^\n$]*\bE/NJ\b")


def line_number(text: str, offset: int) -> int:
    return text.count("\n", 0, offset) + 1


def excerpt(line: str, limit: int = 180) -> str:
    clean = " ".join(line.strip().split())
    return clean if len(clean) <= limit else clean[: limit - 1] + "…"


def normalize_doi(value: str) -> str:
    doi = value.strip().lower()
    doi = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", doi)
    return doi.rstrip(".,; ")


def parse_bibliography(text: str) -> tuple[list[dict[str, Any]], dict[str, list[dict[str, Any]]]]:
    matches = list(ENTRY_RE.finditer(text))
    entries: list[dict[str, Any]] = []
    doi_entries: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for index, match in enumerate(matches):
        kind = match.group(1).lower()
        if kind in {"comment", "preamble", "string"}:
            continue
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        block = text[match.start() : end]
        doi_match = DOI_RE.search(block)
        doi = normalize_doi((doi_match.group(1) or doi_match.group(2))) if doi_match else None
        entry = {
            "key": match.group(2),
            "key_folded": match.group(2).casefold(),
            "type": kind,
            "line": line_number(text, match.start()),
            "doi": doi,
        }
        entries.append(entry)
        if doi:
            doi_entries[doi].append(entry)
    return entries, doi_entries


def finding(
    code: str,
    severity: str,
    source: str,
    line: int,
    message: str,
    text: str = "",
) -> dict[str, Any]:
    item: dict[str, Any] = {
        "code": code,
        "severity": severity,
        "source": source,
        "line": line,
        "message": message,
    }
    if text:
        item["excerpt"] = excerpt(text)
    return item


def audit_text(manuscript: str, bibliography: str) -> dict[str, Any]:
    findings: list[dict[str, Any]] = []
    entries, doi_entries = parse_bibliography(bibliography)

    by_key: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for entry in entries:
        by_key[entry["key_folded"]].append(entry)
    for folded, copies in sorted(by_key.items()):
        if len(copies) > 1:
            findings.append(
                finding(
                    "DUPLICATE_BIBTEX_KEY",
                    "error",
                    "bibliography",
                    copies[1]["line"],
                    f"BibTeX key 重复：{', '.join(item['key'] for item in copies)}",
                )
            )
    for doi, copies in sorted(doi_entries.items()):
        if len(copies) > 1:
            findings.append(
                finding(
                    "DUPLICATE_DOI",
                    "error",
                    "bibliography",
                    copies[1]["line"],
                    f"DOI {doi} 被多个条目使用：{', '.join(item['key'] for item in copies)}",
                )
            )

    citation_occurrences: dict[str, list[int]] = defaultdict(list)
    for match in CITATION_RE.finditer(manuscript):
        citation_occurrences[match.group(1)].append(line_number(manuscript, match.start()))
    defined_exact = {entry["key"] for entry in entries}
    for key in sorted(citation_occurrences):
        if key not in defined_exact:
            findings.append(
                finding(
                    "MISSING_CITATION_KEY",
                    "error",
                    "manuscript",
                    citation_occurrences[key][0],
                    f"引用 key 未在 bibliography 中定义：{key}",
                )
            )

    for number, line in enumerate(manuscript.splitlines(), start=1):
        if "TODO" in line and PRODUCTION_TODO_RE.search(line):
            findings.append(
                finding(
                    "PRODUCTION_RESULT_TODO",
                    "error",
                    "manuscript",
                    number,
                    "正式生产结果仍有 TODO",
                    line,
                )
            )
        if INTERNAL_SOURCE_RE.search(line):
            findings.append(
                finding(
                    "INTERNAL_SOURCE_PATH",
                    "error",
                    "manuscript",
                    number,
                    "图注或正文仍含内部 Source 路径",
                    line,
                )
            )
        if RHO_S_RE.search(line):
            findings.append(
                finding(
                    "RHO_SYMBOL_MIXED",
                    "warning",
                    "manuscript",
                    number,
                    "本稿约定使用 \\rho_{\\rm cell}，此处仍使用 \\rho_s",
                    line,
                )
            )
        if BETA_NUMERIC_RE.search(line):
            findings.append(
                finding(
                    "BETA_J_MISSING",
                    "warning",
                    "manuscript",
                    number,
                    "数值化逆温或固定比例应写成 \\beta J",
                    line,
                )
            )
        if AMBIGUOUS_ENERGY_RE.search(line):
            findings.append(
                finding(
                    "ENERGY_NORMALIZATION_AMBIGUOUS",
                    "warning",
                    "manuscript",
                    number,
                    "能量归一化建议写成 E/(NJ)",
                    line,
                )
            )

        for sentence in re.split(r"(?<=[.!?])\s+", line):
            if NEGATION_RE.search(sentence):
                continue
            static_claim = (
                STATIC_METHOD_RE.search(sentence)
                and DYNAMIC_TARGET_RE.search(sentence)
                and PROOF_VERB_RE.search(sentence)
            )
            if static_claim or LIEB_FLAT_RE.search(sentence):
                findings.append(
                    finding(
                        "BANNED_DYNAMIC_CLAIM",
                        "error",
                        "manuscript",
                        number,
                        "静态 SSE/QMC 或 Lieb--Mattis 定理不能证明 magnon 平带或色散",
                        sentence,
                    )
                )

    severity_counts = {
        level: sum(item["severity"] == level for item in findings)
        for level in ("error", "warning")
    }
    used_keys = sorted(citation_occurrences)
    duplicate_dois = {
        doi: [entry["key"] for entry in copies]
        for doi, copies in sorted(doi_entries.items())
        if len(copies) > 1
    }
    duplicate_keys = [
        [entry["key"] for entry in copies]
        for _, copies in sorted(by_key.items())
        if len(copies) > 1
    ]
    return {
        "schema_version": 1,
        "ok": not findings,
        "summary": {"findings": len(findings), **severity_counts},
        "citations": {
            "used": used_keys,
            "defined": sorted(defined_exact),
            "missing": sorted(set(used_keys) - defined_exact),
        },
        "bibliography": {
            "entry_count": len(entries),
            "duplicate_keys": duplicate_keys,
            "duplicate_dois": duplicate_dois,
        },
        "findings": findings,
    }


def audit_tex_text(tex: str) -> list[dict[str, Any]]:
    """检查投稿源文件中的关键记号。"""
    findings: list[dict[str, Any]] = []
    checks = (
        (RHO_S_RE, "RHO_SYMBOL_MIXED", "本稿约定使用 \\rho_{\\rm cell}，此处仍使用 \\rho_s"),
        (BETA_NUMERIC_RE, "BETA_J_MISSING", "数值化逆温或固定比例应写成 \\beta J"),
        (AMBIGUOUS_ENERGY_RE, "ENERGY_NORMALIZATION_AMBIGUOUS", "能量归一化建议写成 E/(NJ)"),
    )
    for number, line in enumerate(tex.splitlines(), start=1):
        for pattern, code, message in checks:
            if pattern.search(line):
                findings.append(finding(code, "warning", "manuscript_tex", number, message, line))
    return findings


def audit_files(
    manuscript_path: Path,
    bibliography_path: Path,
    tex_path: Path | None = None,
) -> dict[str, Any]:
    result = audit_text(
        manuscript_path.read_text(encoding="utf-8"),
        bibliography_path.read_text(encoding="utf-8"),
    )
    if tex_path is not None and tex_path.exists():
        result["findings"].extend(audit_tex_text(tex_path.read_text(encoding="utf-8")))
        result["summary"] = {
            "findings": len(result["findings"]),
            "error": sum(item["severity"] == "error" for item in result["findings"]),
            "warning": sum(item["severity"] == "warning" for item in result["findings"]),
        }
        result["ok"] = not result["findings"]
    result["manuscript"] = str(manuscript_path.resolve())
    result["bibliography_path"] = str(bibliography_path.resolve())
    result["tex_path"] = str(tex_path.resolve()) if tex_path is not None else None
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manuscript", type=Path, default=DEFAULT_MANUSCRIPT)
    parser.add_argument("--tex", type=Path, default=DEFAULT_TEX)
    parser.add_argument("--bibliography", type=Path, default=DEFAULT_BIBLIOGRAPHY)
    parser.add_argument("--output", type=Path, help="可选 JSON 输出路径")
    args = parser.parse_args()

    result = audit_files(args.manuscript, args.bibliography, args.tex)
    payload = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    else:
        sys.stdout.write(payload)
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
