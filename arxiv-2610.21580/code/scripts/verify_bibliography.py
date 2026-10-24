#!/usr/bin/env python3
"""通过 Crossref 核验稿件引用并生成独立 BibTeX。"""

from __future__ import annotations

import argparse
import concurrent.futures
import difflib
import html
import json
import re
import time
import unicodedata
import urllib.parse
import urllib.request
import urllib.error
from pathlib import Path


def cited_keys(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    result: list[str] = []
    for group in re.findall(r"\\cite\w*\{([^}]+)\}", text):
        for key in group.split(","):
            key = key.strip()
            if key and key not in result:
                result.append(key)
    return result


def parse_entries(path: Path) -> dict[str, dict[str, str]]:
    text = path.read_text(encoding="utf-8")
    starts = list(re.finditer(r"(?m)^@(\w+)\{([^,]+),", text))
    entries: dict[str, dict[str, str]] = {}
    for index, match in enumerate(starts):
        end = starts[index + 1].start() if index + 1 < len(starts) else len(text)
        block = text[match.start() : end].strip()
        fields = {}
        for field_match in re.finditer(r"(?m)^\s*(\w+)\s*=\s*\{(.*)\},?\s*$", block):
            fields[field_match.group(1).lower()] = field_match.group(2)
        entries[match.group(2).strip()] = {
            "type": match.group(1),
            "block": block,
            **fields,
        }
    return entries


def normalized(value: str) -> str:
    value = html.unescape(value)
    value = re.sub(r"<[^>]+>", " ", value)
    value = re.sub(r"\\[\"'`^~=.uvHckbdtr]\s*", "", value)
    value = re.sub(r"\\[A-Za-z]+", "", value)
    value = value.replace("{", "").replace("}", "")
    value = unicodedata.normalize("NFKD", value)
    value = "".join(character for character in value if not unicodedata.combining(character))
    return re.sub(r"[^a-z0-9]+", " ", value.lower()).strip()


def crossref(doi: str) -> dict:
    encoded = urllib.parse.quote(doi, safe="/")
    request = urllib.request.Request(
        f"https://doi.org/{encoded}",
        headers={
            "Accept": "application/vnd.citationstyles.csl+json",
            "User-Agent": "QMC-manuscript-verifier/1.0",
        },
    )
    last_error: Exception | None = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(request, timeout=15) as response:
                payload = json.load(response)
            return payload.get("message", payload)
        except (urllib.error.URLError, TimeoutError, ConnectionError) as error:
            last_error = error
            time.sleep(1.0 * (attempt + 1))
    raise RuntimeError(str(last_error))


def first_year(record: dict) -> int | None:
    for field in ("published-print", "published-online", "issued", "created"):
        parts = record.get(field, {}).get("date-parts", [])
        if parts and parts[0]:
            return int(parts[0][0])
    return None


def verify_one(key: str, entry: dict[str, str] | None) -> dict:
    if entry is None:
        return {"key": key, "status": "missing_entry"}
    doi = entry.get("doi", "").strip()
    if not doi:
        return {"key": key, "status": "manual_needed", "reason": "missing DOI"}
    try:
        found = crossref(doi)
    except Exception as error:
        return {"key": key, "doi": doi, "status": "not_found", "reason": str(error)}
    crossref_title = (found.get("title") or [""])
    if isinstance(crossref_title, list):
        crossref_title = crossref_title[0] if crossref_title else ""
    crossref_journal = found.get("container-title") or ""
    if isinstance(crossref_journal, list):
        crossref_journal = crossref_journal[0] if crossref_journal else ""
    crossref_year = first_year(found)
    title_score = difflib.SequenceMatcher(
        None, normalized(entry.get("title", "")), normalized(crossref_title)
    ).ratio()
    manuscript_year = int(entry["year"]) if entry.get("year", "").isdigit() else None
    first_author = normalized(entry.get("author", "").split(" and ")[0].split(",")[0])
    crossref_first = normalized((found.get("author") or [{}])[0].get("family", ""))
    checks = {
        "doi_resolves": True,
        "title_match": title_score >= 0.9,
        "year_match": manuscript_year == crossref_year,
        "first_author_match": bool(first_author and first_author == crossref_first),
    }
    status = "verified" if all(checks.values()) else "mismatch"
    return {
        "key": key,
        "doi": doi,
        "status": status,
        "checks": checks,
        "title_similarity": title_score,
        "manuscript": {
            "title": entry.get("title"),
            "journal": entry.get("journal"),
            "year": manuscript_year,
            "first_author": first_author,
        },
        "crossref": {
            "title": crossref_title,
            "journal": crossref_journal,
            "year": crossref_year,
            "first_author": crossref_first,
            "type": found.get("type"),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manuscript", type=Path)
    parser.add_argument("bibliography", type=Path)
    parser.add_argument("output_bib", type=Path)
    parser.add_argument("report", type=Path)
    args = parser.parse_args()

    keys = cited_keys(args.manuscript)
    entries = parse_entries(args.bibliography)
    cache = {}
    report_json = args.report.with_suffix(".json")
    if report_json.exists():
        previous = json.loads(report_json.read_text(encoding="utf-8"))
        cache = {
            record["key"]: record
            for record in previous.get("records", [])
            if record.get("status") == "verified"
        }
    pending = [key for key in keys if key not in cache]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        checked = dict(
            zip(pending, executor.map(lambda key: verify_one(key, entries.get(key)), pending))
        )
    records = [cache.get(key, checked.get(key)) for key in keys]

    args.output_bib.parent.mkdir(parents=True, exist_ok=True)
    verified_blocks = [entries[key]["block"] for key in keys if key in entries]
    args.output_bib.write_text(
        "% Crossref-verified references cited by manuscript_prb.tex\n\n"
        + "\n\n".join(verified_blocks)
        + "\n",
        encoding="utf-8",
    )
    counts = {
        status: sum(record["status"] == status for record in records)
        for status in ("verified", "mismatch", "not_found", "manual_needed", "missing_entry")
    }
    payload = {"format": "citation-verification-v1", "counts": counts, "records": records}
    args.report.with_suffix(".json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    lines = ["# Citation verification", "", f"Total cited entries: {len(records)}", ""]
    lines.append("| Key | Status | DOI | Notes |")
    lines.append("|---|---|---|---|")
    for record in records:
        notes = ""
        if record["status"] == "mismatch":
            failed = [name for name, passed in record["checks"].items() if not passed]
            notes = ", ".join(failed)
        elif record["status"] not in ("verified",):
            notes = record.get("reason", "")
        lines.append(
            f"| {record['key']} | {record['status']} | {record.get('doi', '')} | {notes} |"
        )
    args.report.with_suffix(".md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(counts, ensure_ascii=False))


if __name__ == "__main__":
    main()
