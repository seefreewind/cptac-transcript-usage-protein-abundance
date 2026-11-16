#!/usr/bin/env python3
"""Acquire the LinkedOmicsKB v1.2 matrices with resumable multi-connection downloads.

Same S3 URLs as Phase 0. A file is renamed out of .partial only after the byte
count matches Content-Length, the header parses, and a whole-file scan finds
no NUL bytes. A full-size partial that still contains NULs is deleted rather
than resumed. Phase 0 files are not touched.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HEADS = json.loads((ROOT / "phase0/raw/api/s3_head_selected.json").read_text())
ACQ = ROOT / "phase0_5" / "acquisition"
MAN = ROOT / "phase0_5" / "manifests"
LOG = ROOT / "phase0_5" / "logs"
for p in (ACQ, MAN, LOG):
    p.mkdir(parents=True, exist_ok=True)

CANCERS = ["BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC"]
SUFFIXES = (
    "RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
    "RNAseq_isoform_FPKM_log2_Tumor.txt",
    "proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
    "phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt",
    "WES_CNV_gene_ratio_log2.txt",
    "_meta.txt",
    "_phenotype.txt",
)


def utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def layer_of(name: str) -> str:
    if name.endswith("_meta.txt"):
        return "clinical_meta"
    if name.endswith("_phenotype.txt"):
        return "phenotype"
    if "RNAseq_gene_" in name:
        return "rna_gene_tumor"
    if "RNAseq_isoform_" in name:
        return "rna_isoform_tumor"
    if "proteomics_gene_" in name:
        return "protein_tumor"
    if "phospho_site_" in name:
        return "phospho_tumor"
    if "CNV_gene_" in name:
        return "cnv"
    return "other"


def items() -> list[dict]:
    out = []
    for h in HEADS:
        cancer = h["url"].rstrip("/").split("/")[-2]
        if cancer not in CANCERS:
            continue
        if any(h["name"].endswith(s) for s in SUFFIXES):
            out.append(h)
    out.sort(key=lambda h: int(h["bytes"]))
    return out


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def valid_format(path: Path, layer: str) -> tuple[bool, str]:
    with path.open("rb") as f:
        blob = f.read(65536)
    if not blob:
        return False, "empty"
    low = blob[:300].lower()
    if b"<html" in low or b"accessdenied" in low or low.startswith(b"<?xml"):
        return False, "not_a_matrix"
    if b"\x00" in blob:
        return False, "nul_byte"
    line = blob.split(b"\n", 1)[0]
    if b"\t" not in line:
        return False, "header_not_tab_separated"
    first = line.split(b"\t", 1)[0]
    if layer in {"rna_gene_tumor", "rna_isoform_tumor", "protein_tumor", "phospho_tumor", "cnv", "phenotype"}:
        if first != b"idx":
            return False, f"header_token_{first[:40]!r}"
    if layer == "clinical_meta" and b"case" not in line.lower():
        return False, "meta_header_unrecognized"
    return True, "ok"


def head(url: str) -> dict:
    p = subprocess.run(
        ["curl", "-sI", "-L", "--http1.1", "--max-time", "60", "-A", "phase0.5-acquisition", "-w", "\n__HTTP__%{http_code}", url],
        capture_output=True,
    )
    text = p.stdout.decode("latin1", "replace")
    status = ""
    if "__HTTP__" in text:
        text, status = text.rsplit("__HTTP__", 1)
    headers = {}
    for line in text.splitlines():
        if ":" in line:
            k, v = line.split(":", 1)
            headers[k.strip().lower()] = v.strip()
    return {
        "http_status": status.strip(),
        "content_length": headers.get("content-length", ""),
        "etag": headers.get("etag", ""),
        "last_modified": headers.get("last-modified", ""),
        "curl_rc": p.returncode,
    }


def log_event(rec: dict) -> None:
    rec = {"timestamp_utc": utc(), **rec}
    with (LOG / "download_events.jsonl").open("a") as f:
        f.write(json.dumps(rec, ensure_ascii=False) + "\n")


def promote(partial: Path, final: Path, layer: str, expected: int) -> tuple[bool, str, str]:
    if partial.stat().st_size != expected:
        return False, "size_mismatch", ""
    ok, why = valid_format(partial, layer)
    if not ok:
        return False, why, ""
    with partial.open("rb") as handle:
        while True:
            blob = handle.read(1 << 20)
            if not blob:
                break
            if b"\x00" in blob:
                return False, "embedded_nul", ""
    digest = sha256(partial)
    os.replace(partial, final)
    control = Path(str(partial) + ".aria2")
    if control.exists():
        control.unlink()
    return True, "ok", digest


def discard_partial(partial: Path) -> None:
    if partial.exists():
        partial.unlink()
    control = Path(str(partial) + ".aria2")
    if control.exists():
        control.unlink()


def download_one(item: dict) -> dict:
    url = item["url"]
    cancer = url.rstrip("/").split("/")[-2]
    name = item["name"]
    layer = layer_of(name)
    folder = ACQ / cancer
    folder.mkdir(parents=True, exist_ok=True)
    final = folder / name
    partial = folder / (name + ".partial")
    rec = {"cancer": cancer, "layer": layer, "file_name": name, "source_url": url}
    if final.exists() and final.stat().st_size == int(item["bytes"]):
        digest = sha256(final)
        rec.update({"ok": True, "action": "already_final", "bytes": final.stat().st_size, "sha256": digest, "http_status": "200", "content_length": item["bytes"], "etag": item.get("etag", ""), "last_modified": item.get("last_modified", "")})
        log_event(rec)
        print(f"SKIP {cancer} {name}", flush=True)
        return rec
    info = head(url)
    rec.update(info)
    if info["http_status"] != "200" or not str(info["content_length"]).isdigit():
        rec.update({"ok": False, "failure": "head_failed"})
        log_event(rec)
        print(f"FAIL head {name}", flush=True)
        return rec
    expected = int(info["content_length"])
    # A full-size partial can still be a sparse file. Promote only when the
    # NUL scan passes. Otherwise delete it; resuming that control file skips
    # the holes.
    if partial.exists() and partial.stat().st_size == expected:
        ok, why, digest = promote(partial, final, layer, expected)
        if ok:
            rec.update({"ok": True, "failure": "", "sha256": digest, "bytes": expected, "action": "promote_partial"})
            log_event(rec)
            print(f"OK promote {name}", flush=True)
            return rec
        discard_partial(partial)
        rec_fail = dict(rec)
        rec_fail.update({"ok": False, "failure": why, "bytes": expected, "action": "discard_sparse_partial"})
        log_event(rec_fail)
        print(f"DISCARD {name} {why}", flush=True)
    cmd = [
        "aria2c", "-x", "16", "-s", "16", "-k", "1M",
        "--continue=true", "--auto-file-renaming=false", "--allow-overwrite=true",
        "--file-allocation=none", "--max-tries=30", "--retry-wait=3",
        "--timeout=60", "--connect-timeout=20",
        "--console-log-level=warn", "--summary-interval=60",
        "-d", str(folder), "-o", name + ".partial", url,
    ]
    last = rec
    for attempt in (1, 2):
        print(f"GET {cancer} {name} {expected} attempt {attempt}", flush=True)
        with (LOG / "aria2c.log").open("ab") as aria_log:
            aria_log.write(f"\n# {utc()} {name} attempt {attempt}\n".encode())
            proc = subprocess.run(cmd, stdout=aria_log, stderr=subprocess.STDOUT)
        got = partial.stat().st_size if partial.exists() else 0
        if got == expected:
            ok, why, digest = promote(partial, final, layer, expected)
            if ok:
                rec.update({"ok": True, "failure": "", "sha256": digest, "bytes": expected, "action": "download", "aria2_rc": proc.returncode, "attempt": attempt})
                log_event(rec)
                print(f"OK {name}", flush=True)
                return rec
            discard_partial(partial)
            last = dict(rec)
            last.update({"ok": False, "failure": why, "bytes": expected, "action": "download", "aria2_rc": proc.returncode, "attempt": attempt})
            log_event(last)
            print(f"FAIL {name} {why} attempt {attempt}", flush=True)
            continue
        last = dict(rec)
        last.update({"ok": False, "failure": "incomplete_or_aria2_error", "aria2_rc": proc.returncode, "bytes": got, "expected_bytes": expected, "attempt": attempt})
        log_event(last)
        print(f"FAIL {name} rc={proc.returncode} bytes={got}/{expected} attempt {attempt}", flush=True)
        if got == 0:
            break
    return last


def write_tables() -> None:
    events = []
    path = LOG / "download_events.jsonl"
    if path.exists():
        for line in path.read_text().splitlines():
            if line.strip():
                events.append(json.loads(line))
    latest = {}
    for ev in events:
        if ev.get("file_name") and ev.get("source_url", "").startswith("https://cptac-pancancer-data"):
            latest[ev["file_name"]] = ev
    rows = []
    for item in items():
        cancer = item["url"].rstrip("/").split("/")[-2]
        name = item["name"]
        final = ACQ / cancer / name
        ev = latest.get(name, {})
        exists = final.exists()
        digest = ev.get("sha256") or (sha256(final) if exists else "")
        rows.append({
            "cancer": cancer,
            "layer": layer_of(name),
            "file_name": name,
            "source_url": item["url"],
            "source_freeze": "LinkedOmicsKB Release 1 / data_freeze_v1.2_reorganized",
            "retrieval_timestamp_utc": ev.get("timestamp_utc", ""),
            "http_status": ev.get("http_status", ""),
            "content_length": ev.get("content_length") or ev.get("expected_bytes") or item["bytes"],
            "etag": ev.get("etag") or item.get("etag", ""),
            "last_modified": ev.get("last_modified") or item.get("last_modified", ""),
            "local_path": str(final.relative_to(ROOT)) if exists else "",
            "bytes": final.stat().st_size if exists else 0,
            "sha256": digest,
            "validated": "yes" if exists and digest else "no",
            "failure": ev.get("failure", ""),
        })
    fields = list(rows[0].keys())
    with (MAN / "acquisition_manifest.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    log_fields = ["timestamp_utc", "action", "cancer", "layer", "file_name", "source_url", "http_status", "content_length", "etag", "last_modified", "bytes", "sha256", "ok", "failure", "aria2_rc"]
    with (MAN / "download_log.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=log_fields, extrasaction="ignore")
        w.writeheader()
        for ev in events:
            if str(ev.get("source_url", "")).startswith("https://cptac-pancancer-data"):
                w.writerow(ev)
    with (MAN / "download_failures.csv").open("w", newline="") as f:
        cols = ["cancer", "layer", "file_name", "source_url", "http_status", "content_length", "bytes", "failure"]
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            if r["validated"] != "yes":
                w.writerow({k: r.get(k, "") for k in cols})
    with (MAN / "file_integrity_report.csv").open("w", newline="") as f:
        cols = ["cancer", "layer", "file_name", "bytes", "content_length", "size_matches_content_length", "sha256", "format_check", "validated"]
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            fmt = ""
            match = ""
            if r["local_path"]:
                ok, why = valid_format(ROOT / r["local_path"], r["layer"])
                fmt = why
            if r["content_length"] != "" and r["bytes"]:
                match = "yes" if str(r["bytes"]) == str(r["content_length"]) else "no"
            w.writerow({
                "cancer": r["cancer"], "layer": r["layer"], "file_name": r["file_name"],
                "bytes": r["bytes"], "content_length": r["content_length"],
                "size_matches_content_length": match, "sha256": r["sha256"],
                "format_check": fmt, "validated": r["validated"],
            })
    lines = [f"{r['sha256']}  {r['local_path']}" for r in rows if r["sha256"] and r["local_path"]]
    text = "\n".join(lines) + ("\n" if lines else "")
    # preserve a SEPEP section if already appended
    existing = ROOT / "phase0_5" / "SHA256SUMS"
    sepep = ""
    if existing.exists():
        body = existing.read_text()
        if "\n# SEPEP" in body:
            sepep = body[body.index("\n# SEPEP"):]
            if not sepep.endswith("\n"):
                sepep += "\n"
    existing.write_text(text + sepep.lstrip("\n") if False else text)
    if sepep:
        existing.write_text(text + sepep[1:] if sepep.startswith("\n") else text + sepep)
    n = sum(r["validated"] == "yes" for r in rows)
    print(f"validated {n} / {len(rows)}", flush=True)


def main() -> None:
    queue = items()
    print(f"queue {len(queue)}", flush=True)
    for item in queue:
        download_one(item)
        write_tables()
    # second pass for anything still missing
    missing = []
    for item in queue:
        cancer = item["url"].rstrip("/").split("/")[-2]
        if not (ACQ / cancer / item["name"]).exists():
            missing.append(item)
    if missing:
        print(f"second pass {len(missing)}", flush=True)
        for item in missing:
            download_one(item)
            write_tables()


if __name__ == "__main__":
    main()
