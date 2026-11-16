#!/usr/bin/env python3
"""Resumable acquisition of the LinkedOmicsKB v1.2 open matrices.

Downloads to *.partial and renames only after Content-Length and format checks.
Does not touch Phase 0 files.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HEADS = json.loads((ROOT / "phase0/raw/api/s3_head_selected.json").read_text())
ACQ = ROOT / "phase0_5" / "acquisition"
LOG = ROOT / "phase0_5" / "logs"
MAN = ROOT / "phase0_5" / "manifests"
for p in (ACQ, LOG, MAN):
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
EVENT = LOG / "download_events.jsonl"
LOCK = threading.Lock()


def utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def log(event: dict) -> None:
    event["timestamp_utc"] = utc()
    with LOCK:
        with EVENT.open("a") as f:
            f.write(json.dumps(event, ensure_ascii=False) + "\n")


def select() -> list[dict]:
    out = []
    for h in HEADS:
        name = h["name"]
        cancer = h["url"].rstrip("/").split("/")[-2]
        if cancer not in CANCERS:
            continue
        if not any(name.endswith(s) for s in SUFFIXES):
            continue
        out.append(h)
    return out


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


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def head(url: str) -> dict:
    cmd = [
        "curl", "-sI", "-L", "--max-time", "60",
        "-A", "phase0.5-acquisition",
        "-w", "\n__HTTP__%{http_code}",
        url,
    ]
    p = subprocess.run(cmd, capture_output=True)
    text = p.stdout.decode("latin1", "replace")
    status = ""
    if "__HTTP__" in text:
        text, status = text.rsplit("__HTTP__", 1)
        status = status.strip()
    headers = {}
    for line in text.splitlines():
        if ":" in line:
            k, v = line.split(":", 1)
            headers[k.strip().lower()] = v.strip()
    return {
        "curl_rc": p.returncode,
        "http_status": status,
        "content_length": headers.get("content-length", ""),
        "etag": headers.get("etag", ""),
        "last_modified": headers.get("last-modified", ""),
        "content_type": headers.get("content-type", ""),
    }


def valid_format(path: Path, layer: str) -> tuple[bool, str]:
    with path.open("rb") as f:
        blob = f.read(65536)
    if not blob:
        return False, "empty"
    low = blob[:300].lower()
    if b"<html" in low or b"accessdenied" in low or b"<?xml" in low:
        return False, "not_a_matrix"
    if b"\x00" in blob:
        return False, "nul_byte"
    line = blob.split(b"\n", 1)[0]
    if b"\t" not in line:
        return False, "header_not_tab_separated"
    first = line.split(b"\t", 1)[0]
    if layer in {"rna_gene_tumor", "rna_isoform_tumor", "protein_tumor", "phospho_tumor", "cnv"}:
        if first != b"idx":
            return False, f"header_token_{first[:40]!r}"
    if layer == "clinical_meta":
        if b"case" not in first.lower() and b"sample" not in line.lower():
            return False, "meta_header_unrecognized"
    if layer == "phenotype":
        # Samples are rows. The header is idx plus phenotype feature names.
        if first != b"idx":
            return False, "phenotype_header_unrecognized"
    return True, "ok"


def fetch(item: dict) -> dict:
    url = item["url"]
    cancer = url.rstrip("/").split("/")[-2]
    name = item["name"]
    layer = layer_of(name)
    final = ACQ / cancer / name
    partial = ACQ / cancer / (name + ".partial")
    final.parent.mkdir(parents=True, exist_ok=True)
    rec = {
        "cancer": cancer,
        "layer": layer,
        "file_name": name,
        "source_url": url,
        "local_partial": str(partial.relative_to(ROOT)),
        "local_final": str(final.relative_to(ROOT)),
    }
    if final.exists() and final.stat().st_size > 0:
        rec.update({"action": "already_final", "bytes": final.stat().st_size, "ok": True})
        log(rec)
        return rec
    info = head(url)
    rec.update(info)
    log({**rec, "action": "head"})
    if info["http_status"] != "200" or not info["content_length"].isdigit():
        rec.update({"ok": False, "failure": "head_failed"})
        log({**rec, "action": "fail"})
        return rec
    expected = int(info["content_length"])
    if partial.exists() and partial.stat().st_size == expected:
        ok, why = valid_format(partial, layer)
        if ok:
            digest = sha256(partial)
            os.replace(partial, final)
            rec.update({"ok": True, "sha256": digest, "failure": "", "action": "promote_complete_partial", "bytes_on_disk": expected, "expected_bytes": expected})
            log(rec)
            return rec
        rec.update({"ok": False, "failure": why, "action": "complete_partial_bad_format"})
        log(rec)
        return rec
    cmd = [
        "curl", "-L", "--fail",
        "--retry", "8", "--retry-delay", "2", "--retry-all-errors",
        "-C", "-",
        "--connect-timeout", "30",
        "--speed-limit", "2000", "--speed-time", "90",
        "-A", "phase0.5-acquisition",
        "-o", str(partial),
        "-w", "%{http_code}",
        url,
    ]
    started = time.time()
    p = subprocess.run(cmd, capture_output=True)
    http = p.stdout.decode().strip()
    got = partial.stat().st_size if partial.exists() else 0
    rec.update({
        "action": "download",
        "download_http_status": http,
        "curl_rc": p.returncode,
        "bytes_on_disk": got,
        "expected_bytes": expected,
        "elapsed_sec": round(time.time() - started, 1),
    })
    if p.returncode != 0 or got != expected:
        rec.update({"ok": False, "failure": "incomplete_or_curl_error", "stderr_tail": p.stderr.decode("latin1", "replace")[-400:]})
        log(rec)
        return rec
    ok, why = valid_format(partial, layer)
    if not ok:
        rec.update({"ok": False, "failure": why})
        log(rec)
        return rec
    digest = sha256(partial)
    os.replace(partial, final)
    rec.update({"ok": True, "sha256": digest, "failure": ""})
    log(rec)
    return rec


def write_tables(results: list[dict]) -> None:
    # reload events for a complete log even across reruns
    events = []
    if EVENT.exists():
        for line in EVENT.read_text().splitlines():
            if line.strip():
                events.append(json.loads(line))
    finals = {}
    for ev in events:
        if ev.get("action") == "download" and ev.get("ok") and ev.get("sha256"):
            finals[ev["local_final"]] = ev
        if ev.get("action") == "already_final":
            finals.setdefault(ev["local_final"], ev)
    # hash any final that was skipped this run
    rows = []
    for item in select():
        cancer = item["url"].rstrip("/").split("/")[-2]
        name = item["name"]
        final = ACQ / cancer / name
        layer = layer_of(name)
        rel = str(final.relative_to(ROOT))
        ev = finals.get(rel, {})
        exists = final.exists()
        digest = ev.get("sha256") or (sha256(final) if exists else "")
        size = final.stat().st_size if exists else ev.get("bytes_on_disk", 0)
        rows.append({
            "cancer": cancer,
            "layer": layer,
            "file_name": name,
            "source_url": item["url"],
            "source_freeze": "LinkedOmicsKB Release 1 / data_freeze_v1.2_reorganized",
            "retrieval_timestamp_utc": ev.get("timestamp_utc", ""),
            "http_status": ev.get("download_http_status") or ev.get("http_status") or "",
            "content_length": ev.get("expected_bytes") or ev.get("content_length") or "",
            "etag": ev.get("etag", ""),
            "last_modified": ev.get("last_modified", ""),
            "local_path": rel if exists else "",
            "bytes": size if exists else 0,
            "sha256": digest,
            "validated": "yes" if exists and digest else "no",
        })
    fields = list(rows[0].keys())
    with (MAN / "acquisition_manifest.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    log_fields = ["timestamp_utc", "action", "cancer", "layer", "file_name", "source_url", "http_status", "download_http_status", "content_length", "etag", "bytes_on_disk", "expected_bytes", "curl_rc", "ok", "failure", "sha256", "elapsed_sec"]
    with (MAN / "download_log.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=log_fields, extrasaction="ignore")
        w.writeheader()
        for ev in events:
            w.writerow(ev)
    fails = [r for r in rows if r["validated"] != "yes"]
    with (MAN / "download_failures.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["cancer", "layer", "file_name", "source_url", "http_status", "content_length", "bytes"])
        w.writeheader()
        for r in fails:
            w.writerow({k: r.get(k, "") for k in w.fieldnames})
    with (MAN / "file_integrity_report.csv").open("w", newline="") as f:
        cols = ["cancer", "layer", "file_name", "bytes", "content_length", "size_matches_content_length", "sha256", "format_check", "validated"]
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            path = ROOT / r["local_path"] if r["local_path"] else None
            fmt = ""
            if path and path.exists():
                ok, why = valid_format(path, r["layer"])
                fmt = why if ok else why
            match = ""
            if r["content_length"] != "" and r["bytes"]:
                match = "yes" if str(r["bytes"]) == str(r["content_length"]) else "no"
            w.writerow({
                "cancer": r["cancer"],
                "layer": r["layer"],
                "file_name": r["file_name"],
                "bytes": r["bytes"],
                "content_length": r["content_length"],
                "size_matches_content_length": match,
                "sha256": r["sha256"],
                "format_check": fmt,
                "validated": r["validated"],
            })
    lines = []
    for r in rows:
        if r["sha256"] and r["local_path"]:
            lines.append(f"{r['sha256']}  {r['local_path']}")
    (ROOT / "phase0_5" / "SHA256SUMS").write_text("\n".join(lines) + ("\n" if lines else ""))
    print(f"validated {sum(r['validated']=='yes' for r in rows)} / {len(rows)}")


def main() -> None:
    items = select()
    print(f"queue {len(items)}", flush=True)
    # small files first so metadata is available if the large transfer is interrupted
    items.sort(key=lambda h: int(h["bytes"]))
    results = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futs = [pool.submit(fetch, item) for item in items]
        for fut in as_completed(futs):
            rec = fut.result()
            results.append(rec)
            state = "OK" if rec.get("ok") else "FAIL"
            print(f"{state} {rec.get('cancer')} {rec.get('file_name')} {rec.get('failure','')}", flush=True)
    # one retry pass for failures, sequentially, allowing resume
    failed = [item for item in items if not (ACQ / item["url"].rstrip("/").split("/")[-2] / item["name"]).exists()]
    if failed:
        print(f"retry {len(failed)}", flush=True)
        for item in failed:
            rec = fetch(item)
            results.append(rec)
            state = "OK" if rec.get("ok") else "FAIL"
            print(f"RETRY {state} {item['name']} {rec.get('failure','')}", flush=True)
    write_tables(results)


if __name__ == "__main__":
    main()
