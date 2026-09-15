#!/usr/bin/env python3
"""Extract figure references and legend-defined panels from the final PTC DOCX.

Uses the Python standard library. Never modifies the source DOCX. Paragraph
numbers are one-based w:p order under w:body in word/document.xml, including
empty paragraphs. Includes inserted w:t text; excludes deleted w:delText.
"""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
REF_RE = re.compile(r"\b(?P<prefix>(?:(?:Extended Data|Supplementary)\s+)?Fig(?:ure)?s?\.?)\s+(?P<number>\d+)(?P<panels>[a-z](?:\s*[,–-]\s*[a-z])*)?")
LEGEND_RE = re.compile(r"^(?P<figure>(?:(?:Extended Data|Supplementary)\s+)?Fig\.\s+\d+)\s*[:|]")
# New panel descriptions begin after a full stop, never at an inline 'In c,d'.
PANEL_RE = re.compile(r"(?:(?<=\.)\s+|^)(?P<labels>[a-z](?:[–-][a-z]|,[a-z])*)[,]\s+")


def expand_panels(value):
    if not value:
        return []
    result = []
    for term in re.split(r",\s*", value):
        term = term.strip()
        bounds = re.split(r"\s*[–-]\s*", term)
        result.extend([chr(n) for n in range(ord(bounds[0]), ord(bounds[1]) + 1)] if len(bounds) == 2 else bounds)
    return list(dict.fromkeys(result))


def write_tsv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: json.dumps(row[field], ensure_ascii=False) if isinstance(row.get(field), (list, dict)) else row.get(field, "") for field in fields})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("docx", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    source = args.docx.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    with ZipFile(source) as package:
        xml_bytes = package.read("word/document.xml")
        root = ET.fromstring(xml_bytes)
        auxiliary_parts = []
        for name in ["word/footnotes.xml", "word/endnotes.xml", "word/comments.xml"]:
            if name in package.namelist():
                aux = ET.fromstring(package.read(name))
                aux_text = " ".join(e.text or "" for e in aux.findall(".//w:t", NS))
                auxiliary_parts.append({"part": name, "figure_references": [m.group(0) for m in REF_RE.finditer(aux_text)], "use": "Excluded from manuscript inventory; annotations are not main-text or legend evidence."})
    paragraphs = []
    references = []
    legends = {}
    for index, para in enumerate(root.findall(".//w:body//w:p", NS), 1):
        value = "".join(e.text or "" for e in para.findall(".//w:t", NS))
        paragraph = {"paragraph_id": f"p{index:04d}", "paragraph_index": index, "text": value}
        paragraphs.append(paragraph)
        head = LEGEND_RE.match(value)
        role = "legend" if head else "body"
        for match in REF_RE.finditer(value):
            figure_id = f"{match.group('prefix')} {match.group('number')}"
            references.append({"reference_id": f"ref{len(references)+1:03d}", "paragraph_id": paragraph["paragraph_id"], "role": role, "exact_reference": match.group(0), "start_character_0based": match.start(), "end_character_exclusive": match.end(), "figure_id": figure_id, "raw_panel_expression": match.group("panels") or "", "expanded_panels": expand_panels(match.group("panels")), "paragraph_text": value})
        if not head:
            continue
        figure = head.group("figure")
        markers = list(PANEL_RE.finditer(value))
        descriptions = []
        for j, marker in enumerate(markers):
            descriptions.append({"exact_panel_expression": marker.group("labels"), "panels": expand_panels(marker.group("labels")), "start_character_0based": marker.start("labels"), "end_character_exclusive": markers[j+1].start() if j+1 < len(markers) else len(value), "text": value[marker.start("labels"):markers[j+1].start() if j+1 < len(markers) else len(value)].strip()})
        if figure in legends:
            raise ValueError(f"Duplicate legend heading for {figure}")
        legends[figure] = {"figure_id": figure, "paragraph_id": paragraph["paragraph_id"], "paragraph_index": index, "legend_heading_exact": head.group(0), "legend_text": value, "panels": sorted(set(p for block in descriptions for p in block["panels"])), "panel_description_blocks": descriptions}

    inventory = []
    for figure, legend in legends.items():
        category = "extended_data" if figure.startswith("Extended") else "supplementary" if figure.startswith("Supplementary") else "main"
        body_refs = [r for r in references if r["figure_id"] == figure and r["role"] == "body"]
        for panel in legend["panels"] or [""]:
            exact_id = figure + panel
            description_blocks = [d for d in legend["panel_description_blocks"] if panel in d["panels"]]
            inventory.append({"inventory_id": f"ms{len(inventory)+1:03d}", "figure_id": figure, "panel_id": panel, "manuscript_figure_panel_id": exact_id, "canonical_figure_panel_id": figure.replace("Fig.", "Figure") + (panel.upper() if panel else ""), "figure_category": category, "panel_basis": "legend_defined" if panel else "legend_unpanelled", "legend_paragraph_id": legend["paragraph_id"], "legend_paragraph_index": legend["paragraph_index"], "exact_legend_panel_expressions": [d["exact_panel_expression"] for d in description_blocks], "legend_snippet": " | ".join(d["text"] for d in description_blocks) if panel else legend["legend_text"], "legend_text": legend["legend_text"], "explicit_panel_body_references": [r["reference_id"] for r in body_refs if panel and panel in r["expanded_panels"]], "figure_level_body_references": [r["reference_id"] for r in body_refs if not r["expanded_panels"]], "all_body_reference_paragraphs": sorted(set(r["paragraph_id"] for r in body_refs)), "keynote_mapping_status": "Not assessed by manuscript extractor", "source_docx_sha256": source_sha})

    body_only = [r for r in references if r["role"] == "body" and r["figure_id"] not in legends]
    undefined_panels = [r for r in references if r["role"] == "body" and r["figure_id"] in legends and set(r["expanded_panels"]) - set(legends[r["figure_id"]]["panels"])]
    metadata = {"source_docx": str(source), "source_docx_sha256": source_sha, "document_xml_sha256": hashlib.sha256(xml_bytes).hexdigest(), "paragraph_count_including_empty": len(paragraphs), "tracked_insertions": len(root.findall(".//w:ins", NS)), "tracked_deletions": len(root.findall(".//w:del", NS)), "extraction_policy": "Current w:t text including insertions. Body paragraph order includes empty paragraphs. w:delText, comments, headers, footers and auxiliary notes are not manuscript inventory evidence.", "auxiliary_parts_checked": auxiliary_parts, "figure_count": len(legends), "leaf_entry_count": len(inventory), "labelled_panel_count": sum(bool(r["panel_id"]) for r in inventory), "unpanelled_figure_count": sum(not r["panel_id"] for r in inventory), "counts_by_category": {category: {"figures": sum((figure.startswith("Extended") if category == "extended_data" else figure.startswith("Supplementary") if category == "supplementary" else figure.startswith("Fig.")) for figure in legends), "leaf_entries": sum(row["figure_category"] == category for row in inventory)} for category in ["main", "extended_data", "supplementary"]}, "body_only_references_without_legend": body_only, "body_references_to_undefined_panels": undefined_panels, "legends_without_body_citation": [figure for figure in legends if not any(r["figure_id"] == figure and r["role"] == "body" for r in references)]}
    payload = {"metadata": metadata, "inventory": inventory, "legends": list(legends.values()), "references": references}
    (output / "manuscript_inventory.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (output / "manuscript_paragraphs.json").write_text(json.dumps(paragraphs, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (output / "manuscript_paragraphs.txt").write_text("\n".join(p["paragraph_id"] + "\t" + p["text"] for p in paragraphs) + "\n", encoding="utf-8")
    write_tsv(output / "manuscript_figure_panel_inventory.tsv", inventory, list(inventory[0]))
    write_tsv(output / "manuscript_figure_references.tsv", references, list(references[0]))
    (output / "inventory_summary.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
