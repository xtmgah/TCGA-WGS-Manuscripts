# Provenance and confidence

`figure_inventory.tsv` is the authoritative DOCX-derived set: 13 figures, 44 labelled
panels and two unpanelled figures. IDs preserve the manuscript's lower-case labels;
there is no invented S prefix. `figure_provenance.tsv` is the single joined master.

- `keynote_panel_map.tsv`: manual visual verification plus IWA slide/object/media IDs.
- `keynote_media_objects.tsv`: every original and preview reference; previews are
  explicitly flagged and excluded from figure counts.
- `keynote_native_objects.tsv`: native artwork/text geometry, including Fig.4e.
- `embedded_vs_local_files.tsv`: embedded and archived original hashes and version comparisons.
- `source_catalog.tsv`: 27 historical source-file hashes and line counts.
- `source_code_evidence.tsv`: literal export evidence; dynamic and related source
  blocks are identified separately in the master table.
- `validation.tsv`: executed/compared status and limitations, not inferred from an
  existing PDF or an HTML rendering success.
- `unresolved_figures.tsv`: every non-reproduced panel, with reasons and next evidence needed.

Confidence of the Keynote mapping and confidence of generating-code identification
are separate columns. A filename can be recovered with high confidence while its
R generator remains absent. Original source paths are logical archival identifiers;
code execution uses only packaged scripts and configured private inputs.

## Rebuilding the inventory

With the source manuscript available privately:

```sh
python3 tools/extract_manuscript_inventory.py "$PTC_MANUSCRIPT_DOCX" --output-dir "$PTC_AUDIT_DIR"
python3 tools/trace_source_exports.py "$PTC_HISTORICAL_SOURCE_DIR" "$PTC_EXPORT_EVIDENCE_TSV"
```

The first tool reads DOCX XML, preserving paragraph order including empty paragraphs,
inserted current text and exact figure references. The second provides candidate
source evidence; it does not assign provenance automatically. Follow
[Keynote_Extraction.md](../tools/Keynote_Extraction.md) for IWA decoding and visual checks.

No original manuscript, Keynote, historical input dataset or patient-specific PDF
is distributed. Five TCGA sample identifiers appear in recovered original filenames
because exact source filenames are part of the requested audit; no associated
measurements or row-level datasets are included.
