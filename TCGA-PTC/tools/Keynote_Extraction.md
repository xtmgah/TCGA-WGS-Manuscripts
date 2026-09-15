# Keynote figure provenance

The source deck contains 14 actual slides: a title and 13 figure pages. The committed provenance covers every one of the 46 manuscript panel/figure entries. It records 44 imported PDF objects, 42 unique original filenames, and one native schematic. Original manuscript files, embedded figures, patient-level data, rendered audit images and compiled protobuf modules are not included in this repository.

## Published evidence

- [Panel map](../provenance/keynote_panel_map.tsv): figure/panel, slide, object ID, data ID, exact original filename, source hash, geometry and mapping notes.
- [Media objects](../provenance/keynote_media_objects.tsv): all 44 original PDF references and 44 associated generated previews. Preview rows are marked `is_preview=True` and are not independent figures.
- [Native objects](../provenance/keynote_native_objects.tsv): native shapes, groups, panel labels and other text on the figure slides. Fig. 4e is formed by groups `1317724` and `1317726` and their descendants.
- [Original-file comparisons](../provenance/embedded_vs_local_files.tsv): exact-basename matches against the private source archive, SHA-256, page dimensions, software metadata and comparison metrics. `local_candidate` is a **logical path within the private source archive**, not a required path inside this repository.

Fig. 6a/b share one image object, `tp_age_clock_mutation_acc.pdf`. The author confirmed that a/b were added manually to this intentional composite. Extended Data Fig. 2a/b and Supplementary Fig. 3a/b reuse the same two PDF data IDs in separate image objects. These repetitions remain explicit in the map.

## How the mapping was recovered

1. Hash the original `.key` file. Work on an inspection copy when opening it in Keynote.
2. Treat the original `.key` file as a ZIP archive and read its `Index/*.iwa` members. ZIP entry order is not slide order.
3. Decode each IWA four-byte chunk header (one type byte and three little-endian length bytes). Decompress type-zero chunks with Snappy. Parse the resulting length-delimited `ArchiveInfo` and `MessageInfo` protobuf records.
4. Decode `Index/Metadata.iwa` as `TSP.PackageMetadata`. Preserve `DataInfo.preferred_file_name` verbatim. `DataInfo.file_name` is the internal ZIP member name and may contain an added numeric suffix; it is not the original filename.
5. Decode `Index/Document.iwa`, follow the document's `show` reference, and obtain the ordered presentation nodes from `ShowArchive.slideTree`. In this Keynote version the tree uses repeated node references in field 2; older files use a root node in field 1. Exclude theme/master nodes.
6. Traverse each actual slide's `drawables`, including nested groups. For each `TSD.ImageArchive` (type 3005), retain the object ID, image geometry and every data reference. Prefer the original PDF and mark the generated `-small.png` reference as a preview. Preserve repeated image objects.
7. Resolve native text from `TSWP.ShapeInfoArchive` (type 2011), its `containedStorage`, and placeholder objects (type 7). Use panel-letter coordinates and image geometry as evidence.
8. Open only the inspection copy in Keynote, export a native PDF, and visually check **every page**. Match figure titles, panel letters, plotted content and geometry. The assignment manifest must not automatically pair images and letters solely by array order or proximity.
9. Extract each original `Data/` PDF privately. Compare its SHA-256 against same-basename source PDFs. When bytes differ, compare page size, metadata and rendered plots; a filename match alone is insufficient. Many PDFs contain outlined glyphs and have no extractable text.
10. Assert that every manuscript inventory entry is mapped and every used PDF object is represented. Rehash the original Keynote file to confirm it was not modified.

## Parser and schema dependency

The reusable [rebuild_keynote_inventory.py](rebuild_keynote_inventory.py) accepts
the private source deck, compiled schema directory and a private output directory:

```sh
python tools/rebuild_keynote_inventory.py "$PTC_KEYNOTE" --schema-dir "$PTC_PROTO_MODULES" --output-dir "$PTC_KEYNOTE_AUDIT"
```

It reads the deck without modification and recreates the raw slide/media/native
object inventories. The reviewed panel assignments are preserved in the committed
`keynote_panel_map.tsv`; visual review is still needed after any manuscript/deck change.


The author-side audit scripts are `keynote/parse_keynote.py`, `keynote/build_panel_map.py` and `keynote/compare_pdf_versions.py` within the private `ptc_reproducibility_audit` working directory. `build_panel_map.py` holds the explicit object-ID assignments reviewed against the native PDF. These author-side scripts are an evidence-generation record, not a dependency for running the cleaned figure R modules in this repository.

The schema definitions came from [obriensp/iWorkFileFormat](https://github.com/obriensp/iWorkFileFormat/tree/8575e441beaaaa56f480fdd91721f5bb06d07d43), commit `8575e441beaaaa56f480fdd91721f5bb06d07d43`, under `iWorkFileInspector/iWorkFileInspector/Messages/Proto/`. For current `protoc`, the two legacy required extension declarations in `TSWPArchives.proto` and `TSCHArchives.proto` were changed to optional in the private schema copy before generating Python modules. No figure file was changed by that schema compatibility adjustment.

The extraction environment used Python with `python-snappy`, `protobuf`, `PyMuPDF`, `Pillow` and `numpy`, plus Keynote Creator Studio for the independent native PDF export. An equivalent implementation can follow the field-level procedure above and regenerate the same IDs, filenames and hashes from the exact source deck.

## Reproduction boundaries

Recovered filenames establish the figure-to-asset link. They do not prove that a historical script generates the final plotted values. For example, the final Fig. 2d PDF contains FDR values while the same-name historical source PDF contains raw P values; the final Fig. 4c PDF is pooled while its same-name source PDF was overwritten with a cohort-stratified plot. Those differences are retained in the source audit.

Fig. 4e has no historical R generator because it was authored in Keynote. [Fig4_schematic.R](../Rscripts/Fig4_schematic.R) recreates its shapes and labels and explicitly returns `recreated native schematic`. This is a documented reconstruction of a schematic, not a claim to have recovered original analysis code.
