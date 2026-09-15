# Validation summary — 15 September 2026

## Completeness

- 13 manuscript figures, 44 labelled panels and two unpanelled figures: **46 entries**.
- 44 imported PDF objects, 42 distinct recovered original filenames, plus native artwork.
- **30 reproduced** (including the explicit native schematic recreation), **9 source code partially identified**, **7 unresolved**.
- All 46 entries have master provenance, a cleaned module or explicit blocker record, and an R Markdown/HTML section.
- Figure 6a/b share one intentional composite; author confirmation is preserved.

## Execution and plot validation

All 13 figure modules were executed in fresh R processes using the available private
input snapshots. Related historical variants were generated where final export code
was missing and remain classified as partial. Missing generators are explicit status
records, not silently skipped panels. Each run records per-panel output/error status,
full local execution objects and session information. The orchestration entry point
was independently checked after the full run.

Validation included visual comparison to native-Keynote embedded originals and
numerical checks where relevant. Fig1b's 57 colored cohort points matched under canvas
scaling with maximum residual 0.00253 PDF points. Fig5a/f model estimates and P values
matched the historical code exactly; Fig5e maximum coefficient difference was
8.53e-14. Fig6 reproduced the 52.2444-year breakpoint (95% CI 45.7802–58.7087),
pre/post slopes 4.6879/15.783 and adjusted Chornobyl slope 4.581181.

Differences due to canvas size, font handling, labels and software rendering were
assessed separately from data/model correspondence. See the panel-specific validation
notes and methodology review for substantive discrepancies. Matching an existing PDF
hash alone was never used as proof of code execution.

## Documentation and publication checks

All R files parse. All 13 Rmd pages render through rmarkdown/Pandoc to self-contained
HTML; source/validation content and relative links were checked. Public HTML contains
no original or reproduced participant plot payloads. Separate local preview reports
were rendered from the same Rmd with generated plots embedded.

The independent validator checks manuscript/Keynote/master/Rmd completeness, real
cleaned-function line locations, every original image object's assignment, output
statuses, local Markdown/HTML links, forbidden private file types and runtime paths.
Its results are in provenance/package_qa.json and provenance/html_content_qa.json.

The reusable DOCX extractor and Keynote inventory parser were tested against the
private originals; the Keynote parser reproduced the slide order, media inventory and
text-object inventory byte-for-byte. All 27 historical R source hashes, source DOCX and
source Keynote hashes remained unchanged. Existing UCEC repository files were unchanged.

## Limits

The original historical R session was not fully frozen. The package records installed
versions and starts from documented analysis intermediates; it does not independently
re-run all raw-read processing or recover missing final plotting scripts. Public HTML
was checked structurally and for working relative dependencies; no interactive browser
behavior is asserted beyond the standard rmarkdown output. Scientific/legend conflicts
remain for author review, as detailed in provenance/methodology_review.md.
