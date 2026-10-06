# Input contract and distribution

The active package has **105 compressed inputs totaling 16.09 MB**. No raw FASTQ, BAM/CRAM or large upstream workspace is included. Git is the project's only public data location.

| Location | Files | Compressed size | Distribution |
|---|---:|---:|---|
| `processed/files/` | 70 | 15.33 MB | Reviewed aggregate figure-support results and project styles, included in Git |
| `external/files/` | 35 | 0.763 MB | Minimal individual/feature-level inputs plus one reference with unverified provenance; request from the authors; ignored by Git |

A Git checkout does not contain every input needed for recalculation. The authors retain request-only data on Box and/or Biowulf; neither is publicly accessible for this project. There is no public external download/accession to configure. Follow [external/README.md](external/README.md) to request and install the matching private bundle.

`manifest.tsv` is the active input contract. `id` is stable; `logical_path` is the destination inside the disposable output workspace; `path` locates its compressed bytes. `tier` determines Git versus request-only resolution. `sha256` and `compressed_sha256` identify the uncompressed and gzip bytes. `rows` and JSON-encoded `columns` describe the exact table schema. `redistribution` records the reviewed route. All inputs are required for a complete or single-figure build because the current helpers share statistics and geometry.

`redistribution_audit.tsv` reviews all 169 former inputs, including 64 retired files. `column_audit.tsv` records field-level reductions, and `input_review.json` defines their routes and transformations. See [the distribution review](../docs/DATA_DISTRIBUTION.md) for evidence, limitations and request handling. Numbering gaps reflect retired inputs, not missing data.

The Git set contains cohort summaries/model realizations, not individual participant records. Model realization IDs and genomic event labels are not patient identifiers. Detailed individual/segment measurements remain private even when a table has no explicit barcode. Full upstream ordering/timing fits, DESeq2 and GSEA analyses remain frozen; tractable tests and current AC2 summaries are regenerated from their retained inputs.

Missing files or checksum mismatches stop execution before fitting or drawing. Decompression and generated individual-level tables stay inside ignored `outputs/`. The release contract rejects unreviewed files and checks Git data for participant identifiers and individual-record columns. This technical check does not confer source-data permissions.

For a Git-only checkout:

```sh
python3 functions/python/run.py --check-public-inputs
```

After obtaining the private inputs, use `Rscript Rscripts/reproduce.R --check-inputs` and then run the normal build. `functions/python/package.py` creates separate Git/code and **private author-request** archives under ignored `.local/bundles/`; it never uploads data or invokes Git. The public archive includes these audits and only the 70 reviewed Git inputs.
