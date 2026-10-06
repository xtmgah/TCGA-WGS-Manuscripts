# Request-only reproduction inputs

The current data revision is `main-figures-1-5-data-review-v1`. Full figure/statistical reproduction requires **35 request-only files**, about **0.763 MB compressed**, in `files/`. These files are excluded from Git and from the Git/code archive.

Contact the corresponding authors listed in the manuscript, identifying the repository/data revision, intended use and input IDs needed. The authors retain these data on **Box and/or Biowulf**, neither of which is a public download location for this project. They determine whether and how inputs can be supplied under applicable source/institutional permissions; a request does not automatically authorize access. There is no public external accession/DOI or alternate hosting service planned.

The local private archive is `.local/bundles/panglioma-request-data.tar.gz`. It is for the authors' approved transfer workflow, **not for Git publication**. This package does not upload it to Box/Biowulf, provide credentials, grant account access or download data automatically.

After receiving an authorized copy, extract from the repository root:

```sh
mkdir -p data/external
tar -xzf /path/to/panglioma-request-data.tar.gz -C data/external
Rscript Rscripts/reproduce.R --check-inputs
Rscript Rscripts/reproduce.R --figures 1-5 --verify --render-reports
```

Alternatively extract elsewhere and pass `--data-dir /path/to/extracted-directory` (the directory containing `files/`). The active manifest verifies the exact compressed and uncompressed checksums. The earlier 142-file bundle is not interchangeable with this minimized revision.

Most files contain the minimal measured values needed for the tests, figures and cohort checks. The single CIC interval is kept here because its source annotation release/licensing provenance was not established. `request_*` is the chosen distribution route, not a blanket statement that every field is formally controlled-access. See [the file review](../redistribution_audit.tsv) and [distribution documentation](../../docs/DATA_DISTRIBUTION.md).
