# Fresh R/Python environment validation

This milestone tests a newly installed **macOS arm64** environment. R 4.5.1, Python 3.12.13, Pandoc 3.12, the R packages, Python packages and their native libraries were downloaded into a new, isolated directory. It uses the existing macOS operating system and registered Roboto Condensed fonts. It is not a new operating-system installation. A second-platform test was deferred at the author's request; no container was built or tested.

The final build passed in **101.21 seconds**, excluding installation and HTML reports. All **141 numerical checks (1,064 entries)**, **30 AC2 table checks**, **14 Python tests**, **12 statistical R tests** and **two entrypoint-isolation checks** pass. All five figures retain the previous reviewed fingerprints at 150 and 300 dpi; all six HTML reports render. The executable Figure 1 R Markdown entrypoint also passes from a new font cache and reproduces those pixels.

## Isolation and evidence

The test uses a fresh extraction of the Git/code archive plus the 35 authorized request-only inputs. The fresh extraction has no convenience runtime configuration or copied Python dependency directory. R and Python package-location checks require every installed package to resolve inside the new environment, with Python's user site disabled.

The build runs under a macOS sandbox that denies the original analysis project, the primary reproduction package, the previous R framework, Homebrew, the prior bundled Python runtime and manuscript reference PDFs. The new environment and isolated package are allowed. Verification runs separately with access to the isolated reference copies. The full build, numerical regression tests, input/distribution checks, PDF tests and six HTML reports are exercised in this environment.

The aggregate results are recorded in [environment_checks.json](environment_checks.json). [environment/validated-macos/](../environment/validated-macos/) records the installed R/Python packages, Conda package builds and download checksums, R session information, native graphics/BLAS libraries and Pandoc version. These are environment inventories, not a claim that every future dependency resolution or another operating system has been tested.

The test environment and detailed logs remain under ignored `.local/platform_validation/`. No original analysis files or system R/Python installation were modified. Request-only data remain excluded from Git.

## Compatibility fixes found by this test

1. **Cairo blank-space fonts.** Cairo 1.18.4 in the fresh environment represents some spaces in separate empty Type3 fonts. The initial Figure 2 layout stopped because a title was extracted as separate words. `functions/python/cairo_spaces.py` restores only proven-empty space glyphs to adjacent embedded TrueType font subsets, preserving glyph advances and drawn outlines. The composition step applies this to newly generated panels. Regression checks compare rendered pixels at 150 and 300 dpi and extracted labels; painted Type3 glyphs are never treated as spaces. Existing PDFs without this representation are unchanged.
2. **Explicit Python selection.** `PANGLIOMA_PYTHON` now selects that interpreter's installed packages without automatically adding the convenience `.local/python_deps` directory. An explicit caller-provided `PYTHONPATH` remains the caller's choice; the clean test unsets it. The default convenience runtime still works. Two focused R checks cover both paths and environment restoration.
3. **Dependency provenance.** The original `hrbrthemes` 0.9.2 came from GitHub commit `d3fd02949fc201c6db616ccaffbb9858aec6fd2b`, rather than a CRAN release. The fresh environment uses CRAN 0.9.3; the `theme_ipsum_rc` function used here has identical code in both downloaded sources (only a documentation comment differs). The source URL and checksum are recorded in [hrbrthemes-source.json](../environment/hrbrthemes-source.json). GenomicRanges was removed from the direct installation requirements because the minimized CIC coordinate table no longer needs it. Historical snapshots retain their original records.
4. **Font-cache independence.** The direct R Markdown rebuild exposed a previously unrecorded dependency on system-installed Roboto Condensed 3.008. Depending on cache order, Matplotlib could instead choose the bundled older regular face and fail the legend-layout check. The exact additional fonts are now bundled unchanged with their Apache 2.0 license notice and checksums. Python selects only the reviewed bundled faces for this family; a test introduces a competing font entry and verifies all four selected styles. The macOS environment check also verifies the system-selected regular, bold and italic font files against the bundled variable-font files. The original static fonts remain for PDF overlays and the plotting styles that used them.

No statistical method, tolerance, correction family, input value or visual-review threshold was changed.

## Repeat the tested setup on macOS arm64

Use [Micromamba's official installation instructions](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) if it is not already available. Run the commands below from this standalone repository. Register the **variable regular and variable italic** `.ttf` files in `functions/fonts/system/` with macOS Font Book for R/Cairo. `environment/fonts.conf` also directs Fontconfig to the supplied system-font set. Python loads its exact font files directly and does not depend on their being installed system-wide.

```sh
micromamba create -y -p "$PWD/.local/clean-env" \
  --strict-channel-priority -f environment/clean-macos.yml

export PATH="$PWD/.local/clean-env/bin:$PATH"
export PANGLIOMA_PYTHON="$PWD/.local/clean-env/bin/python"
export PANGLIOMA_RSCRIPT="$PWD/.local/clean-env/bin/Rscript"
export PYTHONNOUSERSITE=1
export FONTCONFIG_FILE="$PWD/environment/fonts.conf"
export XDG_CACHE_HOME="$PWD/.local/cache"
unset PYTHONPATH PYTHONHOME R_HOME R_LIBS
export R_LIBS_USER="$PWD/.local/clean-env/lib/R/library"
export R_LIBS_SITE="$R_LIBS_USER"
export R_ENVIRON_USER=/dev/null
export R_PROFILE_USER=/dev/null

python -m pip install --no-cache-dir -r requirements.txt
mkdir -p .local/downloads
curl --fail --location \
  https://cloud.r-project.org/src/contrib/hrbrthemes_0.9.3.tar.gz \
  -o .local/downloads/hrbrthemes_0.9.3.tar.gz
python -c 'import hashlib,pathlib; p=pathlib.Path(".local/downloads/hrbrthemes_0.9.3.tar.gz"); assert hashlib.sha256(p.read_bytes()).hexdigest()=="931d7f65524d090dfff3950e4dd0c850a8bf8bd6b406d65e8683d4d93742729f"'
R CMD INSTALL --no-multiarch .local/downloads/hrbrthemes_0.9.3.tar.gz

Rscript --vanilla Rscripts/check_environment.R
python -m pip check
Rscript --vanilla Rscripts/reproduce.R --figures 1-5 --verify --render-reports
Rscript --vanilla tests/test_statistics.R
Rscript --vanilla tests/test_entrypoint.R
python -m unittest discover -s tests -p 'test_*.py'
```

If CRAN has moved the theme source, use the archived URL recorded in `hrbrthemes-source.json` and require the same checksum. The environment checker compares against the original direct-package snapshot, so its message about `hrbrthemes` 0.9.3 differing from 0.9.2 is expected for this recipe. `requirements.txt` pins direct Python dependencies; the complete installed dependency inventory is retained under `environment/validated-macos/`.

The figure command requires all 105 active inputs. Follow [the author-request instructions](../data/external/README.md) for the 35 inputs omitted from Git. The environment setup does not retrieve or publish those data.
