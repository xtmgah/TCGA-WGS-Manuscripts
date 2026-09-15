#!/usr/bin/env python3
"""Independent, standard-library-only PTC completeness and publication audit.

Run from any directory: python tools/validate_package.py [--root PTC_DIRECTORY]
JSON is printed to stdout. Exit 1 means a required audit check failed.
Historical source paths in provenance are logical metadata, not runtime inputs.
"""
from __future__ import annotations
import argparse
import csv
import html.parser
import json
import re
import sys
from collections import Counter
from pathlib import Path
from urllib.parse import unquote, urlsplit

EXPECTED_LEAVES = 46
EXPECTED_PDF_OBJECTS = 44
EXPECTED_FILENAMES = 42
ALLOWED_STATUSES = {
    'reproduced', 'code identified but reproduction blocked',
    'source code partially identified', 'unresolved'}
FIGURE = re.compile(r'(?:(Extended\s+Data|Supplementary)\s+)?'
                    r'(?:Figure|Fig\.?)\s*S?(\d+)(?:\s*([a-zA-Z])\b)?', re.I)
ABSOLUTE_PATH = re.compile(r'(?<![\w./])/(?:Users|Volumes|data|scratch|tmp|home)/[^\s\"\'`<>)]*')
KNOWN_USERNAME = re.compile(r'\b(?:zhangt8)\b', re.I)
CODE_EXTENSIONS = {'.r', '.py', '.sh', '.bash', '.js'}
FORBIDDEN_ENDINGS = ('.rdata', '.rda', '.rds', '.key', '.docx', '.bam', '.cram',
                     '.fastq', '.fastq.gz', '.fq', '.fq.gz', '.vcf', '.vcf.gz')


def canonical(value):
    m = FIGURE.search(str(value))
    if not m:
        return str(value).strip()
    prefix = m.group(1)
    prefix = ('Extended Data ' if prefix.lower().startswith('extended') else 'Supplementary ') if prefix else ''
    return f'{prefix}Fig. {int(m.group(2))}{(m.group(3) or "").lower()}'


def row_id(row):
    for key in ('manuscript_id', 'manuscript_figure_panel_id', 'canonical_figure_panel_id'):
        if row.get(key):
            return canonical(row[key])
    return canonical(row.get('figure_id', '') + row.get('panel_id', ''))


def load_tsv(root, rel, failures):
    path = root / rel
    if not path.is_file():
        failures.append({'check': 'required_file', 'path': rel, 'reason': 'missing'})
        return []
    with path.open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def source_code(text, suffix):
    if suffix.lower() == '.rmd':
        chunks = re.findall(r'^```\{(?:r|python|bash|sh)\b[^\n]*\}\s*\n(.*?)^```\s*$', text, re.M | re.S)
        text = '\n'.join(chunks)
    # Comments are source-history descriptions, not executable path use.
    return '\n'.join(line for line in text.splitlines() if not re.match(r'^\s*#', line))


class LinkParser(html.parser.HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
    def handle_starttag(self, tag, attrs):
        for key, value in attrs:
            if key in ('href', 'src') and value:
                self.links.append(value)


def document_links(path, text):
    if path.suffix.lower() == '.html':
        parser = LinkParser()
        parser.feed(text)
        return parser.links
    # Omit fenced source code before parsing Markdown links.
    prose = re.sub(r'^```[^\n]*\n.*?^```\s*$', '', text, flags=re.M | re.S)
    return [match[1].strip('<>') for match in re.findall(r'(!?\[[^\]]*\])\(([^\s)]+)(?:\s+[^)]*)?\)', prose)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    args = ap.parse_args()
    root = args.root.resolve()
    failures, warnings = [], []
    inventory = load_tsv(root, 'provenance/figure_inventory.tsv', failures)
    master = load_tsv(root, 'provenance/figure_provenance.tsv', failures)
    keynote = load_tsv(root, 'provenance/keynote_panel_map.tsv', failures)
    media = load_tsv(root, 'provenance/keynote_media_objects.tsv', failures)
    tables = {'manuscript_inventory': inventory, 'master_provenance': master, 'keynote_panel_map': keynote}
    sets = {}
    for name, rows in tables.items():
        ids = [row_id(row) for row in rows]
        sets[name] = set(ids)
        duplicates = {k: v for k, v in Counter(ids).items() if v > 1}
        if len(ids) != EXPECTED_LEAVES or len(set(ids)) != EXPECTED_LEAVES or duplicates:
            failures.append({'check': 'leaf_count', 'table': name, 'rows': len(ids),
                             'unique': len(set(ids)), 'duplicates': duplicates})
    baseline = sets['manuscript_inventory']
    for name, ids in sets.items():
        if ids != baseline:
            failures.append({'check': 'panel_set_equality', 'table': name,
                             'missing': sorted(baseline-ids), 'extra': sorted(ids-baseline)})
    mandatory = {'figure_id', 'panel_id', 'keynote_slide', 'original_plot_filename',
                 'source_script', 'source_code_location', 'input_files', 'dependencies',
                 'output_file', 'reproducibility_status', 'notes'}
    if master:
        missing = sorted(mandatory-set(master[0]))
        if missing:
            failures.append({'check': 'master_required_columns', 'missing': missing})
    status_counts = Counter(row.get('reproducibility_status', '') for row in master)
    invalid_status = sorted(set(status_counts)-ALLOWED_STATUSES)
    if invalid_status:
        failures.append({'check': 'provenance_status', 'invalid': invalid_status})
    imported = [row for row in keynote if row.get('asset_type') == 'imported PDF']
    pdf_objects = {row.get('keynote_object_id', '') for row in imported}
    filenames = {row.get('original_plot_filename', '') for row in imported}
    originals = [row for row in media if str(row.get('is_preview', '')).lower() == 'false']
    media_ids = {row.get('object_id', '') for row in originals}
    if len(pdf_objects) != EXPECTED_PDF_OBJECTS or len(filenames) != EXPECTED_FILENAMES:
        failures.append({'check': 'original_asset_counts', 'pdf_objects': len(pdf_objects),
                         'unique_filenames': len(filenames)})
    if pdf_objects != media_ids:
        failures.append({'check': 'all_media_objects_assigned', 'unassigned': sorted(media_ids-pdf_objects),
                         'unknown_assignment': sorted(pdf_objects-media_ids)})
    kmap = {row_id(row): row for row in keynote}
    a, b = kmap.get('Fig. 6a', {}), kmap.get('Fig. 6b', {})
    composite_ok = bool(a and b and a.get('keynote_object_id') == b.get('keynote_object_id') and
                        a.get('original_plot_filename') == b.get('original_plot_filename'))
    if not composite_ok:
        failures.append({'check': 'intentional_fig6_composite', 'reason': 'a/b must share one original object'})
    native = kmap.get('Fig. 4e', {})
    if native.get('asset_type') != 'native Keynote schematic' or native.get('original_plot_filename'):
        failures.append({'check': 'native_fig4e', 'reason': 'must retain native authoring and no invented filename'})
    # Verify the integrated master has not drifted from its reviewed Keynote map.
    for row in master:
        kid = row_id(row)
        krow = kmap.get(kid, {})
        for field in ('keynote_slide', 'original_plot_filename'):
            if krow and str(row.get(field, '')) != str(krow.get(field, '')):
                failures.append({'check': 'master_keynote_agreement', 'panel': kid, 'field': field,
                                 'master': row.get(field), 'keynote': krow.get(field)})
    rmd_files = sorted((root/'Rscripts').glob('*.Rmd'))
    section_ids, section_locations, section_bodies = set(), {}, {}
    for path in rmd_files:
        document = path.read_text(encoding='utf-8')
        lines = document.splitlines()
        headings = [i for i, line in enumerate(lines) if re.match(r'^#{1,6}\s+', line)]
        for h, index in enumerate(headings):
            line = lines[index]
            if not re.match(r'^#{1,6}\s+', line):
                continue
            for match in FIGURE.finditer(line):
                pid = canonical(match.group())
                if pid in baseline:
                    section_ids.add(pid)
                    section_locations.setdefault(pid, []).append(str(path.relative_to(root)))
                    end = headings[h+1] if h+1 < len(headings) else len(lines)
                    section_bodies.setdefault(pid, []).append((path, '\n'.join(lines[index:end])))
    if section_ids != baseline:
        failures.append({'check': 'rmd_panel_sections', 'missing': sorted(baseline-section_ids),
                         'extra': sorted(section_ids-baseline)})
    for row in master:
        pid = row_id(row)
        cleaned = row.get('cleaned_script', '')
        path = root/cleaned
        if not cleaned or not path.is_file():
            failures.append({'check': 'cleaned_code_file', 'panel': pid, 'path': cleaned, 'reason': 'missing'})
            continue
        script = path.read_text(encoding='utf-8')
        location = row.get('cleaned_code_location', '')
        cited = re.search(r':(\d+)\s*\(([^)]+)\)', location)
        if cited:
            number, function = int(cited.group(1)), cited.group(2)
            definitions = [i for i, line in enumerate(script.splitlines(), 1)
                           if re.search(r'\b'+re.escape(function)+r'\s*<-\s*function', line)]
            if number not in definitions:
                failures.append({'check': 'cleaned_code_locator', 'panel': pid, 'location': location,
                                 'actual_function_lines': definitions})
        linked = False
        for rmd, body in section_bodies.get(pid, []):
            for target in document_links(rmd, body):
                parsed = urlsplit(target)
                if not parsed.scheme and (rmd.parent/unquote(parsed.path)).resolve() == path.resolve():
                    linked = True
        if not linked:
            failures.append({'check': 'panel_to_cleaned_source_link', 'panel': pid, 'expected_script': cleaned,
                             'reason': 'panel section does not link its recorded cleaned script'})
        if (row.get('reproducibility_status') == 'source code partially identified' and
            not any(row.get(k, '').strip() for k in ('source_script', 'source_code_location', 'upstream_scripts'))):
            warnings.append({'check': 'partial_code_evidence', 'panel': pid,
                             'reason': 'partial-code classification has no concrete source or upstream code location'})
    broken_links, runtime_findings, metadata_findings, prohibited = [], [], [], []
    text_extensions = {'.md', '.rmd', '.html', '.r', '.py', '.sh', '.json', '.tsv', '.csv', '.txt', '.yml', '.yaml'}
    for path in sorted(root.rglob('*')):
        if not path.is_file():
            continue
        rel = str(path.relative_to(root))
        low = path.name.lower()
        if low.endswith(FORBIDDEN_ENDINGS) or (re.search(r'(?:_scna|_dpclust|tcga-[a-z0-9]{2}-[a-z0-9]{4})', low)
                                               and path.suffix.lower() in {'.pdf', '.png', '.jpg', '.svg'}):
            prohibited.append(rel)
        if path.suffix.lower() not in text_extensions:
            continue
        try:
            text = path.read_text(encoding='utf-8')
        except UnicodeError:
            continue
        # The auditor necessarily contains the patterns it audits; no runtime
        # file paths or user data are used by those regex constants.
        if path.resolve() != Path(__file__).resolve():
            code = source_code(text, path.suffix) if path.suffix.lower() in CODE_EXTENSIONS | {'.rmd'} else ''
            for pattern, label in ((ABSOLUTE_PATH, 'absolute_personal_or_cluster_path'), (KNOWN_USERNAME, 'personal_username')):
                for m in pattern.finditer(code):
                    runtime_findings.append({'path': rel, 'kind': label, 'match': m.group()})
                if not code:
                    hits = sorted(set(pattern.findall(text)))
                    if hits:
                        metadata_findings.append({'path': rel, 'kind': label, 'matches': hits[:10],
                                                  'classification': 'source metadata or documentation; review, not runtime use'})
        if path.suffix.lower() in {'.md', '.rmd', '.html'}:
            for target in document_links(path, text):
                parsed = urlsplit(target)
                if parsed.scheme or parsed.netloc or not parsed.path or target.startswith('#'):
                    continue
                decoded = unquote(parsed.path)
                dest = (path.parent/decoded).resolve()
                if not dest.exists():
                    broken_links.append({'document': rel, 'target': target, 'reason': 'target does not exist'})
                elif not dest.is_relative_to(root):
                    broken_links.append({'document': rel, 'target': target, 'reason': 'link leaves PTC directory'})
        if rel.startswith('data/') and path.suffix.lower() in {'.tsv', '.csv'}:
            header = text.splitlines()[0] if text.splitlines() else ''
            if re.search(r'\b(?:Tumor_Barcode|Analysis_Barcode|patient_id|sample_id)\b', header, re.I):
                prohibited.append(rel+' (sample-level data columns)')
    if broken_links:
        failures.append({'check': 'relative_document_links', 'findings': broken_links})
    if runtime_findings:
        failures.append({'check': 'runtime_portability', 'findings': runtime_findings})
    if prohibited:
        failures.append({'check': 'private_data_publication', 'files': prohibited})
    if metadata_findings:
        warnings.append({'check': 'metadata_path_review', 'findings': metadata_findings})
    report = {
        'passed': not failures,
        'counts': {'manuscript_leaves': len(baseline), 'provenance_rows': len(master),
                   'keynote_rows': len(keynote), 'rmd_files': len(rmd_files),
                   'html_files': len(list((root/'Rscripts').glob('*.html'))),
                   'rmd_leaf_sections': len(section_ids), 'pdf_objects': len(pdf_objects),
                   'unique_original_filenames': len(filenames),
                   'preview_references': len(media)-len(originals)},
        'status_counts': dict(status_counts),
        'intentional_fig6_shared_object_verified': composite_ok,
        'native_fig4e_verified': native.get('asset_type') == 'native Keynote schematic',
        'rmd_section_locations': section_locations,
        'failures': failures, 'warnings': warnings,
        'limits': ['This audit checks completeness, links, file types and runtime path strings; it does not prove statistical equivalence.',
                   'Historical source-script paths in provenance are logical metadata, so they are not resolved as runtime files.',
                   'Status counts are derived from the current provenance table; successful process exit alone is not reproduction evidence.']}
    print(json.dumps(report, indent=2))
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
