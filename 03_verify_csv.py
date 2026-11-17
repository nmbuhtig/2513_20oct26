#!/usr/bin/env python3




from __future__ import annotations
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE/'results'


def rows(name):
    path=RESULTS/name
    if not path.is_file():
        raise AssertionError(f'Missing {path}. Run 01_generate_csv.py first.')
    with path.open(newline='',encoding='utf-8') as f:
        return list(csv.DictReader(f))


def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1_048_576),b''):
            h.update(block)
    return h.hexdigest()


def verify():
    checks=[]
    def check(condition,label):
        if not condition:
            raise AssertionError(label)
        checks.append('PASS: '+label)
        print(checks[-1],flush=True)
    baseline=rows('first_zero_baseline.csv')
    sweep=rows('parameter_sweep.csv')
    scan=rows('interference_scan.csv')
    summary=rows('interference_cutoff_summary.csv')
    manifest_path=RESULTS/'run_manifest.json'
    check(manifest_path.exists(),'Reproducibility manifest exists')
    meta=json.loads(manifest_path.read_text(encoding='utf-8'))
    check(meta.get('suite_version')=='3.0','Manifest identifies standalone v3')
    check(len(baseline)==1 and baseline[0]['status']=='computed',
          'First-zero baseline evaluated')
    b=baseline[0]
    check(int(b['prime_power_terms'])==149235,'Baseline enumerates 149,235 prime-power terms')
    check(abs(float(b['rho_re'])-0.499999987757224)<1e-7
          and abs(float(b['rho_im'])-14.134725169080086)<1e-7,
          'First-zero benchmark agrees with historical float64 regression')
    check(float(b['posthoc_error'])<1e-6,'Baseline post-hoc error under 1e-6')
    profile=meta['profile']
    check(len(sweep)==(54 if profile=='quick' else 72),
          'Sweep has expected number of quick/full rows')
    check({r['case'] for r in sweep}=={'first','second'},
          'Sweep has both fixed sampling cases')
    settings=meta['interference_parameters']
    cutoff_list=[int(n) for n in settings['N']]
    n_steps=int(settings['steps'])
    check(len(scan)==len(cutoff_list)*n_steps,
          'Interference scan covers every N-by-t combination')
    check({int(r['N']) for r in scan}==set(cutoff_list),
          'Interference scan has precisely requested N values')
    grouped=defaultdict(list)
    for r in scan:
        grouped[int(r['N'])].append(float(r['t']))
    check(all(len(tt)==n_steps and len(set(tt))==n_steps for tt in grouped.values()),
          'Every cutoff has a complete, nonduplicated sampling grid')
    target_grid=sorted(next(iter(grouped.values())))
    check(all(all(abs(a-b)<1e-11 for a,b in zip(target_grid,sorted(tt)))
              for tt in grouped.values()),
          'All N cutoffs use exactly the same sampling grid')
    check(len(summary)==len(cutoff_list) and
          {int(row['N']) for row in summary}==set(cutoff_list),
          'Cutoff summary has one row for each N')
    check(all(r['status']=='computed' for r in baseline+sweep+scan),
          'All supplied numerical estimates completed')
    check(all(math.isfinite(float(r['rho_re'])) and math.isfinite(float(r['rho_im']))
              for r in baseline+sweep+scan),
          'Reconstructed positions are finite numbers')
    check(all(float(r['condition_kappa0'])>=1-1e-10
              and float(r['condition_kappa1'])>=1-1e-10
              for r in baseline+sweep+scan),
          'Cancellation condition indicators satisfy kappa >= 1')
    check(all(r['reference_role']=='posthoc_only'
              for r in baseline+sweep+scan+summary),
          'All CSV outputs label known-zero information as post-hoc')
    check(all((0<float(r['rho_re'])<1)==(r['estimated_in_critical_strip']=='True')
              for r in baseline+sweep+scan),
          'Critical-strip classification agrees with reconstructed coordinates')
    check(all(r['estimated_in_critical_strip']=='True' or
              'estimate_outside_critical_strip' in r['numerical_warnings']
              for r in baseline+sweep+scan),
          'Every out-of-strip estimate has a PRIME-DATA warning')
    check(all(r['reference_selection_rule'] in
              ('prespecified_case_reference','nearest_reference_to_sampling_ordinate')
              for r in baseline+sweep+scan),
          'Reference-selection method is explicitly recorded')
    check(all(int(r['warning_count'])==int(r['data_warning_count'])+int(r['reference_warning_count'])
              for r in baseline+sweep+scan),
          'Data-only and reference-based warning counts are separate')
    check(all(int(r['evaluated'])+int(r['failed'])==int(r['sample_points'])
              for r in summary),
          'Summary row totals are consistent')
    check(all(int(r['out_of_strip'])==sum(not (0<float(x['rho_re'])<1)
                for x in scan if int(x['N'])==int(r['N']))
              for r in summary),
          'Per-cutoff out-of-strip counts match the underlying data')
    check(all(int(r['sample_points'])==n_steps for r in summary),
          'Summary rows have the requested sampling resolution')

    gaps=defaultdict(list)
    for r in scan:
        gaps[round(float(r['t']),12)].append(float(r['partial_gap']))
    check(all(max(vals)-min(vals)<1e-10 for vals in gaps.values()),
          'Finite-list gap is invariant under prime-cutoff changes')
    hashes=meta['csv_sha256']
    check(set(hashes)=={
          'first_zero_baseline.csv','parameter_sweep.csv',
          'interference_scan.csv','interference_cutoff_summary.csv'},
          'Manifest tracks all four CSV datasets')
    check(all(sha256(RESULTS/name)==digest for name,digest in hashes.items()),
          'All four CSV hashes match the run manifest')
    report=RESULTS/'VERIFICATION_REPORT.txt'
    report.write_text('STANDALONE GAUSSIAN PRIME SUITE v3 -- SOFTWARE VERIFICATION\n\n'+
          '\n'.join(checks)+'\n\nWARNING: these are software and regression checks, '+
          'NOT independent interval certificates or analytic proofs.\n',encoding='utf-8')
    print(f'COMPLETED: {len(checks)} software checks passed; report: {report}',flush=True)

if __name__=='__main__':
    verify()
