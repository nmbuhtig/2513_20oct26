#!/usr/bin/env python3






from __future__ import annotations
import argparse
import csv
import hashlib
import subprocess
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
CHECKERS=HERE/'fixed_checkers'
NAMES=[
 'zero_window_certificate_checker_v2.py',
 'prime_only_isolation_checker.py',
 'Phase_IV_extended_zero_free_checker.py',
 'Phase_IV_nonasymptotic_bound_checker.py',
 'Phase_V_sharpened_analytic_bound.py',
 'Phase_VI_H4_checker.py',
 'Phase_VI_optimized_positive_mass_checker.py',
]


def main():
    parser=argparse.ArgumentParser(description='Repeat original fixed arithmetic verification scripts.')
    parser.add_argument('--timeout',type=int,default=600,
                        help='Timeout per checker in seconds, default 600')
    args=parser.parse_args()
    logs=HERE/'logs';logs.mkdir(exist_ok=True)
    results=HERE/'results';results.mkdir(exist_ok=True)
    rows=[]
    for name in NAMES:
        path=CHECKERS/name
        if not path.exists():
            raise SystemExit(f'Missing {path}. Restore the package before running.')
        print(f'RUN: {name}',flush=True)
        h=hashlib.sha256(path.read_bytes()).hexdigest()
        try:
            p=subprocess.run([sys.executable,str(path)],cwd=CHECKERS,
                             capture_output=True,text=True,encoding='utf-8',errors='replace',
                             timeout=args.timeout)
            status='PASS' if p.returncode==0 else 'FAIL'
            content=p.stdout+'\nSTDERR:\n'+p.stderr
        except subprocess.TimeoutExpired:
            status='TIMEOUT'
            content='Timed out. No numerical claim reproduced by this run.\n'
        (logs/(name+'.log.txt')).write_text(content,encoding='utf-8')
        rows.append(dict(checker=name,status=status,script_sha256=h,
                         qualification='fixed arithmetic only; analytic assumptions not independently verified'))
        print('  '+status,flush=True)
    with (results/'fixed_certificate_results.csv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    if any(row['status']!='PASS' for row in rows):
        raise SystemExit('One or more checks failed/timed out. Inspect logs/*.log.txt.')
    print('All seven inherited fixed arithmetic checkers reproduced.\n'
          'No new parameter sweep or analytic proof is certified by this.',flush=True)

if __name__=='__main__':main()
