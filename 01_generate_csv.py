#!/usr/bin/env python3













from __future__ import annotations

import argparse
import cmath
import csv
import hashlib
import json
import math
import platform
import statistics
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent



REFERENCE_GAMMAS = (
    14.134725141734693790457252,
    21.022039638771554992628480,
    25.010857580145688763213791,
    30.424876125859513210311898,
    32.935061587739189690662369,
    37.586178158825671257217764,
    40.918719012147495187398127,
    43.327073280914999519496122,
    48.005150881167159727942473,
    49.773832477672302181916785,
    52.970321477714460644147297,
    56.446247697063394804367760,
)
CASE_ORDINATES = {'first': 14.0, 'second': 21.0}
FLOAT_EPS = 2.220446049250313e-16


def make_prime_powers(max_n: int):
    pass                                                                  
    if not 2 <= max_n <= 15_000_000:
        raise ValueError('Prime cutoff must be between 2 and 15,000,000.')
    sieve = bytearray(b'\x01') * (max_n + 1)
    sieve[:2] = b'\x00\x00'
    for p in range(2, math.isqrt(max_n) + 1):
        if sieve[p]:
            sieve[p*p:max_n+1:p] = b'\x00' * ((max_n-p*p)//p + 1)
    entries = []
    for p in range(2, max_n + 1):
        if not sieve[p]:
            continue
        lam = math.log(p)
        power = p
        while power <= max_n:
            entries.append((power, math.log(power), lam))
            if power > max_n // p:
                break
            power *= p
    entries.sort(key=lambda entry: entry[0])
    return entries


def half_erfc(z: float) -> float:
    return 0.5 * math.erfc(z)


def weight_increment(z0: float, z1: float) -> float:
    pass                                             
    if z0 < 0 and z1 < 0:
        return 0.5 * (math.erfc(-z0) - math.erfc(-z1))
    return 0.5 * (math.erfc(z1) - math.erfc(z0))


def prime_differences(entries, nmax, x, sigma, t, alpha, h):
    pass                                                                      
    if not (nmax >= 2 and x > 1 and 0 < sigma < 1 and alpha > 0 and h > 0):
        raise ValueError('Require N>=2, X>1, 0<sigma<1, alpha>0, h>0.')
    y = math.log(x)
    y0, y1, y2 = y, y+h, y+2*h
    a0 = 1.0 / (2*alpha*math.sqrt(y0))
    a1 = 1.0 / (2*alpha*math.sqrt(y1))
    a2 = 1.0 / (2*alpha*math.sqrt(y2))
    re0, im0, re1, im1, abs0, abs1 = [], [], [], [], [], []
    count = 0
    for n, ln, lam in entries:
        if n > nmax:
            break
        count += 1
        z0, z1, z2 = (ln-y0)*a0, (ln-y1)*a1, (ln-y2)*a2
        g0, g1 = weight_increment(z0,z1),weight_increment(z1,z2)
        amp = -lam * math.exp(-sigma*ln)
        phase = t*ln
        c, ss = math.cos(phase), math.sin(phase)
        v0 = amp*g0
        v1 = amp*g1
        re0.append(v0*c)
        im0.append(-v0*ss)
        re1.append(v1*c)
        im1.append(-v1*ss)
        abs0.append(abs(v0))
        abs1.append(abs(v1))
    d0 = complex(math.fsum(re0), math.fsum(im0))
    d1 = complex(math.fsum(re1), math.fsum(im1))
    a0s, a1s = math.fsum(abs0), math.fsum(abs1)
    k0 = a0s/abs(d0) if d0 != 0 else float('inf')
    k1 = a1s/abs(d1) if d1 != 0 else float('inf')
    return dict(d0=d0, d1=d1, terms=count, kappa0=k0, kappa1=k1,
                numerical_roundoff_heuristic=max(k0,k1)*FLOAT_EPS)


def reconstruct(d0: complex, d1: complex, s: complex, alpha: float, h: float):
    pass                                              
    if d0 == 0 or d1 == 0:
        raise ArithmeticError('D0 or D1 equals zero; estimator is undefined.')
    ratio = d1/d0
    w = cmath.log(ratio)/h
    radical = 1 + 4*alpha*alpha*w
    root = cmath.sqrt(radical)
    rho = s + (root-1)/(2*alpha*alpha)
    return dict(ratio=ratio, w=w, radical=radical, root=root, rho=rho)


def benchmark(entries, *, case, n, x, alpha, h, sigma, t):
    pass                                                                        

                                                                             
                                                                             
       
    row = dict(case=case, N=n, X=x, alpha=alpha, h=h, sigma=sigma,
               t=t, precision='float64_exploratory', reference_role='posthoc_only')
    started = time.perf_counter()
    try:
        signal = prime_differences(entries, n, x, sigma, t, alpha, h)
        result = reconstruct(signal['d0'], signal['d1'], complex(sigma, t), alpha, h)
        rho, ratio = result['rho'], result['ratio']
        data_warnings, ref_warnings = [], []
        in_strip = (0.0 < rho.real < 1.0)
        expected_right = rho.real > sigma
        if not in_strip:
            data_warnings.append('estimate_outside_critical_strip')
        if not expected_right:
            data_warnings.append('estimate_not_right_of_sampling_line')
        if max(signal['kappa0'], signal['kappa1']) >= 1000:
            data_warnings.append('large_cancellation_condition_number')
        if signal['numerical_roundoff_heuristic'] > 1e-8:
            data_warnings.append('large_float64_roundoff_heuristic')
        ratio_phase = math.atan2(ratio.imag, ratio.real)
        branch_margin = math.pi - abs(ratio_phase)
        if branch_margin < 0.15:
            data_warnings.append('ratio_near_principal_log_cut')
        radical = result['radical']
        if radical.real < 0 and abs(radical.imag) < 0.10:
            data_warnings.append('radical_near_principal_sqrt_cut')
        if abs(result['root']) < 0.25:
            data_warnings.append('quadratic_inverse_ill_conditioned')
        row.update(
            status='computed', prime_power_terms=signal['terms'],
            estimate_seconds=round(time.perf_counter()-started, 6),
            d0_re=signal['d0'].real, d0_im=signal['d0'].imag,
            d1_re=signal['d1'].real, d1_im=signal['d1'].imag,
            abs_d0=abs(signal['d0']), abs_d1=abs(signal['d1']),
            condition_kappa0=signal['kappa0'], condition_kappa1=signal['kappa1'],
            float_roundoff_heuristic=signal['numerical_roundoff_heuristic'],
            ratio_re=ratio.real, ratio_im=ratio.imag, ratio_arg=ratio_phase,
            principal_log_cut_margin=branch_margin,
            w_re=result['w'].real, w_im=result['w'].imag,
            radical_re=radical.real, radical_im=radical.imag,
            rho_re=rho.real, rho_im=rho.imag,
            estimated_in_critical_strip=in_strip,
            estimated_right_of_sampling_line=expected_right,
        )


        spectral = []
        for gam in REFERENCE_GAMMAS:
            diff = (0.5+1j*gam) - complex(sigma, t)
            spectral.append(diff + alpha*alpha*diff*diff)
        rates = [w.real for w in spectral]
        ordering = sorted(range(len(rates)), key=lambda k: rates[k], reverse=True)
        gap = rates[ordering[0]] - rates[ordering[1]]
        row.update(partial_gap=gap,
                   partial_dominant_index=ordering[0]+1,
                   partial_dominant_rate=rates[ordering[0]],
                   spectral_reference_count=len(REFERENCE_GAMMAS))
        if gap < 0.25:
            ref_warnings.append('small_finite_reference_spectral_gap')
        if case in CASE_ORDINATES:
            idx = 1 if case == 'first' else 2
            reference_rule = 'prespecified_case_reference'
        else:
            idx = min(range(len(REFERENCE_GAMMAS)),
                      key=lambda j: abs(t - REFERENCE_GAMMAS[j])) + 1
            reference_rule = 'nearest_reference_to_sampling_ordinate'
        reference = complex(0.5, REFERENCE_GAMMAS[idx-1])
        alias = abs(h * spectral[idx-1].imag) >= math.pi
        row.update(reference_re=reference.real, reference_im=reference.imag,
                   reference_index=idx, reference_selection_rule=reference_rule,
                   posthoc_error=abs(rho-reference),
                   target_q=rates[idx-1], target_omega=spectral[idx-1].imag,
                   target_is_partial_dominant=(idx-1 == ordering[0]),
                   target_alias_diagnostic=alias)
        if alias:
            ref_warnings.append('reference_zero_sampling_alias_possible')
        if not row['target_is_partial_dominant']:
            ref_warnings.append('reference_target_not_finite_list_dominant')


        row.update(numerical_warnings=';'.join(data_warnings),
                   reference_only_warnings=';'.join(ref_warnings),
                   warnings=';'.join(data_warnings+ref_warnings),
                   warning_count=len(data_warnings)+len(ref_warnings),
                   data_warning_count=len(data_warnings),
                   reference_warning_count=len(ref_warnings))
    except (ValueError, OverflowError, ArithmeticError, ZeroDivisionError) as exc:
        row.update(status='failed', warnings=f'{type(exc).__name__}: {exc}',
                   numerical_warnings=f'{type(exc).__name__}: {exc}',
                   reference_only_warnings='', warning_count=1,
                   data_warning_count=1, reference_warning_count=0,
                   estimate_seconds=round(time.perf_counter()-started, 6))
    return row


def summarize_interference(rows):
    pass                                                                        
    summary = []
    for n in sorted({int(row['N']) for row in rows}):
        rr = [row for row in rows if int(row['N']) == n]
        good = [row for row in rr if row['status'] == 'computed']
        errs = sorted(float(row['posthoc_error']) for row in good)
        def percentile(q):
            if not errs: return ''
            pos=(len(errs)-1)*q
            lo=int(math.floor(pos)); hi=int(math.ceil(pos))
            return errs[lo] + (pos-lo)*(errs[hi]-errs[lo])
        summary.append(dict(
            N=n, sample_points=len(rr), evaluated=len(good), failed=len(rr)-len(good),
            out_of_strip=sum(not row['estimated_in_critical_strip'] for row in good),
            not_right_of_sampling_line=sum(not row['estimated_right_of_sampling_line'] for row in good),
            high_cancellation=sum('large_cancellation_condition_number' in row['numerical_warnings'] for row in good),
            any_data_warnings=sum(int(row['data_warning_count'])>0 for row in good),
            small_finite_reference_gap=sum(float(row['partial_gap'])<0.25 for row in good),
            median_posthoc_error=statistics.median(errs) if errs else '',
            p90_posthoc_error=percentile(.90),
            max_posthoc_error=max(errs) if errs else '',
            fraction_posthoc_error_under_1=sum(e<1 for e in errs)/len(errs) if errs else '',
            reference_role='posthoc_only',
        ))
    return summary

def save_csv(rows, path):
    path.parent.mkdir(parents=True,exist_ok=True)
    fields=list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    hasher=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1_048_576), b''):
            hasher.update(chunk)
    return hasher.hexdigest()


def comma_ints(text):
    try:
        result = sorted({int(v.strip()) for v in text.split(',') if v.strip()})
    except ValueError as exc:
        raise argparse.ArgumentTypeError('Expected comma-separated integer cutoffs.') from exc
    if not result or min(result)<2 or max(result)>15_000_000:
        raise argparse.ArgumentTypeError('N cutoffs must be in [2, 15000000].')
    return result


def main():
    parser = argparse.ArgumentParser(description='Create exploratory CSV data; Python standard library only.')
    parser.add_argument('--profile', choices=['quick','full'], default='quick')
    parser.add_argument('--output', default='results', help='Output directory (created automatically)')
    parser.add_argument('--interference-n-values', type=comma_ints, default=None,
                        help='Comma-separated cutoffs, e.g. 100000,300000,2000000; same grid at each N')
    parser.add_argument('--interference-steps', type=int, default=None,
                        help='Number of equally spaced t values (default: 21 quick, 31 full)')
    parser.add_argument('--interference-t-min', type=float, default=13.0)
    parser.add_argument('--interference-t-max', type=float, default=23.0)
    parser.add_argument('--interference-x', type=float, default=1000.0)
    parser.add_argument('--interference-alpha', type=float, default=.30)
    parser.add_argument('--interference-h', type=float, default=.5)
    parser.add_argument('--interference-sigma', type=float, default=.25)
    args = parser.parse_args()
    ns = [30_000, 100_000, 300_000] if args.profile=='quick' else [30_000,100_000,300_000,2_000_000]
    n_inter = args.interference_n_values
    if n_inter is None:
        n_inter = [100_000,300_000] if args.profile=='quick' else [100_000,300_000,2_000_000]
    steps = args.interference_steps or (21 if args.profile=='quick' else 31)
    if steps < 2 or steps > 2001:
        parser.error('--interference-steps must be between 2 and 2001')
    if not (math.isfinite(args.interference_t_min) and math.isfinite(args.interference_t_max)
            and args.interference_t_min < args.interference_t_max):
        parser.error('Interference t-min must be smaller than t-max; both must be finite.')
    if not (args.interference_x > 1 and args.interference_alpha > 0 and args.interference_h > 0
            and 0 < args.interference_sigma < 1):
        parser.error('Interference X>1, alpha>0, h>0, sigma in (0,1) required.')
    ts = [args.interference_t_min+(args.interference_t_max-args.interference_t_min)*j/(steps-1)
          for j in range(steps)]
    out = HERE / args.output
    out.mkdir(parents=True, exist_ok=True)
    maxn = max(2_000_000, *ns, *n_inter)
    print(f'Preparing prime powers through {maxn:,}...', flush=True)
    start = time.perf_counter()
    entries = make_prime_powers(maxn)
    print(f'Generated {len(entries):,} prime-power terms in {time.perf_counter()-start:.1f}s',flush=True)
    baseline = benchmark(entries,case='first',n=2_000_000,x=1000.,alpha=.25,
                         h=.5,sigma=.25,t=14.)
    save_csv([baseline],out/'first_zero_baseline.csv')
    print('Baseline: rho_hat = {:.15f} + {:.15f}i; posthoc error = {:.3g}'.format(
        baseline['rho_re'],baseline['rho_im'],baseline['posthoc_error']),flush=True)
    rows=[]
    for case,t in CASE_ORDINATES.items():
        for n in ns:
            for x in [300.,1000.,3000.]:
                for alpha in [.20,.25,.30]:
                    rows.append(benchmark(entries,case=case,n=n,x=x,alpha=alpha,
                                          h=.5,sigma=.25,t=t))
        print(f'Finished sweep case: {case}',flush=True)
    save_csv(rows,out/'parameter_sweep.csv')
    print(f'Wrote {len(rows)} sweep rows',flush=True)
    interfer=[]
    for n in n_inter:
        for i,t in enumerate(ts):
            interfer.append(benchmark(entries,case='interference',n=n,
                           x=args.interference_x, alpha=args.interference_alpha,
                           h=args.interference_h,sigma=args.interference_sigma,t=t))
            if (i+1)%10==0 or i+1==len(ts):
                print(f'Interference N={n:,}: {i+1}/{len(ts)}',flush=True)
    save_csv(interfer,out/'interference_scan.csv')
    summary=summarize_interference(interfer)
    save_csv(summary,out/'interference_cutoff_summary.csv')
    print(f'Wrote {len(interfer)} interference rows; {len(summary)} cutoff summary rows',flush=True)
    created=[out/name for name in ('first_zero_baseline.csv','parameter_sweep.csv',
                                    'interference_scan.csv','interference_cutoff_summary.csv')]
    manifest={
        'suite_version':'3.0',
        'description':'exploratory float64 PRIME-only estimator; reference zeros POST-HOC only',
        'profile':args.profile,
        'interference_parameters':dict(N=n_inter, t_min=args.interference_t_min,
                           t_max=args.interference_t_max,steps=steps,
                           X=args.interference_x,alpha=args.interference_alpha,
                           h=args.interference_h,sigma=args.interference_sigma),
        'generated_utc':datetime.now(timezone.utc).isoformat(),
        'python':sys.version.split()[0], 'platform':platform.platform(),
        'script_sha256':sha256(Path(__file__)),
        'csv_sha256':{path.name:sha256(path) for path in created},
        'reference_gammas':list(REFERENCE_GAMMAS),
        'limitations':[
            'No sweep or interference row is an interval certificate or proof of a zero.',
            'Status computed means evaluation succeeded, NOT that the estimate is a valid zero.',
            'Numerical warnings are based on prime-only quantities; reference warnings use known zeros.',
            'Nearest reference is chosen by sampling t, not by reconstructed coordinate.',
            'A small finite-reference gap is NOT an unconditional zero dominance or isolation proof.',
            'Reference-based parameter optimization uses a known-zero oracle and is not blind recovery.',
            'A fixed cutoff N does not support extrapolation as X increases.',
        ]}
    (out/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(f'COMPLETED. CSV files in: {out}',flush=True)

if __name__=='__main__':
    main()
