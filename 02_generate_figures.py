#!/usr/bin/env python3






from __future__ import annotations

import argparse
import csv
import html
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
W, H = 600, 380
COLORS = {'axis':'#29364A','grid':'#D9E0E9','text':'#18263B',
          'first':'#2E6999','second':'#C46A40','line':'#2E6999',
          'muted':'#65748A','white':'#FFFFFF'}
CUTOFF_COLORS=['#2E6999','#C46A40','#4B8E71','#7654A3','#AF8B34','#6D7482']


def rgb(value):
    t=value.lstrip('#')
    return tuple(int(t[i:i+2],16)/255 for i in (0,2,4))


class Drawing:
    def __init__(self): self.ops=[]
    def line(self,x1,y1,x2,y2,color,width=1): self.ops.append(('line',x1,y1,x2,y2,color,width))
    def poly(self,points,color,width=2): self.ops.append(('poly',points,color,width))
    def circle(self,x,y,r,fill): self.ops.append(('circle',x,y,r,fill))
    def text(self,x,y,s,size=11,color=None,align='left',rotation=0):
        self.ops.append(('text',x,y,str(s),size,color or COLORS['text'],align,rotation))
    def rect(self,x,y,w,h,fill):self.ops.append(('rect',x,y,w,h,fill))

    def write_svg(self,path):
        lines=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img">']
        for op in self.ops:
            if op[0]=='line':
                _,x1,y1,x2,y2,c,lw=op
                lines.append(f'<line x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}" stroke="{c}" stroke-width="{lw}"/>')
            elif op[0]=='poly':
                _,pts,c,lw=op
                if pts:
                    p=' '.join(f'{x:.3f},{y:.3f}' for x,y in pts)
                    lines.append(f'<polyline points="{p}" fill="none" stroke="{c}" stroke-width="{lw}" stroke-linecap="round" stroke-linejoin="round"/>')
            elif op[0]=='circle':
                _,x,y,r,c=op
                lines.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="{r:.2f}" fill="{c}"/>')
            elif op[0]=='rect':
                _,x,y,w,h,c=op
                lines.append(f'<rect x="{x:.3f}" y="{y:.3f}" width="{w:.3f}" height="{h:.3f}" fill="{c}"/>')
            elif op[0]=='text':
                _,x,y,s,size,c,align,rotation=op
                anc={'left':'start','center':'middle','right':'end'}[align]
                transform=f' transform="rotate({rotation} {x:.3f} {y:.3f})"' if rotation else ''
                lines.append(f'<text x="{x:.3f}" y="{y:.3f}" font-family="Arial, Helvetica, sans-serif" font-size="{size}" fill="{c}" text-anchor="{anc}"{transform}>{html.escape(s)}</text>')
        lines.append('</svg>')
        path.write_text('\n'.join(lines)+'\n',encoding='utf-8')

    def write_pdf(self,path):
        def f(n):return f'{n:.4f}'
        def set_stroke(c):return ' '.join(f(i) for i in rgb(c))+' RG'
        def set_fill(c):return ' '.join(f(i) for i in rgb(c))+' rg'
        def yy(y):return H-y
        def esc(s):
            return s.encode('ascii','replace').decode('ascii').replace('\\','\\\\').replace('(','\\(').replace(')','\\)')
        out=['1 J 1 j']
        for op in self.ops:
            if op[0]=='line':
                _,x1,y1,x2,y2,c,lw=op
                out.extend([set_stroke(c),f'{f(lw)} w {f(x1)} {f(yy(y1))} m {f(x2)} {f(yy(y2))} l S'])
            elif op[0]=='poly':
                _,pts,c,lw=op
                if not pts: continue
                cmd=f'{f(pts[0][0])} {f(yy(pts[0][1]))} m '
                cmd+=' '.join(f'{f(x)} {f(yy(y))} l' for x,y in pts[1:])
                out.extend([set_stroke(c),f'{f(lw)} w '+cmd+' S'])
            elif op[0]=='rect':
                _,x,y,w,h,c=op
                out.extend([set_fill(c),f'{f(x)} {f(yy(y+h))} {f(w)} {f(h)} re f'])
            elif op[0]=='circle':
                _,x,y,r,c=op
                k=.552284749831*r

                cy=yy(y)
                cmd=(f'{f(x+r)} {f(cy)} m '
                     f'{f(x+r)} {f(cy+k)} {f(x+k)} {f(cy+r)} {f(x)} {f(cy+r)} c '
                     f'{f(x-k)} {f(cy+r)} {f(x-r)} {f(cy+k)} {f(x-r)} {f(cy)} c '
                     f'{f(x-r)} {f(cy-k)} {f(x-k)} {f(cy-r)} {f(x)} {f(cy-r)} c '
                     f'{f(x+k)} {f(cy-r)} {f(x+r)} {f(cy-k)} {f(x+r)} {f(cy)} c f')
                out.extend([set_fill(c),cmd])
            elif op[0]=='text':
                _,x,y,s,size,c,align,rotation=op

                w=sum((.29 if ch in 'il.,:;!|1 ' else .56 if ch.isalpha() else .52)*size for ch in s)
                if rotation==0:
                    if align=='center':x-=w/2
                    elif align=='right':x-=w
                    matrix=f'1 0 0 1 {f(x)} {f(yy(y))}'
                elif rotation==-90:

                    start_y=yy(y)+(w/2 if align=='center' else w if align=='right' else 0)
                    matrix=f'0 -1 1 0 {f(x)} {f(start_y)}'
                else:matrix=f'1 0 0 1 {f(x)} {f(yy(y))}'
                out.extend([set_fill(c),f'BT /F1 {f(size)} Tf {matrix} Tm ({esc(s)}) Tj ET'])
        content=('\n'.join(out)+'\n').encode('ascii')
        objects=[
            b'<< /Type /Catalog /Pages 2 0 R >>',
            b'<< /Type /Pages /Count 1 /Kids [3 0 R] >>',
            (f'<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {W} {H}] /Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>').encode(),
            b'<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>',
            (f'<< /Length {len(content)} >>\nstream\n').encode()+content+b'endstream',
        ]
        pdf=bytearray(b'%PDF-1.4\n%\xe2\xe3\xcf\xd3\n')
        offsets=[0]
        for i,obj in enumerate(objects,1):
            offsets.append(len(pdf))
            pdf.extend(f'{i} 0 obj\n'.encode()+obj+b'\nendobj\n')
        off=len(pdf)
        pdf.extend(f'xref\n0 {len(objects)+1}\n0000000000 65535 f \n'.encode())
        for pos in offsets[1:]:pdf.extend(f'{pos:010d} 00000 n \n'.encode())
        pdf.extend((f'trailer\n<< /Size {len(objects)+1} /Root 1 0 R >>\nstartxref\n{off}\n%%EOF\n').encode())
        path.write_bytes(pdf)


def read_rows(path):
    if not path.exists():
        raise FileNotFoundError(f'Missing {path}. Execute 01_generate_csv.py first.')
    with path.open(newline='',encoding='utf-8') as stream:
        return list(csv.DictReader(stream))


def good(row,key):
    try:
        value=float(row[key]);return value if math.isfinite(value) else None
    except (ValueError,TypeError,KeyError):return None


def ticks(a,b,logscale):
    if logscale:
        lo=math.floor(math.log10(a));hi=math.ceil(math.log10(b))
        vals=[(10.0**k) for k in range(lo,hi+1) if a<=10.0**k<=b]
        if len(vals)<3:
            vals=sorted(set(vals + [2*10.0**k for k in range(lo,hi+1) if a<=2*10.0**k<=b] +
                            [5*10.0**k for k in range(lo,hi+1) if a<=5*10.0**k<=b]))
        if len(vals)>7:vals=vals[::math.ceil(len(vals)/7)]
        return vals
    step=(b-a)/5
    return [a+i*step for i in range(6)]


def tick_label(v):
    if abs(v)>=10000 or (abs(v)<.001 and v !=0): return f'{v:.0e}'
    if abs(v)>=100:return f'{v:.0f}'
    if abs(v)>=10:return f'{v:.1f}'.rstrip('0').rstrip('.')
    return f'{v:.3g}'


def make_chart(title,xlab,ylab,series,*,logx=False,logy=False,floor_zero=False,note='Exploratory; not an interval certificate.'):
    validseries=[]
    for label,points,color in series:
        clean=[(float(x),float(y)) for x,y in points if x is not None and y is not None
               and math.isfinite(float(x)) and math.isfinite(float(y))
               and (not logx or x>0) and (not logy or y>0)]
        if clean:validseries.append((label,sorted(clean),color))
    if not validseries: return None
    xvals=[x for _,pp,_ in validseries for x,y in pp]
    yvals=[y for _,pp,_ in validseries for x,y in pp]
    def bounds(vals,islog):
        a,b=min(vals),max(vals)
        if islog:
            if a==b:return a/1.4,b*1.4
            pad=(math.log10(b)-math.log10(a))*.10
            return 10**(math.log10(a)-pad),10**(math.log10(b)+pad)
        if a==b:
            diff=max(abs(a)*.05,1.)
            return a-diff,b+diff
        pad=(b-a)*.08
        return a-pad,b+pad
    xmin,xmax=bounds(xvals,logx)
    ymin,ymax=bounds(yvals,logy)
    if floor_zero and not logy:
        ymin=0.0
        ymax=max(max(yvals)*1.10,0.1)
    L,R,T,B=86,566,(100 if len(validseries)>3 else 84),300
    tr=lambda v,lg:math.log10(v) if lg else v
    def px(x):return L+(tr(x,logx)-tr(xmin,logx))/(tr(xmax,logx)-tr(xmin,logx))*(R-L)
    def py(y):return B-(tr(y,logy)-tr(ymin,logy))/(tr(ymax,logy)-tr(ymin,logy))*(B-T)
    d=Drawing();d.rect(0,0,W,H,COLORS['white'])
    d.text(W/2,29,title,16,align='center')

    nseries=len(validseries)
    columns=min(3,nseries)
    col_step=490/columns
    for i,(label,points,color) in enumerate(validseries):
        x=90+(i%columns)*col_step
        y=53+(i//columns)*19
        d.line(x,y,x+18,y,color,2.2);d.circle(x+9,y,2.5,color)
        d.text(x+23,y+4,label,9.5)
    for x in ticks(xmin,xmax,logx):
        X=px(x);d.line(X,T,X,B,COLORS['grid'],.7)
        d.text(X,B+18,tick_label(x),9,COLORS['muted'],align='center')
    for y in ticks(ymin,ymax,logy):
        Y=py(y);d.line(L,Y,R,Y,COLORS['grid'],.7)
        d.text(L-9,Y+3,tick_label(y),9,COLORS['muted'],align='right')
    d.line(L,T,L,B,COLORS['axis'],1.2)
    d.line(L,B,R,B,COLORS['axis'],1.2)
    for label,pp,color in validseries:
        xy=[(px(x),py(y)) for x,y in pp]
        if len(xy)>1:d.poly(xy,color,2.4)
        for x,y in xy:d.circle(x,y,2.9,color)
    d.text((L+R)/2,339,xlab,11,align='center')
    d.text(23,(T+B)/2,ylab,11,align='center',rotation=-90)
    d.text(W/2,368,note,8.7,COLORS['muted'],align='center')
    return d


def run():
    ap=argparse.ArgumentParser(description='Write vector SVG and PDF figures; standard library only.')
    ap.add_argument('--input',default='results')
    ap.add_argument('--output',default='figures')
    a=ap.parse_args()
    indir=HERE/a.input
    out=HERE/a.output
    out.mkdir(parents=True,exist_ok=True)
    sweep=read_rows(indir/'parameter_sweep.csv')
    scan=read_rows(indir/'interference_scan.csv')
    summary=read_rows(indir/'interference_cutoff_summary.csv')
    cutoffs=sorted({int(r['N']) for r in scan})
    if not 1<=len(cutoffs)<=6:
        raise ValueError('Figure generator supports 1-6 interference cutoffs per chart.')
    maxn=max(int(r['N']) for r in sweep)
    configs=[
        ('cutoff_sensitivity','Prime-cutoff sensitivity','N (largest included prime power)',
         'Post-hoc complex error',
         lambda case:[(good(r,'N'),good(r,'posthoc_error')) for r in sweep if r['case']==case
                      and float(r['X'])==1000 and float(r['alpha'])==.25],True,True),
        ('x_sensitivity','Fixed-N sensitivity to X','X','Post-hoc complex error',
         lambda case:[(good(r,'X'),good(r,'posthoc_error')) for r in sweep if r['case']==case
                      and int(r['N'])==maxn and float(r['alpha'])==.25],True,True),
        ('alpha_sensitivity','Smoothing parameter comparison','alpha','Post-hoc complex error',
         lambda case:[(good(r,'alpha'),good(r,'posthoc_error')) for r in sweep if r['case']==case
                      and int(r['N'])==maxn and float(r['X'])==1000],False,True),
        ('runtime_sensitivity','Measured runtime versus cutoff','N (largest included prime power)',
         'Estimator time (seconds)',
         lambda case:[(good(r,'N'),good(r,'estimate_seconds')) for r in sweep if r['case']==case
                      and float(r['X'])==1000 and float(r['alpha'])==.25],True,True),
    ]
    produced=[]
    def emit(name,fig):
        if fig:
            fig.write_svg(out/(name+'.svg'))
            fig.write_pdf(out/(name+'.pdf'))
            produced.append(name)
    for name,title,xlab,ylab,getpoints,lx,ly in configs:
        series=[(case,getpoints(case),COLORS[case]) for case in ('first','second')]
        label_note=('Machine-dependent timing; record Windows/Python versions.'
                    if name=='runtime_sensitivity' else
                    'Post-hoc reference comparison; exploratory float64 arithmetic.')
        fig=make_chart(title,xlab,ylab,series,logx=lx,logy=ly,note=label_note)
        emit(name,fig)


    for name,title,ycol,ylab,logy,note in [
       ('interference_error','Interference error at matched prime cutoffs','posthoc_error',
        'Error to reference near sampling t',True,'Reference chosen by input t; no certified zero tracking.'),
       ('condition_number','Cancellation at matched prime cutoffs','condition_kappa0',
        'Sum |terms| / |D0|',True,'A cancellation indicator, not a rigorous error estimate.'),
    ]:
        series=[]
        for j,n in enumerate(cutoffs):
            data=[(good(r,'t'),good(r,ycol)) for r in scan if int(r['N'])==n]
            series.append((f'N={n:,}',data,CUTOFF_COLORS[j]))
        emit(name,make_chart(title,'Sampling ordinate t',ylab,series,
                             logx=False,logy=logy,note=note))



    nref=cutoffs[0]
    partial=[(good(r,'t'),good(r,'partial_gap')) for r in scan if int(r['N'])==nref]
    emit('partial_spectral_gap',make_chart(
        'Finite-reference-list spectral separation','Sampling ordinate t',
        'Top-two exponent-rate gap',[(f'12 reference zeros',partial,COLORS['line'])],
        logx=False,logy=False,
        note='Reference-list diagnostic only; not global spectral dominance.'))


    summary_series=[
        ('outside strip',[(good(r,'N'),float(r['out_of_strip'])/float(r['evaluated']))
                          for r in summary if float(r['evaluated'])>0],COLORS['second']),
        ('any data warning',[(good(r,'N'),float(r['any_data_warnings'])/float(r['evaluated']))
                             for r in summary if float(r['evaluated'])>0],COLORS['first']),
    ]
    emit('interference_warning_fractions',make_chart(
        'Invalid estimates and warning rates by cutoff',
        'N (largest included prime power)','Fraction of evaluated samples',
        summary_series,logx=True,logy=False,floor_zero=True,
        note='Warnings are diagnostics, not a proof of valid recovery.'))
    med=[(good(r,'N'),good(r,'median_posthoc_error')) for r in summary]
    p90=[(good(r,'N'),good(r,'p90_posthoc_error')) for r in summary]
    emit('interference_error_summary',make_chart(
        'Distribution of interference error by cutoff',
        'N (largest included prime power)','Post-hoc complex error',
        [('median',med,COLORS['first']),('90th percentile',p90,COLORS['second'])],
        logx=True,logy=True,
        note='Same sampling t grid at each cutoff; exploratory float64 error.'))

    with (out/'FIGURE_INDEX.txt').open('w',encoding='utf-8') as f:
        f.write('All figures are vector SVG + PDF, exploratory float64 only.\n')
        f.write('No chart certifies zero existence, uniqueness or the analytic error bound.\n')
        f.write('Finite known zeros used only AFTER prime-only estimator evaluation.\n')
        f.write(f'Interference cutoff curves: {cutoffs}\n\n')
        for n in produced:
            f.write(f'{n}.svg    {n}.pdf\n')
    print(f'WROTE {len(produced)} SVG and {len(produced)} PDF vector figures to {out}')

if __name__=='__main__':
    run()
