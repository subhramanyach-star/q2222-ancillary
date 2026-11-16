"""Build TABLE.md from out/verified.json (standalone verifier output) + out/vn_verified.json + BFF lower bounds."""
import json, math
V=json.load(open('out/verified.json')); R={r['id']:r for r in V}
VS=['0.999','0.99','0.98','0.95','0.9']
def cell(id):
    r=R.get(id)
    if r is None: return 'n/a'
    lo,hi=r['H_lo'],r['H_hi']
    if hi is None: return 'FAIL'
    lo=max(lo,0.0)
    if hi<5e-7: return '0 (exact)'
    if hi-lo<1e-6: return '%.6f'%((lo+hi)/2) if round(lo,6)==round(hi,6) else '[%.6f, %.6f]'%(lo,hi)
    return '[%.6f, %.6f]'%(lo,hi)
def npa(id): r=R.get(id); return '%.3f'%r['npa3_float'] if r and r['npa3_float'] is not None else '–'
def lv(id): r=R.get(id); return ','.join(sorted(set(r['levels']))) if r and r['levels'] else '–'
out=['# Rigorous single-round min-entropy intervals (bits) replacing the NPA level-3 numbers of prl_v13','',
'Each cell: exact value lies in [H_lo, H_hi]; a single number means both ends agree to 6 decimals (interval width < 10⁻⁶).',
'H_lo: exact-rational NPA certificate (level in brackets). H_hi: explicit qubit attack. Paper value (NPA-3 float) below each row.','']
def block(title,rows):
    out.append('### '+title); out.append('')
    out.append('| witness | '+' | '.join('v = '+v for v in VS)+' |'); out.append('|---|'+'---|'*len(VS))
    for name,pat in rows:
        out.append('| '+name+' | '+' | '.join(cell(pat%v) for v in VS)+' |')
        if 'conesup' not in pat: out.append('| *(paper NPA-3)* | '+' | '.join(npa(pat%v) for v in VS)+' |')
        out.append('| *(NPA level used)* | '+' | '.join(lv(pat%v) for v in VS)+' |')
    out.append('')
for th,tag in [('θ = π/8','pi_8'),('θ = 27.76° (Liu et al.)','27.76deg'),('θ = 3π/16 (ancillary)','3pi_16')]:
    block(f'{th}: H_min(A₁|X=1,E)',[('tilted CHSH g_θ',f'th{tag}_tilted_v%s'),('paper "best in cone" (its weights)',f'th{tag}_paperbest_v%s'),
        ('**sup over the whole hexagon C_θ**',f'th{tag}_conesup_v%s'),('CHSH',f'th{tag}_chsh_v%s'),('full statistics',f'th{tag}_full_v%s')])
for dl,tag in [('φ = π/6','pi_6'),('φ = π/12 (ancillary)','pi_12')]:
    for gl,gt in [('H_min(A₀B₀|00,E)','AB'),('H_min(A₀|X=0,E)','A0')]:
        block(f'WBC {dl}: {gl}',[('I_φ (centre)',f'Id{tag}_{gt}_centre_v%s'),('half-way to W₁',f'Id{tag}_{gt}_half_v%s'),('vertex W₁',f'Id{tag}_{gt}_W1_v%s'),
            ('CHSH',f'Id{tag}_{gt}_chsh_v%s'),('full statistics',f'Id{tag}_{gt}_full_v%s')])
for gl,gt in [('H_min(A₀|X=0,E)','A0'),('H_min(A₀B₀|00,E)','AB')]:
    block(f'Tsirelson point P_T: {gl}',[('CHSH (centre)',f'PT_{gt}_CHSH_v%s'),('half-way to g₀',f'PT_{gt}_half-way_v%s'),('vertex g₀',f'PT_{gt}_vertex_v%s'),('full statistics',f'PT_{gt}_full_v%s')])
# dense curves
DN=['0.999','0.995','0.99','0.985','0.98','0.97','0.96','0.95','0.94','0.92','0.9','0.88','0.86','0.84']
out.append('### Fig. 2 dense curves (certified points)'); out.append('')
out.append('| v | (a) tilted | (a) CHSH | (a) full | (b) I_φ centre | (b) CHSH | (b) full |'); out.append('|---|---|---|---|---|---|---|')
for v in DN:
    out.append('| %s | '%v+' | '.join(cell(p%v) for p in ['dense_thpi_8_tilted_v%s','dense_thpi_8_chsh_v%s','dense_thpi_8_full_v%s','dense_Idpi_6_centre_v%s','dense_Idpi_6_chsh_v%s','dense_Idpi_6_full_v%s'])+' |')
out.append('')
# von Neumann
try:
    VN={r['id']:r for r in json.load(open('out/vn_verified.json'))}; B=json.load(open('out/bff_lower.json'))
    out.append('### von Neumann H(A₁|X=1,E): BFF certified lower bound (../rates_v11) and explicit-attack upper bound (this work)'); out.append('')
    out.append('| θ | v | tilted [BFF, attack] | sup over hexagon [BFF best found, attack for all g ∈ C_θ] | full [BFF, attack] |'); out.append('|---|---|---|---|---|')
    for th,tag in [('pi/8','pi_8'),('27.76deg','27.76deg')]:
        for v in ['0.999','0.99','0.98','0.95']:
            def c(kind,bk):
                b=B.get(f'{th}|{float(v)}|{bk}'); a=VN.get(f'vn_{tag}_{kind}_v{v}')
                bs='%.4f'%max(b,0) if b is not None else '–'; as_='%.4f'%a['H_hi'] if a and a['H_hi'] is not None else '–'
                return f'[{bs}, {as_}]'
            out.append(f'| {th} | {v} | '+c('tilted','tilted')+' | '+c('cone','cone')+' | '+c('full','full')+' |')
    out.append('')
except FileNotFoundError: pass
open('TABLE.md','w').write('\n'.join(out)); print('\n'.join(out[:40]))
