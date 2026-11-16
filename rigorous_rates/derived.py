"""Rigorous intervals for the derived numbers quoted in prl_v13 (ratios, percentages), from out/verified.json."""
import json
R={r['id']:r for r in json.load(open('out/verified.json'))}
def I(id): r=R[id]; return (max(r['H_lo'],0.0),r['H_hi'])
def ratio(a,b): (al,ah),(bl,bh)=I(a),I(b); return (al/bh, ah/bl)
def diff(a,b): (al,ah),(bl,bh)=I(a),I(b); return (al-bh, ah-bl)
out=[]
def add(label,paper,iv,fmt='%.4f'): out.append((label,paper,'['+fmt%iv[0]+', '+fmt%iv[1]+']'))
for v,p in [('0.999','3.7'),('0.99','2.1'),('0.98','1.06 ("6%")')]:
    add(f'pi/8 v={v}: sup_hexagon / CHSH',p,ratio(f'thpi_8_conesup_v{v}',f'thpi_8_chsh_v{v}'))
    add(f'pi/8 v={v}: paper-best / CHSH',p,ratio(f'thpi_8_paperbest_v{v}',f'thpi_8_chsh_v{v}'))
for v,p in [('0.999','1.072'),('0.99','1.142'),('0.98','1.115')]:
    add(f'pi/8 v={v}: sup_hexagon / tilted',p,ratio(f'thpi_8_conesup_v{v}',f'thpi_8_tilted_v{v}'))
    add(f'pi/8 v={v}: paper-best / tilted',p,ratio(f'thpi_8_paperbest_v{v}',f'thpi_8_tilted_v{v}'))
add('pi/8 v=0.999: sup_hexagon / full','0.997',ratio('thpi_8_conesup_v0.999','thpi_8_full_v0.999'),'%.5f')
add('pi/8 v=0.999: tilted - CHSH (bits)','0.55',diff('thpi_8_tilted_v0.999','thpi_8_chsh_v0.999'))
add('pi/8 v=0.999: sup_hexagon - tilted (bits)','0.06',diff('thpi_8_conesup_v0.999','thpi_8_tilted_v0.999'))
add('pi/8 v=0.999: sup_hexagon - CHSH (bits)','0.61',diff('thpi_8_conesup_v0.999','thpi_8_chsh_v0.999'))
add('pi/8 v=0.999: sup_hexagon - paper search (bits)','–',diff('thpi_8_conesup_v0.999','thpi_8_paperbest_v0.999'),'%.2e')
for v,p in [('0.999','2.5'),('0.99','1.9'),('0.98','1.5')]:
    add(f'27.76deg v={v}: sup_hexagon / CHSH',p,ratio(f'th27.76deg_conesup_v{v}',f'th27.76deg_chsh_v{v}'))
add('27.76deg v=0.999: tilted / CHSH','2.4',ratio('th27.76deg_tilted_v0.999','th27.76deg_chsh_v0.999'))
for v in ['0.999','0.99']:
    add(f'WBC pi/6 AB v={v}: centre - half-way (bits, >0 required by Cor.1(b))','>0',diff(f'Idpi_6_AB_centre_v{v}',f'Idpi_6_AB_half_v{v}'))
for v in ['0.999','0.99','0.98','0.95','0.9']:
    add(f'P_T A0 v={v}: CHSH - half-way (bits, >0 by Cor. chshopt)','>0',diff(f'PT_A0_CHSH_v{v}',f'PT_A0_half-way_v{v}'))
s=['| quantity | paper | rigorous interval |','|---|---|---|']+['| %s | %s | %s |'%t for t in out]
open('DERIVED.md','w').write('# Derived numbers of prl_v13, rigorous\n\n'+'\n'.join(s)+'\n'); print('\n'.join(s))
