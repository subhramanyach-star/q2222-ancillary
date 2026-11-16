import os
import json
VS=['0.999','0.99','0.98','0.95','0.9']
DENSE=['0.999','0.995','0.99','0.985','0.98','0.97','0.96','0.95','0.94','0.92','0.9','0.88','0.86','0.84']
D=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','rates_numerics')+'/'
from fractions import Fraction
def lamstr(l):
    q=[Fraction(x).limit_denominator(10**9) for x in l]; s=sum(q); return 'lam:'+','.join(str(x/s) for x in q)
C=[]
def add(id,kind,spec,guess,v,ref=None): C.append({'id':id,'kind':kind,'spec':spec,'guess':guess,'v':v,'npa3':ref})
for th,fn in [('pi/8','rates_th1.json'),('27.76deg','rates_th3.json'),('3pi/16','rates_th2.json')]:
    d=json.load(open(D+fn)); R=d['H_min(A1|E)']; tag=th.replace('/','_')
    for j,v in enumerate(VS):
        add(f'th{tag}_tilted_v{v}','single',{'fam':'theta','theta':th,'g':'tilted'},['A',1],v,R['tilted CHSH g_theta'][j])
        add(f'th{tag}_chsh_v{v}','single',{'fam':'theta','theta':th,'g':'chsh'},['A',1],v,R['CHSH'][j])
        add(f'th{tag}_paperbest_v{v}','single',{'fam':'theta','theta':th,'g':lamstr(d['best'][str(float(v))]['lam'])},['A',1],v,R['best in hexagon (search)'][j])
        add(f'th{tag}_full_v{v}','full',{'fam':'theta','theta':th},['A',1],v,R['full statistics'][j])
        add(f'th{tag}_conesup_v{v}','cone',{'fam':'theta','theta':th},['A',1],v,None)
for dl,fn in [('pi/6','rates_Id2.json'),('pi/12','rates_Id1.json')]:
    d=json.load(open(D+fn)); tag=dl.replace('/','_')
    for lab,guess in [('H_min(AB|00,E)',['AB',0,0]),('H_min(A0|E)',['A',0])]:
        R=d[lab]; gt='AB' if guess[0]=='AB' else 'A0'
        for j,v in enumerate(VS):
            for w,k in [('centre','I_delta (centre)'),('chsh','CHSH'),('half','half-way to W1'),('W1','vertex W1')]:
                add(f'Id{tag}_{gt}_{w}_v{v}','single',{'fam':'Id','delta':dl,'g':w},guess,v,R[k][j])
            add(f'Id{tag}_{gt}_full_v{v}','full',{'fam':'Id','delta':dl},guess,v,R['full statistics'][j])
d=json.load(open(D+'rates_PT.json'))
for lab,guess,gt in [('H_min(A0|E)',['A',0],'A0'),('H_min(AB|00,E)',['AB',0,0],'AB')]:
    R=d[lab]
    for j,v in enumerate(VS):
        for w in ['CHSH (centre)','half-way to g0','vertex g0']:
            add(f'PT_{gt}_{w.split()[0]}_v{v}','single',{'fam':'PT','g':w},guess,v,R[w][j])
        add(f'PT_{gt}_full_v{v}','full',{'fam':'PT'},guess,v,R['full statistics'][j])
dn=json.load(open(D+'dense.json'))
for j,v in enumerate(DENSE):
    for w,k in [('tilted','g_theta'),('chsh','CHSH')]:
        add(f'dense_thpi_8_{w}_v{v}','single',{'fam':'theta','theta':'pi/8','g':w},['A',1],v,dn['th1'][k][j])
    add(f'dense_thpi_8_full_v{v}','full',{'fam':'theta','theta':'pi/8'},['A',1],v,dn['th1']['full'][j])
    for w,k in [('centre','I_delta'),('chsh','CHSH')]:
        add(f'dense_Idpi_6_{w}_v{v}','single',{'fam':'Id','delta':'pi/6','g':w},['AB',0,0],v,dn['Id2'][k][j])
    add(f'dense_Idpi_6_full_v{v}','full',{'fam':'Id','delta':'pi/6'},['AB',0,0],v,dn['Id2']['full'][j])
# priority order
pri=lambda c:(0 if c['id'].startswith('thpi_8') else 1 if c['id'].startswith('th27') else 2 if c['id'].startswith('Idpi_6_AB') else 3 if c['id'].startswith('PT') else 4 if c['id'].startswith('dense') else 5)
C.sort(key=pri)
json.dump(C,open('cases.json','w'),indent=0); print(len(C))
