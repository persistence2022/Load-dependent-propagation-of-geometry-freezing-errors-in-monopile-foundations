                                                                                 

                                                                           
                                                                        
   
from pathlib import Path
import sys,csv,json,importlib.metadata as md
ROOT=Path(__file__).resolve().parent

import numpy as np
import openpile
from openpile.utils import py_curves,mt_curves,Hb_curves,Mb_curves

out=ROOT/'results';out.mkdir(exist_ok=True)
records=[]
for model in ('cowden','dunkirk'):
 for D,L in ((5.,10.),(10.,60.)):
  for ratio in (.05,.5,1.):
   z=ratio*L
   for dr in ((.45,.75,.90) if model=='dunkirk' else (0.,)):
    common=dict(X=z,G0=100000.,D=D,output_length=41)
    if model=='cowden':common['Su']=125.
    else:common.update(sig=100.,Dr=dr*100,L=L)
    for component,module in [('p',py_curves),('m',mt_curves),('H',Hb_curves),('M',Mb_curves)]:
     kw=common.copy()
     if model=='cowden' and component in ('H','M'):kw['L']=L
     if model=='dunkirk' and component=='m':kw['p']=2000.
     fn=getattr(module,'cowden_clay' if model=='cowden' else 'dunkirk_sand')
     x,y=fn(**kw)
     for xi,yi in zip(x,y):
      records.append(dict(model=model,D_m=D,L_m=L,z_m=z,Dr_fraction=dr,Su_Pa=125000.,sigma_Pa=100000.,G0_Pa=1e8,p_for_m_N_per_m=2e6,component=component,x_SI=float(xi),openpile_reaction_SI=float(yi)*1000))
with (out/'openpile_reference_curves.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,list(records[0]));w.writeheader();w.writerows(records)
versions={name:md.version(name) for name in ['openpile','numpy','numba','llvmlite','scipy','pandas','pydantic']}
(out/'openpile_actual_environment.json').write_text(json.dumps({'python':sys.version,'executable':sys.executable,'versions':versions,'scope':'Installed official OpenPile curve functions, Numba JIT executed; no system Python packages changed','points':len(records)},indent=2),encoding='utf-8')
print(versions);print('Official OpenPile points generated:',len(records))
