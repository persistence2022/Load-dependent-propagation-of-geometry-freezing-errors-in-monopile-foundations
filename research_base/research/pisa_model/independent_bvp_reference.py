                                                                                 
                                                                            
                                                                              
   
from pathlib import Path
import sys,csv,json
ROOT=Path(__file__).resolve().parent

import numpy as np
from scipy.integrate import solve_bvp

def curve(x,xu,k,n,yu):
    x=np.asarray(x);u=np.abs(x);s=np.sign(x)
    if n==0:return s*np.minimum(k*u,yu)
    aa=1-2*n;b=2*n*u/xu-(1-n)*(1+u*k/yu);c=(1-n)*u*k/yu-n*(u/xu)**2
    q=2*c/(-b+np.sqrt(np.maximum(0,b*b-4*aa*c)))
    return s*np.where(u>=xu,yu,q*yu)

def solve_case(L,targets):
    D=10.;t=.091;E=200e9;nu=.3;kappa=.5;Dr=.75;h=50.;ref=1e7
    EI=E*np.pi*D**3*t/8;GA=E/(2*(1+nu))*np.pi*D*t*kappa
    e=.629
    def soil(z):
        sig=10090*z;G=875*101300/(.3+.7*e**2)*np.sqrt(.6*sig/101300)
        return sig,G
    def reactions(z,w,theta):
        sig,G=soil(z);safe=np.where(sig>0,sig,1);v=w*G/(D*safe);th=theta*G/safe
        p=sig*D*curve(v,146.1-92.11*Dr,8.731-.6982*Dr-.9178*z/D,.917+.06193*Dr,.3667+25.89*Dr+(.3375-8.9*Dr)*z/L)
        mu=.2605+(-.1989+.2019*Dr)*z/L
        m=np.abs(p)*D*np.sign(th)*np.minimum(17*np.abs(th),mu)
        return p,m
    def toe(w,theta):
        sig,G=soil(L);r=L/D
        hb=sig*D**2*curve(w*G/(D*sig),.515+2.883*Dr+(.1695-.7018*Dr)*r,6.505-2.985*Dr+(-.007969-.4299*Dr)*r,.09978+.7974*Dr+(.004994-.07005*Dr)*r,.09952+.7996*Dr+(.03988-.1606*Dr)*r)
        mb=sig*D**3*curve(theta*G/sig,44.89,.3515,.3+.4986*Dr,.09981+.371*Dr+(.01998-.09041*Dr)*r)
        return hb,mb
    def ode(x,y):
        w,theta,mhat,vhat=y;p,m=reactions(x*L,w*D,theta)
        return np.vstack((L/D*(theta+ref*vhat/GA),ref*L*L*mhat/EI,(m-ref*vhat)/ref,L*p/ref))
    mesh=np.linspace(0,1,121);guess=np.zeros((4,len(mesh)));result=[];sol=None
    for target in targets:
        def bc(ya,yb):
            hb,mb=toe(yb[0]*D,yb[1]);return np.array([ya[0]-target/D,ya[2]+h/L*ya[3],yb[3]+hb/ref,yb[2]+mb/(ref*L)])
        if sol is None:
            guess[0]=target/D*(1-2*mesh);guess[1]=-2*target/L;guess[2]=.02;guess[3]=-.01
        else:mesh=sol.x;guess=sol.y*(target/previous)
        sol=solve_bvp(ode,bc,mesh,guess,tol=2e-7,bc_tol=1e-9,max_nodes=20000)
        assert sol.success,(L,target,sol.message)
        H=-sol.y[3,0]*ref;result.append(dict(case_id='C1' if L==20 else 'C4',L_m=L,head_w_m=target,H_N=H,head_theta_rad=sol.y[1,0],collocation_nodes=len(sol.x),max_rms_residual=float(max(sol.rms_residuals)),boundary_residual=float(max(abs(bc(sol.y[:,0],sol.y[:,-1]))))))
        previous=target
    return result

records=[]
for L in [20.,60.]:records+=solve_case(L,[.001,.01,.1,.5,1.])
with (ROOT/'results/independent_bvp_reference.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,list(records[0]));w.writeheader();w.writerows(records)
print(json.dumps(records,indent=2))
