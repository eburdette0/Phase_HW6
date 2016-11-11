# -*- coding: utf-8 -*-
"""
Geol 2460 HW6
11-14-16
@author: Eric Burdette

"""

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

Tcrit=1457+273 #critical temperature from diagram is ~1457C
Omega=2*8.314*Tcrit #regular solution parameter can be solved for and used for rest of predicitons
Tlist=np.array([1450,1440,1425,1400,1300,1350,1200,1100,1000,800,600])+273


def Gmix(Xb,T,omega): # regular solution delta G of mixing
    return 8.314*T*((1-Xb)*np.log(1-Xb)+Xb*np.log(Xb))+omega*Xb*(1-Xb)
    
def spinodal(omega,T): # spinodal points are roots of second derivative
    if T<=Tcrit:
        #polynomial is -R*T+2*omega*Xb-2*omega*x^2=0
        roots=np.roots([1,-1,8.314*T/(2*omega)])
        return roots[0],roots[1] 
    else:
        return None

def Gmix_deriv(Xb,T,omega): # First derivative=0 colse to either end-member is solvus
    return 8.314*T*np.log(Xb/(1-Xb))+omega*(1-2*Xb)



#plot mixing free energies
XbList=np.linspace(0,1,100)
plt.plot(XbList,Gmix(XbList,1773,Omega), label='1773')
for T in Tlist:
    GmixT=Gmix(XbList,T,Omega)
    plt.plot(XbList,GmixT, label=str(T))
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.title('del_Gmix for the En-Di system'); plt.xlabel('X_Di'); plt.ylabel('del_Gmix (J)')
plt.show()

#plot and output solvus/spinodal points for a given temperature
Tspin=[]
xlist_spin=[]
xlist_solv=[]
print('T (K), left sovlus, left spinodal, right spinodal, right solvus')
for T in Tlist:
    sp2,sp1=spinodal(Omega,T)
    solv1=fsolve(Gmix_deriv,x0=0+1e-7,args=(T,Omega))[0]
    solv2=fsolve(Gmix_deriv,x0=1-1e-7,args=(T,Omega))[0]
    Tspin.append(T);xlist_spin.append(sp1); xlist_solv.append(solv1)
    Tspin.append(T);xlist_spin.append(sp2); xlist_solv.append(solv2)
    print('{:04.0f}'.format(T),','+4*' ','{:01.4f}'.format(solv1),','+6*' ','{:01.4f}'.format(sp1),','+7*' ','{:01.4f}'.format(sp2),','+6*' ','{:01.4f}'.format(solv2))
Tspin.append(Tcrit);xlist_spin.append(0.5); xlist_solv.append(0.5)
plt.plot(xlist_spin,Tspin,'x')
plt.plot(xlist_solv,Tspin, 'o')
plt.title('Phase Boundaries for the En-Di system'); plt.xlabel('X_Di'); plt.ylabel('T (K)')

#comments
print(' ')
print('C) Pyronexe compositions are the solvus -> d(del_Gmix)/dXb=0')
print('   Results are as calculated above -> 71%/29%, and 29%/71%')