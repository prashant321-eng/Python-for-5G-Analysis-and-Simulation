import numpy as np
import numpy.random as nr
import numpy.linalg as nl
import matplotlib.pyplot as plt
import m5G

ITER = 20000;
K = 10; # number of users
Mv = np.arange(20,520,30); # number of BS antennas
Eu_dB = 10;  Eu = 10**(Eu_dB/10);
rate_MRC = np.zeros(len(Mv)) ;
bound_MRC = np.zeros(len(Mv));
rate_ZF = np.zeros(len(Mv));

beta = m5G.Dmatrix(K);
sqrtD = np.diag(np.sqrt(beta));



for it in range(ITER):
    print(ITER);
    for mx in range(len(Mv)):
        M=Mv[mx];
        #pu=Eu; #no power scaling
        pu=Eu/M; #with power scaling
        H=(nr.normal(0.0,1.0,(M,K))+1j*nr.normal(0.0,1.0,(M,K)))/np.sqrt(2);
        G=np.matmul(H,sqrtD);
        g0=G[:,0];
        MRCbf = g0/nl.norm(g0);
        nr_MRC=pu*nl.norm(g0)**2
        nr_bound_MRC=pu*M*beta[0];
        dr_bound_MRC=1;
        mu_int=np.matmul(m5G.H(MRCbf),G[:,1:]);
        dr_MRC=1+pu*nl.norm(mu_int)**2;
        dr_bound_MRC=dr_bound_MRC+pu*np.sum(beta[1:]);
        rate_MRC[mx]=rate_MRC[mx]+np.log2(1+nr_MRC/dr_MRC);
        bound_MRC[mx]=bound_MRC[mx]+np.log2(1+nr_bound_MRC/dr_bound_MRC);
        
        GG=np.matmul(m5G.H(G),G);
        nr_ZF=pu;invGG=nl.inv(GG);
        dr_ZF=np.real(invGG[0,0]);
        rate_ZF[mx]=rate_ZF[mx]+np.log2(1+nr_ZF/dr_ZF);
        
rate_MRC = rate_MRC/ITER;
bound_MRC = bound_MRC/ITER;
rate_ZF = rate_ZF/ITER;


plt.plot(Mv, rate_MRC,'g-');
plt.plot(Mv, bound_MRC,'rs');
plt.plot(Mv, rate_ZF,'b^-');
plt.grid(1,which='both')
plt.legend(["MRC", "MRC Bound", "ZF"], loc ="lower right");
plt.suptitle('RATE for MRC and ZF Receivers')
plt.ylabel('RATE ')
plt.xlabel('Number of antennas M') 
plt.savefig('RATE for MRC and ZF Receivers for perfect CSI with scaling.png', dpi = 300)
