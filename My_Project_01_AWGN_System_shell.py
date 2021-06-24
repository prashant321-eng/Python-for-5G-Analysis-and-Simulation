import numpy as np
import matplotlib.pyplot as plt
import numpy.random as nr
from scipy.stats import norm

blockLength = 10000000;
SNRdB = np.arange(1.0,14.1);
BER = np.zeros(SNRdB.size);
SNR = 10**(SNRdB/10);
Sym = 2*nr.randint(2,size=blockLength)-1;



noise = nr.normal(0.0,1.0,blockLength);
for K in range(SNRdB.size):
    TxBits = np.sqrt(SNR[K])*Sym;
    RxBits=TxBits+noise;
    DecBits=2*(RxBits>0)-1;
    BER[K]=BER[K]+np.sum(DecBits!=Sym);

   
BER = BER/blockLength;
plt.yscale('log')
plt.plot(SNRdB, BER,'g-');
plt.plot(SNRdB, 1-norm.cdf(np.sqrt(SNR)),'ro');
plt.grid(1,which='both')
plt.suptitle('BER for AWGN')
plt.xlabel('SNR (dB)')
plt.ylabel('BER') 
plt.savefig('BER for AWGN.png', dpi = 300)