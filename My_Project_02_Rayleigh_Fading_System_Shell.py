import numpy as np
import matplotlib.pyplot as plt
import numpy.random as nr

blockLength = 10000;
nBlocks = 10000;
SNRdB = np.arange(1.0,45.0,3.0);
BER = np.zeros(len(SNRdB));
BERt = np.zeros(len(SNRdB));
SNR = 10**(SNRdB/10);



for blk in range(nBlocks):
    h=(nr.normal(0.0,1.0,blockLength)+1j*nr.normal(0.0,1.0,blockLength))/np.sqrt(2);
    noise=nr.normal(0.0,1.0,blockLength)+1j*nr.normal(0.0,1.0,blockLength);
    Sym=2*nr.randint(2,size=blockLength)-1;
    for K in range(len(SNRdB)):
        TxBits = np.sqrt(SNR[K])*Sym;
        RxBits=h*TxBits+noise;
        DecBits=2*(np.real(np.conj(h)*RxBits)>0)-1;
        BER[K]=BER[K]+np.sum(DecBits!=Sym);
    
BER = BER/blockLength/nBlocks;
BERt = 1/2*(1-np.sqrt(SNR/(2+SNR)));
BERa = 1/2/SNR;
plt.yscale('log')
plt.plot(SNRdB, BER,'g-');
plt.plot(SNRdB, BERt,'ro');
plt.plot(SNRdB, BERa,'bs');
plt.grid(1,which='both')
plt.suptitle('BER for Rayleigh Fading Channel')
plt.legend(["Simulation", "Theory","Approx"], loc ="lower left");
plt.xlabel('SNR (dB)')
plt.ylabel('BER') 
plt.savefig('BER for Rayleigh Fading Channel.png', dpi = 300)