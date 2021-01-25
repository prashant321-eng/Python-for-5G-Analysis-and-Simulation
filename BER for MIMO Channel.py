# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 20:48:51 2020

@author: pc
"""
import numpy as np
import matplotlib.pyplot as plt
import numpy.random as nr
import numpy.linalg as nl

blockLength = 1000;
nBlocks = 1000;
Nr = 2;
Nt = 2;
SNRdB = np.arange(1.0,25.0,1.0);
BER = np.zeros(len(SNRdB));
BERt = np.zeros(len(SNRdB));
SNR = 10**(SNRdB/10);

for blk in range(nBlocks):
    H=(nr.normal(0.0,1.0,(Nr,Nt))+1j*nr.normal(0.0,1.0,(Nr,Nt)))/np.sqrt(2);
    noise=nr.normal(0.0,1.0,(Nr,blockLength))+1j*nr.normal(0.0,1.0,(Nr,blockLength));
    Sym=2*nr.randint(2,size=(Nt,blockLength))-1;
    for K in range(len(SNRdB)):
        TxBits = np.sqrt(SNR[K])*Sym;
        RxBits=np.matmul(H,TxBits)+noise;
        ZFout=np.matmul(nl.pinv(H),RxBits);
        DecBits=2*(np.real(ZFout)>0)-1;
        BER[K]=BER[K]+np.sum(DecBits!=Sym);
    
BER = BER/blockLength/nBlocks/Nt;
BERt = 1/2*(1-np.sqrt(SNR/(2+SNR)));
plt.yscale('log')
plt.plot(SNRdB, BER,'g-');
plt.plot(SNRdB, BERt,'ro');
plt.grid(1,which='both')
plt.suptitle('BER for MIMO Channel')
plt.legend(["Simulation", "Theory"], loc ="lower left");
plt.xlabel('SNR (dB)')
plt.ylabel('BER') 



blockLength = 1000;
nBlocks = 1000;
Nr = 3;
Nt = 2;
SNRdB = np.arange(1.0,25.0,1.0);
BER = np.zeros(len(SNRdB));
BERt = np.zeros(len(SNRdB));
SNR = 10**(SNRdB/10);

for blk in range(nBlocks):
    H=(nr.normal(0.0,1.0,(Nr,Nt))+1j*nr.normal(0.0,1.0,(Nr,Nt)))/np.sqrt(2);
    noise=nr.normal(0.0,1.0,(Nr,blockLength))+1j*nr.normal(0.0,1.0,(Nr,blockLength));
    Sym=2*nr.randint(2,size=(Nt,blockLength))-1;
    for K in range(len(SNRdB)):
        TxBits = np.sqrt(SNR[K])*Sym;
        RxBits=np.matmul(H,TxBits)+noise;
        ZFout=np.matmul(nl.pinv(H),RxBits);
        DecBits=2*(np.real(ZFout)>0)-1;
        BER[K]=BER[K]+np.sum(DecBits!=Sym);
    
BER = BER/blockLength/nBlocks/Nt;
BERt = 1/2*(1-np.sqrt(SNR/(2+SNR)));
plt.yscale('log')
plt.plot(SNRdB, BER,'g-');
plt.plot(SNRdB, BERt,'ro');
plt.grid(1,which='both')
plt.suptitle('BER for MIMO Channel')
plt.legend(["Simulation", "Theory"], loc ="lower left");
plt.xlabel('SNR (dB)')
plt.ylabel('BER') 

blockLength = 1000;
nBlocks = 1000;
Nr = 4;
Nt = 4;
SNRdB = np.arange(1.0,25.0,1.0);
BER = np.zeros(len(SNRdB));
BERt = np.zeros(len(SNRdB));
SNR = 10**(SNRdB/10);

for blk in range(nBlocks):
    H=(nr.normal(0.0,1.0,(Nr,Nt))+1j*nr.normal(0.0,1.0,(Nr,Nt)))/np.sqrt(2);
    noise=nr.normal(0.0,1.0,(Nr,blockLength))+1j*nr.normal(0.0,1.0,(Nr,blockLength));
    Sym=2*nr.randint(2,size=(Nt,blockLength))-1;
    for K in range(len(SNRdB)):
        TxBits = np.sqrt(SNR[K])*Sym;
        RxBits=np.matmul(H,TxBits)+noise;
        ZFout=np.matmul(nl.pinv(H),RxBits);
        DecBits=2*(np.real(ZFout)>0)-1;
        BER[K]=BER[K]+np.sum(DecBits!=Sym);
    
BER = BER/blockLength/nBlocks/Nt;
BERt = 1/2*(1-np.sqrt(SNR/(2+SNR)));
plt.yscale('log')
plt.plot(SNRdB, BER,'g-');
plt.plot(SNRdB, BERt,'ro');
plt.grid(1,which='both')
plt.suptitle('BER for MIMO Channel')
plt.legend(["Simulation", "Theory"], loc ="lower left");
plt.xlabel('SNR (dB)')
plt.ylabel('BER') 

blockLength = 1000;
nBlocks = 1000;
Nr = 2;
Nt = 1;
SNRdB = np.arange(1.0,25.0,1.0);
BER = np.zeros(len(SNRdB));
BERt = np.zeros(len(SNRdB));
SNR = 10**(SNRdB/10);

for blk in range(nBlocks):
    H=(nr.normal(0.0,1.0,(Nr,Nt))+1j*nr.normal(0.0,1.0,(Nr,Nt)))/np.sqrt(2);
    noise=nr.normal(0.0,1.0,(Nr,blockLength))+1j*nr.normal(0.0,1.0,(Nr,blockLength));
    Sym=2*nr.randint(2,size=(Nt,blockLength))-1;
    for K in range(len(SNRdB)):
        TxBits = np.sqrt(SNR[K])*Sym;
        RxBits=np.matmul(H,TxBits)+noise;
        ZFout=np.matmul(nl.pinv(H),RxBits);
        DecBits=2*(np.real(ZFout)>0)-1;
        BER[K]=BER[K]+np.sum(DecBits!=Sym);