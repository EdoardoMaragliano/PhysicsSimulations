import numpy as np
import matplotlib.pyplot as plt
import random as rnd
from   math import *

N0   = 10000
tau1 = 3125
tau2 = 184.5
tau3 = 0.746


dt   = 100
Nt   = 100

##per ogni atomo calcolo la sequenza di decadimento
t1   = np.zeros(N0)     ##istante in cui l'atomo j-esimo decade Pb211->Bi211
t2   = np.zeros(N0)     ##tempo in cui decade Bi211->Po211
t3   = np.zeros(N0)     ##tempo in cui decade Po211->Pb207


for i in range(0,N0):
    t1[i]=-tau1*log(rnd.random());
    t2[i]=t1[i]-tau2*log(rnd.random());
    t3[i]=t2[i]-tau3*log(rnd.random());

n1   = np.zeros(Nt+1)   ##Numero nuclei Pb211 al variare del tempo
n2   = np.zeros(Nt+1)
n3   = np.zeros(Nt+1)
n4   = np.zeros(Nt+1)

# Decay

                     
t    = np.linspace(0,dt*Nt,Nt+1)        ##divido l'intervallo temporale in  Nt segmenti

n1[0]=N0        ##parto con N0 nuclei di piombo
n2[0]=0
n3[0]=0
n4[0]=0


for i in range(1,Nt+1):     ##faccio avanzare il contatore temporale

    tn1 = 0
    tn2 = 0
    tn3 = 0
    tn4 = 0

                                        ##per ogni campionamento temporale i
    for j in range(1,N0):               ##cont    
        
        if(t[i]<t1[j]):                  ##se il tempo è inferiore al tempo di decadimento 1 del nucleo j
            tn1=tn1+1                    ##il nucleo j è di tipo 1   
        elif(t[i]<t2[j]):         
            tn2=tn2+1
        elif(t[i]<+t3[j]):
            tn3=tn3+1
        elif(t[i]>t3[j]):    
            tn4=tn4+1
    n1[i]=tn1
    n2[i]=tn2
    n3[i]=tn3
    n4[i]=tn4

plt.plot(t,n1,'k', label='211Pb')
plt.plot(t,n2,'c', label='211Bi')
plt.plot(t,n3,'g', label='211Po')
plt.plot(t,n4,'r', label='207Pb')
plt.legend()
plt.title('211Pb Decay')
plt.show()
