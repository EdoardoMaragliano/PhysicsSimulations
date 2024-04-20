import numpy as np 
import matplotlib.pyplot as plt
from   ROOT import *

d   = 3
m   = np.array([1,2,3])
Cov = np.array([ [2   , 1.1,   2],
                 [1.1 , 3  , 0.3],
                 [2   , 0.3,   4] ])

h2 = TH2D("h2","",100,-10,10,100,-10,10)
for i in range(0,10000):
    x = np.random.multivariate_normal(mean=m,cov=Cov,size=3)
    h2.Fill(x[0][0],x[0][2])

h2.Draw()
print(h2.GetCorrelationFactor())
gApplication.Run(True)
