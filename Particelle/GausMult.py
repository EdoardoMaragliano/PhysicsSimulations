import numpy as np 
import matplotlib.pyplot as plt
from   ROOT import *

d   = 3
m   = np.array([1,2,3]).reshape(3,1)
cov = np.array([ [2   , 1.1,   2],
                 [1.1 , 3  , 0.3],
                 [2   , 0.3,   4] ])

L   = np.linalg.cholesky(cov)

h2 = TH2D("h2","",100,-10,10,100,-10,10)
for i in range(0,100000):
    u = np.random.normal(loc=0, scale=1, size=d).reshape(3,1)
    x = m + np.dot(L,u)
    h2.Fill(x[0],x[2])

h2.Draw()
print(h2.GetCorrelationFactor())
gApplication.Run(True)
