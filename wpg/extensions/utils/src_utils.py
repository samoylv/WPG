# -*- coding: utf-8 -*
import numpy as np
wav = 6.7
per = 7.5
Ee = 3
gamma = Ee/0.511
N = 201.2/7.5




theta = 1/(gamma*np.sqrt(N))

a = -(gamma**2)*(theta**2)
b = (wav*Ee**2)/(1.306*per)
ksq = (a+b-1)*2
k = np.sqrt(ksq)
print(k)