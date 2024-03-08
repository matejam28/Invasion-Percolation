import numpy as np
from gstools import SRF, Gaussian
import os



def porToPerm(por):
	log10Perm=por*10
	return 10**log10Perm

mu = 0.20
variance = 0.001
nn=100

x = y = z = range(nn)
model = Gaussian(dim=3, var=variance, len_scale=10)
srf = SRF(model, seed=1)
field = srf.structured([x, y, z])
field=field+mu
poroFile=open('PORO.INC','w')
poroFile.write('PORO \n')
permFile=open('PERM.INC','w')
permFile.write('PERMX \n')
for por in field:
	for poro in por:
		for poros in poro:
			poroFile.write(str(poro)+'\n')
			permFile.write(str(porToPerm(poro))+'\n')
poroFile.write('/ \n')
permFile.write('/ \n')
poroFile.close()
permFile.close()