import numpy as np

def pderiv(ar,dx=1.,ax=0,order=2,smth=None):
   """
      pderiv gives the first partial derivative
      of a periodic array along a given axis.

      Inputs:
         ar - The input array
         dx - Grid spacing, defaults to 1.
         ax - Axis along which to take the derivative
         order - Order of accuracy, (1,2) defaults to 2

      Output:
         dar - The derivative array
   """
	if smth is not None:
		ar = gf(ar,sigma=smth)
	if order == 1:
		dar = (np.roll(ar,-1,axis=ax)-ar)/dx
	elif order == 2:
		dar = (np.roll(ar,-1,axis=ax)-np.roll(ar,1,axis=ax))/(2*dx)

	return dar

def divergenced(Yx,Yy,Yz,dx):
	divyx=pderiv(Yx,dx,ax=0,order=1)
	divyy=pderiv(Yy,dx,ax=1,order=1)
	divyz=pderiv(Yy,dx,ax=2,order=1)
	divy=divyx+divyy+divyz
	lyxx,lyyy,lyzz=np.shape(divy)
	divy_1=np.delete(divy, np.s_[lyxx-1], axis=0)
	divy_2=np.delete(divy_1, np.s_[lyyy-1], axis=1)
	divy_3=np.delete(divy_2, np.s_[lyzz-1], axis=2)
	return divy_3


def omnisolidangle(dsdt,dy,dh,lagzz,dx):
##THIS ASSUMES THE LAG ALONG X,Y, AND Z ARE SAME ARRAYS ##
	lent=len(lagzz)
	dSdT=np.zeros(lent)
	divY=np.zeros(lent)
	divH=np.zeros(lent)
	lag=np.zeros(lent)
	size=dsdt.shape[0]
	Z,Y,X=np.ogrid[:size,:size,:size]
	dis=np.sqrt(X**2+Y**2+Z**2)

### CONSTRUCTING THE AVERAGE OVER THE ANNULUS ###
	for k in range(lent):
		a=lagzz[k]
		mask=(((a-0.5)<=dis)&(dis<=a+0.5))
		lagx,lagy,lagz=mask.nonzero()
		dSp=0.
		divYp=0.
		divHp=0.
		for j in range(len(lagx)):
			dSp=dSp+dsdt[lagx[j],lagy[j],lagz[j]]
			divYp=divYp+dy[lagx[j],lagy[j],lagz[j]]
			divHp=divHp+dh[lagx[j],lagy[j],lagz[j]]

		dSdT[k]=dSp/(len(lagx))
		divY[k]=divYp/(len(lagx))
		divH[k]=divHp/(len(lagx))
		lag[k]=a*dx
	return lag,dSdT,divY,divH

