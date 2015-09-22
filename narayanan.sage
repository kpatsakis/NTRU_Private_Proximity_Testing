#curve y**2=x**3+a*x+b
import random
import math
from time import time
import hashlib

def genprime(BITS):
	r=1
	lb=2**(BITS-1)
	ub=2*lb
	while is_prime(r)==False:
		r=random.randrange(lb,ub)
	return r

SEC=128

if SEC==128:
	# Curve25519 y^2=x^3+ax^2+x
	a=486662
	p=2^255-19
	F = Zmod(p)
	E = EllipticCurve([F(0), F(a), F(0), F(1), F(0)])
	x0=9
	y0=E.lift_x(x0)[1]
elif SEC==192:
	#M-383 	y^2 = x^3+2065150x^2+x
	#Aranha-Barreto-Pereira-Ricardini
	p = 2^383 - 187
	F = Zmod(p)
	E = EllipticCurve([F(0), F(2065150), F(0), F(1), F(0)])
	x0=12
	y0=E.lift_x(x0)[1]
else:
	#M-511 y^2 = x^3+530438x^2+x
	#Aranha-Barreto-Pereira-Ricardini
	p = 2^511 - 187
	F = Zmod(p)
	E = EllipticCurve([F(0), F(530438), F(0), F(1), F(0)])
	x0=5
	y0=E.lift_x(x0)[1]

P = E([x0,y0])
k=random.randrange(0,p)
H=k*P
tA=0

#Alice
ts=time()
la=random.randrange(0,2**32)
r=random.randrange(0,p)
C1=r*P
C2=(r+la)*H
tA= time()-ts

#Bob
ts=time()
lb=la
s=random.randrange(0,p)
t=random.randrange(0,p)
C1=s*C1+t*P
C2=s*C2+(t-s*lb)*H
tB= time()-ts

#Alice
ts=time()
res=C2+(-k)*C1
tA+= time()-ts

print tA,tB,tA+tB
