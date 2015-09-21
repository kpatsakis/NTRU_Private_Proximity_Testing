from random import randint, shuffle
from time import time

R.<x> = QQ['x'];

#set the parameters
p=3
q=2048
SEC=128
if SEC==128:
	n=439
	D1=9;D2=8;D3=5
	Dg=146;Dm=112
elif SEC==192:
	n=593
	D1=10;D2=10;D3=8
	Dg=197;Dm=158
else:
	n=743
	D1=11;D2=11;D3=15
	Dg=247;Dm=204

#returns a polynomial of degree NN with coefficients in [-B,B]
def sample(NN, B):
    s=[randint(-B, B) for i in range(NN)]
    poly=R(s)
    return poly

def sample2(NN, B,o,mo):
    oo=[1 for i in range(o)]
    mmoo=[-1 for i in range(mo)]
    zz=[0 for i in range(NN-o-mo)]
    s=oo+mmoo+zz
    shuffle(s)
    poly=R(s)
    return poly

def pn(NN,d1,d2,d3):
	a=randint(2*d1, NN-d2)
	p1=sample2(a, 1, d1, d1)
	p2=sample2(NN-a, 1, d2, d2)
	p3=sample2(NN, 1, d3, d3)
	poly=p1*p2+p3
	return poly

#reduces the coefficients of polynomial f modulo pp
#and "fixes" them so that they belong to [-pp/2,pp/2]
def modCoeffs(f, pp):
    clist=f.list()
    p2=pp/2
    for i in range(len(clist)):
        clist[i] = clist[i]%pp
        if clist[i]>p2:
            clist[i]-=pp
    return R(clist)

def inv_poly_mod2(poly,NNN):
	k=0;b=1;c=0;
	f=poly;g=x^NNN-1
	f=modCoeffs(f, 2)
	res=False
	while True:
		while f(0)==0 and not f.is_zero():
			f=f.shift(-1)
			c=c*x
			c=modCoeffs(c, 2)
			k+=1
		if f.is_one():
			e=NNN-k
			while e<0:
				e+=NNN
			retval= x^e*b
			res=True
			break
		elif f.degree()==-1 or f.is_zero() or f==0:
			break
		if f.degree()<g.degree():
			f,g=g,f
			b,c=c,b
		f=f+g
		b=b+c
		f=modCoeffs(f, 2)
		c=modCoeffs(c, 2)
	if res:
		retval=retval%(x^NNN-1)
		retval=modCoeffs(retval, 2)
		return True, retval
	else:
		return False,0

def inv_poly_mod3(poly,NNN):
	k=0;b=1;c=0;
	f=poly;g=x^NNN-1
	res=False
	while True:
		while f(0)==0 and not f.is_zero():
			f=f.shift(-1)
			c=c*x
			c=modCoeffs(c, 3)
			k+=1
		if f.is_one():
			e=NNN-k
			while e<0:
				e+=NNN
			retval= x^e*b
			res=True
			break
		elif (-f).is_one():
			e=NNN-k
			while e<0:
				e+=NNN
			retval= -x^e*b
			res=True
			break
		elif f.degree()==-1 or f.is_zero() or f==0:
			break
		if f.degree()<g.degree():
			f,g=g,f
			b,c=c,b
		if f(0)==g(0):
			f=f-g
			b=b-c
		else:
			f=f+g
			b=b+c
		f=modCoeffs(f, 3)
		c=modCoeffs(c, 3)
	if res:
		retval=retval%(x^NNN-1)
		retval=modCoeffs(retval, 3)
		return True, retval
	else:
		return False,0

def inv_poly_mod_prime_pow(poly,r,NNN):
	res,b=inv_poly_mod2(poly,NNN)
	if res:
		qr=2
		while qr<2^r:
			qr=qr^2
			b=b*(2-poly*b)
			b=b%(x^NNN-1)
			b=modCoeffs(b, 32)
		return True,b
	else:
		return False,0

def gen_priv_key(deg):
	res=False
	pw=log(q,2)
	while (res==False):
		poly=pn(deg,D1,D2,D3)
		poly=1+3*poly
		ppInv=inv_poly_mod3(poly, deg)[1]
		res,pqInv=inv_poly_mod_prime_pow(poly,pw, deg)
	return poly,ppInv,pqInv

pw=log(q,2)
#Generate the key pair
f,fp,fq=gen_priv_key(n)
res=False
while res==False:
	g=sample2(n, 1, Dg+1, Dg)
	res,gInv=inv_poly_mod_prime_pow(g,pw, n)
h=p*g*fq
h=h%(x^n-1)
h=modCoeffs(h,q)
hinv=gInv*f
hinv=hinv%(x^n-1)
hinv=modCoeffs(hinv,q)

tA=0

#alice
res=False
ts=time()
#generate an invertible "noise" polynomial
while res==False:
	s=sample(n-1, q)
	s=q*s+1
	res,sInv=inv_poly_mod_prime_pow(s,pw, n)
c=s*h+m
c=c%(x^n-1)
c= modCoeffs(c,q)
tA=time()-ts

#bob
ts=time()
r=sample2(n, 1, Dm, Dm)
c=r*(c-m)
c=c%(x^n-1)
tB=time()-ts

#alice
ts=time()
c=c*sInv*gInv*f
c=c%(x^n-1)
c=modCoeffs(c,q)
tA+=time()-ts

print "Alice:",tA,"Bob",tB, "\nTotal:",tA+tB
