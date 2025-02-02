
# This file was *autogenerated* from the file attack.sage
from sage.all_cmdline import *   # import sage library

_sage_const_3 = Integer(3); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_126 = Integer(126); _sage_const_5 = Integer(5); _sage_const_125 = Integer(125); _sage_const_254 = Integer(254); _sage_const_65537 = Integer(65537); _sage_const_962947420735983927056946215901134429196419130606213075415963491270 = Integer(962947420735983927056946215901134429196419130606213075415963491270); _sage_const_4 = Integer(4); _sage_const_253 = Integer(253); _sage_const_13 = Integer(13); _sage_const_11 = Integer(11); _sage_const_10 = Integer(10); _sage_const_17 = Integer(17); _sage_const_1000 = Integer(1000); _sage_const_23 = Integer(23); _sage_const_7 = Integer(7); _sage_const_19 = Integer(19); _sage_const_2305567963945518424753102147331756070 = Integer(2305567963945518424753102147331756070)#Authors: Daniel J. Bernstein and Tanja Lange
#Origin:
from sage.doctest.util import Timer
# from sage.all import var

t = Timer()
L=list()
pnum=list()
expo=list()
expo.append(_sage_const_125 )
expo.append(_sage_const_253 )
pnum.append(_sage_const_126 )
pnum.append(_sage_const_254 )
L.append(_sage_const_2305567963945518424753102147331756070 )
L.append(_sage_const_962947420735983927056946215901134429196419130606213075415963491270 )
# inpt=input('Enter no. of bits of N ')
inpt="256"
if(inpt=='256'):
  L=L[_sage_const_0 ]
  pnum=pnum[_sage_const_0 ]
  expo=expo[_sage_const_0 ]
else:
  L=L[_sage_const_1 ]
  pnum=pnum[_sage_const_1 ]
  expo=expo[_sage_const_1 ]

# L *= 701 # if 701 is included
print 'factors are', factor(L)
g = Mod(_sage_const_65537 ,L)
print 'pnum is ',expo
pmin = _sage_const_3 *_sage_const_2 **pnum
pmax = _sage_const_4 *_sage_const_2 **pnum

t.start()
x = randrange(L)
u = lift(g**randrange(L))
while True:
  p = u + randrange(ceil(pmin/L),floor(pmax/L)) * L
  if p.is_prime(): print 'p',p ;break
print 'time for first prime',t.stop().cputime
print 'u',u
t.start()
u = lift(g**randrange(L))
while True:
  q = u + randrange(ceil(pmin/L),floor(pmax/L)) * L
  if q.is_prime(): print 'q',q ;break
print 'time for second prime',t.stop().cputime

n = p * q
print 'public key',n

smooth = _sage_const_2 **_sage_const_7 *_sage_const_3 **_sage_const_3 *_sage_const_5 **_sage_const_2 *_sage_const_7 *_sage_const_11 *_sage_const_13 *_sage_const_17 *_sage_const_19 *_sage_const_23 
print 'smooth',smooth
def smoothorder(l):
  return smooth % Mod(g,l).multiplicative_order() == _sage_const_0 
v = prod(l for l,e in factor(L) if smoothorder(l))
g = Mod(g,v)
order2 = g.multiplicative_order()
print 'order2',order2
print 'order is',order
# print 'v ',v
# print 'L/v',L/v,factor(L/v)
# print 'factors are',factor(v)
# p = int("c93808a77a3e282b36a80d274d4d2d8d0497815cef52e12adedb1167aa3a2027938c6ac1e39adf74aa881934dee63e0291a93d040f06b19381f99991cc50b91075ec586f023d86a6ef4cb72853f155942b54cab979176ab290922cebe33199b1ab3bf5d035b933153f6b77914c17402026e6f668469f16c59528d9dbbc423b49",16)
u = Mod(p,v)
print 'u',u
u=int(u)
print 'new u is ',u
# u = 117186365649536878271017440711195586067867603549008793949993871455663672155382024627758506969061973692250761281716333265235764727839244690128219862454077941750864265653024596713

print 'p residue class',(p-u)/v
print 'Start loop'
print 'order2 ',order2 
for i in range(_sage_const_1000 ):
  print i
  for j in range(i*order2,(i+_sage_const_1 )*order2):
  # print i
    x = lift(g**j)
    if(x==u):
      print 'found',x
      break
print 'End loop'

t.start()

H = _sage_const_10  + _sage_const_2 **expo // v
u += floor((_sage_const_7 *_sage_const_2 **expo) // v) * v

w = lift(_sage_const_1 /Mod(v,n))

R = QQ['x']; (x,) = R._first_ngens(1)
f = (w*u+H*x)/n
g = H*x

k = _sage_const_3 
m = _sage_const_7 
print 'multiplicity',k
print 'lattice rank',m

basis = [f**j for j in range(_sage_const_0 ,k)] + [f**k*g**j for j in range(m-k)]
basis = [b*n**k for b in basis]
basis = [b.change_ring(ZZ) for b in basis]

M = matrix(m)
for i in range(m):
  M[i] = basis[i].coefficients(sparse=False) + [_sage_const_0 ]*(m-_sage_const_1 -i)
print 'time for creating matrix',t.stop().cputime

t.start()
# print '---',M
M = M.LLL()
print 'time for basis reduction',t.stop().cputime

Q = sum(z*(x/H)**i for i,z in enumerate(M[_sage_const_0 ]))
print(Q.roots())
# polyn=sum([x^2,-2*x,1])
# print 'poly roots',polyn.roots()
for r,multiplicity in Q.roots():
  print 'root is',r
  if u+v*r > _sage_const_0 :
    g = gcd(n,u+v*r)
    if g > _sage_const_1 :
      print 'successful factorization',[g,n/g]
      print 'public key',n

