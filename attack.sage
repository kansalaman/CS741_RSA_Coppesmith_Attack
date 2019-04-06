#Authors: Daniel J. Bernstein and Tanja Lange
#Origin:
from sage.doctest.util import Timer
# from sage.all import var

t = Timer()
L=list()
pnum=list()
expo=list()
expo.append(125)
expo.append(253)
pnum.append(126)
pnum.append(254)
L.append(2305567963945518424753102147331756070)
L.append(962947420735983927056946215901134429196419130606213075415963491270)
# inpt=input('Enter no. of bits of N ')
inpt="256"
if(inpt=='256'):
  L=L[0]
  pnum=pnum[0]
  expo=expo[0]
else:
  L=L[1]
  pnum=pnum[1]
  expo=expo[1]

# L *= 701 # if 701 is included
print 'factors are', factor(L)
g = Mod(65537,L)

pmin = 3*2**pnum
pmax = 4*2**pnum

t.start()
x = randrange(L)
u = lift(g^randrange(L))
while True:
  p = u + randrange(ceil(pmin/L),floor(pmax/L)) * L
  if p.is_prime(): print 'p',p ;break
print 'time for first prime',t.stop().cputime
print 'u',u
t.start()
u = lift(g^randrange(L))
while True:
  q = u + randrange(ceil(pmin/L),floor(pmax/L)) * L
  if q.is_prime(): print 'q',q ;break
print 'time for second prime',t.stop().cputime

n = p * q
print 'public key',n

smooth = 2^7*3^3*5^2*7*11*13*17*19*23
print 'smooth',smooth
def smoothorder(l):
  return smooth % Mod(g,l).multiplicative_order() == 0
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
for i in range(1000):
  print i
  for j in range(i*order2,(i+1)*order2):
  # print i
    x = lift(g^j)
    if(x==u):
      print 'found',x
      break
print 'End loop'

t.start()

H = 10 + 2**expo // v
u += floor((7*2**expo) // v) * v

w = lift(1/Mod(v,n))

R.<x> = QQ[]
f = (w*u+H*x)/n
g = H*x

k = 3
m = 7
print 'multiplicity',k
print 'lattice rank',m

basis = [f^j for j in range(0,k)] + [f^k*g^j for j in range(m-k)]
basis = [b*n^k for b in basis]
basis = [b.change_ring(ZZ) for b in basis]

M = matrix(m)
for i in range(m):
  M[i] = basis[i].coefficients(sparse=False) + [0]*(m-1-i)
print 'time for creating matrix',t.stop().cputime

t.start()
# print '---',M
M = M.LLL()
print 'time for basis reduction',t.stop().cputime

Q = sum(z*(x/H)^i for i,z in enumerate(M[0]))
print(Q.roots())
# polyn=sum([x^2,-2*x,1])
# print 'poly roots',polyn.roots()
for r,multiplicity in Q.roots():
  print 'root is',r
  if u+v*r > 0:
    g = gcd(n,u+v*r)
    if g > 1:
      print 'successful factorization',[g,n/g]
      print 'public key',n