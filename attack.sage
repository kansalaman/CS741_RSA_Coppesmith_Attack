#Authors: Daniel J. Bernstein and Tanja Lange
#Origin:
from sage.doctest.util import Timer
# from sage.all import var

t = Timer()

L = 962947420735983927056946215901134429196419130606213075415963491270
# L *= 701 # if 701 is included
print 'factors are', factor(L)
g = Mod(65537,L)

pmin = 3*2**254
pmax = 4*2**254

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
# print 'kdshfdglkglgtkl',log(p*q)
# n=19920726839070815994998091093341971896358011164873513502252068165898094300547540292355376641783342604056684238872899431321736308213487545942014350950815359041644853865113261438495021685898037777675661465630271412797099783047138041003690144583180000881213443421779202855681661055259092513193388583127872981726859546745771237224182179841241761252421057745479598947549253560288134249980235300420902573907076048126110226066015101444470870185069166904393225123156478224554196871078974124493506357070691906706621812669657437944813699547511830671563781996836773347584294392459515072831363041501254474779301029958385193495057
# n = int("13965a4642b9d88f507ac4223792f8edd0670f302a490523e2f2b616f7a284e849d5fc55193aa76dc4499ac194b410ebdd0bf82448120e062f69769c6384d3d701b60b3ab06adc1bb2f296ad726aa0ce7ba03b2288c02884cba3d341b8c8e8600f492b52fc7d500567ea4e50f3fa16559d7b3cf0db31ab0c751d233b61e5d2c5ab8cc4c4f6f99e3ae5c62196921821b4e6f1ff1fea0aa8015743070015707a2f795ce4d203b423c9f93f65527b5eb6d4f0c8ffb61d747382c2e56969fe23eaaf471ef79e6bc272bb0aa2ff39f649ab6327e627f4dca5b2a4c75269d561d64e039108d25aa35e6638048f488344f74d0799490096e9f47014f7fa9c15717d5b11f",16)
print 'public key',n

smooth = 2^7*3^3*5^2*7*11*13*17*19*23
print 'smooth',smooth
def smoothorder(l):
  return smooth % Mod(g,l).multiplicative_order() == 0
# print 'order1',g.multiplicative_order()
# print 'g is',g.multiplicative_order(), lift(g^g.multiplicative_order())
v = prod(l for l,e in factor(L) if smoothorder(l))
g = Mod(g,v)
order2 = g.multiplicative_order()
print 'order2',order2
order2 = order2/10000
# print 'v ',v
# print 'L/v',L/v,factor(L/v)
# print 'factors are',factor(v)
# p = int("c93808a77a3e282b36a80d274d4d2d8d0497815cef52e12adedb1167aa3a2027938c6ac1e39adf74aa881934dee63e0291a93d040f06b19381f99991cc50b91075ec586f023d86a6ef4cb72853f155942b54cab979176ab290922cebe33199b1ab3bf5d035b933153f6b77914c17402026e6f668469f16c59528d9dbbc423b49",16)
u = Mod(p,v)
print 'u',u
u=int(u)
# u = 117186365649536878271017440711195586067867603549008793949993871455663672155382024627758506969061973692250761281716333265235764727839244690128219862454077941750864265653024596713

# u = 2
# u=
# print 'p residue class',(p-u)/v
# print 'Start loop'
# for i in range(10000):
#   print i
#   for j in range(order2):
#   # print i
#     x = lift(g^i)
#     if(x==u):
#       print 'found',x
#       break
# print 'End loop'

t.start()

H = 10 + 2**253 // v
u += floor((7*2**253) // v) * v

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
    if g > 1: print 'successful factorization',[g,n/g]