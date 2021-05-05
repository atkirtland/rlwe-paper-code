import numpy as np

def a(g,pm):
    return 8*(g/2)-pm*kronecker(g/2,3)
def b(g,pm):
    return 4*(g/2)-pm*kronecker(g/2,3)

primes = [i for i in range(5, 500) if Integer(i).is_prime()]
gaps = [primes[i+1] - primes[i] for i in range(len(primes)-1)]
plusAGaps = [a(g,1) for g in gaps]
minusAGaps = [a(g,-1) for g in gaps]
plusAGapsPartials = np.cumsum(plusAGaps)
plusAGapsPartials = np.insert(plusAGapsPartials,0,0)
minusAGapsPartials = np.cumsum(minusAGaps)
minusAGapsPartials = np.insert(minusAGapsPartials,0,0)
plus2Coeffs = -(22+1+plusAGapsPartials)
minus2Coeffs = -(22-1+minusAGapsPartials)

plusBGaps = [b(g,1) for g in gaps]
minusBGaps = [b(g,-1) for g in gaps]
plusBGapsPartials = np.cumsum(plusBGaps)
plusBGapsPartials = np.insert(plusBGapsPartials,0,0)
minusBGapsPartials = np.cumsum(minusBGaps)
minusBGapsPartials = np.insert(minusBGapsPartials,0,0)
plus1Coeffs = (13+2+plusBGapsPartials)
minus1Coeffs = (13-2+minusBGapsPartials)

alphas = [i for i in range(500) if Integer(i+4).is_prime()]

R.<x> = ZZ[]
polys = [((x-3*primes[i])^alphas[i] 
          *(x-primes[i])^alphas[i]
          *(x^3+plus2Coeffs[i]*x^2+2*primes[i]*plus1Coeffs[i]*x-3^2*primes[i]^2)
          *(x^3+minus2Coeffs[i]*x^2+2*primes[i]*minus1Coeffs[i]*x-3*primes[i]^2)
         )
        for i in range(len(primes))]
        
A = [3*i for i in range(5, 500) if Integer(i).is_prime()]
polys2 = []
for z in A:
   z = Integer(z)
   factors = z.factor()
   s = len(factors)
   v = [(-1)^s]*euler_phi(z)
   for pair in factors:
       factor = pair[0]
       mult = pair[1]
       for i in range(euler_phi(z)):
           if factor.divides(i):
               v[i] *= -euler_phi(factor)
   mhm = matrix.toeplitz(v, v[1:], ZZ)
   polys2.append(mhm.charpoly())

for i in range(len(polys2)):
    print(polys[i] == polys2[i])
