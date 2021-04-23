from networkx import read_graph6
from itertools import product, permutations, combinations
from math import gcd, isqrt
from sympy import factorint, isprime

def stvorec(n):
        if n < 0: return False
        s = isqrt(n)
        return s*s == n

def dopln(S,T):
	if not square(2*T - S*S): return set()
	res = isqrt(2*T - S*S)
	if S%2 != res%2: return set()
	return {[(S + res)//2, (S - res)//2],
                [(S - res)//2, (S + res)//2]}

def pridaj(kluc,prvok,slovnik):
        if kluc not in slovnik: slovnik[kluc] = [prvok]
        else: slovnik[kluc].append(prvok)

def magicky_stvorec_3x3_5_stvorcov(u1,v1,u2,v2):
        p = ((u1*u1 + v1*v1)*(u2*u2 + 2*u2*v2 - v2*v2))**2
        q = ((u1*u1 + 2*u1*v1 - v1*v1)*(u2*u2 + v2*v2))**2
        r = ((- u1*u1 + 2*u1*v1 + v1*v1)*(u2*u2 + v2*v2))**2
        s = ((u1*u1 + v1*v1)*(- u2*u2 + 2*u2*v2 + v2*v2))**2
        t = ((u1*u1 + v1*v1)*(u2*u2 + v2*v2))**2
        x1 = 3*t - p - q
        x2 = 3*t - p - r
        x3 = 3*t - q - s
        x4 = 3*t - r - s
        y1 = x1
        y2 = r + s - p
        y3 = p + q - s
        y4 = x4

def magicky_stvorec_3x3_6_stvorcov(x):
        x2 = x*x
        x4 = x2*x2
        x6 = x4*x2
        x8 = x6*x4
        x10 = x8*x2
        
        v1 = (x2 - 2)*(8*x2 - 1)*(x6 - 6*x4 - 2)
        v2 = (x2 - 2)*(8*x8 - x6 + 30*x4 - 40*x2 + 2)
        v3 = (x2 - 2)*(8*x8 - 25*x6 + 18*x4 - 28*x2 + 2)
        v4 = (4*x10 - 31*x8 + 76*x6 + 76*x4 - 31*x2 + 4)//2
        v5 = (4*x10 + 17*x8 + 4*x6 + 4*x4 + 17*x2 + 4)//2
        v6 = (4*x10 + 65*x8 - 68*x6 - 68*x4 + 65*x2 + 4)//2

def bimagicky_stvorec_5x5(h):
        trojice = dict()
        komplet = dict()
        for a in range(h):
                for b in range(a+1,h):
                        for c in range(b+1,h):
                                index = a*a+b*b+c*c
                                if index not in trojice: trojice[index] = [(a,b,c)]
                                else: trojice[index].append((a,b,c))
                                index2 = (a+b+c,index)
                                if index2 not in komplet: komplet[index2] = [(a,b,c)]
                                else: komplet[index2].append((a,b,c))
                                
        for i,j in trojice.items():
                for j1 in range(len(j)):
                        for j2 in range(j1):
                                if len(set(j[j1] + j[j2])) == 6:
                                        x = (sum(j[j1]) - sum(j[j2]))//2
                                        if abs(x) not in j[j1] and abs(x) not in j[j2]:
                                                for p in range(h):
                                                        if p not in j[j1] and p not in j[j2] and p != abs(x):
                                                                index = (sum(j[j1]) - x + p, i - x*x + p*p)
                                                                if index in komplet:
                                                                        for k in komplet[index]:
                                                                                if len(set(k + j[j1] + j[j2] + (p,-x,x))) == 12:
                                                                                        print(j[j1][0],".",".",".",j[j2][0])
                                                                                        print(".",j[j1][1],".",j[j2][1],".")
                                                                                        print(".",".",p,".",".")
                                                                                        print(".",j[j2][2],".",j[j1][2],".")
                                                                                        print(x,k[0],k[1],k[2],-x)
                                                                                        print()
                                

def vrcholovo_bimagicky_graf_test(n):
	G = []
	G = read_graph6("graph" + str(n) + "c.g6")
	for i in range(len(G)):
		susedia = []
		for k in range(n): susedia.append(set())
		for e in G[i].edges():
			susedia[e[0]].add(e[1])
			susedia[e[1]].add(e[0])
		ok = True
		for v1 in range(n):
			for v2 in range(v1+1,n):
				x = len(susedia[v1].difference(susedia[v2]))
				y = len(susedia[v2].difference(susedia[v1]))
				if (x*y == 0 and x+y > 0) or x == 1 or y == 1 or (x == 2 and y == 2):
					ok = False
					break
			if not ok: break
		if ok: print(G[i].edges())

def bimagicky_obdlznik_3xN(hranica):
        trojice = dict()
        for a in range(2,hranica):
                for b in range(a+1,hranica):
                        for c in range(b+1,hranica):
                                s = a+b+c
                                t = a*a+b*b+c*c
                                v = 2*t - (s-1)*(s-1) - 2
                                if stvorec(v):
                                        w = isqrt(v)
                                        if s%2 != w%2:
                                                x1 = (s-1+w)//2
                                                x2 = (s-1-w)//2
                                                if min(x1,x2) > 1:
                                                        if x1 not in {a,b,c} and x2 not in {a,b,c}:
                                                                index = (min(x1,x2),max(x1,x2))
                                                                if index not in trojice: trojice[index] = [(a,b,c)]
                                                                else: trojice[index].append((a,b,c))
        maximum = 0
        for i,j in trojice.items(): maximum = max(maximum,len(j))
        for n in range(4,maximum+1):
                for x,trojica in trojice.items():
                        for C in combinations(trojica,r=n-1):
                                prvky = []
                                for c in C: prvky += c
                                if len(set(prvky).union({1,x[0],x[1]})) == 3*n:
                                        moznosti = []
                                        for i in range(len(C)):
                                                moznosti.append([])
                                                for P in permutations(C[i]): moznosti[i].append(P)
                                        for PP in product([y for y in range(6)],repeat=n-1):
                                                obdlznik = [[1,x[0],x[1]]]
                                                for ii in range(len(PP)): obdlznik.append(moznosti[ii][PP[ii]])
                                                S = set()
                                                T = set()
                                                for ii in range(3):
                                                        s = 0
                                                        t = 0
                                                        for jj in range(len(obdlznik)):
                                                             s += obdlznik[jj][ii]
                                                             t += obdlznik[jj][ii]*obdlznik[jj][ii]
                                                        S.add(s)
                                                        T.add(t)
                                                if len(S) < 3 and len(T) < 3:
                                                        if len(S) == len(T) == 1: print("bimagicky obdlznik 3 x",n,"|",obdlznik,S,T)
                                                        else: print("ciastocny bimagicky obdlznik 3 x",n,"|",obdlznik,S,T)

                                
def bimagicky_obdlznik_sucet_3xN(sucet):
        trojice = dict()
        for a in range(2,sucet//3+1):
                for b in range(a+1,(sucet-a)//2+1):
                        c = sucet-a-b
                        s = sucet
                        t = a*a+b*b+c*c
                        v = 2*t - (s-1)*(s-1) - 2
                        if stvorec(v):
                                w = isqrt(v)
                                if s%2 != w%2:
                                        x1 = (s-1+w)//2
                                        x2 = (s-1-w)//2
                                        if min(x1,x2) > 1:
                                                if x1 not in {a,b,c} and x2 not in {a,b,c}:
                                                        index = (min(x1,x2),max(x1,x2))
                                                        if index not in trojice: trojice[index] = [(a,b,c)]
                                                        else: trojice[index].append((a,b,c))
        maximum = 0
        for i,j in trojice.items(): maximum = max(maximum,len(j))
        for n in range(4,maximum+1):
                for x,trojica in trojice.items():
                        for C in combinations(trojica,r=n-1):
                                prvky = []
                                for c in C: prvky += c
                                if len(set(prvky).union({1,x[0],x[1]})) == 3*n:
                                        moznosti = []
                                        for i in range(len(C)):
                                                moznosti.append([])
                                                for P in permutations(C[i]): moznosti[i].append(P)
                                        for PP in product([y for y in range(6)],repeat=n-1):
                                                obdlznik = [[1,x[0],x[1]]]
                                                for ii in range(len(PP)): obdlznik.append(moznosti[ii][PP[ii]])
                                                S = set()
                                                T = set()
                                                for ii in range(3):
                                                        s = 0
                                                        t = 0
                                                        for jj in range(len(obdlznik)):
                                                             s += obdlznik[jj][ii]
                                                             t += obdlznik[jj][ii]*obdlznik[jj][ii]
                                                        S.add(s)
                                                        T.add(t)
                                                if len(S) == 1 and len(T) < 3:
                                                        if len(S) == len(T) == 1: print("bimagicky obdlznik 3 x",n,"|",obdlznik,S,T)
                                                        else: print("ciastocny bimagicky obdlznik 3 x",n,"|",obdlznik,S,T)       
        

def multiplikativny_obdlznik_3xN(hranica):
        trojice = dict()
        for a in range(1,hranica):
                if not isprime(a) or 4*a < hranica:
                        for b in range(a+1,hranica):
                                if not isprime(b) or 4*b < hranica:
                                        for c in range(b+1,hranica):
                                                if not isprime(c) or 4*c < hranica:
                                                        index = (a+b+c,a*b*c)
                                                        if index not in trojice: trojice[index] = [(a,b,c)]
                                                        else: trojice[index].append((a,b,c))
        maximum = 0
        for i,j in trojice.items(): maximum = max(maximum,len(j))
        for n in range(4,maximum+1):
                for i,j in trojice.items():
                        for C in combinations(j,r=n):
                                prvky = []
                                for c in C: prvky += c
                                if len(set(prvky)) == 3*n:
                                        moznosti = []
                                        for i in range(1,len(C)):
                                                moznosti.append([])
                                                for P in permutations(C[i]): moznosti[i-1].append(P)
                                        for PP in product([y for y in range(6)],repeat=n-1):
                                                obdlznik = [C[0]]
                                                for ii in range(len(PP)): obdlznik.append(moznosti[ii][PP[ii]])
                                                S = set()
                                                T = set()
                                                for ii in range(3):
                                                        s = 0
                                                        t = 1
                                                        for jj in range(len(obdlznik)):
                                                             s += obdlznik[jj][ii]
                                                             t *= obdlznik[jj][ii]
                                                        S.add(s)
                                                        T.add(t)
                                                if len(T) == 1:
                                                        if len(S) == len(T) == 1: print("multiplikativny magicky obdlznik 3 x",n,"|",obdlznik,S,T)
                                                        else: print("ciastocny multiplikativny magicky obdlznik 3 x",n,"|",obdlznik,S,T)
                                                        
                                                        
def multiplikativny_obdlznik_prvocisla_3xN(hranica,prvocisla):
        trojice = dict()
        vyhovuju = set()
        for x in range(1,hranica):
                ok = True
                for f in factorint(x):
                        if f not in prvocisla:
                                ok = False
                                break
                if ok: vyhovuju.add(x)
                
        for a in range(1,hranica):
                if a in vyhovuju:
                        for b in range(a+1,hranica):
                                 if b in vyhovuju:
                                        for c in range(b+1,hranica):
                                                 if c in vyhovuju:
                                                        index = (a+b+c,a*b*c)
                                                        if index not in trojice: trojice[index] = [(a,b,c)]
                                                        else: trojice[index].append((a,b,c))
        maximum = 0
        for i,j in trojice.items(): maximum = max(maximum,len(j))
        for n in range(4,maximum+1):
                for i,j in trojice.items():
                        for C in combinations(j,r=n):
                                prvky = []
                                for c in C: prvky += c
                                if len(set(prvky)) == 3*n:
                                        moznosti = []
                                        for i in range(1,len(C)):
                                                moznosti.append([])
                                                for P in permutations(C[i]): moznosti[i-1].append(P)
                                        for PP in product([y for y in range(6)],repeat=n-1):
                                                obdlznik = [C[0]]
                                                for ii in range(len(PP)): obdlznik.append(moznosti[ii][PP[ii]])
                                                S = set()
                                                T = set()
                                                for ii in range(3):
                                                        s = 0
                                                        t = 1
                                                        for jj in range(len(obdlznik)):
                                                             s += obdlznik[jj][ii]
                                                             t *= obdlznik[jj][ii]
                                                        S.add(s)
                                                        T.add(t)
                                                if len(T) == 1:
                                                        if len(S) == len(T) == 1: print("multiplikativny magicky obdlznik 3 x",n,"|",obdlznik,S,T)
                                                        else: print("ciastocny multiplikativny magicky obdlznik 3 x",n,"|",obdlznik,S,T)

def multiplikativny_obdlznik_sucet_3xN(sucet):
        trojice = dict()
        for a in range(1,sucet//3+1):
                for b in range(a+1,(sucet-a)//2+1):
                        c = sucet-a-b
                        index = a*b*c
                        if index not in trojice: trojice[index] = [(a,b,c)]
                        else: trojice[index].append((a,b,c))
        maximum = 0
        for i,j in trojice.items(): maximum = max(maximum,len(j))
        for n in range(4,maximum+1):
                for i,j in trojice.items():
                        for C in combinations(j,r=n):
                                prvky = []
                                for c in C: prvky += c
                                if len(set(prvky)) == 3*n:
                                        moznosti = []
                                        for i in range(1,len(C)):
                                                moznosti.append([])
                                                for P in permutations(C[i]): moznosti[i-1].append(P)
                                        for PP in product([y for y in range(6)],repeat=n-1):
                                                obdlznik = [C[0]]
                                                for ii in range(len(PP)): obdlznik.append(moznosti[ii][PP[ii]])
                                                S = set()
                                                T = set()
                                                for ii in range(3):
                                                        s = 0
                                                        t = 1
                                                        for jj in range(len(obdlznik)):
                                                             s += obdlznik[jj][ii]
                                                             t *= obdlznik[jj][ii]
                                                        S.add(s)
                                                        T.add(t)
                                                if len(T) == 1:
                                                        if len(S) == len(T) == 1: print("multiplikativny magicky obdlznik 3 x",n,"|",obdlznik,S,T)
                                                        else: print("ciastocny multiplikativny magicky obdlznik 3 x",n,"|",obdlznik,S,T)


def multiplikativny_obdlznik_sucet_4xN(sucet):
        stvorice = dict()
        for a in range(1,sucet//4+1):
                if not isprime(a) or 4*a < sucet:
                        for b in range(a+1,(sucet-a)//3+1):
                                if not isprime(b) or 4*b < sucet:
                                        for c in range(b+1,(sucet-a-b)//2+1):
                                                if not isprime(c) or 4*c < sucet:
                                                        d = sucet-a-b-c
                                                        index = a*b*c*d
                                                        if index not in stvorice: stvorice[index] = [(a,b,c,d)]
                                                        else: stvorice[index].append((a,b,c,d))
        maximum = 0
        for i,j in stvorice.items(): maximum = max(maximum,len(j))
        for n in range(4,maximum+1):
                for i,j in stvorice.items():
                        for C in combinations(j,r=n):
                                prvky = []
                                for c in C: prvky += c
                                if len(set(prvky)) == 4*n:
                                        moznosti = []
                                        for i in range(1,len(C)):
                                                moznosti.append([])
                                                for P in permutations(C[i]): moznosti[i-1].append(P)
                                        for PP in product([y for y in range(24)],repeat=n-1):
                                                obdlznik = [C[0]]
                                                for ii in range(len(PP)): obdlznik.append(moznosti[ii][PP[ii]])
                                                S = set()
                                                T = set()
                                                for ii in range(4):
                                                        s = 0
                                                        t = 1
                                                        for jj in range(len(obdlznik)):
                                                             s += obdlznik[jj][ii]
                                                             t *= obdlznik[jj][ii]
                                                        S.add(s)
                                                        T.add(t)
                                                if len(T) == 1:
                                                        if len(S) == len(T) == 1: print("multiplikativny magicky obdlznik 4 x",n,"|",obdlznik,S,T)
                                                        else: print("ciastocny multiplikativny magicky obdlznik 4 x",n,"|",obdlznik,S,T)


def vrcholovo_bimagicky_kompletny(i,j):
	if i > j: return ohodnot(j,i)
	if i <= 1 or i == j == 2: return "nie je mozne ohodnotit"
	if i == 2:
		H1 = [j*(j-1)//2 + 1, j*(j-1)*(3*j*j - 7*j + 14)//24]
		H2 = [k for k in range(1,j)]
		H2.append(j*(j-1)*(3*j*j - 7*j + 14)//24 + 1)
		return [H1,H2]
	if i == 3:
		H1 = [1, j*(j+1)//2 - 1, j*(j+1)*(3*j*j - 7*j - 14)//24 + 1]
		H2 = [k for k in range(2,j+1)]
		H2.append(j*(j+1)*(3*j*j - 7*j - 14)//24 + 2)
		return [H1,H2]
	if (i,j) == (4,4): return [[1,4,6,7], [2,3,5,8]]
	if (i,j) == (4,5): return [[2,12,13,15], [1,4,8,10,19]]
	H = vrcholovo_bimagicky_kompletny(i-2,j-3)
	m = max(max(H[0]),max(H[1])) + 1
	H[0] += [4*m, 5*m]
	H[1] += [m, 2*m, 6*m]
	return H

def vrcholovo_superbimagicky_kompletny(n):
        if n < 7: return "nie je mozne ohodnotit"
        if n%4 in {1,2}: return "nie je mozne ohodnotit"
        if n == 7: return [[1,2,4,7], [3,5,6]]
        if n == 8: return [[1,4,6,7], [2,3,5,8]]
        if n == 11: return [[1,3,4,5,9,11],[2,6,7,8,10]]
        if n == 12: return [[1,3,7,8,9,11],[2,4,5,6,10,12]]
        H = vrcholovo_superbimagicky_kompletny(n-8)
        for x in range(1,9):
                if x in {1,4,6,7}: H[0].append(n-8+x)
                else: H[1].append(n-8+x)
        return H

def vrcholovo_multiplikativny_magicky_kompletny(i,j):
        if i > j: return ohodnot(j,i)
        if i <= 1 or i == j == 2: return "nie je mozne ohodnotit"
        if (i,j) == (2,3): return [[5,12],[1,6,10]]
        if (i,j) == (2,4): return [[9,16],[1,2,4,18]]
        if i == 2:
                faktorial = 1
                for f in range(1,j): faktorial *= f
                H1 = [f + 1, f * (f + 1 - j*(j-1)//2)]
                H2 = [k for k in range(1,j)]
                H2.append((f + 1)*(f + 1 - j*(j-1)//2))
                return [H1,H2]
        if i == 3:
                faktorial = 1
                for f in range(1,j+1): faktorial *= f
                H1 = [1, f + 1, f * (f + 3 - j*(j+1)//2)]
                H2 = [k for k in range(2,j+1)]
                H2.append((f + 1)*(f+3 - j*(j+1)//2))
                return [H1,H2]
        if (i,j) == (4,4): return [[1,5,6,12], [2,3,4,15]]
        if (i,j) == (4,5): return [[2,10,20,27], [1,3,6,24,25]]
        H = vrcholovo_multiplikativny_magicky_kompletny(i-2,j-3)
        x = max(max(H[0]),max(H[1])) + 1
        y = max(max(H[0]),max(H[1])) + 2
        H[0] += [2*x*y, 2*x*y - x - y]
        H[1] += [2*(2*x*y - x - y), x, y]
        return H


