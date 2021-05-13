from networkx import read_graph6, enumerate_all_cliques
from itertools import product, permutations, combinations
from math import gcd, inf, isqrt
from sympy import factorint, isprime, divisors
from datetime import datetime
from random import choice, randint

def jeDruhouMocninou(n):
        if n < 0:
                return False
        return isqrt(n)*isqrt(n) == n

def rozne(P):
        return len(P) == len(set(P))

def pridaj(kluc,prvok,asociativnePole):
        if kluc not in asociativnePole:
                asociativnePole[kluc] = [prvok]
        else:
                asociativnePole[kluc].append(prvok)

def bimagickySucet(pole):
        odpoved = 0
        for prvok in pole:
                odpoved += prvok*prvok
        return odpoved

def sucin(pole):
        odpoved = 1
        for prvok in pole:
                odpoved *= prvok
        return odpoved

def asponSedemStvorcov(stvorec):
        if rozne(stvorec[0] + stvorec[1] + stvorec[2]):
                pocetStvorcov = 0
                for riadok in stvorec:
                     for prvok in riadok:
                             if jeDruhouMocninou(prvok):
                                     pocetStvorcov += 1
                if pocetStvorcov >= 7:
                        return True
                return False
        return False

def magickyStvorec3x3sPiatimiStvorcami(u1,v1,u2,v2):
        p = ((u1*u1 + v1*v1)*(u2*u2 + 2*u2*v2 - v2*v2))**2
        q = ((u1*u1 + 2*u1*v1 - v1*v1)*(u2*u2 + v2*v2))**2
        r = ((- u1*u1 + 2*u1*v1 + v1*v1)*(u2*u2 + v2*v2))**2
        s = ((u1*u1 + v1*v1)*(- u2*u2 + 2*u2*v2 + v2*v2))**2
        t = ((u1*u1 + v1*v1)*(u2*u2 + v2*v2))**2

        for stvorec in (([p,3*t-p-q,q],
                         [3*t-p-r,t,3*t-q-s],
                         [r,3*t-r-s,s]),
                        
                        ([2*(r+s),4*p,2*(q+s)],
                         [4*q,4*t,4*r],
                         [2*(p+r),4*s,2*(p+q)]),
                        
                        ([p,q,3*t-p-q],
                         [r+s-p,t,p+q-s],
                         [3*t-r-s,r,s]),

                        ([p,r,3*t-p-r],
                         [q+s-p,t,p+r-s],
                         [3*t-q-s,q,s])):
                
                if asponSedemStvorcov(stvorec):
                        print(stvorec)               

def magickyStvorec3x3soSiestimiStvorcami(x):
        X = [1]
        for _ in range(10):
                X.append(X[-1]*x)

        stvorec1 = ([0,0,0],
                    [0,0,0],
                    [0,0,0])
        
        stvorec1[0][0] = (2*X[5] + 4*X[3] - 7*X[1])**2
        stvorec1[0][1] = (X[2] - 2)*(8*X[2] - 1)*(X[6] - 6*X[4] - 2)
        stvorec1[0][2] = (5*X[4] - 2*X[2] + 2)**2
        stvorec1[1][0] = (X[4] + 8*X[2] - 2)**2
        stvorec1[1][1] = (2*X[5] - 2*X[3] + 5*X[1])**2
        stvorec1[1][2] = (X[2] - 2)*(8*X[8] - X[6] + 30*X[4] - 40*X[2] + 2)
        stvorec1[2][0] = (X[2] - 2)*(8*X[8] - 25*X[6] + 18*X[4] - 28*X[2] + 2)
        stvorec1[2][1] = (7*X[4] - 4*X[2] - 2)**2
        stvorec1[2][2] = (2*X[5] - 8*X[3] - X[1])**2

        stvorec2 = ([0,0,0],
                    [0,0,0],
                    [0,0,0])
        
        stvorec2[0][0] = stvorec1[1][0]
        stvorec2[0][1] = stvorec1[0][0]
        stvorec2[0][2] = (4*X[10] - 31*X[8] + 76*X[6] + 76*X[4] - 31*X[2] + 4)//2
        stvorec2[1][0] = stvorec1[2][2]
        stvorec2[1][1] = (4*X[10] + 17*X[8] + 4*X[6] + 4*X[4] + 17*X[2] + 4)//2
        stvorec2[1][2] = stvorec1[2][1]
        stvorec2[2][0] = (4*X[10] + 65*X[8] - 68*X[6] - 68*X[4] + 65*X[2] + 4)//2
        stvorec2[2][1] = stvorec1[0][2]
        stvorec2[2][2] = stvorec1[1][1]

        for stvorec in (stvorec1,stvorec2):
                if asponSedemStvorcov(stvorec):
                        print(stvorec)
        
def bimagickyStvorec5x5(h,zleSucty=3):
        T = dict()
        x = 1
        y = h
        spolu = 0

        for a in range(0,isqrt((y+2)//3)+1):
                for b in range(a+1,isqrt((y-a*a+1)//2)+1):
                        for c in range(max(b+1,isqrt(max(x-a*a-b*b,0))),isqrt(y-a*a-b*b)+1):
                                n = a*a + b*b + c*c
                                if x <= n and n < y:
                                        if n not in T: T[n] = [[a,b,c]]
                                        else: T[n].append([a,b,c])
                    
        for i,j in T.items():
                if len(j) >= 3: spolu += 4*len(j)*(len(j)-1)*(len(j)-2)//6

        counter = 0


        for KK,j in T.items():
                if len(j) >= 3:
                        K = 4*KK
                        for i3,i2,i1 in combinations([i for i in range(len(j))],r=3):
                                a = j[i1][0]
                                b = j[i1][1]
                                c = j[i1][2]
                                dd = j[i2][0]
                                ee = j[i2][1]
                                ff = j[i2][2]
                                gg = j[i3][0]
                                hh = j[i3][1]
                                ii = j[i3][2]
                                for d,e,f in ((dd,ee,ff),(-dd,-ee,-ff)):
                                        for g,h,i in ((gg,hh,ii),(-gg,-hh,-ii)):
                                                if counter%10000 == 0:
                                                        print()
                                                        print(counter,"z",spolu,"|",str(datetime.now()))
                                                        print()
                                                counter += 1
                                                hodnoty = dict()

                                                for Vabc in range(4):
                                                        V1 = [a+b-c,a-b+c,-a+b+c,-a-b-c][Vabc]
                                                        for Vdef in range(4):
                                                                V2 = [d+e-f,d-e+f,-d+e+f,-d-e-f][Vdef]
                                                                for Vghi in range(4):
                                                                        V3 = [g+h-i,g-h+i,-g+h+i,-g-h-i][Vghi]
                                                                        n = 4*(V1*V2 + V1*V3 + V2*V3 + V1*V1 + V2*V2 + V3*V3 - K//2)
                                                                        output = (Vabc,Vdef,Vghi,V1,V2,V3)

                                                                        if n == 0 or n%4 == 2: pass
                                                                        else:
                                                                                for D1 in divisors(abs(n),True):
                                                                                        D2 = abs(n)//D1
                                                                                        if D1%2 == D2%2 and D1 >= D2:
                                                                                                mm = ((D1+D2)//2,(-D1-D2)//2)
                                                                                                ll = ((D1-D2)//2,(D2-D1)//2)
                                                                                                if n < 0:
                                                                                                        mm,ll = ll,mm
                                                                                                for M in mm:
                                                                                                        for L in ll:
                                                                                                                S = -V1-V2-V3+M
                                                                                                                if (S-V1-V2-V3)%2 == L%2:
                                                                                                                        x1 = min((S-V1-V2-V3+L)//2,(S-V1-V2-V3-L)//2)
                                                                                                                        x2 = max((S-V1-V2-V3+L)//2,(S-V1-V2-V3-L)//2)
                                                                                                                        if S not in hodnoty: hodnoty[S] = {output + (x1,x2)}
                                                                                                                        else: hodnoty[S].add(output + (x1,x2))
                                                                                                                        if M == 0: break
                                                for s,l in hodnoty.items():
                                                        if len(l) >= 4:
                                                                
                                                                H = []
                                                                for ll in l: H.append(ll)
                                                                H.sort()
                                                                h0 = -1
                                                                h1 = -1
                                                                h2 = -1
                                                                h3 = -1
                                                                
                                                                for hh in range(len(H)):
                                                                        if H[hh][0] == 0 and h0 == -1: h0 = hh
                                                                        elif H[hh][0] == 1 and h0 != -1 and h1 == -1: h1 = hh
                                                                        elif H[hh][0] == 2 and h1 != -1 and h2 == -1: h2 = hh
                                                                        elif H[hh][0] == 3 and h2 != -1 and h3 == -1: h3 = hh

                                                                if h0 != -1 and h1 != -1 and h2 != -1 and h3 != -1:
                                                                        for k0 in range(h0,h1):
                                                                                for k1 in range(h1,h2):
                                                                                        for k2 in range(h2,h3):
                                                                                                for k3 in range(h3,len(H)):
                                                                                                        if len({H[k0][1],H[k1][1],H[k2][1],H[k3][1]}) == 4:
                                                                                                                if len({H[k0][2],H[k1][2],H[k2][2],H[k3][2]}) == 4:
                                                                                                                        pouzite = {s}.union(set(H[k0][3:8]).union(set(H[k1][3:8]).union(set(H[k2][3:8]).union(set(H[k3][3:8])))))
                                                                                                                        if len(pouzite) == 21:
                                                                                                                                s2 = K + s*s
                                                                                                                                R0 = [[H[k0][3],H[k0][4],H[k0][5]],[H[k0][6],H[k0][7]]]
                                                                                                                                R1 = [[H[k1][3],H[k1][4],H[k1][5]],[H[k1][6],H[k1][7]]]
                                                                                                                                R2 = [[H[k2][3],H[k2][4],H[k2][5]],[H[k2][6],H[k2][7]]]
                                                                                                                                R3 = [[H[k3][3],H[k3][4],H[k3][5]],[H[k3][6],H[k3][7]]]
                                                                                                                                for P in [[1,0,2],[0,1,2],[0,2,1]]:
                                                                                                                                        for R in [[R0,R1,R2,R3],[R0,R1,R3,R2],[R0,R2,R3,R1]]:
                                                                                                                                                for x0,y0 in [[R[0][1][0],R[0][1][1]],[R[0][1][1],R[0][1][0]]]:
                                                                                                                                                        for x1,y1 in [[R[1][1][0],R[1][1][1]],[R[1][1][1],R[1][1][0]]]:
                                                                                                                                                                for x2,y2 in [[R[2][1][0],R[2][1][1]],[R[2][1][1],R[2][1][0]]]:
                                                                                                                                                                        for x3,y3 in [[R[3][1][0],R[3][1][1]],[R[3][1][1],R[3][1][0]]]:
                                                                                                                                                                                stvorec = [[R[0][0][P[0]],x0,R[0][0][P[1]],y0,R[0][0][P[2]]],
                                                                                                                                                                                           [x1,R[1][0][P[0]],R[1][0][P[1]],R[1][0][P[2]],y1],
                                                                                                                                                                                           [-1,-1,s,-1,-1],
                                                                                                                                                                                           [x2,R[2][0][P[2]],R[2][0][P[1]],R[2][0][P[0]],y2],
                                                                                                                                                                                           [R[3][0][P[2]],x3,R[3][0][P[1]],y3,R[3][0][P[0]]]]
                                                                                                                                                                                
                                                                                                                                                                                for stlpec in {0,1,3,4}:
                                                                                                                                                                                        sucetStlpec = 0
                                                                                                                                                                                        for riadok in {0,1,3,4}:
                                                                                                                                                                                                sucetStlpec += stvorec[riadok][stlpec]
                                                                                                                                                                                        stvorec[2][stlpec] = s - sucetStlpec

                                                                                                                                                                                if len(pouzite.union({stvorec[2][0],stvorec[2][1],stvorec[2][3],stvorec[2][4]})) == 25:
                                                                                                                                                                                        nespravne = 5
                                                                                                                                                                                        for stlpec in {0,1,3,4}:
                                                                                                                                                                                                bimagickySucetStlpec = 0
                                                                                                                                                                                                for riadok in range(5):
                                                                                                                                                                                                        bimagickySucetStlpec += stvorec[riadok][stlpec]**2
                                                                                                                                                                                                if bimagickySucetStlpec == s2:
                                                                                                                                                                                                        nespravne -= 1

                                                                                                                                                                                        bimagickySucetRiadok = 0
                                                                                                                                                                                        riadok = 2
                                                                                                                                                                                        for stlpec in range(5):
                                                                                                                                                                                                bimagickySucetRiadok += stvorec[riadok][stlpec]**2
                                                                                                                                                                                        if bimagickySucetRiadok == s2:
                                                                                                                                                                                                nespravne -= 1
                                                                                                                                                                                        
                                                                                                                                                                                        if nespravne <= zleSucty:
                                                                                                                                                                                                print("magicky sucet:",s)
                                                                                                                                                                                                print("bimagicky sucet:",s2)
                                                                                                                                                                                                print("pocet nespravnych suctov:",nespravne)
                                                                                                                                                                                                for riadok in stvorec: print(riadok)
                                                                                                                                                                                                print()



def variacneRozpatie(stvorec):
        velkost = len(stvorec)
        prvky = []
        for riadok in stvorec:
                prvky += riadok
        if not rozne(prvky):
                return (inf,inf)
        
        sucty = []
        sucetDiagonala1 = 0
        sucetDiagonala2 = 0
        
        for i in range(velkost):
                sucetDiagonala1 += stvorec[i][i]
                sucetDiagonala2 += stvorec[i][velkost-i-1]
                
                sucetRiadok = 0
                sucetStlpec = 0
                for j in range(velkost):
                        sucetRiadok += stvorec[i][j]
                        sucetStlpec += stvorec[j][i]
                sucty.append(sucetRiadok)
                sucty.append(sucetStlpec)

        sucty.append(sucetDiagonala1)
        sucty.append(sucetDiagonala2)
        return (max(sucty) - min(sucty),len(set(sucty)))


def suVzorkyDisjunktne(vzorky):
        odpoved = True
        pocetVzoriek = len(vzorky)
        if pocetVzoriek > 0:
                velkostVzorky = len(vzorky[0])
                for index in range(velkostVzorky):
                        vybratePolicka = []
                        for vzorka in vzorky:
                                vybratePolicka.append(vzorka[index])
                        if not rozne(vybratePolicka):
                                odpoved = False
                                break
        return odpoved

def generujStvorcoveVzorky(vsetkyVzorky,vybraneVzorky=[],vybraneIndexy=[]):
        pocetVsetkych = len(vsetkyVzorky)
        velkostVzorky = len(vsetkyVzorky[0])
        
        while True:
                pocetVybranych = len(vybraneVzorky)
                if suVzorkyDisjunktne(vybraneVzorky) and pocetVybranych < velkostVzorky:
                        if pocetVybranych == 0:
                                vybraneVzorky.append(vsetkyVzorky[0])
                                vybraneIndexy.append(0)
                        else:
                                vybraneVzorky.append(vsetkyVzorky[vybraneIndexy[-1]+1])
                                vybraneIndexy.append(vybraneIndexy[-1]+1)
                        pocetVybranych += 1
                else:
                        if suVzorkyDisjunktne(vybraneVzorky):
                                yield vybraneVzorky

                        while pocetVybranych > 0:
                                poslednyIndex = vybraneIndexy.pop()
                                vybraneVzorky.pop()
                                pocetVybranych -= 1
                                if poslednyIndex < pocetVsetkych - (velkostVzorky - pocetVybranych):
                                        vybraneVzorky.append(vsetkyVzorky[poslednyIndex+1])
                                        vybraneIndexy.append(poslednyIndex+1)
                                        pocetVybranych += 1
                                        break
                if pocetVybranych == 0:
                        break
                

def prenasobStvorec(hodnotyVzoriek,pouziteVzorky):
        velkost = len(pouziteVzorky[0])
        stvorec = []
        for riadok in range(velkost):
                stvorec.append([])
                for _ in range(velkost):
                        stvorec[riadok].append(1)

        for index in range(len(hodnotyVzoriek)):
                for riadok in range(velkost):
                        stvorec[riadok][pouziteVzorky[index][riadok]] *= hodnotyVzoriek[index]

        return stvorec

def generujVzorky(velkost):
        for vzorka in permutations((i for i in range(velkost))):
                diagonala1 = 0
                diagonala2 = 0
                for v in range(len(vzorka)):
                        if vzorka[v] == v: diagonala1 += 1
                        if vzorka[v] == velkost-v-1: diagonala2 += 1
                if diagonala1 == diagonala2 == 1:
                        yield vzorka


def vymenDveHodnoty(hodnotyVzoriek):
        for index1,index2 in permutations([i for i in range(len(hodnotyVzoriek))],2):
                pamat1,pamat2 = hodnotyVzoriek[index1],hodnotyVzoriek[index2]
                hodnotyVzoriek[index1],hodnotyVzoriek[index2] = pamat2,pamat1
                yield hodnotyVzoriek
                hodnotyVzoriek[index1],hodnotyVzoriek[index2] = pamat1,pamat2


def zmenHodnotu(hodnotyVzoriek,h):
        for index in range(len(hodnotyVzoriek)):
                pamat = hodnotyVzoriek[index]
                for novaHodnota in range(1,h+1):
                        hodnotyVzoriek[index] = novaHodnota
                        yield hodnotyVzoriek
                hodnotyVzoriek[index] = pamat

def pripocitajRozsah(hodnotyVzoriek,velkost,rozsah):
        for indexStvorcovejVzorky in range(0,len(hodnotyVzoriek),velkost):
                for pripocitat in product([i for i in range(-rozsah,rozsah+1)],repeat=velkost):
                        pamat = []
                        for hodnotaVzorky in hodnotyVzoriek:
                                pamat.append(hodnotaVzorky)

                        for i in range(len(pripocitat)):
                                hodnotyVzoriek[indexStvorcovejVzorky+i] += pripocitat[i]

                        yield hodnotyVzoriek
                                
                        for i in range(len(pripocitat)):
                                hodnotyVzoriek[indexStvorcovejVzorky+i] = pamat[indexStvorcovejVzorky+i]  
                                
def multiplikativnyMagickyStvorec(h,p):
        velkost = 6
        vzorky = []
        stvorcoveVzorky = []
        
        for vzorka in generujVzorky(velkost):
                vzorky.append(tuple(vzorka))
        
        for stvorcovaVzorka in generujStvorcoveVzorky(vzorky):
                stvorcoveVzorky.append(tuple(stvorcovaVzorka))

        rozsah = 1        
        stavNajlepsi = (inf,inf)
        optimum = (0,0)
        
        while True:
                hodnotyVzoriek = []
                pouziteVzorky = []
                for i in range(p):
                        hodnotyVzoriek += [randint(1,h) for _ in range(velkost)]
                        pouziteVzorky += choice(stvorcoveVzorky)
                        
                stvorec = prenasobStvorec(hodnotyVzoriek,pouziteVzorky)
                stav = variacneRozpatie(stvorec)

                while stav > optimum:
                        stavZaciatok = stav
                        if stav != (inf,inf):
                                for operacia in {vymenDveHodnoty(hodnotyVzoriek),
                                                 zmenHodnotu(hodnotyVzoriek,h),
                                                 pripocitajRozsah(hodnotyVzoriek,velkost,rozsah)}:
                                        for hodnotyVzoriekNove in operacia:
                                                stvorecNovy = prenasobStvorec(hodnotyVzoriekNove,pouziteVzorky)
                                                stavNovy = variacneRozpatie(stvorecNovy)
                                                if stavNovy < stav:
                                                        stav = stavNovy
                                                        stvorec = [riadok for riadok in stvorecNovy]
                                                        hodnotyVzoriek = [hodnota for hodnota in hodnotyVzoriekNove]

                        if stav == stavZaciatok:
                                if stav < stavNajlepsi:
                                        stavNajlepsi = stav
                                        print("rozpatie suctov:",stavNajlepsi[0])
                                        print("pocet roznych suctov:",stavNajlepsi[1])
                                        for riadok in stvorec:
                                                print(riadok)
                                        print()
                                break
                                

def vrcholovoBimagickyGrafTest(n):
	G = []
	G = read_graph6("graph" + str(n) + "c.g6")
	
	for i in range(len(G)):
		susedia = []
		for _ in range(n): susedia.append(set())
		for e in G[i].edges():
			susedia[e[0]].add(e[1])
			susedia[e[1]].add(e[0])
			
		vyhovuje = True
		for v1 in range(n):
			for v2 in range(v1+1,n):
				x = len(susedia[v1].difference(susedia[v2]))
				y = len(susedia[v2].difference(susedia[v1]))
				if (x*y == 0 and x+y > 0) or x == 1 or y == 1 or (x == 2 and y == 2):
					vyhovuje = False
					break
			if not vyhovuje: break
		if vyhovuje: print(G[i].edges())


def najdiJednotkovuTrojicu(a,b,c):
        s = a+b+c
        t = a*a+b*b+c*c
        vyraz = 2*t - (s-1)*(s-1) - 2
        if jeDruhouMocninou(vyraz):
                n = isqrt(vyraz)
                if s%2 != n%2:
                        x1 = (s-1+n)//2
                        x2 = (s-1-n)//2
                        if min(x1,x2) > 1:
                                if rozne({1,x1,x2,a,b,c}):
                                        return ((1,min(x1,x2),max(x1,x2)),(a,b,c))

def generujBimagickeTrojice(hranica):
        for a,b,c in combinations([i for i in range(2,hranica+1)],r=3):
                vysledok = najdiJednotkovuTrojicu(a,b,c)
                if vysledok is not None:
                        yield vysledok

def generujTrojiceSoSuctom(sucet):
        for a in range(2,sucet//3+1):
                for b in range(a+1,(sucet-a)//2+1):
                        c = sucet-a-b
                        yield (a,b,c)
                        
def generujBimagickeTrojiceSoSuctom(sucet):
        for a,b,c in generujTrojiceSoSuctom(sucet):
                vysledok = najdiJednotkovuTrojicu(a,b,c)
                if vysledok is not None:
                        yield vysledok

def generujStvoriceSoSuctom(sucet):
        for a in range(2,sucet//4+1):
                for b in range(a+1,(sucet-a)//3+1):
                        for c in range(b+1,(sucet-a-b)//2+1):
                                d = sucet-a-b-c
                                yield (a,b,c,d)

def generujNticeSoSuctom(n,sucet):
        if n*(n+1)//2 > sucet:
                return
        
        Ntica = [i for i in range(1,n)]
        Ntica.append(sucet - n*(n-1)//2)
        
        while True:
                if Ntica[n-1] - Ntica[n-2] > 1:
                        Ntica[n-1] -= 1
                        Ntica[n-2] += 1
                else:
                        posledny = n-2
                        while True:
                                if posledny <= 0:
                                        return
                                if Ntica[posledny] - Ntica[posledny-1] > 1:
                                        Ntica[posledny-1] += 1
                                        for dalsi in range(posledny,n-1):
                                                Ntica[dalsi] = Ntica[dalsi-1] + 1
                                        sucetOkremPosledneho = sum(Ntica) - Ntica[n-1]
                                        Ntica[n-1] = sucet - sucetOkremPosledneho
                                        break
                                else:
                                        posledny -= 1 
                yield Ntica
        
def bimagickyObdlznik3xN(n,h):
        trojice = dict()

        for trojicaSJednotkou,trojicaBezJednotky in generujBimagickeTrojice(h):
                pridaj(trojicaSJednotkou,trojicaBezJednotky,trojice)

        for x,trojica in trojice.items():
                for C in combinations(trojica,r=n-1):
                        prvky = []
                        for c in C: prvky += c
                        if len(set(prvky).union(set(x))) == 3*n:
                                moznosti = []
                                for i in range(len(C)):
                                        moznosti.append([])
                                        for P in permutations(C[i]): moznosti[i].append(P)
                                for PP in product([y for y in range(6)],repeat=n-1):
                                        obdlznik = [x]
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
                                        if len(S) == 1:
                                                if len(S) == len(T) == 1: print("bimagicky obdlznik 3 x",n,"|",obdlznik,S,T)
                                                elif len(T) < 3: print("ciastocny bimagicky obdlznik 3 x",n,"|",obdlznik,S,T)    
        
def bimagickyObdlznikSucet3xN(n,s):
        trojice = dict()

        for trojicaSJednotkou,trojicaBezJednotky in generujBimagickeTrojiceSoSuctom(s):
                pridaj(trojicaSJednotkou,trojicaBezJednotky,trojice)
                
        for x,trojica in trojice.items():
                for C in combinations(trojica,r=n-1):
                        prvky = []
                        for c in C: prvky += c
                        if len(set(prvky).union(set(x))) == 3*n:
                                moznosti = []
                                for i in range(len(C)):
                                        moznosti.append([])
                                        for P in permutations(C[i]): moznosti[i].append(P)
                                for PP in product([y for y in range(6)],repeat=n-1):
                                        obdlznik = [x]
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


def bimagickyObdlznikSucet4xN(n,s):
        trojice = dict()
        
        for trojica in generujTrojiceSoSuctom(s-1):
                magickySucetTrojice = sum(trojica)
                bimagickySucetTrojice = bimagickySucet(trojica)
                pridaj((magickySucetTrojice,bimagickySucetTrojice),trojica,trojice)

        stvorice = dict()

        for a,b,c,d in generujStvoriceSoSuctom(s):
                t = a*a+b*b+c*c+d*d
                if (s-1,t-1) in trojice:
                        for x,y,z in trojice[(s-1,t-1)]:
                                if rozne({a,b,c,d,1,x,y,z}):
                                        pridaj((1,x,y,z),(a,b,c,d),stvorice)

        for x,stvorica in stvorice.items():
                for C in combinations(stvorica,r=n-1):
                        prvky = []
                        for c in C: prvky += c
                        if len(set(prvky).union(set(x))) == 4*n:
                                moznosti = []
                                for i in range(len(C)):
                                        moznosti.append([])
                                        for P in permutations(C[i]): moznosti[i].append(P)
                                for PP in product([y for y in range(24)],repeat=n-1):
                                        obdlznik = [x]
                                        for ii in range(len(PP)): obdlznik.append(moznosti[ii][PP[ii]])
                                        S = set()
                                        T = set()
                                        for ii in range(4):
                                                s = 0
                                                t = 0
                                                for jj in range(len(obdlznik)):
                                                     s += obdlznik[jj][ii]
                                                     t += obdlznik[jj][ii]*obdlznik[jj][ii]
                                                S.add(s)
                                                T.add(t)
                                        if len(S) == 1 and len(T) < 3:
                                                if len(S) == len(T) == 1: print("bimagicky obdlznik 4 x",n,"|",obdlznik,S,T)
                                                else: print("ciastocny bimagicky obdlznik 4 x",n,"|",obdlznik,S,T)  

def generujVyhovujuceMultiplikativneMagickeNticeSoSuctom(n,sucet):
        vyhovuju = set()
        for cislo in range(1,sucet+1):
                if not isprime(cislo) or n*cislo <= sucet:
                        vyhovuju.add(cislo)
                        
        for Ntica in generujNticeSoSuctom(n,sucet):
                vyhovuje = True
                for prvok in Ntica:
                        if prvok not in vyhovuju:
                                vyhovuje = False
                                break
                if vyhovuje:
                        p = sucin(Ntica)
                        yield (p,tuple(Ntica))

def multiplikativnyMagickyObdlznikSucet3xN(n,sucet):
        trojice = dict()

        for sucinTrojice,trojica in generujVyhovujuceMultiplikativneMagickeNticeSoSuctom(3,sucet):
                pridaj(sucinTrojice,trojica,trojice)

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


def multiplikativnyMagickyObdlznikSucet4xN(n,sucet):
        stvorice = dict()

        for sucinStvorice,stvorica in generujVyhovujuceMultiplikativneMagickeNticeSoSuctom(4,sucet):
                pridaj(sucinStvorice,stvorica,stvorice)
                                                        
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
                                        if len(S) < 3:
                                                if len(S) == len(T) == 1: print("multiplikativny magicky obdlznik 4 x",n,"|",obdlznik,S,T)
                                                else: print("ciastocny multiplikativny magicky obdlznik 4 x",n,"|",obdlznik,S,T)

def vrcholovoBimagickyKompletny(i,j):
        if i > j:
                return vrcholovoBimagickyKompletny(j,i)
        if i <= 1 or i == j == 2:
                return "nie je mozne ohodnotit"
        
        if i == 2:
                H1 = [j*(j-1)//2 + 1, j*(j-1)*(3*j*j - 7*j + 14)//24]
                H2 = [k for k in range(1,j)]
                H2.append(j*(j-1)*(3*j*j - 7*j + 14)//24 + 1)
                return (H1,H2)
        
        if i == 3:
                H1 = [1, j*(j+1)//2 - 1, j*(j+1)*(3*j*j - j - 14)//24 + 1]
                H2 = [k for k in range(2,j+1)]
                H2.append(j*(j+1)*(3*j*j - j - 14)//24 + 2)
                return (H1,H2)
        
        if (i,j) == (4,4):
                return ([1,4,6,7], [2,3,5,8])
        if (i,j) == (4,5):
                return ([2,12,13,15], [1,4,8,10,19])
        
        H = vrcholovoBimagickyKompletny(i-2,j-3)
        m = max(max(H[0]),max(H[1])) + 1
        H = (H[0] + [4*m, 5*m], H[1] + [m, 2*m, 6*m])
        return H
                      

def vrcholovoSuperbimagickyKompletny(n):
        if n < 7:
                return "nie je mozne ohodnotit"
        if n%4 in {1,2}:
                return "nie je mozne ohodnotit"
        if n == 7:
                return ([1,2,4,7], [3,5,6])
        if n == 8:
                return ([1,4,6,7], [2,3,5,8])
        if n == 11:
                return ([1,3,4,5,9,11], [2,6,7,8,10])
        if n == 12:
                return ([1,3,7,8,9,11], [2,4,5,6,10,12])
        
        H = vrcholovoSuperbimagickyKompletny(n-8)
        for x in range(1,9):
                if x in {1,4,6,7}: H[0].append(n-8+x)
                else: H[1].append(n-8+x)
        return H


def vrcholovoMultiplikativnyMagickyKompletny(i,j):
        if i > j:
                return vrcholovoMultiplikativnyMagickyKompletny(j,i)
        if i <= 1 or i == j == 2:
                return "nie je mozne ohodnotit"
        if (i,j) == (2,3):
                return ([5,12], [1,6,10])
        if (i,j) == (2,4):
                return ([9,16], [1,2,4,18])
        
        if i == 2:
                faktorial = 1
                for f in range(1,j): faktorial *= f
                H1 = [faktorial + 1, faktorial * (faktorial + 1 - j*(j-1)//2)]
                H2 = [k for k in range(1,j)]
                H2.append((faktorial + 1)*(faktorial + 1 - j*(j-1)//2))
                return (H1,H2)
        
        if i == 3:
                faktorial = 1
                for f in range(1,j+1): faktorial *= f
                H1 = [1, faktorial + 1, faktorial * (faktorial + 3 - j*(j+1)//2)]
                H2 = [k for k in range(2,j+1)]
                H2.append((faktorial + 1)*(faktorial + 3 - j*(j+1)//2))
                return (H1,H2)
        
        if (i,j) == (4,4):
                return ([1,5,6,12], [2,3,4,15])
        if (i,j) == (4,5):
                return ([2,10,20,27], [1,3,6,24,25])
        
        H = vrcholovoMultiplikativnyMagickyKompletny(i-2,j-3)
        x = max(max(H[0]),max(H[1])) + 1
        y = max(max(H[0]),max(H[1])) + 2
        H = (H[0] + [2*x*y, 2*x*y - x - y], H[1] + [2*(2*x*y - x - y), x, y])
        return H
