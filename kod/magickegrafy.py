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

def vypisRiesenie(nazovUtvaru,utvar,dodatocneUdaje):
        print(nazovUtvaru)
        dlzkaNajvacsiehoPrvku = 0
        for castUtvaru in utvar:
                for prvok in castUtvaru:
                       dlzkaNajvacsiehoPrvku = max(dlzkaNajvacsiehoPrvku, len(str(prvok)))
                       
        for castUtvaru in utvar:
                vypis = ""
                for index in range(len(castUtvaru)):
                        vypis += str(castUtvaru[index]) + (dlzkaNajvacsiehoPrvku - len(str(castUtvaru[index])) + 1)*" "
                print(vypis)        
                        
        for dodatocnyUdaj in dodatocneUdaje:
                print(dodatocnyUdaj)
        print()

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

#Algoritmus 4.1
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
                        vypisRiesenie("magicky stvorec velkosti 3 x 3 s aspon 7 druhymi mocninami",
                                      stvorec,
                                      set())               

#Algoritmus 4.2
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
                        vypisRiesenie("magicky stvorec velkosti 3 x 3 s aspon 7 druhymi mocninami",
                                      stvorec,
                                      set())

#Algoritmus 4.3        
def bimagickyStvorec5x5(h,zleSucty=3):
        trojice = dict()
        x = 1
        y = h
        spolu = 0

        for a in range(0,isqrt((y+2)//3)+1):
                for b in range(a+1,isqrt((y-a*a+1)//2)+1):
                        for c in range(max(b+1,isqrt(max(x-a*a-b*b,0))),isqrt(y-a*a-b*b)+1):
                                n = a*a + b*b + c*c
                                if x <= n and n < y:
                                        pridaj(n,[a,b,c],trojice)

        for KK,j in trojice.items():
                if len(j) >= 3:
                        K = 4*KK
                        for i3,i2,i1 in combinations([i for i in range(len(j))],r=3):
                                a,b,c = j[i1][0],j[i1][1],j[i1][2]
                                dd,ee,ff = j[i2][0],j[i2][1],j[i2][2]
                                gg,hh,ii = j[i3][0],j[i3][1],j[i3][2]
                                for d,e,f in ((dd,ee,ff),(-dd,-ee,-ff)):
                                        for g,h,i in ((gg,hh,ii),(-gg,-hh,-ii)):
                                                mozneProstrednePrvky = dict()
                                                diagonala1 = (a+b-c,a-b+c,-a+b+c,-a-b-c)
                                                prostrednyStlpec = (d+e-f,d-e+f,-d+e+f,-d-e-f)
                                                diagonala2 = (g+h-i,g-h+i,-g+h+i,-g-h-i)
                                                pouzite = diagonala1 + prostrednyStlpec + diagonala2
                                                
                                                if rozne(pouzite):
                                                        for pi in range(4):
                                                                p = diagonala1[pi]
                                                                for qi in range(4):
                                                                        q = prostrednyStlpec[qi]
                                                                        for ri in range(4):
                                                                                r = diagonala2[ri]
                                                                                vyraz = 4*(p*q + p*r + q*r + p**2 + q**2 + r**2) - 2*K
                                                                                udaje = (pi,qi,ri,p,q,r)

                                                                                if vyraz != 0:
                                                                                        for D1 in divisors(abs(vyraz),True):
                                                                                                D2 = abs(vyraz)//D1
                                                                                                if D1%2 == D2%2 and D1 >= D2:
                                                                                                        sucetRovnic = ((D1+D2)//2,(-D1-D2)//2)
                                                                                                        rozdielRovnic = ((D1-D2)//2,(D2-D1)//2)
                                                                                                        if vyraz < 0:
                                                                                                                sucetRovnic,rozdielRovnic = rozdielRovnic,sucetRovnic
                                                                                                        for m in sucetRovnic:
                                                                                                                for n in rozdielRovnic:
                                                                                                                        s = -p-q-r+m
                                                                                                                        if (s-p-q-r)%2 == n%2:
                                                                                                                                x1 = (s-(p+q+r)+n)//2
                                                                                                                                x2 = (s-(p+q+r)-n)//2
                                                                                                                                x = min(x1,x2)
                                                                                                                                y = max(x1,x2)
                                                                                                                                if rozne(pouzite+(x,y,s)):
                                                                                                                                        pridaj(s,udaje+(x,y),mozneProstrednePrvky)
                                                                                                                                if m == 0: break
                                                for s,vsetkyUdaje in mozneProstrednePrvky.items():
                                                        if len(set(vsetkyUdaje)) >= 4:
                                                                H = []
                                                                for udaj in set(vsetkyUdaje):
                                                                        H.append(udaj)
                                                                H.sort()
                                                                prvyUdajZacinajuci0 = -1
                                                                prvyUdajZacinajuci1 = -1
                                                                prvyUdajZacinajuci2 = -1
                                                                prvyUdajZacinajuci3 = -1
                                                                
                                                                for hh in range(len(H)):
                                                                        if H[hh][0] == 0 and prvyUdajZacinajuci0 == -1:
                                                                                prvyUdajZacinajuci0 = hh
                                                                        elif H[hh][0] == 1 and prvyUdajZacinajuci0 != -1 and prvyUdajZacinajuci1 == -1:
                                                                                prvyUdajZacinajuci1 = hh
                                                                        elif H[hh][0] == 2 and prvyUdajZacinajuci1 != -1 and prvyUdajZacinajuci2 == -1:
                                                                                prvyUdajZacinajuci2 = hh
                                                                        elif H[hh][0] == 3 and prvyUdajZacinajuci2 != -1 and prvyUdajZacinajuci3 == -1:
                                                                                prvyUdajZacinajuci3 = hh

                                                                if prvyUdajZacinajuci0 != -1 and prvyUdajZacinajuci1 != -1 and prvyUdajZacinajuci2 != -1 and prvyUdajZacinajuci3 != -1:
                                                                        for udajZacinajuci0 in range(prvyUdajZacinajuci0,prvyUdajZacinajuci1):
                                                                                for udajZacinajuci1 in range(prvyUdajZacinajuci1,prvyUdajZacinajuci2):
                                                                                        for udajZacinajuci2 in range(prvyUdajZacinajuci2,prvyUdajZacinajuci3):
                                                                                                for udajZacinajuci3 in range(prvyUdajZacinajuci3,len(H)):
                                                                                                        if len({H[udajZacinajuci0][1],H[udajZacinajuci1][1],H[udajZacinajuci2][1],H[udajZacinajuci3][1]}) == 4:
                                                                                                                if len({H[udajZacinajuci0][2],H[udajZacinajuci1][2],H[udajZacinajuci2][2],H[udajZacinajuci3][2]}) == 4:
                                                                                                                        pouzite = [s] + list(H[udajZacinajuci0][3:8]) + list(H[udajZacinajuci1][3:8]) + list(H[udajZacinajuci2][3:8]) + list(H[udajZacinajuci3][3:8])
                                                                                                                        if rozne(pouzite):
                                                                                                                                s2 = K + s*s
                                                                                                                                R0 = [[H[udajZacinajuci0][3],H[udajZacinajuci0][4],H[udajZacinajuci0][5]],[H[udajZacinajuci0][6],H[udajZacinajuci0][7]]]
                                                                                                                                R1 = [[H[udajZacinajuci1][3],H[udajZacinajuci1][4],H[udajZacinajuci1][5]],[H[udajZacinajuci1][6],H[udajZacinajuci1][7]]]
                                                                                                                                R2 = [[H[udajZacinajuci2][3],H[udajZacinajuci2][4],H[udajZacinajuci2][5]],[H[udajZacinajuci2][6],H[udajZacinajuci2][7]]]
                                                                                                                                R3 = [[H[udajZacinajuci3][3],H[udajZacinajuci3][4],H[udajZacinajuci3][5]],[H[udajZacinajuci3][6],H[udajZacinajuci3][7]]]
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

                                                                                                                                                                                if rozne(pouzite + [stvorec[2][0],stvorec[2][1],stvorec[2][3],stvorec[2][4]]):
                                                                                                                                                                                        nespravneBimagickeSucty = 5
                                                                                                                                                                                        for stlpec in {0,1,3,4}:
                                                                                                                                                                                                bimagickySucetStlpec = 0
                                                                                                                                                                                                for riadok in range(5):
                                                                                                                                                                                                        bimagickySucetStlpec += stvorec[riadok][stlpec]**2
                                                                                                                                                                                                if bimagickySucetStlpec == s2:
                                                                                                                                                                                                        nespravneBimagickeSucty -= 1

                                                                                                                                                                                        bimagickySucetRiadok = 0
                                                                                                                                                                                        riadok = 2
                                                                                                                                                                                        for stlpec in range(5):
                                                                                                                                                                                                bimagickySucetRiadok += stvorec[riadok][stlpec]**2
                                                                                                                                                                                        if bimagickySucetRiadok == s2:
                                                                                                                                                                                                nespravneBimagickeSucty -= 1
                                                                                                                                                                                        
                                                                                                                                                                                        if nespravneBimagickeSucty <= zleSucty:
                                                                                                                                                                                                if nespravneBimagickeSucty == 0:
                                                                                                                                                                                                        vypisRiesenie("bimagicky stvorec velkosti 5 x 5 so zapornymi prvkami",
                                                                                                                                                                                                                      stvorec,
                                                                                                                                                                                                                      set())
                                                                                                                                                                                                else:
                                                                                                                                                                                                        vypisRiesenie("magicky stvorec velkosti 5 x 5 so zapornymi prvkami",
                                                                                                                                                                                                                      stvorec,
                                                                                                                                                                                                                      {"pocet nespravnych bimagickych suctov: " + str(nespravneBimagickeSucty)})



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

#Algoritmus 4.4                                
def multiplikativnyMagickyStvorec(p,h):
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
                                        vypisRiesenie("multiplikativny stvorec velkosti 6 x 6",
                                                      stvorec,
                                                      ("rozpatie suctov: " + str(stavNajlepsi[0]),
                                                       "pocet roznych suctov: " + str(stavNajlepsi[1])))
                                break


#Algoritmus 4.5
def vrcholovoBimagickyGrafTest(n):
        G = []
        G = read_graph6("graph" + str(n) + "c.g6")

        for i in range(len(G)):
                susedia = []
                for _ in range(n):
                        susedia.append(set())
                        
                for hrana in G[i].edges():
                        susedia[hrana[0]].add(hrana[1])
                        susedia[hrana[1]].add(hrana[0])

                vyhovuje = True
                for v1,v2 in combinations((i for i in range(n)),r=2):
                        x = len(susedia[v1].difference(susedia[v2]))
                        y = len(susedia[v2].difference(susedia[v1]))
                        if (x*y == 0 and x+y > 0) or x == 1 or y == 1 or (x == 2 and y == 2):
                                vyhovuje = False
                                break
                        
                        if not vyhovuje:
                                break

                if vyhovuje:
                        vypisRiesenie("potencialne vrcholovo bimagicky graf s nasledovnymi hranami",
                                      G[i].edges(),
                                      {})


#Algoritmus 4.6
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
                      
#Algoritmus 4.7
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

#Algoritmus 4.8
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

def generujBimagickeObdlzniky(n,vsetkyNtice):
        for NticaSJednotkou,NticeBezJednotky in vsetkyNtice.items():
                m = len(NticaSJednotkou)
                
                for C in combinations(NticeBezJednotky,r=n-1):
                        zoznamNticBezJednotky = []
                        for Ntica in C:
                                zoznamNticBezJednotky += Ntica
                                
                        if rozne(zoznamNticBezJednotky + list(NticaSJednotkou)):
                                permutacie = []
                                for i in range(len(C)):
                                        permutacie.append([])
                                        for P in permutations(C[i]):
                                                permutacie[i].append(P)
                                                
                                pocetPermutacii = 1
                                for i in range(1,m+1):
                                        pocetPermutacii *= i
                                      
                                for vyberPermutacii in product([y for y in range(pocetPermutacii)],repeat=n-1):
                                        obdlznik = [NticaSJednotkou]
                                        for ii in range(len(vyberPermutacii)):
                                                obdlznik.append(permutacie[ii][vyberPermutacii[ii]])
                                        yield obdlznik
                                        

def transpozicia(obdlznik):
        m = len(obdlznik[0])
        n = len(obdlznik)
        novyObdlznik = []
        
        for stlpec in range(m):
                novyObdlznik.append([])
                for riadok in range(n):
                        novyObdlznik[stlpec].append(obdlznik[riadok][stlpec])
        return novyObdlznik


def overBimagickyObdlznik(obdlznik,ciastocneRiesenia):
        m = len(obdlznik[0])
        n = len(obdlznik)
        
        magickeSucty = set()
        bimagickeSucty = set()
        
        for riadok in range(m):
                magickySucet = 0
                bimagickySucet = 0
                for stlpec in range(n):
                        magickySucet += obdlznik[stlpec][riadok]
                        bimagickySucet += obdlznik[stlpec][riadok]**2
                magickeSucty.add(magickySucet)
                bimagickeSucty.add(bimagickySucet)
                
        if len(magickeSucty) == 1 and len(bimagickeSucty) < 3:
                if len(magickeSucty) == len(bimagickeSucty) == 1:
                        vypisRiesenie("bimagicky obdlznik velkosti " + str(m) + " x " + str(n),
                                      transpozicia(obdlznik),
                                      set())
                elif ciastocneRiesenia:
                        vypisRiesenie("ciastocny bimagicky obdlznik velkosti " + str(m) + " x " + str(n),
                                      transpozicia(obdlznik),
                                      {"bimagicke sucty v riadkoch: " + str(bimagickeSucty)})
                        

#Algoritmus 4.9
def bimagickyObdlznikSucet3xN(n,s,ciastocneRiesenia=True):
        trojice = dict()

        for trojicaSJednotkou,trojicaBezJednotky in generujBimagickeTrojiceSoSuctom(s):
                pridaj(trojicaSJednotkou,trojicaBezJednotky,trojice)

        for obdlznik in generujBimagickeObdlzniky(n,trojice):
                overBimagickyObdlznik(obdlznik,ciastocneRiesenia)
                

#Algoritmus 4.10
def bimagickyObdlznik3xN(n,h,ciastocneRiesenia=True):
        trojice = dict()

        for trojicaSJednotkou,trojicaBezJednotky in generujBimagickeTrojice(h):
                pridaj(trojicaSJednotkou,trojicaBezJednotky,trojice)

        for obdlznik in generujBimagickeObdlzniky(n,trojice):
                overBimagickyObdlznik(obdlznik,ciastocneRiesenia)


#Modifikacia algoritmu 4.9 pre vacsi obdlznik
def bimagickyObdlznikSucet4xN(n,s,ciastocneRiesenia=True):
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

        for obdlznik in generujBimagickeObdlzniky(n,stvorice):
                overBimagickyObdlznik(obdlznik,ciastocneRiesenia)
                

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

def generujMultiplikativneMagickeObdlzniky(n,vsetkyNtice):
        for sucinNtice,Ntice in vsetkyNtice.items():
                m = len(Ntice[0])
                for C in combinations(Ntice,r=n):
                        zoznamNtic = []
                        for Ntica in C:
                                zoznamNtic += Ntica
                                
                        if rozne(zoznamNtic):
                                permutacie = []
                                for i in range(1,len(C)):
                                        permutacie.append([])
                                        for P in permutations(C[i]):
                                                permutacie[i-1].append(P)

                                pocetPermutacii = 1
                                for i in range(1,m+1):
                                        pocetPermutacii *= i
                                                
                                for vyberPermutacii in product([y for y in range(pocetPermutacii)],repeat=n-1):
                                        obdlznik = [C[0]]
                                        for ii in range(len(vyberPermutacii)):
                                                obdlznik.append(permutacie[ii][vyberPermutacii[ii]])
                                        yield obdlznik


def overMultiplikativnyMagickyObdlznik(obdlznik,ciastocneRiesenia):
        m = len(obdlznik[0])
        n = len(obdlznik)
        
        magickeSucty = set()
        multiplikativneSuciny = set()
        
        for riadok in range(m):
                magickySucet = 0
                multiplikativnySucin = 1
                for stlpec in range(n):
                     magickySucet += obdlznik[stlpec][riadok]
                     multiplikativnySucin *= obdlznik[stlpec][riadok]
                magickeSucty.add(magickySucet)
                multiplikativneSuciny.add(multiplikativnySucin)
                
        if len(multiplikativneSuciny) == 1:
                if len(magickeSucty) == len(multiplikativneSuciny) == 1:
                        vypisRiesenie("multiplikativny magicky obdlznik velkosti " + str(m) + " x " + str(n),
                                      transpozicia(obdlznik),
                                      set())
                elif ciastocneRiesenia:
                        vypisRiesenie("ciastocny multiplikativny magicky obdlznik velkosti " + str(m) + " x " + str(n),
                                      transpozicia(obdlznik),
                                      {"magicke sucty v riadkoch: " + str(magickeSucty)})

#Algoritmus 4.11
def multiplikativnyMagickyObdlznikSucet3xN(n,sucet,ciastocneRiesenia=True):
        trojice = dict()

        for sucinTrojice,trojica in generujVyhovujuceMultiplikativneMagickeNticeSoSuctom(3,sucet):
                pridaj(sucinTrojice,trojica,trojice)

        for obdlznik in generujMultiplikativneMagickeObdlzniky(n,trojice):
                overMultiplikativnyMagickyObdlznik(obdlznik,ciastocneRiesenia)


#Modifikacia algoritmu 4.11 pre vacsi obdlznik
def multiplikativnyMagickyObdlznikSucet4xN(n,sucet,ciastocneRiesenia=True):
        stvorice = dict()

        for sucinStvorice,stvorica in generujVyhovujuceMultiplikativneMagickeNticeSoSuctom(4,sucet):
                pridaj(sucinStvorice,stvorica,stvorice)

        for obdlznik in generujMultiplikativneMagickeObdlzniky(n,stvorice):
                overMultiplikativnyMagickyObdlznik(obdlznik,ciastocneRiesenia)
