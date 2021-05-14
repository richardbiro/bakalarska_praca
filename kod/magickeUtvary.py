from itertools import combinations, permutations, product 
from math import gcd, inf, isqrt
from networkx import read_graph6
from random import choice, randint
from sympy import divisors, isprime


def jeDruhouMocninou(n):
        if n < 0:
                return False
        return isqrt(n) * isqrt(n) == n


def rozne(P):
        return len(P) == len(set(P))


def pridaj(kluc, prvok, asociativnePole):
        if kluc not in asociativnePole:
                asociativnePole[kluc] = [prvok]
        else:
                asociativnePole[kluc].append(prvok)


def bimagickySucet(pole):
        odpoved = 0
        for prvok in pole:
                odpoved += prvok * prvok
        return odpoved


def sucin(pole):
        odpoved = 1
        for prvok in pole:
                odpoved *= prvok
        return odpoved


def vypisRiesenie(nazovUtvaru, utvar, dodatocneUdaje):
        print(nazovUtvaru)
        
        dlzkaNajvacsiehoPrvku = 0
        for castUtvaru in utvar:
                for prvok in castUtvaru:
                       dlzkaNajvacsiehoPrvku = max(dlzkaNajvacsiehoPrvku, len(str(prvok)))
                       
        for castUtvaru in utvar:
                vypis = ""
                for index in range(len(castUtvaru)):
                        vypis += (dlzkaNajvacsiehoPrvku - len(str(castUtvaru[index])) + 1) * " " + str(castUtvaru[index])
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
def magickyStvorec3x3sPiatimiStvorcami(u1, v1, u2, v2):
        p = ((u1 * u1 + v1 * v1) * (u2 * u2 + 2 * u2 * v2 - v2 * v2))**2
        q = ((u1 * u1 + 2 * u1 * v1 - v1 * v1) * (u2 * u2 + v2 * v2))**2
        r = ((- u1 * u1 + 2 * u1 * v1 + v1 * v1) * (u2 * u2 + v2 * v2))**2
        s = ((u1 * u1 + v1 * v1) * (- u2 * u2 + 2 * u2 * v2 + v2 * v2))**2
        t = ((u1 * u1 + v1 * v1) * (u2 * u2 + v2 * v2))**2

        for stvorec in (([            p, 3 * t - p - q, q            ],
                         [3 * t - p - r,             t, 3 * t - q - s],
                         [            r, 3 * t - r - s, s            ]),
                        
                        ([2 * (r + s), 4 * p, 2 * (q + s)],
                         [      4 * q, 4 * t, 4 * r      ],
                         [2 * (p + r), 4 * s, 2 * (p + q)]),
                        
                        ([            p, q, 3 * t - p - q],
                         [    r + s - p, t, p + q - s    ],
                         [3 * t - r - s, r, s            ]),

                        ([            p, r, 3 * t - p - r],
                         [    q + s - p, t, p + r - s    ],
                         [3 * t - q - s, q, s            ])):
                
                if asponSedemStvorcov(stvorec):
                        vypisRiesenie("magicky stvorec velkosti 3 x 3 s aspon 7 druhymi mocninami",
                                      stvorec,
                                      set())               

#Algoritmus 4.2
def magickyStvorec3x3soSiestimiStvorcami(x):
        X = [1]
        for _ in range(10):
                X.append(X[-1] * x)

        x1 = 8 * X[8] - 49 * X[6] + 6 * X[4] - 16 * X[2] + 2
        x2 = 8 * X[8] - X[6] + 30 * X[4] - 40 * X[2] + 2
        x3 = 8 * X[8] - 25 * X[6] + 18 * X[4] - 28 * X[2] + 2

        stvorec1 = ([0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0])

        stvorec1[0][0] = (2 * X[5] + 4 * X[3] - 7 * X[1])**2
        stvorec1[0][1] = x1 * (X[2] - 2)
        stvorec1[0][2] = (5 * X[4] - 2 * X[2] + 2)**2
        stvorec1[1][0] = (X[4] + 8 * X[2] - 2)**2
        stvorec1[1][1] = (2 * X[5] - 2 * X[3] + 5 * X[1])**2
        stvorec1[1][2] = x2 * (X[2] - 2)
        stvorec1[2][0] = x3 * (X[2] - 2)
        stvorec1[2][1] = (7 * X[4] - 4 * X[2] - 2)**2
        stvorec1[2][2] = (2 * X[5] - 8 * X[3] - X[1])**2

        stvorec2 = ([0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0])
        
        stvorec2[0][0] = stvorec1[0][2]
        stvorec2[0][1] = stvorec1[0][0]
        stvorec2[0][2] = (4 * X[10] - 31 * X[8] + 76 * X[6] + 76 * X[4] - 31 * X[2] + 4)//2
        stvorec2[1][0] = stvorec1[2][2]
        stvorec2[1][1] = (4 * X[10] + 17 * X[8] + 4 * X[6] + 4 * X[4] + 17 * X[2] + 4)//2
        stvorec2[1][2] = stvorec1[2][1]
        stvorec2[2][0] = (4 * X[10] + 65 * X[8] - 68 * X[6] - 68 * X[4] + 65 * X[2] + 4)//2
        stvorec2[2][1] = stvorec1[1][0]
        stvorec2[2][2] = stvorec1[1][1]

        for stvorec in (stvorec1,stvorec2):
                if asponSedemStvorcov(stvorec):
                        vypisRiesenie("magicky stvorec velkosti 3 x 3 s aspon 7 druhymi mocninami",
                                      stvorec,
                                      set())


def generujTrojiceSBimagickymSuctom(dolnaHranica, hornaHranica):
        for a in range(0, isqrt((hornaHranica + 2)//3) + 1):
                for b in range(a + 1, isqrt((hornaHranica - a * a + 1)//2) + 1):
                        for c in range(max(b + 1, isqrt(max(dolnaHranica - a * a - b * b, 0))), isqrt(hornaHranica - a * a - b * b) + 1):
                                yield a, b, c


def najdiRieseniaFaktorizaciou(vyraz, p, q, r):
        for D1 in divisors(abs(vyraz), True):
                D2 = abs(vyraz)//D1
                if D1 % 2 == D2 % 2 and D1 >= D2:
                        sucetRovnic = ((D1+D2)//2, (-D1-D2)//2)
                        rozdielRovnic = ((D1-D2)//2, (D2-D1)//2)
                        
                        if vyraz < 0:
                                sucetRovnic, rozdielRovnic = rozdielRovnic, sucetRovnic
                                
                        for m in sucetRovnic:
                                for n in rozdielRovnic:
                                        s = - p - q - r + m
                                        if (s - p - q - r) % 2 == n % 2:
                                                x1 = (s - (p + q + r) + n)//2
                                                x2 = (s - (p + q + r) - n)//2
                                                x = min(x1, x2)
                                                y = max(x1, x2)
                                                yield x, y, s
                                                if m == 0:
                                                        break


def generujPrvkyBimagickehoStvorca(mozneProstrednePrvky):
        for s, zoznamUdajov in mozneProstrednePrvky.items():
                if len(set(zoznamUdajov)) >= 4:
                        U = []
                        for udaj in set(zoznamUdajov):
                                U.append(udaj)
                        U.sort()
                        
                        prvyUdajZacinajuci0 = -1
                        prvyUdajZacinajuci1 = -1
                        prvyUdajZacinajuci2 = -1
                        prvyUdajZacinajuci3 = -1
                        
                        for Ui in range(len(U)):
                                if U[Ui][0] == 0 and prvyUdajZacinajuci0 == -1:
                                        prvyUdajZacinajuci0 = Ui
                                elif U[Ui][0] == 1 and prvyUdajZacinajuci0 != -1 and prvyUdajZacinajuci1 == -1:
                                        prvyUdajZacinajuci1 = Ui
                                elif U[Ui][0] == 2 and prvyUdajZacinajuci1 != -1 and prvyUdajZacinajuci2 == -1:
                                        prvyUdajZacinajuci2 = Ui
                                elif U[Ui][0] == 3 and prvyUdajZacinajuci2 != -1 and prvyUdajZacinajuci3 == -1:
                                        prvyUdajZacinajuci3 = Ui

                        if -1 not in {prvyUdajZacinajuci0, prvyUdajZacinajuci1, prvyUdajZacinajuci2, prvyUdajZacinajuci3}:
                                for udajZacinajuci0 in range(prvyUdajZacinajuci0, prvyUdajZacinajuci1):
                                        for udajZacinajuci1 in range(prvyUdajZacinajuci1, prvyUdajZacinajuci2):
                                                for udajZacinajuci2 in range(prvyUdajZacinajuci2, prvyUdajZacinajuci3):
                                                        for udajZacinajuci3 in range(prvyUdajZacinajuci3, len(U)):
                                                                pi1, qi1, ri1, p1, q1, r1, X1, Y1 = U[udajZacinajuci0]
                                                                pi2, qi2, ri2, p2, q2, r2, X2, Y2 = U[udajZacinajuci1]
                                                                pi3, qi3, ri3, p3, q3, r3, X3, Y3 = U[udajZacinajuci2]
                                                                pi4 ,qi4, ri4, p4, q4, r4, X4, Y4 = U[udajZacinajuci3]
                                                                
                                                                if {pi1, pi2, pi3, pi4} == {qi1, qi2, qi3, qi4} == {ri1, ri2, ri3, ri4} == {0, 1, 2, 3}:
                                                                        prvkyStvorca = (p1, q1, r1, X1, Y1,
                                                                                        p2, q2, r2, X2, Y2,
                                                                                        p3, q3, r3, X3, Y3,
                                                                                        p4, q4, r4, X4, Y4, s)
                                                                        
                                                                        if rozne(prvkyStvorca):
                                                                                udaje1 = ((p1, q1, r1), (X1, Y1))
                                                                                udaje2 = ((p2, q2, r2), (X2, Y2))
                                                                                udaje4 = ((p3, q3, r3), (X3, Y3))
                                                                                udaje5 = ((p4, q4, r4), (X4, Y4))
                                                                                yield s, udaje1, udaje2, udaje4, udaje5


def generujBimagickeStvorce(s, udaje1, udaje2, udaje4, udaje5):
        for poradie in ((1, 0, 2), (0, 1, 2), (0, 2, 1)):
                for riadok1, riadok2, riadok4, riadok5 in ((udaje1, udaje2, udaje4, udaje5),
                                                           (udaje1, udaje2, udaje5, udaje4),
                                                           (udaje1, udaje4, udaje5, udaje2)):
                        
                        p1, q1, r1 = riadok1[0][poradie[0]], riadok1[0][poradie[1]], riadok1[0][poradie[2]]
                        p2, q2, r2 = riadok2[0][poradie[0]], riadok2[0][poradie[1]], riadok2[0][poradie[2]]
                        p3, q3, r3 = riadok4[0][poradie[0]], riadok4[0][poradie[1]], riadok4[0][poradie[2]]
                        p4, q4, r4 = riadok5[0][poradie[0]], riadok5[0][poradie[1]], riadok5[0][poradie[2]]
                        
                        for x1, y1 in permutations((riadok1[1][0], riadok1[1][1])):
                                for x2, y2 in permutations((riadok2[1][0], riadok2[1][1])):
                                        for x3, y3 in permutations((riadok4[1][0], riadok4[1][1])):
                                                for x4, y4 in permutations((riadok5[1][0], riadok5[1][1])):
                                                        
                                                        stvorec = ([p1, x1, q1, y1, r1],
                                                                   [x2, p2, q2, r2, y2],
                                                                   [-1, -1,  s, -1, -1],
                                                                   [x3, r3, q3, p3, y3],
                                                                   [r4, x4, q4, y4, p4])
                                                        
                                                        for stlpec in {0, 1, 3, 4}:
                                                                sucetStlpec = 0
                                                                for riadok in {0, 1, 3, 4}:
                                                                        sucetStlpec += stvorec[riadok][stlpec]
                                                                stvorec[2][stlpec] = s - sucetStlpec

                                                        yield stvorec


def overBimagickyStvorec(stvorec, K, zleSucty):
        if rozne(stvorec[0] + stvorec[1] + stvorec[2] + stvorec[3] + stvorec[4]):
                s = stvorec[2][2]
                nespravneBimagickeSucty = 5
                
                for stlpec in {0, 1, 3, 4}:
                        bimagickySucetStlpec = 0
                        for riadok in range(5):
                                bimagickySucetStlpec += stvorec[riadok][stlpec]**2
                                
                        if bimagickySucetStlpec == K + s * s:
                                nespravneBimagickeSucty -= 1

                bimagickySucetRiadok = 0
                riadok = 2
                for stlpec in range(5):
                        bimagickySucetRiadok += stvorec[riadok][stlpec]**2
                        
                if bimagickySucetRiadok == K + s * s:
                        nespravneBimagickeSucty -= 1
                                                        
                if nespravneBimagickeSucty <= zleSucty:
                        if nespravneBimagickeSucty == 0:
                                vypisRiesenie("bimagicky stvorec velkosti 5 x 5 so zapornymi prvkami",
                                              stvorec,
                                              set())
                        else:
                                vypisRiesenie("ciastocny bimagicky stvorec velkosti 5 x 5 so zapornymi prvkami",
                                              stvorec,
                                              {"pocet nespravnych bimagickych suctov: " + str(nespravneBimagickeSucty)})


#Algoritmus 4.3        
def bimagickyStvorec5x5(h, zleSucty=3):
        trojice = dict()
        dolnaHranica = 1
        hornaHranica = h

        for a, b, c in generujTrojiceSBimagickymSuctom(dolnaHranica, hornaHranica):
                bimagickySucetTrojice = a * a + b * b + c * c
                if dolnaHranica <= bimagickySucetTrojice and bimagickySucetTrojice < hornaHranica:
                        pridaj(bimagickySucetTrojice, (a, b, c), trojice)

        for bimagickySucetTrojice, vsetkyTrojice in trojice.items():
                if len(vsetkyTrojice) >= 3:
                        K = 4 * bimagickySucetTrojice
                        for i3, i2, i1 in combinations([i for i in range(len(vsetkyTrojice))], r = 3):
                                a, b, c = vsetkyTrojice[i1][0], vsetkyTrojice[i1][1], vsetkyTrojice[i1][2]
                                dd, ee, ff = vsetkyTrojice[i2][0], vsetkyTrojice[i2][1], vsetkyTrojice[i2][2]
                                gg, hh, ii = vsetkyTrojice[i3][0], vsetkyTrojice[i3][1], vsetkyTrojice[i3][2]
                                
                                for d, e, f in ((dd, ee, ff), (-dd, -ee, -ff)):
                                        for g,h,i in ((gg, hh, ii), (-gg, -hh, -ii)):
                                                mozneProstrednePrvky = dict()
                                                diagonala1 = (a + b - c, a - b + c, - a + b + c, - a - b - c)
                                                prostrednyStlpec = (d + e - f, d - e + f, - d + e + f, - d - e - f)
                                                diagonala2 = (g + h - i, g - h + i, - g + h + i, - g - h - i)
                                                pouzite = diagonala1 + prostrednyStlpec + diagonala2
                                                
                                                if rozne(pouzite):
                                                        for pi in range(4):
                                                                p = diagonala1[pi]
                                                                for qi in range(4):
                                                                        q = prostrednyStlpec[qi]
                                                                        for ri in range(4):
                                                                                r = diagonala2[ri]
                                                                                udaje = (pi, qi, ri, p, q , r)
                                                                                
                                                                                vyraz = 4 * (p * q + p * r + q * r + p**2 + q**2 + r**2) - 2 * K
                                                                                if vyraz != 0:
                                                                                        for x, y, s in najdiRieseniaFaktorizaciou(vyraz, p, q, r):
                                                                                                if rozne(pouzite + (x, y, s)):
                                                                                                        pridaj(s, udaje + (x, y), mozneProstrednePrvky)
                                                                                                        
                                                for s, udaje1, udaje2, udaje4, udaje5 in generujPrvkyBimagickehoStvorca(mozneProstrednePrvky):
                                                        for stvorec in generujBimagickeStvorce(s, udaje1, udaje2, udaje4, udaje5):
                                                                overBimagickyStvorec(stvorec, K, zleSucty)


def variacneRozpatieAPocetSuctov(stvorec):
        velkost = len(stvorec)
        prvky = []
        for riadok in stvorec:
                prvky += riadok
        if not rozne(prvky):
                return (inf, inf)
        
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
        return (max(sucty) - min(sucty), len(set(sucty)))


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


def generujStvorcoveVzorky(vsetkyVzorky, vybraneVzorky=[], vybraneIndexy=[]):
        pocetVsetkych = len(vsetkyVzorky)
        velkostVzorky = len(vsetkyVzorky[0])
        
        while True:
                pocetVybranych = len(vybraneVzorky)
                if suVzorkyDisjunktne(vybraneVzorky) and pocetVybranych < velkostVzorky:
                        if pocetVybranych == 0:
                                vybraneVzorky.append(vsetkyVzorky[0])
                                vybraneIndexy.append(0)
                        else:
                                vybraneVzorky.append(vsetkyVzorky[vybraneIndexy[-1] + 1])
                                vybraneIndexy.append(vybraneIndexy[-1] + 1)
                                
                        pocetVybranych += 1
                        
                else:
                        if suVzorkyDisjunktne(vybraneVzorky):
                                yield vybraneVzorky

                        while pocetVybranych > 0:
                                poslednyIndex = vybraneIndexy.pop()
                                vybraneVzorky.pop()
                                pocetVybranych -= 1
                                if poslednyIndex < pocetVsetkych - (velkostVzorky - pocetVybranych):
                                        vybraneVzorky.append(vsetkyVzorky[poslednyIndex + 1])
                                        vybraneIndexy.append(poslednyIndex + 1)
                                        pocetVybranych += 1
                                        break
                if pocetVybranych == 0:
                        break
                

def prenasobStvorec(hodnotyVzoriek, pouziteVzorky):
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
                        if vzorka[v] == velkost - v - 1: diagonala2 += 1
                        
                if diagonala1 == diagonala2 == 1:
                        yield vzorka


def vymenDveHodnoty(hodnotyVzoriek):
        for index1, index2 in permutations([i for i in range(len(hodnotyVzoriek))], 2):
                pamat1, pamat2 = hodnotyVzoriek[index1], hodnotyVzoriek[index2]
                hodnotyVzoriek[index1], hodnotyVzoriek[index2] = pamat2, pamat1
                yield hodnotyVzoriek
                hodnotyVzoriek[index1], hodnotyVzoriek[index2] = pamat1, pamat2


def zmenHodnotu(hodnotyVzoriek, h):
        for index in range(len(hodnotyVzoriek)):
                pamat = hodnotyVzoriek[index]
                
                for novaHodnota in range(1, h + 1):
                        hodnotyVzoriek[index] = novaHodnota
                        yield hodnotyVzoriek
                        
                hodnotyVzoriek[index] = pamat


def pripocitajRozsah(hodnotyVzoriek, velkost, rozsah):
        for indexStvorcovejVzorky in range(0, len(hodnotyVzoriek), velkost):
                for pripocitat in product([i for i in range(-rozsah, rozsah + 1)], repeat=velkost):
                        pamat = []
                        for hodnotaVzorky in hodnotyVzoriek:
                                pamat.append(hodnotaVzorky)

                        for i in range(len(pripocitat)):
                                hodnotyVzoriek[indexStvorcovejVzorky + i] += pripocitat[i]

                        yield hodnotyVzoriek
                                
                        for i in range(len(pripocitat)):
                                hodnotyVzoriek[indexStvorcovejVzorky + i] = pamat[indexStvorcovejVzorky + i]  


#Algoritmus 4.4                                
def multiplikativnyMagickyStvorec6x6(p,h):
        velkost = 6
        vzorky = []
        stvorcoveVzorky = []
        
        for vzorka in generujVzorky(velkost):
                vzorky.append(tuple(vzorka))
        
        for stvorcovaVzorka in generujStvorcoveVzorky(vzorky):
                stvorcoveVzorky.append(tuple(stvorcovaVzorka))

        rozsah = 1        
        stavNajlepsi = (inf, inf)
        optimum = (0, 0)
        
        while True:
                hodnotyVzoriek = []
                pouziteVzorky = []
                for i in range(p):
                        hodnotyVzoriek += [randint(1, h) for _ in range(velkost)]
                        pouziteVzorky += choice(stvorcoveVzorky)
                        
                stvorec = prenasobStvorec(hodnotyVzoriek,pouziteVzorky)
                stav = variacneRozpatieAPocetSuctov(stvorec)

                while stav > optimum:
                        stavZaciatok = stav
                        
                        if stav != (inf, inf):
                                for operacia in {vymenDveHodnoty(hodnotyVzoriek),
                                                 zmenHodnotu(hodnotyVzoriek, h),
                                                 pripocitajRozsah(hodnotyVzoriek, velkost, rozsah)}:
                                        
                                        for hodnotyVzoriekNove in operacia:
                                                stvorecNovy = prenasobStvorec(hodnotyVzoriekNove, pouziteVzorky)
                                                stavNovy = variacneRozpatieAPocetSuctov(stvorecNovy)
                                                if stavNovy < stav:
                                                        stav = stavNovy
                                                        stvorec = [riadok for riadok in stvorecNovy]
                                                        hodnotyVzoriek = [hodnota for hodnota in hodnotyVzoriekNove]

                        if stav == stavZaciatok:
                                if stav < stavNajlepsi:
                                        stavNajlepsi = stav
                                        vypisRiesenie("multiplikativny stvorec velkosti 6 x 6",
                                                      stvorec,
                                                      ("variacne rozpatie magickych suctov: " + str(stavNajlepsi[0]),
                                                       "pocet roznych magickych suctov: " + str(stavNajlepsi[1])))
                                break


#Algoritmus 4.5
def vrcholovoBimagickyGrafTest(n):
        G = []
        G = read_graph6("zoznam grafov/graph" + str(n) + "c.g6")

        for i in range(len(G)):
                susedia = []
                for _ in range(n):
                        susedia.append(set())
                        
                for hrana in G[i].edges():
                        susedia[hrana[0]].add(hrana[1])
                        susedia[hrana[1]].add(hrana[0])

                vyhovuje = True
                
                for v1,v2 in combinations((i for i in range(n)), r = 2):
                        x = len(susedia[v1].difference(susedia[v2]))
                        y = len(susedia[v2].difference(susedia[v1]))
                        if (x * y == 0 and x + y > 0) or x == 1 or y == 1 or (x == 2 and y == 2):
                                vyhovuje = False
                                break
                        
                        if not vyhovuje:
                                break

                if vyhovuje:
                        vypisRiesenie("potencialne vrcholovo bimagicky graf s nasledovnymi hranami",
                                      G[i].edges(),
                                      {"pocet vrcholov: " + str(n)})


#Algoritmus 4.6
def vrcholovoBimagickyKompletnyGraf(i, j):
        if i > j:
                H = vrcholovoBimagickyKompletnyGraf(j, i)
                if H is not None:
                        return H[::-1]
                return
        
        if i <= 1 or i == j == 2:
                return
        
        if i == 2:
                H1 = [j * (j - 1)//2 + 1, j * (j - 1) * (3 * j * j - 7 * j + 14)//24]
                H2 = [k for k in range(1, j)]
                H2.append(j * (j - 1) * (3 * j * j - 7 * j + 14)//24 + 1)
                return (H1,H2)
        
        if i == 3:
                H1 = [1, j * (j + 1)//2 - 1, j * (j + 1) * (3 * j * j - j - 14)//24 + 1]
                H2 = [k for k in range(2, j + 1)]
                H2.append(j * (j + 1) * (3 * j * j - j - 14)//24 + 2)
                return (H1,H2)
        
        if (i, j) == (4, 4):
                return ([1, 4, 6, 7], [2, 3, 5, 8])
        if (i, j) == (4, 5):
                return ([2, 12, 13, 15], [1, 4, 8, 10, 19])
        
        H = vrcholovoBimagickyKompletnyGraf(i - 2, j - 3)
        m = max(max(H[0]), max(H[1])) + 1
        H = (H[0] + [4 * m, 5 * m], H[1] + [m, 2 * m, 6 * m])
        return H

               
#Algoritmus 4.7
def vrcholovoSuperbimagickyKompletnyGraf(n):
        if n < 7:
                return
        if n % 4 in {1, 2}:
                return
        if n == 7:
                return ([1, 2, 4, 7], [3, 5, 6])
        if n == 8:
                return ([1, 4, 6, 7], [2, 3, 5, 8])
        if n == 11:
                return ([1, 3, 4, 5, 9, 11], [2, 6, 7, 8, 10])
        if n == 12:
                return ([1, 3, 7, 8, 9, 11], [2, 4, 5, 6, 10, 12])
        
        H = vrcholovoSuperbimagickyKompletnyGraf(n - 8)
        for x in range(1, 9):
                if x in {1, 4, 6, 7}:
                        H[0].append(n - 8 + x)
                else:
                        H[1].append(n - 8 + x)
        return H


#Algoritmus 4.8
def vrcholovoMultiplikativnyMagickyKompletnyGraf(i, j):
        if i > j:
                H = vrcholovoMultiplikativnyMagickyKompletnyGraf(j, i)
                if H is not None:
                        return H[::-1]
                return
        
        if i <= 1 or i == j == 2:
                return

        if (i, j) == (2, 3):
                return ([5, 12], [1, 6, 10])
        if (i, j) == (2, 4):
                return ([9, 16], [1, 2, 4, 18])
        
        if i == 2:
                faktorial = 1
                for f in range(1, j): faktorial *= f
                H1 = [faktorial + 1, faktorial * (faktorial + 1 - j * (j - 1)//2)]
                H2 = [k for k in range(1,j)]
                H2.append((faktorial + 1) * (faktorial + 1 - j * (j - 1)//2))
                return (H1, H2)
        
        if i == 3:
                faktorial = 1
                for f in range(1, j + 1): faktorial *= f
                H1 = [1, faktorial + 1, faktorial * (faktorial + 3 - j * (j + 1)//2)]
                H2 = [k for k in range(2, j + 1)]
                H2.append((faktorial + 1) * (faktorial + 3 - j * (j + 1)//2))
                return (H1, H2)
        
        if (i,j) == (4, 4):
                return ([1, 5, 6, 12], [2, 3, 4, 15])
        if (i,j) == (4, 5):
                return ([2, 10, 20, 27], [1, 3, 6, 24, 25])
        
        H = vrcholovoMultiplikativnyMagickyKompletnyGraf(i - 2, j - 3)
        x = max(max(H[0]), max(H[1])) + 1
        y = max(max(H[0]), max(H[1])) + 2
        H = (H[0] + [2 * x * y, 2 * x * y - x - y], H[1] + [2 * (2 * x * y - x - y), x, y])
        return H


def najdiJednotkovuTrojicu(a, b, c):
        s = a + b + c
        t = a * a + b * b + c * c
        vyraz = 2 * t - (s - 1) * (s - 1) - 2
        if jeDruhouMocninou(vyraz):
                n = isqrt(vyraz)
                if s % 2 != n % 2:
                        x1 = (s - 1 + n)//2
                        x2 = (s - 1 - n)//2
                        if min(x1, x2) > 1:
                                if rozne({1, x1, x2, a, b, c}):
                                        return ((1, min(x1, x2), max(x1, x2)), (a, b, c))


def generujBimagickeTrojice(hranica):
        for a, b, c in combinations([i for i in range(2, hranica + 1)],r = 3):
                vysledok = najdiJednotkovuTrojicu(a, b, c)
                if vysledok is not None:
                        yield vysledok


def generujTrojiceSoSuctom(sucet):
        for a in range(2, sucet//3 + 1):
                for b in range(a + 1,(sucet - a)//2 + 1):
                        c = sucet - a - b
                        yield (a, b, c)

                        
def generujBimagickeTrojiceSoSuctom(sucet):
        for a, b, c in generujTrojiceSoSuctom(sucet):
                vysledok = najdiJednotkovuTrojicu(a, b, c)
                if vysledok is not None:
                        yield vysledok


def generujStvoriceSoSuctom(sucet):
        for a in range(2, sucet//4 + 1):
                for b in range(a + 1, (sucet - a)//3 + 1):
                        for c in range(b + 1,(sucet - a - b)//2 + 1):
                                d = sucet - a - b - c
                                yield (a, b, c, d)


def generujNticeSoSuctom(n,sucet):
        if n * (n + 1)//2 > sucet:
                return
        
        Ntica = [i for i in range(1, n)]
        Ntica.append(sucet - n * (n - 1)//2)
        
        while True:
                if Ntica[n - 1] - Ntica[n - 2] > 1:
                        Ntica[n - 1] -= 1
                        Ntica[n - 2] += 1
                else:
                        posledny = n - 2
                        while True:
                                if posledny <= 0:
                                        return
                                
                                if Ntica[posledny] - Ntica[posledny - 1] > 1:
                                        Ntica[posledny - 1] += 1
                                        for dalsi in range(posledny, n - 1):
                                                Ntica[dalsi] = Ntica[dalsi - 1] + 1
                                        sucetOkremPosledneho = sum(Ntica) - Ntica[n - 1]
                                        Ntica[n - 1] = sucet - sucetOkremPosledneho
                                        break
                                
                                else:
                                        posledny -= 1 
                yield Ntica


def generujBimagickeObdlzniky(n, vsetkyNtice):
        for NticaSJednotkou, NticeBezJednotky in vsetkyNtice.items():
                m = len(NticaSJednotkou)
                
                for C in combinations(NticeBezJednotky, r = n - 1):
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
                                for i in range(1, m + 1):
                                        pocetPermutacii *= i
                                      
                                for vyberPermutacii in product([y for y in range(pocetPermutacii)], repeat = n - 1):
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
def bimagickyObdlznik3xNSoSuctom(n, s, ciastocneRiesenia=True):
        trojice = dict()

        for trojicaSJednotkou, trojicaBezJednotky in generujBimagickeTrojiceSoSuctom(s):
                pridaj(trojicaSJednotkou, trojicaBezJednotky, trojice)

        for obdlznik in generujBimagickeObdlzniky(n, trojice):
                overBimagickyObdlznik(obdlznik, ciastocneRiesenia)
                

#Algoritmus 4.10
def bimagickyObdlznik3xN(n, h, ciastocneRiesenia=True):
        trojice = dict()

        for trojicaSJednotkou, trojicaBezJednotky in generujBimagickeTrojice(h):
                pridaj(trojicaSJednotkou, trojicaBezJednotky, trojice)

        for obdlznik in generujBimagickeObdlzniky(n, trojice):
                overBimagickyObdlznik(obdlznik, ciastocneRiesenia)


#Modifikacia algoritmu 4.9 pre vacsi obdlznik
def bimagickyObdlznik4xNSoSuctom(n, s, ciastocneRiesenia=True):
        trojice = dict()
        
        for trojica in generujTrojiceSoSuctom(s - 1):
                magickySucetTrojice = sum(trojica)
                bimagickySucetTrojice = bimagickySucet(trojica)
                pridaj((magickySucetTrojice, bimagickySucetTrojice), trojica, trojice)

        stvorice = dict()

        for a, b, c, d in generujStvoriceSoSuctom(s):
                t = a * a + b * b + c * c + d * d
                if (s - 1, t - 1) in trojice:
                        for x, y, z in trojice[(s - 1, t - 1)]:
                                if rozne({a, b, c, d, 1, x, y, z}):
                                        pridaj((1, x, y, z), (a, b, c, d), stvorice)

        for obdlznik in generujBimagickeObdlzniky(n, stvorice):
                overBimagickyObdlznik(obdlznik, ciastocneRiesenia)
                

def generujVyhovujuceMultiplikativneMagickeNticeSoSuctom(n, sucet):
        vyhovuju = set()
        for cislo in range(1, sucet + 1):
                if not isprime(cislo) or n * cislo <= sucet:
                        vyhovuju.add(cislo)
                        
        for Ntica in generujNticeSoSuctom(n, sucet):
                vyhovuje = True
                for prvok in Ntica:
                        if prvok not in vyhovuju:
                                vyhovuje = False
                                break
                if vyhovuje:
                        p = sucin(Ntica)
                        yield (p, tuple(Ntica))


def generujMultiplikativneMagickeObdlzniky(n, vsetkyNtice):
        for sucinNtice, Ntice in vsetkyNtice.items():
                m = len(Ntice[0])
                for C in combinations(Ntice, r = n):
                        zoznamNtic = []
                        for Ntica in C:
                                zoznamNtic += Ntica
                                
                        if rozne(zoznamNtic):
                                permutacie = []
                                for i in range(1, len(C)):
                                        permutacie.append([])
                                        for P in permutations(C[i]):
                                                permutacie[i - 1].append(P)

                                pocetPermutacii = 1
                                for i in range(1,m + 1):
                                        pocetPermutacii *= i
                                                
                                for vyberPermutacii in product([y for y in range(pocetPermutacii)],repeat = n - 1):
                                        obdlznik = [C[0]]
                                        for ii in range(len(vyberPermutacii)):
                                                obdlznik.append(permutacie[ii][vyberPermutacii[ii]])
                                        yield obdlznik


def overMultiplikativnyMagickyObdlznik(obdlznik, ciastocneRiesenia):
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
def multiplikativnyMagickyObdlznik3xNSoSuctom(n, sucet, ciastocneRiesenia = True):
        trojice = dict()

        for sucinTrojice,trojica in generujVyhovujuceMultiplikativneMagickeNticeSoSuctom(3, sucet):
                pridaj(sucinTrojice, trojica, trojice)

        for obdlznik in generujMultiplikativneMagickeObdlzniky(n, trojice):
                overMultiplikativnyMagickyObdlznik(obdlznik, ciastocneRiesenia)


#Modifikacia algoritmu 4.11 pre vacsi obdlznik
def multiplikativnyMagickyObdlznik4xNSoSuctom(n, sucet, ciastocneRiesenia = True):
        stvorice = dict()

        for sucinStvorice, stvorica in generujVyhovujuceMultiplikativneMagickeNticeSoSuctom(4, sucet):
                pridaj(sucinStvorice, stvorica, stvorice)

        for obdlznik in generujMultiplikativneMagickeObdlzniky(n, stvorice):
                overMultiplikativnyMagickyObdlznik(obdlznik, ciastocneRiesenia)
