from networkx import read_graph6, enumerate_all_cliques
from itertools import product, permutations, combinations
from math import gcd, isqrt
from sympy import factorint, isprime, divisors
from datetime import datetime
import random

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

def sedem_stvorcov(X):
        if len(set(X[0] + X[1] + X[2])) == 9:
                pocet = 0
                for i in X:
                     for j in i:
                             if stvorec(j): pocet += 1
                if pocet >= 7: return True
                return False
        return False

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

        v = [[p,x1,q],
             [x2,t,x3],
             [r,x4,s]]

        if sedem_stvorcov(v): print(v)
        
        w = [[p,q,y1],
             [y2,t,y3],
             [y4,r,s]]

        if sedem_stvorcov(w): print(w)

        ww = [[p,r,x2],
              [q+s-p,t,p+r-s],
              [x3,q,s]]

        if sedem_stvorcov(ww): print(ww)
        
        if p%2 == q%2 == r%2 == s%2:
                www = [[(r+s)//2,p,(q+s)//2],
                      [q,t,r],
                      [(p+r)//2,s,(p+q)//2]]
                if sedem_stvorcov(www): print(www)
                

def magicky_stvorec_3x3_6_stvorcov(x):
        X = [1]
        for _ in range(10): X.append(X[-1]*x)

        v = [[0,0,0],
             [0,0,0],
             [0,0,0]]
        v[0][0] = (2*X[5] + 4*X[3] - 7*X[1])**2
        v[0][1] = (X[2] - 2)*(8*X[2] - 1)*(X[6] - 6*X[4] - 2)
        v[0][2] = (5*X[4] - 2*X[2] + 2)**2
        v[1][0] = (X[4] + 8*X[2] - 2)**2
        v[1][1] = (2*X[5] - 2*X[3] + 5*X[1])**2
        v[1][2] = (X[2] - 2)*(8*X[8] - X[6] + 30*X[4] - 40*X[2] + 2)
        v[2][0] = (X[2] - 2)*(8*X[8] - 25*X[6] + 18*X[4] - 28*X[2] + 2)
        v[2][1] = (7*X[4] - 4*X[2] - 2)**2
        v[2][2] = (2*X[5] - 8*X[3] - X[1])**2
        if sedem_stvorcov(v): print(v)

        w = [[0,0,0],
             [0,0,0],
             [0,0,0]]
        w[0][0] = v[1][0]
        w[0][1] = v[0][0]
        w[0][2] = (4*X[10] - 31*X[8] + 76*X[6] + 76*X[4] - 31*X[2] + 4)//2
        w[1][0] = v[2][2]
        w[1][1] = (4*X[10] + 17*X[8] + 4*X[6] + 4*X[4] + 17*X[2] + 4)//2
        w[1][2] = v[2][1]
        w[2][0] = (4*X[10] + 65*X[8] - 68*X[6] - 68*X[4] + 65*X[2] + 4)//2
        w[2][1] = v[0][2]
        w[2][2] = v[1][1]
        if sedem_stvorcov(w): print(w)
        
def bimagicky_stvorec_5x5(h):
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
                        for i1 in range(2,len(j)):
                                a = j[i1][0]
                                b = j[i1][1]
                                c = j[i1][2]
                                for i2 in range(1,i1):
                                        dd = j[i2][0]
                                        ee = j[i2][1]
                                        ff = j[i2][2]
                                        for i3 in range(i2):
                                                gg = j[i3][0]
                                                hh = j[i3][1]
                                                ii = j[i3][2]
                                                for d,e,f in [[dd,ee,ff],[-dd,-ee,-ff]]:
                                                        for g,h,i in [[gg,hh,ii],[-gg,-hh,-ii]]:
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
                                                                                        n = 2*((V1+V2+V3)*(V1+V2+V3) + V1*V1 + V2*V2 + V3*V3 - K)
                                                                                        output = str(Vabc) + " " + str(Vdef) + " " + str(Vghi) + " " + str(V1) + " " + str(V2) + " " + str(V3) + " "

                                                                                        if n == 0 or n%4 == 2: pass
                                                                                        elif n > 0:
                                                                                                for D1 in divisors(n,True):
                                                                                                        D2 = n//D1
                                                                                                        if D1%2 == D2%2 and D1 >= D2:
                                                                                                                for M in [(D1+D2)//2,(-D1-D2)//2]:
                                                                                                                        for L in [(D1-D2)//2,(D2-D1)//2]:
                                                                                                                                S = -V1-V2-V3+M
                                                                                                                                if (S-V1-V2-V3)%2 == L%2:
                                                                                                                                        x1 = min((S-V1-V2-V3+L)//2,(S-V1-V2-V3-L)//2)
                                                                                                                                        x2 = max((S-V1-V2-V3+L)//2,(S-V1-V2-V3-L)//2)
                                                                                                                                        if S not in hodnoty: hodnoty[S] = {output + str(x1) + " " + str(x2)}
                                                                                                                                        else: hodnoty[S].add(output + str(x1) + " " + str(x2))
                                                                                                                                        if M == 0: break
                                                                                        else:
                                                                                                for D1 in divisors(-n,True):
                                                                                                        D2 = (-n)//D1
                                                                                                        if D1%2 == D2%2 and D1 >= D2:
                                                                                                                for L in [(D1+D2)//2,(-D1-D2)//2]:
                                                                                                                        for M in [(D1-D2)//2,(D2-D1)//2]:
                                                                                                                                S = -V1-V2-V3+M
                                                                                                                                if (S-V1-V2-V3)%2 == L%2:
                                                                                                                                        x1 = min((S-V1-V2-V3+L)//2,(S-V1-V2-V3-L)//2)
                                                                                                                                        x2 = max((S-V1-V2-V3+L)//2,(S-V1-V2-V3-L)//2)
                                                                                                                                        if S not in hodnoty: hodnoty[S] = {output + str(x1) + " " + str(x2)}
                                                                                                                                        else: hodnoty[S].add(output + str(x1) + " " + str(x2))
                                                                                                                                        if M == 0: break
                                                                for s,l in hodnoty.items():
                                                                        if len(l) >= 4:
                                                                                
                                                                                H = []
                                                                                for ll in l: H.append(list(map(int,ll.split(' '))))
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
                                                                                                                                                                                                stvorec[2][0] = s - stvorec[0][0] - stvorec[1][0] - stvorec[3][0] - stvorec[4][0]
                                                                                                                                                                                                stvorec[2][1] = s - stvorec[0][1] - stvorec[1][1] - stvorec[3][1] - stvorec[4][1]
                                                                                                                                                                                                stvorec[2][3] = s - stvorec[0][3] - stvorec[1][3] - stvorec[3][3] - stvorec[4][3]
                                                                                                                                                                                                stvorec[2][4] = s - stvorec[0][4] - stvorec[1][4] - stvorec[3][4] - stvorec[4][4]
                                                                                                                                                                                                if len(pouzite.union({stvorec[2][0],stvorec[2][1],stvorec[2][3],stvorec[2][4]})) == 25:
                                                                                                                                                                                                        nespravne = 5
                                                                                                                                                                                                        if stvorec[0][0]*stvorec[0][0] + stvorec[1][0]*stvorec[1][0] + stvorec[2][0]*stvorec[2][0] + stvorec[3][0]*stvorec[3][0] + stvorec[4][0]*stvorec[4][0] == s2: nespravne -= 1
                                                                                                                                                                                                        if stvorec[0][1]*stvorec[0][1] + stvorec[1][1]*stvorec[1][1] + stvorec[2][1]*stvorec[2][1] + stvorec[3][1]*stvorec[3][1] + stvorec[4][1]*stvorec[4][1] == s2: nespravne -= 1
                                                                                                                                                                                                        if stvorec[0][3]*stvorec[0][3] + stvorec[1][3]*stvorec[1][3] + stvorec[2][3]*stvorec[2][3] + stvorec[3][3]*stvorec[3][3] + stvorec[4][3]*stvorec[4][3] == s2: nespravne -= 1
                                                                                                                                                                                                        if stvorec[0][4]*stvorec[0][4] + stvorec[1][4]*stvorec[1][4] + stvorec[2][4]*stvorec[2][4] + stvorec[3][4]*stvorec[3][4] + stvorec[4][4]*stvorec[4][4] == s2: nespravne -= 1
                                                                                                                                                                                                        if stvorec[2][0]*stvorec[2][0] + stvorec[2][1]*stvorec[2][1] + stvorec[2][2]*stvorec[2][2] + stvorec[2][3]*stvorec[2][3] + stvorec[2][4]*stvorec[2][4] == s2: nespravne -= 1
                                                                                                                                                                                                        if nespravne <= 3:
                                                                                                                                                                                                                print("S =",s)
                                                                                                                                                                                                                print("S2 =",s2)
                                                                                                                                                                                                                print(nespravne,"nespravnych")
                                                                                                                                                                                                                for ss in stvorec: print(ss)
                                                                                                                                                                                                                print()



def func(X):
    if len(set(X[0] + X[1] + X[2] + X[3] + X[4] + X[5])) < 36: return [10**40,10**40]

    sucty = []
    for i in range(6):
        sucty.append(X[i][0] + X[i][1] + X[i][2] + X[i][3] + X[i][4] + X[i][5])
        sucty.append(X[0][i] + X[1][i] + X[2][i] + X[3][i] + X[4][i] + X[5][i])

    sucty.append(X[0][0] + X[1][1] + X[2][2] + X[3][3] + X[4][4] + X[5][5])
    sucty.append(X[0][5] + X[1][4] + X[2][3] + X[3][2] + X[4][1] + X[5][0])
    return [max(sucty) - min(sucty),len(set(sucty))]


def vzorka(hodnoty,pouzite):
    A = []
    for i in range(6):
        A.append([])
        for j in range(6):
            A[i].append(1)
    
    for index in range(len(hodnoty)):
        for i in range(6): A[i][pouzite[index][i]] *= hodnoty[index]

    return A

                                
def multiplikativny_magicky_stvorec(h):
        vzorky = []
        hodnoty = []

        for a in range(6):
            for b in range(6):
                for c in range(6):
                    for d in range(6):
                        for e in range(6):
                            for f in range(6):
                                if len(set([a,b,c,d,e,f])) == 6:
                                    P = []
                                    if a == 0: P.append(a)
                                    if b == 1: P.append(b)
                                    if c == 2: P.append(c)
                                    if d == 3: P.append(d)
                                    if e == 4: P.append(e)
                                    if f == 5: P.append(f)
                                    if len(P) == 1:
                                        P = []
                                        if a == 5: P.append(a)
                                        if b == 4: P.append(b)
                                        if c == 3: P.append(c)
                                        if d == 2: P.append(d)
                                        if e == 1: P.append(e)
                                        if f == 0: P.append(f)
                                        if len(P) == 1: vzorky.append([a,b,c,d,e,f])

        extrem = 10**40
        minimum = [extrem,extrem]
        tolerancia = [0,1]
        parametre = len(vzorky)
        
        stvorce = []
        it = len(vzorky)//6
        for i in range(0,it):
            stvorec = [vzorky[i]]
            ok = True
            for j in range(it,2*it):
                stvorec = [vzorky[i],vzorky[j]]
                for jj in range(1):
                    for jjj in range(6):
                        if stvorec[jj][jjj] == vzorky[j][jjj]: ok = False
                if not ok:
                    ok = True
                else:
                    for k in range(2*it,3*it):
                        stvorec = [vzorky[i],vzorky[j],vzorky[k]]
                        for kk in range(2):
                            for kkk in range(6):
                                if stvorec[kk][kkk] == vzorky[k][kkk]: ok = False
                        if not ok:
                            ok = True
                        else:
                            for l in range(3*it,4*it):
                                stvorec = [vzorky[i],vzorky[j],vzorky[k],vzorky[l]]
                                for ll in range(3):
                                    for lll in range(6):
                                        if stvorec[ll][lll] == vzorky[l][lll]: ok = False
                                if not ok:
                                    ok = True
                                else:
                                    for m in range(4*it,5*it):
                                        stvorec = [vzorky[i],vzorky[j],vzorky[k],vzorky[l],vzorky[m]]
                                        for mm in range(4):
                                            for mmm in range(6):
                                                if stvorec[mm][mmm] == vzorky[m][mmm]: ok = False
                                        if not ok:
                                            ok = True
                                        else:
                                            for n in range(5*it,6*it):
                                                stvorec = [vzorky[i],vzorky[j],vzorky[k],vzorky[l],vzorky[m],vzorky[n]]
                                                for nn in range(5):
                                                    for nnn in range(6):
                                                        if stvorec[nn][nnn] == vzorky[n][nnn]: ok = False
                                                if not ok:
                                                    ok = True
                                                else:
                                                    stvorce.append(stvorec)

        parametre = 3*6
        rozsah = 3
        skok = 1
        while True:
                hodnoty = []
                pouzite = []
                for i in range(parametre//6):
                        d = random.randint(1,h)
                        moznosti = [random.randint(1,h) for i in range(d,d+6)]
                        random.shuffle(moznosti)
                        hodnoty += moznosti
                        pouzite += random.choice(stvorce)
                X = vzorka(hodnoty,pouzite)
                vvvX = []
                for i in range(parametre): vvvX.append(hodnoty[i])
                v = func(X)

                while v > [0,0]:
                        vvv = v
                        if v != [extrem,extrem]:
                                for index1 in range(parametre):
                                        for index2 in range(index1+1,parametre):
                                                pamat1 = hodnoty[index1]
                                                pamat2 = hodnoty[index2]
                                                hodnoty[index1] = pamat2
                                                hodnoty[index2] = pamat1

                                                Y = vzorka(hodnoty,pouzite)

                                                if func(Y) < vvv:
                                                        vvv = func(Y)
                                                        vvvX = []
                                                        for i in range(parametre): vvvX.append(hodnoty[i])
                                                        
                                                hodnoty[index1] = pamat1
                                                hodnoty[index2] = pamat2
                                if vvv == v:
                                        for index in range(parametre):
                                                for cislo in range(1,h+1):
                                                        pamat = hodnoty[index]
                                                        hodnoty[index] = cislo
                                                        Y = vzorka(hodnoty,pouzite)
                                                        if func(Y) < vvv:
                                                                vvv = func(Y)
                                                                vvvX = []
                                                                for i in range(parametre): vvvX.append(hodnoty[i])
                                                        hodnoty[index] = pamat
                                if vvv == v:
                                        for x in range(0,parametre,6):
                                                for i in range(rozsah**6):
                                                        pamat = []
                                                        for j in range(parametre): pamat.append(hodnoty[j])
                                                        k = i
                                                        l = 0
                                                            
                                                        while k > 0:
                                                                hodnoty[x+l] += (k%rozsah) - rozsah//2
                                                                if hodnoty[x+l] < 1: hodnoty[x+l] = 1
                                                                k = k//rozsah
                                                                l += 1
                                                        Y = vzorka(hodnoty,pouzite)
                                                        if func(Y) < vvv:
                                                                vvv = func(Y)
                                                                vvvX = []
                                                                for j in range(parametre): vvvX.append(hodnoty[j])
                                                        for j in range(parametre): hodnoty[j] = pamat[j]   
                                

                                if vvv == v and v[0] < 100:
                                        for index1 in range(parametre):
                                                for cislo1 in range(1,h+1):
                                                        for index2 in range(index1+1,parametre):
                                                                for cislo2 in range(1,h+1):
                                                                        pamat1 = hodnoty[index1]
                                                                        pamat2 = hodnoty[index2]
                                                                        hodnoty[index1] = cislo1
                                                                        hodnoty[index2] = cislo2

                                                                        Y = vzorka(hodnoty,pouzite)

                                                                        if func(Y) < vvv:
                                                                                vvv = func(Y)
                                                                                vvvX = []
                                                                                for i in range(parametre): vvvX.append(hodnoty[i])
                                                                                
                                                                        hodnoty[index1] = pamat1
                                                                        hodnoty[index2] = pamat2
                        if vvv == v:
                                if v < minimum:
                                        minimum = v
                                        print("pocet roznych suctov:",vvv[1])
                                        print("rozpatie suctov:",vvv[0])
                                        XX = vzorka(hodnoty,pouzite)
                                        for i in XX: print(i)
                                        print()
                                break
                        else:
                                hodnoty = []
                                for i in range(len(vvvX)): hodnoty.append(vvvX[i])
                                vzorka(hodnoty,pouzite)
                                v = vvv
                                

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

def bimagicky_obdlznik_3xN_nulovy(h):
        dvojice = dict()
        for a in range(h+1):
                for b in range(-a+1,a,a%2+1):
                        t = a*a + b*b + (-a-b)*(-a-b)
                        if t not in dvojice: dvojice[t] = [(a,b)]
                        else: dvojice[t].append((a,b))
        maximum = 0
        for i,j in dvojice.items(): maximum = max(maximum,len(j))
        for n in range(4,maximum):
                for sucet,j in dvojice.items():
                        for C in combinations(j,r=n):
                                prvky = []
                                for c in C: prvky += [c[0],c[1],-c[0]-c[1]]
                                if len(set(prvky)) == 3*n:
                                        moznosti = []
                                        for i in range(1,len(C)):
                                                moznosti.append([])
                                                for P in permutations((C[i][0],C[i][1],-C[i][0]-C[i][1])): moznosti[i-1].append(P)
                                        for PP in product([y for y in range(6)],repeat=n-1):
                                                obdlznik = [[C[0][0],C[0][1],-C[0][0]-C[0][1]]]
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


def bimagicky_obdlznik_sucet_4xN(h):
        trojice = dict()
        for a in range(2,h//3 + 1):
                for b in range(a+1,(h-a)//2 + 1):
                        for c in range(b+1,h-a-b + 1):
                                index = (a+b+c,a*a+b*b+c*c)
                                if index not in trojice: trojice[index] = [(a,b,c)]
                                else: trojice[index].append((a,b,c))
                                
        for sucet in range(1,h):
                print(sucet)
                stvorice = dict()
                for a in range(2,sucet//4+1):
                        for b in range(a+1,(sucet-a)//3+1):
                                for c in range(b+1,(sucet-a-b)//2+1):
                                        d = sucet-a-b-c
                                        s = sucet
                                        t = a*a+b*b+c*c+d*d
                                        if (s-1,t-1) in trojice:
                                                for x,y,z in trojice[(s-1,t-1)]:
                                                        if len({1,a,b,c,d,x,y,z}) == 8:
                                                                index = (x,y,z)
                                                                if index not in stvorice: stvorice[index] = [(a,b,c,d)]
                                                                else: stvorice[index].append((a,b,c,d))
                                                
                maximum = 0
                for i,j in stvorice.items(): maximum = max(maximum,len(j))
                for n in range(5,maximum+1):
                        for x,stvorica in stvorice.items():
                                for C in combinations(stvorica,r=n-1):
                                        prvky = []
                                        for c in C: prvky += c
                                        if len(set(prvky).union({1,x[0],x[1],x[2]})) == 4*n:
                                                moznosti = []
                                                for i in range(len(C)):
                                                        moznosti.append([])
                                                        for P in permutations(C[i]): moznosti[i].append(P)
                                                for PP in product([y for y in range(24)],repeat=n-1):
                                                        obdlznik = [[1,x[0],x[1],x[2]]]
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
        for n in range(5,maximum+1):
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


