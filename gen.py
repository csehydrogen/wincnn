# -*- coding: utf-8 -*-
from wincnn import cookToomFilter
from sympy import pprint, Rational, IndexedBase, Matrix, simplify

def genMeta(r, m, G, genDst, genSrc):
    # P1 P2 P3 P6 P5 P4 P7 P8 P9 order
    ijs = []
    for i in xrange(0, 1):
        for j in xrange(0, 1):
            ijs.append(i * m + j)
    for i in xrange(0, 1):
        for j in xrange(1, m - 1):
            ijs.append(i * m + j)
    for i in xrange(0, 1):
        for j in xrange(m - 1, m):
            ijs.append(i * m + j)
    for i in xrange(1, m - 1):
        for j in xrange(m - 1, m):
            ijs.append(i * m + j)
    for i in xrange(1, m - 1):
        for j in xrange(1, m - 1):
            ijs.append(i * m + j)
    for i in xrange(1, m - 1):
        for j in xrange(0, 1):
            ijs.append(i * m + j)
    for i in xrange(m - 1, m):
        for j in xrange(0, 1):
            ijs.append(i * m + j)
    for i in xrange(m - 1, m):
        for j in xrange(1, m - 1):
            ijs.append(i * m + j)
    for i in xrange(m - 1, m):
        for j in xrange(m - 1, m):
            ijs.append(i * m + j)

    # 매트릭스 원소별, 계수 별, src Index 별 저장
    gi = IndexedBase('g');
    g = Matrix(r, r, lambda i, j : gi[i * r + j])
    U = G * g * G.T
    dbs = []
    for ij in ijs:
        i, j = ij // m, ij % m
        u = U[i,j]
        db = {}
        dbs.append(db)
        for k in xrange(r * r):
            c = u.coeff(gi[k])
            if c == 0:
                continue
            ac = abs(c)
            if not ac in db:
                db[ac] = []
            db[ac].append((1 if c == ac else -1, k))

    # cc가 캐시 역할, gs는 캐시에 언제 대입하는지 저장
    line = 0
    cc = []
    gs = []
    for i in xrange(m * m):
        cc.append(set())
        gs.append(set())

    for z in xrange(1024):
        # x + y 꼴의 빈도를 센다
        stat = {}
        for ij in xrange(m * m):
            db = dbs[ij]
            for ac in db:
                ks = db[ac]
                l = len(ks)
                for i in xrange(l):
                    if ks[i][0] == 0:
                        continue
                    for j in xrange(i + 1, l):
                        if ks[j][0] == 0:
                            continue
                        t = (ks[i], ks[j])
                        if not t in stat:
                            stat[t] = [0, []]
                        stat[t][0] += 1
                        stat[t][1].append(ij)
        # 없으면 종료
        if len(stat) == 0:
            print("stat is empty : %d" % z)
            break
        p = max([[k, v] for k, v in stat.items()], key = lambda x : x[1][0])
        # 있지만 단 한번 나오는 x + y라면 굳이 치환할 필요가 없다.
        if p[1][0] <= 1:
            print("useless substitution : %d" % z)
            break
        x, y = p[0]
        start, end = p[1][1][0], p[1][1][-1]
        # life 동안 빈 캐시가 있는지 확인, 없으면 추가
        for l in xrange(line + 1):
            flag = True
            if l < line:
                for i in xrange(start, end + 1):
                    if l in cc[i]:
                        flag = False
                        break
            else:
                line += 1
            if flag:
                gs[start].add((l, x, y))
                for i in xrange(start, end + 1):
                    cc[i].add(l)
                break
        # 치환
        for ij in p[1][1]:
            db = dbs[ij]
            for ac in db:
                ks = db[ac]
                if x in ks and y in ks:
                    ks.remove(x)
                    ks.remove(y)
                    ks.append((0, l))

    def sign(x):
        return "+" if x > 0 else "-"
    def genStmt(ij):
        # 캐시 대입 코드 생성
        for it in gs[ij]:
            s = "t[%d] = %s %s %s %s;" % (it[0], sign(it[1][0]), genSrc(it[1][1]), sign(it[2][0]), genSrc(it[2][1]))
            print(s)
        db = dbs[ij]
        # dst 대입 코드 생성
        s = ""
        s += genDst(ijs[ij] // m, ijs[ij] % m)
        for k in db:
            v = db[k]
            s += " + %d./%d. * (" % (k.p, k.q)
            for it in v:
                if it[0] == 0: # 캐시에 있는 경우
                    s += " + t[%d]" % (it[1])
                else:
                    s += " %s %s" % (sign(it[0]), genSrc(it[1]))
            s += ")"
        s += ";"
        print(s)

    print("Dtype t[%d];" % (line))
    for ij in xrange(m * m):
        genStmt(ij)

def genKernel(a, n, r):
    AT, G, BT, f = cookToomFilter(a, n, r)
    m = n + r - 1
    print("winoWeight")
    genMeta(r, m, G,
            lambda i, j : "dst[gIdx + %d * gap] =" % (i * m + j),
            lambda k : "src[kIdx + %d]" % k)
    print("winoSrc")
    genMeta(m, m, BT,
            lambda i, j : "dst[bIdx + %d * gap] =" % (i * m + j),
            lambda k : "src[sIdx + %d * dataW + %d]" % (k // m, k % m))
    print("winoDst")
    genMeta(m, n, AT,
            lambda i, j : "dst[rIdx + %d * outW + %d] =" % (i, j),
            lambda k : "src[mIdx + %d * gap]" % k)

# F(2x2, 3x3)
genKernel((0, 1, -1), 2, 3)
# F(4x4, 3x3)
genKernel((0, 1, -1, 2, -2), 4, 3)
# F(6x6, 3x3)
genKernel((0, 1, -1, 2, -2, Rational(1, 2), -Rational(1, 2)), 6, 3)
