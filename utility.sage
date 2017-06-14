#!/usr/bin/env python

from itertools import combinations_with_replacement

def convolve(A,B):
    rdict = {}

    for p, v in A:
        for q, w in B:
            k = tuple(vector(q)+vector(p))
            if k not in rdict:
                rdict[k] = 0
            rdict[k] += v*w

    return [(vector(k),rdict[k]) for k in rdict if rdict[k] != 0]


def filter_to_zdomain(filter, L = None):
    """
    filter: a list of ((lattice_site),weight) tuples
    L: the generating matrix for the lattice
    
    returns: A symbolic expression in z. Note that this
    expression is not within the laurent polynomial ring 
    (it technically is, but needs to be cast as such)
    """
    if len(filter) == 0 or len(filter[0]) == 0:
        return 0
        
    s = len(filter[0][0])
    
    if L is not None and L.det() == 0:
        return 0
        
    if L is not None:
        Li = L.inverse()
        filter = [(Li*vector(v), w) for (v,w) in filter]
        has_frac = lambda x: any([_ != floor(_) for _ in x])
        if any([has_frac(v) for (v,w) in filter]):
            return 0
            
    # 
    Z = [var('z_%d' % i) for i in xrange(s)]
    
    return sum([
        c*prod([Z[i]**(-v) for i,v in enumerate(v)]) 
        for v,c in filter
    ])

def zdomain_to_filter(Zdom, s, L=None):
    Z = [var('z_%d' % i) for i in xrange(s)]
    r = []
    for t in SR(Zdom).expand().operands():
        rd = {Z[j]:1 for j in xrange(s)}
        c = t.subs(rd)
        t = t/c
        z = []
        for i in xrange(s):
            v = var('v')
            rd = {Z[j]:1 for j in xrange(s) if j != i }
            tp = t.subs({Z[i]:v})
            tp = tp.subs(rd)
            sgn = -1
            if tp.subs({v:2}) < 1:
                sgn = 1
                tp = 1/tp
            ops = tp.operands()
            p = 0
            if len(ops) == 0:
                if len(tp.variables()) > 0:
                    p = 1
                else:
                    p = 0
            elif len(ops) == 2:
                p = ops[1]
            else:
                print "XXX"
            z += [p*sgn]
        if L is not None:
            r += [(tuple(L*vector(z)), c)]
        else:
            r += [(tuple(z), c)]
    rd = {}
    for (t, c) in r:
        if t not in rd:
            rd[t] = 0
        rd[t] += c
    return [(k, rd[k]) for k in rd]

    
def filter_dft(filter):
    if len(filter) == 0 or len(filter[0]) == 0:
        return 0
        
    s = len(filter[0][0])
    w = vector([var('w_%d' % i) for i in xrange(s)])
    
    expression = sum([c*exp(-I*(w.dot_product(vector(v)))) for v,c in filter])
    return expression.simplify_rectform().simplify()


def extend_along(n, ww, dir):
    p = vector([0,0,0])
    f = []
    for w in [(-1)^i * binomial(n,i)/(2^n) for i in xrange(n + 1)]:
        f += [(tuple(p), w*ww)]
        p += dir
    return f 

def expl(expression):
    return [exp(_) for _ in expression]

def create_ring_filter(rings=2):
    orn = rings
    rings = ceil(rings)
    e0 = vector([ 0, 0, 0])
    e1 = vector([-1, 1, 1])
    e2 = vector([ 1,-1, 1])
    e3 = vector([ 1, 1,-1])
    e4 = vector([-1,-1,-1])
    
    s = set([tuple(sum([vector(v) for v in k])) for k in list(combinations_with_replacement([e0,e1,e2,e3,e4,-e1,-e2,-e3,-e4], rings))])
    
    rd = {}
    for p in s:
        k = vector(p).dot_product(vector(p))
        if k not in rd:
            rd[k] = []
        rd[k] += [p]
    verts =  sorted(rd.keys())[:orn]
    s = sum([rd[k] for k in verts], [])
    return [(vector(s),var('f_%d' %i)) for i,s in enumerate(s)]
    

def signed_permutation_matrices(n = 3):
    if n <= 1:
        yield matrix([[ 1]])
        yield matrix([[-1]])
    else:
        # Construct matrix out of sub matrices
        for i in xrange(n):
            for A in signed_permutation_matrices(n-1):
                P = A[0:i]
                Q = A[i:n]
                P = matrix(P.rows() + [vector([0]*(n-1))] + Q.rows())
                Q = matrix([1  if _ == i else 0 for _ in xrange(n) ]).transpose()
                yield matrix(block_matrix([Q,P], ncols=2).rows())
                Q[i,0] = -1
                yield matrix(block_matrix([Q,P], ncols=2).rows())
                
    
def impose_symmetry(filter):
    eclass = {}
    
    for _ in signed_permutation_matrices(3):
        for s,v in filter:
            s = tuple(_*s)
            
            if s not in eclass:
                eclass[s] = set()
            eclass[s] = eclass[s].union(set([v]))
    rdict = {}
    cid = 0
    for s in eclass:
        k = tuple(sorted([str(_) for _ in set(eclass[s])]))
        if k not in rdict:
            rdict[k] = var('c_%d' % cid)
            cid += 1
    rfilt = []    
    for s,v in filter:
        s = tuple(s)
        k = tuple(sorted([str(_) for _ in set(eclass[s])]))
        
        rfilt += [(s,rdict[k])]
    return (rfilt, [var('c_%d' % i) for i in xrange(cid)])


def find_dual_lp(m_0, sring=2):
    for ring in xrange(sring, 100):
        print "creating ring.."
        f = create_ring_filter(ring)
        
        print "giving it symmetry..."
        (f,cvars) = impose_symmetry(f)
        
        print "making filter.."
        f = CompactFilterElement(f)
        F = (m_0*f).flip_rw(pi_[0], L) + (m_0*f).flip_rw(pi_[1], L) + (m_0*f).flip_rw(pi_[2], L) + (m_0*f).flip_rw(pi_[3], L) + (m_0*f).flip_rw(pi_[4], L) + (m_0*f).flip_rw(pi_[5], L) + (m_0*f).flip_rw(pi_[6], L) + (m_0*f).flip_rw(pi_[7], L) 
        
        print "Solving..."
        print cvars
        soln = solve([v ==0 for _,v in F.filter if tuple(_) != (0,0,0) ] +  [v == 1 for _,v in F.filter if tuple(_) == (0,0,0)], cvars )
        
        if len(soln) > 0:
            print "Got solution"
            break
        else:
            print "No solution, increasing filter ring"
            
    soln = soln[0]
    rdict = {}
    for expression in soln:
        variable = expression.lhs()
        value = expression.rhs()
        rdict[variable] = value
        
    return CFE([(s,v.subs(rdict)) for s,v in f.filter])
