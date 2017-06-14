
class CompactFilterElement(RingElement):
    def __init__(self, F, parent = SR):
        if len(F) == 0:
            self.filter = [((0,0,0), 0)]
        else:
            rdict = {tuple(k): v for k,v in F}
            self.filter = [(vector(k),parent(rdict[k])) for k in rdict if rdict[k] != 0]
        RingElement.__init__(self, parent)

    def plot(self):
        return points([(k[0],k[1], k[2]) for k,v in self.filter])

    def to_polynomial(self):
        z_ = [var('z_%d' % i) for i in xrange(len(self.filter[0][0]))]

        return sum([ v*prod([_**__ for (_,__) in zip(z_, p)])  for p,v in self.filter])

    def dft(self):
        w = len(self.filter[0])
        return filter_dft(self.filter)

    def _add_(self, other):
        rdict = {tuple(k):v for k,v in self.filter}
        for k,v in other.filter:
            k = tuple(k)
            if k not in rdict:
                rdict[k] = 0
            rdict[k] += v
        filter = [(vector(k),self.parent()(rdict[k])) for k in rdict]
        return CompactFilterElement(filter, self.parent())

    def _sub_(self, other):
        rdict = {tuple(k):v for k,v in self.filter}
        for k,v in other.filter:
            k = tuple(k)
            if k not in rdict:
                rdict[k] = 0
            rdict[k] -= v
        filter = [(vector(k),self.parent()(rdict[k])) for k in rdict]
        return CompactFilterElement(filter, self.parent())

    def _mul_(self,other):
        if type(other) == sage.symbolic.expression.Expression or type(other)==sage.rings.integer.Integer:
            A = [(p,self.parent()(v*other)) for p,v in self.filter]
            return CompactFilterElement(A, self.parent())
        A = self.filter
        B = other.filter
        filter = convolve(A,B)
        return CompactFilterElement(filter, self.parent())

    def _neg_(self):
        rdict = {tuple(k): v for k,v in self.filter}
        return CompactFilterElement([(vector(k),self.parent()(-rdict[k])) for k in rdict if rdict[k] != 0], self.parent())

    def flip(self):
        rdict = {tuple(k): v for k,v in self.filter}
        return CompactFilterElement([(-vector(k),self.parent()(rdict[k])) for k in rdict if rdict[k] != 0], self.parent())

    def flip_rw(self, pie, L):
        rdict = {} # {tuple(vector(k)): v for k,v in self.filter}
        Li = L.inverse()

        for k,v in self.filter:
            n = vector(k)
            rdict[tuple(-n)] = exp(I*pie.dot_product(Li*n))*v
        return CompactFilterElement([(vector(k),self.parent()(rdict[k])) for k in rdict if rdict[k] != 0], self.parent())

CFE = CompactFilterElement
