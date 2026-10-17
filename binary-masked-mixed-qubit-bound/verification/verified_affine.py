"""Outward Arb affine arithmetic for a fresh 2MQC certificate replay.

Semantics: c + sum(a[name]*epsilon[name]) + r, where every epsilon is in
[-1,1] and r is an Arb ball. Arb's second constructor argument is a radius,
NEVER an upper endpoint. Products include linear/remainder cross terms.
"""
from dataclasses import dataclass
from decimal import Decimal, localcontext

from flint import arb, ctx

ctx.prec = 256


def to_arb(x):
    return x if isinstance(x, arb) else arb(str(x))


def abs_upper(x):
    """An outward upper bound, itself retained as an exact dyadic Arb point."""
    return abs(x).upper()


def symmetric(radius):
    return arb(0, abs_upper(to_arb(radius)))


@dataclass
class Aff:
    c: arb
    a: dict
    r: arb

    def linear_radius(self):
        return sum((abs_upper(x) for x in self.a.values()), arb(0))

    def interval(self):
        return self.c + self.r + symmetric(self.linear_radius())

    def _co(self, x):
        return x if isinstance(x, Aff) else Aff(to_arb(x), {}, arb(0))

    def __add__(self, x):
        x = self._co(x)
        keys = self.a.keys() | x.a.keys()
        return Aff(self.c+x.c,
                   {k:self.a.get(k,arb(0))+x.a.get(k,arb(0)) for k in keys},
                   self.r+x.r)

    __radd__ = __add__

    def __neg__(self):
        return Aff(-self.c, {k:-v for k,v in self.a.items()}, -self.r)

    def __sub__(self, x):
        return self + -self._co(x)

    def __rsub__(self, x):
        return self._co(x) - self

    def __mul__(self, x):
        x = self._co(x)
        keys = self.a.keys() | x.a.keys()
        linear = {k:self.a.get(k,arb(0))*x.c+x.a.get(k,arb(0))*self.c for k in keys}
        ra = self.linear_radius()+abs_upper(self.r)
        rb = x.linear_radius()+abs_upper(x.r)
        remainder = x.c*self.r+self.c*x.r+symmetric(ra*rb)
        return Aff(self.c*x.c, linear, remainder)

    __rmul__ = __mul__

    def unary2(self, f, df, d2f):
        X = self.interval()
        deriv = df(self.c)
        radius = self.linear_radius()+abs_upper(self.r)
        remainder = deriv*self.r+symmetric(abs_upper(d2f(X))*radius*radius/2)
        return Aff(f(self.c), {k:deriv*v for k,v in self.a.items()}, remainder)

    def inv(self):
        return self.unary2(lambda x:1/x, lambda x:-1/(x*x), lambda x:2/(x*x*x))

    def __truediv__(self, x):
        return self * self._co(x).inv()

    def __rtruediv__(self, x):
        return self._co(x) * self.inv()

    def sqrt(self):
        return self.unary2(lambda x:x.sqrt(),lambda x:1/(2*x.sqrt()),lambda x:-1/(4*x*x.sqrt()))

    def log(self):
        return self.unary2(lambda x:x.log(),lambda x:1/x,lambda x:-1/(x*x))

    def atanh(self):
        return self.unary2(lambda x:x.atanh(),lambda x:1/(1-x*x),lambda x:2*x/(1-x*x)**2)


@dataclass
class D:
    v: Aff
    c: Aff
    u: Aff

    def _d(self,x):
        return x if isinstance(x,D) else D(aff_const(x),aff_const(0),aff_const(0))

    def __add__(self,x):
        x=self._d(x)
        return D(self.v+x.v,self.c+x.c,self.u+x.u)

    __radd__=__add__

    def __neg__(self):
        return D(-self.v,-self.c,-self.u)

    def __sub__(self,x):
        return self+-self._d(x)

    def __rsub__(self,x):
        return self._d(x)-self

    def __mul__(self,x):
        x=self._d(x)
        return D(self.v*x.v,self.c*x.v+self.v*x.c,self.u*x.v+self.v*x.u)

    __rmul__=__mul__

    def inv(self):
        vi=self.v.inv()
        return D(vi,-self.c*vi*vi,-self.u*vi*vi)

    def __truediv__(self,x):
        return self*self._d(x).inv()

    def __rtruediv__(self,x):
        return self._d(x)*self.inv()

    def unary2(self,f,df,d2f,d3f):
        value=self.v.unary2(f,df,d2f)
        deriv=self.v.unary2(df,d2f,d3f)
        return D(value,deriv*self.c,deriv*self.u)

    def sqrt(self):
        return self.unary2(lambda x:x.sqrt(),lambda x:1/(2*x.sqrt()),
                           lambda x:-1/(4*x*x.sqrt()),lambda x:3/(8*x*x*x.sqrt()))

    def log(self):
        return self.unary2(lambda x:x.log(),lambda x:1/x,lambda x:-1/(x*x),lambda x:2/(x*x*x))

    def atanh(self):
        return self.unary2(lambda x:x.atanh(),lambda x:1/(1-x*x),
                           lambda x:2*x/(1-x*x)**2,lambda x:2*(1+3*x*x)/(1-x*x)**3)


def aff_const(x):
    return Aff(to_arb(x),{},arb(0))


def var(c,r,name):
    return Aff(to_arb(c),{name:to_arb(r)},arb(0))


def A(x):
    z=(1-x).sqrt()
    return z/z.atanh()


def evaluate(C,u):
    d=1-C
    rho=((1+C)*u-d)/(2*u)
    x=u*u
    y=rho*x
    ax,ay=A(x),A(y)
    M=1+rho-2*ay/ax
    lam=(d*d-4*d*M)/(d*d+4*C*M)
    q=d+C*lam
    K=(3+C*C-(1-C*C)*u)/4
    # Legacy logarithmic residual log(R/H). On a zero its gradient is
    # -1/H times that of the manuscript's canonical residual H-R.
    F=(u*K/q).log()+A(lam).log()-ay.log()
    distance=K*lam-d*d*(1+u)/4
    J=distance.c*F.u-distance.u*F.c
    return {"F":F,"J":J,"a":distance,"lam":lam,"rho":rho,"M":M,"q":q,"u":u}


def eval_tube(c0,c1,center,slope,margin,side=None):
    with localcontext() as dc:
        dc.prec=120
        c0,c1,center,slope,margin=map(lambda x:Decimal(str(x)),(c0,c1,center,slope,margin))
        cm,cr=(c0+c1)/2,(c1-c0)/2
        Caff=var(cm,cr,"C")
        coeff={"C":to_arb(slope*cr)}
        if side is None:
            coeff["E"]=to_arb(margin)
        else:
            center += Decimal(side)*margin
        Uaff=Aff(to_arb(center),coeff,arb(0))
    C=D(Caff,aff_const(1),aff_const(0))
    u=D(Uaff,aff_const(0),aff_const(1))
    return evaluate(C,u)


def interval(x):
    if isinstance(x,D):
        x=x.v
    return x.interval() if isinstance(x,Aff) else x
