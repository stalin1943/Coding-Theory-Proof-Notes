"""Fresh conditional interval cover for the canonical stationary equation.

Only necessary physical-root constraints are used to contract value intervals;
derivative intervals always enclose the original analytic expressions at those
points. Each accepted box either excludes a physical root, excludes a multiple
root, or proves the required orientation. A nonempty pending list is NOT a
certificate of the stated global conclusion.
"""
import argparse
from dataclasses import dataclass
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import time

from flint import arb, fmpq, ctx

ctx.prec=192


def point(x):
    if isinstance(x,arb):
        return x
    x=Q(x)
    return arb(fmpq(x.numerator,x.denominator))


def ball(lo,hi):
    lo,hi=point(lo),point(hi)
    return arb((lo+hi)/2,abs((hi-lo)/2).upper())


def clip(x,lo,hi):
    lo,hi=point(lo),point(hi)
    if x.upper()<lo or x.lower()>hi:
        raise Empty()
    return x.intersection(ball(lo,hi))


class Empty(Exception):
    pass


@dataclass
class Jet:
    v: arb
    c: arb=arb(0)
    t: arb=arb(0)
    cc: arb=arb(0)
    ct: arb=arb(0)
    tt: arb=arb(0)

    def co(self,x):
        return x if isinstance(x,Jet) else Jet(point(x))

    def __add__(self,x):
        x=self.co(x)
        return Jet(self.v+x.v,self.c+x.c,self.t+x.t,
                   self.cc+x.cc,self.ct+x.ct,self.tt+x.tt)

    __radd__=__add__

    def __neg__(self):
        return Jet(-self.v,-self.c,-self.t,-self.cc,-self.ct,-self.tt)

    def __sub__(self,x):
        return self+-self.co(x)

    def __rsub__(self,x):
        return self.co(x)-self

    def __mul__(self,x):
        x=self.co(x)
        return Jet(self.v*x.v,self.c*x.v+self.v*x.c,self.t*x.v+self.v*x.t,
                   self.cc*x.v+2*self.c*x.c+self.v*x.cc,
                   self.ct*x.v+self.c*x.t+self.t*x.c+self.v*x.ct,
                   self.tt*x.v+2*self.t*x.t+self.v*x.tt)

    __rmul__=__mul__

    def chain(self,f,fp,fpp):
        return Jet(f,fp*self.c,fp*self.t,
                   fpp*self.c*self.c+fp*self.cc,
                   fpp*self.c*self.t+fp*self.ct,
                   fpp*self.t*self.t+fp*self.tt)

    def inv(self):
        return self.chain(1/self.v,-1/(self.v*self.v),2/(self.v*self.v*self.v))

    def __truediv__(self,x):
        return self*self.co(x).inv()

    def __rtruediv__(self,x):
        return self.co(x)*self.inv()

    def restrict(self,lo,hi):
        return Jet(clip(self.v,lo,hi),self.c,self.t,self.cc,self.ct,self.tt)


def g_components(z):
    """g and two derivatives; regular series includes z=1 exactly."""
    if z.lower()>point(Q(3,4)):
        r=1-z
        R=abs(r).upper()
        N=28
        v=arb(0); p=arb(0); pp=arb(0)
        for j in range(N-1,-1,-1):
            pp=pp*r+2*p
            p=p*r+v
            v=v*r+point(Q(1,2*j+1))
        tail=R**N/((2*N+1)*(1-R))
        tailp=(N*R**(N-1)/(1-R)+R**N/(1-R)**2)/(2*N+1)
        tailpp=(N*(N-1)*R**(N-2)/(1-R)+2*N*R**(N-1)/(1-R)**2+
                2*R**N/(1-R)**3)/(2*N+1)
        return v+arb(0,tail.upper()),-p+arb(0,tailp.upper()),pp+arb(0,tailpp.upper())
    s=(1-z).sqrt()
    v=s.atanh()/s
    p=(v-1/z)/(2*(1-z))
    pp=(p+1/(z*z))/(2*(1-z))+(v-1/z)/(2*(1-z)**2)
    return v,p,pp


def g(z):
    return z.chain(*g_components(z.v))


def evaluate(box,physical=True):
    c0,c1,t0,t1=box
    C=Jet(ball(c0,c1),arb(1))
    t=Jet(ball(t0,t1),arb(0),arb(1))
    d=1-C
    u=((d+2*C*t)/(1+C)).restrict((1-c1)/(1+c1),1)
    rho=(C*t/u).restrict(0,Q(3,20) if physical else c1)
    x=(u*u).restrict(((1-c1)/(1+c1))**2,1)
    if physical:
        if rho.v.lower()>=x.v.upper()/4:
            raise Empty()
        rho=rho.restrict(0,x.v.upper()/4)
    y=rho*x if physical else C*t*u
    # On the low strip every physical stationary root has g(x)<2 and
    # g(y)<16/3; g(z)>=.5*log(1/z) implies y>2^-16.
    if physical:
        y=y.restrict(Q(1,2**16),Q(3,20))
    gx,gy=g(x),g(y)
    M=1+rho-2*gx/gy
    if physical:
        M=M.restrict(0,d.v.upper()/4)
    elif not (M.v.is_finite() and M.v.lower()>0 and M.v.upper()<d.v.lower()/4):
        raise ValueError("The entire rectangle is not inside the analytic lambda domain")
    lam=d*(d-4*M)/(d*d+4*C*M)
    if physical:
        lam=lam.restrict(x.v.lower(),1)
    elif not (lam.v.is_finite() and lam.v.lower()>0 and lam.v.upper()<1):
        raise ValueError("Natural lambda enclosure not interior")
    qq=d+C*lam
    K=(3+C*C-(1-C*C)*u)/4
    F=qq*g(lam)-u*K*gy
    dist=K*lam-d*d*(1+u)/4
    return F,dist


def nz(x):
    return x.is_finite() and (x.lower()>0 or x.upper()<0)


def classify(box,mode):
    try:
        F,a=evaluate(box)
    except Empty:
        return "physical_exclusion"
    except (ValueError,ZeroDivisionError):
        return None
    if nz(F.v):
        return "residual_exclusion"
    # A centered mean-value enclosure is used only when the original
    # analytic formula is defined throughout the WHOLE convex rectangle.
    # Conditional root contractions alone do not justify a segment argument.
    try:
        whole,whole_a=evaluate(box,physical=False)
        c0,c1,t0,t1=box
        cm,tm=(c0+c1)/2,(t0+t1)/2
        center,center_a=evaluate((cm,cm,tm,tm),physical=False)
        dc=ball(-(c1-c0)/2,(c1-c0)/2)
        dt=ball(-(t1-t0)/2,(t1-t0)/2)
        centered=center.v+whole.c*dc+whole.t*dt
        if nz(centered):
            return "centered_residual_exclusion"
        F.c=F.c.intersection(center.c+whole.cc*dc+whole.ct*dt)
        F.t=F.t.intersection(center.t+whole.ct*dc+whole.tt*dt)
        a.c=a.c.intersection(center_a.c+whole_a.cc*dc+whole_a.ct*dt)
        a.t=a.t.intersection(center_a.t+whole_a.ct*dc+whole_a.tt*dt)
    except (ValueError,Empty,ZeroDivisionError):
        pass
    if mode=="fold":
        if nz(F.t):
            return "simple_fiber"
        orientation=F.c*F.tt
        if orientation.is_finite() and orientation.upper()<0:
            return "positive_birth_orientation"
    else:
        J=a.c*F.t-a.t*F.c
        if nz(J):
            return "noncritical_distance"
        Jc=a.cc*F.t+a.c*F.ct-a.ct*F.c-a.t*F.cc
        Jt=a.ct*F.t+a.c*F.tt-a.tt*F.c-a.t*F.ct
        N=Jc*F.t-Jt*F.c
        if F.t.is_finite() and F.t.lower()>0 and N.is_finite() and N.upper()<0:
            return "upper_strict_distance_maximum"
    return None


def encode(box):
    return [str(x) for x in box]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--mode",choices=("fold","distance"),default="fold")
    parser.add_argument("--max-nodes",type=int,default=20000)
    parser.add_argument("--resume",action="store_true")
    parser.add_argument("--upgrade-v1",action="store_true")
    args=parser.parse_args()
    out=Path(__file__).with_name("fresh_low_"+args.mode+"_cover.json")
    sourcehash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.resume:
        previous=json.loads(out.read_text())
        history=previous.get("engine_history",[])
        if previous["engine_sha256"]!=sourcehash:
            old=Path(__file__).with_name("cover_2mqc_folds_v1.py")
            assert args.upgrade_v1 and hashlib.sha256(old.read_bytes()).hexdigest()==previous["engine_sha256"]
            history.append(previous["engine_sha256"])
        else:
            assert not args.upgrade_v1
        pending=[tuple(map(Q,x)) for x in previous["pending"]]
        leaves=previous["leaves"]
        nodes=previous["nodes"]
    else:
        assert not out.exists(), "Use --resume or a separately named workspace; preserve existing evidence."
        pending=[(Q(1,60),Q(53,100),Q(1,2**16),Q(1))]
        leaves=[]
        nodes=0
        history=[]
    startnodes=nodes
    start=time.perf_counter()
    while pending and nodes-startnodes<args.max_nodes:
        box=pending.pop()
        kind=classify(box,args.mode)
        nodes+=1
        if kind:
            leaves.append({"box":encode(box),"conclusion":kind})
        else:
            c0,c1,t0,t1=box
            if (c1-c0)/Q(77,150)>t1-t0:
                mid=(c0+c1)/2
                pending.extend([(mid,c1,t0,t1),(c0,mid,t0,t1)])
            else:
                mid=(t0+t1)/2
                pending.extend([(c0,c1,mid,t1),(c0,c1,t0,mid)])
        if nodes%5000==0:
            print(json.dumps({"nodes":nodes,"leaves":len(leaves),"pending":len(pending),
                              "seconds":round(time.perf_counter()-start,2)}),flush=True)
    result={"mode":args.mode,"engine_sha256":sourcehash,"nodes":nodes,"leaves":leaves,
            "engine_history":history,
            "pending":[encode(x) for x in pending],"complete_cover":not pending,
            "domain":["1/60","53/100","1/65536","1"],
            "analytic_inputs":["projective stationary equivalence","quarter barrier",
                               "lambda > u^2 at every stationary root","g(y)<16/3 on the low strip"],
            "global_theorem_verified":False,"elapsed_last_run":time.perf_counter()-start}
    out.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k not in ("leaves","pending")}),flush=True)
    print("remaining_boxes",len(pending),flush=True)


if __name__=="__main__":
    main()
