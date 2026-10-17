"""Conditional interval cover in the implicit cubic (C,theta) chart.

For C>=53/100 and 0<=theta<=1, (1+C)u-theta*u^3/2=1-C
has a unique root in (0,1/3). Endpoint monotonicity encloses that root
over each rectangle. Floating-point Newton values only propose brackets;
both bracket signs are then proved with Arb before any enclosure is used.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import time
from flint import arb,ctx
import cover_2mqc_folds as low

ctx.prec=224
Jet=low.Jet
point=low.point
ball=low.ball
g=low.g
Empty=low.Empty


def log(z):
    return z.chain(z.v.log(),1/z.v,-1/(z.v*z.v))


def sqrt(z):
    r=z.v.sqrt()
    return z.chain(r,1/(2*r),-1/(4*z.v*r))


@lru_cache(maxsize=100000)
def scalar_u_bracket(c,theta):
    cf,tf=float(c),float(theta)
    proposal=(1-cf)/(1+cf)
    for _ in range(8):
        proposal-=((1+cf)*proposal-tf*proposal**3/2-(1-cf))/(1+cf-1.5*tf*proposal**2)
    C,T=point(c),point(theta)
    radius=point(Q(1,10**16))
    center=arb(str(proposal))
    f=lambda u:(1+C)*u-T*u*u*u/2-(1-C)
    for _ in range(20):
        U=arb(center,radius.upper())
        lo,hi=U.lower(),U.upper()
        if lo>0 and hi<point(Q(1,3)) and f(lo).upper()<0 and f(hi).lower()>0:
            return lo,hi
        radius*=2
    raise RuntimeError(("Implicit cubic endpoint not enclosed",str(c),str(theta)))


def implicit_u(box):
    c0,c1,t0,t1=box
    # du/dC<0 and du/dtheta>0 on the entire physical cubic root branch.
    lower=scalar_u_bracket(c1,t0)[0]
    upper=scalar_u_bracket(c0,t1)[1]
    u=ball(lower,upper)
    C,T=ball(c0,c1),ball(t0,t1)
    A=1+C-point(Q(3,2))*T*u*u
    assert A.lower()>0
    uc=-(u+1)/A
    ut=u*u*u/(2*A)
    ucc=-(2*uc-3*T*u*uc*uc)/A
    uct=-(ut-point(Q(3,2))*u*u*uc-3*T*u*uc*ut)/A
    utt=(3*u*u*ut+3*T*u*ut*ut)/A
    return Jet(u,uc,ut,ucc,uct,utt)


def entropy_remainder(z):
    # Exact cancellation-free g(z)-log(4/z)/2 for small positive z.
    r=sqrt(1-z)
    return z*log(4/z)/(2*r*(1+r))+log(1-z/(2*(1+r)))/r


def evaluate(box,physical=True):
    c0,c1,t0,t1=box
    C=Jet(ball(c0,c1),arb(1))
    theta=Jet(ball(t0,t1),arb(0),arb(1))
    u=implicit_u(box)
    d=1-C
    x=u*u
    rho=theta*x/4
    y=rho*x
    gy=g(y)
    difference=-log(theta)/2+entropy_remainder(y)-2*entropy_remainder(x)
    M=rho+difference/gy
    if physical:
        M=M.restrict(0,d.v.upper()/4)
    elif not (M.v.is_finite() and M.v.lower()>0 and M.v.upper()<d.v.lower()/4):
        raise ValueError("Whole rectangle not inside the analytic lambda domain")
    lam=d*(d-4*M)/(d*d+4*C*M)
    if physical:
        lam=lam.restrict(x.v.lower(),1)
    elif not (lam.v.is_finite() and lam.v.lower()>0 and lam.v.upper()<1):
        raise ValueError("Whole lambda interval is not interior")
    qq=d+C*lam
    K=(3+C*C-(1-C*C)*u)/4
    F=qq*g(lam)-u*K*gy
    # Exact stable rearrangement prevents lower-branch distance cancellation.
    den=2-d-2*rho
    base=d*d*rho*(2*d*rho-3*d+2*rho*rho-6*rho+4)/(den*den*den)
    distance=base+K*(lam-u*u)
    return F,distance


def classify(box,mode):
    try:
        F,a=evaluate(box)
    except Empty:
        return "physical_exclusion"
    except (ValueError,ZeroDivisionError):
        return None
    if low.nz(F.v):
        return "residual_exclusion"
    try:
        whole,whole_a=evaluate(box,physical=False)
        c0,c1,t0,t1=box
        cm,tm=(c0+c1)/2,(t0+t1)/2
        center,center_a=evaluate((cm,cm,tm,tm),physical=False)
        dc=ball(-(c1-c0)/2,(c1-c0)/2)
        dt=ball(-(t1-t0)/2,(t1-t0)/2)
        centered=center.v+whole.c*dc+whole.t*dt
        if low.nz(centered):
            return "centered_residual_exclusion"
        F.c=F.c.intersection(center.c+whole.cc*dc+whole.ct*dt)
        F.t=F.t.intersection(center.t+whole.ct*dc+whole.tt*dt)
        a.c=a.c.intersection(center_a.c+whole_a.cc*dc+whole_a.ct*dt)
        a.t=a.t.intersection(center_a.t+whole_a.ct*dc+whole_a.tt*dt)
    except (ValueError,Empty,ZeroDivisionError):
        pass
    if mode=="fold":
        if low.nz(F.t):
            return "simple_fiber"
    else:
        J=a.c*F.t-a.t*F.c
        if low.nz(J):
            return "noncritical_distance"
    return None


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--mode",choices=("fold","distance"),default="fold")
    parser.add_argument("--max-nodes",type=int,default=100000)
    parser.add_argument("--resume",action="store_true")
    args=parser.parse_args()
    c0=Q(53,100)
    c1=Q(24,25) if args.mode=="fold" else Q(4,5)
    t0,t1=Q(1,20),Q(1)
    # Strong uniform separation: E(z)=g(z)-log(4/z)/2 increases from 0.
    # u<1/3 gives E(x)<2/25, and g(x)<5. At a physical level
    # D=g(y)-2g(x)=(M-rho)g(y)<470/353 because M<47/400.
    # Hence -.5 log(theta) <470/353+4/25 <.5 log(20).
    assert low.g_components(point(Q(1,49)**2))[0].upper()<5
    assert entropy_remainder(Jet(point(Q(1,9)))).v.upper()<point(Q(2,25))
    assert point(Q(470,353)+Q(4,25))<arb(20).log()/2
    out=Path(__file__).with_name("fresh_high_conditioned_"+args.mode+"_cover.json")
    hashes={Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            Path(low.__file__).name:hashlib.sha256(Path(low.__file__).read_bytes()).hexdigest()}
    if args.resume:
        previous=json.loads(out.read_text())
        assert previous["engine_hashes"]==hashes
        pending=[tuple(map(Q,b)) for b in previous["pending"]]
        leaves=previous["leaves"];nodes=previous["nodes"]
    else:
        assert not out.exists(), "Preserve existing evidence: use --resume"
        pending=[(c0,c1,t0,t1)];leaves=[];nodes=0
    firstnodes=nodes
    start=time.perf_counter()
    while pending and nodes-firstnodes<args.max_nodes:
        box=pending.pop()
        kind=classify(box,args.mode)
        nodes+=1
        if kind:
            leaves.append({"box":low.encode(box),"conclusion":kind})
        else:
            a,b,c,d=box
            if (b-a)/(c1-c0)>d-c:
                mid=(a+b)/2
                pending.extend([(mid,b,c,d),(a,mid,c,d)])
            else:
                mid=(c+d)/2
                pending.extend([(a,b,mid,d),(a,b,c,mid)])
        if nodes%5000==0:
            print(json.dumps({"nodes":nodes,"leaves":len(leaves),"pending":len(pending),
                              "seconds":round(time.perf_counter()-start,2)}),flush=True)
    result={"mode":args.mode,"engine_hashes":hashes,"nodes":nodes,"leaves":leaves,
            "pending":[low.encode(b) for b in pending],"complete_cover":not pending,
            "domain":[str(c0),str(c1),str(t0),str(t1)],
            "analytic_inputs":["stationary equivalence","quarter barrier",
                               "lambda > u^2 at a stationary root","implicit cubic is increasing on [0,1/3]",
                               "entropy-remainder separation theta>1/20, proved by the checked constants"],
            "global_theorem_verified":False,"elapsed_last_run":time.perf_counter()-start}
    out.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k not in ("leaves","pending")}),flush=True)
    print("remaining_boxes",len(pending),flush=True)


if __name__=="__main__":
    main()
