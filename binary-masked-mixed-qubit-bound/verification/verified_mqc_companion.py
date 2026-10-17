"""Correlated MQC face brackets using the corrected affine arithmetic."""
import math
from decimal import Decimal
from flint import arb
import verified_affine as ar


def regular_A(x):
    X=x.interval() if isinstance(x,ar.Aff) else x
    if X.lower()>arb(".75"):
        t=1-x
        total=1
        power=1
        N=24
        for j in range(1,N+1):
            power*=t
            total+=power/(2*j+1)
        radius=ar.abs_upper(t.interval() if isinstance(t,ar.Aff) else t)
        tail=radius**(N+1)/((2*N+3)*(1-radius))
        total+=ar.Aff(arb(0),{},ar.symmetric(tail)) if isinstance(total,ar.Aff) else ar.symmetric(tail)
        return 1/total
    z=(1-x).sqrt()
    return z/z.atanh()


def H(x):
    X=x.interval() if isinstance(x,ar.Aff) else x
    if X.lower()>arb(".75"):
        t=1-x
        total=arb(2).log()
        power=1
        N=24
        for j in range(1,N+1):
            power*=t
            total-=power/(2*j*(2*j-1))
        radius=ar.abs_upper(t.interval() if isinstance(t,ar.Aff) else t)
        tail=radius**(N+1)/((2*N+2)*(2*N+1)*(1-radius))
        total+=ar.Aff(arb(0),{},ar.symmetric(tail)) if isinstance(total,ar.Aff) else ar.symmetric(tail)
        return total
    r=(1-x).sqrt()
    p=x/(2*(1+r))
    return -(p*p.log()+(1-p)*(1-p).log())


def float_g(x):
    if x>=1-1e-7:
        z=1-x
        return 1+z/3+z*z/5+z**3/7
    r=math.sqrt(1-x)
    return .5*math.log((1+r)**2/x)/r


def float_A(x):
    return 0.0 if x<=0 else 1/float_g(x)


def float_gp(x):
    if x>=1-1e-5:
        z=1-x
        return -1/3-2*z/5-3*z*z/7-4*z**3/9
    return (float_g(x)-1/x)/(2*(1-x))


def face_center(a):
    lo,hi=1-a,1.0
    for _ in range(65):
        w=(lo+hi)/2
        r=1-(1-a)/w
        value=2*float_A(w*w*r)-(1+r)*float_A(w*w)
        if value>0:
            hi=w
        else:
            lo=w
    w=(lo+hi)/2
    z=w*(w-1+a)
    p=2*w-1+a
    Fa=-float_g(z)-w*p*float_gp(z)
    Fw=2*float_g(w*w)+4*w*w*float_gp(w*w)-2*float_g(z)-p*p*float_gp(z)
    return w,-Fa/Fw


def face_equation(w,a):
    r=1-(1-a)/w
    return 2*regular_A(r*w*w)-(1+r)*regular_A(w*w)


def interval(x):
    return ar.interval(x)


def entropy_interval(x):
    """Monotone range bound, including the regular endpoint x=1."""
    lo=max(x.lower(),arb(0))
    hi=min(x.upper(),arb(1))
    if lo>hi:
        raise ValueError("Entropy domain empty")
    def point_h(z):
        if z==0:
            return arb(0)
        if z==1:
            return arb(2).log()
        r=(1-z).sqrt()
        p=z/(2*(1+r))
        return -p*p.log()-(1-p)*(1-p).log()
    l,u=point_h(lo).lower(),point_h(hi).upper()
    return arb((l+u)/2,((u-l)/2).upper())


def companion(tube,Caff=None):
    a_raw=tube["a"].v
    coeff=dict(a_raw.a)
    if ar.abs_upper(a_raw.r)!=0:
        coeff["Arem"]=ar.abs_upper(a_raw.r)
    a=ar.Aff(a_raw.c,coeff,arb(0))
    ai=a.interval()
    su,lam,rho=tube["u"].v,tube["lam"].v,tube["rho"].v
    x=su*su
    y=rho*x
    strict=arb(2).log()+H(x)-H(lam)-H(y)
    C=(1-su+2*su*rho)/(1+su) if Caff is None else Caff
    K=(3+C*C-(1-C*C)*su)/4
    if ai.upper()<arb(".14883"):
        gap=strict-(arb(2).log()-H(a))
        latent=1/regular_A(a)-1/(K*regular_A(lam))
        return {"mode":"endpoint","gap":gap.interval(),"latent":latent.interval(),
                "a":ai,"strict":strict.interval()}
    if ai.lower()<=arb(".14884"):
        # A common near-endpoint bracket handles the switch without
        # pretending its optimizer is analytic through the endpoint.
        w0=ar.aff_const(".999")
        if not face_equation(w0,a).interval().upper()<0:
            raise ValueError("Face-switch lower bracket not certified")
        # Unique face minimization implies w in [.999,1], including w=1.
        w=arb(".9995",".0005")
        w=w.intersection(arb(".5",".5"))
        z=w*(w-1+ai)
        lower_face=entropy_interval(w*w)-entropy_interval(z)
        gap=strict.interval()-lower_face
        return {"mode":"switch_rectangle","gap":gap,"latent":None,
                "a":ai,"strict":strict.interval(),"face_lower":lower_face}
    center,slope=face_center(float(a.c))
    Wbase=ar.Aff(arb(str(center)),{name:arb(str(slope))*value for name,value in a.a.items()},arb(0))
    radius=max(arb("1e-14"),ar.abs_upper(a.interval()-a.c)*arb(".0001"))
    for _ in range(40):
        Wminus=Wbase-radius
        Wplus=Wbase+radius
        if Wminus.interval().lower()>0 and Wplus.interval().upper()<1 and (Wminus-1+a).interval().lower()>0:
            hm=face_equation(Wminus,a).interval()
            hp=face_equation(Wplus,a).interval()
            if hm.is_finite() and hp.is_finite() and hm.upper()<0 and hp.lower()>0:
                break
        radius*=arb("1.5")
    else:
        raise ValueError("Correlated face side signs not certified")
    wcoeff=dict(Wbase.a)
    wcoeff["FaceNoise"]=radius
    w=ar.Aff(Wbase.c,wcoeff,arb(0))
    z=w*(w-1+a)
    gap=strict-(H(w*w)-H(z))
    latent=w/regular_A(z)-1/(K*regular_A(lam))
    return {"mode":"interior","gap":gap.interval(),"latent":latent.interval(),
            "a":ai,"strict":strict.interval(),"face_side_signs":[hm,hp],
            "w_center":str(center),"w_radius":str(radius),"w_slope":str(slope)}


# The endpoint branch switch is enclosed, not assumed from decimal lookup.
assert (2*regular_A(arb(".14883"))-1-arb(".14883")).upper()<0
assert (2*regular_A(arb(".14884"))-1-arb(".14884")).lower()>0
