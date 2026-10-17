"""Fresh interval contraction proofs for birth, spinodal and coexistence.

mpmath only locates rational-decimal centers and an approximate inverse.
Existence and uniqueness are proved by an outward Arb Jacobian enclosure,
strict Krawczyk inclusion and a contraction norm less than one.
"""
import json
from pathlib import Path
import mpmath as mp
from flint import arb,ctx
import cover_2mqc_folds as cv

ctx.prec=384
mp.mp.dps=110


def g_components(z):
    # All local event boxes are separated from z=1. The direct analytic
    # formula avoids the fixed-order series tail used for global cells.
    r=(1-z).sqrt()
    value=r.atanh()/r
    first=(value-1/z)/(2*(1-z))
    second=(first+1/(z*z))/(2*(1-z))+(value-1/z)/(2*(1-z)**2)
    return value,first,second


def g(z):
    return z.chain(*g_components(z.v))


def geom(C,u):
    d=1-C
    rho=((1+C)*u-d)/(2*u)
    x=u*u;y=rho*x
    M=1+rho-2*g(x)/g(y)
    lam=d*(d-4*M)/(d*d+4*C*M)
    qq=d+C*lam
    K=(3+C*C-(1-C*C)*u)/4
    a=K*lam-d*d*(1+u)/4
    F=qq*g(lam)-u*K*g(y)
    return {"rho":rho,"x":x,"y":y,"M":M,"lam":lam,"q":qq,"K":K,"a":a,"F":F}


def H(z):
    r=(1-z).sqrt()
    p=z/(2*(1+r))
    return -p*p.log()-(1-p)*(1-p).log()


def system(values,kind):
    C=cv.Jet(values[0],arb(1))
    u=cv.Jet(values[1],arb(0),arb(1))
    d=geom(C,u)
    F,a=d["F"],d["a"]
    if kind=="birth":
        return [F.v,F.t],[[F.c,F.t],[F.ct,F.tt]],d
    if kind=="spinodal":
        J=a.c*F.t-a.t*F.c
        Jc=a.cc*F.t+a.c*F.ct-a.ct*F.c-a.t*F.cc
        Ju=a.ct*F.t+a.c*F.tt-a.tt*F.c-a.t*F.ct
        return [F.v,J],[[F.c,F.t],[Jc,Ju]],d
    w=values[2]
    z=w*(w-1+a.v)
    gw,gwp,_=g_components(w*w)
    gz,gzp,_=g_components(z)
    face=2*w*gw-(2*w-1+a.v)*gz
    face_a=-gz-w*(2*w-1+a.v)*gzp
    face_w=2*gw+4*w*w*gwp-2*gz-(2*w-1+a.v)**2*gzp
    gap=arb(2).log()+H(d["x"].v)-H(d["lam"].v)-H(d["y"].v)-H(w*w)+H(z)
    gx=g_components(d["x"].v)[0]
    gl=g_components(d["lam"].v)[0]
    gy=g_components(d["y"].v)[0]
    gc=(gx*d["x"].c-gl*d["lam"].c-gy*d["y"].c+gz*w*a.c)/2
    gu=(gx*d["x"].t-gl*d["lam"].t-gy*d["y"].t+gz*w*a.t)/2
    gapw=-w*gw+(2*w-1+a.v)*gz/2
    d["face_w"]=w;d["face_image"]=z
    d["slope_difference"]=(w*gz-gl/d["K"].v)/arb(2).log()
    return [F.v,face,gap],[[F.c,F.t,arb(0)],[face_a*a.c,face_a*a.t,face_w],[gc,gu,gapw]],d


def mg(z):
    r=mp.sqrt(1-z)
    return mp.atanh(r)/r


def mh(z):
    p=z/(2*(1+mp.sqrt(1-z)))
    return -p*mp.log(p)-(1-p)*mp.log(1-p)


def mpgeom(C,u):
    d=1-C
    rho=((1+C)*u-d)/(2*u);x=u*u;y=rho*x
    M=1+rho-2*mg(x)/mg(y)
    lam=d*(d-4*M)/(d*d+4*C*M)
    qq=d+C*lam;K=(3+C*C-(1-C*C)*u)/4
    return qq*mg(lam)-u*K*mg(y),K*lam-d*d*(1+u)/4,x,y,lam,K


def mpfunc(kind,*values):
    C,u=values[:2]
    F,a,x,y,lam,K=mpgeom(C,u)
    if kind=="birth":
        return F,mp.diff(lambda z:mpgeom(C,z)[0],u)
    if kind=="spinodal":
        Fc=mp.diff(lambda z:mpgeom(z,u)[0],C)
        Fu=mp.diff(lambda z:mpgeom(C,z)[0],u)
        ac=mp.diff(lambda z:mpgeom(z,u)[1],C)
        au=mp.diff(lambda z:mpgeom(C,z)[1],u)
        return F,ac*Fu-au*Fc
    w=values[2];z=w*(w-1+a)
    return F,2*w*mg(w*w)-(2*w-1+a)*mg(z),mp.log(2)+mh(x)-mh(lam)-mh(y)-mh(w*w)+mh(z)


def label(x):
    return x.str(65,more=True)


def certify(kind,guess):
    n=len(guess)
    f=lambda *x:mpfunc(kind,*x)
    root=mp.findroot(f,tuple(mp.mpf(x) for x in guess),tol=mp.mpf("1e-95"),maxsteps=60)
    strings=[mp.nstr(root[i],100) for i in range(n)]
    center=[arb(x) for x in strings]
    radius=arb("1e-24")
    box=[arb(x,radius) for x in center]
    f0,_,_=system(center,kind)
    _,J,data=system(box,kind)
    for key in ("x","y","rho","lam","a"):
        v=data[key].v
        assert v.is_finite() and v.lower()>0 and v.upper()<1,(kind,key,v)
    assert data["lam"].v.lower()>data["x"].v.upper()
    assert data["M"].v.lower()>0 and data["M"].v.upper()<(1-box[0]).lower()/4
    if kind=="coexistence":
        assert box[2].upper()<1 and box[2].lower()>(1-data["a"].v).upper()
        assert data["slope_difference"].lower()>0
    # Analytic Jacobian is also checked against an independently differentiated
    # high-precision point formula. This is diagnostic, not the interval proof.
    _,Jcenter,_=system(center,kind)
    mpcenter=[mp.mpf(x) for x in strings]
    numerical=mp.matrix(n,n)
    for i in range(n):
        for j in range(n):
            def component(z):
                args=list(mpcenter);args[j]=z
                return f(*args)[i]
            numerical[i,j]=mp.diff(component,mpcenter[j])
            difference=Jcenter[i][j]-arb(mp.nstr(numerical[i,j],100))
            assert abs(difference).upper()<arb("1e-60"),(kind,i,j,difference)
    inverse=numerical**-1
    Y=[[arb(mp.nstr(inverse[i,j],95)) for j in range(n)] for i in range(n)]
    R=[[arb(int(i==j))-sum((Y[i][k]*J[k][j] for k in range(n)),arb(0))
        for j in range(n)] for i in range(n)]
    contraction=max(sum((abs(R[i][j]).upper() for j in range(n)),arb(0)) for i in range(n))
    assert contraction.upper()<1
    Kraw=[center[i]-sum((Y[i][j]*f0[j] for j in range(n)),arb(0))+
          sum((R[i][j]*(box[j]-center[j]) for j in range(n)),arb(0)) for i in range(n)]
    assert all(Kraw[i].lower()>box[i].lower() and Kraw[i].upper()<box[i].upper() for i in range(n))
    report={"event":kind,"center_decimal":strings,"radius":"1e-24",
            "boxes":[label(x) for x in box],"Krawczyk_image":[label(x) for x in Kraw],
            "contraction_infinity_norm":label(contraction),"strict_inclusion":True,
            "distance_a":label(data["a"].v),"relative_distance":label(data["a"].v/2),
            "lambda":label(data["lam"].v),"global_identification_verified":False}
    if kind=="birth":
        report["F_C"]=label(J[0][0]);report["F_uu"]=label(J[1][1])
        assert J[0][0].upper()<0 and J[1][1].lower()>0
        birthJ=data["a"].c*data["F"].t-data["a"].t*data["F"].c
        assert birthJ.lower()>0
        report["distance_J"]=label(birthJ)
    if kind=="spinodal":
        assert data["F"].t.lower()>0
        report["F_u"]=label(data["F"].t)
    if kind=="coexistence":
        report["slope_difference_d_delta"]=label(data["slope_difference"])
    print(json.dumps(report,indent=2),flush=True)
    return report


def main():
    records=[]
    for kind,guess in (("birth",(".252377",".65")),("spinodal",(".419046",".43")),
                       ("coexistence",(".527346",".3167",".65"))):
        records.append(certify(kind,guess))
        Path(__file__).with_name("fresh_2mqc_event_certificates.json").write_text(
            json.dumps({"events":records,"all_local_contractions_verified":len(records)==3},indent=2),encoding="utf-8")
    result={"events":records,"all_local_contractions_verified":True,
            "scope":"Local existence and uniqueness only; global branch identification and exclusivity are separate proofs."}
    Path(__file__).with_name("fresh_2mqc_event_certificates.json").write_text(json.dumps(result,indent=2),encoding="utf-8")


if __name__=="__main__":
    main()
