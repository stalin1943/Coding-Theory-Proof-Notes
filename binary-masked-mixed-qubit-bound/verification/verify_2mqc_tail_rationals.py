"""Exact rational Bernstein certificates for the analytic endpoint proof."""
import json
from math import comb
from pathlib import Path
import sympy as s

d,u,gamma=s.symbols("d u gamma")
R=s.Rational
C=1-d
U=d*(2+d+d*d)/4
uL=d/(2-d)
K=lambda z:(3+C*C-(1-C*C)*z)/4
md=d*(1-d)/(4*(2-d))
z=U*U
t=(d+d*d)/2
B0=2*K(U)/((2-d)**2*(1+U*U/4-md))
Rc=B0*(1-t)
r0=2*(t-t*t/2+t**3/3)+z+z*z/(2*(1-z))
beta=R(5,4)
qbar=1+beta*C
mbeta=d*(1-beta*d)/(4*qbar)
sx=1-U*U/2-U**4/2
Pb=2*U*K(uL)/(d*qbar*(1-mbeta)*sx)
Pe=U*U/(4-U*U-U**4)
Dg=d*d*(1+d/(2-d))*(3-d)/4
Gs=Dg/2-4*Pe
m2=d*(1-d*d)/(4*(1+C*d))
B2=2*K(U)/((2-d)*(1+C*d)*(1+U*U/4-m2))
N=8
Tu=64+8*u**3-8*u**4+u**6
PN=4*(2*N+1)*u*(1+u)*(8+u*u)-(2*N+1)*u*u*Tu*sum(
    (1-u*u)**j/R(2*j+1) for j in range(N))-Tu*(1-u*u)**N

rho=d*d*gamma
ug=d/(2-d-2*rho)
x=ug*ug
Q1=lambda w:(d+C*w)*(1+rho)-d*d*(1-w)/4
Q1prime=C*(1+rho)+d*d/4
phi_margin=R(9,100)*d-(Q1(x)-2*ug*K(ug)+Q1prime*x/4)/(2*ug*K(ug))


def coefficients(poly,variables,box):
    dummy=s.symbols("z0:"+str(len(variables)))
    transformed=s.Poly(s.expand(poly.subs({
        v:a+(b-a)*z for v,z,(a,b) in zip(variables,dummy,box)},simultaneous=True)),*dummy,domain=s.QQ)
    degrees=transformed.degree_list()
    powers=transformed.as_dict()
    if len(variables)==1:
        n=degrees[0]
        return [sum(value*R(comb(k,i[0]),comb(n,i[0])) for i,value in powers.items() if i[0]<=k)
                for k in range(n+1)]
    n,m=degrees
    return [sum(value*R(comb(k,i),comb(n,i))*R(comb(l,j),comb(m,j))
                for (i,j),value in powers.items() if i<=k and j<=l)
            for k in range(n+1) for l in range(m+1)]


def certify(poly,variables,box):
    pending=[(box,0)]
    leaves=[]
    while pending:
        current,depth=pending.pop()
        vals=coefficients(poly,variables,current)
        if min(vals)>0:
            leaves.append({"box":[[str(a),str(b)] for a,b in current],
                           "minimum_Bernstein_coefficient":str(min(vals)),
                           "coefficients":[str(v) for v in vals]})
            continue
        if depth>=16:
            raise AssertionError(("Unresolved exact polynomial sign",poly,current,min(vals)))
        axis=depth%len(variables)
        left=list(current);right=list(current)
        a,b=current[axis];mid=(a+b)/2
        left[axis]=(a,mid);right[axis]=(mid,b)
        pending.extend([(tuple(right),depth+1),(tuple(left),depth+1)])
    return leaves


def rational_certificate(name,expression,variables,box):
    num,den=map(s.expand,s.fraction(s.cancel(expression)))
    middle={v:(a+b)/2 for v,(a,b) in zip(variables,box)}
    if den.subs(middle)<0:
        num,den=-num,-den
    removed=0
    if d in variables:
        while num.subs(d,0)==0:
            num=s.cancel(num/d).expand()
            removed+=1
    num_cover=certify(num,variables,box)
    den_cover=certify(den,variables,box)
    print(name,"verified",len(num_cover),len(den_cover),flush=True)
    return {"name":name,"variables":[str(v) for v in variables],
            "removed_positive_endpoint_power":removed,"numerator":str(num),"denominator":str(den),
            "numerator_cover":num_cover,"denominator_cover":den_cover}


tests=[
    ("lambda_one_entropy_majorant",PN,(u,),((R(1,5),R(1)),)),
    ("cubic_upper_support",(2-d)*U-U**3/2-d,(d,),((R(0),R(1,5)),)),
    ("central_slope",2*Rc-1,(d,),((R(0),R(1,20)),)),
    ("central_residual",(2*Rc-1)*R(69,50)+d/2-Rc*r0,(d,),((R(0),R(1,20)),)),
    ("upper_bracket_slope",1-2*Pb,(d,),((R(0),R(1,20)),)),
    ("upper_bracket_intercept",-((2*Pb-1)*R(299,100)+Pb*(R(139,50)-d)
                               -(R(29,25)-beta*d/(1-beta*d))),(d,),((R(0),R(1,20)),)),
    ("upper_gap_slope",Gs,(d,),((R(0),R(1,20)),)),
    ("upper_gap_intercept",Gs*R(299,100)+Dg/2*(R(29,25)-beta*d/(1-beta*d))
                         -R(189,25)*Pe,(d,),((R(0),R(1,20)),)),
    ("lower_central_factor",B2-(1-3*d/4),(d,),((R(0),R(1,5)),)),
    ("lower_root_displacement_geometry",phi_margin,(d,gamma),((R(0),R(1,5)),(R(0),R(49,625)))),
]
results=[rational_certificate(*test) for test in tests]
report={"passed":True,"arithmetic":"exact SymPy rationals; Bernstein convex-hull bounds",
        "scope":"Only the explicitly listed rational inequalities in the analytic endpoint arguments",
        "certificates":results,"global_theorem_verified":False}
Path(__file__).with_name("2mqc_tail_rational_certificates.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
print("All",len(results),"rational inequalities certified.",flush=True)
