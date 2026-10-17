"""Whole-interval entropy and exact constant checks for the endpoint proof."""
import json
from fractions import Fraction as Q
from pathlib import Path
from flint import arb,ctx
import cover_2mqc_folds as cv

ctx.prec=256
p=cv.point
def g(z):
    r=(1-z).sqrt()
    return r.atanh()/r
def J0(U):
    x=U*U
    G=((2*arb(2).sqrt()+x)**2/(arb(2).sqrt()*(2-x))).log()/3
    return G-U*(6+x)*(2+U)/(3*(8+x))
pending=[(Q(0),Q(1))]
leaves=[]
while pending:
    a,b=pending.pop()
    value=J0(cv.ball(a,b))
    if value.is_finite() and value.lower()>0:
        leaves.append({"interval":[str(a),str(b)],"J0_enclosure":value.str(45,more=True)})
    else:
        assert b-a>Q(1,2**20)
        mid=(a+b)/2
        pending.extend([(mid,b),(a,mid)])

checks=[]
def above(name,left,right):
    assert left>right,(name,left,right)
    checks.append({"name":name,"left":str(left),"right":str(right)})
above("rho0_upper",p(Q(40,23)),g(p(Q(3,20))).upper())
above("rho0_lower",g(p(Q(1,16))).lower(),p(Q(32,17)))
above("boundary_small_u",p(Q(5,2)),g(p(Q(1,25))).upper())
above("upper_tail_x_entropy",g(p(Q(1,1000))).lower(),p(3))
above("larger_tail_x_entropy",g(p(Q(1,75))).lower(),p(Q(5,2)))
above("d_g_d_squared_small",p(Q(1,4)),(p(Q(1,20))*g(p(Q(1,400)))).upper())
above("d_g_d_squared_large",p(Q(1,2)),(p(Q(1,5))*g(p(Q(1,25)))).upper())
above("d_g_x_large",p(Q(5,8)),(p(Q(1,5))*g(p(Q(1,100)))).upper())
above("lower_central_entropy_increment",arb(Q(25,14).numerator).log()-arb(Q(25,14).denominator).log(),p(Q(57,100)))
above("lower_displacement_entropy_increment",(p(Q(250,563))*p(Q(5,4)).log()).lower(),p(Q(9,100)))
above("lower_separator_entropy_increment",(p(Q(125,278))*p(Q(50,13)).log()).lower(),p(Q(59,100)))
for zmax,beta in ((Q(1,4000000),Q(1,4)),(Q(1,25000),Q(7,15))):
    z=p(zmax)
    upper=z.root(4)*(4/z).log()/(2*(1-z).sqrt())
    above("scaled_entropy_"+str(zmax),p(beta),upper.upper())
    above("elasticity_remainder_"+str(zmax),p(1/beta-2),(z.root(4)**3).upper())
above("g_y_above_8",arb(16000000).log()/2,p(8))

Ds=Q(3897,2000);Ks=Q(19,20);us=Q(3,100);qs=Q(79,80)
Dl=Q(89,50);Kl=Q(89,100);ul=Q(14,125);ql=Q(19,20);rl=Q(49,15625)
above("small_T_over_V",Q(1,40),2*us*((Q(1,4)+Q(1,4000)/Ds)/qs+1/(6*Ds)))
above("small_U_over_V",Q(1,50),us/(Ks*Ds))
above("large_T_over_V",Q(1,10),2*ul*((Q(1,4)+rl/Dl)/ql+1/(5*Dl)))
above("large_U_over_V",Q(2,25),ul/(Kl*Dl))
brace=1-Q(3,3800)-2*Q(40,39)*Q(1,6)-2*Q(40,39)*Q(1,4000)/qs-Q(4,16)*Q(1,50)
above("Gd_brace",brace,Q(16,25))
above("Gd_coefficient",Q(512,821)-Q(103,190)*Q(1,20),Q(59,100))
above("lambda_d_lower",(Q(59,100)-Q(20,39))/(Q(20,39)+Q(2,79)),Q(7,50))
ulambda=Q(40,39)*us*us/(4*Ds*qs)
above("a_derivative_coefficient",Ks-Q(7,8)*Q(1,400)*ulambda,Q(949,1000))
above("a_d_positive",Q(949,1000)*Q(7,50)-Q(103,3200)-Q(103,4000)-Q(21,16000),Q(7,100))
# Bound the correction to u_d at fixed lambda. 0<=partial_d m<=1/4.
bcoef=(1+us)/Ds
correction=us*us/(Ds*Q(39,40))*(2*Q(1,4000)*bcoef/(us*qs)+Q(7,12)*bcoef/us+Q(1,4)/qs)
above("section_u_d_upper",Q(3,5),bcoef+correction)
above("section_u_d_lower",Q(1,2)-correction,Q(0))
result={"passed":True,"J0_interval_cover":leaves,"constant_checks":checks,
        "scope":"Uniform analytic-lemma constants and J0 positivity, not a stand-alone global phase proof."}
Path(__file__).with_name("2mqc_analytic_constant_certificates.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
print("J0 leaves",len(leaves),"constant inequalities",len(checks),"all passed",flush=True)

