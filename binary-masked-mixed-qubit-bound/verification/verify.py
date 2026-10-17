"""Recompute the computational lemmas used in proof.pdf.

Dependencies: python-flint 0.9.0, mpmath 1.3.0, sympy 1.14.0.
Run: python -B verification/verify.py --work PATH_TO_EMPTY_DIRECTORY
The work directory receives logs and reports. No legacy research checkout,
stored sign strings, or previously asserted success flags are used.
The analytic reductions and endpoint arguments remain part of the paper.
"""
import argparse
import contextlib
from decimal import Decimal,localcontext
from fractions import Fraction as Q
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from flint import arb,ctx
import cover_2mqc_folds as low
import cover_2mqc_high as high
import verified_affine as ar
import verified_mqc_companion as face
import certify_2mqc_events as events

ctx.prec=384
ROOT=Path(__file__).resolve().parent


def read(name):
    return json.loads(gzip.decompress((ROOT/"data"/name).read_bytes()))


def root_cover(name,engine,mode):
    data=read(name)
    expected={
        "low_fold.json.gz":("1/60","53/100","1/65536","1"),
        "low_distance.json.gz":("1/60","53/100","1/65536","1"),
        "high_fold.json.gz":("53/100","24/25","1/20","1"),
        "high_distance.json.gz":("53/100","4/5","1/20","1"),
    }
    assert tuple(map(Q,data["domain"]))==tuple(map(Q,expected[name])),(name,"wrong theorem domain")
    terminal={tuple(x) for x in data["boxes"]}
    assert len(terminal)==len(data["boxes"])
    assert data["nodes"]==2*len(terminal)-1
    domain=tuple(map(Q,data["domain"]))
    width=domain[1]-domain[0]
    pending=[domain]
    nodes=checked=0
    while pending:
        box=pending.pop()
        nodes+=1
        assert nodes<=data["nodes"],(name,"uncovered subdivision")
        key=tuple(map(str,box))
        if key in terminal:
            terminal.remove(key)
            assert engine.classify(box,mode) is not None,(name,key,"unresolved function enclosure")
            checked+=1
            if checked%10000==0:
                print(name,"recomputed",checked,flush=True)
            continue
        a,b,c,d=box
        if (b-a)/width>d-c:
            m=(a+b)/2
            pending.extend([(m,b,c,d),(a,m,c,d)])
        else:
            m=(c+d)/2
            pending.extend([(a,b,m,d),(a,b,c,m)])
    assert not terminal and nodes==data["nodes"]
    return {"domain":data["domain"],"nodes":nodes,"terminal_boxes_recomputed":checked}


def endpoint_roots():
    data=read("endpoint.json.gz")
    C=Q(53,100)
    parts=[]
    for index,ends in enumerate(data["brackets"]):
        a,b=map(Q,ends)
        cc=low.point(C);d=1-cc;t=low.ball(a,b)
        u=(d+2*cc*t)/(1+cc)
        x=u*u;rho=cc*t/u;y=rho*x
        M=1+rho-2*low.g_components(x)[0]/low.g_components(y)[0]
        lam=d*(d-4*M)/(d*d+4*cc*M)
        assert u.lower()>0 and u.upper()<1 and rho.lower()>0
        assert M.lower()>0 and M.upper()<d/4
        assert lam.lower()>x.upper() and lam.upper()<1
        F,A=low.evaluate((C,C,a,b))
        Fa,_=low.evaluate((C,C,a,a));Fb,_=low.evaluate((C,C,b,b))
        if index==0:
            assert F.t.upper()<0 and Fa.v.lower()>0 and Fb.v.upper()<0
            J=A.c*F.t-A.t*F.c
            assert J.lower()>0
        else:
            assert F.t.lower()>0 and Fa.v.upper()<0 and Fb.v.lower()>0
        parts.append((a,b))
    for ends in data["excluded"]:
        a,b=map(Q,ends)
        try:
            F,_=low.evaluate((C,C,a,b))
        except low.Empty:
            pass
        else:
            assert low.nz(F.v),("endpoint root exclusion",ends)
        parts.append((a,b))
    parts.sort()
    assert parts[0][0]==Q(1,65536) and parts[-1][1]==1
    assert all(x[1]==y[0] for x,y in zip(parts,parts[1:]))
    lower=1-2*low.g_components(low.point(Q(59,61)**2))[0]/low.g_components(low.point(Q(1,60)))[0]
    assert lower.lower()>arb(".2589403965")
    # Necessary-root cutoffs used in the two rectangle covers.
    assert low.g_components(low.point(Q(9,100)))[0].upper()<2
    assert low.g_components(low.point(Q(1,49)**2))[0].upper()<5
    assert high.entropy_remainder(low.Jet(low.point(Q(1,9)))).v.upper()<low.point(Q(2,25))
    assert low.point(Q(470,353)+Q(4,25))<arb(20).log()/2
    return {"simple_roots_at_C_053":2,"excluded_intervals":len(data["excluded"]),
            "initial_no_root_interval":["0","1/60"],"lower_root_J_positive":True,
            "low_t_and_high_theta_cutoffs_checked":True}


def value_cover():
    cases=read("values.json.gz")
    assert Decimal(cases[0]["C"][0])==Decimal(".41904")
    assert Decimal(cases[-1]["C"][1])==Decimal(".95")
    assert all(Decimal(a["C"][1])==Decimal(b["C"][0]) for a,b in zip(cases,cases[1:]))
    for flag,start,end in (("monotone",".41905",".55"),("negative",".53",".95"),
                           ("distance",".8",".95"),("positive",".41904",".41905")):
        rows=[r for r in cases if r[flag]]
        assert Decimal(rows[0]["C"][0])==Decimal(start)
        assert Decimal(rows[-1]["C"][1])==Decimal(end)
        assert all(Decimal(a["C"][1])==Decimal(b["C"][0]) for a,b in zip(rows,rows[1:]))
    for number,row in enumerate(cases,1):
        args=(*row["C"],row["u_center"],row["slope"],row["u_radius"])
        tube=ar.eval_tube(*args)
        left=ar.eval_tube(*args,side=-1)
        right=ar.eval_tube(*args,side=1)
        assert ar.interval(left["F"]).lower()>0
        assert ar.interval(right["F"]).upper()<0
        assert ar.interval(tube["F"].u).upper()<0
        for key in ("M","rho","lam","a","q","u"):
            assert ar.interval(tube[key]).is_finite() and ar.interval(tube[key]).lower()>0
        for key in ("rho","lam","a","u"):
            assert ar.interval(tube[key]).upper()<1
        with localcontext() as dc:
            dc.prec=120
            a,b=map(Decimal,row["C"])
            C=ar.var((a+b)/2,(b-a)/2,"C")
        comparison=face.companion(tube,C)
        gap=comparison["gap"]
        assert gap.is_finite()
        if row["positive"]:
            assert (gap/arb(2).log()).lower()>arb(".00102288236")
        if row["negative"]:
            assert gap.upper()<0
        if row["monotone"]:
            assert comparison["latent"] is not None and comparison["latent"].lower()>0
        if row["distance"]:
            assert ar.interval(tube["J"]).lower()>0
        if number%1000==0:
            print("value intervals recomputed",number,flush=True)
    return {"intervals_recomputed":len(cases),
            "monotone_intervals":sum(r["monotone"] for r in cases),
            "negative_intervals":sum(r["negative"] for r in cases),
            "upper_distance_intervals":sum(r["distance"] for r in cases)}


def local_events(work):
    data=read("events.json.gz")
    reports=[]
    with (work/"events.log").open("w",encoding="utf-8") as log:
        with contextlib.redirect_stdout(log):
            for item in data:
                reports.append(events.certify(item["event"],item["center_decimal"]))
    for report in reports:
        box=[arb(arb(x),arb(report["radius"])) for x in report["center_decimal"]]
        _,_,d=events.system(box,report["event"])
        if report["event"]=="birth":
            assert box[0].lower()>arb(".25237734804709390") and box[0].upper()<arb(".25237734804709396")
        elif report["event"]=="spinodal":
            assert box[0].lower()>arb(".41904") and box[0].upper()<arb(".41905")
        else:
            assert box[0].lower()>arb(".52734600028746") and box[0].upper()<arb(".52734600028748")
            assert d["F"].t.lower()>0
            delta=d["a"].v/2
            assert delta.lower()>arb(".2002251206412457") and delta.upper()<arb(".2002251206412894")
            assert (1-d["lam"].v).lower()>arb(".3782489593555764")
            jump=d["slope_difference"]
            assert jump.lower()>arb(".0747075697651890") and jump.upper()<arb(".0747075697653768")
    return reports


def main():
    if sys.flags.optimize:
        raise RuntimeError("Run without -O: assertions are part of verification.")
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work",type=Path,required=True)
    args=parser.parse_args()
    work=args.work.resolve()
    if work.exists() and any(work.iterdir()):
        raise ValueError("Use an empty work directory; existing results are preserved.")
    work.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter()
    manifest=json.loads((ROOT/"SHA256SUMS.json").read_text())
    for name,digest in manifest.items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    summary={"integrity_files_checked":len(manifest)}
    # These scripts generate their reports only in the requested work directory.
    env=dict(os.environ)
    env["PYTHONPATH"]=str(ROOT)+os.pathsep+env.get("PYTHONPATH","")
    for script in ("verify_2mqc_tail_rationals.py","verify_2mqc_analytic_constants.py"):
        shutil.copy2(ROOT/script,work/script)
        with (work/(script+".log")).open("w",encoding="utf-8") as log:
            subprocess.run([sys.executable,"-B",str(work/script)],cwd=work,env=env,
                           stdout=log,stderr=subprocess.STDOUT,check=True)
    summary["endpoint"]=endpoint_roots()
    summary["events"]=local_events(work)
    summary["covers"]={}
    for name,engine,mode in (("low_fold.json.gz",low,"fold"),("high_fold.json.gz",high,"fold"),
                              ("low_distance.json.gz",low,"distance"),("high_distance.json.gz",high,"distance")):
        print("Recomputing",name,flush=True)
        summary["covers"][name]=root_cover(name,engine,mode)
    summary["values"]=value_cover()
    summary["all_computational_lemmas_verified"]=True
    summary["scope"]="Combine these computational lemmas with the analytic reductions and endpoint arguments in proof.pdf."
    summary["elapsed_seconds"]=time.perf_counter()-start
    (work/"verification_result.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print(json.dumps({k:v for k,v in summary.items() if k!="events"},indent=2),flush=True)


if __name__=="__main__":
    main()
