"""Exact radial aggregation, noninjectivity, and swept-ring deviation checks."""
from fractions import Fraction as F
import random

def radial_agg(p):return [sum(row) for row in p]

def one_step_2d(p,s,alpha,w):
    mass=alpha*sum(sum(row) for row in p[:s])
    return [[(1-alpha)*v+mass*w[r][t] if r<s else v+mass*w[r][t]
             for t,v in enumerate(row)] for r,row in enumerate(p)]

def one_step_1d(p,s,alpha,g):
    mass=alpha*sum(p[:s])
    return [(1-alpha)*v+mass*g[r] if r<s else v+mass*g[r] for r,v in enumerate(p)]

def radial_selection(p):
    P=radial_agg(p);pref=F(0);D=[]
    for j,v in enumerate(P,1):pref+=v;D.append(pref-F(j,len(P)))
    return D.index(max(D))+1,max(D)

def main():
    rng=random.Random(20261006);R,T=12,5;alpha=F(1,4)
    p=[[F(rng.randint(1,50)) for _ in range(T)] for _ in range(R)]
    total=sum(map(sum,p));p=[[v/total for v in row] for row in p]
    q=[list(reversed(row)) for row in p]
    weights=[[F(rng.randint(0,30)) for _ in range(T)] for _ in range(R)]
    weights[-1]=[F(0)]*T
    total=sum(map(sum,weights));w=[[v/total for v in row] for row in weights]
    g=radial_agg(w);one=radial_agg(p)
    assert p!=q and radial_agg(p)==radial_agg(q)
    steps=0
    for t in range(5):
        s,d=radial_selection(p);sq,dq=radial_selection(q)
        assert s==sq and d==dq
        if d<=F(1,10**12):break
        pnew=one_step_2d(p,s,alpha,w);qnew=one_step_2d(q,s,alpha,w)
        onenew=one_step_1d(one,s,alpha,g)
        assert radial_agg(pnew)==radial_agg(qnew)==onenew
        assert pnew!=qnew
        assert min(v for row in pnew for v in row)>=0 and sum(map(sum,pnew))==1
        for r in range(R):
            if g[r]==0:continue
            nu=[v/g[r] for v in w[r]]
            for theta in range(T):
                before=p[r][theta]-nu[theta]*one[r]
                after=pnew[r][theta]-nu[theta]*onenew[r]
                assert after==(1-alpha)*before if r<s else after==before
        p,q,one=pnew,qnew,onenew;steps+=1
    assert steps>0
    print("PASS: exact radial aggregation over active steps:",steps)
    print("PASS: distinct angular states with same radial trajectories; aggregation is many-to-one.")
    print("PASS: deviations contract on swept rings and stay unchanged on unswept rings.")
    return {"status":"PASS","R":R,"T":T,"verified_active_steps":steps,
            "zero_kernel_ring":R,"many_to_one_demonstrated":True}

if __name__=="__main__":main()
