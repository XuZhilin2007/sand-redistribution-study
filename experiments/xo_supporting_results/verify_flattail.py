"""Exact flat-tail checks: unconditional monotonic floor; conditional conservation."""
from fractions import Fraction as F
import random,time

def canonical_init(K):return [F(2*K-2*i+1,K*K) for i in range(1,K+1)]

def make_flat_tail_law_safe(K,m,rng):
    weights=[F(rng.randint(1,50)) for _ in range(m)]
    total=sum(weights)
    return [w/total for w in weights]+[F(0)]*(K-m)

def run(K,g,p0,alpha,H,tol=F(1,10**12),m=None):
    assert 1<=m<K and 0<alpha<1
    assert len(g)==len(p0)==K and min(g)>=0 and min(p0)>=0
    assert sum(g)==sum(p0)==1 and all(v==0 for v in g[m:])
    p=list(p0);fm0=sum(p[:m]);a0=fm0-F(m,K)
    Mstar=max(sum(p[:j])-F(j,K) for j in range(m+1,K+1))
    dominance=Mstar<=a0
    old_floor=max(a0,Mstar);changed=old_fail=0;previous=fm0;last=-1;smax=0
    initial_d=None;stopped=None
    for t in range(H):
        pref=F(0);D=[]
        for j,v in enumerate(p,1):pref+=v;D.append(pref-F(j,K))
        d=max(D);s=D.index(d)+1;fm=sum(p[:m])
        assert min(p)>=0 and pref==1 and fm>=previous and d>=a0
        if t==0:initial_d=d
        changed+=fm!=fm0;old_fail+=d<old_floor
        if dominance:assert fm==fm0 and s<=m
        if s not in (m-1,m):last=t
        smax=max(smax,s)
        if d<=tol:
            assert a0<=tol
            stopped=t;break
        mass=alpha*sum(p[:s])
        predicted=fm if s<=m else fm+alpha*(sum(p[:s])-fm)
        p=[(1-alpha)*v+mass*g[i] if i<s else v+mass*g[i] for i,v in enumerate(p)]
        assert sum(p[:m])==predicted
        previous=fm
    return {"K":K,"m":m,"a0":str(a0),"Mstar":str(Mstar),"initial_d":str(initial_d),
            "weak_suffix_dominance":dominance,"Fm_changed_states":changed,
            "old_Mstar_floor_failed_states":old_fail,"stopped_at":stopped,
            "last_interior":last,"s_max":smax,"H":H,
            "unconditional_monotone_floor_ok":True}

def main():
    started=time.time();rng=random.Random(20261006);A=[];B=[]
    for _ in range(40):
        K=rng.randint(8,40);m=rng.randint(K//2+1,K-1)
        g=make_flat_tail_law_safe(K,m,rng)
        A.append(run(K,g,canonical_init(K),F(1,4),150,m=m))
    for _ in range(20):
        K=rng.randint(8,30);m=rng.randint(K//2+1,K-1)
        g=make_flat_tail_law_safe(K,m,rng);p=canonical_init(K)
        j=rng.randint(m+1,K);take=min(F(1,10),p[0]/2)
        p[j-1]+=take;p[0]-=take
        B.append(run(K,g,p,F(1,4),150,m=m))
    positive=run(8,[F(1,6)]*6+[F(0)]*2,
                 [F(77,600)]*6+[F(21,100),F(1,50)],F(1,4),100,m=6)
    assert not positive["weak_suffix_dominance"] and F(positive["a0"])>F(1,10**12)
    assert positive["Fm_changed_states"]>0 and positive["stopped_at"] is None
    K=12
    def G(y):
        if y<=F(1,3):return y/2
        if y<=F(2,3):return 3*y/2-F(1,3)
        if y<=F(5,6):return 2*y-F(2,3)
        return F(1)
    xo=[G(F(i,K))-G(F(i-1,K)) for i in range(1,K+1)]
    sanity=run(K,xo,canonical_init(K),F(1,4),2000,m=10)
    assert sanity["weak_suffix_dominance"] and sanity["Fm_changed_states"]==0
    assert F(10**12-1,10**24)<F(1,10**12)
    stats={"canonical_runs":len(A),"perturbed_runs":len(B),
           "perturbed_true_violation":sum(not r["weak_suffix_dominance"] for r in B),
           "perturbed_Fm_changed":sum(r["Fm_changed_states"]>0 for r in B),
           "perturbed_stopped":sum(r["stopped_at"] is not None for r in B),
           "canonical_confined_before_half":sum(r["last_interior"]<75 for r in A),
           "canonical_interior_at_final_state":sum(r["last_interior"]==149 for r in A)}
    print("PASS: unconditional Fm monotonicity and a0 floor in all valid states.")
    print("Finite sample statistics:",stats)
    print("Positive-floor, suffix-violating example:",positive)
    print("No infinite-time boundary-entry claim from the H=150 sample.")
    return {"status":"PASS","statistics":stats,"canonical":A,"perturbed":B,
            "positive_floor_counterexample_to_necessity":positive,"xo_sanity":sanity,
            "elapsed_seconds":round(time.time()-started,3)}

if __name__=="__main__":main()
