"""Restricted forced-l scalar FLOAT proxies after exact setup, plus exact finite initial-state controls."""
from fractions import Fraction as F
import random,time

def G(y):
    if y<=F(1,3):return y/2
    if y<=F(2,3):return 3*y/2-F(1,3)
    if y<=F(5,6):return 2*y-F(2,3)
    return F(1)

def anchor(K):
    m=(5*K+5)//6;h=F(1,K);b=F(m,K)
    return m,h,b,2*b-b*b,b*(1-b),K%6,F(6-K%6,3*K)

def discrepancies(p):
    K=len(p);pref=F(0);D=[]
    for j,v in enumerate(p,1):pref+=v;D.append(pref-F(j,K))
    return D

def update(p,s,g,alpha=F(1,4)):
    mass=alpha*sum(p[:s])
    result=[(1-alpha)*v+mass*g[i] if i<s else v+mass*g[i] for i,v in enumerate(p)]
    assert min(result)>=0 and sum(result)==1
    return result

def full_state_to_upper(K):
    m,h,b,c,a,r,q=anchor(K)
    p=[F(2*K-2*i+1,K*K) for i in range(1,K+1)]
    g=[G(F(i,K))-G(F(i-1,K)) for i in range(1,K+1)]
    for t in range(10000):
        D=discrepancies(p);s=D.index(max(D))+1
        if s==m:
            assert all(v<a for v in D[:m-2]) and p[m-1]>h
            return {"p":p,"g":g,"x":p[m-1],"setup_rounds":t,"a":a,"c":c,"q":q,"h":h,"m":m}
        p=update(p,s,g)
    raise AssertionError("upper state not reached")

def s12_noise(K,eps,H,seed):
    assert K%6<3 and 0<=eps<1
    state=full_state_to_upper(K);rng=random.Random(seed)
    x,a,c,q,h=map(float,(state["x"],state["a"],state["c"],state["q"],state["h"]))
    forced=selected_lower=plateau=0;run=0;runs=[];max_excess=0.0
    for _ in range(H):
        assert x>h
        d=a+max(h-x,0.0)
        plateau+=d==a;max_excess=max(max_excess,d-a)
        noise=rng.random()<eps
        if noise:
            forced+=1;selected_lower+=1;run+=1;x=x+(c-x)*q/4
        else:
            if run:runs.append(run);run=0
            x=0.75*x+c*q/4
    if run:runs.append(run)
    assert plateau==H and selected_lower==forced and max_excess==0
    return {"K":K,"eps":eps,"epsilon_fraction":str(F(str(eps))),"H":H,"seed":seed,
            "alpha":"1/4","arithmetic":"scalar float proxy after exact full-state setup",
            "observable_note":"d is inferred from boundary closure, not measured from a full mass vector",
            "setup_rounds":state["setup_rounds"],
            "forced_count":forced,"forced_share":forced/H,"unforced_share":(H-forced)/H,
            "selected_lower_share":selected_lower/H,"proxy_d_plateau_share":plateau/H,
            "proxy_max_d_excess":max_excess,"forced_run_count":len(runs),
            "max_forced_run":max(runs,default=0),"mean_forced_run":sum(runs)/len(runs) if runs else 0}

def initial_state(K,kind,seed=42):
    if kind=="canonical":p=[F(2*K-2*i+1,K*K) for i in range(1,K+1)]
    elif kind=="uniform":p=[F(1,K)]*K
    elif kind=="near_uniform":
        rng=random.Random(seed);weights=[F(rng.randint(90,110)) for _ in range(K)]
        p=[w/sum(weights) for w in weights]
    elif kind=="suffix_transfer_stop_control":
        p=[F(1,K)]*K;take=min(F(1,20),p[0]/2)
        p[-1]+=take;p[0]-=take
    else:raise ValueError(kind)
    assert min(p)>=0 and sum(p)==1
    return p

def s6_init(K,kind,seed=42,H=1000):
    m,h,b,c,a,r,q=anchor(K);p=initial_state(K,kind,seed)
    g=[G(F(i,K))-G(F(i-1,K)) for i in range(1,K+1)]
    a0=sum(p[:m])-F(m,K);Mstar=max(sum(p[:j])-F(j,K) for j in range(m+1,K+1))
    initial_d=max(discrepancies(p));initial_min=min(p);last=-1;stop=None
    for t in range(H):
        assert min(p)>=0 and sum(p)==1
        D=discrepancies(p);d=max(D);s=D.index(d)+1
        if s not in (m-1,m):last=t
        if d<=F(1,10**12):stop=t;break
        p=update(p,s,g)
    return {"K":K,"kind":kind,"H":H,"stopped_at":stop,"last_interior":last,
            "initial_d":str(initial_d),"initial_min_mass":str(initial_min),
            "a0":str(a0),"Mstar":str(Mstar),"mass_valid":True}

def main():
    started=time.time();noise=[];inits=[]
    for K in (50,200,8,14,26):
        for eps in (0.01,0.05,0.10,0.20):
            row=s12_noise(K,eps,200000,1000*K+int(100*eps));noise.append(row)
            print("Restricted scalar float proxy statistics:",row,flush=True)
    for K in (20,50,100):
        for kind in ("canonical","uniform","near_uniform","suffix_transfer_stop_control"):
            row=s6_init(K,kind);inits.append(row)
            if kind in ("uniform","suffix_transfer_stop_control"):assert row["stopped_at"]==0
            print("Valid finite initial-state case:",row,flush=True)
    print("PASS: scalar proxy d-plateau share=1 for restricted forced-l design; no general robustness claim.")
    return {"status":"PASS","restricted_noise":noise,"finite_initial_states":inits,
            "elapsed_seconds":round(time.time()-started,3)}

if __name__=="__main__":main()
