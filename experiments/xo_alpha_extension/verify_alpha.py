"""Exact finite validation of Theorem A; finite float runs support, not prove, T5."""
from fractions import Fraction
import time

F = Fraction
TOL = F(1, 10**12)

def G(y):
    if y <= F(1,3): return y/2
    if y <= F(2,3): return 3*y/2-F(1,3)
    if y <= F(5,6): return 2*y-F(2,3)
    return F(1)

def anchor(K):
    m=(5*K+5)//6; h=F(1,K); b=F(m,K)
    return m,h,b,2*b-b*b,b*(1-b),K%6,F(6-K%6,3*K)

def classify(s_list,m,K):
    first=next((t for t,s in enumerate(s_list) if s in (m-1,m)),None)
    tail=0
    for s in reversed(s_list):
        if s!=m: break
        tail+=1
    post=set(s_list[first:]) if first is not None else set()
    label="LOCK" if tail>=K else "CYCLE" if first is not None and post=={m-1,m} else "OTHER"
    return {"t_star":first,"R_m":tail,"class":label,"post_set":sorted(post)}

def run_alpha(K,alpha,H,max_full=None):
    assert K>=6 and 0<alpha<1
    m,h,b,c,a,r,q=anchor(K)
    g=[G(F(i,K))-G(F(i-1,K)) for i in range(1,K+1)]
    p=[F(2*K-2*i+1,K*K) for i in range(1,K+1)]
    limit=min(H,max_full) if max_full is not None else H
    bound=int(F(81*K,80)/alpha)+1
    t_star=t_conf=None; expected_x=None; upper_seen=False
    s_list=[]; closure_steps=0
    for t in range(limit):
        pref=F(0); D=[]
        for j,v in enumerate(p,1):
            pref+=v; D.append(pref-F(j,K))
        d=max(D); s=D.index(d)+1; s_list.append(s)
        assert min(p)>=0 and pref==1 and sum(p[:m])==c
        assert s<=m and d>=a>TOL and all(z<a for z in D[m:])
        if t_conf is None and all(z<a for z in D[:m-2]): t_conf=t
        if t_star is None and s in (m-1,m): t_star=t
        if t_star is not None:
            assert s in (m-1,m), ("interior return",K,str(alpha),t,s)
            x=p[m-1]
            if expected_x is not None: assert x==expected_x
            assert s==(m if x>h else m-1)
            assert d==a+max(h-x,F(0))
            expected_x=(1-alpha)*x+alpha*c*q if x>h else x+alpha*(c-x)*q
            closure_steps+=1
            if r<3:
                if upper_seen: assert s==m
                upper_seen=upper_seen or s==m
        mass=alpha*sum(p[:s])
        p=[(1-alpha)*v+mass*g[i] if i<s else v+mass*g[i] for i,v in enumerate(p)]
    assert t_star is not None and t_conf is not None and t_conf<=bound and t_star<=bound
    if r<3: assert upper_seen
    label=classify(s_list,m,K)
    return {"K":K,"alpha":str(alpha),"requested_H":H,"verified_H":limit,
            "t_star":t_star,"t_conf":t_conf,"T_alpha":bound,"closure_states":closure_steps,
            "finite_classifier":label,
            "both_branches_in_verified_window":set(s_list[t_star:])=={m-1,m}}

def finite_float_proxy(K,alpha,N=200000):
    _,h,_,c,_,r,q=anchor(K)
    h,c,q,alpha=map(float,(h,c,q,alpha))
    x=h+alpha*(c-h)*q if r<3 else 0.0
    up=down=0
    for _ in range(N):
        before=x>h
        x=(1-alpha)*x+alpha*c*q if before else x+alpha*(c-x)*q
        after=x>h
        up+=not before and after; down+=before and not after
        if r<3: assert after
    if r>=3: assert up>=2 and down>=2
    return {"K":K,"alpha":str(alpha),"steps":N,"upcrossings":up,"downcrossings":down}

def main():
    started=time.time()
    alphas=[F(1,8),F(1,4),F(3,8),F(1,2),F(3,4)]
    census=list(range(6,54))+[100,200,400]
    matrix=[]
    for alpha in alphas:
        for K in census:
            matrix.append(run_alpha(K,alpha,10*K,max_full=300 if K>60 else None))
        print("exact alpha matrix block PASS:",str(alpha),flush=True)
    extra=[]
    for alpha in [F(1,1000),F(1,100),F(4,27),F(6,31),F(1,3),F(99,100),F(999,1000)]:
        extra.append(run_alpha(6,alpha,400))
    proxies=[]
    for alpha in alphas:
        for K in census: proxies.append(finite_float_proxy(K,alpha))
        print("finite float block completed 51 x 200000 steps:",str(alpha),flush=True)
    counter=run_alpha(9,F(1,8),90)
    assert counter["finite_classifier"]["class"]=="LOCK"
    assert counter["finite_classifier"]["R_m"]==10 and counter["both_branches_in_verified_window"]
    print("Exact finite-label counterexample K=9 alpha=1/8 H=90:",counter["finite_classifier"])
    print("PASS: exact matrix=255; extra K6=7; float runs=255, actual steps=51000000")
    print("Float runs are finite support; universal quantifiers come from the proof.")
    return {"status":"PASS","matrix":matrix,"K6_extra":extra,"finite_float":proxies,
            "finite_counterexample":counter,"elapsed_seconds":round(time.time()-started,3)}

if __name__=="__main__":
    main()
