"""Exact first-step candidate verification; corrected s0/s1/s2 indexing."""
from fractions import Fraction as F
import time

def G(y):
    if y<=F(1,3):return y/2
    if y<=F(2,3):return 3*y/2-F(1,3)
    if y<=F(5,6):return 2*y-F(2,3)
    return F(1)

def d1_exactly(K,j):
    u=F(K//2,K);f0=2*u-u*u;y=F(j,K);d0=y*(1-y)
    return F(3,4)*d0+(f0*G(y)-y)/4 if j<=K//2 else d0-f0*(1-G(y))/4

def integer_score(K,j):
    # Independent common denominator 24 K^3.
    s0=K//2;f0num=s0*(2*K-s0)
    if 3*j<=K:gnum=3*j
    elif 3*j<=2*K:gnum=9*j-2*K
    elif 6*j<=5*K:gnum=12*j-4*K
    else:gnum=6*K
    if j<=s0:return 18*K*j*(K-j)+f0num*gnum-6*j*K*K
    return 24*K*j*(K-j)-f0num*(6*K-gnum)

def argmax_D1(K):
    bestj=1;best=integer_score(K,1)
    for j in range(2,K+1):
        value=integer_score(K,j)
        if value>best:best,bestj=value,j
    return bestj

def formula_candidates(K):
    u=F(K//2,K);f0=2*u-u*u
    points=[K*(1+3*f0/8)/2,K*(1+f0/2)/2,F(2*K,3)]
    result=set()
    for point in points:
        lo=point.numerator//point.denominator
        hi=-(-point.numerator//point.denominator)
        result.update(j for j in (lo,hi) if K//2<j<=K)
    return sorted(result)

def formula_argmax(K):
    return min(formula_candidates(K),key=lambda j:(-d1_exactly(K,j),j))

def first_selections(K):
    p=[F(2*K-2*i+1,K*K) for i in range(1,K+1)]
    g=[G(F(i,K))-G(F(i-1,K)) for i in range(1,K+1)];ss=[]
    for t in range(3):
        pref=F(0);D=[]
        for j,v in enumerate(p,1):pref+=v;D.append(pref-F(j,K))
        s=D.index(max(D))+1;ss.append(s)
        if t==2:break
        mass=sum(p[:s])/4
        p=[F(3,4)*v+mass*g[i] if i<s else v+mass*g[i] for i,v in enumerate(p)]
        assert min(p)>=0 and sum(p)==1
    return ss

def s2_of(K):return first_selections(K)[2]

def main():
    started=time.time()
    for K in range(6,4001):
        assert argmax_D1(K)==formula_argmax(K),K
        assert len(formula_candidates(K))<=6
    records=[]
    for K in range(6,61):
        s0,s1,s2=first_selections(K)
        assert s0==K//2 and s1==formula_argmax(K)
        records.append({"K":K,"s0":s0,"s1":s1,"s2":s2,
                        "s1_over_K":str(F(s1,K)),"s2_over_s1":str(F(s2,s1))})
    assert first_selections(50)==[25,32,41] and s2_of(50)==41
    f0=F(3,4);y2=(1+3*f0/8)/2;y3=(1+f0/2)/2
    v2=(1+3*f0/8)**2/4-f0/3;v3=(1+f0/2)**2/4-5*f0/12
    assert y2==F(41,64) and y3==F(11,16) and v2-v3==F(1,4096)
    print("PASS: first-step candidates vs exact full-bin argmax K=6..4000: 3995/3995")
    print("K=50 corrected indices: s0=25, s1=32, s2=41; s1/K=16/25; s2/s1=41/32")
    print("Corrected finite second-step records (no asymptotic claim):",records)
    return {"status":"PASS","first_step_checked":3995,"second_step_records":records,
            "elapsed_seconds":round(time.time()-started,3)}

if __name__=="__main__":main()
