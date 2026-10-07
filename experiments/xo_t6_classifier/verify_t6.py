"""Exact finite checks corresponding to the repaired T6 proof, alpha=1/4."""
from fractions import Fraction
import time
F=Fraction

def G(y):
    if y<=F(1,3): return y/2
    if y<=F(2,3): return 3*y/2-F(1,3)
    if y<=F(5,6): return 2*y-F(2,3)
    return F(1)

def anchor(K):
    m=(5*K+5)//6;h=F(1,K);b=F(m,K)
    return m,h,b,2*b-b*b,b*(1-b),K%6,F(6-K%6,3*K)

def power_bound(factor,ratio,strict=False):
    power=F(1)
    for n in range(1,10001):
        power*=factor
        if power<ratio if strict else power<=ratio: return n
    raise AssertionError("bound not reached")

def maxA_bound(K):
    _,h,_,c,_,_,q=anchor(K);L=c*q
    if L>=h:return None
    return power_bound(F(3,4),(h-L)/(h-F(3,4)*L))

def maxB_post_bound(K):
    _,h,_,c,_,_,q=anchor(K);L=c*q
    if L>=h:return None
    return power_bound(1-q/4,(c-h)/(c-(3*h+L)/4))

def maxB_global_bound(K):
    _,h,_,c,_,_,q=anchor(K)
    return power_bound(1-q/4,(c-h)/c,strict=True)

def scalar_classifier(K,H=None,entry_cap=100000):
    assert K>=6
    if H is None:H=10*K
    assert H>0
    m,h,b,c,a,r,q=anchor(K)
    p=[F(2*K-2*i+1,K*K) for i in range(1,K+1)]
    g=[G(F(i,K))-G(F(i-1,K)) for i in range(1,K+1)]
    ss=[];xs=[];first=None
    for t in range(min(H,entry_cap)):
        pref=F(0);D=[]
        for j,v in enumerate(p,1):pref+=v;D.append(pref-F(j,K))
        s=D.index(max(D))+1
        if s in (m-1,m):first=t;x=p[m-1];break
        ss.append(s);xs.append(p[m-1])
        mass=sum(p[:s])/4
        p=[F(3,4)*v+mass*g[i] if i<s else v+mass*g[i] for i,v in enumerate(p)]
        assert min(p)>=0 and sum(p)==1 and sum(p[:m])==c
    if first is None:
        assert len(ss)==H,"entry cap shorter than requested window"
    else:
        for t in range(first,H):
            ss.append(m if x>h else m-1);xs.append(x)
            x=F(3,4)*x+c*q/4 if x>h else x+(c-x)*q/4
    assert len(ss)==H
    tail=0
    for s in reversed(ss):
        if s!=m:break
        tail+=1
    post=set(ss[first:]) if first is not None else set()
    label="LOCK" if tail>=K else "CYCLE" if post=={m-1,m} else "OTHER"
    switches=sum(ss[t]==m and ss[t+1]==m-1 for t in range(H-1))
    rebounds=sum(ss[t]==m and xs[t+1]<h for t in range(H-1))
    return {"K":K,"H":H,"t_star":first,"R_m":tail,"class":label,
            "A_to_B_switches":switches,"strict_rebounds":rebounds,
            "s_list":ss,"x_list":xs}

def dwell_summary(result):
    K=result["K"];m=anchor(K)[0];first=result["t_star"]
    post=result["s_list"][first:];runs=[]
    for s in post:
        if runs and runs[-1][0]==s:runs[-1][1]+=1
        else:runs.append([s,1])
    A=[n for j,(s,n) in enumerate(runs) if s==m and j>0 and runs[j-1][0]==m-1]
    B=[n for j,(s,n) in enumerate(runs) if s==m-1 and j>0 and runs[j-1][0]==m]
    return {"K":K,"initial_branch":runs[0][0],"initial_dwell":runs[0][1],
            "maxA_after_B":max(A,default=0),"maxB_after_A":max(B,default=0),
            "boundA_after_B":maxA_bound(K),"boundB_after_A":maxB_post_bound(K)}

def main():
    started=time.time()
    checks={
      "lambda max r3":anchor(9)[3]==F(80,81),
      "lambda max r4":F(2,3)*anchor(10)[3]<=F(80,81),
      "lambda max r5":anchor(11)[3]/3<=F(80,81),
      "ratio at lambda max":(1-F(80,81))/(1-F(3,4)*F(80,81))==F(1,21),
      "11-step post-lower A bound":F(3,4)**11<=F(1,21),
      "initial upper 3K base":81*9*F(3,4)**27<1,
      "initial upper monotone ratio":F(10,9)*F(27,64)==F(15,32)<1,
      "lower bound <=3K at min cycle K":F(72,29)*9+1<=27,
      "lower bound slope":F(72,29)<3,
      "LOCK corrected count at K6":F(2011,580)*6-2>=6,
      "LOCK margin slope":F(2011,580)>1,
      "CYCLE corrected chain at K9":F(141,20)*9+12<89,
      "CYCLE chain slope":F(141,20)<10,
      "small A10":maxA_bound(10)==2,
      "small A11":maxA_bound(11)==1,
    }
    for name,ok in checks.items():
        assert ok,name
        print("PASS closed-form check:",name,flush=True)
    records=[];dwells=[]
    for K in range(6,301):
        result=scalar_classifier(K)
        assert result["class"]==("LOCK" if K%6<3 else "CYCLE"),(K,result["class"])
        if K==9:assert result["R_m"]==3 and result["t_star"]==2
        if K==6:assert result["R_m"]==58 and result["t_star"]==1
        if K%6>=3:assert result["A_to_B_switches"]>0 and result["R_m"]<K
        if K in (9,10,11,15,16,17,21,22,23):
            row=dwell_summary(result)
            assert row["maxA_after_B"]<=row["boundA_after_B"]
            assert row["maxB_after_A"]<=row["boundB_after_A"]
            dwells.append(row);print("Computed dwell:",row,flush=True)
        records.append({key:v for key,v in result.items() if key not in ("s_list","x_list")})
    print("PASS: classifier examples K=6..300: 295/295; all checks assert on failure.")
    return {"status":"PASS","closed_form_checks":checks,"classifications":records,
            "dwell_examples":dwells,"elapsed_seconds":round(time.time()-started,3)}

if __name__=="__main__":main()
