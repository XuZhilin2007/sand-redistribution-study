"""Audit-only reconstruction. No project modules imported."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json

def weights(k):
    segments = [(F(0), F(1,3), F(1,2)), (F(1,3), F(2,3), F(3,2)), (F(2,3), F(5,6), F(2))]
    return [sum((max(F(0), min(F(j,k),b)-max(F(j-1,k),a))*v for a,b,v in segments),F(0)) for j in range(1,k+1)]

def trajectory(k,a,word=None,n=4):
    p=[F(2*k-2*j+1,k*k) for j in range(1,k+1)]; g=weights(k)
    states=[p]; selectors=[]; ties=[]
    for t in range(n):
        d=[sum(p[:j])-F(j,k) for j in range(1,k+1)]
        winners=[j+1 for j,x in enumerate(d) if x==max(d)]
        s=winners[0] if word is None else word[t]
        mass=sum(p[:s]); p=[x+a*(mass*g[j]-(x if j<s else 0)) for j,x in enumerate(p)]
        assert sum(p)==1
        states.append(p); selectors.append(s); ties.append(winners)
    return [k*sum(x*x for x in p)-1 for p in states],selectors,ties

def add(a,b):
    return [(a[j] if j<len(a) else 0)+(b[j] if j<len(b) else 0) for j in range(max(len(a),len(b)))]

def times(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c

def interpolate(k,word,t):
    # Degree <= 2*len(word) follows inductively from the affine-in-alpha update.
    n=2*len(word); nodes=[F(j,n+2) for j in range(n+1)]
    vals=[]
    for a in nodes:
        u,_,_=trajectory(k,a,word,len(word)); vals.append(u[t]-u[2])
    result=[F(0)]; basis=[F(1)]
    for order in range(n+1):
        result=add(result,[vals[0]*x for x in basis])
        vals=[(vals[j+1]-vals[j])/(nodes[j+order+1]-nodes[j]) for j in range(len(vals)-1)]
        basis=times(basis,[-nodes[order],F(1)])
    while result and result[-1]==0: result.pop()
    return result

def evaluate(p,x):
    return sum(v*x**j for j,v in enumerate(p))

def bernstein(p,lo,hi):
    n=len(p)-1
    shifted=[sum(p[i]*comb(i,j)*lo**(i-j)*(hi-lo)**j for i in range(j,n+1)) for j in range(n+1)]
    return [sum(shifted[j]*F(comb(k,j),comb(n,j)) for j in range(k+1)) for k in range(n+1)]

polys={'lower':[-37296000,272381200,-613816360,555506236,-240679790,42889525],
       'upper':[-1136,5068,-4361,1164],
       'a6':[-2080,17183,-38636,32735,-13378,2328]}
cases=[('lower',10,[5,7,8],3,45000000),('upper',6,[3,4],1,-1296),('a6',6,[3,4,5],3,2592)]
records={}
for name,k,word,t,scale in cases:
    result=interpolate(k,word,t)
    assert result==[F(0)]+[F(c,scale) for c in polys[name]]
    lo,hi={'lower':(F(6,25),F(1,4)), 'upper':(F(19,100),F(1,3)), 'a6':(F(19,100),F(1,5))}[name]
    derivative=[j*v for j,v in enumerate(polys[name])][1:]
    coefficients=bernstein(derivative,lo,hi)
    assert min(coefficients)>0
    assert evaluate(polys[name],lo)<0<evaluate(polys[name],hi)
    for _ in range(120):
        mid=(lo+hi)/2
        if evaluate(polys[name],mid)<0: lo=mid
        else: hi=mid
    records[name]={'loss_coefficients':list(map(str,result)), 'positive_derivative_bernstein':list(map(str,coefficients)), 'bracket':[str(lo),str(hi)]}

witnesses=[]
for k,a,t,expected in [(13,F(13,40),1,F(-31191,22497280)),(10,F(9,10),4,F(-4697222139,20000000000000)),(9,F(1,5),3,F(-83303117,5978711250)),(6,F(7,24),1,F(-161,11943936))]:
    u,word,ties=trajectory(k,a)
    assert u[t]-u[2]==expected
    if k==10: assert all(u[j]>u[2] for j in (0,1,3))
    witnesses.append({'K':k,'alpha':str(a),'word':word,'ties':ties,'losses':list(map(str,u)),'gap':str(expected)})
u,_,_=trajectory(13,F(13,40),[6,10],2)
v,_,_=trajectory(13,F(13,40),None,2)
assert u[2]-v[2]==F(-155929,16224000)
report={'status':'PASS','method':'Independent kernel intersection lengths, exact interpolation with analytic degree bound, Bernstein derivative coefficients, exact full-vector updates; no project imports','polynomials':records,'witnesses':witnesses,'K13_jump':str(u[2]-v[2])}
Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS: three polynomial identities and root uniqueness; four exact witnesses; K13 jump')
