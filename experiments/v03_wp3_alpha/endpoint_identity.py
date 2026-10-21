"""Exact endpoint identity for K=10 (Theorem G-ep): with lambda = 1-a,
beta = 1 - a/15, mu = beta/lambda, H = (1+lam)(1+lam beta) + lam^2 beta^2,
P(a) = 17a^4 - 99a^3 - 2448a^2 + 8230a - 1875,

    delta_K(a) - delta^+(lambda, mu, 2/5)  ==  lambda * P(a) / (10500 * H)

as a POLYNOMIAL IDENTITY in a (verified by rational polynomial arithmetic:
10500*H*(delta_K*D - N) - lambda*P*D is the ZERO polynomial, where
delta^+ = N/D with the LN Theorem-1 formula, S = lambda^2 for p/q = 2/5).
Also verifies P(0) < 0 < P(1) and P' >= 3037 on [0,1] (exact grid + bound).
"""
from fractions import Fraction as F


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
            for i in range(n)]


def pmul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def psub(a, b):
    return padd(a, [-x for x in b])


def pscale(a, s):
    return [x * s for x in a]


def peval(a, x):
    r = F(0)
    for c in reversed(a):
        r = r * x + c
    return r


def ptrim(a):
    out = list(a)
    while out and out[-1] == 0:
        out.pop()
    return out


one = [F(1)]
lam = [F(1), F(-1)]          # lambda = 1 - a
beta = [F(1), F(-1, 15)]     # beta = 1 - a/15
eta = F(89, 140)

# delta^+ = N/D with mu = beta/lam and S = lam^2 (p/q = 2/5, floors verified):
#   N = (1-lam)(1 + mu*lam^2 + lam^4 mu^2 (1-lam))
#     = (1-lam)(1 + lam*beta + lam^2 beta^2 (1-lam))
#   D = 1 + (mu-1) lam^2 + lam^4 mu (mu - lam mu - 1)
#     = 1 + (beta-lam) lam + lam^2 beta (beta - lam beta - lam)
N1 = padd(one, pmul(lam, beta))
N2 = padd(one, pscale(lam, -1))
N3 = pmul(pmul(lam, lam), pmul(pmul(beta, beta), N2))
N = pmul(N2, padd(N1, N3))
Dinner = psub(psub(beta, pmul(lam, beta)), lam)
D = padd(padd(one, pmul(psub(beta, lam), lam)),
         pmul(pmul(lam, lam), pmul(beta, Dinner)))

delta_K = padd(one, pscale(lam, -eta))          # 1 - eta*lambda
lhs_num = psub(pmul(delta_K, D), N)             # (delta_K - delta^+) * D
P = [F(-1875), F(8230), F(-2448), F(-99), F(17)]
H = padd(pmul(padd(one, lam), padd(one, pmul(lam, beta))),
         pmul(pmul(lam, lam), pmul(beta, beta)))
residual = psub(pscale(pmul(H, lhs_num), 10500), pmul(pmul(lam, P), D))
assert not ptrim(residual), "identity fails"
assert peval(P, 0) < 0 < peval(P, 1)
Pd = [F(8230), F(-4896), F(-297), F(68)]
grid_min = min(peval(Pd, F(i, 100)) for i in range(101))
assert grid_min >= 3037, grid_min
# bracket check
assert peval(P, F(245996, 10 ** 6)) < 0 < peval(P, F(245997, 10 ** 6))
# numeric sanity at alpha = 1/4
a = F(1, 4)
lhs = peval(lhs_num, a) / peval(D, a)
rhs = peval(pscale(pmul(lam, P), 1), a) / (10500 * peval(H, a))
assert lhs == rhs
print("ENDPOINT IDENTITY PASS (zero residual; P(0)<0<P(1); "
      "P'>=3037 on [0,1] grid; bracket exact; alpha=1/4 sanity exact)")
