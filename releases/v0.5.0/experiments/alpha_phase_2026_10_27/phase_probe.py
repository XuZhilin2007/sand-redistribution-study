"""Exploration only: finite K/time phase census, separate from v0.4.0.

Standard library, full vectors until first boundary selection, then the exact
two-coordinate representation evaluated in float. No float PASS is a proof.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from math import lcm
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent


def cdf(x):
    if x <= Q(1, 3):
        return x / 2
    if x <= Q(2, 3):
        return 3 * x / 2 - Q(1, 3)
    if x <= Q(5, 6):
        return 2 * x - Q(2, 3)
    return Q(1)


@lru_cache(maxsize=None)
def model(k, exact=False):
    h = Q(1, k)
    p = [Q(2 * k - 2 * i + 1, k * k) for i in range(1, k + 1)]
    g = [cdf(i * h) - cdf((i - 1) * h) for i in range(1, k + 1)]
    return (p, g, h) if exact else ([float(x) for x in p], [float(x) for x in g], float(h))


@lru_cache(maxsize=None)
def first_step_lines(k):
    """Integer coefficients of D1, so odd-K initial ties are never rounded."""
    p, g, h = model(k, True)
    denominator = lcm(*(v.denominator for v in g))
    selected = k // 2
    mass = selected * (2 * k - selected)
    cumulative = Q(0)
    lines = []
    for j, gj in enumerate(g, 1):
        cumulative += gj
        base = j * (k - j)
        slope = mass * int(cumulative * denominator)
        slope -= j * (2 * k - j) * denominator if j <= selected else mass * denominator
        lines.append((base * denominator, slope))
    return lines


def select(p, h):
    acc, best, selected = 0, -1, 0
    for j, x in enumerate(p, 1):
        acc += x
        d = acc - j * h
        if d > best:
            best, selected = d, j
    return selected


def step(p, g, alpha, s):
    mass = alpha * sum(p[:s])
    return [(1 - alpha) * x + mass * y if i < s else x + mass * y
            for i, (x, y) in enumerate(zip(p, g))]


def loss(p, h):
    return len(p) * sum((x - h) ** 2 for x in p)


def replay(k, alpha, horizon):
    p, g, h = model(k, exact=isinstance(alpha, Q))
    states, words = [p], []
    for _ in range(horizon):
        s = select(p, h)
        p = step(p, g, alpha, s)
        words.append(s)
        states.append(p)
    return states, words, [loss(p, h) for p in states]


def probe(k, alpha, horizon=400):
    rational_alpha = Q(alpha)
    alpha = float(alpha)
    p, g, h = model(k)
    m = (5 * k + 5) // 6
    b = m * h
    c = b * (2 - b)
    q = g[m - 1]
    rho = 1 - alpha
    values = [loss(p, h)]
    words, t0 = [], None
    for t in range(horizon + 1):
        if t == 0:
            s = k // 2
        elif t == 1:
            aa, ab = rational_alpha.numerator, rational_alpha.denominator
            lines = first_step_lines(k)
            s = 1 + max(range(k), key=lambda i: (ab * lines[i][0] + aa * lines[i][1], -i))
        else:
            s = select(p, h)
        if t >= 2 and (min(values[:2]) < values[2] - 1e-12 or
                       (t >= 3 and values[-1] < values[2] - 1e-12)):
            return False, min(range(len(values)), key=values.__getitem__), words[:5], t0
        if s in (m - 1, m) and t >= 2:
            t0 = t
            x = p[m - 1]
            gamma = (c - x) / (1 - q)
            e = [p[i] - gamma * g[i] for i in range(m - 1)]
            base = [c * y / (1 - q) for y in g[:m - 1]] + [0.] + p[m:]
            direction = [-y / (1 - q) for y in g[:m - 1]] + [1.] + [0.] * (k - m)
            e += [0.] * (k - m + 1)
            a = k * sum(y * y for y in direction)
            bb = 2 * k * sum((y - h) * v for y, v in zip(base, direction))
            cc = loss(base, h)
            dd = k * sum(v * w for v, w in zip(direction, e))
            ff = k * sum((y - h) * w for y, w in zip(base, e))
            ee = k * sum(w * w for w in e)
            z = 1.
            for j in range(t, horizon):
                st = m if x > h else m - 1
                x = rho * x + alpha * c * q if st == m else (1 - alpha * q) * x + alpha * c * q
                z *= rho
                val = a * x * x + bb * x + cc + 2 * z * (dd * x + ff) + ee * z * z
                values.append(val)
                if len(words) < 5:
                    words.append(st)
                if val < values[2] - 1e-12:
                    return False, j + 1, words[:5], t0
            return True, 2, words[:5], t0
        if t == horizon:
            break
        p = step(p, g, alpha, s)
        words.append(s)
        values.append(loss(p, h))
    return True, min(range(len(values)), key=values.__getitem__), words[:5], t0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-k', type=int, default=120)
    parser.add_argument('--denominator', type=int, default=1000)
    parser.add_argument('--horizon', type=int, default=400)
    parser.add_argument('--output', type=Path, default=HERE / 'finite_phase.json')
    args = parser.parse_args()
    started = time.monotonic()
    rows = []
    for k in range(6, args.max_k + 1):
        good = []
        witnesses = {}
        for j in range(1, args.denominator):
            ok, t, word, entry = probe(k, Q(j, args.denominator), args.horizon)
            if ok:
                good.append(j)
            elif t not in witnesses:
                witnesses[t] = {'alpha': f'{j}/{args.denominator}', 'word': word, 'entry': entry}
        runs = []
        for j in good:
            if not runs or j > runs[-1][1] + 1:
                runs.append([j, j])
            else:
                runs[-1][1] = j
        rows.append({'K': k, 'sampled_good_runs_numerators': runs, 'first_failure_by_time': witnesses})
        if k <= 20 or len(runs) != 1 or k % 25 == 0:
            print(k, runs, 'failure times', sorted(witnesses), flush=True)
    report = {'evidence': 'FLOAT FINITE K FINITE TIME ONLY', 'denominator': args.denominator,
              'prefix_selection': 's0 and s1 exact; later selectors float; near ties require exact replay',
              'horizon': args.horizon, 'K_domain': [6, args.max_k], 'rows': rows,
              'elapsed_seconds': time.monotonic() - started,
              'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('wrote', args.output, 'seconds', report['elapsed_seconds'])


if __name__ == '__main__':
    main()
