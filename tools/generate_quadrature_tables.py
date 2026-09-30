#!/usr/bin/env python3
"""Generate every node and weight table embedded in the LNIT headers.

Tables are computed with 80 digits and written as long double literals with 36
significant digits, enough for the IEEE quad type (33 digits). A long double
literal converted to Scalar is correctly rounded for float, double and every
long double format, so all the instantiations get their full precision.

Rules:
  - closed Clenshaw-Curtis with n+1 nodes x_k = cos(k pi / n), n even:
        w_k = (c_k / n) (1 - 2 sum_{j=1}^{n/2} b_j cos(2 j k pi / n) / (4 j^2 - 1)),
    c_0 = c_n = 1, c_k = 2 otherwise, b_j = 1 except b_{n/2} = 1/2;
  - Gauss-Legendre with 15 nodes;
  - error-estimation companions, interpolatory on a subset of the nodes of their
    main rule (weights solve sum_i w_i x_i^k = int_{-1}^{1} x^k for k < number of nodes):
    14 and 6 nodes of Gauss-Legendre 15, 10 nodes of Clenshaw-Curtis 13;
  - Gauss-Laguerre with 33 nodes and Gauss-Hermite with 66 nodes, whose tables hold
    the weights multiplied by exp(x) and exp(x^2), because the integrand is passed
    without its weight function.

Usage:
    python3 tools/generate_quadrature_tables.py            # print the C++ arrays
    python3 tools/generate_quadrature_tables.py --update   # rewrite the headers in place

Requires mpmath.
"""
import argparse
import pathlib
import re
import sys

from mpmath import mp, mpf, cos, exp, pi, sqrt, factorial, gamma, nstr, lu_solve, matrix, fabs

mp.dps = 80
DIGITS = 36
PER_LINE = 3

INCLUDE = pathlib.Path(__file__).resolve().parent.parent / "include" / "LNIT"
ADAPTIVE = INCLUDE / "AdaptiveQuadratures"


# --- rules -------------------------------------------------------------------

def clenshaw_curtis(n):
    """nodes (decreasing, from 1 to -1) and weights of the closed rule with n+1 nodes"""
    nodes = [cos(k * pi / n) for k in range(n + 1)]
    nodes[n // 2] = mpf(0)
    weights = []
    for k in range(n + 1):
        s = mpf(0)
        for j in range(1, n // 2 + 1):
            b = mpf(1) / 2 if 2 * j == n else mpf(1)
            s += b / (4 * j * j - 1) * cos(2 * j * k * pi / n)
        c = 1 if k % n == 0 else 2
        weights.append(mpf(c) / n * (1 - 2 * s))
    return nodes, weights


def legendre(n, x):
    """P_n(x) and P_n'(x)"""
    p0, p1 = mpf(1), x
    for k in range(2, n + 1):
        p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
    return p1, n * (x * p1 - p0) / (x * x - 1)


def gauss_legendre(n):
    """nodes (increasing) and weights"""
    nodes, weights = [], []
    for i in range(n):
        x = -cos(pi * (i + mpf(3) / 4) / (n + mpf(1) / 2))
        for _ in range(100):
            p, dp = legendre(n, x)
            dx = p / dp
            x -= dx
            if fabs(dx) < mpf(10) ** (-mp.dps + 5):
                break
        _, dp = legendre(n, x)
        nodes.append(x)
        weights.append(2 / ((1 - x * x) * dp * dp))
    if n % 2 == 1:
        nodes[n // 2] = mpf(0)
    return nodes, weights


def interpolatory(nodes, subset):
    """weights of the interpolatory rule on nodes[i], i in subset; zero elsewhere"""
    xs = [nodes[i] for i in subset]
    m = len(xs)
    A = matrix(m, m)
    b = matrix(m, 1)
    for k in range(m):
        for j, x in enumerate(xs):
            A[k, j] = x ** k
        b[k] = mpf(0) if k % 2 else mpf(2) / (k + 1)
    w = lu_solve(A, b)
    weights = [mpf(0)] * len(nodes)
    for j, i in enumerate(subset):
        weights[i] = w[j]
    return weights


def newton_roots(poly, guesses):
    """refine the roots of poly (returns value and derivative) from double guesses"""
    roots = []
    for x in guesses:
        x = mpf(x)
        for _ in range(200):
            p, dp = poly(x)
            dx = p / dp
            x -= dx
            if fabs(dx) <= fabs(x) * mpf(10) ** (-mp.dps + 5):
                break
        roots.append(x)
    return roots


def laguerre(n, x):
    """L_n(x) and L_n'(x)"""
    l0, l1 = mpf(1), 1 - x
    for k in range(2, n + 1):
        l0, l1 = l1, ((2 * k - 1 - x) * l1 - (k - 1) * l0) / k
    return l1, n * (l1 - l0) / x


def gauss_laguerre(n, guesses):
    """nodes and weights times exp(x)"""
    nodes = newton_roots(lambda x: laguerre(n, x), guesses)
    weights = []
    for x in nodes:
        ln1, _ = laguerre(n + 1, x)
        weights.append(x / ((n + 1) ** 2 * ln1 * ln1) * exp(x))
    return nodes, weights


def hermite(n, x):
    """physicists' H_n(x) and H_n'(x) = 2 n H_{n-1}(x)"""
    h0, h1 = mpf(1), 2 * x
    for k in range(2, n + 1):
        h0, h1 = h1, 2 * x * h1 - 2 * (k - 1) * h0
    return h1, 2 * n * h0


def gauss_hermite(n, guesses):
    """nodes and weights times exp(x^2)"""
    nodes = newton_roots(lambda x: hermite(n, x), guesses)
    weights = []
    for x in nodes:
        hn1, _ = hermite(n - 1, x)
        weights.append(2 ** (n - 1) * factorial(n) * sqrt(pi) / (n * n * hn1 * hn1) * exp(x * x))
    return nodes, weights


# --- C++ output --------------------------------------------------------------

ARRAY = re.compile(r"\tstatic constexpr (?:inline )?std::array<Scalar, (\d+)>\s*(\w+)\s*=?\s*\{.*?\};", re.S)


def literal(x):
    if x == 0:
        return "Scalar{}"
    return f"Scalar({nstr(x, DIGITS, strip_zeros=False, min_fixed=-5, max_fixed=5)}L)"


def cxx_array(name, values):
    items = [literal(v) + ("," if i + 1 < len(values) else "};") for i, v in enumerate(values)]
    width = max(len(s) for s in items)
    lines = [f"\tstatic constexpr std::array<Scalar, {len(values)}> {name} = {{"]
    for start in range(0, len(items), PER_LINE):
        chunk = items[start:start + PER_LINE]
        lines.append("\t\t" + " ".join(s.ljust(width) for s in chunk).rstrip())
    return "\n".join(lines)


def current_values(path, name):
    """the double values of a table as currently written in a header"""
    for match in ARRAY.finditer(path.read_text()):
        if match.group(2) == name:
            body = match.group(0).split("{", 1)[1]
            out = []
            for item in re.findall(r"(-?)\s*Scalar\s*[({]\s*([-+0-9.eE]*)\s*L?\s*[)}]", body):
                sign, number = item
                value = mpf(number) if number else mpf(0)
                out.append(-value if sign else value)
            return out
    sys.exit(f"{path}: no table {name}")


def update_header(path, name, text):
    source = path.read_text()
    matches = [m for m in ARRAY.finditer(source) if m.group(2) == name]
    if len(matches) != 1:
        sys.exit(f"{path}: expected exactly one definition of {name}")
    m = matches[0]
    path.write_text(source[:m.start()] + text + source[m.end():])


# --- exactness checks (independent of the previous tables) -------------------

def legendre_exactness(nodes, weights, degree):
    """largest error of the rule on x^k over [-1, 1], k <= degree"""
    err = mpf(0)
    for k in range(degree + 1):
        exact = mpf(0) if k % 2 else mpf(2) / (k + 1)
        err = max(err, fabs(sum(w * x ** k for x, w in zip(nodes, weights)) - exact))
    return err


def laguerre_exactness(nodes, weights, degree):
    """largest relative error on int_0^inf x^k exp(-x) = k!, k <= degree"""
    return max(fabs(sum(w * exp(-x) * x ** k for x, w in zip(nodes, weights)) / factorial(k) - 1)
               for k in range(degree + 1))


def hermite_exactness(nodes, weights, degree):
    """largest relative error on int x^k exp(-x^2) = Gamma((k+1)/2) (even k), 0 (odd k), k <= degree;
    for odd k, relative to the sum of the absolute values of the terms"""
    err = mpf(0)
    for k in range(degree + 1):
        terms = [w * exp(-x * x) * x ** k for x, w in zip(nodes, weights)]
        s = sum(terms)
        err = max(err, fabs(s) / sum(fabs(t) for t in terms) if k % 2 else fabs(s / gamma(mpf(k + 1) / 2) - 1))
    return err


def check_exactness():
    tolerance = mpf(10) ** -60
    checks = []
    for n in (8, 12, 16, 32):
        x, w = clenshaw_curtis(n)
        checks.append((f"Clenshaw-Curtis {n + 1}", legendre_exactness(x, w, n + 1)))
    x13, _ = clenshaw_curtis(12)
    sub = [i for i in range(13) if i not in (0, 6, 12)]
    checks.append(("Clenshaw-Curtis 13, 10-node companion", legendre_exactness(x13, interpolatory(x13, sub), 9)))
    x15, w15 = gauss_legendre(15)
    checks.append(("Gauss-Legendre 15", legendre_exactness(x15, w15, 29)))
    checks.append(("Gauss-Legendre 15, 14-node companion", legendre_exactness(x15, interpolatory(x15, [i for i in range(15) if i != 7]), 13)))
    checks.append(("Gauss-Legendre 15, 6-node companion", legendre_exactness(x15, interpolatory(x15, [1, 3, 5, 9, 11, 13]), 5)))
    return checks, tolerance


def tables():
    cc9 = clenshaw_curtis(8)
    cc13 = clenshaw_curtis(12)
    cc17 = clenshaw_curtis(16)
    cc33 = clenshaw_curtis(32)
    gl15 = gauss_legendre(15)
    cc_h = ADAPTIVE / "ClenshawCurtisAdaptiveQuadrature.hpp"
    hyb_h = ADAPTIVE / "ClenshawCurtisHybridAdaptiveQuadrature.hpp"
    gl_h = ADAPTIVE / "GaussLegendreAdaptiveQuadrature.hpp"
    glcc_h = ADAPTIVE / "GLCCAdaptiveQuadrature.hpp"
    lag_h = INCLUDE / "GaussLaguerreQuadrature.hpp"
    her_h = INCLUDE / "GaussHermiteQuadrature.hpp"
    lag = gauss_laguerre(33, current_values(lag_h, "s_xi"))
    her = gauss_hermite(66, current_values(her_h, "s_xi"))
    return [
        (cc_h, "s_wi09", cc9[1]),
        (cc_h, "s_wi17", cc17[1]),
        (cc_h, "s_wi33", cc33[1]),
        (cc_h, "s_xi", cc33[0]),
        (hyb_h, "s_wi", cc13[1]),
        (hyb_h, "s_alternateWi", interpolatory(cc13[0], [i for i in range(13) if i not in (0, 6, 12)])),
        (hyb_h, "s_xi", cc13[0]),
        (gl_h, "s_xi", gl15[0]),
        (gl_h, "s_wi15", gl15[1]),
        (gl_h, "s_wi14", interpolatory(gl15[0], [i for i in range(15) if i != 7])),
        (gl_h, "s_wi06", interpolatory(gl15[0], [1, 3, 5, 9, 11, 13])),
        (glcc_h, "s_xi_gl", gl15[0]),
        (glcc_h, "s_wi_gl", gl15[1]),
        (glcc_h, "s_xi_cc", cc33[0]),
        (glcc_h, "s_wi_cc", cc33[1]),
        (lag_h, "s_wi", lag[1]),
        (lag_h, "s_xi", lag[0]),
        (her_h, "s_wi", her[1]),
        (her_h, "s_xi", her[0]),
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--update", action="store_true", help="rewrite the headers in place")
    args = parser.parse_args()

    checks, tolerance = check_exactness()
    lag = gauss_laguerre(33, current_values(INCLUDE / "GaussLaguerreQuadrature.hpp", "s_xi"))
    her = gauss_hermite(66, current_values(INCLUDE / "GaussHermiteQuadrature.hpp", "s_xi"))
    checks.append(("Gauss-Laguerre 33", laguerre_exactness(*lag, 65)))
    checks.append(("Gauss-Hermite 66", hermite_exactness(*her, 131)))
    for rule, err in checks:
        if not err < tolerance:
            sys.exit(f"{rule}: not exact to its degree (error {nstr(err, 3)})")
        print(f"// {rule}: exact to its degree (error {nstr(err, 3)})")

    for path, name, values in tables():
        # The new tables must agree with the current ones to about their double precision:
        # this checks that each table is the rule the header implements. The margin is
        # loose because some previous tables were computed in double precision.
        old = current_values(path, name)
        if len(old) != len(values):
            sys.exit(f"{path.name} {name}: {len(old)} values in the header, {len(values)} generated")
        scale = max(fabs(v) for v in values)
        diff = max(fabs(a - b) for a, b in zip(old, values)) / scale
        if diff > mpf("1e-13"):
            sys.exit(f"{path.name} {name}: generated table differs from the header by {nstr(diff, 3)}")
        text = cxx_array(name, values)
        if args.update:
            update_header(path, name, text)
            print(f"updated {path.name}: {name} (previous table off by {nstr(diff, 3)})")
        else:
            print(f"// {path.name} (previous table off by {nstr(diff, 3)})")
            print(text)
            print()


if __name__ == "__main__":
    main()
