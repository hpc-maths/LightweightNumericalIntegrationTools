#!/usr/bin/env python3
"""Generate the Clenshaw-Curtis weight tables embedded in the LNIT headers.

Closed Clenshaw-Curtis rule with n+1 nodes x_k = cos(k pi / n), n even:

    w_k = (c_k / n) * (1 - 2 * sum_{j=1}^{n/2} b_j cos(2 j k pi / n) / (4 j^2 - 1))

with c_0 = c_n = 1, c_k = 2 otherwise, and b_j = 1 except b_{n/2} = 1/2.
The 1/2 factor on the last Chebyshev mode is what the tables shipped before
this script were missing (they were exact only up to degree n-1).

Usage:
    python3 tools/generate_clenshaw_curtis_tables.py            # print the C++ arrays
    python3 tools/generate_clenshaw_curtis_tables.py --update   # rewrite the headers in place

Requires mpmath. Values are computed with 50 digits and printed with 20
significant digits, like the other tables of the library.
"""
import argparse
import pathlib
import re
import sys

from mpmath import mp, mpf, cos, pi, nstr

mp.dps = 50
DIGITS = 20
PER_LINE = 5

ROOT = pathlib.Path(__file__).resolve().parent.parent / "include" / "LNIT" / "AdaptiveQuadratures"

# (header, array name, n) ; the rule has n+1 nodes
TABLES = [
    ("ClenshawCurtisAdaptiveQuadrature.hpp",       "s_wi09", 8),
    ("ClenshawCurtisAdaptiveQuadrature.hpp",       "s_wi17", 16),
    ("ClenshawCurtisAdaptiveQuadrature.hpp",       "s_wi33", 32),
    ("ClenshawCurtisHybridAdaptiveQuadrature.hpp", "s_wi",   12),
    ("GLCCAdaptiveQuadrature.hpp",                 "s_wi_cc", 32),
]


def clenshaw_curtis_weights(n):
    weights = []
    for k in range(n + 1):
        s = mpf(0)
        for j in range(1, n // 2 + 1):
            b = mpf(1) / 2 if 2 * j == n else mpf(1)
            s += b / (4 * j * j - 1) * cos(2 * j * k * pi / n)
        c = 1 if k % n == 0 else 2
        weights.append(mpf(c) / n * (1 - 2 * s))
    return weights


def literal(x):
    return nstr(x, DIGITS, strip_zeros=False)


def cxx_array(name, weights):
    n = len(weights)
    width = max(len(literal(w)) for w in weights)
    lines = [f"\tstatic constexpr std::array<Scalar, {n}> {name} = {{"]
    for start in range(0, n, PER_LINE):
        chunk = weights[start:start + PER_LINE]
        items = ", ".join(f"Scalar({literal(w).ljust(width)})" for w in chunk)
        end = "};" if start + PER_LINE >= n else ", "
        lines.append(f"\t\t{items}{end}")
    return "\n".join(lines)


def update_header(path, name, text):
    source = path.read_text()
    pattern = re.compile(
        r"\tstatic constexpr std::array<Scalar, \d+>\s+" + re.escape(name) + r"\s*=\s*\{.*?\};",
        re.S,
    )
    if len(pattern.findall(source)) != 1:
        sys.exit(f"{path}: expected exactly one definition of {name}")
    path.write_text(pattern.sub(lambda _: text, source))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--update", action="store_true", help="rewrite the headers in place")
    args = parser.parse_args()

    for header, name, n in TABLES:
        weights = clenshaw_curtis_weights(n)
        assert abs(sum(weights) - 2) < mpf(10) ** -45
        text = cxx_array(name, weights)
        if args.update:
            update_header(ROOT / header, name, text)
            print(f"updated {header}: {name}")
        else:
            print(f"// {header}")
            print(text)
            print()


if __name__ == "__main__":
    main()
