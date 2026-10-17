#!/usr/bin/env python3
"""Verify every final obstruction directly on the two binary rank supports.
No projective-plane lookup, graph normal form, or quotient-syndrome routine.
Only Python standard-library exact integer arithmetic is used.
"""
import argparse
import collections
import json
from pathlib import Path


def multiply(a: int, b: int, modulus: int) -> int:
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a & 64:
            a ^= modulus
    return r


def rank(columns: list[int]) -> int:
    rows = {}
    for x in columns:
        while x:
            p = x.bit_length() - 1
            if p in rows:
                x ^= rows[p]
            else:
                rows[p] = x
                break
    return len(rows)


def verify(path: Path, modulus: int, independent: bool = False) -> dict:
    mt = [[multiply(a, b, modulus) for b in range(64)] for a in range(64)]
    def evaluate(h, v):
        return (mt[h & 63][v & 63]
                ^ mt[(h >> 6) & 63][(v >> 6) & 63]
                ^ mt[h >> 12][v >> 12])
    inv = [0] + [next(b for b in range(1, 64) if mt[a][b] == 1) for a in range(1, 64)]
    hist = collections.Counter()
    example = None
    n = 0
    for lineno, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip():
            continue
        a = list(map(int, line.split()))
        expected = 15 if independent else 17
        assert len(a) == expected, (path, lineno, len(a))
        if independent:
            orbit, k, param, result, heavy, h0, h1 = a[:7]
            basis = a[7:]
        else:
            orbit, k, core_heavy, x, y, result, heavy, h0, h1 = a[:9]
            basis = a[9:]
        assert result == 3 and rank(basis) == 8, (path, lineno)
        c0 = [evaluate(h0, v) for v in basis]
        c1 = [evaluate(h1, v) for v in basis]
        r0, r1 = rank(c0), rank(c1)
        joint = rank([u | (v << 6) for u, v in zip(c0, c1)])
        assert (r0, r1, joint) == (4, 4, 8), (path, lineno, r0, r1, joint)
        # Compute point weights without a projective chart.
        space = [0]
        for v in basis:
            space += [w ^ v for w in space[:]]
        assert len(set(space)) == 256
        assert basis[:2] == [1, 64] and any(v >> 12 for v in basis)
        counts = collections.Counter()
        for v in space[1:]:
            xyz = (v & 63, (v >> 6) & 63, v >> 12)
            scale = inv[next(x for x in xyz if x)]
            key = tuple(mt[scale][x] for x in xyz)
            counts[key] += 1
        assert max(counts.values()) <= 3
        assert sum(n == 3 for n in counts.values()) == heavy, (path, lineno, heavy)
        hist[heavy] += 1
        if example is None or (heavy == 0 and example['heavy_points'] != 0):
            example = dict(orbit=orbit, projection_type=k, heavy_points=heavy,
                           basis=basis, messages=[h0, h1], codewords=[c0, c1],
                           rank_support_dimensions=[r0, r1], intersection_dimension=0)
        n += 1
    return dict(file=path.name, field_modulus=modulus, verified=n,
                obstruction='two rank-four codewords with disjoint binary rank supports',
                heavy_point_histogram=dict(sorted(hist.items())), example=example)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('file', type=Path)
    p.add_argument('--modulus', type=int, choices=(91, 67), default=91)
    p.add_argument('--independent', action='store_true')
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    result = verify(args.file, args.modulus, args.independent)
    text = json.dumps(result, ensure_ascii=False, indent=2)
    print(text)
    if args.output:
        args.output.write_text(text + '\n')

if __name__ == '__main__':
    main()
