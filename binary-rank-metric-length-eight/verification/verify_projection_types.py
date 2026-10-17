#!/usr/bin/env python3
"""Check all 651 binary planes in F64 and all eleven projection orbits."""
from collections import Counter
from verify_witnesses import multiply


def check() -> dict:
    m = [[multiply(a, b, 91) for b in range(64)] for a in range(64)]
    inv = [0] + [next(b for b in range(1, 64) if m[a][b] == 1) for a in range(1, 64)]
    def key(plane):
        return min(tuple(sorted(m[inv[a]][b] for b in plane)) for a in plane)
    planes = [(a, b, a ^ b) for a in range(1, 64)
              for b in range(a + 1, 64) if b < (a ^ b)]
    assert len(planes) == 651
    counts = Counter(key(p) for p in planes)
    reps = [14, 2, 4, 6, 8, 12, 16, 20, 22, 28, 36]
    selected = [key((1, b, 1 ^ b)) for b in reps]
    assert len(set(selected)) == 11 and set(selected) == set(counts)
    assert sorted(counts.values()) == [21] + [63] * 10
    tr = []
    for a in range(64):
        x, t = a, 0
        for _ in range(6):
            t ^= x
            x = m[x][x]
        assert t in (0, 1)
        tr.append(t)
    projections = [tuple(x for x in range(64) if not tr[x] and not tr[m[b][x]])
                   for b in reps]
    assert all(len(z) == 16 for z in projections)
    all_z = set()
    orbit_sizes = []
    for z in projections:
        orbit = {tuple(sorted(m[c][x] for x in z)) for c in range(1, 64)}
        assert not all_z.intersection(orbit)
        all_z.update(orbit)
        orbit_sizes.append(len(orbit))
    assert len(all_z) == 651 and orbit_sizes == [21] + [63] * 10
    return dict(binary_planes=651, scalar_classes=11,
                projection_orbit_sizes=orbit_sizes, covered_four_spaces=len(all_z))

if __name__ == '__main__':
    print(check())
