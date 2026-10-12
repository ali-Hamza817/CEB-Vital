"""Independent reconstruction of ViTAL from

  S. J. Ali, J. Iqbal, F. Khan, Y.-K. Moon, "Vigorous technique for augmented
  lifetime in WSNs", Scientific Reports 15:15400 (2025),
  https://doi.org/10.1038/s41598-025-96830-w  (CC BY 4.0)

built only from the paper text. Every unstated choice is listed in
PREREGISTRATION.md (IDs P1..P18 are cited in comments below).

Protocols
  LEACH      Heinzelman et al. 2000 (exact)
  EE-LEACH*  reconstruction (P13): T(n) * Eres/E0
  I-LEACH*   reconstruction (P13): Eres >= mean alive energy;
             T(n) * (neigh/mean neigh) * (mean dist/dist)
  4F         Eq. (20)-(22) CH selection only
  ViTAL      HSA deployment + CDS-restricted 4F selection (P10)
"""
import numpy as np
from scipy.spatial import cKDTree

EELEC, EFS, EMP, EDA = 50e-9, 10e-12, 0.0013e-12, 5e-9   # P3, P4
D0 = np.sqrt(EFS / EMP)
SIDE, BS0, E0, P, RTX, RS = 100.0, (5.0, 5.0), 0.5, 0.1, 15.0, 5.0  # P1, P5
PT_DBM, PL0, NEXP, SIGMA, NOISE_DBM = 0.0, 40.0, 3.0, 4.0, -100.0  # P8
EPOCH = int(round(1 / P))


def tx(L, d):
    """First-order radio model: energy to send L bits over distance d."""
    return L * EELEC + L * np.where(d < D0, EFS * d**2, EMP * d**4)


# ---------------------------------------------------------------- deployment
_GX, _GY = np.meshgrid(np.arange(SIDE) + 0.5, np.arange(SIDE) + 0.5)
_GRID = np.c_[_GX.ravel(), _GY.ravel()]
_OFF = np.array([(i, j) for i in range(-6, 7) for j in range(-6, 7)])


def coverage(xy):
    """Binary disk model (Eq. 2) on 1 m demand cells: covered fraction cr (Eq. 8)."""
    cnt = np.zeros((int(SIDE), int(SIDE)), np.int32)
    c = np.floor(xy).astype(int)
    cells = (c[:, None, :] + _OFF[None]).reshape(-1, 2)
    owner = np.repeat(np.arange(len(xy)), len(_OFF))
    ok = (cells >= 0).all(1) & (cells < SIDE).all(1)
    cells, owner = cells[ok], owner[ok]
    hit = np.hypot(cells[:, 0] + 0.5 - xy[owner, 0], cells[:, 1] + 0.5 - xy[owner, 1]) <= RS
    np.add.at(cnt, (cells[hit, 1], cells[hit, 0]), 1)
    return (cnt > 0).mean()


def objective(xy):
    """Eq. (9) with s fixed (1/s and K are constants): cr * dmin (Eq. 9a)."""
    dmin = cKDTree(xy).query(xy, k=2)[0][:, 1].min()
    return coverage(xy) * dmin / np.hypot(SIDE - RS, SIDE - RS)


def deploy_random(N, rng):
    return rng.uniform(0, SIDE, (N, 2))


def deploy_hsa(N, rng, hms=10, hmcr=0.9, par=0.3, bw=5.0, iters=2000):
    """Harmony search (Algorithm 1, ref. 22), P11. Bounds Eq. (10)-(13) with Ru=0."""
    lo, hi = RS, SIDE - RS
    hm = rng.uniform(lo, hi, (hms, N, 2))
    fit = np.array([objective(h) for h in hm])
    for _ in range(iters):
        pick = rng.integers(0, hms, N)
        new = hm[pick, np.arange(N)]                         # memory consideration
        fresh = rng.random(N) > hmcr
        new[fresh] = rng.uniform(lo, hi, (fresh.sum(), 2))   # random choice
        adj = ~fresh & (rng.random(N) < par)                 # pitch adjustment, Eq. (19)
        new[adj] += rng.uniform(-1, 1, (adj.sum(), 2)) * bw
        np.clip(new, lo, hi, out=new)
        f = objective(new)
        w = fit.argmin()
        if f > fit[w]:
            hm[w], fit[w] = new, f
    return hm[fit.argmax()]


def deploy(kind, N, seed):
    rng = np.random.default_rng([seed, N, 1])
    return deploy_random(N, rng) if kind == "random" else deploy_hsa(N, rng)


# ---------------------------------------------------------------- CDS (P10)
def adjacency(xy):
    A = np.zeros((len(xy), len(xy)), bool)
    p = cKDTree(xy).query_pairs(RTX, output_type="ndarray")
    A[p[:, 0], p[:, 1]] = A[p[:, 1], p[:, 0]] = True
    return A


def greedy_cds(xy, alive):
    """Guha-Khuller greedy connected dominating set over alive nodes (each component)."""
    idx = np.flatnonzero(alive)
    A = adjacency(xy[idx])
    white = np.ones(len(idx), bool)
    gray = np.zeros(len(idx), bool)
    black = np.zeros(len(idx), bool)
    deg = A.sum(1)
    while white.any():
        gain = np.where(gray, (A & white).sum(1), -1)
        if gain.max() > 0:
            v = gain.argmax()
        else:                                   # start a new component
            v = np.flatnonzero(white)[deg[white].argmax()]
        black[v], white[v], gray[v] = True, False, False
        nb = A[v] & white
        white[nb], gray[nb] = False, True
    cds = np.zeros(len(xy), bool)
    cds[idx[black]] = True
    return cds


# ---------------------------------------------------------------- simulation
def simulate(proto, N, seed, deploy_kind="random", L=800, w=(0.25, 0.25, 0.25, 0.25),
             dist_inverted=False, mobile=False, rmax=30000, xy=None):
    """One run. Returns FND/HND/LND, packets to BS, residual energy traces."""
    if xy is None:
        xy = deploy(deploy_kind, N, seed)
    rng = np.random.default_rng([seed, N, 2])
    shadow = rng.normal(0, SIGMA, N)                       # static shadowing, P8
    nbr = cKDTree(xy).query_pairs(RTX, output_type="ndarray")
    E = np.full(N, E0)
    G = np.zeros(N, bool)
    bs = np.array(BS0)
    wp, speed, pause = bs.copy(), 0.0, 0                    # random waypoint state, P14
    cds, cds_alive = None, -1
    dead_at = np.full(N, -1)
    pk = 0
    no_ch = 0
    res = {}
    for r in range(rmax):
        alive = E > 0
        newly = ~alive & (dead_at < 0)
        dead_at[newly] = r
        if r in (500, 1000):
            res[r] = np.clip(E, 0, None).mean()
        na = alive.sum()
        if na == 0:
            break
        if mobile:
            if pause > 0:
                pause -= 1
            else:
                step = wp - bs
                dist = np.hypot(*step)
                if dist <= speed:                       # arrive, pause 10 s, new leg
                    bs, pause = wp.copy(), 10
                    wp, speed = rng.uniform(0, SIDE, 2), rng.uniform(1, 10)
                else:
                    bs = bs + step / dist * speed
        dbs = np.hypot(xy[:, 0] - bs[0], xy[:, 1] - bs[1])

        k = r % EPOCH
        if k == 0:
            G[:] = False
        T = np.where(alive & ~G, P / (1 - P * k), 0.0)

        if proto in ("I-LEACH*", "4F", "ViTAL"):
            a = alive.astype(int)
            neigh = np.bincount(nbr[:, 0], a[nbr[:, 1]], N) + np.bincount(nbr[:, 1], a[nbr[:, 0]], N)
        if proto == "EE-LEACH*":
            thr = T * E / E0
        elif proto == "I-LEACH*":
            m_n = max(neigh[alive].mean(), 1e-9)
            thr = T * (neigh / m_n) * (dbs[alive].mean() / np.maximum(dbs, 1e-9))
            thr[E < E[alive].mean()] = 0
        elif proto in ("4F", "ViTAL"):
            snr = np.maximum(PT_DBM - PL0 - 10 * NEXP * np.log10(np.maximum(dbs, 1)) - shadow - NOISE_DBM, 0)
            dterm = dbs[alive].mean() / np.maximum(dbs, 1e-9) if dist_inverted else dbs / max(dbs[alive].mean(), 1e-9)
            nv = (w[0] * E / E[alive].sum() + w[1] * neigh / N + w[2] * dterm
                  + w[3] * snr / max(snr[alive].mean(), 1e-9))        # Eq. (20)
            thr = nv * T                                               # Eq. (22), P9
            if proto == "ViTAL":
                if na != cds_alive:
                    cds, cds_alive = greedy_cds(xy, alive), na
                thr = np.where(cds, thr, 0)
        else:
            thr = T

        ch = alive & (rng.random(N) <= thr)
        G |= ch
        chi = np.flatnonzero(ch)
        mem = np.flatnonzero(alive & ~ch)
        if len(chi):
            E[chi] -= tx(L, dbs[chi]) + L * EDA
            dch = np.hypot(xy[mem, None, 0] - xy[None, chi, 0], xy[mem, None, 1] - xy[None, chi, 1])
            j = dch.argmin(1)
            dmin = dch[np.arange(len(mem)), j]
            direct = dbs[mem] <= dmin
            E[mem] -= tx(L, np.where(direct, dbs[mem], dmin))
            np.subtract.at(E, chi[j[~direct]], (EELEC + EDA) * L)
            pk += len(chi) + direct.sum()                              # P15
        else:
            E[mem] -= tx(L, dbs[mem])
            pk += len(mem)
            no_ch += 1

    dead_at[(dead_at < 0) & (E <= 0)] = r
    d = np.sort(dead_at[dead_at >= 0])

    def nth(k):
        return int(d[k - 1]) if len(d) >= k else None

    return dict(FND=nth(1), HND=nth(N // 2), LND=nth(N), packets=int(pk), rounds_without_ch=no_ch,
                res500=float(res.get(500, 0.0)), res1000=float(res.get(1000, 0.0)))


if __name__ == "__main__":
    import time
    assert np.isclose(D0, 87.7058, atol=1e-3)
    sq = np.array([[x, y] for x in np.arange(5, 100, 10) for y in np.arange(5, 100, 10)], float)
    assert abs(coverage(sq) - 100 * np.pi * 25 / 1e4) < 0.02           # 100 disjoint disks
    from scipy.sparse.csgraph import connected_components
    xy = deploy_random(100, np.random.default_rng(0))
    c = greedy_cds(xy, np.ones(100, bool))
    A = adjacency(xy)
    assert (c | (A & c).any(1)).all()                                   # dominating
    ncomp, lab = connected_components(A, directed=False)
    for k in range(ncomp):                                              # connected per component
        m = c & (lab == k)
        assert connected_components(A[np.ix_(m, m)], directed=False)[0] == 1
    t = time.time()
    for p in ("LEACH", "EE-LEACH*", "I-LEACH*", "4F", "ViTAL"):
        o = simulate(p, 100, 0)
        assert o["FND"] <= o["HND"] <= o["LND"]
        print(p, o)
    print("self-check OK  %.1fs" % (time.time() - t))
