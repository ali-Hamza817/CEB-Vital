"""Simulator for the CEB-ViTAL study (continuation of ViTAL, Ali et al., Sci. Rep. 15:15400, 2025).

Derived from the EB-ViTAL v1 simulator. With side=100, bs0=(5,5) and sched=None every protocol
reproduces the earlier runs stored in results/regress_*.json bit-for-bit (checked in run_experiments.py).

Protocols (proto)
  LEACH, EE-LEACH*, I-LEACH*, 4F, ViTAL      as in the reconstruction
  DEEC                                         Qing, Zhu & Wang (2006), homogeneous form
  EB-ViTAL                                     M1 + M2 + M3 (eb flags switch components)
Activity scheduling (sched), applies to any protocol
  None     every alive node senses and reports each round
  "cov"    CovSched: coverage-preserving, energy-weighted greedy active set (PREREGISTRATION_v2.md)
  "rand"   RandSleep: same number of active nodes as CovSched, chosen uniformly at random
CEB-ViTAL = EB-ViTAL + sched="cov" on the HSA deployment.
"""
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree

EELEC, EFS, EMP, EDA = 50e-9, 10e-12, 0.0013e-12, 5e-9
D0 = np.sqrt(EFS / EMP)
SIDE, BS0, E0, P, RTX, RS = 100.0, (5.0, 5.0), 0.5, 0.1, 15.0, 5.0
PT_DBM, PL0, NEXP, SIGMA, NOISE_DBM = 0.0, 40.0, 3.0, 4.0, -100.0
EPOCH = int(round(1 / P))
L_CTRL = 200            # bits of schedule broadcast received by every alive node per (re)schedule
_OFF = np.array([(i, j) for i in range(-6, 7) for j in range(-6, 7)])


def tx(L, d):
    """First-order radio model: energy to send L bits over distance d."""
    return L * EELEC + L * np.where(d < D0, EFS * d**2, EMP * d**4)


def mround(x):
    return np.floor(x + 0.5)


# ---------------------------------------------------------------- coverage (binary disk, 1 m cells)
def _stamps(xy, side):
    s = int(side)
    c = np.floor(xy).astype(int)
    cells = (c[:, None, :] + _OFF[None]).reshape(-1, 2)
    owner = np.repeat(np.arange(len(xy)), len(_OFF))
    ok = (cells >= 0).all(1) & (cells < s).all(1)
    cells, owner = cells[ok], owner[ok]
    hit = np.hypot(cells[:, 0] + 0.5 - xy[owner, 0], cells[:, 1] + 0.5 - xy[owner, 1]) <= RS
    return owner[hit], cells[hit, 1] * s + cells[hit, 0]


def coverage(xy, side=SIDE):
    """Covered fraction of the field (ViTAL Eq. 2, 8)."""
    _, cell = _stamps(xy, side)
    return np.unique(cell).size / int(side) ** 2


def incidence(xy, side):
    """Sparse node x cell matrix: 1 if the cell centre is within Rs of the node."""
    own, cell = _stamps(xy, side)
    return csr_matrix((np.ones(len(own)), (own, cell)), shape=(len(xy), int(side) ** 2))


def covered_frac(A, active):
    return np.count_nonzero(np.asarray(A[np.flatnonzero(active)].sum(0)).ravel()) / A.shape[1] if active.any() else 0.0


def greedy_cover(A, cand, E, tau):
    """CovSched: pick nodes maximising (new cells covered) x (E/Emax) until tau of the
    cells coverable by all candidates are covered."""
    idx = np.flatnonzero(cand)
    reach = np.asarray(A[idx].sum(0)).ravel() > 0
    target = tau * reach.sum()
    unc = reach.astype(float)
    wgt = np.where(cand, E / E[cand].max(), 0.0)
    chosen = np.zeros(A.shape[0], bool)
    got = 0.0
    while got < target:
        g = A @ unc
        g[chosen] = 0
        i = int(np.argmax(g * wgt))
        if g[i] <= 0 or wgt[i] <= 0:
            break
        chosen[i] = True
        got += g[i]
        unc[A.indices[A.indptr[i]:A.indptr[i + 1]]] = 0
    return chosen


# ---------------------------------------------------------------- deployment
def objective(xy, side=SIDE):
    """ViTAL Eq. (9) with s fixed: cr * dmin (Eq. 9a)."""
    dmin = cKDTree(xy).query(xy, k=2)[0][:, 1].min()
    return coverage(xy, side) * dmin / np.hypot(side - RS, side - RS)


def deploy_random(N, rng, side=SIDE):
    return rng.uniform(0, side, (N, 2))


def deploy_hsa(N, rng, side=SIDE, hms=10, hmcr=0.9, par=0.3, bw=5.0, iters=2000):
    """Harmony search (ViTAL Algorithm 1), bounds Eq. (10)-(13) with Ru = 0."""
    lo, hi = RS, side - RS
    hm = rng.uniform(lo, hi, (hms, N, 2))
    fit = np.array([objective(h, side) for h in hm])
    for _ in range(iters):
        pick = rng.integers(0, hms, N)
        new = hm[pick, np.arange(N)]
        fresh = rng.random(N) > hmcr
        new[fresh] = rng.uniform(lo, hi, (fresh.sum(), 2))
        adj = ~fresh & (rng.random(N) < par)
        new[adj] += rng.uniform(-1, 1, (adj.sum(), 2)) * bw
        np.clip(new, lo, hi, out=new)
        f = objective(new, side)
        w = fit.argmin()
        if f > fit[w]:
            hm[w], fit[w] = new, f
    return hm[fit.argmax()]


def deploy(kind, N, seed, side=SIDE):
    rng = np.random.default_rng([seed, N, 1] if side == SIDE else [seed, N, 1, int(side)])
    return deploy_random(N, rng, side) if kind == "random" else deploy_hsa(N, rng, side)


# ---------------------------------------------------------------- CDS (ViTAL phase 2)
def adjacency(xy):
    A = np.zeros((len(xy), len(xy)), bool)
    p = cKDTree(xy).query_pairs(RTX, output_type="ndarray")
    A[p[:, 0], p[:, 1]] = A[p[:, 1], p[:, 0]] = True
    return A


def greedy_cds(xy, alive):
    idx = np.flatnonzero(alive)
    A = adjacency(xy[idx])
    white = np.ones(len(idx), bool)
    gray = np.zeros(len(idx), bool)
    black = np.zeros(len(idx), bool)
    deg = A.sum(1)
    while white.any():
        gain = np.where(gray, (A & white).sum(1), -1)
        v = gain.argmax() if gain.max() > 0 else np.flatnonzero(white)[deg[white].argmax()]
        black[v], white[v], gray[v] = True, False, False
        nb = A[v] & white
        white[nb], gray[nb] = False, True
    cds = np.zeros(len(xy), bool)
    cds[idx[black]] = True
    return cds


# ---------------------------------------------------------------- EB-ViTAL M2, M3
def k_opt(n_alive, d_alive, side=SIDE):
    k = np.sqrt(n_alive / (2 * np.pi)) * D0 * side / max(d_alive.mean() ** 2, 1e-12)
    return int(round(np.clip(k, 1, n_alive)))


def elect_separated(xy, score, alive, k, side=SIDE):
    r2 = side ** 2 / (np.pi * k)
    chosen = []
    for i in np.argsort(-score, kind="stable"):
        if len(chosen) == k or not alive[i]:
            break
        if not chosen or (((xy[chosen] - xy[i]) ** 2).sum(1) >= r2).all():
            chosen.append(i)
    ch = np.zeros(len(xy), bool)
    ch[chosen] = True
    return ch


def route_multihop(E, chi, xy, dbs, L, ebar):
    k = len(chi)
    p = xy[chi]
    C = np.zeros((k + 1, k + 1))
    C[:k, :k] = tx(L, np.hypot(p[:, None, 0] - p[None, :, 0], p[:, None, 1] - p[None, :, 1])) \
        + L * EELEC * (ebar / np.maximum(E[chi], 1e-12))[None, :]
    np.fill_diagonal(C, 0)
    C[:k, k] = tx(L, dbs[chi])
    _, pred = dijkstra(C.T, directed=True, indices=k, return_predecessors=True)
    spent = np.zeros(k)
    for i in range(k):
        cur = i
        while cur != k:
            nxt = pred[cur]
            spent[cur] += C[cur, k] if nxt == k else tx(L, np.hypot(*(p[cur] - p[nxt])))
            if nxt != k:
                spent[nxt] += L * EELEC
            cur = nxt
    return spent


# ---------------------------------------------------------------- simulation
def simulate(proto, N, seed, deploy_kind="random", L=800, w=(0.25, 0.25, 0.25, 0.25),
             dist_inverted=False, mobile=False, rmax=30000, xy=None, eb=None,
             side=SIDE, bs0=BS0, sched=None, tau=0.95, r_sched=10, cov_trace=None, sleep_frac=0.0, snap_round=None, e_trace=None):
    """One run. Returns lifetime, coverage and delivery metrics. cov_trace: optional list,
    receives the per-round covered fraction (plots only). sleep_frac: post-hoc robustness check;
    every sleeping alive node pays sleep_frac * L * Eelec per round (0 in all pre-registered runs).
    snap_round: if set, the output has "snapshot" = network state at that round (figures only).
    e_trace: optional list, receives (round, mean E, std E of alive nodes, n_alive, n_sensing)."""
    eb = {"m1": True, "m2": True, "m3": True} | (eb or {})
    if xy is None:
        xy = deploy(deploy_kind, N, seed, side)
    rng = np.random.default_rng([seed, N, 2])
    rng_sleep = np.random.default_rng([seed, N, 3])
    shadow = rng.normal(0, SIGMA, N)
    nbr = cKDTree(xy).query_pairs(RTX, output_type="ndarray")
    Ainc = incidence(xy, side)
    E = np.full(N, E0)
    G = np.zeros(N, bool)
    bs = np.array(bs0, float)
    wp, speed, pause = bs.copy(), 0.0, 0
    cds, cds_alive = None, -1
    dead_at = np.full(N, -1)
    pk = no_ch = readings = chs = 0
    snap = None
    res = {}
    active = np.ones(N, bool)
    c0 = covered_frac(Ainc, active)
    cov, cov_key = c0, None
    cov_int, t90, t50, ctrl, act_sum = 0.0, None, None, 0.0, 0
    for r in range(rmax):
        alive = E > 0
        newly = ~alive & (dead_at < 0)
        dead_at[newly] = r
        if r in (500, 1000):
            res[r] = np.clip(E, 0, None).mean()
        if not alive.any():
            break
        if sched is None:
            active = alive
        elif r % r_sched == 0 or (active & ~alive).any():          # (re)schedule
            chosen = greedy_cover(Ainc, alive, E, tau)
            if sched == "rand":                                       # same count, random members
                pick = rng_sleep.choice(np.flatnonzero(alive), int(chosen.sum()), replace=False)
                chosen = np.zeros(N, bool)
                chosen[pick] = True
            active = chosen & alive
            E[alive] -= L_CTRL * EELEC
            ctrl += alive.sum() * L_CTRL * EELEC
            alive = E > 0
            active &= alive
        part = active & alive                                       # nodes that sense and report this round
        if sleep_frac:
            E[alive & ~part] -= sleep_frac * L * EELEC
        na = part.sum()
        key = part.tobytes()
        if key != cov_key:
            cov, cov_key = covered_frac(Ainc, part), key
        cov_int += cov
        if cov_trace is not None:
            cov_trace.append(cov)
        if e_trace is not None:
            ea = E[alive]
            if ea.size:
                e_trace.append((r, float(ea.mean()), float(ea.std()), int(ea.size), int(na)))
        act_sum += na
        if t90 is None and cov < 0.9 * c0:
            t90 = r
        if t50 is None and cov < 0.5 * c0:
            t50 = r
        if na == 0:
            continue
        if mobile:
            if pause > 0:
                pause -= 1
            else:
                step = wp - bs
                dist = np.hypot(*step)
                if dist <= speed:
                    bs, pause = wp.copy(), 10
                    wp, speed = rng.uniform(0, side, 2), rng.uniform(1, 10)
                else:
                    bs = bs + step / dist * speed
        dbs = np.hypot(xy[:, 0] - bs[0], xy[:, 1] - bs[1])

        k = r % EPOCH
        if k == 0 and proto != "DEEC":
            G[:] = False
        T = np.where(part & ~G, P / (1 - P * k), 0.0)

        readings += na
        if proto in ("I-LEACH*", "4F", "ViTAL", "EB-ViTAL"):
            a = part.astype(int)
            neigh = np.bincount(nbr[:, 0], a[nbr[:, 1]], N) + np.bincount(nbr[:, 1], a[nbr[:, 0]], N)
        if proto == "EE-LEACH*":
            thr = T * E / E0
        elif proto == "I-LEACH*":
            m_n = max(neigh[part].mean(), 1e-9)
            thr = T * (neigh / m_n) * (dbs[part].mean() / np.maximum(dbs, 1e-9))
            thr[E < E[part].mean()] = 0
        elif proto in ("4F", "ViTAL"):
            snr = np.maximum(PT_DBM - PL0 - 10 * NEXP * np.log10(np.maximum(dbs, 1)) - shadow - NOISE_DBM, 0)
            dterm = dbs[part].mean() / np.maximum(dbs, 1e-9) if dist_inverted else dbs / max(dbs[part].mean(), 1e-9)
            nv = (w[0] * E / E[part].sum() + w[1] * neigh / N + w[2] * dterm
                  + w[3] * snr / max(snr[part].mean(), 1e-9))
            thr = nv * T
            if proto == "ViTAL":
                if na != cds_alive:
                    cds, cds_alive = greedy_cds(xy, part), na
                thr = np.where(cds, thr, 0)
        elif proto == "EB-ViTAL":
            snr = np.maximum(PT_DBM - PL0 - 10 * NEXP * np.log10(np.maximum(dbs, 1)) - shadow - NOISE_DBM, 0)
            if eb["m1"]:
                w1 = 0.25 + 0.5 * (1 - E[part].mean() / E0)
                wo = (1 - w1) / 3
                score = w1 * E / E[part].max() + wo * (neigh / max(neigh[part].max(), 1)
                                                       + 1 - dbs / max(dbs[part].max(), 1e-9)
                                                       + snr / max(snr[part].max(), 1e-9))
            else:
                score = 0.25 * (E / E[part].sum() + neigh / N + dbs / max(dbs[part].mean(), 1e-9)
                                + snr / max(snr[part].mean(), 1e-9))
            score = np.where(part, score, -np.inf)
            thr = np.where(part, score, 0) * T
        elif proto == "DEEC":
            pi = P * np.maximum(E, 0) / E[part].mean()
            ep = mround(1 / np.maximum(pi, 1e-12))
            kk = r % ep
            G[kk == 0] = False
            thr = np.where(part & ~G & (pi > 0), pi / (1 - pi * kk), 0.0)
        else:
            thr = T

        rnd = rng.random(N)
        if proto == "EB-ViTAL" and eb["m2"]:
            ch = elect_separated(xy, score, part, k_opt(na, dbs[part], side), side)
        else:
            ch = part & (rnd <= thr)
        G |= ch
        chi = np.flatnonzero(ch)
        mem = np.flatnonzero(part & ~ch)
        chs += len(chi)
        if len(chi):
            if proto == "EB-ViTAL" and eb["m3"]:
                E[chi] -= route_multihop(E, chi, xy, dbs, L, E[part].mean()) + L * EDA
            else:
                E[chi] -= tx(L, dbs[chi]) + L * EDA
            dch = np.hypot(xy[mem, None, 0] - xy[None, chi, 0], xy[mem, None, 1] - xy[None, chi, 1])
            j = dch.argmin(1)
            dmin = dch[np.arange(len(mem)), j]
            direct = dbs[mem] <= dmin
            E[mem] -= tx(L, np.where(direct, dbs[mem], dmin))
            np.subtract.at(E, chi[j[~direct]], (EELEC + EDA) * L)
            pk += len(chi) + direct.sum()
        else:
            E[mem] -= tx(L, dbs[mem])
            pk += len(mem)
            no_ch += 1
        if r == snap_round:
            head = np.where(direct, -1, chi[j]) if len(chi) else np.full(len(mem), -1)
            snap = dict(xy=xy.copy(), E=E.copy(), alive=alive.copy(), part=part.copy(), ch=chi.copy(),
                        mem=mem.copy(), head=head, bs=bs.copy(), round=r)

    dead_at[(dead_at < 0) & (E <= 0)] = r
    d = np.sort(dead_at[dead_at >= 0])

    def nth(k):
        return int(d[k - 1]) if len(d) >= k else None

    return dict(FND=nth(1), HND=nth(N // 2), LND=nth(N), packets=int(pk), rounds_without_ch=no_ch,
                readings=int(readings), mean_ch=chs / max(r, 1), deaths=d.tolist(),
                res500=float(res.get(500, 0.0)), res1000=float(res.get(1000, 0.0)),
                C0=float(c0), T90=t90, T50=t50, cov_int=float(cov_int), ctrl_energy=float(ctrl),
                mean_active=act_sum / max(r, 1), rounds=int(r)) | ({"snapshot": snap} if snap is not None else {})


if __name__ == "__main__":
    import time
    rng = np.random.default_rng(0)
    sq = np.array([[x, y] for x in np.arange(5, 100, 10) for y in np.arange(5, 100, 10)], float)
    assert abs(coverage(sq) - 100 * np.pi * 25 / 1e4) < 0.02
    A = incidence(sq, 100)
    assert abs(covered_frac(A, np.ones(100, bool)) - coverage(sq)) < 1e-12
    # greedy cover: duplicated nodes are redundant, the higher-energy copy is chosen
    dup = np.r_[sq, sq]
    Ad = incidence(dup, 100)
    Ed = np.r_[np.full(100, 0.2), np.full(100, 0.5)]
    ch = greedy_cover(Ad, np.ones(200, bool), Ed, 1.0)
    assert ch.sum() == 100 and ch[100:].all()
    ch = greedy_cover(Ad, np.ones(200, bool), Ed, 0.5)
    assert 50 <= ch.sum() <= 51
    t = time.time()
    for p in ("LEACH", "DEEC", "ViTAL", "EB-ViTAL"):
        for s in (None, "cov", "rand"):
            o = simulate(p, 100, 0, sched=s)
            assert o["FND"] <= o["HND"] <= o["LND"] and 0 < o["C0"] <= 1
            print(f"{p:9s} sched={s!s:5s} FND={o['FND']:6d} LND={o['LND']:6d} T90={o['T90']} readings={o['readings']} active={o['mean_active']:.1f}")
    o = simulate("EB-ViTAL", 100, 0, side=200, bs0=(5, 5), sched="cov", rmax=60000)
    print("side=200", {k: o[k] for k in ("FND", "LND", "T90", "C0")})
    print("self-check OK  %.1fs" % (time.time() - t))
