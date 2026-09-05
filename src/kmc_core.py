"""
Rejection-free lattice kinetic Monte Carlo for a self renewing copper oxide surface.

Method
------
Bortz-Kalos-Lebowitz (n-fold way) on a periodic square lattice. Every elementary
process is exactly enumerated, so no move is ever rejected and the time advance

    dt = -ln(u) / R_total

is exact for the master equation defined by the process list.

Sampling is O(1) per event. Single site processes are drawn from per class site
pools held as swap-with-last arrays with an index map. Two site processes are
drawn from per class bond pools built over the 2 N lattice bonds, so a pair is
picked uniformly among genuinely eligible pairs rather than by rejection on a
neighbour.

Physical picture
----------------
Each site carries a substrate phase and an adsorbate.

    phase     CU     bare metallic copper, freshly exposed, catalytically inert
              CUO    copper oxide, the active phase, generated in situ
              CUSO4  copper sulphate, catalytically dead

    ads       EMPTY  vacant active site
              O      atomic oxygen from dissociated O2
              SO2    molecularly adsorbed sulphur dioxide
              SO3    adsorbed sulphur trioxide, the branch point

The branch point is the whole model. An adsorbed SO3 either desorbs, which
produces recoverable product, or it reacts into the lattice as sulphate, which
kills the site. The competition between those two rates, and the rate at which
dead sites are recovered, decides whether the surface acts as a catalyst or a sorbent.
"""

import math
import numpy as np

# ---------------------------------------------------------------- state codes

CU, CUO, CUSO4 = 0, 1, 2
EMPTY, O, SO2, SO3 = 0, 1, 2, 3

# single site process ids
P_OXIDISE_CU   = 0   # CU + 1/2 O2 -> CUO                     in situ catalyst genesis
P_ADS_SO2      = 1   # SO2(g) + *  -> SO2*
P_DES_SO2      = 2   # SO2*        -> SO2(g) + *
P_DES_SO3      = 3   # SO3*        -> SO3(g) + *              PRODUCT
P_SULFATE      = 4   # SO3* + CuO  -> CuSO4                   POISON
P_DESULFATE    = 5   # CuSO4       -> CuO + SO3(g)
P_SPALL_SULF   = 6   # CuSO4       -> CU (fresh)              renewal, stressed oxide
P_SPALL_OX     = 7   # CuO         -> CU (fresh)              renewal, erosion
P_ADS_SO3      = 8   # SO3(g) + *  -> SO3*                    reverse of P_DES_SO3
N_SITE_PROC    = 9

# two site (bond) process ids
B_ADS_O2       = 0   # O2(g) + 2*      -> 2 O*
B_DES_O2       = 1   # 2 O*            -> O2(g) + 2*
B_LH           = 2   # SO2* + O*       -> SO3* + *            rate determining
N_BOND_PROC    = 3


class Pool:
    """O(1) insert, remove and uniform draw over a set of integer ids."""

    __slots__ = ("items", "pos", "n")

    def __init__(self, capacity):
        self.items = np.full(capacity, -1, dtype=np.int64)
        self.pos = np.full(capacity, -1, dtype=np.int64)
        self.n = 0

    def add(self, i):
        if self.pos[i] >= 0:
            return
        self.items[self.n] = i
        self.pos[i] = self.n
        self.n += 1

    def remove(self, i):
        p = self.pos[i]
        if p < 0:
            return
        last = self.n - 1
        j = self.items[last]
        self.items[p] = j
        self.pos[j] = p
        self.items[last] = -1
        self.pos[i] = -1
        self.n = last

    def draw(self, rng):
        return int(self.items[rng.integers(self.n)])


class Lattice:
    """Periodic square lattice with precomputed neighbour and bond topology."""

    def __init__(self, L):
        self.L = L
        self.N = L * L
        n = self.N
        # four neighbours per site
        self.nbr = np.empty((n, 4), dtype=np.int64)
        for r in range(L):
            for c in range(L):
                i = r * L + c
                self.nbr[i, 0] = r * L + (c + 1) % L          # right
                self.nbr[i, 1] = r * L + (c - 1) % L          # left
                self.nbr[i, 2] = ((r + 1) % L) * L + c        # down
                self.nbr[i, 3] = ((r - 1) % L) * L + c        # up
        # bonds: 2 per site, right and down. bond b joins bond_ends[b]
        self.n_bonds = 2 * n
        self.bond_ends = np.empty((self.n_bonds, 2), dtype=np.int64)
        for i in range(n):
            self.bond_ends[2 * i] = (i, self.nbr[i, 0])
            self.bond_ends[2 * i + 1] = (i, self.nbr[i, 2])
        # bonds touching a given site (4 of them)
        self.site_bonds = np.empty((n, 4), dtype=np.int64)
        for i in range(n):
            r, c = divmod(i, L)
            left = r * L + (c - 1) % L
            up = ((r - 1) % L) * L + c
            self.site_bonds[i, 0] = 2 * i          # own right bond
            self.site_bonds[i, 1] = 2 * i + 1      # own down bond
            self.site_bonds[i, 2] = 2 * left       # left neighbour's right bond
            self.site_bonds[i, 3] = 2 * up + 1     # up neighbour's down bond


class KMC:
    def __init__(self, lattice, rates, seed=0, init_phase=CU):
        """
        rates: dict mapping process id -> rate constant in s^-1, using keys
               ("site", pid) and ("bond", pid).
        """
        self.lat = lattice
        self.rng = np.random.default_rng(seed)
        n = lattice.N

        self.phase = np.full(n, init_phase, dtype=np.int8)
        self.ads = np.full(n, EMPTY, dtype=np.int8)

        self.k_site = np.array([rates[("site", p)] for p in range(N_SITE_PROC)], dtype=float)
        self.k_bond = np.array([rates[("bond", p)] for p in range(N_BOND_PROC)], dtype=float)

        self.site_pools = [Pool(n) for _ in range(N_SITE_PROC)]
        self.bond_pools = [Pool(lattice.n_bonds) for _ in range(N_BOND_PROC)]

        self.t = 0.0
        self.n_events = 0
        # tallies
        self.n_so3_des = 0          # SO3 desorption events (gross product release)
        self.n_so3_readsorbed = 0   # SO3 readsorption events from the gas
        self.n_so3_from_desulf = 0  # SO3 released by sulphate decomposition
        self.n_so2_adsorbed = 0
        self.n_so2_desorbed = 0
        self.n_sulfation = 0
        self.n_desulfation = 0
        self.n_spall_sulf = 0
        self.n_spall_ox = 0
        self.n_cu_oxidised = 0      # metallic copper consumed into oxide

        for i in range(n):
            self._refresh_site(i)
        for b in range(lattice.n_bonds):
            self._refresh_bond(b)

    # ------------------------------------------------------------ classification

    def _refresh_site(self, i):
        ph = self.phase[i]
        ad = self.ads[i]
        pools = self.site_pools

        want = [False] * N_SITE_PROC
        if ph == CU:
            want[P_OXIDISE_CU] = True
        elif ph == CUO:
            if ad == EMPTY:
                want[P_ADS_SO2] = True
                want[P_ADS_SO3] = True
            elif ad == SO2:
                want[P_DES_SO2] = True
            elif ad == SO3:
                want[P_DES_SO3] = True
                want[P_SULFATE] = True
            want[P_SPALL_OX] = True
        elif ph == CUSO4:
            want[P_DESULFATE] = True
            want[P_SPALL_SULF] = True

        for p in range(N_SITE_PROC):
            if want[p]:
                pools[p].add(i)
            else:
                pools[p].remove(i)

    def _refresh_bond(self, b):
        i, j = self.lat.bond_ends[b]
        pi, pj = self.phase[i], self.phase[j]
        ai, aj = self.ads[i], self.ads[j]
        pools = self.bond_pools

        o2_ads = (pi == CUO and ai == EMPTY and pj == CUO and aj == EMPTY)
        o2_des = (ai == O and aj == O)
        lh = ((ai == SO2 and aj == O) or (ai == O and aj == SO2))

        for p, ok in ((B_ADS_O2, o2_ads), (B_DES_O2, o2_des), (B_LH, lh)):
            if ok:
                pools[p].add(b)
            else:
                pools[p].remove(b)

    def _touch(self, i):
        """Reclassify a site and every bond it participates in."""
        self._refresh_site(i)
        for b in self.lat.site_bonds[i]:
            self._refresh_bond(int(b))

    # ------------------------------------------------------------------- stepping

    def step(self):
        ks, kb = self.k_site, self.k_bond
        sp, bp = self.site_pools, self.bond_pools

        r_site = [ks[p] * sp[p].n for p in range(N_SITE_PROC)]
        r_bond = [kb[p] * bp[p].n for p in range(N_BOND_PROC)]
        R = sum(r_site) + sum(r_bond)
        if R <= 0.0:
            return False

        u = self.rng.random()
        self.t += -math.log(1.0 - u) / R

        x = self.rng.random() * R
        acc = 0.0
        chosen = None
        for p in range(N_SITE_PROC):
            acc += r_site[p]
            if x < acc:
                chosen = ("s", p)
                break
        if chosen is None:
            for p in range(N_BOND_PROC):
                acc += r_bond[p]
                if x < acc:
                    chosen = ("b", p)
                    break
        if chosen is None:
            chosen = ("b", N_BOND_PROC - 1)

        kind, p = chosen
        if kind == "s":
            self._do_site(p, sp[p].draw(self.rng))
        else:
            self._do_bond(p, bp[p].draw(self.rng))

        self.n_events += 1
        return True

    def _do_site(self, p, i):
        if p == P_OXIDISE_CU:
            self.phase[i] = CUO
            self.ads[i] = EMPTY
            self.n_cu_oxidised += 1
        elif p == P_ADS_SO2:
            self.ads[i] = SO2
            self.n_so2_adsorbed += 1
        elif p == P_DES_SO2:
            self.ads[i] = EMPTY
            self.n_so2_desorbed += 1
        elif p == P_DES_SO3:
            self.ads[i] = EMPTY
            self.n_so3_des += 1
        elif p == P_ADS_SO3:
            self.ads[i] = SO3
            self.n_so3_readsorbed += 1
        elif p == P_SULFATE:
            self.phase[i] = CUSO4
            self.ads[i] = EMPTY
            self.n_sulfation += 1
        elif p == P_DESULFATE:
            self.phase[i] = CUO
            self.ads[i] = EMPTY
            self.n_desulfation += 1
            self.n_so3_from_desulf += 1
        elif p == P_SPALL_SULF:
            self.phase[i] = CU
            self.ads[i] = EMPTY
            self.n_spall_sulf += 1
        elif p == P_SPALL_OX:
            self.phase[i] = CU
            self.ads[i] = EMPTY
            self.n_spall_ox += 1
        self._touch(i)

    def _do_bond(self, p, b):
        i, j = self.lat.bond_ends[b]
        i, j = int(i), int(j)
        if p == B_ADS_O2:
            self.ads[i] = O
            self.ads[j] = O
        elif p == B_DES_O2:
            self.ads[i] = EMPTY
            self.ads[j] = EMPTY
        elif p == B_LH:
            if self.ads[i] == SO2:
                self.ads[i] = SO3
                self.ads[j] = EMPTY
            else:
                self.ads[j] = SO3
                self.ads[i] = EMPTY
        self._touch(i)
        self._touch(j)

    # -------------------------------------------------------------------- readout

    @property
    def n_so3_net(self):
        """Net SO3 delivered to the gas: desorption minus readsorption, plus the
        SO3 released when sulphate decomposes."""
        return self.n_so3_des - self.n_so3_readsorbed + self.n_so3_from_desulf

    def coverages(self):
        n = self.lat.N
        return {
            "f_CU":    float(np.count_nonzero(self.phase == CU)) / n,
            "f_CUO":   float(np.count_nonzero(self.phase == CUO)) / n,
            "f_CUSO4": float(np.count_nonzero(self.phase == CUSO4)) / n,
            "th_O":    float(np.count_nonzero(self.ads == O)) / n,
            "th_SO2":  float(np.count_nonzero(self.ads == SO2)) / n,
            "th_SO3":  float(np.count_nonzero(self.ads == SO3)) / n,
            "th_free": float(np.count_nonzero((self.ads == EMPTY) & (self.phase == CUO))) / n,
        }
