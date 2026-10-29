#######################################
#   DENSE PACKING
# Disks touching each other in a domain
# given by a signed distance field (SDF:
# positive inside, the distance to the
# boundary). Used by
# initial_packing_position.py.
#######################################
"""
A grid of cells of size h (a fraction of the smallest radius) covers the
domain and keeps, for each cell, its clearance: the distance from its
center to the nearest obstacle (the boundary, or a disk already placed),
and the index of that disk (-1 for the boundary). A disk of radius r
fits wherever the clearance is at least r, and touches an obstacle where
the clearance is r. The cells are stored by blocks of B x B cells, only
for the blocks that reach into the domain, and the largest clearance of
each block, and of each superblock of S x S blocks, is kept to find
quickly where a disk fits. Placing a disk lowers the clearance of the
cells that it is now the nearest obstacle of.

Two ways of choosing the radii:

pack_sizes  given radii (e.g. from a distribution), placed from the
            largest: each one goes in the tightest place where it fits
            (the block whose largest clearance is the smallest that
            fits), against its nearest obstacle, then rolls on it until
            it touches a second one. The radii that do not fit anywhere
            are left out.
pack_fill   the radii fill the gaps: the largest gap first, with the
            largest disk it holds (touching the obstacles around it).
            While there are gaps larger than r_max (open areas), the
            radii are drawn uniformly between r_max / 2 and r_max and
            placed as in pack_sizes (always r_max would tile the open
            areas with a lattice of equal disks, a pile-up at r_max in
            the floe size distribution). The smallest radius r_min, with
            r_max = ratio * r_min, is set so that the gaps left are all
            smaller than r_min when the n disks are placed: the densest
            packing of n disks for this ratio.

Needs numpy, scipy and numba.
"""

import time

import numba as nb
import numpy as np
from scipy.ndimage import distance_transform_edt, map_coordinates

B = 16      # cells per side of a block
S = 16      # blocks per side of a superblock
K = 4       # cells per smallest radius


# ----------------------------------------------------------------------
#   domain
# ----------------------------------------------------------------------

def read_grid(path):
    """grid file of the masks (header: nx ny dx x_min y_max projection,
    then the rows from the top), returns the values with row 0 at the
    bottom, dx and the lower left corner"""
    with open(path) as f:
        head = f.readline().split()
        values = np.fromstring(f.read(), dtype=np.float64, sep=" ")
    nx, ny, dx = int(head[0]), int(head[1]), float(head[2])
    x0, y0 = float(head[3]), float(head[4]) - ny * dx
    return values.reshape(ny, nx)[::-1], dx, x0, y0


class Domain:
    """SDF of the domain at the centers of the cells of a grid,
    (x0 + (ix + 1/2) dx, y0 + (iy + 1/2) dx), row 0 at the bottom, as the
    masks of GODAR (pkg/mask/mask_io.f90)"""

    def __init__(self, sdf, dx, x0, y0):
        self.dx, self.x0, self.y0 = float(dx), float(x0), float(y0)
        # the SDF can be > 0 at the edge of the grid (water there, the
        # distance transform does not see the edge), and the bilinear
        # SDF keeps the values of the edge beyond it; GODAR removes the
        # particles that leave the grid: the domain stops one cell before
        # the edge
        ny, nx = sdf.shape
        ex = (np.minimum(np.arange(nx), nx - 1 - np.arange(nx)) - 0.5) * dx
        ey = (np.minimum(np.arange(ny), ny - 1 - np.arange(ny)) - 0.5) * dx
        self.sdf = np.ascontiguousarray(
            np.minimum(sdf, np.minimum(ex[None, :], ey[:, None])),
            dtype=np.float64)

    @classmethod
    def rectangle(cls, lx, ly):
        """[0, lx] x [0, ly]"""
        dx = min(lx, ly) / 100
        x = (np.arange(int(np.ceil(lx / dx)) + 2) - 0.5) * dx
        y = (np.arange(int(np.ceil(ly / dx)) + 2) - 0.5) * dx
        sdf = np.minimum(np.minimum(x, lx - x)[None, :],
                         np.minimum(y, ly - y)[:, None])
        return cls(sdf, dx, -dx, -dx)

    def inside(self, grid, dx, x0, y0):
        """the part of the domain where the 0/1 grid (row 0 at the
        bottom) is 1, e.g. the sea ice extent"""
        # signed distance to the edge of the 1 cells, as build_sdf.py
        sdf = np.where(grid == 1,
                       distance_transform_edt(grid == 1) * dx - dx / 2,
                       -(distance_transform_edt(grid == 0) * dx - dx / 2))
        if (sdf.shape, dx, x0, y0) != (self.sdf.shape, self.dx, self.x0,
                                       self.y0):
            # on the grid of the domain
            ny, nx = self.sdf.shape
            fx = (self.x0 + (np.arange(nx) + 0.5) * self.dx - x0) / dx - 0.5
            fy = (self.y0 + (np.arange(ny) + 0.5) * self.dx - y0) / dx - 0.5
            sdf = map_coordinates(sdf, np.meshgrid(fy, fx, indexing="ij"),
                                  order=1, mode="nearest")
        return Domain(np.minimum(self.sdf, sdf), self.dx, self.x0, self.y0)

    def area(self):
        return float((self.sdf > 0).sum()) * self.dx**2


# ----------------------------------------------------------------------
#   SDF
# ----------------------------------------------------------------------

@nb.njit(cache=True)
def bilinear(a, fx, fy):
    """a at the fractional cell index (fx, fy), as sdf_at in GODAR"""
    ny, nx = a.shape
    ix = min(max(int(np.floor(fx)), 0), nx - 2)
    iy = min(max(int(np.floor(fy)), 0), ny - 2)
    tx = min(max(fx - ix, 0.0), 1.0)
    ty = min(max(fy - iy, 0.0), 1.0)
    return ((1 - tx) * (1 - ty) * a[iy, ix] + tx * (1 - ty) * a[iy, ix + 1]
            + (1 - tx) * ty * a[iy + 1, ix] + tx * ty * a[iy + 1, ix + 1])


@nb.njit(cache=True)
def sdf_at(sdf, dx, x0, y0, x, y):
    return bilinear(sdf, (x - x0) / dx - 0.5, (y - y0) / dx - 0.5)


@nb.njit(cache=True)
def sdf_grad(sdf, dx, x0, y0, x, y):
    """gradient of the SDF at (x, y)"""
    fx = (x - x0) / dx - 0.5
    fy = (y - y0) / dx - 0.5
    gx = (bilinear(sdf, fx + 0.25, fy) - bilinear(sdf, fx - 0.25, fy)) / (
        0.5 * dx)
    gy = (bilinear(sdf, fx, fy + 0.25) - bilinear(sdf, fx, fy - 0.25)) / (
        0.5 * dx)
    return gx, gy


# ----------------------------------------------------------------------
#   clearance grid
# ----------------------------------------------------------------------

@nb.njit(cache=True, parallel=True)
def init_slots(sdf, dx, x0, y0, h, gx0, gy0, nby, nbx):
    """slot of each block that reaches into the domain, -1 for the others"""
    near = np.zeros((nby, nbx), np.bool_)
    half = 0.75 * B * h         # more than half the diagonal of a block
    for by in nb.prange(nby):
        for bx in range(nbx):
            near[by, bx] = sdf_at(sdf, dx, x0, y0, gx0 + (bx + 0.5) * B * h,
                                  gy0 + (by + 0.5) * B * h) > -half - dx
    slot = np.full((nby, nbx), -1, np.int32)
    k = 0
    for by in range(nby):
        for bx in range(nbx):
            if near[by, bx]:
                slot[by, bx] = k
                k += 1
    return slot, k


@nb.njit(cache=True, parallel=True)
def init_clearance(sdf, dx, x0, y0, h, gx0, gy0, slot, nslot):
    nby, nbx = slot.shape
    D = np.empty((nslot, B, B), np.float32)
    own = np.empty((nslot, B, B), np.int32)
    bmax = np.full((nby, nbx), -np.inf, np.float32)
    for by in nb.prange(nby):
        for bx in range(nbx):
            s = slot[by, bx]
            if s < 0:
                continue
            m = -np.inf
            for j in range(B):
                y = gy0 + (by * B + j + 0.5) * h
                for i in range(B):
                    v = sdf_at(sdf, dx, x0, y0,
                               gx0 + (bx * B + i + 0.5) * h, y)
                    D[s, j, i] = v
                    own[s, j, i] = -1
                    m = max(m, v)
            bmax[by, bx] = m
    return D, own, bmax


@nb.njit(cache=True)
def block_max(D, s):
    m = -np.inf
    for j in range(B):
        for i in range(B):
            m = max(m, D[s, j, i])
    return m


@nb.njit(cache=True)
def super_max(bmax, sy, sx):
    nby, nbx = bmax.shape
    m = -np.inf
    for by in range(sy * S, min((sy + 1) * S, nby)):
        for bx in range(sx * S, min((sx + 1) * S, nbx)):
            m = max(m, bmax[by, bx])
    return m


@nb.njit(cache=True)
def init_super(bmax):
    nby, nbx = bmax.shape
    smax = np.empty(((nby + S - 1) // S, (nbx + S - 1) // S), np.float32)
    for sy in range(smax.shape[0]):
        for sx in range(smax.shape[1]):
            smax[sy, sx] = super_max(bmax, sy, sx)
    return smax


class Grid:
    """clearance grid of a domain, cells of size h"""

    def __init__(self, domain, h):
        rows = np.where((domain.sdf > 0).any(axis=1))[0]
        cols = np.where((domain.sdf > 0).any(axis=0))[0]
        if len(rows) == 0:
            raise ValueError("empty domain")
        d = domain.dx
        self.h = h
        self.gx0 = domain.x0 + (cols[0] - 1) * d
        self.gy0 = domain.y0 + (rows[0] - 1) * d
        nbx = int(np.ceil((cols[-1] - cols[0] + 3) * d / (B * h)))
        nby = int(np.ceil((rows[-1] - rows[0] + 3) * d / (B * h)))
        self.slot, nslot = init_slots(domain.sdf, d, domain.x0, domain.y0, h,
                                      self.gx0, self.gy0, nby, nbx)
        self.D, self.own, self.bmax = init_clearance(
            domain.sdf, d, domain.x0, domain.y0, h, self.gx0, self.gy0,
            self.slot, nslot)
        self.smax = init_super(self.bmax)
        self.area = float((self.D > 0).sum()) * h * h

    def args(self):
        return (self.D, self.own, self.slot, self.bmax, self.smax, self.h,
                self.gx0, self.gy0)


# ----------------------------------------------------------------------
#   where a disk fits
# ----------------------------------------------------------------------

@nb.njit(cache=True)
def tightest(bmax, smax, r):
    """block whose largest clearance is the smallest >= r, (-1, -1) if
    none"""
    nby, nbx = bmax.shape
    best, bby, bbx = np.inf, -1, -1
    for sy in range(smax.shape[0]):
        for sx in range(smax.shape[1]):
            if smax[sy, sx] < r:
                continue
            for by in range(sy * S, min((sy + 1) * S, nby)):
                for bx in range(sx * S, min((sx + 1) * S, nbx)):
                    v = bmax[by, bx]
                    if v >= r and v < best:
                        best, bby, bbx = v, by, bx
    return bby, bbx


@nb.njit(cache=True)
def largest(bmax, smax):
    """block with the largest clearance, and that clearance"""
    nby, nbx = bmax.shape
    m, msy, msx = -np.inf, 0, 0
    for sy in range(smax.shape[0]):
        for sx in range(smax.shape[1]):
            if smax[sy, sx] > m:
                m, msy, msx = smax[sy, sx], sy, sx
    bby, bbx = msy * S, msx * S
    for by in range(msy * S, min((msy + 1) * S, nby)):
        for bx in range(msx * S, min((msx + 1) * S, nbx)):
            if bmax[by, bx] == m:
                bby, bbx = by, bx
    return bby, bbx, m


@nb.njit(cache=True)
def nearby(own, slot, h, gx0, gy0, px, py, w, seen, stamp, cand):
    """disks nearest to the cells within w of (px, py), in cand; returns
    their number"""
    nby, nbx = slot.shape
    ix0 = max(0, int((px - w - gx0) / h))
    ix1 = min(nbx * B - 1, int((px + w - gx0) / h))
    iy0 = max(0, int((py - w - gy0) / h))
    iy1 = min(nby * B - 1, int((py + w - gy0) / h))
    nc = 0
    for iy in range(iy0, iy1 + 1):
        for ix in range(ix0, ix1 + 1):
            s = slot[iy // B, ix // B]
            if s < 0:
                continue
            c = own[s, iy % B, ix % B]
            if c >= 0 and seen[c] != stamp:
                seen[c] = stamp
                cand[nc] = c
                nc += 1
    return nc


@nb.njit(cache=True)
def clearance(px, py, cand, nc, X, Y, R, sdf, dx, x0, y0):
    """distance from (px, py) to the nearest obstacle among the boundary
    and the disks of cand"""
    f = sdf_at(sdf, dx, x0, y0, px, py)
    for k in range(nc):
        c = cand[k]
        f = min(f, np.sqrt((px - X[c]) ** 2 + (py - Y[c]) ** 2) - R[c])
    return f


@nb.njit(cache=True)
def fits(qx, qy, r, cand, nc, X, Y, R, sdf, dx, x0, y0):
    """a disk of radius r at (qx, qy) overlaps neither the boundary nor
    the disks of cand (touching them is fine)"""
    tol = 1e-9 * r
    if sdf_at(sdf, dx, x0, y0, qx, qy) < r - tol:
        return False
    for k in range(nc):
        c = cand[k]
        if (qx - X[c]) ** 2 + (qy - Y[c]) ** 2 < (R[c] + r - tol) ** 2:
            return False
    return True


@nb.njit(cache=True)
def two_contacts(px, py, r, cand, nc, X, Y, R, sdf, dx, x0, y0):
    """position closest to (px, py), within r of it, where the disk
    touches two obstacles (two disks, or a disk and the boundary) without
    overlapping any; (px, py) if there is none"""
    bx, by, bd = px, py, r * r
    # the disks it can touch within r of (px, py)
    close = np.empty(nc, np.int64)
    ncl = 0
    for k in range(nc):
        c = cand[k]
        if np.sqrt((px - X[c]) ** 2 + (py - Y[c]) ** 2) - R[c] <= 2.01 * r:
            close[ncl] = c
            ncl += 1
    boundary = sdf_at(sdf, dx, x0, y0, px, py) <= 2.5 * r
    for k in range(ncl):
        a = close[k]
        d1 = R[a] + r
        # touching a and another disk: intersections of two circles
        for l in range(k + 1, ncl):
            b = close[l]
            d2 = R[b] + r
            ex, ey = X[b] - X[a], Y[b] - Y[a]
            d = np.sqrt(ex * ex + ey * ey)
            if d == 0 or d > d1 + d2 or d < abs(d1 - d2):
                continue
            al = (d1 * d1 - d2 * d2 + d * d) / (2 * d)
            hh = np.sqrt(max(d1 * d1 - al * al, 0.0))
            mx, my = X[a] + al * ex / d, Y[a] + al * ey / d
            for sgn in (-1.0, 1.0):
                qx, qy = mx - sgn * hh * ey / d, my + sgn * hh * ex / d
                dq = (qx - px) ** 2 + (qy - py) ** 2
                if dq < bd and fits(qx, qy, r, cand, nc, X, Y, R,
                                    sdf, dx, x0, y0):
                    bx, by, bd = qx, qy, dq
        # touching a and the boundary: SDF = r on the circle around a
        if not boundary:
            continue
        ns = 64
        prev = 0.0
        for s in range(ns + 1):
            t = 2 * np.pi * s / ns
            f = sdf_at(sdf, dx, x0, y0, X[a] + d1 * np.cos(t),
                       Y[a] + d1 * np.sin(t)) - r
            if s > 0 and (f >= 0) != (prev >= 0):
                lo, hi, flo = t - 2 * np.pi / ns, t, prev
                for _ in range(40):
                    mid = 0.5 * (lo + hi)
                    fm = sdf_at(sdf, dx, x0, y0, X[a] + d1 * np.cos(mid),
                                Y[a] + d1 * np.sin(mid)) - r
                    if (fm >= 0) == (flo >= 0):
                        lo, flo = mid, fm
                    else:
                        hi = mid
                t_in = lo if flo >= 0 else hi       # the side inside
                qx, qy = X[a] + d1 * np.cos(t_in), Y[a] + d1 * np.sin(t_in)
                dq = (qx - px) ** 2 + (qy - py) ** 2
                if dq < bd and fits(qx, qy, r, cand, nc, X, Y, R,
                                    sdf, dx, x0, y0):
                    bx, by, bd = qx, qy, dq
            prev = f
    return bx, by


@nb.njit(cache=True)
def against(D, own, slot, h, gx0, gy0, by, bx, r, X, Y, R, seen, stamp,
            cand, sdf, dx, x0, y0):
    """position of a disk of radius r in block (by, bx), which has room
    for it: at the cell where it just fits, moved against the nearest
    obstacle, then rolled on it until it touches a second one"""
    s = slot[by, bx]
    cmin, cj, ci = np.inf, 0, 0
    for j in range(B):
        for i in range(B):
            v = D[s, j, i]
            if v >= r and v < cmin:
                cmin, cj, ci = v, j, i
    cx = gx0 + (bx * B + ci + 0.5) * h
    cy = gy0 + (by * B + cj + 0.5) * h
    nc = nearby(own, slot, h, gx0, gy0, cx, cy, 2 * r + cmin + 2 * h,
                seen, stamp, cand)

    # the clearance at the cell is cmin >= r: moving by cmin - r toward
    # the nearest obstacle keeps it >= r from all the others
    o = own[s, cj, ci]
    if o >= 0:
        ux, uy = X[o] - cx, Y[o] - cy
    else:
        gx, gy = sdf_grad(sdf, dx, x0, y0, cx, cy)
        ux, uy = -gx, -gy
    norm = np.sqrt(ux * ux + uy * uy)
    px, py = cx, cy
    if norm > 0:
        px += (cmin - r) * ux / norm
        py += (cmin - r) * uy / norm
    # the bilinear SDF is not exactly a distance: step back if needed
    for _ in range(5):
        f = sdf_at(sdf, dx, x0, y0, px, py) - r
        if f >= 0:
            break
        gx, gy = sdf_grad(sdf, dx, x0, y0, px, py)
        g2 = gx * gx + gy * gy
        if g2 == 0:
            break
        px -= 1.01 * f * gx / g2
        py -= 1.01 * f * gy / g2
    if not fits(px, py, r, cand, nc, X, Y, R, sdf, dx, x0, y0):
        px, py = cx, cy
    return two_contacts(px, py, r, cand, nc, X, Y, R, sdf, dx, x0, y0)


@nb.njit(cache=True)
def inscribed(px, py, step, cand, nc, X, Y, R, sdf, dx, x0, y0):
    """largest disk around (px, py): local maximum of the clearance
    (compass search from (px, py)), returns its center and radius"""
    f = clearance(px, py, cand, nc, X, Y, R, sdf, dx, x0, y0)
    for _ in range(200):
        bx, by, bf = px, py, f
        for d in range(8):
            a = d * np.pi / 4
            qx, qy = px + step * np.cos(a), py + step * np.sin(a)
            fq = clearance(qx, qy, cand, nc, X, Y, R, sdf, dx, x0, y0)
            if fq > bf:
                bx, by, bf = qx, qy, fq
        if bf > f:
            px, py, f = bx, by, bf
        else:
            step *= 0.5
            if step < 1e-7 * f:
                break
    return px, py, f


@nb.njit(cache=True)
def add_disk(D, own, slot, bmax, smax, h, gx0, gy0, px, py, r, m, dirty):
    """lower the clearance of the cells nearest to the new disk m"""
    nby, nbx = slot.shape
    gmax = -np.inf
    for sy in range(smax.shape[0]):
        for sx in range(smax.shape[1]):
            gmax = max(gmax, smax[sy, sx])
    # beyond r + gmax, the disk is farther than the obstacles already there
    reach = r + gmax
    by0 = max(0, int(np.floor((py - reach - gy0) / (B * h))))
    by1 = min(nby - 1, int(np.floor((py + reach - gy0) / (B * h))))
    bx0 = max(0, int(np.floor((px - reach - gx0) / (B * h))))
    bx1 = min(nbx - 1, int(np.floor((px + reach - gx0) / (B * h))))
    for by in range(by0, by1 + 1):
        ylo = gy0 + by * B * h
        ey = max(ylo - py, 0.0, py - ylo - B * h)
        for bx in range(bx0, bx1 + 1):
            s = slot[by, bx]
            if s < 0:
                continue
            xlo = gx0 + bx * B * h
            ex = max(xlo - px, 0.0, px - xlo - B * h)
            if bmax[by, bx] <= np.sqrt(ex * ex + ey * ey) - r:
                continue
            changed = False
            for j in range(B):
                qy = ylo + (j + 0.5) * h - py
                for i in range(B):
                    qx = xlo + (i + 0.5) * h - px
                    d = np.sqrt(qx * qx + qy * qy) - r
                    if d < D[s, j, i]:
                        D[s, j, i] = d
                        own[s, j, i] = m
                        changed = True
            if changed:
                bmax[by, bx] = block_max(D, s)
                dirty[by // S, bx // S] = True
    for sy in range(by0 // S, by1 // S + 1):
        for sx in range(bx0 // S, bx1 // S + 1):
            if dirty[sy, sx]:
                smax[sy, sx] = super_max(bmax, sy, sx)
                dirty[sy, sx] = False


# ----------------------------------------------------------------------
#   packings
# ----------------------------------------------------------------------

@nb.njit(cache=True)
def sizes_kernel(D, own, slot, bmax, smax, h, gx0, gy0, radii,
                 sdf, dx, x0, y0):
    n = radii.size
    X, Y, R = np.empty(n), np.empty(n), np.empty(n)
    seen = np.full(n, -1, np.int64)
    cand = np.empty(n, np.int64)
    dirty = np.zeros(smax.shape, np.bool_)
    m = 0
    for k in range(n):
        r = radii[k]
        by, bx = tightest(bmax, smax, r)
        if by < 0:
            continue
        px, py = against(D, own, slot, h, gx0, gy0, by, bx, r, X, Y, R,
                         seen, k, cand, sdf, dx, x0, y0)
        X[m], Y[m], R[m] = px, py, r
        add_disk(D, own, slot, bmax, smax, h, gx0, gy0, px, py, r, m, dirty)
        m += 1
    return X[:m], Y[:m], R[:m]


@nb.njit(cache=True)
def fill_kernel(D, own, slot, bmax, smax, h, gx0, gy0, r_min, r_max, nmax,
                seed, sdf, dx, x0, y0):
    np.random.seed(seed)
    X, Y, R = np.empty(nmax), np.empty(nmax), np.empty(nmax)
    seen = np.full(nmax, -1, np.int64)      # stamps m and nmax + m
    cand = np.empty(nmax, np.int64)
    dirty = np.zeros(smax.shape, np.bool_)
    m = 0
    while m < nmax:
        by, bx, gmax = largest(bmax, smax)
        if gmax < r_min:
            break
        if gmax >= r_max:
            # open area: a random radius, so that the disks do not form
            # a lattice of equal disks
            r = np.random.uniform(0.5 * r_max, r_max)
            by, bx = tightest(bmax, smax, r)
            px, py = against(D, own, slot, h, gx0, gy0, by, bx, r, X, Y, R,
                             seen, m, cand, sdf, dx, x0, y0)
        else:
            # the largest disk in the largest gap
            s = slot[by, bx]
            cmax, cj, ci = -np.inf, 0, 0
            for j in range(B):
                for i in range(B):
                    if D[s, j, i] > cmax:
                        cmax, cj, ci = D[s, j, i], j, i
            cx = gx0 + (bx * B + ci + 0.5) * h
            cy = gy0 + (by * B + cj + 0.5) * h
            nc = nearby(own, slot, h, gx0, gy0, cx, cy, cmax + 3 * h,
                        seen, m, cand)
            px, py, r = inscribed(cx, cy, 0.5 * h, cand, nc, X, Y, R,
                                  sdf, dx, x0, y0)
            # in case the search went far: the disks around the center
            nc = nearby(own, slot, h, gx0, gy0, px, py, r + 2 * h,
                        seen, nmax + m, cand)
            r = min(r, clearance(px, py, cand, nc, X, Y, R,
                                 sdf, dx, x0, y0))
            r = min(r * (1 - 1e-9), r_max)
            if r <= 0.5 * r_min:
                break                   # not expected: stop rather than
                                        # place a wrong disk
        X[m], Y[m], R[m] = px, py, r
        add_disk(D, own, slot, bmax, smax, h, gx0, gy0, px, py, r, m, dirty)
        m += 1
    return X[:m], Y[:m], R[:m]


def pack_sizes(domain, radii):
    """places the radii (from the largest) touching each other; returns
    the positions (n, 2) and the radii of the disks that fit, and the
    area of the domain"""
    radii = np.sort(np.asarray(radii, dtype=np.float64))[::-1]
    grid = Grid(domain, radii.min() / K)
    X, Y, R = sizes_kernel(*grid.args(), radii, domain.sdf, domain.dx,
                           domain.x0, domain.y0)
    return np.c_[X, Y], R, grid.area


def pack_fill(domain, n, ratio, seed=0, verbose=True):
    """n disks filling the gaps, radii between r_min and ratio * r_min;
    returns the positions (n, 2), the radii and the area of the domain.
    seed: of the random radii in the open areas.

    The largest gaps are filled first, and a run stops at n disks: it
    holds the n largest disks that fit, the smallest radius is >= r_min.
    r_min is the largest value for which n disks fit (searched by
    bisection: the number of disks does not change smoothly with
    r_min)"""
    # first guess on the large side
    r_min = np.sqrt(domain.area() / (np.pi * n)) / 4
    lo, hi, best = None, None, None
    for _ in range(40):
        t = time.time()
        grid = Grid(domain, r_min / K)
        X, Y, R = fill_kernel(*grid.args(), r_min, ratio * r_min, n, seed,
                              domain.sdf, domain.dx, domain.x0, domain.y0)
        if verbose:
            print(f"  r_min = {r_min:.1f} m: {len(R)} disks, smallest "
                  f"{R.min():.1f} m, {100 * np.pi * np.sum(R**2) / grid.area:.2f}"
                  f" % of the domain ({time.time() - t:.1f} s)")
        if len(R) == n:
            lo, best = r_min, (np.c_[X, Y], R, grid.area)
            if R.min() <= 1.03 * r_min:
                break                   # the gaps left are about < r_min
        else:
            hi = r_min
        if lo and hi:
            if hi < 1.005 * lo:
                break
            r_min = np.sqrt(lo * hi)
        elif hi:
            r_min *= np.sqrt(len(R) / n)
        else:
            r_min = np.sqrt(r_min * R.min())
    if best is None:
        raise RuntimeError(f"pack_fill: never {n} disks")
    return best
