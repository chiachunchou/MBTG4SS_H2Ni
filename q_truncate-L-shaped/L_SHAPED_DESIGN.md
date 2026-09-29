# L-shaped MBTG package design

This package is branched from `q_truncate-target`.  The present directory is
the implementation workspace for combining a fixed L-shaped computational
domain with the moving-boundary truncated-grid (MBTG) method.  The first
version implements the fixed orientation shown below and keeps it disabled by
default for backward-compatible rectangular runs.

Input parameters:

```ini
isLShape=false
l_shape_x1_cut=...
l_shape_x2_cut=...
```

The two cut values use the same physical units as `xi1/xf1` and `xi2/xf2`.
When `isLShape=true`, the lower strip is `x2 < l_shape_x2_cut` over the full
`x1` range, and the upper block is `x2 >= l_shape_x2_cut` with
`x1 < l_shape_x1_cut`.

## What the attached papers contribute

The 1993 H2+OH paper supplies the fixed L-shaped idea: a lower horizontal
asymptotic strip is joined to an upper interaction block on the left.  Its
configuration-space picture is a union of two rectangles, with an absorbing
zone at the outgoing edge.  It is a domain-reduction concept, not a direct
set of parameters for the present H2/Ni calculation.

The MBTG4SS paper supplies the current two-dimensional gas-surface benchmark:
uniform rectangular grid, fourth-order Runge-Kutta propagation, low-density
and small-gradient truncation, boundary extrapolation, barrier probabilities,
local saddle-point fluxes, and wall-clock comparisons.  Its reported
truncation settings `(3,2,2,1)`, `(5,4,4,3)`, and `(7,6,6,5)` are algorithm
parameters and should not be confused with L-shaped geometry parameters.

The two papers therefore imply a two-level domain:

1. `DomainMask`: a fixed geometric L-shaped mask that defines where the
   Hamiltonian problem exists;
2. `TAMask`: the moving MBTG active mask, which is a subset of
   `DomainMask` and changes with time.

## Proposed geometry

The first implementation will support the orientation shown in the attached
L-shaped figure:

```
upper block:  x1 in [left_min, left_max], x2 in [lower_cut, upper_max]
lower strip:  x1 in [left_min, right_max], x2 in [lower_min, lower_cut]
```

The mask is the union of these two rectangles.  Geometry should be specified
in physical coordinates and converted to integer grid indices once during
initialization.  This avoids silently changing the shape when `h1`, `h2`, or
the global origin changes.  The configuration should include an explicit
orientation/shape name so a later rotated or multi-branch shape cannot be
mistaken for this one.

## Changes required in `Scatter2d`

### 1. Build and validate the fixed mask

Add a `DomainMask` buffer with one entry per row-major grid point.  Validate
that the two rectangles are inside the rectangular bounding box, overlap at
the intended cut line, and leave enough points for the finite-difference
stencil.  Set `F`, `FF`, and `PF` to zero outside the mask immediately after
reading `IniCndt.dat` and `PES.dat`.

The current `TAMask` is only allocated for truncated runs.  The L-shaped
branch must keep the two concepts separate: a full-grid L-shaped run needs
`DomainMask` but no moving truncation, while an L-shaped MBTG run uses both
masks.

### 2. Make every propagation loop geometry-aware

All four RK4 stages, normalization, truncation tests, boundary detection,
extrapolation targets, density output, and observable reductions must reject
points for which `DomainMask[idx]` is false.  A target point may be promoted
only when it is inside the L-shaped domain.

At the re-entrant corner and along the two internal L boundaries, the first
version should use Dirichlet zero values outside the domain.  A point whose
five-point/second-order stencil reaches outside the mask must use zero
neighbor values and must not cause an extrapolation target to cross the
missing arm.  This makes the boundary condition explicit and avoids reading
uninitialized values across the cut-out.

### 3. Preserve the physics outputs

Keep the existing original-compatible `Flux_*.dat` and `Prob_*.dat` files.
For the MBTG4SS comparison, map them to the paper's dissociation and three
diffusion barriers using the existing `idx_abs_*`, `idx_dsc_*`, and `idx_df*`
settings.  The L-shaped geometry must not change the definitions of these
observables; it only changes which grid points participate in the sums and
propagation.

### 4. Add diagnostics before optimizing

Log and optionally write:

- total rectangular grid points;
- points in the fixed L-shaped domain and their fraction of the rectangle;
- active MBTG points and their fraction of the L-shaped domain;
- RK4 time, L-mask/geometry overhead, truncation overhead, and I/O time;
- norm and the number of points removed by the geometry versus MBTG.

The distinction between geometry reduction and adaptive reduction is needed
to answer how much extra time the L-shaped simplification saves.

## Required benchmark matrix

Run all cases with the same compiler, MPI/OpenMP placement, input files,
time step, final time, output period, and number of repetitions:

| Case | Geometry | MBTG | Purpose |
|---|---|---|---|
| FG-rect | rectangular | off | original reference |
| FG-L | L-shaped | off | error caused by fixed-domain reduction |
| TG-rect | rectangular | on | current target baseline |
| TG-L | L-shaped | on | combined method |

Report both speedups, `FG-rect/FG-L` and `TG-rect/TG-L`, and the total
speedup `FG-rect/TG-L`.  Compare `TG-L` against `FG-L` for truncation error,
and compare `FG-L` against `FG-rect` for the physical error introduced by the
L-shaped domain.  A visually plausible wave packet is not sufficient: check
the norm, relative L2 error, all six probability files, and all four local
flux files.

## Implementation order

1. Add input parameters and a disabled-by-default geometry parser.
2. Implement `DomainMask` and a full-grid L-shaped run (`FG-L`).
3. Verify norm and observables against the rectangular full-grid run.
4. Restrict MBTG initialization, extrapolation, and promotion to
   `DomainMask` (`TG-L`).
5. Add timing/error diagnostics and run the four-case benchmark matrix.
6. Only after the diagnostics are stable, optimize mask traversal and data
   locality.

The numerical values from the H2+OH figure (`R1`, `R2`, `R3`, `R4`, `r1`,
`r2`, `r3`) should remain documented as a geometry example, not hard-coded
defaults for the H2/Ni package.  The current H2/Ni PES, initial wave packet,
barrier coordinates, and physical units must determine the actual production
geometry.
