# Scatter2d migration

`Scatter2d` now follows the target package's execution model:

- the 2-D fields use contiguous row-major buffers and `i1 * W1 + i2` indexing;
- truncated RK4 work is restricted to the TA bounding box and guarded by `TAMask`;
- `TB`, `TBL`, `ExFF`, and the extrapolation loop retain the target MBTG flow;
- the custom `IniCndt.dat` and `PES.dat` inputs are preserved;
- the original physical observables are written as `Flux_d.dat`, `Flux_1.dat`,
  `Flux_2.dat`, `Flux_3.dat`, `Prob_r_a.dat`, `Prob_a_d.dat`, `Prob_d_1.dat`,
  `Prob_1_2.dat`, `Prob_2_3.dat`, and `Prob_3_p.dat`.

The port also restores the original cutoff rule: a point is removed only when its
probability and both first differences are below their thresholds. Flat-index
boundary checks are performed before neighbor access, and all array ownership uses
matching `new[]`/`delete[]` pairs.

The old compile-time model-specific potential functions are not carried over;
the package uses the file-defined PES.  The original four saddle-point fluxes
and six conditional-probability regions are retained for the physical analysis.

## Vista regression inputs

The supplied `tests/0.fg` and `tests/1.tg_7665` inputs use the same 801 x 1801
grid and the same `IniCndt.dat`/`PES.dat`.  Their initial truncated state is
`TA=7893, TB=332`; this is a useful initialization regression check.  Vista's
`IniCndt.dat` stores two values per row, while `PES.dat` stores one real value
per row; the reader accepts both forms and accepts subnormal values such as
`1e-308`.

The original `idx_abs_*`, `idx_dsc_*`, and `idx_df*` settings are restored and
drive the `Flux_*.dat` and `Prob_*.dat` observables.  The target-only
`idx_des_*` and `idx_dis_*` settings apply only to the optional theoretical
remaining-probability re-normalization when `isReNorm=true`; for the original
physical analysis, leave `isReNorm=false` (the default).  `tests/2.tg_5443`
does not contain its own `IniCndt.dat`, and its log starts with an invalidly
tiny normalization, so it is retained as historical output rather than a
standalone regression input.
