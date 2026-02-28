# Discoveries

## CLASS build: `-ffast-math` causes HyRec failure (2026-02-25)

Building CLASS (v2.9.4, dmeff branch) with `gcc-14 -O4 -ffast-math` causes a
floating-point precision issue in HyRec.  During the sigma8 shooting,
HyRec's recombination integration overshoots to a tiny negative redshift
(`z ≈ -4.8e-13`), which triggers the strict bounds check in
`background_tau_of_z` (line 201 of `source/background.c`):

```
out of range: z=-4.792174e-13 < z_min=0.000000e+00
```

**Fix:** rebuild without `-ffast-math`:

```bash
cd ../class_dmeff_rui_used
make clean && make class CC=gcc-14 OPTFLAG="-O3"
```

This keeps full optimization (`-O3`) and OpenMP but avoids the
floating-point reordering that causes the rounding error.


## CLASS grid edge: DAO residual beyond `P_k_max` (2026-02-26)

CLASS outputs a grid point slightly beyond `P_k_max_h/Mpc = 200`
(at k ≈ 203 h/Mpc).  At this edge point, DAO oscillations can produce a
non-negligible residual (e.g. ratio = 0.058 vs T_WDM = 0.017) that
causes the envelope condition to fail even though all points within the
intended range are fine.

**Fix:** restrict the envelope check to `k <= 200 h/Mpc`.


## Paper erratum: k_hm values in Table 1 (2026-02-25)

Two k_hm entries in Table 1 of COZMIC I (arXiv:2410.03635) differ from
values computed directly from the CLASS transfer function output files:

1. **n=2, m=1e-2, envelope (σ=7.1e-24):** paper says k_hm = 21.8 Mpc⁻¹,
   but the actual first crossing of T = d_dmeff/d_cdm = 0.5 is at 11.70.
   The DAO oscillation causes the signed ratio to go negative and |T|
   rises back above 0.5; T² crosses 0.25 upward at k ≈ 21.86.  The paper
   likely reports this second (upward) crossing of T² = 0.25 rather than
   the first (downward) crossing.

2. **n=2, m=1e-2, halfmode (σ=1.3e-25):** paper says k_hm = 58.0 Mpc⁻¹,
   CLASS output gives 56.92 (~1.9% discrepancy).  Minor rounding or
   interpolation difference.


## Paper erratum: spurious h⁻¹ factor in Eq. 4 (2026-02-28)

Equation 4 of COZMIC I (arXiv:2410.03635) includes a factor of h⁻¹ that
should not be there.


## "Midpoint" sigma exceeds envelope for m=1e-4 (2026-02-25)

The estimated IDM upper bounds from Section 6.4 (stored as "midpoint" in
sim-table.dat) fall outside the halfmode–envelope bracket for m=1e-4:

- **n=2, m=1e-4:** midpoint 2.91e-27 > envelope 2.80e-27 (slightly above)
- **n=4, m=1e-4:** midpoint 1.05e-25 > envelope 3.40e-26 (factor ~3 above)

The other four masses (1e-2, 1 GeV for both n=2 and n=4) have midpoint
values correctly between halfmode and envelope.  The m=1e-4 cases have
the smallest DAOs, so the abundance-matched interpolation (which maps
IDM subhalo counts to effective WDM masses) can overshoot the envelope
cross section.
