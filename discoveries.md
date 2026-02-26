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
