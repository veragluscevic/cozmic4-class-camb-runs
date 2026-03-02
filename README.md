# cozmic4-class-camb-runs

## COZMIC1-template-files/

Reference files from COZMIC I (arXiv:2410.03635).

- **CDM_class_sync_tk.dat** — Raw CLASS synchronous-gauge CDM transfer function output for COZMIC I cosmology (h=0.7, Omega_m=0.286, Omega_b=0.049, sigma_8=0.82, n_s=0.96). Used as the CDM reference by `plot_transfer.py` and `analyze_pk.py`.
- **test_CDM.dat** — Same CDM run converted to CAMB 13-column format via `class_to_camb.py`. Used as a validation reference by `test_conversion.py`.
- **minimal_syncronous.ini** — CLASS input template for synchronous-gauge runs. Used by `prepare-transfers.sh`, `generate_all_inis_from_sim-table.py`, and `analyze_pk.py` to generate per-model ini files.
- **minimal_newtonian.ini** — CLASS input template for Newtonian-gauge runs.
- **data_tk/** — Original CAMB-formatted transfer functions from Rui. Used by `test_conversion.py` as the reference CAMB output for validating `class_to_camb.py` conversions.
