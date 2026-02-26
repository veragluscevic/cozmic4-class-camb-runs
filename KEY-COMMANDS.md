# Key Commands

## 1. Generate CLASS ini files

```bash
python generate_inis.py
```

Reads `sim-table.dat` and templates (`minimal_syncronous.ini`, `minimal_newtonian.ini`) to create ini files in `inis/`.

## 2. Run CLASS

```bash
bash run_class.sh
```

Runs CLASS on all ini files in `inis/`, saves output to `output/`.

## 3. Convert CLASS to CAMB format

```bash
bash convert_all.sh
```

Converts all CLASS output (sync + newt + background) to CAMB 13-column format in `converted/`.

Single file:

```bash
python class_to_camb.py output/n2_1e-2GeV_1.3e-25_sync_tk.dat output/n2_1e-2GeV_1.3e-25_newt_tk.dat output/n2_1e-2GeV_1.3e-25_sync_background.dat -o converted/camb_n2_1e-2GeV_1.3e-25_tk.dat
```

## 4. Plot transfer functions

```bash
python plot_transfer.py                                # defaults: n=2, m=1e-2
python plot_transfer.py -n 4 -m 1e-4                   # specific (n, m)
python plot_transfer.py -m all --no-midpoint            # all masses, skip midpoint
python plot_transfer.py -m all --no-midpoint -s         # save to plots/
python plot_transfer.py --mwdm 6.5                      # change WDM reference mass
```

## 5. Analyze transfer functions (compute k_hm)

```bash
python analyze_pk.py                                    # defaults: n=2, m=1e-2
python analyze_pk.py -m all -n 2                        # all masses for n=2
python analyze_pk.py -m all -n 4 --no-midpoint          # skip midpoint entries
python analyze_pk.py -m all --save-khm                  # write k_hm to sim-table.dat
python analyze_pk.py -m 1e-2 --plot -s                  # also plot and save figure
python analyze_pk.py -m all --save-khm --plot -s        # do everything
```

## 6. Compare CLASS vs CAMB output

```bash
python plot_test.py --class-file <class_file> --camb-file <camb_file> -i <column_index>
python plot_test.py --class-file <class_file> --camb-file <camb_file> -i 1 -s
```
