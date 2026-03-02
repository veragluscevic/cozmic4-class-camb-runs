import numpy as np
from scipy.interpolate import Akima1DInterpolator
import matplotlib.pyplot as plt
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--index', type=int, default=1)
parser.add_argument('--class-file', default='transfers/camb_n2_1e-2GeV_7.1e-24_tk.dat')
parser.add_argument('--camb-file', default='COZMIC1-files/data_tk/idm_n2_1e-2GeV_envelope_z99_Tk.dat')
parser.add_argument('-o', '--output', type=str, default=None,
                    help='Save plot to the specified path (e.g. plots/test_conversion.png).')
args = parser.parse_args()

index = args.index
class_file = args.class_file
camb_file = args.camb_file

data_class = np.loadtxt(class_file)
data_camb = np.loadtxt(camb_file)
data_cdm = np.loadtxt('COZMIC1-files/test_CDM.dat')

k_camb = data_camb[:,0]
tk_camb = data_camb[:, index]

k_class = data_class[:,0]
tk_class = data_class[:, index]

k_cdm = data_cdm[:,0]
tk_cdm = data_cdm[:, index]

# Hybrid interpolation: use log-log for positive data, Akima for oscillatory/negative
# Clip k to valid range to prevent extrapolation artifacts at boundaries
k_min, k_max = k_camb.min(), k_camb.max()
k_clipped = np.clip(k_class, k_min, k_max)

# Check if all values are positive (use log-log interpolation if so)
if np.all(tk_camb > 0):
    # Log-log interpolation: best for smooth positive transfer functions
    tk_camb_interp = np.exp(np.interp(np.log(k_clipped), np.log(k_camb), np.log(tk_camb)))
else:
    # Akima interpolation: handles negative values and oscillatory behavior
    sort_idx = np.argsort(k_camb)
    k_camb_sorted = k_camb[sort_idx]
    tk_camb_sorted = tk_camb[sort_idx]
    akima_interp = Akima1DInterpolator(k_camb_sorted, tk_camb_sorted)
    tk_camb_interp = akima_interp(k_clipped)

# Interpolate CDM reference onto the same k_clipped grid
if np.all(tk_cdm > 0):
    tk_cdm_interp = np.exp(np.interp(np.log(k_clipped), np.log(k_cdm), np.log(tk_cdm)))
else:
    sort_idx = np.argsort(k_cdm)
    tk_cdm_interp = Akima1DInterpolator(k_cdm[sort_idx], tk_cdm[sort_idx])(k_clipped)


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle(f'{class_file}  (index={index})')

ax1.semilogx(k_clipped, (np.abs(tk_camb_interp) - (np.abs(tk_class)))**2 / tk_cdm_interp**2, '.', color='k')
#ax1.set_ylim(0.7,1.3)
ax1.set_xlabel('k')
ax1.set_ylabel('$(|T_{camb}| - |T_{class}|)^2 / T_{cdm}^2$')
ax1.legend()

ax2.loglog(k_camb, (tk_camb), '.', color='b', label='Rui')
ax2.loglog(k_class, np.abs(tk_class), '.', color='r', label='CLASS')
ax2.loglog(k_cdm, (tk_cdm), ':', color='gray', linewidth=0.8, label='CDM')
ax2.set_xlim(10,200)
ax2.set_xlabel('k')
ax2.set_ylabel('$T_k^2$')
ax2.legend()

plt.tight_layout()

max_rel_diff = np.max((np.abs(tk_camb_interp) - np.abs(tk_class))**2 / tk_cdm_interp**2)
if max_rel_diff < 0.01:
    print(f"YAY! TEST for i={index} PASSED! (max relative difference: {max_rel_diff:.4f})")
else:
    print(f"HMM. TEST for i={index} DOESN'T LOOK GREAT. CHECK PLOTS. (max relative difference: {max_rel_diff:.4f})")

if args.output:
    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    plt.savefig(args.output, dpi=300)
    print(f'Saved to {args.output}')

plt.show()

