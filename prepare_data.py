#!/usr/bin/env python3
import json, os
import numpy as np
from astropy.io import fits

# 1. Load the table-based full-sky map from EXT=1
hdul = fits.open('data/convergence_map.fits')
tbl  = hdul[1].data              # BinTableHDU with one column of vectors
colname = tbl.dtype.names[0]     # grab the (only) column name
# Stack into a 2D numpy array: shape = (nrows, vector_length)
data = np.vstack(tbl[colname])
hdul.close()

# 2. Normalize to [0,1]
dmin, dmax = np.nanmin(data), np.nanmax(data)
norm = (data - dmin) / (dmax - dmin)
norm = np.nan_to_num(norm)

# 3. Sample a 100×100 grid for the visual
N = 100
ys = np.linspace(0, norm.shape[0] - 1, N).astype(int)
xs = np.linspace(0, norm.shape[1] - 1, N).astype(int)

particles = [
    {"x": float(ix) / (norm.shape[1] - 1),
     "y": float(iy) / (norm.shape[0] - 1),
     "density": float(norm[iy, ix])}
    for iy in ys for ix in xs
]

# 4. Extract the middle row for the sound data
mid = norm.shape[0] // 2
sounddata = [float(v) for v in norm[mid, :]]

# 5. Write out JSON files
os.makedirs('data', exist_ok=True)
with open('data/particles.json', 'w') as f:
    json.dump(particles, f, indent=2)
with open('data/sounddata.json', 'w') as f:
    json.dump(sounddata, f, indent=2)

print(f"✔ particles.json ({len(particles)} points) and sounddata.json ({len(sounddata)} samples) created.")
