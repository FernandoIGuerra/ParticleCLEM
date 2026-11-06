"""Write a table of the format of make_profile.py from a plotfile of the 1-D run
(Detonation1D), so the case starts from that solution instead of ZND. The front is the
rightmost point above twice the ambient pressure, as in Detonation1D/det1d_long_velocity.ipynb.
Columns, CGS: distance behind the front [cm], rho [g/cm3], T [K], u_lab [cm/s], Y_k in the
order of mechanism.H. usage: python3 make_profile_1d.py <plotfile> <output.dat>"""
import sys
import numpy as np
import yt

yt.set_log_level(50)
P_FRONT = 2.0 * 6.67e4   # dyn/cm2
Y_LINE = 0.0705          # cm; off the cell faces of all levels
ds = yt.load(sys.argv[1])
ray = ds.ortho_ray(0, (Y_LINE, 0.0))
i = np.argsort(ray["index", "x"].d)
x = ray["index", "x"].d[i]
p = ray["boxlib", "pressure"].d[i]
above = np.nonzero(p >= P_FRONT)[0]
assert above[-1] < len(x) - 1, "the shock has left the domain in this plotfile"
k = above[-1]
xf = x[k] + (P_FRONT - p[k]) / (p[k + 1] - p[k]) * (x[k + 1] - x[k])

order = ["H2", "H", "AR", "N2", "HE", "O", "OH", "HCO", "HO2", "H2O", "CO", "O2", "H2O2", "CO2"]
assert sorted(order) == sorted(f[1][2:-1] for f in ds.field_list if f[1].startswith("Y(")), \
    "species differ from Davis/mechanism.H"
get = lambda name: ray["boxlib", name].d[i]
cols = [get("density"), get("Temp"), get("x_velocity")] + [get(f"Y({s})") for s in order]
keep = np.nonzero((x < xf) & (x >= xf - 30.0))[0][::-1]     # ascending in the distance behind
with open(sys.argv[2], "w") as f:
    f.write(f"# 1-D run {sys.argv[1]}, t = {float(ds.current_time) * 1e6:.2f} us, front x = {xf:.4f} cm; "
            "x[cm] rho[g/cm3] T[K] u_lab[cm/s] " + " ".join(order) + "\n")
    f.write(f"{len(keep)}\n")
    for j in keep:
        f.write(" ".join(f"{v:.12e}" for v in [xf - x[j]] + [c[j] for c in cols]) + "\n")
print(f"t = {float(ds.current_time) * 1e6:.2f} us, front at x = {xf:.3f} cm, {len(keep)} points, "
      f"max distance {xf - x[keep[-1]]:.2f} cm")
