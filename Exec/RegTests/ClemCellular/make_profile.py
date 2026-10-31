"""Write the ZND table the case starts from, from Validation/detonation/znd.npz (znd.py,
Cantera, Davis mechanism, 2H2 + O2 + 7Ar, 6.67 kPa, 298 K, D_CJ = 1616.02 m/s).
Columns, CGS: distance behind the shock [cm], rho [g/cm3], T [K], u_lab [cm/s], Y_k in the
order of mechanism.H. usage: python3 make_profile.py <znd.npz> <output.dat>"""
import sys
import numpy as np

z = np.load(sys.argv[1])
x = z["x"] * 100.0
keep = x <= 30.0
order = ["H2", "H", "AR", "N2", "HE", "O", "OH", "HCO", "HO2", "H2O", "CO", "O2", "H2O2", "CO2"]
assert list(z["names"]) == order, "species order differs from Davis/mechanism.H"
with open(sys.argv[2], "w") as f:
    f.write(f"# ZND, D_CJ = {float(z['D']) * 100.0:.2f} cm/s; x[cm] rho[g/cm3] T[K] u_lab[cm/s] " + " ".join(order) + "\n")
    f.write(f"{int(keep.sum())}\n")
    for i in np.where(keep)[0]:
        row = [x[i], z["rho"][i] * 1e-3, z["T"][i], z["u"][i] * 100.0] + list(z["Y"][:, i])
        f.write(" ".join(f"{v:.12e}" for v in row) + "\n")
