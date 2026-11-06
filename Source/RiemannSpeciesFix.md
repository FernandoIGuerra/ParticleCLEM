# Species flux of the Godunov Riemann solver: upwinding by the face velocity

Notes of 2026-10-26. A bug in `Source/Riemann.H` (`riemann`, the approximate Godunov solver used by the MOL hydro flux, `MOL.cpp`) drove the mass fractions out of $[0,1]$ in the 2D CrossJet case and stopped pure PeleC with NaNs. This note gives the symptom, the diagnosis, the mechanism, the fix, and why the fix does not add numerical diffusion.

## 1. Symptom

2D CrossJet (`Exec/RegTests/CrossJet2D`, H$_2$ slot jet in a Mach 1.5 crossflow, two AMR levels, finest 0.078 mm), pure PeleC, no LES, no chemistry:

* from about 10 µs, cells with $Y_{H_2} > 1$ and $Y_{N_2} < 0$ appear in fast cold H$_2$ pockets downstream of the injector ($u \approx 1500$ m/s, local Mach about 1.2);
* the overshoot travels with the pocket (about 1100 m/s) and grows: $\max Y_{H_2}$ = 1.01 at 12 µs, 1.13 at 14 µs, 1.97 at 16 µs, with $\min Y_{N_2} = -0.74$;
* $\sum_k Y_k = 1$ throughout: the species are not lost, they are carried in negative and over-unity amounts;
* the time step collapses (from $10^{-8}$ to $1.7\times10^{-9}$ s) and the run stops at 75 µs: *"NaNs detected in the Solution"*.

## 2. Diagnosis

Runs on the half channel to 16 µs, maximum $Y_{H_2}$ over the domain:

| Test | Setting | $\max Y_{H_2}$ at 16 µs |
|---|---|---|
| default | second-order MOL, Godunov flux, AMR, diffusion on | 1.97 |
| reconstruction | `pelec.mol_iorder = 1` (no slopes, first order) | 1.97 |
| species diffusion | `pelec.diffuse_spec = 0`, `pelec.diffuse_enth = 0` | 1.92 |
| AMR interpolation | `pelec.state_interp_order = 0` (piecewise constant) | 1.98 |
| no AMR | one uniform level at the finest resolution (800 × 192) | 1.91 |
| subcycling | `amr.subcycling_mode = None` | 1.96 |
| reflux | `pelec.do_reflux = 0` | 1.97 |
| time step | `pelec.cfl = 0.05` (half) | 3.14 |
| flux | `pelec.use_laxf_flux = 1` (Lax-Friedrichs / Rusanov) | 1.002 |
| **fix** | **Godunov, species upwinded by the face velocity** | **1.0008** |

So it is neither the AMR machinery, nor the reconstruction, nor the diffusion, nor a time-step (positivity) limit: a first-order scheme on a single uniform level overshoots, and a smaller time step makes it worse. Only the choice of flux matters, which points to the Riemann solver.

## 3. How the solver builds the fluxes

For a face with left and right states $(\rho, u, p, Y_k)_{L,R}$ (`Riemann.H`, lines 40 to 186):

1. **Star state** (two-shock linearisation):
$$
u^* = \frac{W_L u_L + W_R u_R + p_L - p_R}{W_L + W_R}, \qquad
p^* = \frac{W_R p_L + W_L p_R + W_L W_R (u_L - u_R)}{W_L + W_R}.
$$
2. **Upwind side by the sign of $u^*$**: the state $o$ is $L$ if $u^* > 0$, $R$ otherwise:
$$
(\rho_o, u_o, p_o, Y_{k,o}) = \begin{cases}(\rho_L, u_L, p_L, Y_{k,L}) & u^* > 0\\ (\rho_R, u_R, p_R, Y_{k,R}) & u^* \le 0\end{cases}
$$
   and the star density on that side, $\rho^* = \rho_o + (p^* - p_o)/c_o^2$, with the same composition: $(\rho Y_k)^* = \rho^* Y_{k,o}$.
3. **Godunov state at the face**: if the wave on the $o$ side is a rarefaction, the face state is interpolated through the fan with a weight $f \in [0,1]$; supersonic cases take the $o$ state or the star state outright:
$$
u_{gd} = f\,u^* + (1-f)\,u_o, \qquad (\rho Y_k)_{gd} = f\,(\rho Y_k)^* + (1-f)\,(\rho Y_k)_o = \rho_{gd}\, Y_{k,o}.
$$
4. **Fluxes**:
$$
F_\rho = \rho_{gd}\,u_{gd}, \qquad F_{\rho Y_k} = (\rho Y_k)_{gd}\, u_{gd} = F_\rho\, Y_{k,o}.
$$

The composition carried by the species flux is $Y_{k,o}$, the side chosen by the **sign of $u^*$**. The direction of the mass flux is the **sign of $u_{gd}$**.

## 4. The mechanism

$u_o$ is the velocity of the side chosen by $u^*$, and it can have the opposite sign to $u^*$. For example, the left state moves left ($u_L < 0$) but the pressure jump makes $u^* > 0$. Inside a rarefaction fan $u_{gd} = f u^* + (1-f) u_o$ then changes sign with $f$: this is a **transonic rarefaction**, where the fan crosses $u = 0$. With $u^* > 0$ and $u_{gd} < 0$:

* the mass flows from right to left ($F_\rho < 0$),
* but it carries the left composition $Y_{k,L}$.

The right cell therefore loses $|F_\rho|\,Y_{k,L}\,\Delta t/\Delta x$ of species it may not contain. If the right cell is pure H$_2$ and the left one is air, it loses nitrogen it does not have:

$$
(\rho Y_{N_2})_R^{\,new} = \underbrace{(\rho Y_{N_2})_R}_{=0} - |F_\rho|\,Y_{N_2,L}\,\frac{\Delta t}{\Delta x} < 0 ,
$$

and since $\sum_k \rho Y_k = \rho$ is still conserved, $Y_{H_2} > 1$ in that cell. Nothing bounds the error: it grows step after step where cold supersonic H$_2$ meets air, as in the pockets shed by the jet. A smaller time step does not help because the defect is in the direction of the flux, not in its size relative to the cell content.

## 5. The fix

In an exact Godunov solution the species are passive scalars: they are discontinuous only across the contact and travel with the fluid, so the composition at a face is that of the side **the mass comes from**, i.e. upwind of the face velocity:

$$
F_{\rho Y_k} = F_\rho \; Y_{k,\mathrm{up}}, \qquad
Y_{k,\mathrm{up}} = \begin{cases} Y_{k,L} & u_{gd} > 0 \\ Y_{k,R} & u_{gd} < 0 \\ (\rho Y_k)_{gd}/\rho_{gd} & u_{gd} = 0 \end{cases}
$$

```cpp
// Riemann.H, end of riemann()
uflx_rho = rgd * qint_iu;
for (int n = 0; n < NUM_SPECIES; n++) {
  const amrex::Real Y_up =
    (qint_iu > 0.0) ? spl[n] : ((qint_iu < 0.0) ? spr[n] : rspgd[n] / rgd);
  uflx_rhoY[n] = uflx_rho * Y_up;
}
```

`spl`, `spr` are the reconstructed face states of `MOL.cpp`, normalised so that $\sum_k Y_k = 1$, hence $\sum_k F_{\rho Y_k} = F_\rho$ exactly.

With this flux every species update is a convex combination of upwind mass fractions, so $Y_k$ stays in $[0, 1]$ up to the second-order slopes: in the test above $\max Y_{H_2} = 1.0008$, the time step stays at $7.3\times10^{-9}$ s, and the half channel runs past 16 µs with no growth.

## 6. Why this is not a diffusive (upwind) scheme

* **The scheme is still Godunov.** $u^*$, $p^*$, the star and fan states, and the mass, momentum and energy fluxes are computed exactly as before. Shock positions and strengths, contact and rarefaction structure come from those and are unchanged. The second-order MUSCL reconstruction (`mol_iorder = 2`) is unchanged.
* **Only the composition riding on the mass flux changes, and only on some faces.** Where $\mathrm{sign}(u^*) = \mathrm{sign}(u_{gd})$, which is almost every face, $Y_{k,\mathrm{up}} = Y_{k,o}$ and the new flux is identical to the old one to the last bit. The two differ only at transonic rarefaction faces, where $u_{gd}$ is near zero and the mass flux is small.
* **It is the Godunov flux of a passive scalar.** Upwinding the scalar by the mass flux is what the exact Riemann solution prescribes for $Y_k$. The old flux was not more accurate there, only inconsistent.
* **The dissipative alternative is Lax-Friedrichs.** `pelec.use_laxf_flux = 1` also keeps $Y_k$ bounded (1.002), but it adds a dissipation proportional to the largest wave speed to every variable on every face, shocks included. That is the heavy diffusion to avoid, and it is not used.

## 7. Open points

1. **How many faces are affected** has not been counted. A short run comparing old and new species fluxes face by face would give the fraction of faces and of the mass flux concerned.
2. **The energy flux** still uses the internal energy of the Godunov state built with the $u^*$-side composition, $\rho e_{gd} = \rho_{gd}\, e(\rho_{gd}, Y_{k,o}, p_{gd})$. At the same transonic faces it carries the internal energy of a composition that is not the one moving. It did not cause visible trouble once the species were fixed, but it is the same inconsistency and could be treated the same way.
3. **The fix is in the shared `Source/Riemann.H`**, so it applies to every PeleC build, with or without CLEM. The CLEM runs made before 2026-10-26 used the old flux. With the LES species frozen they did not use the LES species fluxes, but their mass and energy fluxes come from the same solver; rebuild before any comparison with the quasi-DNS.
4. Worth reporting upstream to PeleC, with the 2D CrossJet case as reproducer.

## 8. Reproduce

```bash
cd Exec/RegTests/CrossJet2D
make -j16 CLEM_MODEL=FALSE                     # pure PeleC, 2D
mpiexec -n 32 ./PeleC2d.gnu.MPI.ex inputs-2d-amr2.inp stop_time=1.6e-5
# then the maximum of Y(H2) in plt* at 10, 12, 14, 16 us
```

Without the fix: 1.001, 1.012, 1.13, 1.97. With it: 1.0002, 1.0007, 1.0008, 1.0001.
