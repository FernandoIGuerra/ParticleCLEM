# Energy coupling between the LES and the CLEM lines

Notes of 2026-10-15. Why the CLEMNew energy write-back fails under AMR, why classical PeleC does not, and the two changes that make the 2D CrossJet case run: the line energy anchor and the removal of the species-enthalpy diffusion flux from the LES. Companion of `multiresolution/EnergyCouplingProblem.pdf` (the 3D diagnosis).

Notation: $U = [\rho, \rho\mathbf{u}, \rho E, \rho Y_k]$, $E = e + k$, $k = \tfrac12|\mathbf{u}|^2$. Superscript $n$ is the start of the step, $^{**}$ the end of the RK2 (MOL) LES step, $n+1$ the end of the full step. On the lines, element $p$ of cell $c$ has mass $m_p$, volume $v_p$, $e_p$, $Y_{k,p}$, $T_p$, $P_p$. The conserved filter of a line is

$$
[\rho\phi]_c = \frac{1}{V_c}\sum_{p\in c} m_p \phi_p .
$$

---

## 1. Why classical PeleC works

PeleC solves, in every cell, one conservative system with one set of face fluxes:

$$
\partial_t(\rho Y_k) + \nabla\cdot(\rho Y_k\mathbf{u}) = -\nabla\cdot\mathbf{J}_k ,
$$

$$
\partial_t(\rho E) + \nabla\cdot\big((\rho E + P)\mathbf{u}\big) = \nabla\cdot(\boldsymbol{\tau}\cdot\mathbf{u}) + \nabla\cdot(\lambda\nabla T) - \nabla\cdot\Big(\sum_k h_k \mathbf{J}_k\Big) .
$$

In `Diffterm.H` (`SpeciesEnergyFlux`) every species diffusion flux $\mathbf{J}_k$ is computed together with the energy it carries, $h_k\mathbf{J}_k$, on the same face, from the same face state. Three consequences:

1. **Energy and composition travel together.** Whatever moves H2 (advection, diffusion) moves the energy of that H2 in the same step, through the same face. The temperature recovered from $(\rho, e, Y_k)$ is therefore consistent: a cell never holds the energy of an H2-air mixture with the composition of air.
2. **The irreversible terms are in the energy equation.** Viscous heating $\nabla\cdot(\boldsymbol{\tau}\cdot\mathbf{u})$, conduction $\nabla\cdot(\lambda\nabla T)$ and the numerical dissipation of the hydro scheme at shocks and unresolved expansions all land in $\rho E$, and nothing overwrites $\rho E$ afterwards.
3. **Pressure, density and energy are one state.** $P = P(\rho, e, Y_k)$ in every cell, so there is no second description of the same gas to drift away from.

The wake behind the injector (a pressure drop by a factor of about 5 across the downstream edge of the slot in 2D, reverse flow at up to Mach 0.9 along the wall) is real and resolved; PeleC handles it with this system.

CLEM splits that system in two, the LES ($\rho$, $\rho\mathbf{u}$, $\rho E$) and the lines ($Y_k$, $T$), and every failure below is one of the three properties above being broken by the split.

---

## 2. The CLEM step

With `clem::setup::freeze_les_species = true` the LES does not advance $\rho Y_k$: the species part of the MOL source is zeroed in both RK2 stages,

$$
U^{**} = \mathrm{RK2}(U^n), \qquad (\rho Y_k)^{**} = (\rho Y_k)^n .
$$

Then `PeleC::AdvanceClem` (`advance_clem.cpp`):

1. **Forward coupling (pressure).** Every element goes isentropically to the LES pressure $P^{**} = P(Q^{**})$:
$$
\rho_p \leftarrow \rho_p\Big(\frac{P^{**}}{P_p}\Big)^{1/\gamma_p},\qquad
T_p \leftarrow T_p\Big(\frac{P^{**}}{P_p}\Big)^{(\gamma_p-1)/\gamma_p},\qquad
P_p \leftarrow P^{**}.
$$
2. **Line processes** $\mathcal{L}$ (`DiffusionStirringReactor`): diffusion along the line now, stirring and reaction later.
3. **Splice**: the lines follow the LES face mass fluxes.
4. **Regrid** to $n_{lem}$ equal-mass elements (conserves $m$, $v$ and $\sum m_p e_p$).
5. **Volume renormalisation**: isentropic rescale by $f = V_c / \sum_p v_p$.
6. **Backward coupling (write-back)** of $\rho Y_k$ and $\rho e$.
7. **Energy anchor** of the lines (section 5), new.

Species write-back (unchanged):

$$
(\rho Y_k)^{n+1} = \Big[(\rho Y_k)^n + [\rho Y_k]_c^{\text{post}} - [\rho Y_k]_c^{\text{pre}}\Big]\,
\frac{\rho^{**}}{\sum_j \big[(\rho Y_j)^n + [\rho Y_j]_c^{\text{post}} - [\rho Y_j]_c^{\text{pre}}\big]} ,
$$

with pre/post taken at the start and the end of the whole CLEM step.

---

## 3. Energy write-back, first form (CLEMNew, `EnergyCoupling::LinesOnly`)

$$
(\rho e)^{n+1} = \underbrace{(\rho E)^n - k^n}_{(\rho e)^n} + \Delta[\rho e], \qquad
\Delta[\rho e] = [\rho e]_c^{\text{post}} - [\rho e]_c^{\text{pre}} \ \text{(whole CLEM step)},
$$

$$
(\rho E)^{n+1} = (\rho e)^{n+1} + k^{n+1}.
$$

$(\rho E)^{**}$ does not appear: the LES energy update of the step is replaced by the line energy. The lines carry only transport (splice) and isentropic $P\,dV$ work (steps 1 and 5), so property 2 of section 1 is lost: viscous heating, conduction and shock dissipation are discarded every step.

**Failure (3D, two levels).** In the wake cell at the rim of the injector the lines lose about $2.4\times10^{8}$ erg/g (about 2.4 K) per substep relative to the LES energy equation. Lower energy gives lower $P$, more expansion, more isentropic cooling: the cell slides from 1000 K to 163 K and from 1 atm to 0.1 atm, then the LES returns NaN. Stable at zero and one refinement level (the CLEMNew resolution), unstable at 0.078 mm where the injector lip is resolved.

---

## 4. Energy write-back, second form (`EnergyCoupling::LesPlusSubgrid`, default)

The LES keeps what its energy equation did; the lines add only what the LES does not contain:

$$
(\rho e)^{n+1} = (\rho E)^{**} - k^{**} + \Delta[\rho e]_{\mathcal{L}}, \qquad
\Delta[\rho e]_{\mathcal{L}} = [\rho e]_c\big|_{\text{after }\mathcal{L}} - [\rho e]_c\big|_{\text{before }\mathcal{L}},
$$

$$
(\rho E)^{n+1} = (\rho e)^{n+1} + k^{**}.
$$

The bracket surrounds the line processes only (step 2), not the splice nor the pressure coupling, so transport and $P\,dV$ work are counted once, by the LES. At the coarse/fine sync the energy moved by the reflux is the LES's own (PeleC refluxes $\rho E$); the line reflux no longer touches $\rho e$.

This restores property 2 for the LES, but two new failures appeared.

### 4a. The lines drift (fixed by the anchor, section 5)

The lines still evolve their energy by transport and isentropic work only, and nothing ties them back to the LES. In 2D the line-vs-LES temperature gap grew to 560 K, then to $10^4$ K, the volume defect to 97 %; diffusion on those lines stopped conserving energy, $\Delta[\rho e]_{\mathcal{L}}$ became huge, the LES jumped by 1000 K, then 900 000 K, and a NaN element reached the regrid (*"a new element did not get its equal share of the mass"*).

### 4b. Energy and composition separate (fixed by section 6)

$\rho E$ is transported by the LES, $Y_k$ by the lines. Property 1 of section 1 is lost wherever the two transports differ. The largest difference was not numerical but a missing pairing, see section 6.

---

## 5. The line energy anchor (`AnchorLineEnergy`, `setup::anchor_line_energy = true`)

After the write-back, the LES holds the cell energy $e_c^{n+1} = (\rho e)^{n+1}/\rho^{n+1}$. Every element of the line gets the same shift

$$
\delta e_c = e_c^{n+1} - \frac{\sum_p m_p e_p}{\sum_p m_p}, \qquad e_p \leftarrow e_p + \delta e_c ,
$$

then $T_p = T(\rho_p, e_p, Y_{k,p})$ and $P_p = P(\rho_p, T_p, Y_{k,p})$ at the element's own density.

* The mass-mean energy of the line equals the LES energy of the cell, exactly, at the end of every step.
* The differences between the elements (the sub-grid structure built by diffusion, and later by stirring and reaction) are kept: the shift is uniform.
* The LES is not changed, so the LES energy stays conservative.
* $|\delta e_c|/c_v$ is the energy the lines miss in one step (irreversible heating, conduction). It is reported as *"line energy anchored to the LES, max shift"*: in the 2D run 2 to 3 K per step on level 2, 15 to 30 K on levels 0 and 1.

Together with the forward coupling this closes the loop in both directions: pressure from the LES to the lines at the start of the step, energy from the LES to the lines at the end, the sub-grid energy increment $\Delta[\rho e]_{\mathcal{L}}$ from the lines to the LES. The two descriptions now agree at $t^{n+1}$ in mass, energy and composition.

A gap between the mass-mean line temperature and the LES temperature remains (up to about 25 K on level 2) where a line holds both H2 and air elements: $\sum_p m_p T_p/\sum_p m_p \neq T(\bar\rho, \bar e, \bar Y)$ because $T(e, Y)$ is nonlinear in $Y$ (the heat capacity of H2 is about 14 times that of air). This is a property of the sub-grid structure, not a drift.

---

## 6. The species-enthalpy flux (`Diffterm.H`, `pc_les_species_diffuse`)

**The bug.** With the species frozen, `Diffusion.cpp` zeroes the species diffusion fluxes $\mathbf{J}_k$ after they are computed, but `SpeciesEnergyFlux` had already added their energy to the energy flux. The LES therefore solved

$$
\partial_t(\rho Y_k) + \dots = 0 \qquad\text{(species frozen)},
$$

$$
\partial_t(\rho E) + \dots = \dots - \nabla\cdot\Big(\sum_k h_k \mathbf{J}_k\Big) \qquad\text{(their energy still diffusing)}.
$$

The energy of H2 diffused into the air while the H2 stayed where the lines put it. In a cell whose LES energy is that of an H2-air mixture but whose composition (from the lines) is nearly pure air, $T(e, Y)$ is thousands of kelvin, because H2 holds about 14 times more energy per kelvin than air. In the 2D run this hot spot sat at the wall, at the H2/air front of the wake, exactly where pure PeleC has its diffuse H2 tail: 1666 K at 12 µs, 3449 K at 18 µs, then broken lines.

**The fix.** When the CLEM owns the species, $\mathbf{J}_k$ and $h_k\mathbf{J}_k$ are dropped together:

$$
\mathbf{J}_k = 0 \quad\text{and}\quad \sum_k h_k\mathbf{J}_k = 0 \qquad\text{in the LES, when } \texttt{freeze\_les\_species},
$$

```cpp
// Diffterm.H
constexpr bool pc_les_species_diffuse()   // false when clem::setup::freeze_les_species
...
if (update && pc_les_species_diffuse()) {
  FluxTypes::SpeciesEnergyFluxType()(iv, ivm, dxinv[dir], coef, q, flx);
} else if (!pc_les_species_diffuse()) {
  for (int ns = 0; ns < NUM_SPECIES; ++ns) flx(iv, UFS + ns) = 0.0;
}
```

in both `pc_diffusion_flux` and `pc_diffusion_flux_eb`. Pure PeleC builds (no `PELE_USE_CLEM`) are unchanged. Species diffusion and the energy it carries now both belong to the lines, as in PeleC they both belong to the same face flux.

---

## 7. Results (2D CrossJet, two levels, finest 0.078 mm, impulsive jet)

`Exec/RegTests/CrossJet2D`, input `inputs-2d-amr2.inp`, runs on the Desktop.

| Run | End | Hottest LES cell | Notes |
|---|---|---|---|
| Pure PeleC | NaN at 75 µs (step 7599) | 1344 to 1402 K (upstream separation) | NaN at x = 5.96 cm, not in the wake, section 8 |
| CLEM, `LesPlusSubgrid` | 5.7 µs | &ndash; | lines drift (4a) |
| CLEM, + anchor | 20 µs | 3449 K at 18 µs, wall, H2 front | species-enthalpy flux (4b, 6) |
| CLEM, + anchor + section 6 | completed step 7600 (76.4 µs), dt still 1e-8 s | 1344 K at 30 µs, same cell as pure PeleC; 1419 K at the end | $0 \le Y_{H_2} \le 1$ throughout; max line vs LES gap over the run 27 / 58 / 178 K (levels 0 / 1 / 2); max anchor shift 52 / 40 / 8 K |

Along the wall (`CrossJet2D_CLEM_fix/wall_line_vs_pure.png`, 10 to 74 µs) pressure, temperature, density, velocity and Mach number of CLEM and pure PeleC agree closely: the upstream separation shock, the drop across the slot, the reverse flow in the wake, and the passage of the H2 front at 10 µs. The differences are in the H2 fraction just behind the slot at early time (0.92 against 0.99) and in the diffuse H2 tail behind the front, which CLEM does not have (section 8, point 1), and in the details of the unsteady wake downstream of 4.7 cm at 74 µs.

The CLEM run passes the step at which pure PeleC stops: its species come from the lines, which keep $0 \le Y_k \le 1$ by construction, while pure PeleC's LES species overshoot (section 8, point 2).

---

## 8. Open points

1. **Inter-cell species diffusion is missing in CLEM.** The lines diffuse only within their own cell; PeleC also diffuses between cells. Behind the slot the CLEM wall H2 fraction is 0.87 to 0.95 against 0.99 in PeleC, and the H2 front is sharper (no diffuse tail). At 0.078 mm in a slow wall wake this is a resolved process, not a sub-grid one. A passive test (LES = pure PeleC, lines carried without feedback, line H2 compared with LES H2 cell by cell) would measure it.
2. **Pure PeleC itself breaks species bounds.** $Y_{H_2}$ reaches 1.97 and $Y_{N_2}$ $-0.74$ on level 2 at 16 µs, again 1.46 at 72 µs, in fast H2 pockets at x = 5.57 to 5.79 cm; the time step collapses and the run stops with NaN at step 7599. Piecewise-constant coarse/fine interpolation (`pelec.state_interp_order = 0`) does not remove it, so it comes from the fine-level solver. Not investigated further yet.
3. **The 3D case** has not been rerun with sections 5 and 6.
4. **Energy conservation of the lines.** The anchor injects $\sum_c \rho_c V_c\,\delta e_c$ into the lines each step; the LES is conservative, the lines are not by themselves. With reaction on the lines, $\delta e_c$ must not absorb the heat release: reaction is inside $\mathcal{L}$, so its energy goes to the LES through $\Delta[\rho e]_{\mathcal{L}}$ and the anchor then matches it.

## Switches (`setup.H`)

| Switch | Default | Meaning |
|---|---|---|
| `freeze_les_species` | `true` | lines own the species; with section 6 also their diffusion enthalpy flux |
| `energy_coupling` | `LesPlusSubgrid` | `LinesOnly` = CLEMNew; `Off` = diagnostic, LES energy only |
| `anchor_line_energy` | `true` | section 5 |
