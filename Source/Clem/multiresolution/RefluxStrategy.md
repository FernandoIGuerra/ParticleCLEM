# Particle CLEM in PeleC: state of the code and the open reflux question

Notes written on 2026-10-12 to continue the discussion in another session or on another device. They are self-contained: Part 1 says what exists, Part 2 states the open problem (the reflux of the lines at coarse/fine interfaces), Part 3 lists what still has to be decided.

Repository: `https://github.com/FernandoIGuerra/ParticleCLEM.git`, branch `particle-clem` (a PeleC `development` clone plus `Source/Clem`).

---

## Part 1. What exists

### 1.1 The model in one paragraph

A combustion subgrid model for LES of any chemical mechanism. Every LES cell of **every AMR level**, covered or not, holds a one-dimensional LEM line of `n_lem = 12` AMReX particles (elements). The LES solves $\rho$, $\rho\mathbf{u}$, $\rho E$; the lines solve $Y_k$ and $T$ in a Lagrangian formulation. The lines are transported between cells by the LES mass flux (splicing), processed along the line (diffusion; stirring and reaction to come), and returned to the LES through a conserved-delta write-back.

Fixed choices: MOL time integrator only (`pelec.do_mol = 1`, `pelec.mol_iters = 1`), CVODE as chemistry integrator, CPU + MPI first.

### 1.2 Element state

Per particle, primitive: `rho, mass, vol, T, eint, press, phi, Y_0 … Y_{N-1}` and two integers `lem_index` (position along the line, 0 is the head) and `entry_side` (bookkeeping of the splice). `phi` is a dummy passive scalar for tests.

### 1.3 One step of one level (`PeleC::AdvanceClem`)

Called in place of `react_state` inside `do_mol_advance`, once per level per (sub)step, with that level's `dt`:

0. Halo of the level above filled from the lines of this level at the old time.
1. **ForwardCoupling**: every element goes isentropically from its own pressure to the LES cell pressure.
2. **DiffusionStirringReactor**: explicit sub-step loop along each line (diffusion only today).
3. **Splice**: lines follow the stored LES mass flux, in the order of the policy.
4. **Regrid**: every line back to `n_lem` equal-mass elements.
5. **BackwardCoupling**: volume renormalisation, then $(\rho Y_k)^{n+1} = (\rho Y_k)^n + \Delta[\rho Y_k]$ projected onto the LES density, and $(\rho e)^{n+1} = (\rho e)^n + \Delta[\rho e]$.

The step ends with the regrid, so **between steps every line of every level is regular** (12 equal-mass elements, `lem_index` 0..11). The halo and the average-down rely on this.

In `PeleC::post_timestep` of a level with a finer level above it (`PeleC::PostTimestepClem`): the mass reflux of Part 4, then **AverageDown**.

The LES does not advance the species: the species part of the MOL source is zeroed in both RK2 stages (`freeze_les_species`).

### 1.4 Files (`Source/Clem`)

| File | Role |
|---|---|
| `setup.H` | every compile-time choice: `n_lem`, models, policies, switches |
| `particle_index.H`, `particle_container.{H,cpp}`, `particle_plot.cpp` | storage, cell map, seeding, filter, plotfiles |
| `initialize_clem.cpp`, `advance_clem.cpp` | `PeleC::InitializeClem`, `PeleC::AdvanceClem` |
| `mass_flux.{H,cpp}` | per-level face mass flux of the MOL step, and its check |
| `coupling/` | `Coupling` base, `WriteBackCoupling` |
| `diffusion/` | `Diffusion` base, `MixtureAveragedDiffusion` |
| `reactor/react.cpp` | the sub-step loop `DiffusionStirringReactor` |
| `splicing/splicer.{H,cpp}` | the splice |
| `reindexing/` | ordering policies $\Sigma^+$, $\Sigma^-$, alternating |
| `regridding/` | `Regridder` base, `EqualMassRegridder` |
| `multiresolution/` | `HaloManager`, `AverageDown`, `PostTimestepClem`, this note |
| `line_scratch.H`, `les_comps.H`, `eos_util.H` | shared helpers |

Modules are CRTP classes selected in `setup.H`. Naming: PascalCase for functions and classes, snake_case for variables and file names.

### 1.5 The pieces that matter for the reflux

**Mass flux.** `clem_mass_flux[3]` is a per-level `PeleC` member, face-centred, one component: $\tfrac12(F^n + F^{n+1})$ with $F = \rho^* u^* A$, the RK2 average of the MOL stages. It closes the LES mass balance of every level and substep to round-off.

**Splice** (`clem::Splice`, one level at a time). For each cell, the outflux faces are served in the order $\mathcal{O}^{\omega}$, each cutting its segment from the **head** of what is left; the element across a cut is split and both pieces keep their state. Arriving blocks are appended at the **tail**, ranked by the same $\mathcal{O}^{\omega}$ (by size of the transfer). A face whose neighbour is not on the level (domain boundary **or coarse/fine boundary**) is a boundary of the lines: outflux leaves, influx enters as one block with the state of the **LES ghost cell**. The halo is not used by the splice yet.

**Halo** (`HaloManager`, one per fine level). Each fine cell position of the first layer around the fine grids that lies on the coarser level holds an exact copy of the line of its coarse parent $C$, same elements, same order, same state, with mass and volume divided by $r^d$ (8 in 3D).

**Average-down** (`clem::AverageDown`). Element $k$ of a covered cell is the average of element $k$ of its $r^d$ children: mass and volume summed, $e$, $Y_k$, $\phi$ mass-weighted, $T$ recovered from $(\rho, e, Y)$, pressure from the EOS. The regridder uses the same rule (average $e$, recover $T$).

### 1.6 What the tests showed

Reduced CrossJet grid (base 200×24×8, H2 jet in Mach 1.5 crossflow, BurkeDryer), `amr.regrid_int = -1`.

| Quantity | Result |
|---|---|
| Line mass vs LES cell mass, single level, 40 steps | $10^{-15}$ |
| Line mass vs LES, finest of three levels | $10^{-15}$ |
| H2 mass on the lines vs in the LES, single level | agree to 0.05–0.1 % |
| Covered cells directly under the finest level, mass, $Y_{H_2}$, $\rho e$ vs LES | $10^{-15}$ |
| Uncovered cells of level 1 next to level 2, line mass vs LES | drift, $5\times10^{-7}$ after 10 coarse steps |
| Covered cells of level 0, line mass vs LES | $10^{-7}$ (inherited from the level-1 drift) |

The drift of the uncovered coarse cells next to a fine level is the missing reflux: PeleC replaces the coarse flux on the interface by the fine ones, the lines do not.

### 1.7 Not yet in

Reaction, stirring, restart, input-file parameters (everything is a constant in `setup.H`). The species part of the reflux is deliberately not corrected (Part 4).

---

## Part 2. The open problem: reflux of the lines

### 2.1 Setting

One coarse step of level $\ell$ and the $r$ substeps of level $\ell+1$, in the order PeleC runs them (`Amr::timeStep`):

1. advance level $\ell$ over $\Delta t_c$ (its `AdvanceClem`),
2. advance level $\ell+1$ over $\Delta t_f$, $r$ times,
3. `post_timestep` of $\ell+1$, then of $\ell$ (LES reflux, LES average-down, then the lines).

A coarse cell $C$ next to the fine level exchanges through its interface face $F$ twice, with different amounts: once in its own splice with the coarse transfer $\mu_F$ (a provisional exchange with the covered cell behind the face), and through the fine cells, over their substeps, with the fine transfers $\mu^s_q$ on the $r^{d-1}$ sub-faces $q$ of $F$. The fine level is the truth, as in the LES reflux.

### 2.2 The reflux of AMREX_CLEM, and why it does not carry over

The prototype (`amrex-tutorials/ExampleCodes/Particles/AMREX_CLEM`, one passive scalar) does this per coarse cell $C$:

- register the mass and $m\phi$ that crossed each interface face on both levels, giving the mismatches $dm$ and $d\Phi$;
- mass: cut $|dm|$ from the head of $C$, or put it back at the head with the state of the head element;
- scalar: shift every element of $C$ by the same $d\phi$ so that $\sum m\phi$ matches;
- clip to $[0,1]$ and renormalise; regrid.

Three failures, all worse with $N$ species plus energy:

**Penetration.** The fine level takes its inflow from the head of a ghost copy of $C$ that holds $m_C/r^d$. Over one coarse step a fine cell takes $X = \sum_s \mu^s_q$ from its ghost, a fraction $r^d X / m_C$ of that ghost's line. The $r^{d-1}$ fine cells on the face together remove $r^{d-1} X$ from $C$, a fraction $r^{d-1} X / m_C$ of the real line. The fine level therefore draws from a depth $r$ times larger than what actually leaves $C$ (twice, for $r=2$), and draws the same depth range $r^{d-1}$ times over. What the fine level received is not what $C$ lost.

**Non-conservation of $\delta\phi_k$.** Because of that mismatch the uniform shift is needed, and the shift can push a mass fraction outside $[0,1]$. The clip and renormalisation that follow destroy the conservation the shift was meant to restore. With one scalar this is a small reported loss. With 13 species, trace radicals go negative almost every time, the renormalisation moves the major species, and a shifted energy can give a temperature the element cannot have.

**Structure.** The shift spreads the correction over all the elements of $C$, although the exchange happened at its head. A front crossing the interface is smeared over the whole coarse cell.

### 2.3 Candidate strategies

**A. Keep the reflux as a correction (register, cut/add, shift, clip).** Known to conserve mass; fails for the scalars as above. Listed for completeness.

**B. The coarse line itself is shared (section 10 of `MultiLevelSplicingTheory.md`, the "fine-driven interface").** The transfers are partitioned into two disjoint groups: $\mathcal{G}_f$ (fine faces and interface sub-faces) and $\mathcal{G}_c$ (coarse faces that are not on the interface). Each transfer is executed exactly once per coarse step. The fine splice treats $C$ as an ordinary member at sub-face granularity: $C$ serves its outflux sub-faces from its own head, in the order $\mathcal{O}^{\omega}$, and appends the blocks arriving from the fine cells at its tail. The coarse splice skips the interface faces.

- No ghost, no register, no shift, no clip.
- Only cuts and concatenations, so every $Y_k$ and $e$ is conserved exactly and every value is that of a real element: positivity and $\sum_k Y_k = 1$ hold by construction.
- Penetration is the true one: the sub-faces take consecutive ranges from the same head.
- Mass of $C$ after the step equals the refluxed LES mass, because the refluxed LES is exactly "own flux on non-interface faces, fine fluxes on the interface".

**C. One ghost shared by the fine cells of a face.** The fine cells facing $C$ draw consecutively from a single copy of $C$'s head segment instead of each owning a $1/r^d$ clone. This fixes the depth and the duplicated range, but $C$ and the copy still disagree when fine and coarse transfers differ, so a correction of $C$ remains.

### 2.4 What strategy B costs inside PeleC

- **Order.** PeleC advances the coarse level first, then the fine substeps: the reverse of section 10.4 of the notes. $C$ would do its coarse faces, and its regrid, before it takes part in the interface exchanges.
- **Communication.** $C$ lives on the coarse grids, possibly on another rank, and each fine substep must cut from it and append to it.
- **A second write-back.** $C$ changes after its own step was already written to the LES. The cells next to the interface need a write-back at the sync, and PeleC's own reflux of species and energy there must be switched off (today PeleC still sends those fluxes to its flux registers).
- **An extra regrid.** $C$ is irregular after the fine substeps and must be regridded at the sync, before the average-down. The energy change of that regrid has to be inside the write-back of the sync.
- **The role of the halo.** With a shared line the halo is no longer the donor. It could remain as a read-only prediction of $C$'s head for the fine level, or go.
- **Admissibility.** $C$ must hold enough mass for its coarse outflux and for all the interface sub-face outfluxes of the $r$ substeps.

---

## Part 3. Open questions

1. **What was the "shared" strategy?** The last message of the session ended at "a potential strategy was a shared…". B and C above are the two readings; it may be neither.
2. **If B:** how does $C$ take part in the fine splice across ranks; is the halo kept; where exactly does the second write-back go; how are species and energy removed from the LES reflux.
3. **If a correction is kept:** what replaces the uniform shift for $Y_k$ and for the thermal variable so that the result stays physical and conserved.
4. **Ordering inside $C$.** In B the face schedule of $C$ is split in time (coarse faces first, interface sub-faces after), so it is no longer one $\mathcal{O}^{\omega}$ over all its faces. Is that acceptable for the three policies.
5. **Level-0 energy mismatch.** Covered cells of level 0 differ from the LES $\rho e$ by up to 0.3 % while level 1 is exact. Hypothesis, not verified: after its own average-down PeleC recomputes the internal energy of the level-1 covered cells from total minus kinetic energy, so the LES values that level 0 is compared against have moved.
6. **The $8\times10^{-8}$ between stored and EOS energy** at the 1000 K crossflow, attributed (not confirmed) to the switch between the two temperature ranges of the thermodynamic fits.
7. **After the reflux:** lines following AMR regrids, restart, then reaction and stirring.

---

## Part 4. Decision of 2026-10-13: mass-only reflux (implemented)

Fixing $\delta Y_k$ is not trivial and a shift can create negative species, so the first version corrects **only the mass** and leaves the composition as it falls. The purpose is to see what the unfixed $\delta Y_k$ does in a complex flow.

**Rule.** The mismatch on an interface face, $\delta\mu_f = \sum_s \mu^s_{\text{fine}} - \mu_{\text{coarse}}$, is treated as a transfer still to be made and handled as a splice between $C$ and the coarse cell $N$ that covers the fine cells behind the face:

- $C$ gave too little: $\delta\mu_f$ is cut from the **head** of $C$ and dropped (the fine cells already hold it).
- $C$ needs more: $\delta\mu_f$ is taken from the **head** of $N$ and appended at the **tail** of $C$. $N$ is only copied from.
- Several interface faces: served and ranked in the order of the policy, as in the splice.
- A mismatch larger than the line it comes from aborts the run.

**Order at the sync** (`PeleC::PostTimestepClem`): mass reflux, regrid of the touched cells only, conserved-delta write-back in those cells, then the average-down. The reflux therefore reads $N$ as the coarse splice left it, before the average-down overwrites it.

**Where $\delta\mu$ comes from.** Each level accumulates `dt * clem_mass_flux` over its substeps in `clem_mass_transfer`; at the sync it is summed onto the coarse faces and the coarse transfer subtracted.

**The LES.** The species fluxes are removed from PeleC's flux registers, so PeleC refluxes $\rho$, $\rho\mathbf{u}$ and $\rho E$ only. In the touched cells $\rho Y_k \mathrel{+}= \Delta[\rho Y_k]_{\text{reflux}}$, projected onto the refluxed density, and $\rho e \mathrel{+}= \Delta[\rho e]_{\text{reflux}}$, both from the lines alone. $\rho E$ is the one PeleC refluxed.

**Measured** (reduced CrossJet, three levels, 10 coarse steps): after every sync the line mass equals the LES cell mass to $10^{-15}$ on every level; the projection is $10^{-16}$, so the mass the lines exchanged is exactly the mass PeleC refluxed. The mismatches were small, $|dm|/m \le 4\times10^{-8}$, and hydrogen had not reached an interface, so **the effect of the unfixed $\delta Y_k$ has not been observed yet**.

**Known caveat.** After the sync PeleC recomputes the internal energy from $\rho E$ minus the kinetic energy whenever that is not tiny (`dual_energy_eta2`), so in practice the temperature of the touched cells follows PeleC's energy reflux and the $\Delta[\rho e]_{\text{reflux}}$ written to $\rho e$ is overwritten.

Files: `multiresolution/refluxer.{H,cpp}`, `multiresolution/post_timestep.cpp`; the regridder takes an optional mask of cells.

Part 2 and Part 3 remain as the record of the alternatives if the mass-only reflux proves insufficient.

---

## Part 5. Lines following an AMR regrid (implemented 2026-10-13)

**Rule.** A new fine cell gets the clone of the line of its coarse parent (same elements, same order, same state, mass and volume divided by $r^d$), exactly as the halo. A fine cell that leaves its level has its line dropped; the coarse cell under it already holds the average. A cell that stays is refiled under the new grids. Then every covered cell is averaged down again.

**The LES must match.** In a new fine cell the LES takes the state of its parent unchanged, **piecewise constant for every component** (`PeleC::SetNewCellsFromCoarse`, called where a level is rebuilt). A linear interpolation gives the children different densities, which identical clones cannot follow. The ghost-cell fills during the steps keep PeleC's usual interpolation.

**Where.** `clem::AdaptToGrids` (`multiresolution/grid_adapter.{H,cpp}`), called from `PeleC::post_regrid` of the base level. It does nothing when the grids are those the lines are already filed under, which is most regrids.

**Measured** (reduced CrossJet, two levels, 40 steps, 4 ranks, refinement following the density gradient): one regrid grew the fine level from 2,048 to 3,072 cells; 1,024 new cells were cloned; right after it the line mass matched the LES to $2\times10^{-15}$ on both levels, and the run continued with the same accuracy. **Not exercised yet:** fine cells being removed, a level appearing or disappearing, three levels with regridding.

---

## Appendix. Reference material

| What | Where |
|---|---|
| Single-level reference with full physics | `FPCE/CLEMNew/PeleC/Source/Clem`, notes in `CLEMNew/PeleC/InformationImplementation/CouplingLESwithCLEM.md` |
| AMR strategy prototype, passive scalar | `FPCE/amrex/Tutorials/amrex-tutorials/ExampleCodes/Particles/AMREX_CLEM/Source` (`FillHaloParticle`, `RefluxParticle`, `AvgDownParticle`, `ClemAmr`) |
| Splicing theory | `…/Particles/ClemAdvection/OperatorClem/AdvectionSplicingTheory.md` |
| Multilevel theory, incl. section 10 | `…/Particles/ClemAdvection/OperatorClem/MultiLevelSplicingTheory.md` |

Build and run the test case:

```bash
cd Exec/RegTests/CrossJet
make -j8                      # CLEM_MODEL = TRUE in the GNUmakefile
mpiexec -n 4 ./PeleC3d.gnu.MPI.ex inputs-case1-clem-init.inp \
  amr.n_cell="200 24 8" max_step=10 amr.regrid_int=-1 \
  pelec.do_react=1 pelec.chem_integrator=ReactorCvode cvode.solve_type=denseAJ_direct
```

`pelec.do_react = 1` is needed because the CLEM step sits in PeleC's reaction branch. The full base grid (400×48×16) needs roughly 3 GB for the particles alone.
