/**
 * probe.cpp: a diagnostic of one LES cell and its line through the CLEM step.
 *
 *   PeleC::ProbeClem  prints, at the stage named by `tag`, the LES state of
 *                     the probed cells and the filter of their lines
 *
 * Off unless the input file gives the cells, all on one level:
 *
 *   clem.probe_level     = 2
 *   clem.probe_cells     = 544 0 32   543 0 32      (i j k of each cell)
 *   clem.probe_from_step = 70                       (coarse step, default 0)
 */
#include <AMReX_ParmParse.H>

#include "PeleC.H"
#include "IndexDefines.H"

#include "particle_container.H"

void
PeleC::ProbeClem(const char* tag)
{
  static int probe_level = -1;
  static int from_step   = 0;
  static amrex::Vector<int> cells;
  static bool read = false;
  if (!read) {
    amrex::ParmParse pp("clem");
    pp.query("probe_level", probe_level);
    pp.query("probe_from_step", from_step);
    pp.queryarr("probe_cells", cells);
    read = true;
  }
  if (level != probe_level || parent->levelSteps(0) < from_step || ClemPC == nullptr) {
    return;
  }

  auto eos                       = pele::physics::PhysicsType::eos();
  const amrex::MultiFab& S_new   = get_new_data(State_Type);
  const amrex::MultiFab& S_old   = get_old_data(State_Type);
  const amrex::Real cell_vol     = AMREX_D_TERM(geom.CellSize(0), *geom.CellSize(1), *geom.CellSize(2));
  auto& pc                       = ClemPC->Particles();
  if (!ClemPC->IsCellMapValid(level)) {
    ClemPC->RebuildCellMap(level);
  }

  const auto print_les = [&](const char* name, const amrex::Array4<const amrex::Real>& s,
                             const amrex::IntVect& iv) {
    const amrex::Real rho = s(iv, URHO);
    amrex::Real Y[NUM_SPECIES];
    amrex::Real rho_Y_sum = 0.0;
    for (int k = 0; k < NUM_SPECIES; ++k) {
      Y[k] = s(iv, UFS + k) / rho;
      rho_Y_sum += s(iv, UFS + k);
    }
    const amrex::Real kinetic =
      0.5 * (s(iv, UMX) * s(iv, UMX) + s(iv, UMY) * s(iv, UMY) + s(iv, UMZ) * s(iv, UMZ)) / rho;
    amrex::Real press = 0.0;
    eos.RTY2P(rho, s(iv, UTEMP), Y, press);
    amrex::AllPrint() << "PROBE " << name << " rho " << rho << " sumY " << rho_Y_sum / rho
                      << " Y_H2 " << Y[H2_ID] << " T " << s(iv, UTEMP) << " e_int "
                      << s(iv, UEINT) / rho << " e_tot-ke " << (s(iv, UEDEN) - kinetic) / rho
                      << " ke " << kinetic / rho << " p " << press << '\n';
  };

  for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.validbox();
    for (int c = 0; c + AMREX_SPACEDIM <= static_cast<int>(cells.size()); c += AMREX_SPACEDIM) {
      const amrex::IntVect iv(AMREX_D_DECL(cells[c], cells[c + 1], cells[c + 2]));
      if (!bx.contains(iv)) {
        continue;
      }
      amrex::AllPrint() << "PROBE ---- step " << parent->levelSteps(level) << " level " << level
                        << " cell " << iv << " : " << tag << '\n';
      print_les("  LES old ", S_old.const_array(mfi), iv);
      print_les("  LES new ", S_new.const_array(mfi), iv);

      const auto& tile = pc.ParticlesAt(level, mfi.index(), mfi.LocalTileIndex());
      const auto* particles = tile.GetArrayOfStructs().data();
      amrex::Real mass = 0.0, vol = 0.0, m_e = 0.0, m_h2 = 0.0, m_T = 0.0;
      amrex::Real T_min = 1.0e30, T_max = -1.0e30, p_min = 1.0e30, p_max = -1.0e30;
      int n_el = 0;
      for (const int n : ClemPC->ParticlesInCell(level, mfi.index(), iv)) {
        const auto& p       = particles[n];
        const amrex::Real m = p.rdata(clem::RealData::mass);
        mass += m;
        vol += p.rdata(clem::RealData::vol);
        m_e += m * p.rdata(clem::RealData::eint);
        m_T += m * p.rdata(clem::RealData::T);
        m_h2 += m * p.rdata(clem::RealData::Y0 + H2_ID);
        T_min = amrex::min(T_min, p.rdata(clem::RealData::T));
        T_max = amrex::max(T_max, p.rdata(clem::RealData::T));
        p_min = amrex::min(p_min, p.rdata(clem::RealData::press));
        p_max = amrex::max(p_max, p.rdata(clem::RealData::press));
        ++n_el;
      }
      amrex::AllPrint() << "PROBE   line    n " << n_el << " rho " << mass / cell_vol
                        << " vol/V " << vol / cell_vol << " Y_H2 " << m_h2 / mass << " T "
                        << T_min << " / " << m_T / mass << " / " << T_max << " e_int "
                        << m_e / mass << " p " << p_min << " / " << p_max << '\n';
    }
  }
}
