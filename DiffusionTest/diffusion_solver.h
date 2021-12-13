#ifndef DIFFUSION_SOLVER_H
#define DIFFUSION_SOLVER_H

#include "ChiPhysics/SolverBase/chi_solver.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/UnknownManager/unknown_manager.h"
#include "ChiMath/PETScUtils/petsc_utils.h"

#include "FiniteElementSpace/fe_space.h"

//################################################################### Class def
/**Test diffusion solver for different methods.*/
class DiffusionSolver : public chi_physics::Solver
{
protected:
  chi_math::finite_element::SpatialDiscretizationPtr m_sdm = nullptr;
  chi_mesh::MeshContinuumPtr m_grid = nullptr;
  chi_math::UnknownManager   m_dof_handler;

  Mat A = nullptr;
  Vec x = nullptr;
  Vec b = nullptr;

public:
  DiffusionSolver();

  void SetSpatialDiscretization(
    const chi_math::finite_element::SpatialDiscretizationPtr& sdm)
  {
    m_sdm = sdm;
  }

  void Initialize() override;
  void Execute() override {}

  virtual ~DiffusionSolver() = default;
};

#endif //DIFFUSION_SOLVER_H