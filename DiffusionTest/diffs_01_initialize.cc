#include "diffusion_solver.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include "FiniteElementSpace/fe_space.h"

#include "chi_log.h"

//###################################################################
/**Initialization routine.*/
void DiffusionSolver::Initialize()
{
  auto& chi_log = ChiLog::GetInstance();
  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log() << "Initializing diffusion solver";
  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Get grid
  auto handler = chi_mesh::GetCurrentHandler();
  m_grid = handler->GetGrid();

  //============================================= Initialize spatial
  //                                              discretization
  const auto sdm_string_option = basic_options("SDM").StringValue();

  if (m_sdm == nullptr)
    throw std::logic_error("Diffusion solver: Spatial Discretization not set.");

  //============================================= Initialize unknown manager
  m_dof_handler.Clear();
  m_dof_handler.AddUnknown(chi_math::UnknownType::SCALAR);
  m_dof_handler.SetUnknownTextName(0, "Temperature");

  //============================================= Get number of DOFs
  int64_t num_local_dofs = m_sdm->NumLocalDOFs(m_dof_handler);
  int64_t num_globl_dofs = m_sdm->NumGlobalDOFs(m_dof_handler);

  //============================================= Make PETSc items
  if (x != nullptr) VecDestroy(&x);
  if (b != nullptr) VecDestroy(&b);
  if (A != nullptr) MatDestroy(&A);

//  if (sdm_string_option == "PWLC")
  {
    using namespace chi_math::PETScUtils;
    x = CreateVector(num_local_dofs, num_globl_dofs);
    b = CreateVector(num_local_dofs, num_globl_dofs);

    const auto sparsity_pattern = m_sdm->BuildCFEMSparsityPattern(m_dof_handler);
    const auto& nnz_in_diag = sparsity_pattern.first;
    const auto& nnz_off_diag = sparsity_pattern.second;

    A = CreateSquareMatrix(num_local_dofs, num_globl_dofs);

    InitMatrixSparsity(A, nnz_in_diag, nnz_off_diag);
  }

  //============================================= Assemble system
  typedef std::vector<double> Row;
  typedef std::vector<Row> Matrix;

  for (const auto& cell : m_grid->local_cells)
  {
    const size_t num_nodes  = m_sdm->GetCellNumNodes(cell);
    const auto node_indices = m_sdm->CellNodeIndices(cell);

    const auto dof_global_indices = m_sdm->MapDOFsGlobal(cell, m_dof_handler, 0, 0);

    Matrix              cell_matrix;
    std::vector<double> cell_rhs;

    cell_matrix.assign(num_nodes, Row(num_nodes, 0.0));
    cell_rhs.assign(num_nodes, 0.0);

    for (const size_t i : node_indices)
      for (const size_t j : node_indices)
        cell_matrix[i][j] += 1.0;

    for (const size_t i : node_indices)
      cell_rhs[i] += 1.0;

    std::vector<double> cell_matrix_cont(num_nodes * num_nodes, 0.0);
    size_t n = 0;
    for (const size_t i : node_indices)
      for (const size_t j : node_indices)
        cell_matrix_cont[n++] = cell_matrix[i][j];

    MatSetValues(A,
                 static_cast<PetscInt>(num_nodes), dof_global_indices.data(),
                 static_cast<PetscInt>(num_nodes), dof_global_indices.data(),
                 cell_matrix_cont.data(), ADD_VALUES);

    VecSetValues(b,
                 static_cast<PetscInt>(num_nodes), dof_global_indices.data(),
                 cell_rhs.data(), ADD_VALUES);

    VecSetValues(x,
                 static_cast<PetscInt>(num_nodes), dof_global_indices.data(),
                 cell_rhs.data(), ADD_VALUES);
  }

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  VecAssemblyBegin(x);
  VecAssemblyEnd(x);

  MatInfo info;
  MatGetInfo(A,MAT_GLOBAL_SUM,&info);

  chi_log.Log(LOG_0) << "Number of mallocs used = " << info.mallocs
                       << "\nNumber of non-zeros allocated = "
                       << info.nz_allocated
                       << "\nNumber of non-zeros used = "
                       << info.nz_used
                       << "\nNumber of unneeded non-zeros = "
                       << info.nz_unneeded;


  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log() << "Done initializing diffusion solver.";
  MPI_Barrier(MPI_COMM_WORLD);
}