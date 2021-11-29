#include "diffusion_solver.h"

DiffusionSolver::DiffusionSolver() :
  chi_physics::Solver("DiffusionSolver")
{
  basic_options.AddOption("SDM", std::string("PWLC"));
}