#ifndef CHI_MESH_UTILS_H
#define CHI_MESH_UTILS_H

#include "ChiMesh/Cell/cell.h"
#include <vector>

namespace chi_mesh
{
  struct CellVolumeAndFaceAreas
  {
    double volume;
    std::vector<double> face_areas;
  };

  CellVolumeAndFaceAreas
    ComputeCellVolumeAndFaceAreas(const chi_mesh::Cell& cell,
                                  const chi_mesh::MeshContinuum& grid);
}//namespace chi_mesh

#endif //CHI_MESH_UTILS_H
