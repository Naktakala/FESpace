#include "chi_mesh_utils.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

chi_mesh::CellVolumeAndFaceAreas chi_mesh::
  ComputeCellVolumeAndFaceAreas(const chi_mesh::Cell &cell,
                                  const chi_mesh::MeshContinuum& grid)
{
  CellVolumeAndFaceAreas cell_info;

  if (cell.Type() == chi_mesh::CellType::SLAB)
  {
    const auto& v0 = grid.vertices[cell.vertex_ids[0]];
    const auto& v1 = grid.vertices[cell.vertex_ids[1]];

    cell_info.volume = (v1-v0).Norm();
    cell_info.face_areas.reserve(2);
    cell_info.face_areas.push_back(1.0); //First face  unit area
    cell_info.face_areas.push_back(1.0); //Second face unit area
  }
  else if (cell.Type() == chi_mesh::CellType::POLYGON)
  {
    const size_t num_faces = cell.faces.size();

    double volume = 0.0;
    std::vector<double> face_areas;
    face_areas.reserve(num_faces);
    for (const auto& face_edge : cell.faces)
    {
      const auto& v0 = grid.vertices[face_edge.vertex_ids[0]];
      const auto& v1 = grid.vertices[face_edge.vertex_ids[1]];
      const auto& v2 = cell.centroid;

      chi_mesh::Vector3 s01 = v1-v0; //side 01
      chi_mesh::Vector3 s02 = v2-v0; //side 02

      volume += 0.5*(s01.Cross(s02).Norm())*1.0; //Area triangle times unit h

      face_areas.push_back((v1-v0).Norm()); //Length of an edge time unit h
    }//for face_edge

    cell_info.volume = volume;
    cell_info.face_areas = std::move(face_areas);
  }
  else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
  {
    const size_t num_faces = cell.faces.size();

    double volume = 0.0;
    std::vector<double> face_areas;
    face_areas.reserve(num_faces);
    for (const auto& face : cell.faces)
    {
      size_t num_f_verts = face.vertex_ids.size();
      double face_area = 0.0;

      //Loop over sub-facets and form tetrahedra
      for (size_t fv = 0; fv < num_f_verts; ++fv)
      {
        size_t fvplus1 = (fv < (num_f_verts-1))? fv+1 : 0;

        const auto& v0 = grid.vertices[face.vertex_ids[fv     ]];
        const auto& v1 = grid.vertices[face.vertex_ids[fvplus1]];
        const auto& v2 = face.centroid;
        const auto& v3 = cell.centroid;

        chi_mesh::Vector3 s01 = v1-v0; //side 01
        chi_mesh::Vector3 s02 = v2-v0; //side 02
        chi_mesh::Vector3 s03 = v3-v0; //side 02

        face_area += 0.5*(s01.Cross(s02).Norm()); //Area triangle

        chi_mesh::Matrix3x3 J;

        J.SetColJVec(0,s01);
        J.SetColJVec(1,s02);
        J.SetColJVec(2,s03);

        volume += J.Det()/6.0; //Volume of Tetrahedron
      }

      face_areas.push_back(face_area);
    }//for face

    cell_info.volume = volume;
    cell_info.face_areas = face_areas;
  }
  else
    throw std::logic_error(std::string(__FUNCTION__) + "Unsupported cell-type"
                           " encountered.");

  return cell_info;
}