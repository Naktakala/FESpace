#ifndef FESPACETEST_FENODEINFO_H
#define FESPACETEST_FENODEINFO_H

#include <vector>
#include <set>

#include "ChiMesh/chi_mesh.h"
#include "ChiDataTypes/byte_array.h"

namespace chi_math::finite_element
{

enum class NodeType
{
  CORNER   = 1, ///< On the corner of cell and numerous faces
  EDGE     = 2, ///< On the edge of numerous faces
  FACE     = 3, ///< On a face but not on any edges or corners
  INTERNAL = 4, ///< Completely within a cell
};
typedef std::vector<std::set<uint64_t>> IdentifyingInfo;
struct NodeInfo
{
  NodeType type;
  std::set<uint64_t> vertex_id_info;
  std::set<uint64_t> cell_id_info;
  std::set<uint64_t> aux_id_info;

  chi_mesh::Vector3 location;

  NodeInfo(NodeType in_type,
           IdentifyingInfo in_id_info,
           const chi_mesh::Vector3& in_location) :
  type(in_type),
  location(in_location)
  {
    const size_t num_id_info_items = in_id_info.size();
    if (num_id_info_items >= 1) vertex_id_info = in_id_info[0];
    if (num_id_info_items >= 2) cell_id_info = in_id_info[1];
    if (num_id_info_items >= 3) aux_id_info = in_id_info[2];
  }

  bool operator==(const NodeInfo& other) const;

  chi_data_types::ByteArray Serialize() const;
  static NodeInfo DeSerialize(const chi_data_types::ByteArray& raw,
                              size_t& address);
  std::string ToString() const;
};

}// namespace chi_math::finite_element


#endif //FESPACETEST_FENODEINFO_H
