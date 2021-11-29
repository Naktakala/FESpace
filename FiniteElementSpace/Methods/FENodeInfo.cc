#include "FENodeInfo.h"

using namespace chi_math::finite_element;

bool NodeInfo::operator==(const NodeInfo &other) const
{
  if (type != other.type) return false;

  if (vertex_id_info != other.vertex_id_info) return false;
  if (cell_id_info   != other.cell_id_info  ) return false;
  if (aux_id_info    != other.aux_id_info   ) return false;

  return true;
}

chi_data_types::ByteArray NodeInfo::Serialize() const
{
  chi_data_types::ByteArray raw;

  auto WriteIDInfo = [&raw](const std::set<uint64_t>& id_set)
  {
    raw.Write<size_t>(id_set.size());
    for (uint64_t id_value : id_set)
      raw.Write<uint64_t>(id_value);
  };

  raw.Write<NodeType>(type);

  WriteIDInfo(vertex_id_info);
  WriteIDInfo(cell_id_info);
  WriteIDInfo(aux_id_info);

  raw.Write<chi_mesh::Vector3>(location);

  return raw;
}

NodeInfo NodeInfo::DeSerialize(const chi_data_types::ByteArray &raw,
                               size_t &address)
{
  auto ReadIDInfo = [&raw, &address]()
  {
    std::set<uint64_t> id_set;
    const auto num_id_values = raw.Read<size_t>(address, &address);
    for (size_t i=0; i<num_id_values; ++i)
      id_set.insert(raw.Read<uint64_t>(address, &address));

    return id_set;
  };

  const auto type = raw.Read<NodeType>(address,&address);

  std::vector<std::set<uint64_t>> id_info =
    {ReadIDInfo(),  //vertex_id_info
     ReadIDInfo(),  //cell_id_info
     ReadIDInfo()}; //aux_id_info

  auto location = raw.Read<chi_mesh::Vector3>(address, &address);

  return NodeInfo(type, id_info, location);
}

std::string NodeInfo::ToString() const
{
  std::stringstream outstr;

  outstr << "type: " << static_cast<int>(type) << "\n";
  outstr << "id_info: {";

  auto AttachIDInfo = [&outstr](const std::set<uint64_t>& id_set)
  {
    std::vector<uint64_t> vec_set;
    vec_set.reserve(id_set.size());
    for (uint64_t value : id_set)
      vec_set.push_back(value);

    outstr << "{";
    for (uint64_t value : vec_set)
    {
      outstr << value;
      if (value != vec_set.back()) outstr << ",";
    }
    outstr << "}";
  };

  AttachIDInfo(vertex_id_info); outstr << ",";
  AttachIDInfo(cell_id_info);   outstr << ",";
  AttachIDInfo(aux_id_info);


  outstr << "}\n";

  outstr << "location: " << location.PrintStr();

  return outstr.str();
}
