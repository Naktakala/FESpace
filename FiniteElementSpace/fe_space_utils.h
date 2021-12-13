#ifndef FESPACETEST_FE_SPACE_UTILS_H
#define FESPACETEST_FE_SPACE_UTILS_H

#include "Methods/FENodeInfo.h"

namespace chi_math::finite_element
{

  /**Document this.*/
  class NodeListFindManager
  {
  private:
    const std::vector<NodeInfo>& m_node_list;

    std::map<uint64_t, std::set<size_t>> m_vertex_subscriptions;
    std::map<uint64_t, std::set<size_t>> m_cell_subscriptions;
    std::map<uint64_t, std::set<size_t>> m_aux_number_subscriptions;

  public:
    explicit NodeListFindManager(const std::vector<NodeInfo>& node_list) :
      m_node_list(node_list)
    {
      size_t node_id = 0;
      for (const auto& node : node_list)
      {
        for (uint64_t id : node.vertex_id_info)
          m_vertex_subscriptions[id].insert(node_id);
        for (uint64_t id : node.cell_id_info)
          m_cell_subscriptions[id].insert(node_id);
        for (uint64_t id : node.aux_id_info)
          m_aux_number_subscriptions[id].insert(node_id);
        ++node_id;
      }
    }

    std::pair<bool, size_t> FindNode(const NodeInfo& node) const
    {
      for (uint64_t id : node.vertex_id_info)
        if (m_vertex_subscriptions.count(id) > 0)
          for (size_t node_id : m_vertex_subscriptions.at(id))
            if (node == m_node_list[node_id])
              return std::make_pair(true, node_id);

      for (uint64_t id : node.cell_id_info)
        if (m_cell_subscriptions.count(id) > 0)
          for (size_t node_id : m_cell_subscriptions.at(id))
            if (node == m_node_list[node_id])
              return std::make_pair(true, node_id);

      for (uint64_t id : node.aux_id_info)
        if (m_aux_number_subscriptions.count(id) > 0)
          for (size_t node_id : m_aux_number_subscriptions.at(id))
            if (node == m_node_list[node_id])
              return std::make_pair(true, node_id);

      return std::make_pair(false,m_node_list.size());
    }

  };

  /**Document this too.*/
  class NodeInfoListManager
  {
  private:
    std::vector<NodeInfo> m_node_list{};

    std::map<uint64_t, std::set<size_t>> m_vertex_subscriptions;
    std::map<uint64_t, std::set<size_t>> m_cell_subscriptions;
    std::map<uint64_t, std::set<size_t>> m_aux_number_subscriptions;

  public:
    void Reserve(size_t reserve_size)
    {
      m_node_list.reserve(reserve_size);
    }

    void PushBack(const NodeInfo& node)
    {
      const size_t node_id = m_node_list.size();
      m_node_list.push_back(node);

      for (uint64_t id : node.vertex_id_info)
        m_vertex_subscriptions[id].insert(node_id);

      for (uint64_t id : node.cell_id_info)
        m_cell_subscriptions[id].insert(node_id);

      for (uint64_t id : node.aux_id_info)
        m_aux_number_subscriptions[id].insert(node_id);
    }

    std::vector<NodeInfo>& Data() {return m_node_list;}
    const std::vector<NodeInfo>& Data() const {return m_node_list;}

  public:
    std::pair<bool, size_t> FindNode(const NodeInfo& node) const
    {
      for (uint64_t id : node.vertex_id_info)
        if (m_vertex_subscriptions.count(id) > 0)
          for (size_t node_id : m_vertex_subscriptions.at(id))
            if (node == m_node_list[node_id])
              return std::make_pair(true, node_id);

      for (uint64_t id : node.cell_id_info)
        if (m_cell_subscriptions.count(id) > 0)
          for (size_t node_id : m_cell_subscriptions.at(id))
            if (node == m_node_list[node_id])
              return std::make_pair(true, node_id);

      for (uint64_t id : node.aux_id_info)
        if (m_aux_number_subscriptions.count(id) > 0)
          for (size_t node_id : m_aux_number_subscriptions.at(id))
            if (node == m_node_list[node_id])
              return std::make_pair(true, node_id);

      return std::make_pair(false,m_node_list.size());
    }

  };

}//namespace chi_math::finite_element::SpatialDiscretization

#endif //FESPACETEST_FE_SPACE_UTILS_H
