#include "fe_space.h"

#include <algorithm>

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

#include "ChiMPI/chi_mpi_map_all2all.h"

//###################################################################
/**Assembles the nodes for this discretization.*/
void chi_math::finite_element::SpatialDiscretization::
  AssembleNodes(const std::vector<NodeInfo>& node_register_node_info,
                const std::vector<std::pair<NodeInfo,uint64_t>>& ghost_node_register_node_info)
{
  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_ALL) << "Assembling nodes";
  MPI_Barrier(MPI_COMM_WORLD);

  /**Document this.*/
  class NodeFinder
  {
  private:
    const std::vector<NodeInfo>& m_node_list;

    std::map<uint64_t, std::set<size_t>> m_vertex_subscriptions;
    std::map<uint64_t, std::set<size_t>> m_cell_subscriptions;
    std::map<uint64_t, std::set<size_t>> m_aux_number_subscriptions;

  public:
    explicit NodeFinder(const std::vector<NodeInfo>& node_list) :
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

    std::pair<bool, size_t> FindNode(const NodeInfo& node)
    {
      for (uint64_t id : node.vertex_id_info)
        for (size_t node_id : m_vertex_subscriptions[id])
          if (node == m_node_list[node_id])
            return std::make_pair(true, node_id);

      for (uint64_t id : node.cell_id_info)
        for (size_t node_id : m_cell_subscriptions[id])
          if (node == m_node_list[node_id])
            return std::make_pair(true, node_id);

      for (uint64_t id : node.aux_id_info)
        for (size_t node_id : m_aux_number_subscriptions[id])
          if (node == m_node_list[node_id])
            return std::make_pair(true, node_id);

//      size_t node_id = 0;
//      for (const auto& list_node : m_node_list)
//      {
//        if (node == list_node)
//            return std::make_pair(true, node_id);
//        ++node_id;
//      }

      return std::make_pair(false,m_node_list.size());
    }

  };

  /**Document this too.*/
  class DynamicNodeFinder
  {
  private:
    std::vector<NodeInfo> m_node_list;

    std::map<uint64_t, std::set<size_t>> m_vertex_subscriptions;
    std::map<uint64_t, std::set<size_t>> m_cell_subscriptions;
    std::map<uint64_t, std::set<size_t>> m_aux_number_subscriptions;

  public:
    void reserve(size_t reserve_size)
    {
      m_node_list.reserve(reserve_size);
    }

    void push_back(const NodeInfo& node)
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

    std::vector<NodeInfo>& data() {return m_node_list;}

  public:
    std::pair<bool, size_t> FindNode(const NodeInfo& node)
    {
      for (uint64_t id : node.vertex_id_info)
        for (size_t node_id : m_vertex_subscriptions[id])
          if (node == m_node_list[node_id])
            return std::make_pair(true, node_id);

      for (uint64_t id : node.cell_id_info)
        for (size_t node_id : m_cell_subscriptions[id])
          if (node == m_node_list[node_id])
            return std::make_pair(true, node_id);

      for (uint64_t id : node.aux_id_info)
        for (size_t node_id : m_aux_number_subscriptions[id])
          if (node == m_node_list[node_id])
            return std::make_pair(true, node_id);

//      size_t node_id = 0;
//      for (const auto& list_node : m_node_list)
//      {
//        if (node == list_node)
//            return std::make_pair(true, node_id);
//        ++node_id;
//      }

      return std::make_pair(false,m_node_list.size());
    }

  };

  const std::string fname = __FUNCTION__;
  const uint64_t home = static_cast<uint64_t>(chi_mpi.location_id);
  const size_t node_register_size = node_register_node_info.size();

  //============================================= Initialize node register
  // We will later fill this register with a
  // mapping to unique geometrical info
  m_node_register_2_local_unique_map.assign(node_register_size, 0);

  //============================================= Build node register and
  //                                              get locally scoped unique
  //                                              nodes
  // This list will not contain any duplicates.
  DynamicNodeFinder finder_local_scope_unique_nodes;
  {
    size_t node_register_id = 0;
    for (const auto& node : node_register_node_info)
    {
      const auto location_info   = finder_local_scope_unique_nodes.FindNode(node);
      const bool   list_has_node = location_info.first;
      const size_t list_position = location_info.second; //will be end if not found

      if (not list_has_node) finder_local_scope_unique_nodes.push_back(node);

      m_node_register_2_local_unique_map[node_register_id] = list_position;

      ++node_register_id;
    }

  }
  const std::vector<NodeInfo>& local_scope_unique_nodes =
    finder_local_scope_unique_nodes.data();

  //============================================= Make m_local_unique_node_locations
  m_local_unique_node_locations.reserve(local_scope_unique_nodes.size());
  {
    auto& location_list = m_local_unique_node_locations;
    for (const auto& node : local_scope_unique_nodes)
      location_list.push_back(node.location);
  }

//  MPI_Barrier(MPI_COMM_WORLD);
//  chi_log.Log(LOG_ALL) << "Number of local scope unique nodes: "
//                       << local_scope_unique_nodes.size();
//  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Build partition mapped
  //                                              ghost scope nodes
  // This list will not contain any duplicates.
  // If a duplicate is encountered, the node with
  // the lowest partition-id will prevail.
  std::vector<std::pair<NodeInfo,uint64_t>> ghost_scope_nodes_pids;
  {
    for (const auto& master_node_pid_pair : ghost_node_register_node_info)
    {
      const auto&   node = master_node_pid_pair.first;
      const uint64_t pid = master_node_pid_pair.second;

      bool node_already_there = false;
      for (auto& node_pid_pair : ghost_scope_nodes_pids)
        if (node_pid_pair.first == node)
        {
          node_already_there = true;
          if (pid < node_pid_pair.second)
            node_pid_pair = std::make_pair(node, pid);
        }

      if (not node_already_there)
        ghost_scope_nodes_pids.emplace_back(node, pid);
    }
  }

//  MPI_Barrier(MPI_COMM_WORLD);
//  chi_log.Log(LOG_ALL) << "Number of ghost scope nodes: "
//                       << ghost_scope_nodes_pids.size();
//  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Filter local scope and break
  //                                              ties to get true local-nodes
  //                                              and true ghost-nodes
  //
  std::vector<NodeInfo> true_local_nodes;
  std::vector<std::pair<NodeInfo, uint64_t>> true_ghost_nodes;
  {
    const size_t num_lsu_nodes = local_scope_unique_nodes.size();
    std::vector<bool> true_local_flags(num_lsu_nodes,true);

    size_t num_ghosted_lsu_nodes = 0;
    for (const auto& ghost_node_pid_pair : ghost_scope_nodes_pids)
    {
      const auto&         node = ghost_node_pid_pair.first;
      const uint64_t ghost_pid = ghost_node_pid_pair.second;

      const auto location_info   = finder_local_scope_unique_nodes.FindNode(node);
      const bool   list_has_node = location_info.first;
      const size_t list_position = location_info.second; //will be end if not found

      if (list_has_node and home > ghost_pid)
      {
        true_ghost_nodes.emplace_back(node, ghost_pid);
        true_local_flags[list_position] = false;
        ++num_ghosted_lsu_nodes;
      }
    }

    const size_t num_est_true_local_nodes = num_lsu_nodes - num_ghosted_lsu_nodes;
    true_local_nodes.reserve(num_est_true_local_nodes);

    size_t node_id = 0;
    for (const auto& node : local_scope_unique_nodes)
    {
      if (true_local_flags[node_id]) true_local_nodes.push_back(node);
      ++node_id;
    }
  }

//  MPI_Barrier(MPI_COMM_WORLD);
//  chi_log.Log(LOG_ALL) << "Number of true local nodes: "
//                       << true_local_nodes.size()
//                       << " and true ghost nodes "
//                       << true_ghost_nodes.size();
//  MPI_Barrier(MPI_COMM_WORLD);

  // At this point we have a list of the true local nodes and a list of the
  // true ghost-nodes and to which partition they belong.

  // Next we map local node global addresses.

  //============================================= Gather local node counts from
  //                                              all processes
  const size_t num_true_local_nodes = true_local_nodes.size();
  std::vector<size_t> process_local_node_count(chi_mpi.process_count, 0);

  MPI_Allgather(&num_true_local_nodes,           //sendbuf
                1,                               //sendcount
                MPI_UNSIGNED_LONG,               //send datatype
                process_local_node_count.data(), //recvbuf
                1,                               //recvcount
                MPI_UNSIGNED_LONG,               //recv datatype,
                MPI_COMM_WORLD);                 //communicator

  //============================================= Create true local to global
  //                                              node mapping
  std::vector<int64_t> true_local_node_id_to_global_id_map;
  {
    true_local_node_id_to_global_id_map.assign(num_true_local_nodes, 0);
    int64_t running_total_node_count = 0;
    for (uint64_t pid=0; pid<static_cast<uint64_t>(chi_mpi.process_count); ++pid)
    {
      if (pid == home)
      {
        for (size_t i=0; i<num_true_local_nodes; ++i)
          true_local_node_id_to_global_id_map[i] =
            static_cast<int64_t>(running_total_node_count + i);
        break;
      }

      running_total_node_count +=
        static_cast<int64_t>(process_local_node_count[pid]);
    }//for pid
  }

//  MPI_Barrier(MPI_COMM_WORLD);
//  chi_log.Log(LOG_ALL) << "Creating true_local_node_id_to_global_id_map";
//  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Initialize node-register
  //                                              local-id map
  NodeFinder finder_true_local_node(true_local_nodes);
  m_node_register_local_ids.assign(node_register_size, -1);
  {
    uint64_t node_register_id = 0;
    for (const auto& node : node_register_node_info)
    {
      auto location_info = finder_true_local_node.FindNode(node);
      const bool has_node   = location_info.first;
      const size_t local_id = location_info.second;

      if (has_node)
        m_node_register_local_ids[node_register_id] =
          static_cast<int64_t>(local_id);

      ++node_register_id;
    }
  }

  // Now that we have local- and global-addresses for all true
  // local nodes we can start the process of mapping the ghost nodes.

//  MPI_Barrier(MPI_COMM_WORLD);
//  chi_log.Log(LOG_ALL) << "Checkpoint A";
//  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Consolidate ghost nodes
  //                                              onto ghosted partitions
  // This helps us tell each of the neighboring
  // partitions which nodes we want global addresses
  // for
  std::map<uint64_t, std::vector<NodeInfo>> true_ghost_nodes_consolidated;
  {
    for (const auto& ghost_node_pid : true_ghost_nodes)
    {
      const auto& node = ghost_node_pid.first;
      uint64_t    pid  = ghost_node_pid.second;

      auto& node_list = true_ghost_nodes_consolidated[pid];

      node_list.push_back(node);
    }
  }

//  MPI_Barrier(MPI_COMM_WORLD);
//  chi_log.Log(LOG_ALL) << "Checkpoint B";
//  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Serialize consolidated ghost
  //                                              nodes
  // This serialized data will be communicated via
  // an all-to-all process to provide each process
  // with a list of pairs of info. Each pair will
  // a querying partition-id and a list of nodes for
  // which that partition needs global-ids.
  typedef std::map<uint64_t, std::vector<std::byte>> MapPIDSerialData;

  MapPIDSerialData true_ghost_nodes_serialized;
  {
    for (const auto& pid_ghost_list : true_ghost_nodes_consolidated)
    {
      const uint64_t                      pid = pid_ghost_list.first;
      const std::vector<NodeInfo>& ghost_list = pid_ghost_list.second;

      chi_data_types::ByteArray serialized_data;

      serialized_data.Write<size_t>(ghost_list.size());

      for (auto& ghost_node : ghost_list)
        serialized_data.Append(ghost_node.Serialize());

      true_ghost_nodes_serialized[pid] = serialized_data.Data();
    }
  }

//  MPI_Barrier(MPI_COMM_WORLD);
//  chi_log.Log(LOG_ALL) << "Checkpoint C";
//  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Communicate serialized data
  MapPIDSerialData query_pids_nodes_serialized =
    ChiMPI2::MapAllToAll(true_ghost_nodes_serialized, MPI_BYTE);

  //============================================= Convert/DeSerialize query data
  //                                              to usable data
  std::map<uint64_t, std::vector<NodeInfo>> query_pids_nodes;
  {
    for (auto& pid_serial_data : query_pids_nodes_serialized)
    {
      const uint64_t      pid = pid_serial_data.first;
      auto& serial_data = pid_serial_data.second;

      chi_data_types::ByteArray buffer(std::move(serial_data));

      size_t address=0;
      while (address < buffer.Size())
      {
        const size_t num_nodes = buffer.Read<size_t>(address, &address);

        auto& node_list = query_pids_nodes[pid];
        node_list.reserve(num_nodes);

        for (size_t n=0; n<num_nodes; ++n)
        {
          node_list.emplace_back(NodeInfo::DeSerialize(buffer, address));
        }
      }
    }
  }

//  MPI_Barrier(MPI_COMM_WORLD);
//  chi_log.Log(LOG_ALL) << "Checkpoint D";
//  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Map the query-nodes
  // All the query nodes should be true local nodes
  // unless an error occurred upstream of this.
  // We simply have to lookup each query node in
  // the local list, get its local-id and then map
  // its global id.
  std::map<uint64_t, std::vector<int64_t>> query_pids_nodes_mapped;
  {
    for (const auto& pid_node_list : query_pids_nodes)
    {
      const uint64_t    pid = pid_node_list.first;
      const auto& node_list = pid_node_list.second;

      auto& map_list = query_pids_nodes_mapped[pid];
      map_list.reserve(node_list.size());

      for (const auto& node : node_list)
      {
        auto location_info = finder_true_local_node.FindNode(node);
        const bool has_node   = location_info.first;
        const size_t location = location_info.second;

        if (has_node)
          map_list.push_back(true_local_node_id_to_global_id_map[location]);
      }
    }
  }

//  MPI_Barrier(MPI_COMM_WORLD);
//  chi_log.Log(LOG_ALL) << "Checkpoint E";
//  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Communicate back mapped nodes
  // Now that the query-nodes have been mapped we
  // have to send all this information back to the
  // query-partitions.
  std::map<uint64_t, std::vector<int64_t>> ghost_pids_node_mapping =
    ChiMPI2::MapAllToAll(query_pids_nodes_mapped, MPI_LONG_LONG_INT);

//  MPI_Barrier(MPI_COMM_WORLD);
//  chi_log.Log(LOG_ALL) << "Checkpoint FG";
//  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Make global mappings
  {
    m_node_register_global_ids.assign(node_register_size, -1);

    uint64_t node_register_id = 0;
    for (const auto& node : node_register_node_info)
    {
      //=============================== Search true_local_nodes
      {
        const auto location_info = finder_true_local_node.FindNode(node);
        const bool   list_has_node = location_info.first;
        const size_t list_position = location_info.second; //will be end if not found

        if (list_has_node)
        {
          m_node_register_global_ids[node_register_id] =
            true_local_node_id_to_global_id_map[list_position];
          goto next_node_register;
        }
      }

      //=============================== Search consolidated true_ghosts
      for (const auto& pid_node_list : true_ghost_nodes_consolidated)
      {
        const uint64_t    pid = pid_node_list.first;
        const auto& node_list = pid_node_list.second;
        const auto& node_map  = ghost_pids_node_mapping.at(pid);

        for (size_t n=0; n<node_list.size(); ++n)
          if (node_list[n] == node)
          {
            m_node_register_global_ids[node_register_id] = node_map[n];
            goto end_of_ghost_loop;
          }
      }//for each ghost list
      end_of_ghost_loop:;

      next_node_register: ++node_register_id;
    }//for node in node_register
  }

//  MPI_Barrier(MPI_COMM_WORLD);
//  chi_log.Log(LOG_ALL) << "Checkpoint H";
//  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Check global-id map
  size_t num_invalid_global_ids = 0;
  for (int64_t value : m_node_register_global_ids)
    if (value < 0) ++num_invalid_global_ids;

  if (num_invalid_global_ids > 0)
    throw std::logic_error(fname + ": A total of " +
                           std::to_string(num_invalid_global_ids) +
                           " nodes in the node register do not have "
                           "global ids.");

//  MPI_Barrier(MPI_COMM_WORLD);
//  chi_log.Log(LOG_ALL) << "Number of local nodes and cells: "
//                       << num_true_local_nodes << " "
//                       << m_grid->local_cells.size();
//  MPI_Barrier(MPI_COMM_WORLD);

}