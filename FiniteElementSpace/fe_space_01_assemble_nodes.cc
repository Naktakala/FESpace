#include "fe_space.h"

#include <algorithm>

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiMPI/chi_mpi_map_all2all.h"

//###################################################################
/**Assembles the nodes for this discretization.

## _

### Abbreviations
- NR = Node Register.
- LNR = Local Node Register.
- GNR = Ghost Node Register.
- LNLL = Local Node Location List.
- GNLL = Ghost Node Location List.
- TLNR = True Local Node Register.
- IGNR = Intersected Ghost Node Register.
- QNR = Query Node Register.
- IQNR = Intersected Query Node Register.
- FQNR = Filtered Query Node Register.
- pid. Partition-id.

### Gameplan
 The process here starts with a Local Node Register (LNR) and a
 Ghost Node Register (GNR). A LNR is created from only the nodes on the local
 cells and is just a list (`std::vector<NodeInfo>`). A GNR is created from the
 the nodes on the ghost cells, with each ghost cell belonging to a particular
 partition, which requires each node to be accompanied by its partition-id at
 the moment it gets registered.
 Cells register their nodes in isolation from other cells, meaning two cells
 can register two nodes that are technically the same. We cannot determine
 whether a node (in a register) is duplicated by comparing its `Vector3`
 location because floating point precision issues will incur undue reliability
 issues. Instead, duplicates are identified via integer identification
 quantities, i.e. vertex-ids, cell-ids, cell-face node counts, etc.
 Duplicate nodes that lie on a partition interface have their ownership
 determined by the partition-id with the lowest integer-value.

 The steps in this process are as follows:
    - We filter the LNR, creating a filtered LNR (FLNR)
    - We use the FLNR to extract the local node location list (LNLL)
    - We use the LNR and FLNR to define a mapping from an LNR-id to a
      FLNR-id, which is identical to a mapping between LNR-id and LNLL-id.

    - Next we filter the GNR, creating a filtered GNR (FGNR)
    - Similar to the LNR, we use the FGNR to extract a ghost node location
      list (GNLL)
    - We also then use the GNR and FGNR to define a mapping from an LNR-id to a
      GNLL-id.

    - Subtract the FGNR from the FLNR.
      We need to know the true locally owned nodes so that we can get some
      sane way to determine local- and global-ids for any of the nodes (not just
      the true local ones). We therefore create a true local node register (TLNR)
      which is formed from the FLNR by removing/subtracting duplicate nodes in
      the FGNR. The set of duplicate nodes are only subtracted, individually, if
      the ghost node's partition-id is lower than the local partition-id
      (meaning the ghost's partition-id owns that node).

      NOTE: We can only determine the ownership of a node if it is in the LNR.
            This is because we can only break the partition-based ties for nodes
            that lie on the partition interfaces. This is in turn caused by the
            fact that we only have information on the first layer of ghost cells.

    - Populate local node register local-ids. With the TLNR established we
      can map nodes in the LNR to a possible node in the TLNR and therefore
      a local-id.
    - Assign global-ids to all true local nodes.
      Using the TLNR, we can start from partition-0 and start stacking
      the true local nodes for each consecutive location. Partition-0 with true
      local node local-ids  0, 1, ..., N_L0 will map to global-ids
      0, 1, ..., N_L0. Here N_L0 is the number of true local nodes -1 for
      partition-0. In general N_LX is the number of true local nodes -1 for
      partition-X.
      After partition-0 each partition-J will have the global-id of
      its first local node map to the global-id + 1 of the last local node of
      partition-(J-1), and thereafter the global-ids continue sequentially for
      the rest of the local nodes of partition-J. This can practically be done
      by simply all-gathering the num_true_local_nodes (N_L0, N_L1, etc.) of
      each partition. We have a routine that encapsulates this process called
      `CreateGlobalIDMapForLNR`.

    - Create an intersected GNR (IGNR).
      Intersect the set of nodes in the FGNR with the set of nodes in the FLNR.
      This will give us a list of nodes, on the partition-interfaces, which are
      not true local nodes, however, we know exactly to which partition-ids
      they are True Local to.
    - Create consolidated GNR for the IGNR (consolidated IGNR).
      For more efficient communication we consolidate nodes belonging to the
      same partition-ids into a map structure. This essentially allows us to
      easily use an AllToAllv MPI command.
    - Cross-communicate to establish intersected query node register (IQNR).
      By cross-communicating the IGNRs via an AllToAllv command, each location
      then has a map structure of query nodes (with the map-keys indicating
      which partition the query-nodes originate from).
    - Map the intersected query nodes.
      Since all the query nodes are guaranteed to be true local nodes, and since
      we know the global-ids of all the true local nodes, we simply
      search the TLNR for these nodes and provide the global-id.
    - Communicate back the IQNR node mappings. This returns as a consolidated
      IGNR mapping.

    - Make global-ids for all nodes in the LNR.
      At this point a node in the LNR will either be in the TLNR or the IGNR and
      therefore we ought to have a global-id for it.

    - Since all the nodes in the LNR are mapped we can now ask neighboring
      partitions to map nodes that are in the GNR even if they are not in the
      LNR. We start by consolidating the FGNR.
    - Cross-communicate to establish filtered query node register (FQNR).
    - We now look through the LNR and find a global id for each query node.
    - Communicate back the FQNR node mappings. This returns as a consolidated
      FGNR mapping.

    - Finally we assign the global-ids for each node in the GNR.

 */
void chi_math::finite_element::SpatialDiscretization::
  AssembleNodes(const std::vector<NodeInfo>& LNR, //Local Node Register
                const VecNodeInfoIDPair& GNR)     //Ghost Node Register
{
  const std::string fname = __FUNCTION__;
  typedef std::vector<NodeInfo> VecNodeInfo;

  chi_log.Log(LOG_0) << "Assembling nodes";
  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Filter the local node register
  NodeInfoListManager FLNR_manager;
  {
    for (const auto& node : LNR)
    {
      const auto location_info   = FLNR_manager.FindNode(node);
      const bool   list_has_node = location_info.first;

      if (not list_has_node) FLNR_manager.PushBack(node);
    }
  }
  // FLNR = Filtered Local Node Register
  const VecNodeInfo& FLNR = FLNR_manager.Data();

  //============================================= Make Local Node Location List
  // LNLL = Local Node Location List
  m_LNLL.reserve(FLNR.size());
  for (const auto& node : FLNR)
    m_LNLL.push_back(node.location);

  //============================================= Make LNR_id_2_LNLL_id_map
  m_LNR_id_2_LNLL_id_map.assign(LNR.size(), 0);
  {
    size_t node_register_id = 0;
    for (const auto& node : LNR)
    {
      const auto location_info   = FLNR_manager.FindNode(node);
      const bool   list_has_node = location_info.first;
      const size_t list_position = location_info.second; //will be end if not found

      if (list_has_node)
        m_LNR_id_2_LNLL_id_map[node_register_id] = list_position;
      else
        throw std::logic_error(fname + ": Error in LNR_id_2_LNLL_id_map");

      ++node_register_id;
    }
  }

  //============================================= Filter ghost node register
  // Eliminates duplicates. If a duplicate is
  // encountered, the node with the lowest
  // partition-id will prevail.
  //
  // FGNR = Filtered Ghost Node Register
  VecNodeInfoIDPair FGNR = FilterGhostNodeRegister(GNR);

  //============================================= Make Ghost Node Location List
  // GNLL = Ghost Node Location List
  m_GNLL.reserve(FGNR.size());
  for (const auto& node_id_info_pair: FGNR)
    m_GNLL.push_back(node_id_info_pair.first.location);

  //============================================= Make GNR_id_2_GNLL_id_map
  m_GNR_id_2_GNLL_id_map.assign(GNR.size(), 0);
  {
    const size_t FGNR_size = FGNR.size();
    size_t node_register_id = 0;
    for (const auto& node_info_id_pair : GNR)
    {
      const auto& registery_node = node_info_id_pair.first;

      // frid = filtered_register_id
      for (size_t frid=0; frid<FGNR_size; ++frid)
      {
        const auto& filtered_registery_node = FGNR[frid].first;

        if (registery_node == filtered_registery_node)
          m_GNR_id_2_GNLL_id_map[node_register_id] = frid;
      }
      ++node_register_id;
    }
  }

  //============================================= Subtract the filtered_GNR set
  //                                              from the filtered_LNR to
  //                                              determine true local nodes
  // TLNR = True Local Node Register
  VecNodeInfo TLNR = SubtractGNRFromLNR(FLNR_manager, FGNR);

  //============================================= Populate local node-register
  //                                              local-ids
  NodeListFindManager TLNR_finder(TLNR);
  m_LNR_local_ids.assign(LNR.size(), -1);
  {
    uint64_t node_register_id = 0;
    for (const auto& node : LNR)
    {
      auto location_info = TLNR_finder.FindNode(node);
      const bool has_node   = location_info.first;
      const size_t local_id = location_info.second;

      if (has_node)
        m_LNR_local_ids[node_register_id] = static_cast<int64_t>(local_id);

      ++node_register_id;
    }
  }

  //============================================= Assign global-ids to all
  //                                              true local nodes.
  std::vector<int64_t> TLNR_id_to_global_id_map = CreateGlobalIDMapForLNR(TLNR);

  //============================================= Intersect FLNR and FGNR
  VecNodeInfoIDPair IGNR = IntersectGNRWithLNR(FLNR_manager, FGNR);

  //============================================= Consolidate intersected ghost
  //                                              nodes onto ghosted partitions
  // This helps us tell each of the neighboring
  // partitions which nodes we want global
  // addresses for
  std::map<uint64_t, VecNodeInfo> consolidated_IGNR = ConsolidateGNR(IGNR);

  //============================================= Cross-communicate intersected
  //                                              ghost nodes
  LocINodeInfoMap IQNR = SerializeAndCommunicateGNR(consolidated_IGNR);

  //============================================= Map the intersected query-nodes
  // All the query nodes should be true local nodes
  // unless an error occurred upstream of this.
  // We simply have to lookup each query node in
  // the local list, get its local-id and then map
  // its global id.
  std::map<uint64_t, std::vector<int64_t>> IQNR_mapped;
  {
    for (const auto& pid_node_list : IQNR)
    {
      const uint64_t    pid = pid_node_list.first;
      const auto& node_list = pid_node_list.second;

      auto& map_list = IQNR_mapped[pid];
      map_list.reserve(node_list.size());

      for (const auto& node : node_list)
      {
        auto location_info = TLNR_finder.FindNode(node);
        const bool has_node   = location_info.first;
        const size_t location = location_info.second;

        if (has_node)
          map_list.push_back(TLNR_id_to_global_id_map[location]);
      }
    }
  }

  //============================================= Communicate back IGNR nodes
  // Now that the query-nodes have been mapped we
  // have to send all this information back to the
  // query-partitions.
  std::map<uint64_t, std::vector<int64_t>> consolidated_IGNR_mapping =
    ChiMPI2::MapAllToAll(IQNR_mapped, MPI_LONG_LONG_INT);

  //============================================= Make global mappings for LNR
  {
    m_LNR_global_ids.assign(LNR.size(), -1);

    uint64_t node_register_id = 0;
    for (const auto& node : LNR)
    {
      //=============================== Search True Local Nodes
      {
        const auto location_info = TLNR_finder.FindNode(node);
        const bool   list_has_node = location_info.first;
        const size_t list_position = location_info.second; //will be end if not found

        if (list_has_node)
        {
          m_LNR_global_ids[node_register_id] =
            TLNR_id_to_global_id_map[list_position];
          goto next_local_node;
        }
      }

      //=============================== Search consolidated Ghost Node Register
      for (const auto& pid_node_list : consolidated_IGNR)
      {
        const uint64_t    pid = pid_node_list.first;
        const auto& node_list = pid_node_list.second;
        const auto& node_map  = consolidated_IGNR_mapping.at(pid);

        for (size_t n=0; n<node_list.size(); ++n)
          if (node_list[n] == node)
          {
            m_LNR_global_ids[node_register_id] = node_map[n];
            goto next_local_node;
          }
      }//for each ghost list

      next_local_node: ++node_register_id;
    }//for node in LNR
  }

  //============================================= Consolidate filtered ghost
  //                                              nodes onto ghosted partitions
  // This helps us tell each of the neighboring
  // partitions which nodes we want global
  // addresses for
  std::map<uint64_t, VecNodeInfo> consolidated_FGNR = ConsolidateGNR(FGNR);

  //============================================= Cross-communicate filtered
  //                                              ghost nodes
  LocINodeInfoMap FQNR = SerializeAndCommunicateGNR(consolidated_FGNR);

  //============================================= Map the filtered query-nodes
  // The query nodes here will be either true local
  // nodes or intersected ghost nodes. We need to
  // search both these list to find the mapping
  std::map<uint64_t, std::vector<int64_t>> FQNR_mapped;
  {
    NodeListFindManager LNR_finder(LNR);

    for (const auto& fqnr_pid_node_list : FQNR)
    {
      const uint64_t    fqnr_pid = fqnr_pid_node_list.first;
      const auto& fqnr_node_list = fqnr_pid_node_list.second;

      auto& fqnr_map_list = FQNR_mapped[fqnr_pid];
      fqnr_map_list.reserve(fqnr_node_list.size());

      for (const auto& fqnr_node : fqnr_node_list)
      {
        auto location_info = LNR_finder.FindNode(fqnr_node);
        const bool has_node   = location_info.first;
        const size_t LNR_id = location_info.second;

        if (has_node)
          fqnr_map_list.push_back(m_LNR_global_ids[LNR_id]);
        else
          throw std::logic_error(fname + ": Query node mapping failure.");
      }
    }
  }

  //============================================= Communicate back mapped nodes
  // Now that the query-nodes have been mapped we
  // have to send all this information back to the
  // query-partitions.
  std::map<uint64_t, std::vector<int64_t>> consolidated_FGNR_mapping =
    ChiMPI2::MapAllToAll(FQNR_mapped, MPI_LONG_LONG_INT);

  //============================================= Make global mappings for GNR
  {
    m_GNR_global_ids.assign(GNR.size(), -1);

    uint64_t node_register_id = 0;
    for (const auto& node_info_id_pair : GNR)
    {
      const auto& node = node_info_id_pair.first;

      //=============================== Search consolidated Ghost Node Register
      for (const auto& pid_node_list : consolidated_FGNR)
      {
        const uint64_t    pid = pid_node_list.first;
        const auto& node_list = pid_node_list.second;
        const auto& node_map  = consolidated_FGNR_mapping.at(pid);

        for (size_t n=0; n<node_list.size(); ++n)
          if (node_list[n] == node)
          {
            m_GNR_global_ids[node_register_id] = node_map[n];
            goto next_ghost_node;
          }
      }//for each ghost list

      next_ghost_node: ++node_register_id;
    }//for node-id pair in GNR
  }

  //============================================= Check global-id map
  size_t num_invalid_global_ids = 0;
  for (int64_t value : m_LNR_global_ids)
    if (value < 0) ++num_invalid_global_ids;

  for (int64_t value : m_GNR_global_ids)
    if (value < 0) ++num_invalid_global_ids;

  if (num_invalid_global_ids > 0)
    throw std::logic_error(fname + ": A total of " +
                           std::to_string(num_invalid_global_ids) +
                           " nodes in the node register do not have "
                           "global ids.");

}