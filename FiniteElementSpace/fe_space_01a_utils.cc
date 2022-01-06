#include "fe_space.h"

#include "chi_mpi.h"

#include "ChiMPI/chi_mpi_utils.h"

using namespace chi_math::finite_element;

//###################################################################
/**Filters a vector of pairs for duplicates, based on whether the
 * NodeInfo portion of the pair (first part of pair) is a duplicate. If a
 * duplicate is encountered the pair with the lowest pid (second part of pair)
 * will persist.
*/
VecNodeInfoIDPair SpatialDiscretization::
  FilterGhostNodeRegister(const VecNodeInfoIDPair &GNR)
{
  VecNodeInfoIDPair FGNR; // Filtered Ghost Node Register

  for (const auto& gnr_node_pid_pair : GNR)
  {
    const auto&   gnr_node = gnr_node_pid_pair.first;
    const uint64_t gnr_pid = gnr_node_pid_pair.second;

    bool node_already_there = false;
    for (auto& node_pid_pair : FGNR)
      if (node_pid_pair.first == gnr_node)
      {
        node_already_there = true;
        if (gnr_pid < node_pid_pair.second)
          node_pid_pair = std::make_pair(gnr_node, gnr_pid);
      }

    if (not node_already_there)
      FGNR.emplace_back(gnr_node, gnr_pid);
  }

  return FGNR;
}

//###################################################################
/**Make a new list of local nodes (`NodeInfo`) by not adding nodes from the LNR
 * having duplicates in the GNR. A node-pid pair in the GNR is considered a
 * duplicate if the NodeInfo-part of the node-pid pair matches a node in the
 * LNR, and if the pid-part of the same node-pid pair is less than the local
 * partition-id.*/
std::vector<NodeInfo> SpatialDiscretization::
  SubtractGNRFromLNR(const NodeInfoListManager &LNR_manager,
                     const VecNodeInfoIDPair &GNR)
{
  auto& chi_mpi = ChiMPI::GetInstance();

  const auto& LNR = LNR_manager.Data();

  std::vector<NodeInfo> TLNR;

  const size_t LNR_size = LNR.size();
  std::vector<bool> true_local_flags(LNR_size, true);

  size_t num_GNs_in_LNR = 0;
  for (const auto& gnr_node_pid_pair : GNR)
  {
    const auto&   gnr_node = gnr_node_pid_pair.first;
    const uint64_t gnr_pid = gnr_node_pid_pair.second;

    const auto   lnr_find_info     = LNR_manager.FindNode(gnr_node);
    const bool   lnr_list_has_node = lnr_find_info.first;
    const size_t lnr_list_position = lnr_find_info.second; //will be end if not found

    if (lnr_list_has_node and  gnr_pid < chi_mpi.location_id/*home*/)
    {
      true_local_flags[lnr_list_position] = false;
      ++num_GNs_in_LNR;
    }
  }

  //TLN = True Local Nodes
  const size_t estimated_num_TLN = LNR_size - num_GNs_in_LNR;
  TLNR.reserve(estimated_num_TLN);

  size_t node_id = 0;
  for (const auto& node : LNR)
  {
    if (true_local_flags[node_id]) TLNR.push_back(node);
    ++node_id;
  }

  TLNR.shrink_to_fit();
  return TLNR;
}

//###################################################################
/**Make a new list of ghost nodes (`NodeInfo-PID` pair) by adding nodes from
 * the GNR having duplicates in the LNR. A node-pid pair in the GNR is
 * considered a duplicate if the NodeInfo-part of the node-pid pair matches a
 * node in the LNR, and if the pid-part of the same node-pid pair is less than
 * the local partition-id.*/
VecNodeInfoIDPair SpatialDiscretization::
  IntersectGNRWithLNR(const NodeInfoListManager &LNR_manager,
                      const VecNodeInfoIDPair &GNR)
{
  auto& chi_mpi = ChiMPI::GetInstance();

  const auto& LNR = LNR_manager.Data();

  VecNodeInfoIDPair IGNR; // Intersected Ghost Node Register

  IGNR.reserve(GNR.size());

  for (const auto& gnr_node_pid_pair : GNR)
  {
    const auto&   gnr_node = gnr_node_pid_pair.first;
    const uint64_t gnr_pid = gnr_node_pid_pair.second;

    const auto lnr_find_info     = LNR_manager.FindNode(gnr_node);
    const bool lnr_list_has_node = lnr_find_info.first;

    if (lnr_list_has_node and gnr_pid < chi_mpi.location_id/*home*/)
      IGNR.push_back(gnr_node_pid_pair);
  }

  IGNR.shrink_to_fit();
  return IGNR;
}

//###################################################################
/**Turns a GNR, which is a vector of `NodeInfo-PID` pairs, into a map where the
 * keys are the unique pids contained in the GNR and values are vectors storing
 * all the nodes associated with the pid.*/
LocINodeInfoMap SpatialDiscretization::
  ConsolidateGNR(const VecNodeInfoIDPair &GNR)
{
  std::map<uint64_t, std::vector<NodeInfo>> consolidated_GNR;

  for (const auto& gnr_node_pid : GNR)
  {
    const auto& gnr_node = gnr_node_pid.first;
    uint64_t    gnr_pid  = gnr_node_pid.second;

    auto& node_list = consolidated_GNR[gnr_pid]; // Will make new list if needed

    node_list.push_back(gnr_node);
  }

  return consolidated_GNR;
}

//###################################################################
/**Communicates a consolidated GNR to its mapped partitions. And returns
 * a consolidated QNR with the keys being the originating partition-ids.
 * Since `NodeInfo` is not a primitive type we have to convert a
 * consolidated GNR to a serialized version before communicating it via an
 * AllToAll communication.*/
LocINodeInfoMap SpatialDiscretization::
  SerializeAndCommunicateGNR(const LocINodeInfoMap& consolidated_GNR)
{
  //============================================= Serialize consolidated ghost
  //                                              nodes
  // This serialized data will be communicated via
  // an all-to-all process to provide each process
  // with a list of pairs of info. Each pair will
  // a querying partition-id and a list of nodes for
  // which that partition needs global-ids.
  std::map<uint64_t, std::vector<std::byte>> consolidated_GNR_serialized;
  {
    for (const auto& pid_ghost_list : consolidated_GNR)
    {
      const uint64_t                      pid = pid_ghost_list.first;
      const std::vector<NodeInfo>& ghost_list = pid_ghost_list.second;

      chi_data_types::ByteArray serialized_data;

      serialized_data.Write<size_t>(ghost_list.size());

      for (auto& ghost_node : ghost_list)
        serialized_data.Append(ghost_node.Serialize());

      consolidated_GNR_serialized[pid] = serialized_data.Data();
    }
  }

  //============================================= Communicate serialized data
  // All processes participate and will populate
  // a Query Node Register (QNR).
  std::map<uint64_t, std::vector<std::byte>> QNR_serialized =
    chi_mpi_utils::MapAllToAll(consolidated_GNR_serialized, MPI_BYTE);

  //============================================= Convert/DeSerialize query data
  //                                              to usable data
  std::map<uint64_t, std::vector<NodeInfo>> QNR;
  {
    for (auto& pid_serial_data : QNR_serialized)
    {
      const uint64_t      pid = pid_serial_data.first;
      auto& serial_data = pid_serial_data.second;

      chi_data_types::ByteArray buffer(std::move(serial_data));

      size_t address=0;
      while (address < buffer.Size())
      {
        const size_t num_nodes = buffer.Read<size_t>(address, &address);

        auto& node_list = QNR[pid];
        node_list.reserve(num_nodes);

        for (size_t n=0; n<num_nodes; ++n)
          node_list.emplace_back(NodeInfo::DeSerialize(buffer, address));
      }//while no overflow
    }//for key,value in QNR_serialized
  }

  return QNR;
}

//###################################################################
/**Given a list of true local nodes, this routine stacks local nodes
 * from all partitions to determine global ids for all the local nodes.*/
std::vector<int64_t> SpatialDiscretization::
  CreateGlobalIDMapForLNR(const std::vector<NodeInfo> &TLNR,
                          int64_t& num_global_nodes)
{
  auto& chi_mpi = ChiMPI::GetInstance();

  //============================================= Gather local node counts from
  //                                              all processes
  const size_t num_TLN = TLNR.size();
  std::vector<size_t> locI_num_TLN(chi_mpi.process_count, 0);

  MPI_Allgather(&num_TLN,                        //sendbuf
                1,                               //sendcount
                MPI_UNSIGNED_LONG,               //send datatype
                locI_num_TLN.data(),             //recvbuf
                1,                               //recvcount
                MPI_UNSIGNED_LONG,               //recv datatype,
                MPI_COMM_WORLD);                 //communicator

  //============================================= Create true local to global
  //                                              node mapping
  std::vector<int64_t> TLN_id_to_global_id_map;
  {
    auto& mapping = TLN_id_to_global_id_map;
    mapping.assign(num_TLN, -1);

    int64_t running_total_NC = 0; //NC = Node Count
    for (uint64_t pid=0; pid<static_cast<uint64_t>(chi_mpi.process_count); ++pid)
    {
      if (pid == chi_mpi.location_id/*home*/)
      {
        for (size_t i=0; i < num_TLN; ++i)
          mapping[i] = static_cast<int64_t>(running_total_NC + i);
        break;
      }

      running_total_NC += static_cast<int64_t>(locI_num_TLN[pid]);
    }//for pid
  }

  //============================================= Compute number of global nodes
  num_global_nodes = 0;
  for (uint64_t pid=0; pid<static_cast<uint64_t>(chi_mpi.process_count); ++pid)
    num_global_nodes += static_cast<int64_t>(locI_num_TLN[pid]);

  return TLN_id_to_global_id_map;
}


//###################################################################
/**Make a new list of local nodes (`NodeInfo`) by not adding nodes from the LNR
 * having duplicates in the GNR. A node-pid pair in the GNR is considered a
 * duplicate if the NodeInfo-part of the node-pid pair matches a node in the
 * LNR, and if the pid-part of the same node-pid pair is less than the local
 * partition-id.*/
VecNodeInfoIDPair SpatialDiscretization::
  SubtractLNRFromGNR(const NodeListFindManager &LNR_manager,
                     const VecNodeInfoIDPair &GNR)
{
  auto& chi_mpi = ChiMPI::GetInstance();

  VecNodeInfoIDPair SGNR; // Subtracted GNR
  SGNR.reserve(GNR.size());

  for (const auto& gnr_node_info_id_pair : GNR)
  {
    const auto& gnr_node = gnr_node_info_id_pair.first;
    const auto& gnr_pid  = gnr_node_info_id_pair.second;

    const auto   lnr_find_info     = LNR_manager.FindNode(gnr_node);
    const bool   lnr_list_has_node = lnr_find_info.first;

    if (not lnr_list_has_node)
      SGNR.push_back(gnr_node_info_id_pair);
  }

  SGNR.shrink_to_fit();
  return SGNR;
}