#ifndef CHI_MPI_MAP_ALL2ALL_H
#define CHI_MPI_MAP_ALL2ALL_H

#include<map>
#include<vector>

#include "chi_mpi.h"

namespace ChiMPI2
{

template<class T> std::map<uint64_t, std::vector<T>>
  MapAllToAll(const std::map<uint64_t, std::vector<T>>& pid_data_pairs,
              const MPI_Datatype data_mpi_type)
{
  auto& chi_mpi = ChiMPI::GetInstance();

  //============================================= Make sendcounts and
  //                                              senddispls
  std::vector<int> sendcounts(chi_mpi.process_count, 0);
  std::vector<int> senddispls(chi_mpi.process_count, 0);
  {
    size_t accumulated_displ = 0;
    for (const auto& pid_data_pair : pid_data_pairs)
    {
      const uint64_t pid = pid_data_pair.first;
      const auto&  data  = pid_data_pair.second;

      sendcounts[pid] = static_cast<int>(data.size());
      senddispls[pid] = static_cast<int>(accumulated_displ);
      accumulated_displ += data.size();
    }
  }

  //============================================= Communicate sendcounts to
  //                                              get recvcounts
  std::vector<int> recvcounts(chi_mpi.process_count, 0);

  MPI_Alltoall(sendcounts.data(), //sendbuf
               1, MPI_INT,        //sendcount, sendtype
               recvcounts.data(), //recvbuf
               1, MPI_INT,        //recvcount, recvtype
               MPI_COMM_WORLD);   //communicator

  //============================================= Populate recvdispls,
  //                                              sender_pids_set, and
  //                                              total_recv_count
  // All three these quantities are constructed
  // from recvcounts.
  std::vector<int>   recvdispls(chi_mpi.process_count, 0);
  std::set<uint64_t> sender_pids_set; //set of neighbor-partitions sending data
  size_t total_recv_count;
  {
    int displacement=0;
    for (int pid=0; pid < chi_mpi.process_count; ++pid)
    {
      recvdispls[pid] = displacement;
      displacement += recvcounts[pid];

      if (recvcounts[pid] > 0)
        sender_pids_set.insert(static_cast<uint64_t>(pid));
    }//for pid
    total_recv_count = displacement;
  }

  //============================================= Make sendbuf
  // The data for each partition is now loaded
  // into a single buffer
  std::vector<T> sendbuf;
  for (const auto& pid_data_pair : pid_data_pairs)
    sendbuf.insert(sendbuf.end(),
                   pid_data_pair.second.begin(),
                   pid_data_pair.second.end());

  //============================================= Make recvbuf
  std::vector<T> recvbuf(total_recv_count);

  //============================================= Communicate serial data
  MPI_Alltoallv(sendbuf.data(),        //sendbuf
                sendcounts.data(),     //sendcounts
                senddispls.data(),     //senddispls
                data_mpi_type,         //sendtype
                recvbuf.data(),        //recvbuf
                recvcounts.data(),     //recvcounts
                recvdispls.data(),     //recvdispls
                data_mpi_type,         //recvtype
                MPI_COMM_WORLD);       //comm

  std::map<uint64_t, std::vector<T>> output_data;
  {
    for (uint64_t pid : sender_pids_set)
    {
      const int data_count = recvcounts.at(pid);
      const int data_displ = recvdispls.at(pid);

      auto& data = output_data[pid];
      data.resize(data_count);

      for (int i=0; i<data_count; ++i)
        data.at(i) = recvbuf.at(data_displ + i);
    }
  }

  return output_data;
}

}//namespace ChiMPI2

#endif//CHI_MPI_MAP_ALL2ALL_H
