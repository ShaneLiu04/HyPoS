#include "comm/halo_exchanger.hpp"
#include "utils/logger.hpp"
#include <cstring>

namespace hypo {

// ============================================================================
// PointToPointExchanger
// ============================================================================

void PointToPointExchanger::initialize(Subgrid& subgrid) {
    comm_ = subgrid.comm();

    Index faceSizeLR = subgrid.nyLocal() * subgrid.nzLocal() * subgrid.haloWidth();
    Index faceSizeUD = subgrid.nxLocal() * subgrid.nzLocal() * subgrid.haloWidth();
    Index faceSizeBF = subgrid.nxLocal() * subgrid.nyLocal() * subgrid.haloWidth();

    sendBufLeft_.resize(faceSizeLR);
    sendBufRight_.resize(faceSizeLR);
    recvBufLeft_.resize(faceSizeLR);
    recvBufRight_.resize(faceSizeLR);

    sendBufDown_.resize(faceSizeUD);
    sendBufUp_.resize(faceSizeUD);
    recvBufDown_.resize(faceSizeUD);
    recvBufUp_.resize(faceSizeUD);

    requests_.resize(8, MPI_REQUEST_NULL);
}

void PointToPointExchanger::packSend(Subgrid& subgrid, int direction, std::vector<Real>& buf) const {
    Real* u = subgrid.u().data();
    const Index nxT = subgrid.nxTotal();
    const Index nyT = subgrid.nyTotal();
    const Index hw = subgrid.haloWidth();
    Index idx = 0;

    if (direction == 0) { // Left: pack rightmost interior face (i = nxLocal)
        for (Index k = 0; k < subgrid.nzTotal(); ++k) {
            for (Index j = 0; j < subgrid.nyTotal(); ++j) {
                for (Index h = 0; h < static_cast<Index>(hw); ++h) {
                    buf[idx++] = u[k * nyT * nxT + j * nxT + (subgrid.iEnd() - hw + h)];
                }
            }
        }
    } else if (direction == 1) { // Right: pack leftmost interior face (i = halo)
        for (Index k = 0; k < subgrid.nzTotal(); ++k) {
            for (Index j = 0; j < subgrid.nyTotal(); ++j) {
                for (Index h = 0; h < static_cast<Index>(hw); ++h) {
                    buf[idx++] = u[k * nyT * nxT + j * nxT + (subgrid.iBegin() + h)];
                }
            }
        }
    } else if (direction == 2) { // Down: pack top interior face (j = nyLocal)
        for (Index k = 0; k < subgrid.nzTotal(); ++k) {
            for (Index h = 0; h < static_cast<Index>(hw); ++h) {
                for (Index i = 0; i < subgrid.nxTotal(); ++i) {
                    buf[idx++] = u[k * nyT * nxT + (subgrid.jEnd() - hw + h) * nxT + i];
                }
            }
        }
    } else if (direction == 3) { // Up: pack bottom interior face (j = halo)
        for (Index k = 0; k < subgrid.nzTotal(); ++k) {
            for (Index h = 0; h < static_cast<Index>(hw); ++h) {
                for (Index i = 0; i < subgrid.nxTotal(); ++i) {
                    buf[idx++] = u[k * nyT * nxT + (subgrid.jBegin() + h) * nxT + i];
                }
            }
        }
    }
}

void PointToPointExchanger::unpackRecv(Subgrid& subgrid, int direction, const std::vector<Real>& buf) const {
    Real* u = subgrid.u().data();
    const Index nxT = subgrid.nxTotal();
    const Index nyT = subgrid.nyTotal();
    const Index hw = subgrid.haloWidth();
    Index idx = 0;

    if (direction == 0) { // Left halo (i = 0..halo-1)
        for (Index k = 0; k < subgrid.nzTotal(); ++k) {
            for (Index j = 0; j < subgrid.nyTotal(); ++j) {
                for (Index h = 0; h < static_cast<Index>(hw); ++h) {
                    u[k * nyT * nxT + j * nxT + h] = buf[idx++];
                }
            }
        }
    } else if (direction == 1) { // Right halo (i = nxLocal+halo..nxTotal-1)
        for (Index k = 0; k < subgrid.nzTotal(); ++k) {
            for (Index j = 0; j < subgrid.nyTotal(); ++j) {
                for (Index h = 0; h < static_cast<Index>(hw); ++h) {
                    u[k * nyT * nxT + j * nxT + (subgrid.iEnd() + h)] = buf[idx++];
                }
            }
        }
    } else if (direction == 2) { // Bottom halo (j = 0..halo-1)
        for (Index k = 0; k < subgrid.nzTotal(); ++k) {
            for (Index h = 0; h < static_cast<Index>(hw); ++h) {
                for (Index i = 0; i < subgrid.nxTotal(); ++i) {
                    u[k * nyT * nxT + h * nxT + i] = buf[idx++];
                }
            }
        }
    } else if (direction == 3) { // Top halo (j = nyLocal+halo..nyTotal-1)
        for (Index k = 0; k < subgrid.nzTotal(); ++k) {
            for (Index h = 0; h < static_cast<Index>(hw); ++h) {
                for (Index i = 0; i < subgrid.nxTotal(); ++i) {
                    u[k * nyT * nxT + (subgrid.jEnd() + h) * nxT + i] = buf[idx++];
                }
            }
        }
    }
}

void PointToPointExchanger::exchange(Subgrid& subgrid) {
    beginExchange(subgrid);
    endExchange(subgrid);
}

void PointToPointExchanger::beginExchange(Subgrid& subgrid) {
    for (auto& req : requests_) {
        req = MPI_REQUEST_NULL;
    }
    int reqIdx = 0;

    // Pack and send left
    if (subgrid.neighborLeft() != MPI_PROC_NULL) {
        packSend(subgrid, 0, sendBufLeft_);
        MPI_Isend(sendBufLeft_.data(), static_cast<int>(sendBufLeft_.size()), MPI_DOUBLE,
                  subgrid.neighborLeft(), 0, comm_, &requests_[reqIdx]);
    }
    ++reqIdx;
    // Receive left
    if (subgrid.neighborLeft() != MPI_PROC_NULL) {
        MPI_Irecv(recvBufLeft_.data(), static_cast<int>(recvBufLeft_.size()), MPI_DOUBLE,
                  subgrid.neighborLeft(), 1, comm_, &requests_[reqIdx]);
    }
    ++reqIdx;

    // Pack and send right
    if (subgrid.neighborRight() != MPI_PROC_NULL) {
        packSend(subgrid, 1, sendBufRight_);
        MPI_Isend(sendBufRight_.data(), static_cast<int>(sendBufRight_.size()), MPI_DOUBLE,
                  subgrid.neighborRight(), 1, comm_, &requests_[reqIdx]);
    }
    ++reqIdx;
    // Receive right
    if (subgrid.neighborRight() != MPI_PROC_NULL) {
        MPI_Irecv(recvBufRight_.data(), static_cast<int>(recvBufRight_.size()), MPI_DOUBLE,
                  subgrid.neighborRight(), 0, comm_, &requests_[reqIdx]);
    }
    ++reqIdx;

    // Pack and send down
    if (subgrid.neighborDown() != MPI_PROC_NULL) {
        packSend(subgrid, 2, sendBufDown_);
        MPI_Isend(sendBufDown_.data(), static_cast<int>(sendBufDown_.size()), MPI_DOUBLE,
                  subgrid.neighborDown(), 2, comm_, &requests_[reqIdx]);
    }
    ++reqIdx;
    // Receive down
    if (subgrid.neighborDown() != MPI_PROC_NULL) {
        MPI_Irecv(recvBufDown_.data(), static_cast<int>(recvBufDown_.size()), MPI_DOUBLE,
                  subgrid.neighborDown(), 3, comm_, &requests_[reqIdx]);
    }
    ++reqIdx;

    // Pack and send up
    if (subgrid.neighborUp() != MPI_PROC_NULL) {
        packSend(subgrid, 3, sendBufUp_);
        MPI_Isend(sendBufUp_.data(), static_cast<int>(sendBufUp_.size()), MPI_DOUBLE,
                  subgrid.neighborUp(), 3, comm_, &requests_[reqIdx]);
    }
    ++reqIdx;
    // Receive up
    if (subgrid.neighborUp() != MPI_PROC_NULL) {
        MPI_Irecv(recvBufUp_.data(), static_cast<int>(recvBufUp_.size()), MPI_DOUBLE,
                  subgrid.neighborUp(), 2, comm_, &requests_[reqIdx]);
    }
    ++reqIdx;
}

void PointToPointExchanger::endExchange(Subgrid& subgrid) {
    MPI_Waitall(static_cast<int>(requests_.size()), requests_.data(), MPI_STATUSES_IGNORE);

    // Unpack received data
    if (subgrid.neighborLeft() != MPI_PROC_NULL) {
        unpackRecv(subgrid, 0, recvBufLeft_);
    }
    if (subgrid.neighborRight() != MPI_PROC_NULL) {
        unpackRecv(subgrid, 1, recvBufRight_);
    }
    if (subgrid.neighborDown() != MPI_PROC_NULL) {
        unpackRecv(subgrid, 2, recvBufDown_);
    }
    if (subgrid.neighborUp() != MPI_PROC_NULL) {
        unpackRecv(subgrid, 3, recvBufUp_);
    }
}


// ============================================================================
// CollectiveExchanger
// ============================================================================

void CollectiveExchanger::initialize(Subgrid& subgrid) {
    // For simplicity, use the subgrid's comm directly (assumed to be a cartesian comm)
    comm_ = subgrid.comm();

    // Each process sends 4 (2D) faces to its neighbors
    // We use a simple all-to-all approach: pack all 4 faces into one buffer,
    // then use MPI_Allgatherv or MPI_Neighbor_allgather
    // For simplicity and compatibility, we'll use a manual approach with MPI_Alltoallv
    // since MPI_Neighbor_allgather may not be available in all implementations.

    Index faceSizeLR = subgrid.nyLocal() * subgrid.nzLocal() * subgrid.haloWidth();
    Index faceSizeUD = subgrid.nxLocal() * subgrid.nzLocal() * subgrid.haloWidth();
    Index totalSend = 0;
    if (subgrid.neighborLeft() != MPI_PROC_NULL) totalSend += faceSizeLR;
    if (subgrid.neighborRight() != MPI_PROC_NULL) totalSend += faceSizeLR;
    if (subgrid.neighborDown() != MPI_PROC_NULL) totalSend += faceSizeUD;
    if (subgrid.neighborUp() != MPI_PROC_NULL) totalSend += faceSizeUD;

    sendBuf_.resize(totalSend);
    recvBuf_.resize(totalSend);
}

void CollectiveExchanger::exchange(Subgrid& subgrid) {
    // Fallback: use P2P for now since MPI_Neighbor_allgather requires a
    // dedicated cartesian comm with specific topology. In a real implementation,
    // we'd use MPI_Ineighbor_allgatherv with a proper cartesian communicator.
    // For this reference implementation, we delegate to P2P logic inline.

    // Pack all 4 faces into sendBuf_
    Real* u = subgrid.u().data();
    const Index nxT = subgrid.nxTotal();
    const Index nyT = subgrid.nyTotal();
    const Index hw = subgrid.haloWidth();
    Index offset = 0;

    auto pack_face = [&](auto index_fn, Index count) {
        for (Index i = 0; i < count; ++i) {
            sendBuf_[offset++] = u[index_fn(i)];
        }
    };

    // Left face
    if (subgrid.neighborLeft() != MPI_PROC_NULL) {
        for (Index k = 0; k < subgrid.nzTotal(); ++k) {
            for (Index j = 0; j < subgrid.nyTotal(); ++j) {
                for (Index h = 0; h < static_cast<Index>(hw); ++h) {
                    sendBuf_[offset++] = u[k * nyT * nxT + j * nxT + (subgrid.iEnd() - hw + h)];
                }
            }
        }
    }
    // Right face
    if (subgrid.neighborRight() != MPI_PROC_NULL) {
        for (Index k = 0; k < subgrid.nzTotal(); ++k) {
            for (Index j = 0; j < subgrid.nyTotal(); ++j) {
                for (Index h = 0; h < static_cast<Index>(hw); ++h) {
                    sendBuf_[offset++] = u[k * nyT * nxT + j * nxT + (subgrid.iBegin() + h)];
                }
            }
        }
    }
    // Down face
    if (subgrid.neighborDown() != MPI_PROC_NULL) {
        for (Index k = 0; k < subgrid.nzTotal(); ++k) {
            for (Index h = 0; h < static_cast<Index>(hw); ++h) {
                for (Index i = 0; i < subgrid.nxTotal(); ++i) {
                    sendBuf_[offset++] = u[k * nyT * nxT + (subgrid.jEnd() - hw + h) * nxT + i];
                }
            }
        }
    }
    // Up face
    if (subgrid.neighborUp() != MPI_PROC_NULL) {
        for (Index k = 0; k < subgrid.nzTotal(); ++k) {
            for (Index h = 0; h < static_cast<Index>(hw); ++h) {
                for (Index i = 0; i < subgrid.nxTotal(); ++i) {
                    sendBuf_[offset++] = u[k * nyT * nxT + (subgrid.jBegin() + h) * nxT + i];
                }
            }
        }
    }

    // For collective exchanger, we use a simple approach: each process sends its
    // packed faces to the corresponding neighbors using blocking sends/receives.
    // This is not true "collective" in the MPI sense but demonstrates the concept.
    // In production, one would use MPI_Neighbor_allgatherv with a dedicated comm.
    offset = 0;
    std::vector<MPI_Request> reqs;
    auto post_recv = [&](int neighbor, Index count) {
        if (neighbor != MPI_PROC_NULL && count > 0) {
            reqs.push_back(MPI_REQUEST_NULL);
            MPI_Irecv(&recvBuf_[offset], static_cast<int>(count), MPI_DOUBLE, neighbor, 0, comm_, &reqs.back());
            offset += count;
        }
    };
    auto post_send = [&](int neighbor, Index start, Index count) {
        if (neighbor != MPI_PROC_NULL && count > 0) {
            reqs.push_back(MPI_REQUEST_NULL);
            MPI_Isend(&sendBuf_[start], static_cast<int>(count), MPI_DOUBLE, neighbor, 0, comm_, &reqs.back());
        }
    };

    // This is a simplified placeholder. For a complete implementation, the
    // packing/unpacking and communication pattern must match exactly.
    // In this reference version, we delegate to P2P for simplicity.
    PointToPointExchanger p2p;
    p2p.initialize(subgrid);
    p2p.exchange(subgrid);
}

void CollectiveExchanger::beginExchange(Subgrid& subgrid) {
    // Not implemented in reference version; falls back to P2P
    PointToPointExchanger p2p;
    p2p.initialize(subgrid);
    p2p.beginExchange(subgrid);
    request_ = MPI_REQUEST_NULL; // placeholder
}

void CollectiveExchanger::endExchange(Subgrid& subgrid) {
    // Not implemented in reference version
    if (request_ != MPI_REQUEST_NULL) {
        MPI_Wait(&request_, MPI_STATUS_IGNORE);
        request_ = MPI_REQUEST_NULL;
    }
}

} // namespace hypo
