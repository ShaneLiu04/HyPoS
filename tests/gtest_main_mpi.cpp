#include <gtest/gtest.h>
#include <mpi.h>
#include <cstdlib>

int main(int argc, char** argv) {
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}

// AR005 (design §4.1 FP1): process-count probe guarding against MPI
// runtimes that silently degrade multi-rank launches to a single
// process (mpich 4.2.0 on Ubuntu 24.04, see README "Known issues").
// Activated via HYPOS_EXPECT_NP: unset -> SKIP (the unit ctest entry
// runs this same binary as a plain single process); set -> the actual
// MPI_Comm_size must match, otherwise the suite would be "passing"
// with N-1 ranks missing. Unparsable or non-positive values fail
// safe (treated as a mismatch).
TEST(MpiEnvTest, SizeProbeMatchesExpected) {
    const char* env = std::getenv("HYPOS_EXPECT_NP");
    if (env == nullptr) {
        GTEST_SKIP() << "HYPOS_EXPECT_NP not set; probe inactive";
    }
    char* end = nullptr;
    const long expected = std::strtol(env, &end, 10);
    ASSERT_TRUE(end != env && *end == '\0' && expected >= 1)
        << "HYPOS_EXPECT_NP is not a positive integer: '" << env << "'";
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    EXPECT_EQ(size, expected)
        << "MPI launch degraded: expected " << expected
        << " ranks but MPI_Comm_size reports " << size
        << " (this MPI runtime may be silently collapsing multi-rank "
           "launches to a single process)";
}
