#include <gtest/gtest.h>

#include <csim/init_ops.hpp>
#include <csim/memory_ops.hpp>
#include <csim/utility.hpp>

#include "../util/util.hpp"

TEST(MemoryOperationTest, AllocateAndRelease) {
    const UINT n = 10;
    const ITYPE dim = 1ULL << n;
    auto ptr = allocate_quantum_state(dim);
    release_quantum_state(ptr);
}

TEST(MemoryOperationTest, MemoryZeroCheck) {
    const UINT n = 10;
    const ITYPE dim = 1ULL << n;

    auto ptr = allocate_quantum_state(dim);
    initialize_quantum_state(ptr, dim);
    for (ITYPE ind = 0; ind < dim; ++ind) {
        if (ind == 0)
            ASSERT_NEAR(_cabs(ptr[ind] - 1.), 0., eps);
        else
            ASSERT_NEAR(_cabs(ptr[ind]), 0, eps);
    }
    release_quantum_state(ptr);
}

TEST(MemoryOperationTest, HaarRandomState) {
    const UINT n = 10;
    const ITYPE dim = 1ULL << n;
    auto ptr = allocate_quantum_state(dim);
    initialize_Haar_random_state(ptr, dim);
    release_quantum_state(ptr);
}

TEST(MemoryOperationTest, LargeMemory) {
    const UINT n = 20;  // this requires about 8GB
    const ITYPE dim = 1ULL << n;
    auto ptr = allocate_quantum_state(dim);
    initialize_quantum_state(ptr, dim);
    release_quantum_state(ptr);
}
