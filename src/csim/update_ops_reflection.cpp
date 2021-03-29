
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constant.hpp"
#include "stat_ops.hpp"
#include "update_ops.hpp"
#include "utility.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif

void reflection_gate(const CTYPE* reflection_state, CTYPE* state, ITYPE dim) {
    CTYPE coef = state_inner_product(reflection_state, state, dim);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (ITYPE state_index = 0; state_index < dim; ++state_index) {
        state[state_index] =
            2.0 * coef * reflection_state[state_index] - state[state_index];
    }
}
