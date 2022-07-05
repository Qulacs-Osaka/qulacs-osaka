import argparse
from typing import Optional


def split_into_report_chunk(raw_output: str) -> list[str]:
    """Split raw output into reports of each memory leak report.
    They are separated by a blank line.
    """
    chunks = raw_output.split("\n\n")
    # Skip chunk including only new line character.
    chunks = list(filter(lambda chunk: len(chunk) > 0, chunks))
    return chunks


def extract_root_cause_of_memory_leak(report_chunk: str) -> Optional[str]:
    """From memory leak report, extract a file location of initial cause of it in this library.

    # Example
        The output for following report is `src/cppsim/gate_factory.cpp:54` from #8.

        Indirect leak of 16 byte(s) in 1 object(s) allocated from:
        #0 0x7f3a488d4647 in operator new(unsigned long) ../../../../src/libsanitizer/asan/asan_new_delete.cpp:99
        #1 0x55798f60c272 in __gnu_cxx::new_allocator<TargetQubitInfo>::allocate(unsigned long, void const*) /usr/include/c++/10/ext/new_allocator.h:115
        #2 0x55798f60c272 in std::allocator_traits<std::allocator<TargetQubitInfo> >::allocate(std::allocator<TargetQubitInfo>&, unsigned long) /usr/include/c++/10/bits/alloc_traits.h:460
        #3 0x55798f60c272 in std::_Vector_base<TargetQubitInfo, std::allocator<TargetQubitInfo> >::_M_allocate(unsigned long) /usr/include/c++/10/bits/stl_vector.h:346
        #4 0x55798f60c272 in void std::vector<TargetQubitInfo, std::allocator<TargetQubitInfo> >::_M_realloc_insert<TargetQubitInfo>(__gnu_cxx::__normal_iterator<TargetQubitInfo*, std::vector<TargetQubitInfo, std::allocator<TargetQubitInfo> > >, TargetQubitInfo&&) /usr/include/c++/10/bits/vector.tcc:440
        #5 0x55798f5a82f0 in void std::vector<TargetQubitInfo, std::allocator<TargetQubitInfo> >::emplace_back<TargetQubitInfo>(TargetQubitInfo&&) /usr/include/c++/10/bits/vector.tcc:121
        #6 0x55798f5a82f0 in std::vector<TargetQubitInfo, std::allocator<TargetQubitInfo> >::push_back(TargetQubitInfo&&) /usr/include/c++/10/bits/stl_vector.h:1204
        #7 0x55798f5a82f0 in ClsP0Gate::ClsP0Gate(unsigned int) /workspaces/qulacs-osaka/src/cppsim/gate_named_one.hpp:359
        #8 0x55798f5a82f0 in gate::P0(unsigned int) /workspaces/qulacs-osaka/src/cppsim/gate_factory.cpp:54
        #9 0x55798f426df0 in DensityMatrixGeneralGateTest_InstrumentGate_Test::TestBody() /workspaces/qulacs-osaka/test/cppsim/test_noise_dm.cpp:49
        #10 0x55798f78c7dc in void testing::internal::HandleExceptionsInMethodIfSupported<testing::Test, void>(testing::Test*, void (testing::Test::*)(), char const*) (/workspaces/qulacs-osaka/bin/cppsim_test+0x8387dc)
    """
    if "==" in report_chunk:
        # First few lines of stderr are summary, so ignore them.
        return None
    report_lines = report_chunk.split("\n")
    file_names_of_this_library = list(
        filter(lambda name: "qulacs-osaka/src" in name, report_lines)
    )
    if len(file_names_of_this_library) > 0:
        return file_names_of_this_library[-1]
    return None


def main():
    parser = argparse.ArgumentParser(
        "Extract file names of memory leak from the stderr of a program using address sanitizer."
    )
    parser.add_argument("file_name", type=str, help="File name of stderr dump.")
    args = parser.parse_args()

    with open(args.file_name, "r") as f:
        raw_output = f.read()
        chunks = split_into_report_chunk(raw_output)
        initial_causes = [extract_root_cause_of_memory_leak(chunk) for chunk in chunks]
        unique_initial_causes = set(initial_causes)
        for cause in unique_initial_causes:
            if cause is not None:
                print(cause)


if __name__ == "__main__":
    main()
