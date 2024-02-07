#include "sdsl/bit_vectors.hpp"
#include "sdsl/succ_support_v.hpp"
#include "sdsl/succ_support_sd.hpp"
#include "gtest/gtest.h"
#include <string>

using namespace sdsl;
using namespace std;

string test_file;

namespace
{

template<class T>
class succ_support_test : public ::testing::Test { };

using testing::Types;

typedef Types<succ_support_v<1>,
        succ_support_rrr<1, 256>,
        succ_support_rrr<1, 129>,
        succ_support_rrr<1, 192>,
        succ_support_rrr<1, 255>,
        succ_support_rrr<1, 15>,
        succ_support_rrr<1, 31>,
        succ_support_rrr<1, 63>,
        succ_support_rrr<1, 127>,
        succ_support_rrr<1, 128>,
        succ_support_sd<1>,
        succ_support_hyb<1>
        > Implementations;

TYPED_TEST_CASE(succ_support_test, Implementations);

//! Test the select method
TYPED_TEST(succ_support_test, succ_method)
{
    static_assert(sdsl::util::is_regular<TypeParam>::value, "Type is not regular");
    bit_vector bvec;
    ASSERT_TRUE(load_from_file(bvec, test_file));
    typename TypeParam::bit_vector_type bv(bvec);
    TypeParam ss(&bv);
    for (uint64_t j=0; j < bvec.size(); ++j) {
        uint64_t succ;
        for (succ=j; succ < bvec.size(); ++succ) {
            if(bvec[succ]) break;
        }
        if (succ < bvec.size()) {
            ASSERT_EQ(succ, ss.succ(j));
        }
    }
}

}// end namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 2) {
        // LCOV_EXCL_START
        cout << "Usage: " << argv[0] << " FILE " << endl;
        cout << "  Reads a bitvector from FILE and executes tests." << endl;
        return 1;
        // LCOV_EXCL_STOP
    }
    test_file = argv[1];
    return RUN_ALL_TESTS();
}
