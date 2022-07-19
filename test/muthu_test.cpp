/***
BSD 2-Clause License

Copyright (c) 2018, Adrián
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**/


//
// Created by Adrián on 14/7/22.
//

#include <sdsl/wm_int.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <unordered_map>
#include <unordered_set>
#include "gtest/gtest.h"


using namespace sdsl;
using namespace std;

uint64_t n;
uint64_t sigma;

namespace {

    template<class T>
    class muthu_test : public ::testing::Test {
    };

    using testing::Types;

    typedef Types<bit_vector> Implementations;


    uint64_t distinct(int_vector<> &vec, uint64_t lb, uint64_t rb) {
        std::unordered_set<uint64_t> set;
        for (uint64_t i = lb; i <= rb; ++i) {
            set.insert(vec[i]);
        }
        return set.size();
    }


    TYPED_TEST_CASE(muthu_test, Implementations);

    TYPED_TEST(muthu_test, count) {
        int_vector<> data(n + 1);
        data[0] = 0;
        for (uint64_t i = 1; i < data.size(); ++i) {
            data[i] = rand() % sigma;
        }
        int_vector<> muthu(n + 1);
        muthu[0] = 0;
        std::unordered_map<uint64_t, uint64_t> hash_table;
        for (uint64_t i = 1; i < muthu.size(); ++i) {
            auto ht_it = hash_table.find(data[i]);
            if (ht_it == hash_table.end()) {
                muthu[i] = 0;
                hash_table.insert({data[i], i});
            } else {
                muthu[i] = ht_it->second;
                ht_it->second = i;
            }
        }

        wm_int<bit_vector> wt_muthu;
        wm_int<bit_vector> wt_data;
        sdsl::construct_im(wt_muthu, muthu);
        sdsl::construct_im(wt_data, data);

        uint64_t length = 10;
        for (uint64_t i = 1; i + length - 1 < muthu.size(); ++i) {
            auto obtained = wt_muthu.count_range_search_2d(i, i + length - 1, 0, i - 1);
            auto expected = wt_data.count_distinct_values(i, i + length-1);
            std::cout << i << std::endl;
            ASSERT_EQ(expected, obtained);
        }
    }
}


int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 2) {
        // LCOV_EXCL_START
        cout << "Usage: " << argv[0] << " N SIGMA " << endl;
        return 1;
        // LCOV_EXCL_STOP
    }
    n = atoll(argv[1]);
    sigma = atoll(argv[2]);
    return RUN_ALL_TESTS();
}
