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
// Created by Adrián on 16/12/2018.
//

#ifndef RUNS_VECTORS_PREV_SUPPORT_V_HPP
#define RUNS_VECTORS_PREV_SUPPORT_V_HPP

#include "vectors.hpp"
#include "definitions.hpp"
#include "math.hpp"


namespace sdsl {

    template<uint8_t t_b=1>
    class prev_support_v
    {
    private:
        static_assert(t_b == 1u or t_b == 0u, "prev_support_v: bit pattern must be `0`or `1`");
    public:
        typedef sdsl::bit_vector                    bit_vector_type;
        typedef sdsl::int_vector_size_type          size_type;
        enum { bit_pat = t_b };
    private:
        const sdsl::bit_vector* m_v;  //!< Pointer to the previous supported bit_vector
        sdsl::int_vector<> m_super_blocks; //super blocks
        sdsl::int_vector<> m_basic_blocks; //basic blocks

        void copy(const prev_support_v& p){
            m_v = p.m_v;
            m_super_blocks = p.m_super_blocks;
            m_basic_blocks = p.m_basic_blocks;
        }

    public:
        explicit prev_support_v(const sdsl::bit_vector* v = nullptr) {
            set_vector(v);
            if (v == nullptr || v->empty()) {
                m_basic_blocks = sdsl::int_vector<>(1 ,0, 1);
                m_super_blocks = sdsl::int_vector<>(1 ,0, 1);
                return ;
            }
            size_type super_block_size = util::math::ceil_div(v->capacity(), W2);
            size_type basic_block_size = util::math::ceil_div(v->capacity(), W);
            m_super_blocks.resize(super_block_size);   // resize structure for super_blocks
            m_basic_blocks.resize(basic_block_size);   // resize structure for basic_blocks
            size_type i, j=0;
            m_basic_blocks[0] = 0;
            m_super_blocks[0] = m_v->size();
            int64_t pred = v->size();
            const uint64_t* data = m_v->data();
            for(i=0; i < basic_block_size-1; i++){
                if(i >0 && i % W == 0){
                    m_super_blocks[(i / W)] = pred;
                }
                uint64_t d = (t_b) ? sdsl::bits::rev(data[i]) : sdsl::bits::rev(~data[i]);
                if(d){
                    pred = i*W + (W-1) - sdsl::bits::lo(d); //right-most bit pred(B, (i+1)*b)
                }
                size_type bb_value = std::min((i+1)*W-pred, (size_type) (W2));
                m_basic_blocks[i+1] = bb_value;
            }
            //sdsl::util::bit_compress(m_super_blocks);
            //sdsl::util::bit_compress(m_basic_blocks);
        }
        prev_support_v(const prev_support_v& p){
            copy(p);
        };

        prev_support_v(prev_support_v&& p){
            *this = std::move(p);
        };

        prev_support_v& operator=(const prev_support_v& p){
            if (this != &p) {
                copy(p);
            }
            return *this;
        };

        prev_support_v& operator=(prev_support_v&& p){
            if (this != &p) {
                m_super_blocks = std::move(p.m_super_blocks);
                m_basic_blocks = std::move(p.m_basic_blocks);
                m_v = std::move(p.m_v);
            }
            return *this;
        };


        uint64_t clear(const uint64_t bits, size_type idx) const{
            auto displacement = (W-1) - idx;
            return (bits >> displacement) << displacement;
        }

        uint64_t clear_rev(const uint64_t bits, size_type idx) const{
            auto displacement = (W-1) - idx;
            return (bits << displacement) >> displacement;
        }

        uint64_t prev(const size_type idx) const;




        inline uint64_t operator()(size_type idx)const {
            return prev(idx);
        }

        size_type size()const {
            return m_v->size();
        }

        void set_vector(const sdsl::bit_vector* v=nullptr) {
            m_v = v;
        }


        void swap(prev_support_v& ps) {
            if (this != &ps) { // if rs and _this_ are not the same object
                m_basic_blocks.swap(ps.m_basic_blocks);
                m_super_blocks.swap(ps.m_super_blocks);

            }
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr,
                            std::string name="")const {
            size_type written_bytes = 0;
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            written_bytes += m_basic_blocks.serialize(out, child, "basic_blocks");
            written_bytes += m_super_blocks.serialize(out, child, "super_blocks");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream &in){
            m_basic_blocks.load(in);
            m_super_blocks.load(in);
        }
    };


    template <uint8_t t_b>
    inline uint64_t prev_support_v<t_b>::prev(const size_type idx) const {
        assert(m_v != nullptr);
        assert(idx <= m_v->size());
        if(idx < 0) return m_v->size();
        const uint64_t* data = m_v->data();
        uint64_t b_block_idx = (idx >> 6); //basic block index
        uint64_t b_clear = clear_rev(data[b_block_idx], (idx & (W-1)));
        if(b_clear) {//if there are 1s
            return (b_block_idx << 6) + sdsl::bits::hi(b_clear);
        }else if(b_block_idx == 0){
            return m_v->size();
        }else if(m_basic_blocks[b_block_idx] < W2){
            return b_block_idx *W - m_basic_blocks[b_block_idx];
        }else{
            return m_super_blocks[(idx >> 12)];
        }
    }


    template <>
    inline uint64_t prev_support_v<0>::prev(const size_type idx) const {
        assert(m_v != nullptr);
        assert(idx <= m_v->size());
        if(idx < 0) return m_v->size();
        const uint64_t* data = m_v->data();
        uint64_t b_block_idx = (idx >> 6); //basic block index
        uint64_t b_clear = clear_rev(~data[b_block_idx], (idx & (W-1)));
        if(b_clear) {//if there are 1s
            return (b_block_idx << 6) + sdsl::bits::hi(b_clear);
        }else if(b_block_idx == 0){
            return m_v->size();
        }else if(m_basic_blocks[b_block_idx] < W2){
            return b_block_idx *W - m_basic_blocks[b_block_idx];
        }else{
            return m_super_blocks[(idx >> 12)];
        }
    }

}


#endif //RCT_PREV_SUPPORT_V_HPP
