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
// Created by Adrián on 19/10/2019.
//

#ifndef SDSL_SUCC_SUPPORT_SD_HPP
#define SDSL_SUCC_SUPPORT_SD_HPP

#include "vectors.hpp"
#include "sd_vector.hpp"
#include "succ_support_v.hpp"
#include "definitions.hpp"
#include "math.hpp"

namespace sdsl {

    template<uint8_t t_b          = 1,
            class t_hi_bit_vector= bit_vector,
            class t_select_1     = typename t_hi_bit_vector::select_1_type,
            class t_select_0     = typename t_hi_bit_vector::select_0_type>
    class succ_support_sd;
    //! Select data structure for sd_vector
/*! \tparam t_b             Bit pattern.
 *  \tparam t_hi_bit_vector Type of the bitvector used for the unary decoded differences of
 *                          the high part of the positions of the 1s.
 *  \tparam t_select_1      Type of the select structure which is used to select ones in HI.
 *  \tparam t_select_0      Type of the select structure which is used to select zeros in HI.
 */
    template<uint8_t t_b, class t_hi_bit_vector, class t_select_1, class t_select_0>
    class succ_support_sd
    {
    public:
        typedef bit_vector::size_type size_type;
        typedef sdsl::sd_vector<t_hi_bit_vector, t_select_1, t_select_0> bit_vector_type;
        enum { bit_pat = t_b };
        enum { bit_pat_len = (uint8_t)1 };
    private:
        const bit_vector_type* m_v;
        succ_support_v<1> m_succ_high;

        void copy(const succ_support_sd& ss){
            m_v = ss.m_v;
            m_succ_high = ss.m_succ_high;
            if(m_v != nullptr){
                m_succ_high.set_vector(&(m_v->high));
            }
        }
    public:

        succ_support_sd(){};

        explicit succ_support_sd(const bit_vector_type* v)
        {
            set_vector(v);
            if(v != nullptr){
                sdsl::util::init_support(m_succ_high, &(m_v->high));
            }
        }

        //! Returns the position of the i-th occurrence in the bit vector.
        size_type succ(size_type i)const
        {
            size_type high_val = (i >> (m_v->wl));
            size_type val_low = i & bits::lo_set[ m_v->wl ];
            size_type sel_high = m_v->high_0_select(high_val+1);
            size_type rank_low = sel_high - high_val;


            int64_t r = rank_low;
            int64_t s = sel_high;
            do {
                --s; --r;
            }while(s >= 0 and m_v->high[s] and m_v->low[r] >= val_low);

            if(m_v->high[s+1]){
                return m_v->low[r+1] + ((high_val) << m_v->wl);
            }else{
                size_type succ_high = m_succ_high(sel_high);
                if(succ_high < m_v->high.size()){
                    return m_v->low[rank_low] +
                           ((high_val+(succ_high-sel_high)) << m_v->wl);
                }
                return m_v->size();
            }

        }

        size_type operator()(size_type i)const
        {
            return succ(i);
        }

        size_type size()const
        {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=nullptr)
        {
            m_v = v;
        }

        succ_support_sd(const succ_support_sd& p){
            copy(p);
        };

        succ_support_sd(succ_support_sd&& p){
            *this = std::move(p);
        };

        succ_support_sd& operator=(const succ_support_sd& ss)
        {
            if (this != &ss) {
                copy(ss);
            }
            return *this;
        }

        succ_support_sd& operator=(succ_support_sd&& ss)
        {
            if (this != &ss) {
                m_v = std::move(ss.m_v);
                m_succ_high = std::move(ss.m_succ_high);
                if(m_v != nullptr){
                    m_succ_high.set_vector(&(m_v->high));
                }
            }
            return *this;
        }

        void swap(succ_support_sd& ss) {
            std::swap(m_v, ss.m_v);
            m_succ_high.swap(ss.m_succ_high);
            if(m_v != nullptr){
                m_succ_high.set_vector(&(m_v->high));
            }else{
                m_succ_high.set_vector(nullptr);
            }
            if(ss.m_v != nullptr){
                ss.m_succ_high.set_vector(&(ss.m_v->high));
            }else{
                ss.m_succ_high.set_vector(nullptr);
            }
        }

        void load(std::istream& in, const bit_vector_type* v=nullptr)
        {
            set_vector(v);
            if(m_v != nullptr){
                m_succ_high.load(in, &(m_v->high));
            }

        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_succ_high.serialize(out, child, "succ_high");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }
    };


    template<uint8_t t_b          = 1,
            class t_hi_bit_vector= bit_vector,
            class t_select_1     = typename t_hi_bit_vector::select_1_type,
            class t_select_0     = typename t_hi_bit_vector::select_0_type>
    class succ_support_naive_sd;
    //! Select data structure for sd_vector
/*! \tparam t_b             Bit pattern.
 *  \tparam t_hi_bit_vector Type of the bitvector used for the unary decoded differences of
 *                          the high part of the positions of the 1s.
 *  \tparam t_select_1      Type of the select structure which is used to select ones in HI.
 *  \tparam t_select_0      Type of the select structure which is used to select zeros in HI.
 */
    template<uint8_t t_b, class t_hi_bit_vector, class t_select_1, class t_select_0>
    class succ_support_naive_sd
    {
    public:
        typedef bit_vector::size_type size_type;
        typedef sdsl::sd_vector<t_hi_bit_vector, t_select_1, t_select_0> bit_vector_type;
        enum { bit_pat = t_b };
        enum { bit_pat_len = (uint8_t)1 };
    private:
        const bit_vector_type* m_v;

        void copy(const succ_support_naive_sd& ss){
            m_v = ss.m_v;
        }
    public:

        explicit succ_support_naive_sd(const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        //! Returns the position of the i-th occurrence in the bit vector.
        size_type succ(size_type i)const
        {
            size_type high_val = (i >> (m_v->wl));
            size_type val_low = i & bits::lo_set[ m_v->wl ];
            size_type sel_high = m_v->high_0_select(high_val+1);
            size_type rank_low = sel_high - high_val;


            int64_t r = rank_low;
            int64_t s = sel_high;
            do {
                --s; --r;
            }while(s >= 0 and m_v->high[s] and m_v->low[r] >= val_low);

            if(m_v->high[s+1]){
                return m_v->low[r+1] + ((high_val) << m_v->wl);
            }else{
                uint64_t* data = m_v->high.data();
                size_type succ_high = bits::next(data, sel_high, m_v->high.size());
                //size_type succ_high = m_succ_high(sel_high);
                if(succ_high < m_v->high.size()){
                    return m_v->low[rank_low] +
                           ((high_val+(succ_high-sel_high)) << m_v->wl);
                }
                return m_v->size();
            }

        }

        size_type operator()(size_type i)const
        {
            return succ(i);
        }

        size_type size()const
        {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=nullptr)
        {
            m_v = v;
        }

        succ_support_naive_sd& operator=(const succ_support_naive_sd& ss)
        {
            if (this != &ss) {
                copy(ss);
            }
            return *this;
        }

        succ_support_naive_sd& operator=(succ_support_naive_sd&& ss)
        {
            if (this != &ss) {
                m_v = std::move(ss.m_v);
            }
            return *this;
        }

        void swap(succ_support_naive_sd& ss) {
            std::swap(m_v, ss.m_v);
        }

        void load(std::istream& in, const bit_vector_type* v=nullptr)
        {
            set_vector(v);

        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            sdsl::structure_tree::add_size(child, 0);
            return 0;
        }
    };
}

#endif //RUNS_VECTORS_SUCC_SUPPORT_SD_HPP
