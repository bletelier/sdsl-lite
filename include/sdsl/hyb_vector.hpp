/*
 *   Copyright (c) 2014 Juha Karkkainen, Dominik Kempa and Simon J. Puglisi
 *
 *   Permission is hereby granted, free of charge, to any person
 *   obtaining a copy of this software and associated documentation
 *   files (the "Software"), to deal in the Software without
 *   restriction, including without limitation the rights to use,
 *   copy, modify, merge, publish, distribute, sublicense, and/or sell
 *   copies of the Software, and to permit persons to whom the
 *   Software is furnished to do so, subject to the following
 *   conditions:
 *
 *   The above copyright notice and this permission notice shall be
 *   included in all copies or substantial portions of the Software.
 *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *   OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *   HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *   WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *   OTHER DEALINGS IN THE SOFTWARE.
 *
 *   Simon Gog made the following changes:
 *     - replace std::vectors by int_vectors
 *     - add support for rank0
 *     - added naive implementation of method get_int
 *     - TODO: added a naive implementation of select
*/
#ifndef INCLUDED_SDSL_HYB_VECTOR
#define INCLUDED_SDSL_HYB_VECTOR

#include "int_vector.hpp"
#include "util.hpp"
#include "iterators.hpp"
#include "io.hpp"
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iostream>

namespace sdsl
{

// Needed for friend declarations.
    template<uint8_t t_b = 1, uint32_t k_sb_rate = 16> class rank_support_hyb;
    template<uint8_t t_b = 1, uint32_t k_sb_rate = 16> class select_support_hyb;
    template<uint8_t t_b = 1, uint32_t k_sb_rate = 16> class succ_support_hyb;

//! A hybrid-encoded compressed bitvector representation
/*!
 * \tparam k_sblock_rate  Superblock rate (number of blocks inside superblock)
 *
 * References:
 *  - Juha Karkkainen, Dominik Kempa and Simon J. Puglisi.
 *    Hybrid Compression of Bitvectors for the FM-Index.
 *    DCC 2014.
 */
template<uint32_t k_sblock_rate = 16>
class hyb_vector
{
    public:
        typedef bit_vector::size_type size_type;
        typedef bit_vector::value_type value_type;
        typedef bit_vector::difference_type difference_type;
        typedef random_access_const_iterator<hyb_vector> iterator;
        typedef rank_support_hyb<1, k_sblock_rate> rank_1_type;
        typedef rank_support_hyb<0, k_sblock_rate> rank_0_type;
        typedef select_support_hyb<1, k_sblock_rate> select_1_type;
        typedef select_support_hyb<0, k_sblock_rate> select_0_type;
        typedef succ_support_hyb<1, k_sblock_rate> succ_1_type;
        typedef succ_support_hyb<0, k_sblock_rate> succ_0_type;

        friend class rank_support_hyb<1, k_sblock_rate>;
        friend class rank_support_hyb<0, k_sblock_rate>;
        friend class select_support_hyb<1, k_sblock_rate>;
        friend class select_support_hyb<0, k_sblock_rate>;
        friend class succ_support_hyb<1, k_sblock_rate>;
        friend class succ_support_hyb<0, k_sblock_rate>;

    private:
        static const uint32_t k_block_size;
        static const uint32_t k_block_bytes;
        static const uint32_t k_sblock_header_size;
        static const uint32_t k_sblock_size;
        static const uint32_t k_hblock_rate;

        size_type m_size = 0;           // original bitvector size
        int_vector<8> m_trunk;          // body of encoded blocks
        int_vector<8> m_sblock_header;  // sblock headers
        int_vector<64> m_hblock_header; // hblock headers

        void copy(const hyb_vector& hybrid)
        {
            m_size = hybrid.m_size;
            m_trunk = hybrid.m_trunk;
            m_sblock_header = hybrid.m_sblock_header;
            m_hblock_header = hybrid.m_hblock_header;
        }

    public:
        //! Default constructor
        hyb_vector() = default;

        //! Copy constructor
        hyb_vector(const hyb_vector& hybrid)
        {
            copy(hybrid);
        }

        //! Move constructor
        hyb_vector(hyb_vector&& hybrid)
            : m_size(std::move(hybrid.m_size)),
              m_trunk(std::move(hybrid.m_trunk)),
              m_sblock_header(std::move(hybrid.m_sblock_header)),
              m_hblock_header(std::move(hybrid.m_hblock_header)) {}

        //! Constructor
        hyb_vector(const bit_vector& bv)
        {
            m_size = bv.size();

            // Compute the number of blocks.
            size_type n_blocks = (m_size + k_block_size - 1) / k_block_size;
            size_type n_sblocks = (n_blocks + k_sblock_rate - 1) / k_sblock_rate;
            size_type n_hblocks = (n_blocks + k_hblock_rate - 1) / k_hblock_rate;

            size_type trunk_size = 0;

            // runs_lookup[i] = number of runs - 1 in the binary encoding of i.
            int_vector<8> runs_lookup(65536,0);
            runs_lookup[0] = 0;
            for (uint32_t i = 1; i < 65536; ++i) {
                runs_lookup[i] = runs_lookup[i >> 1];
                if (i >= 32768) --runs_lookup[i];
                if ((i & 1) != ((i >> 1) & 1)) ++runs_lookup[i];
            }

            // Compute optimal encoding for each block.
            const uint64_t* bv_ptr = bv.data();
            for (size_type block_id = 0; block_id < n_blocks; ++block_id) {
                size_type block_beg = block_id * k_block_size;
                size_type block_end = block_beg + k_block_size;

                uint32_t ones = 0;
                uint32_t runs = 0;

                if (block_end <= m_size) {
                    // Count the number of ones, fast.
                    const uint64_t* ptr64 = bv_ptr;
                    for (uint8_t i = 0; i < 4; ++i) ones += bits::cnt(*ptr64++);

                    // Count the number of runs, fast.
                    ptr64 = bv_ptr;
                    for (uint8_t i = 0; i < 4; ++i) {
                        // Count changes of bits inside 16-bit words of *ptr64.
                        for (uint8_t j = 0; j < 4; ++j) runs += runs_lookup[((*ptr64)>>(16*j))&0xffff];

                        // Count changes of bits between 16-bit words of *ptr64.
                        for (uint8_t j = 0; j < 3; ++j) runs += ((((*ptr64)>>(16*j+15))&1) ^ (((*ptr64)>>(16*j+16))&1));
                        ++ptr64;
                    }

                    // Count changes of bits between 64-bit words.
                    ptr64 = bv_ptr;
                    for (uint8_t i = 0; i < 3; ++i) {
                        runs += ((((*ptr64)>>63)&1) ^ ((*(ptr64 + 1))&1));
                        ++ptr64;
                    }
                    ++runs;
                } else {
                    // Count number of ones and runs, slow.
                    uint8_t prevbit = 2;
                    for (size_type i = block_beg; i < block_end; ++i) {
                        uint8_t bit = (i < m_size ? bv[i] : 0);
                        if (bit == 1) ++ones;
                        if (bit != prevbit) ++runs;
                        prevbit = bit;
                    }
                }

                // Choose best encoding.
                uint32_t minority_enc_size = std::min(ones, k_block_size - ones);
                uint32_t runs_enc_size = (uint32_t)std::max(0, (int32_t)runs - 2);
                uint32_t best_enc_size = std::min(minority_enc_size, runs_enc_size);
                best_enc_size = std::min(best_enc_size, k_block_bytes);

                // Update the answer.
                trunk_size += best_enc_size;
                bv_ptr += k_block_size / 64;
            }

            // Allocate the memory.
            m_sblock_header = int_vector<8>(n_sblocks * k_sblock_header_size, 0);
            m_hblock_header = int_vector<64>(n_hblocks * 2, 0);
            m_trunk         = int_vector<8>(trunk_size, 0);

            // The actual encoding follows.
            size_type tot_rank = 0;    // stores current rank value
            size_type sblock_ones = 0; // number of 1s inside superblock
            size_type trunk_ptr = 0;

            // Process blocks left to right.
            bv_ptr = bv.data();
            for (size_type block_id = 0; block_id < n_blocks; ++block_id) {
                size_type block_beg = block_id * k_block_size;
                size_type block_end = block_beg + k_block_size;
                size_type sblock_id = block_id / k_sblock_rate;
                size_type hblock_id = block_id / k_hblock_rate;

                // Update hblock header.
                if (!(block_id % k_hblock_rate)) {
                    m_hblock_header[2 * hblock_id] = trunk_ptr;
                    m_hblock_header[2 * hblock_id + 1] = tot_rank;
                }

                // Update sblock header.
                if (!(block_id % k_sblock_rate)) {
                    uint32_t* ptr = (uint32_t*)(((uint8_t*)m_sblock_header.data()) + k_sblock_header_size * sblock_id);
                    *ptr++ = trunk_ptr - m_hblock_header[2 * hblock_id];
                    *ptr = tot_rank - m_hblock_header[2 * hblock_id + 1];

                    // If the sblock is uniform, flip the bit.
                    if (sblock_id && (!sblock_ones || sblock_ones == k_sblock_size)) {
                        ptr = (uint32_t*)(((uint8_t*)m_sblock_header.data()) + k_sblock_header_size * (sblock_id - 1));
                        *ptr |= 0x80000000;
                    }

                    // Reset the number of ones in sblock.
                    sblock_ones = 0;
                }

                uint32_t ones = 0;
                uint32_t runs = 0;

                // Compute the number of 1-bits and runs inside current block.
                if (block_end <= m_size) {
                    // Count the number of ones, fast.
                    const uint64_t* ptr64 = bv_ptr;
                    for (uint8_t i = 0; i < 4; ++i) ones += bits::cnt(*ptr64++);

                    // Count the number of runs, fast.
                    ptr64 = bv_ptr;
                    for (uint8_t i = 0; i < 4; ++i) {
                        for (uint8_t j = 0; j < 4; ++j) runs += runs_lookup[((*ptr64)>>(16*j))&0xffff];
                        for (uint8_t j = 0; j < 3; ++j) runs += ((((*ptr64)>>(16*j+15))&1) ^ (((*ptr64)>>(16*j+16))&1));
                        ++ptr64;
                    }
                    ptr64 = bv_ptr;
                    for (uint8_t i = 0; i < 3; ++i) {
                        runs += ((((*ptr64)>>63)&1) ^ ((*(ptr64 + 1))&1));
                        ++ptr64;
                    }
                    ++runs;
                } else {
                    // Count number of ones and runs, slow.
                    uint8_t prevbit = 2;
                    for (size_type i = block_beg; i < block_end; ++i) {
                        uint8_t bit = (i < m_size ? bv[i] : 0);
                        if (bit == 1) ++ones;
                        if (bit != prevbit) ++runs;
                        prevbit = bit;
                    }
                }
                uint32_t zeros = k_block_size - ones;

                // Store block popcount.
                uint16_t* header_ptr16 = (uint16_t*)(((uint8_t*)m_sblock_header.data()) +
                                                     sblock_id * k_sblock_header_size + 8 + (block_id % k_sblock_rate) * 2);

                (*header_ptr16) = ones;
                if (ones == k_block_size)
                    (*header_ptr16) |= 0x200;

                if (0 < ones && ones < k_block_size) { // non uniform block
                    uint32_t minority_enc_size = std::min(ones, zeros);
                    uint32_t runs_enc_size = (uint32_t)std::max(0, (int32_t)runs - 2);
                    uint32_t best_enc_size = std::min(minority_enc_size, runs_enc_size);

                    if (k_block_bytes <= best_enc_size) {
                        // Use plain encoding.
                        (*header_ptr16) |= (k_block_bytes << 10);

                        // Copy original 256 bits from bv into trunk.
                        if (block_end <= m_size) {
                            for (uint8_t i = 0; i < 4; ++i) {
                                *((uint64_t*)(((uint8_t*)m_trunk.data()) + trunk_ptr)) = *(bv_ptr + i);
                                trunk_ptr += 8;
                            }
                        } else {
                            for (size_type i = block_beg; i < block_end; i += 64) {
                                uint64_t w = 0;
                                for (size_type j = i; j < std::min(i + 64, block_end); ++j) {
                                    uint8_t bit = (j < m_size ? bv[j] : 0);
                                    if (bit) w |= ((uint64_t)1 << (j - i));
                                }
                                *((uint64_t*)(((uint8_t*)m_trunk.data()) + trunk_ptr)) = w;
                                trunk_ptr += 8;
                            }
                        }
                    } else {
                        if (runs_enc_size < minority_enc_size) {

                            // Use runs encoding.
                            (*header_ptr16) |= (runs_enc_size << 10);
                            (*header_ptr16) |= (bv[block_beg] << 9);

                            if (block_end <= m_size) {
                                // Find run ends, fast.
                                uint32_t runid = 0;
                                const uint64_t* ptr64 = bv_ptr;

                                uint64_t w = 0;
                                for (uint8_t i = 0; runid < runs_enc_size && i < 4; ++i) {
                                    // Check if run end aligns with the end of the 64-bit word.
                                    if (i > 0 && (w & 1) != ((*ptr64) & 1))
                                        m_trunk[trunk_ptr + runid++] = 64 * i - 1;

                                    w = (*ptr64++);
                                    for (uint8_t j = 0; runid < runs_enc_size && j < 63; ++j) {
                                        if ((w & 1) != ((w >> 1) & 1))
                                            m_trunk[trunk_ptr + runid++] = j + i * 64;
                                        w >>= 1;
                                    }
                                }
                                trunk_ptr += runid;
                            } else {
                                // Find run ends, slow.
                                uint8_t prevbit = 2;
                                uint32_t runid = 0;

                                for (size_type i = block_beg; runid < runs_enc_size; ++i) {
                                    uint8_t bit = (i < m_size ? bv[i] : 0);
                                    if (bit != prevbit && i != block_beg)
                                        m_trunk[trunk_ptr + runid++] = (i - block_beg - 1);
                                    prevbit = bit;
                                }
                                trunk_ptr += runid;
                            }
                        } else {
                            // Use minority encoding.
                            // Update sblock header.
                            (*header_ptr16) |= (minority_enc_size << 10);
                            if (ones < zeros)(*header_ptr16) |= 0x200;
                            uint8_t keybit = (ones < zeros);

                            // Find positions of 1-bits, fast.
                            if (block_end <= m_size) {
                                const uint64_t* ptr64 = bv_ptr;
                                for (uint8_t i = 0; i < 4; ++i) {
                                    uint64_t w = (*ptr64++);
                                    for (uint8_t j = 0; j < 64; ++j) {
                                        if ((w & 1) == keybit)
                                            m_trunk[trunk_ptr++] = j + 64 * i;
                                        w >>= 1;
                                    }
                                }
                            } else {
                                for (size_type i = block_beg; i < block_end; ++i) {
                                    uint8_t bit = (i < m_size ? bv[i] : 0);

                                    if (bit == keybit)
                                        m_trunk[trunk_ptr++] = i - block_beg;
                                }
                            }
                        }
                    }
                }

                // Update global rank.
                tot_rank += ones;
                sblock_ones += ones;
                bv_ptr += k_block_size / 64;
            }
        }

    private:
        //! Given i returns bv[i - 1].
        value_type access0(size_type i) const
        {
            assert(i > 0);
            assert(i <= m_size);

            size_type block_id = (i - 1) / k_block_size;
            size_type sblock_id = block_id / k_sblock_rate;
            size_type hblock_id = block_id / k_hblock_rate;

            size_type trunk_base = m_hblock_header[2 * hblock_id];

            uint32_t local_i = i - block_id * k_block_size;

            // Read superblock header.
            const uint8_t* header_ptr8 = ((const uint8_t*)m_sblock_header.data()) + (sblock_id * k_sblock_header_size);
            uint32_t* header_ptr32 = (uint32_t*)header_ptr8;
            size_type trunk_ptr = trunk_base + ((*header_ptr32) & 0x3fffffff);
            header_ptr8 += 8;

            uint16_t* header_ptr16 = (uint16_t*)header_ptr8;

            // Uniform superblock optimization.
            if ((*header_ptr32) & 0x80000000)
                return (value_type)((*(header_ptr8 + 1)) & 0x01);

            // Fast forward through preceding blocks in the superblock.
            for (size_type j = sblock_id * k_sblock_rate; j != block_id; ++j) {
                trunk_ptr += ((*header_ptr16) >> 10);     // Update trunk pointer.
                ++header_ptr16;
            }

            const uint8_t* trunk_p = ((const uint8_t*)m_trunk.data()) + trunk_ptr;

            uint32_t encoding_size = ((*header_ptr16) >> 10);
            uint32_t ones = ((*header_ptr16) & 0x1ff);
            uint32_t zeros = k_block_size - ones;

            // Number of runs <= 2.
            uint32_t special_bit = (((*header_ptr16) & 0x200) >> 9);
            if (!encoding_size) {
                uint32_t first_run_length = special_bit * ones + (1 - special_bit) * zeros;
                uint8_t inside_second_run = (first_run_length < local_i);
                return (inside_second_run ^ special_bit);
            }

            // Number of runs > 2.
            if (encoding_size < k_block_bytes) {
                if (std::min(ones, zeros) == encoding_size) {
                    // Minority encoding.
                    uint32_t tot = 0;
                    while (tot < encoding_size && *trunk_p < local_i) {
                        ++trunk_p;
                        ++tot;
                    }
                    uint8_t last_was_majority = ((!tot) || (*(trunk_p - 1) != local_i - 1));
                    return (last_was_majority ^ special_bit);
                }

                // Runs encoding.
                if (special_bit) {
                    uint32_t j = 0;
                    uint32_t acc = 0;

                    int32_t last = -1;
                    while (j + 1 < encoding_size && *(trunk_p + 1) < local_i) {
                        acc += *trunk_p - last; ++trunk_p;
                        last = *trunk_p; ++trunk_p; j += 2;
                    }

                    uint8_t access_i = 0;
                    if (j + 1 >= encoding_size) {
                        if (j < encoding_size) {  // j == encoding_size - 1
                            if (local_i <= (uint32_t)(*trunk_p) + 1) access_i = (((int32_t)local_i - last - 1) > 0);
                            else {
                                acc += (int32_t)(*trunk_p) - last;
                                if (ones - acc <= k_block_size - local_i) access_i = 0;
                                else access_i = 1;
                            }
                        } else {  // j == encoding_size
                            if ((int32_t)(ones - acc) < (int32_t)local_i - last - 1) access_i = 0;
                            else access_i = (((int32_t)local_i - last - 1) > 0);
                        }
                    } else {
                        if ((*trunk_p) < local_i - 1) access_i = 0;
                        else access_i = (((int32_t)local_i - last - 1) > 0);
                    }

                    return access_i;
                } else {
                    uint32_t j = 0;
                    uint32_t acc = 0;
                    int32_t last = -1;

                    while (j + 1 < encoding_size && *(trunk_p + 1) < local_i) {
                        acc += *trunk_p - last; ++trunk_p;
                        last = *trunk_p; ++trunk_p; j += 2;
                    }

                    uint8_t access_i = 0;
                    if (j + 1 >= encoding_size) {
                        if (j < encoding_size) {
                            if (local_i <= (uint32_t)(*trunk_p) + 1) access_i = (((int32_t)local_i - last - 1) == 0);
                            else {
                                acc += (*trunk_p) - last;
                                if (zeros - acc <= k_block_size - local_i) access_i = 1;
                                else access_i = 0;
                            }
                        } else {
                            if ((int32_t)(zeros - acc) < (int32_t)local_i - last - 1) access_i = 1;
                            else access_i = ((local_i - last - 1) == 0);
                        }
                    } else {
                        if ((*trunk_p) < local_i - 1) access_i = 1;
                        else access_i = (((int32_t)local_i - last - 1) == 0);
                    }

                    return access_i;
                }
            } else {
                // plain encoding.
                uint64_t* trunk_ptr64 = (uint64_t*)(((uint8_t*)m_trunk.data()) + trunk_ptr);
                uint32_t bit;
                for (bit = 0; bit + 64 <= local_i; bit += 64) trunk_ptr64++;

                uint8_t access_i = 0;
                if (bit != local_i) access_i = (((*trunk_ptr64) >> (local_i - bit - 1)) & 1);
                else access_i = (((*(trunk_ptr64 - 1)) >> 63) & 1);
                return access_i;
            }
        }

    public:
        //! Swap method
        void swap(hyb_vector& hybrid)
        {
            if (this != &hybrid) {
                std::swap(m_size, hybrid.m_size);
                std::swap(m_trunk, hybrid.m_trunk);
                std::swap(m_sblock_header, hybrid.m_sblock_header);
                std::swap(m_hblock_header, hybrid.m_hblock_header);
            }
        }

        //! Get the integer value of the binary string of length len starting at position idx.
        /*! \param idx Starting index of the binary representation of the integer.
         *  \param len Length of the binary representation of the integer. Default value is 64.
         *  \returns The integer value of the binary string of length len starting at position idx.
         *
         *  \pre idx+len-1 in [0..size()-1]
         *  \pre len in [1..64]
         */
        uint64_t get_int(size_type idx, const uint8_t len=64) const
        {
            uint64_t res = 0;
            for (size_t i=0; i<len; ++i) {
                res <<= 1;
                res |= (*this)[idx+len-1-i];
            }
            return res;
        }

        //! Accessing the i-th element of the original bitvector
        value_type operator[](size_type i) const
        {
            return access0(i + 1);
        }

        //! Assignment operator
        hyb_vector& operator=(const hyb_vector& hybrid)
        {
            if (this != &hybrid)
                copy(hybrid);
            return *this;
        }

        //! Move assignment operator
        hyb_vector& operator=(hyb_vector&& hybrid)
        {
            swap(hybrid);
            return *this;
        }

        //! Returns the size of the original bitvector
        size_type size() const
        {
            return m_size;
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const
        {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "size");
            written_bytes += m_trunk.serialize(out, child, "trunk");
            written_bytes += m_sblock_header.serialize(out, child, "sblock_header");
            written_bytes += m_hblock_header.serialize(out, child, "hblock_header");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream
        void load(std::istream& in)
        {
            read_member(m_size, in);
            m_trunk.load(in);
            m_sblock_header.load(in);
            m_hblock_header.load(in);
        }

        iterator begin() const
        {
            return iterator(this, 0);
        }

        iterator end() const
        {
            return iterator(this, size());
        }
};

template<uint32_t k_sblock_rate> const uint32_t hyb_vector<k_sblock_rate>::k_block_size = 256;
template<uint32_t k_sblock_rate> const uint32_t hyb_vector<k_sblock_rate>::k_block_bytes = 32;
template<uint32_t k_sblock_rate> const uint32_t hyb_vector<k_sblock_rate>::k_sblock_header_size = 8 + 2 * k_sblock_rate;
template<uint32_t k_sblock_rate> const uint32_t hyb_vector<k_sblock_rate>::k_sblock_size = 256 * k_sblock_rate;
template<uint32_t k_sblock_rate> const uint32_t hyb_vector<k_sblock_rate>::k_hblock_rate = (1U << 31) / 256;


template<uint8_t t_bp>
struct rank_result {
    typedef bit_vector::size_type size_type;
    static size_type adapt(size_type res, size_type)
    {
        return res;
    }
};

template<>
struct rank_result<0> {
    typedef bit_vector::size_type size_type;
    static size_type adapt(size_type res, size_type i)
    {
        return i-res;
    }
};

    template<uint8_t t_bp>
    struct succ_trait {
        typedef bit_vector::size_type size_type;
        static bool uniform_bit(size_type block_size, size_type ones, size_type)
        {
            return block_size == ones;
        }

        static bool empty_bit(size_type block_size, size_type, size_type zeros){
            return block_size == zeros;
        }
    };

    template<>
    struct succ_trait<0> {
        typedef bit_vector::size_type size_type;
        static size_type uniform_bit(size_type block_size, size_type, size_type zeros)
        {
            return block_size == zeros;
        }

        static bool empty_bit(size_type block_size, size_type ones, size_type){
            return block_size == ones;
        }
    };


//! Rank_support for the hyb_vector class
/*!
 * \tparam t_b      The bit pattern of size one. (so `0` or `1`)
 * \tparam k_sblock_rate  Superblock rate (number of blocks inside superblock)
 */
// TODO:
template<uint8_t t_b, uint32_t k_sblock_rate>
class rank_support_hyb
{
    public:
        typedef hyb_vector<k_sblock_rate> bit_vector_type;
        typedef typename bit_vector_type::size_type size_type;
        enum { bit_pat = t_b };
        enum { bit_pat_len = (uint8_t)1 };
    private:
        const bit_vector_type* m_v;


    public:
        //! Standard constructor
        explicit rank_support_hyb(const bit_vector_type* v = nullptr)
        {
            set_vector(v);
        }

        //! Answers rank queries
        const size_type rank(size_type i) const
        {
            assert(m_v != nullptr);
            assert(i <= m_v->size());

            // Handle easy case.
            if (i <= 0) return 0;

            size_type block_id = (i - 1) / bit_vector_type::k_block_size;
            size_type sblock_id = block_id / k_sblock_rate;
            size_type hblock_id = block_id / bit_vector_type::k_hblock_rate;

            size_type trunk_base = m_v->m_hblock_header[2 * hblock_id];
            size_type hblock_rank = m_v->m_hblock_header[2 * hblock_id + 1];

            uint32_t local_i = i - block_id * bit_vector_type::k_block_size;

            // Read superblock header.
            const uint8_t* header_ptr8 = ((const uint8_t*)(m_v->m_sblock_header.data())) + (sblock_id * bit_vector_type::k_sblock_header_size);
            uint32_t* header_ptr32 = (uint32_t*)header_ptr8;
            size_type trunk_ptr = trunk_base + ((*header_ptr32) & 0x3fffffff);
            size_type sblock_rank = *(header_ptr32 + 1);
            header_ptr8 += 8;

            uint16_t* header_ptr16 = (uint16_t*)header_ptr8;

            // Uniform superblock optimization.
            if ((*header_ptr32) & 0x80000000) {
                return rank_result<t_b>::adapt(
                        hblock_rank + sblock_rank +
                        ((*(header_ptr8 + 1)) & 0x01) * (i - sblock_id * bit_vector_type::k_sblock_size),
                        i);
            }

            // Fast forward through preceding blocks in the superblock.
            size_type block_rank = 0;
            for (size_type j = sblock_id * k_sblock_rate; j != block_id; ++j) {
                trunk_ptr += ((*header_ptr16) >> 10);     // Update trunk pointer.
                block_rank += ((*header_ptr16) & 0x1ff);  // Add 1s in the block.
                ++header_ptr16;
            }

            const uint8_t* trunk_p = ((uint8_t*)m_v->m_trunk.data()) + trunk_ptr;

            uint32_t encoding_size = ((*header_ptr16) >> 10);
            uint32_t ones = ((*header_ptr16) & 0x1ff);
            uint32_t zeros = bit_vector_type::k_block_size - ones;

            // Number of runs <= 2.
            uint32_t special_bit = (((*header_ptr16) & 0x200) >> 9);
            if (!encoding_size) {
                uint32_t first_run_length = special_bit * ones + (1 - special_bit) * zeros;
                uint32_t local_rank = std::min(local_i, first_run_length);

                return rank_result<t_b>::adapt(
                        hblock_rank + sblock_rank + block_rank +
                        (special_bit * local_rank + (1 - special_bit) * (local_i - local_rank)),
                        i);
            }

            // Number of runs > 2.
            if (encoding_size < bit_vector_type::k_block_bytes) {
                if (std::min(ones, zeros) == encoding_size) {
                    // Minority encoding.
                    uint32_t tot = 0;
                    while (tot < encoding_size && (*trunk_p++) < local_i) ++tot;

                    return rank_result<t_b>::adapt(
                            hblock_rank + sblock_rank + block_rank +
                            special_bit * tot + (1 - special_bit) * (local_i - tot),
                            i);
                }

                // Runs encoding.
                if (special_bit) {
                    uint32_t j = 0;
                    uint32_t acc = 0;

                    int32_t last = -1;
                    while (j + 1 < encoding_size && *(trunk_p + 1) < local_i) {
                        acc += *trunk_p - last; ++trunk_p;
                        last = *trunk_p; ++trunk_p; j += 2;
                    }

                    if (j + 1 >= encoding_size) {
                        if (j < encoding_size) {
                            if (*trunk_p  >= local_i) acc += local_i - last - 1;
                            else {
                                acc += (*trunk_p) - last;
                                acc += (ones - acc) - std::min(ones - acc, bit_vector_type::k_block_size - local_i);
                            }
                        } else acc += std::min(ones - acc, local_i - last - 1);
                    } else acc += std::min((int32_t)(*trunk_p), (int32_t)local_i - 1) - last;

                    return rank_result<t_b>::adapt(hblock_rank + sblock_rank + block_rank + acc,i);
                } else {
                    uint32_t j = 0;
                    uint32_t acc = 0;
                    int32_t last = -1;

                    while (j + 1 < encoding_size && *(trunk_p + 1) < local_i) {
                        acc += *trunk_p - last; ++trunk_p;
                        last = *trunk_p; ++trunk_p; j += 2;
                    }

                    if (j + 1 >= encoding_size) {
                        if (j < encoding_size) {
                            if (*trunk_p  >= local_i) acc += local_i - last - 1;
                            else {
                                acc += (*trunk_p) - last;
                                acc += (zeros - acc) - std::min(zeros - acc, bit_vector_type::k_block_size - local_i);
                            }
                        } else acc += std::min(zeros - acc, local_i - last - 1);
                    } else acc += std::min((int32_t)(*trunk_p), (int32_t)local_i - 1) - last;

                    return rank_result<t_b>::adapt(hblock_rank + sblock_rank + block_rank + (local_i - acc),i);
                }
            } else {
                // plain encoding.
                uint64_t* trunk_ptr64 = (uint64_t*)(((uint8_t*)m_v->m_trunk.data()) + trunk_ptr);
                uint32_t bit;
                for (bit = 0; bit + 64 <= local_i; bit += 64)
                    block_rank += bits::cnt(*trunk_ptr64++);
                if (bit != local_i)
                    block_rank += bits::cnt((*trunk_ptr64) & (((uint64_t)1 << (local_i - bit)) - 1));

                return rank_result<t_b>::adapt(hblock_rank + sblock_rank + block_rank, i);
            }
        }

        //! Shorthand for rank(i)
        const size_type operator()(size_type i) const
        {
            return rank(i);
        }

        //! Return the size of the original vector
        const size_type size() const
        {
            return m_v->size();
        }

        //! Set the supported vector
        void set_vector(const bit_vector_type* v = nullptr)
        {
            m_v = v;
        }

        //! Assignment operator
        rank_support_hyb& operator=(const rank_support_hyb& rs)
        {
            if (this != &rs) {
                set_vector(rs.m_v);
            }
            return *this;
        }

        //! Swap method
        void swap(rank_support_hyb&) {}

        //! Load the data structure from a stream and set the supported vector
        void load(std::istream&, const bit_vector_type* v = nullptr)
        {
            set_vector(v);
        }

        //! Serializes the data structure into a stream
        size_type serialize(std::ostream&, structure_tree_node* v = nullptr, std::string name = "") const
        {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            structure_tree::add_size(child, 0);
            return 0;
        }
};

//! Select support for the hyb_vector class
/*!
 * \tparam t_b            The bit pattern of size one. (so `0` or `1`)
 * \tparam k_sblock_rate  Superblock rate (number of blocks inside superblock)
 * TODO: implement select queries, currently this is dummy class.
 */
template<uint8_t t_b, uint32_t k_sblock_rate>
class select_support_hyb
{
    public:
        typedef hyb_vector<k_sblock_rate> bit_vector_type;
        typedef rank_support_hyb<t_b, k_sblock_rate> rank_type;
        typedef typename bit_vector_type::size_type size_type;
        enum { bit_pat = t_b };
        enum { bit_pat_len = (uint8_t)1 };
    private:
        const rank_type* m_rank;

    public:
        //! Standard constructor
        explicit select_support_hyb(const bit_vector_type* v = nullptr)
        { }

        explicit select_support_hyb(const rank_type* r){
            m_rank = r;
        }

        //! Answers select queries
        //! Answers select queries
        size_type select(size_type i) const
        {
            //fprintf(stderr, "\nzombit_vector: select queries are not currently supported\n");
            //std::exit(EXIT_FAILURE);

            uint64_t l = 0, r = m_rank->size()-1;
            uint64_t mid, cnt;
            while(l < r){
                mid = (l + r) >> 1;
                cnt = m_rank->rank(mid+1);
                if(cnt < i){
                    l = mid + 1;
                }else{
                    r = mid;
                }
            }
            return l;
        }

        //! Shorthand for select(i)
        size_type operator()(size_type i) const
        {
            return select(i);
        }

        //! Return the size of the original vector
        size_type size() const
        {
            return m_rank->size();
        }

        //! Assignment operator
        select_support_hyb& operator=(const select_support_hyb& rs)
        {
            if (this != &rs) {
                m_rank = rs.m_rank;
            }
            return *this;
        }

        //! Swap method
        void swap(select_support_hyb&) {}

        //! Load the data structure from a stream and set the supported vector
        void load(std::istream&, const rank_type* rank)
        {
            m_rank = rank;
        }

        //! Load the data structure from a stream and set the supported vector
        void load(std::istream&, const bit_vector_type* bv)
        {
        }

        void set_vector(const bit_vector_type* bv){

        }

        //! Serializes the data structure into a stream
        size_type serialize(std::ostream&, sdsl::structure_tree_node* v = nullptr, std::string name = "") const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            sdsl::structure_tree::add_size(child, 0);
            return 0;
        }
};

    template<uint8_t t_b, uint32_t k_sblock_rate>
    class succ_support_hyb
    {
    public:
        typedef hyb_vector<k_sblock_rate> bit_vector_type;
        typedef typename bit_vector_type::size_type size_type;
        enum { bit_pat = t_b };
        enum { bit_pat_len = (uint8_t)1 };
    private:
        const bit_vector_type* m_v;

        size_type local_succ(size_type local_i, size_type block_id, const uint8_t* trunk_ptr, uint64_t* trunk_ptr64, uint16_t* header_ptr16) const{

            uint32_t encoding_size = ((*header_ptr16) >> 10);
            uint32_t ones = ((*header_ptr16) & 0x1ff);
            uint32_t zeros = bit_vector_type::k_block_size - ones;

            if(succ_trait<t_b>::uniform_bit(bit_vector_type::k_block_size, ones,  zeros)){
                return local_i;
            }
            if(succ_trait<t_b>::empty_bit(bit_vector_type::k_block_size, ones, zeros)){
                return bit_vector_type::k_block_size;
            }

            // Number of runs <= 2.
            uint32_t special_bit = (((*header_ptr16) & 0x200) >> 9);
            if (!encoding_size) {
                uint32_t first_run_length = special_bit * ones + (1 - special_bit) * zeros;
                if((uint8_t) (special_bit) == t_b){ //Solution in the first run
                    if(local_i < first_run_length) {
                        return local_i;
                    }else{
                        return bit_vector_type::k_block_size;
                    }
                }else{ //Solution in the second run
                    return std::max(local_i, (size_type) first_run_length);
                }
            }

            // Number of runs > 2.
            if (encoding_size < bit_vector_type::k_block_bytes) {
                if (std::min(ones, zeros) == encoding_size) {
                    // Minority encoding.
                    uint8_t keybit = (ones < zeros);
                    uint32_t tot = 0;
                    //Looking for the first keybit >= local_i
                    while (tot < encoding_size && *trunk_ptr < local_i) {
                        ++trunk_ptr;
                        ++tot;
                    }
                    if(keybit == t_b){ //Easy part
                        if(tot == encoding_size) return bit_vector_type::k_block_size;
                        return *trunk_ptr;
                    }else{ //Check the !t_b positions
                        if(tot == encoding_size || local_i != *trunk_ptr) return local_i;
                        while(tot < encoding_size-1 && (*(trunk_ptr+1) - *trunk_ptr) == 1){
                            ++trunk_ptr;
                            ++tot;
                        }
                        return *trunk_ptr+1; //TODO: check creo que se non atopa debería devolver o k_block_size
                    }
                }else{
                    // Runs encoding.

                    if(special_bit == t_b){
                        uint32_t j = 0;
                        uint32_t acc_ones = 0;
                        int32_t last = -1;

                        while (j + 1 < encoding_size && *(trunk_ptr + 1) < local_i) {
                            acc_ones += *trunk_ptr - last; ++trunk_ptr;
                            last = *trunk_ptr; ++trunk_ptr; j += 2;
                        }
                        //The first part is 111111 and the second 00000
                        if(j+1 >= encoding_size){
                            if(j < encoding_size) { // j == encoding_size-1
                                if(local_i <= *(trunk_ptr)){
                                    return local_i;
                                }else{
                                    acc_ones += (int32_t) *trunk_ptr - last;
                                    auto first_one = bit_vector_type::k_block_size - (ones-acc_ones);
                                    if (local_i < first_one) {
                                        return first_one;
                                    }else{
                                        return local_i;
                                    }
                                }
                            }else{ // j == encoding_size
                                auto last_one = last + (ones - acc_ones);
                                if(local_i <= last_one){
                                    return local_i;
                                }else{
                                    return bit_vector_type::k_block_size;
                                }
                            }
                        }else{
                            if(local_i <= *trunk_ptr){
                                return local_i;
                            }else{
                                return *(trunk_ptr+1)+1;
                            }
                        }
                    }else{
                        uint32_t j = 0;
                        uint32_t acc = 0;
                        int32_t last = -1;

                        while (j + 1 < encoding_size && *(trunk_ptr + 1) < local_i) {
                            acc += *trunk_ptr - last; ++trunk_ptr;
                            last = *trunk_ptr; ++trunk_ptr; j += 2;
                        }

                        //The first part is 00000 and the second 11111
                        if(j + 1 >= encoding_size){
                            if (j < encoding_size) {  // j == encoding_size - 1
                                acc += (int32_t) *trunk_ptr - last;
                                auto last_one = bit_vector_type::k_block_size - (zeros - acc)-1;
                                if(local_i > *trunk_ptr && local_i <= last_one){
                                    return local_i;
                                }else if (local_i <= last_one){
                                    return *(trunk_ptr)+1;
                                }else{
                                    return bit_vector_type::k_block_size;
                                }
                            }else{
                                auto first_one = last + (zeros-acc)+1;
                                if (local_i > first_one) {
                                    return local_i;
                                } else {
                                    return first_one;
                                }
                            }
                        }else{
                            if(local_i > *trunk_ptr){
                                return local_i;
                            }else{
                                return *trunk_ptr+1;
                            }
                        }
                    }
                }
            } else {
                // plain encoding.
                return bits::next(trunk_ptr64, local_i);
            }
        }

        size_type local_succ2(size_type local_i, size_type block_id, const uint8_t* trunk_ptr, uint64_t* trunk_ptr64, uint16_t* header_ptr16) const{

            uint32_t encoding_size = ((*header_ptr16) >> 10);
            uint32_t ones = ((*header_ptr16) & 0x1ff);
            uint32_t zeros = bit_vector_type::k_block_size - ones;

            if(succ_trait<t_b>::uniform_bit(bit_vector_type::k_block_size, ones,  zeros)){
                return local_i;
            }
            if(succ_trait<t_b>::empty_bit(bit_vector_type::k_block_size, ones, zeros)){
                return bit_vector_type::k_block_size;
            }

            // Number of runs <= 2.
            uint32_t special_bit = (((*header_ptr16) & 0x200) >> 9);
            if (!encoding_size) {
                uint32_t tot = 0;
                uint32_t first_run_length = special_bit * ones + (1 - special_bit) * zeros;
                if((uint8_t) (special_bit) == t_b){ //Solution in the first run
                    if(local_i < first_run_length) {
                        return local_i;
                    }else{
                        return bit_vector_type::k_block_size;
                    }
                }else{ //Solution in the second run
                    return std::max(local_i, (size_type) first_run_length);
                }
            }

            if (encoding_size < bit_vector_type::k_block_bytes) {
                if (std::min(ones, zeros) == encoding_size) {
                    // Minority encoding.
                    uint8_t keybit = (ones < zeros);
                    uint32_t tot = 0;
                    //Looking for the first t_b greater than local_i
                    while (tot < encoding_size && (*trunk_ptr) < local_i) {
                        ++trunk_ptr;
                        ++tot;
                    }
                    if(keybit == t_b){ //Easy part
                        if(tot == encoding_size) return bit_vector_type::k_block_size;
                        return *trunk_ptr;
                    }else{ //Check the !t_b positions
                        size_type interval_beg = tot ? *(trunk_ptr-1)+1 : 0;
                        size_type interval_end = *trunk_ptr-1;
                        //Check if local_i is contained into the first !t_b interval
                        if(interval_beg <= local_i && local_i <= interval_end){
                            return local_i;
                        }
                        //Continue with the others !t_b intervals
                        ++trunk_ptr;
                        ++tot;
                        while (tot < encoding_size) {
                            //Check if exists an interval of !t_b which ends after local_i
                            interval_beg = *(trunk_ptr-1)+1;
                            interval_end = *trunk_ptr-1;
                            if(interval_beg <= interval_end){
                                return interval_beg;
                            }
                            ++trunk_ptr;
                            ++tot;
                        }
                        //Last interval
                        if(*(trunk_ptr-1) < bit_vector_type::k_block_size-1){
                            return *(trunk_ptr-1)+1;
                        }
                        return bit_vector_type::k_block_size;
                    }
                }else{
                    if(special_bit == t_b){
                        uint32_t j = 0;
                        uint32_t acc_ones = 0;
                        int32_t last = -1;

                        while (j + 1 < encoding_size && *(trunk_ptr + 1) < local_i) {
                            acc_ones += *trunk_ptr - last; ++trunk_ptr;
                            last = *trunk_ptr; ++trunk_ptr; j += 2;
                        }
                        //The first part is 111111 and the second 00000
                        if(j+1 >= encoding_size){
                            if(j < encoding_size) { // j == encoding_size-1
                                if(local_i <= *(trunk_ptr)){
                                    return local_i;
                                }else{
                                    acc_ones += (int32_t) *trunk_ptr - last;
                                    auto first_one = bit_vector_type::k_block_size - (ones-acc_ones);
                                    if (local_i < first_one) {
                                        return first_one;
                                    }else{
                                        return local_i;
                                    }
                                }
                            }else{ // j == encoding_size
                                auto last_one = last + (ones - acc_ones);
                                if(local_i <= last_one){
                                    return local_i;
                                }else{
                                    return bit_vector_type::k_block_size;
                                }
                            }
                        }else{
                            if(local_i <= *trunk_ptr){
                                return local_i;
                            }else{
                                return *(trunk_ptr+1)+1;
                            }
                        }
                    }else{
                        uint32_t j = 0;
                        uint32_t acc = 0;
                        int32_t last = -1;

                        while (j + 1 < encoding_size && *(trunk_ptr + 1) < local_i) {
                            acc += *trunk_ptr - last; ++trunk_ptr;
                            last = *trunk_ptr; ++trunk_ptr; j += 2;
                        }

                        //The first part is 00000 and the second 11111
                        if(j + 1 >= encoding_size){
                            if (j < encoding_size) {  // j == encoding_size - 1
                                acc += (int32_t) *trunk_ptr - last;
                                auto last_one = bit_vector_type::k_block_size - (zeros - acc)-1;
                                if(local_i > *trunk_ptr && local_i <= last_one){
                                    return local_i;
                                }else if (local_i <= last_one){
                                    return *(trunk_ptr)+1;
                                }else{
                                    return bit_vector_type::k_block_size;
                                }
                            }else{
                                auto first_one = last + (zeros-acc)+1;
                                if (local_i > first_one) {
                                    return local_i;
                                } else {
                                    return first_one;
                                }
                            }
                        }else{
                            if(local_i > *trunk_ptr){
                                return local_i;
                            }else{
                                return *trunk_ptr+1;
                            }
                        }
                    }
                }
            } else {
                // plain encoding.
                return bits::next(trunk_ptr64, local_i);
                //if(next >= bit_vector_type::k_block_size) return -1;
                //return next;
            }
        }

        size_type hb_leftmost_greater(size_type beg, size_type end, size_type hblock_rank) const{
            // auto i = sblock_id_beg;
            // auto j = sblock_id_end+1;

            size_type diff = 1, rank;
            while(beg + diff <= end){
                rank = m_v->m_hblock_header[2 * (beg+diff) + 1];
                if(rank > hblock_rank) break;
                diff = diff << 1;
            }

            size_type hb_end = std::min(beg + diff, end);
            size_type hb_beg = beg + (diff >> 1) + 1;
            size_type idx, r = hb_beg;
            while(hb_beg <= hb_end){
                idx = (hb_beg+hb_end)>>1;
                rank = m_v->m_hblock_header[2 * idx + 1];
                if (rank > hblock_rank){
                    r = idx;
                    hb_end = idx-1;
                }else {
                    hb_beg = idx+1;
                }
            }
            return r;
        }

        size_type sb_leftmost_greater(size_type beg, size_type end, size_type sblock_rank) const{
            size_type rank;
            size_type diff = 2, sb_beg = beg+1, sb_end= beg+diff;
            while(sb_end <= end){
                const uint8_t* header_ptr8 = ((const uint8_t*)(m_v->m_sblock_header.data())) +
                                             (sb_end * bit_vector_type::k_sblock_header_size);
                uint32_t* header_ptr32 = (uint32_t*)header_ptr8;
                rank = *(header_ptr32 + 1);
                if(rank > sblock_rank) break;
                diff = diff << 1;
                sb_beg = sb_end;
                sb_end = beg + diff;
            }
            size_type idx;
            size_type r = end; //se non hai, pode haber no ultimo
            sb_end = std::min(sb_end, end);
            while(sb_beg <= sb_end){
                idx = (sb_beg+sb_end)>>1;
                const uint8_t* header_ptr8 = ((const uint8_t*)(m_v->m_sblock_header.data()))
                                             + (idx * bit_vector_type::k_sblock_header_size);
                uint32_t* header_ptr32 = (uint32_t*)header_ptr8;
                rank = *(header_ptr32 + 1);
                if (rank > sblock_rank){
                    r = idx-1;
                    sb_end = idx-1;
                }else {
                    sb_beg = idx+1;
                }
            }
            return r;
        }



    public:
        //! Standard constructor
        explicit succ_support_hyb(const bit_vector_type* v = nullptr)
        {
            set_vector(v);
        }

        //! Answers select queries
        size_type succ(size_type i) const
        {
            size_type n_blocks = (m_v->m_size + bit_vector_type::k_block_size - 1) / bit_vector_type::k_block_size;
            size_type n_sblocks = (n_blocks + k_sblock_rate - 1) / k_sblock_rate;
            size_type n_hblocks = (n_blocks + bit_vector_type::k_hblock_rate - 1) / bit_vector_type::k_hblock_rate;

            size_type block_id = i / bit_vector_type::k_block_size;
            size_type sblock_id = block_id / k_sblock_rate;
            size_type hblock_id = block_id / bit_vector_type::k_hblock_rate;

            size_type trunk_base = m_v->m_hblock_header[2 * hblock_id];
            size_type hblock_rank = m_v->m_hblock_header[2 * hblock_id + 1];

            uint32_t off = i - block_id * bit_vector_type::k_block_size;
            // Read superblock header.
            const uint8_t* header_ptr8 = ((const uint8_t*)(m_v->m_sblock_header.data())) + (sblock_id * bit_vector_type::k_sblock_header_size);
            uint32_t* header_ptr32 = (uint32_t*)header_ptr8;
            size_type trunk_ptr = trunk_base + ((*header_ptr32) & 0x3fffffff);
            size_type sblock_rank = *(header_ptr32 + 1);
            header_ptr8 += 8;
            uint16_t* header_ptr16 = (uint16_t*)header_ptr8;

            // Uniform superblock optimization.
            if (((*header_ptr32) & 0x80000000) && ((*(header_ptr8 + 1)) & 0x01) == t_b) {
                return i;
            }

            // Fast-forward through preceding blocks in the superblock.
            for (size_type j = sblock_id * k_sblock_rate; j != block_id; ++j) {
                trunk_ptr += ((*header_ptr16) >> 10);     // Update trunk pointer.
                ++header_ptr16;
            }
            const uint8_t* trunk_p = ((const uint8_t*) m_v->m_trunk.data()) + trunk_ptr;
            uint64_t* trunk_ptr64 = (uint64_t*)(((uint8_t*) m_v->m_trunk.data()) + trunk_ptr);

            size_type local_i = local_succ(off, block_id, trunk_p, trunk_ptr64, header_ptr16);
            if(local_i < bit_vector_type::k_block_size){
                return block_id * bit_vector_type::k_block_size + local_i;
            }else{
                //Not found. Try with the rest of blocks within the superblock
                off = 0;
                trunk_ptr += ((*header_ptr16) >> 10);     // Update trunk pointer.
                ++header_ptr16;
                ++block_id;
                while (block_id < (sblock_id +1)* k_sblock_rate) {
                    //while (block_id+1 < n_blocks) {
                    if(((*header_ptr16) & 0x1ff) > 0){
                        trunk_p = ((const uint8_t*) m_v->m_trunk.data()) + trunk_ptr;
                        trunk_ptr64 = (uint64_t*)(((uint8_t*) m_v->m_trunk.data()) + trunk_ptr);
                        return block_id * bit_vector_type::k_block_size + local_succ(off, block_id, trunk_p, trunk_ptr64, header_ptr16);
                    }
                    trunk_ptr += ((*header_ptr16) >> 10);     // Update trunk pointer.
                    ++header_ptr16;
                    ++block_id;

                }
                //Not found. Try with the rest of superblocks
                //sblock_id = upper_bound(sblock_id+1, n_sblocks-1, aux_sblock_rank);
                if(sblock_id+1 == n_sblocks) return m_v->size();

                //Getting number of ones up contained within sblock_id (stored at sblock_id+1)
                header_ptr8 = ((const uint8_t*)(m_v->m_sblock_header.data())) + ((sblock_id+1) * bit_vector_type::k_sblock_header_size);
                header_ptr32 = (uint32_t*)header_ptr8;
                sblock_rank = *(header_ptr32 + 1);
                //Getting the last sblock in the current hblock
                size_type end_sb = std::min((hblock_id+1)* (bit_vector_type::k_hblock_rate/k_sblock_rate)-1, n_sblocks-1);

                sblock_id = sb_leftmost_greater(sblock_id, end_sb, sblock_rank);
                block_id = sblock_id * k_sblock_rate;
                header_ptr8 = ((const uint8_t*)(m_v->m_sblock_header.data())) + (sblock_id * bit_vector_type::k_sblock_header_size);
                header_ptr32 = (uint32_t*)header_ptr8;
               // hblock_id = block_id / bit_vector_type::k_hblock_rate;
               // trunk_base = m_v->m_hblock_header[2 * hblock_id];
                trunk_ptr = trunk_base + ((*header_ptr32) & 0x3fffffff);
                header_ptr8 += 8;
                header_ptr16 = (uint16_t*)header_ptr8;
                //Found inside a superblock. Look for the correct block
                while (block_id < (sblock_id +1)* k_sblock_rate) {
                    if(((*header_ptr16) & 0x1ff) > 0){
                        trunk_p = ((const uint8_t*) m_v->m_trunk.data()) + trunk_ptr;
                        trunk_ptr64 = (uint64_t*)(((uint8_t*) m_v->m_trunk.data()) + trunk_ptr);
                        return block_id * bit_vector_type::k_block_size + local_succ(off, block_id, trunk_p, trunk_ptr64, header_ptr16);
                    }
                    trunk_ptr += ((*header_ptr16) >> 10);     // Update trunk pointer.
                    ++header_ptr16;
                    ++block_id;
                }
                //Not found. Try with the rest of hyperblocks
                if(hblock_id + 1 == n_hblocks) return m_v->size();
                hblock_rank = m_v->m_hblock_header[2 * (hblock_id+1) + 1];
                hblock_id = hb_leftmost_greater(hblock_id, n_hblocks-1, hblock_rank);
                if(hblock_id == n_hblocks) return m_v->size();
                sblock_id = (hblock_id * bit_vector_type::k_hblock_rate) / k_sblock_rate;
                end_sb = std::min((hblock_id+1)* (bit_vector_type::k_hblock_rate/k_sblock_rate)-1, n_sblocks-1);
                sblock_id = sb_leftmost_greater(sblock_id, end_sb, 0);
                block_id = sblock_id * k_sblock_rate;
                header_ptr8 = ((const uint8_t*)(m_v->m_sblock_header.data())) + (sblock_id * bit_vector_type::k_sblock_header_size);
                header_ptr32 = (uint32_t*)header_ptr8;
                trunk_base = m_v->m_hblock_header[2 * hblock_id];
                trunk_ptr = trunk_base + ((*header_ptr32) & 0x3fffffff);
                header_ptr8 += 8;
                header_ptr16 = (uint16_t*) header_ptr8;
                //Found inside a superblock. Look for the correct block
                while (block_id < (sblock_id +1)* k_sblock_rate) {
                    if(((*header_ptr16) & 0x1ff) > 0){
                        trunk_p = ((const uint8_t*) m_v->m_trunk.data()) + trunk_ptr;
                        trunk_ptr64 = (uint64_t*)(((uint8_t*) m_v->m_trunk.data()) + trunk_ptr);
                        return block_id * bit_vector_type::k_block_size + local_succ(off, block_id, trunk_p, trunk_ptr64, header_ptr16);
                    }
                    trunk_ptr += ((*header_ptr16) >> 10);     // Update trunk pointer.
                    ++header_ptr16;
                    ++block_id;
                }

            }
            return m_v->size();


        }

        //! Shorthand for select(i)
        const size_type operator()(size_type i) const
        {
            return succ(i);
        }

        //! Return the size of the original vector
        const size_type size() const
        {
            return m_v->size();
        }

        //! Set the supported vector
        void set_vector(const bit_vector_type* v = nullptr)
        {
            m_v = v;
        }

        //! Assignment operator
        succ_support_hyb& operator=(const succ_support_hyb& rs)
        {
            if (this != &rs) {
                set_vector(rs.m_v);
            }
            return *this;
        }

        //! Swap method
        void swap(succ_support_hyb&) {}

        //! Load the data structure from a stream and set the supported vector
        void load(std::istream&, const bit_vector_type* v = nullptr)
        {
            set_vector(v);
        }

        //! Serializes the data structure into a stream
        size_type serialize(std::ostream&, structure_tree_node* v = nullptr, std::string name = "") const
        {
            structure_tree_node* child = structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            structure_tree::add_size(child, 0);
            return 0;
        }
    };

}  // end namespace sdsl

#endif  // INCLUDED_SDSL_HYB_VECTOR
