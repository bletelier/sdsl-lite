
/** 
    A class for arrays that can be initialized in constant time (total!)
    Copyright (C) 2021 Diego Arroyuelo

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
**/


#ifndef INIT_ARRAY
#define INIT_ARRAY

#include <cstdint>

using namespace std;


namespace sdsl {

    template<typename array_t>
    class initializable_array_classic
    {
        uint64_t D; // array size
        array_t init_value;  // value used to initialize the array
        std::vector<array_t> V;  // array
        std::vector<array_t> S;
        std::vector<array_t> U; // auxiliary array

        uint64_t top;

        void copy(const initializable_array_classic &o){
            D = o.D;
            init_value = o.init_value;
            V = o.V;
            S = o.S;
            U = o.U;
        }

    public:

        initializable_array_classic() = default;

        initializable_array_classic(const uint64_t _D, array_t _init_value)
        {
            D = _D;
            top = 0;
            V.reserve(D);
            S.reserve(D);
            U.reserve(D);
            init_value = _init_value;
        };

        void init(const uint64_t _D, array_t _init_value)
        {
            D = _D;
            top = 0;
            V.reserve(D);
            S.reserve(D);
            U.reserve(D);
            init_value = _init_value;
        };

        //! Copy constructor
        initializable_array_classic(const initializable_array_classic& o)
        {
            copy(o);
        }

        //! Move constructor
        initializable_array_classic(initializable_array_classic&& o)
        {
            *this = std::move(o);
        }


        initializable_array_classic &operator=(const initializable_array_classic &o) {
            if (this != &o) {
                copy(o);
            }
            return *this;
        }

        initializable_array_classic &operator=(initializable_array_classic &&o) {
            if (this != &o) {
                D = o.D;
                init_value = o.init_value;
                V = std::move(o.V);
                S = std::move(o.S);
                U = std::move(o.U);
            }
            return *this;
        }

        void swap(initializable_array_classic &o) {
            // m_bp.swap(bp_support.m_bp); use set_vector to set the supported bit_vector
            std::swap(D, o.D);
            std::swap(init_value, o.init_value);
            std::swap(V, o.V);
            std::swap(S, o.S);
            std::swap(U, o.U);
        }


        array_t operator[](const uint64_t i) const
        {
            if (i >= D) {
                cerr << "Caution init_array.hpp: index " << i << " is greater than array size " << D << endl;
            }
            //std::cout << "U[i] = " << U[i] << " top = " << top << " S[U[i]] = " << S[U[i]] << " i = " << i << std::endl;
            if (U[i] < top && S[U[i]] == i)
                return V[i];
            else
                return init_value;
        };

        array_t atPos(const uint64_t i) const
        {
            if (i >= D) {
                cerr << "Caution init_array.hpp function atPos: index " << i << " is greater than array size " << D << endl;
            }

            //std::cout << "U[i] = " << U[i] << " top = " << top << " S[U[i]] = " << S[U[i]] << " i = " << i << std::endl;
            if (U[i] < top && S[U[i]] == i)
                return V[i];
            else
                return init_value;
        };

        array_t& operator[](uint64_t i)
        {
            if (i >= D) {
                cerr << "Caution init_array.hpp operator array_t& [uint64_t i]: index " << i << " is greater than array size " << D << endl;
            }

            //std::cout << "Oh! " << i << std::endl;
            if (U[i] < top && S[U[i]] == i)
                return V[i];
            else {
                U[i] = top;
                S[top++] = i;
                return V[i];
            }
        };



        uint64_t size()
        {
            return D;
        };

        uint64_t size_in_bytes()
        {
            return sizeof(array_t)*D + sizeof(array_t) + 3*sizeof(uint64_t)*D + sizeof(uint64_t) + 3*sizeof(uint64_t*);
        };

    };


    template<typename array_t>
    class initializable_array
    {
        uint64_t n; // array size
        array_t init_value;  // value used to initialize the array

        std::vector<array_t> A;  // array
        initializable_array_classic<uint64_t> W;


        void copy(const initializable_array &o){
            n = o.n;
            init_value = o.init_value;
            A = o.A;
            W = o.W;
        }


    public:

        initializable_array() = default;

        initializable_array(const uint64_t _n, array_t _init_value)
        {
            n = _n;
            A.reserve(n);
            init_value = _init_value;
            W.init((n+63)/64, 0); // = initializable_array_classic<uint64_t>((n+63)/64 + 1, 0);
        };

        array_t operator[](const uint64_t i) const
        {
            //if (i >= D) {
            //    cerr << "Caution init_array.hpp: index " << i << " is greater than array size " << D << endl;
            //}
            //std::cout << "U[i] = " << U[i] << " top = " << top << " S[U[i]] = " << S[U[i]] << " i = " << i << std::endl;
            //cout << "[ ]" << i << " " << n <<std::endl;
            //register uint64_t temp = W.atPos(i/64);
            /*if (!temp)
                return init_value;
            else {*/
            if ((W.atPos(i/64)>>(i%64/*&0x3F*/)) & 1ull) {
                return A[i];
            }
            else
                return init_value;
            /*}*/
        };

        //! Copy constructor
        initializable_array(const initializable_array& o)
        {
            copy(o);
        }

        //! Move constructor
        initializable_array(initializable_array&& o)
        {
            *this = std::move(o);
        }


        initializable_array &operator=(const initializable_array &o) {
            if (this != &o) {
                copy(o);
            }
            return *this;
        }

        initializable_array &operator=(initializable_array &&o) {
            if (this != &o) {
                n = o.n;
                init_value = o.init_value;
                A = std::move(o.A);
                W = std::move(o.W);
            }
            return *this;
        }

        void swap(initializable_array &o) {
            std::swap(n, o.n);
            std::swap(init_value, o.init_value);
            std::swap(A, o.A);
            std::swap(W, o.W);
        }

        array_t atPos(const uint64_t i) const
        {
            //if (i >= D) {
            //    cerr << "Caution init_array.hpp function atPos: index " << i << " is greater than array size " << D << endl;
            //}

            //std::cout << "U[i] = " << U[i] << " top = " << top << " S[U[i]] = " << S[U[i]] << " i = " << i << std::endl;
            //cout << "[ ]" << i << " " << n <<std::endl;
            //register uint64_t temp = W.atPos(i/64);
            /*if (!temp)
                return init_value;
            else {*/
            if ((W.atPos(i/64)>>(i%64/*&0x3F*/)) & 1ull) {
                return A[i];
            }
            else
                return init_value;
            //}
        };

        array_t& operator[](uint64_t i)
        {
            //if (i >= D) {
            //    cerr << "Caution init_array.hpp operator array_t& [uint64_t i]: index " << i << " is greater than array size " << D << endl;
            //}

            //std::cout << "[ ]" << i << " " << n <<std::endl;
            register uint64_t temp = W.atPos(i/64);
            if ((/*W.atPos(i/64)*/temp>>(i%64/*&0x3F*/)) & 1ull) {
                return A[i];
            }
            else {
                W[i/64] = /*W.atPos(i/64)*/ temp | (1ull<<(i%64/*&0x3F*/));
                return A[i];
            }
        };

        // number of elements in the array
        uint64_t size()
        {
            return n;
        };

        uint64_t size_in_bytes()
        {
            return sizeof(array_t)*n + W.size_in_bytes() + sizeof(array_t) + sizeof(array_t *) + sizeof(uint64_t);
        };

    };

}

#endif
