
#ifndef HTC_HPP
#define HTC_HPP

#include <cstdint>
#include <iostream>
#include <omp.h>
#include <cstdio>
#include <sstream>
#include <fstream>

#include "mkmh.hpp"


namespace mkmh{
    class HASHTCounter{

    public:
        HASHTCounter();
        HASHTCounter(uint64_t sz);
        ~HASHTCounter();
        int& operator[](hash_t key);

        void increment(hash_t key);
        void bulk_increment(hash_t* h, int num);

        int& get(hash_t key);

        void get(hash_t key, int& ret);

        int size(void);
        void size(int sz);
        void resize(int sz);
        inline void set(int pos, int val){
            *(counts + pos) = val;
        };

        int* begin(void);

        std::string to_string();
        void write_to_binary(std::string filename);
        void print();
        
    private:
        uint64_t my_size;
        int* counts;
        
};
}

#endif
