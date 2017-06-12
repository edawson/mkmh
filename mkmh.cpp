#include "mkmh.hpp"
#include <string.h>

namespace mkmh{

    using namespace std;

    unordered_map<char, char> RCOMP (
            {
            {'A', 'T'},
            {'a', 'T'},
            {'T','A'},
            {'t', 'a'},
            {'C', 'G'},
            {'c', 'G'},
            {'G', 'C'},
            {'g', 'C'}

            });

    /**
     * TODO: move to_upper and rev_comp and reverse to C-style
     * rip out the rev_comp_helper function
     * move chars to integers for speed
     */
        char rev_arr [26] = {84, 66, 71, 68, 69, 70, 67, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 65,
                       85, 86, 87, 88, 89, 90};
     /*
     * char[26] up_arr = {};
     */
    void reverse_reverse_complement(const char* seq, char* ret, int len){
        for (int i = len - 1; i >=0; i--){
            ret[ len - 1 - i ] = (char) rev_arr[ (int) seq[i] - 65];
        }
    }

    string reverse_complement(string seq){
        stringstream ret;
/**
        for (int i = 0; i < seq.size(); i++){
            ret << RCOMP[seq[i]];

        }
        return ret.str();

        **/

        for (int i = 0; i < seq.length(); i++){
            char c = seq[i];
            switch (c){
                case 'A':
                    ret << "T";
                    break;
                case 'a':
                    ret << "t";
                    break;
                case 'T':
                    ret << "A";
                    break;
                case 't':
                    ret << "a";
                    break;
                case 'C':
                    ret << "G";
                    break;
                case 'c':
                    ret << "g";
                    break;
                case 'G':
                    ret << "C";
                    break;
                case 'g':
                    ret << "c";
                    break;
                     //Handle X, N, Y, all that stuff. 
                default:
                    ret << c;
                    break;
            }
        }
        return ret.str();

    }

    string reverse(string seq){
        string copy = string(seq);
        std::reverse(copy.begin(), copy.end());
        return copy;
    }

    void to_upper(char* seq, int length){
        
        for (int i = 0; i < length; i++){
            char c = seq[i];
            seq[i] = ( (c - 91) > 0 ? c - 32 : c);
        }
    }

    /**
    void to_upper(char* seq, int length){

        for (int i = 0; i < length; i++){
            char c = seq[i];
            switch (c){
                case 'A':
                    seq[i] = 'A';
                    break;
                case 'a':
                    seq[i] = 'A';
                    break;
                case 'T':
                    seq[i] = 'T';
                    break;
                case 't':
                    seq[i] = 'T';
                    break;
                case 'C':
                    seq[i] = 'C';
                    break;
                case 'c':
                    seq[i] = 'C';
                    break;
                case 'G':
                    seq[i] = 'G';
                    break;
                case 'g':
                    seq[i] = 'G';
                    break;
                    // Handle X, N, Y, all that stuff.
                default:
                    seq[i] = seq[i];
                    break;
            }
        }

    }
**/
    string to_upper(string seq){
        stringstream ret;
        std::locale loc;
        //for (int i = 0; i < seq.length(); i++){
        //    ret << std::toupper(seq[i], loc);
        //}
        for (int i = 0; i < seq.length(); i++){
            char c = seq[i];
            switch (c){
                case 'A':
                    ret << "A";
                    break;
                case 'a':
                    ret << "A";
                    break;
                case 'T':
                    ret << "T";
                    break;
                case 't':
                    ret << "T";
                    break;
                case 'C':
                    ret << "C";
                    break;
                case 'c':
                    ret << "C";
                    break;
                case 'G':
                    ret << "G";
                    break;
                case 'g':
                    ret << "G";
                    break;
                    /* Handle X, N, Y, all that stuff. */
                default:
                    ret << std::toupper(c, loc);
                    break;
            }
        }

        //std::transform(seq.begin(), str.end(), str.begin(), ::to_upper);
        return ret.str();
    }

    vector<string> kmer_set(vector<string> kmers){
        set<string> uniqs = set<string> (kmers.begin(), kmers.end());
        vector<string> ret = vector<string> (uniqs.begin(), uniqs.end());
        return ret;
    }

    hash_t kmer_to_integer(string kmer){
        //00 - 'G', 01 - 'A', 10 - 'T', and 11 - 'C'
        int8_t ret[8];
        hash_t conv;
        int count = 0;
        for (int i = 0; i < kmer.length(); i++){
            char c = kmer[i];
            switch(c){
                //01
                case 'A':
                case 'a':
                    //                conv ^= (-x ^ conv) & (1 << n);
                    //                num << 1;
                    break;
                case 'T':
                case 't':
                    //                num << 2;
                    break;
                case 'C':
                case 'c':
                    //                num << 3;
                    break;
                case 'G':
                case 'g':
                    //                num << 4;
                    break;
                default:
                    //                num << 5;
                    break;

            }
        }


        return conv;
    }

    priority_queue<string> kmer_heap(string seq, vector<int> kmer){

        vector<string> base;

        //priority_queue<string> ret(base.begin(), base.end());
        for (auto k : kmer){
            vector<string> outmers(seq.length() - k, "");
#pragma omp parallel for
            for (int i = 0; i < seq.length() - k; i++){
                string forward = seq.substr(i, k);
                string revrev = reverse(reverse_complement(forward));

                //ret.push( (revrev < forward ? revrev : forward) );
                outmers[i] = (revrev < forward ? revrev : forward);
            }
            base.reserve(outmers.size() + base.size());
            base.insert(base.end(), outmers.begin(), outmers.end());
        }

        priority_queue<string> ret(base.begin(), base.end());
        return ret;
    }

    /* Returns the forward and reverse-reverse complement kmers of a sequence */
    vector<string> kmerize(string seq, int k){
        vector<string> ret(seq.length() - k, "");

#pragma omp parallel for
        for (int i = 0; i < seq.length() - k; i++){
            string s = seq.substr(i, k);
            //#pragma omp atomic read
            ret[i] = s;
            //ret.push_back(s);
            //ret.push_back(reverse(reverse_complement(s)));
        }
        return ret;
    }
    vector<string> multi_kmerize(string seq, vector<int> kSizes){
        int i = 0;
        vector<string> ret;
        //ret.reserve(kSizes.size() * 1000);
        for (auto k : kSizes){
            vector<string> kmers = kmerize(seq, k);
            ret.reserve(ret.size() + kmers.size());
            ret.insert(ret.end(), kmers.begin(), kmers.end());

            //for (i = 0; i + k < seq.length(); i++){
            //    ret.push_back(seq.substr(i, i+k));
            //    ret.push_back(reverse(reverse_complement(seq.substr(i, i+k))));
            //}


        }
        return ret;
    }

    /* Returns the forward shingles size k of a sequence */
    vector<string> shingle(string seq, int k){
        int i = 0;
        vector<string> ret;
        for (i = 0; i < seq.length() - k; i++){
            ret.push_back(seq.substr(i, k));
        }
        return ret;
    }

    vector<string> multi_shingle(string seq, vector<int> kSizes){
        int i = 0;
        vector<string> ret;
        for (auto k : kSizes){
            for (i = 0; i + k < seq.length(); i++){
                ret.push_back(seq.substr(i, k));
            }
        }
        return ret;
    }

    // vector<hash_t> preserve_kmer_mh64(string seq, vector<int> kSizes, int hashSize);
    vector<hash_t> allhash_unsorted_64(string& seq, vector<int>& k){
        int seqlen = seq.length();
        return allhash_unsorted_64_fast(seq.c_str(), k);

    }

    vector<hash_t> minhash_64_fast(string seq, vector<int> kmer, int hashSize, bool useBottom){
        vector<hash_t> ret;
        ret.reserve(seq.length() * kmer.size());
        for (auto k : kmer){
            vector<hash_t> tmp = calc_hashes(seq, k);
            ret.insert(ret.end(), tmp.begin(), tmp.end());
        }
        std::sort(ret.begin(), ret.end());

        int hashmax = hashSize < ret.size() ? hashSize : ret.size() - 1 ;

        return useBottom ?
            vector<hash_t> (ret.begin(), ret.begin() + hashmax) :
            vector<hash_t> (ret.rbegin(), ret.rbegin() + hashmax);
    }

    vector<hash_t> calc_hashes(string seq, int k){
        const char* x = seq.c_str();
        return calc_hashes(x, int(seq.length()), k);
    }

    hash_t calc_hash(string seq){
        int k = seq.length();
        char khash[16];
        char rev_rev_khash[16];
        const char* start = seq.c_str();
        char* rev_rev_s = new char[k]; 
        reverse_reverse_complement(start, rev_rev_s, k);
        if (!canonical(rev_rev_s, k)){
            //cerr << "Noncanonical bases found; exluding... " << rr_string << endl;
            //continue;     
        }
 
        MurmurHash3_x64_128(start, seq.size(), 42, khash);
        MurmurHash3_x64_128((const char*)rev_rev_s, seq.size(), 42, rev_rev_khash);

            //hash_t tmp_for = hash_t(khash[2]) << 32 | hash_t(khash[1]);
            //hash_t tmp_rev = hash_t(rev_rev_khash[2]) << 32 | hash_t(rev_rev_khash[1]);
        hash_t tmp_rev = *((hash_t *) rev_rev_khash);
        hash_t tmp_for = *((hash_t *) khash);
        return ( tmp_for < tmp_rev ? tmp_for : tmp_rev );

    }

    hash_t calc_hash(char* seq, int seqlen){
        char khash[16];
        char rev_rev_khash[16];
        const char* start = seq;
        char* rev_rev_s = new char[seqlen];
        reverse_reverse_complement(start, rev_rev_s, seqlen);
        if (!canonical(rev_rev_s, seqlen)){
            cerr << "Noncanonical bases found; exluding... " << rev_rev_s << endl;
            return 0;     
        }
            // need to handle reverse of char*
        MurmurHash3_x64_128(start, seqlen, 42, khash);
        MurmurHash3_x64_128( (const char*) rev_rev_s, seqlen, 42, rev_rev_khash);

            //hash_t tmp_for = hash_t(khash[2]) << 32 | hash_t(khash[1]);
            //hash_t tmp_rev = hash_t(rev_rev_khash[2]) << 32 | hash_t(rev_rev_khash[1]);
        hash_t tmp_rev = *((hash_t *) rev_rev_khash);
        hash_t tmp_for = *((hash_t *) khash);
        return ( tmp_for < tmp_rev ? tmp_for : tmp_rev );

       
    }

    vector<hash_t> calc_hashes(const char* x, int seq_length, int k){
        vector<hash_t> ret;
        ret.reserve(seq_length - k);
        //#pragma omp parallel for
        for (int i = 0; i < seq_length - k; i++){
            char khash[16];
            char rev_rev_khash[16];
            const char* start = x + i;

            char* rev_rev_s = new char[k];
            reverse_reverse_complement(start, rev_rev_s, k);
            if (!canonical(rev_rev_s, k)){
                //cerr << "Noncanonical bases found; exluding... " << rr_string << endl;
                continue;     
            }
            // need to handle reverse of char*
            MurmurHash3_x64_128(start, k, 42, khash);
            MurmurHash3_x64_128((const char*) rev_rev_s, k, 42, rev_rev_khash);

            //hash_t tmp_for = hash_t(khash[2]) << 32 | hash_t(khash[1]);
            //hash_t tmp_rev = hash_t(rev_rev_khash[2]) << 32 | hash_t(rev_rev_khash[1]);
            hash_t tmp_rev = *((hash_t *) rev_rev_khash);
            hash_t tmp_for = *((hash_t *) khash);
#pragma omp critical
            ret.push_back( tmp_for < tmp_rev ? tmp_for : tmp_rev );
        }
        return ret;
    }

    vector<hash_t> allhash_unsorted_64_fast(const char* seq, vector<int>& kmer){
        vector<hash_t> ret;
        ret.reserve(strlen(seq) * kmer.size());
        for (auto k : kmer){
            cerr << "K: " << k << endl;
            vector<hash_t> tmp = calc_hashes(seq, k);

            cerr << tmp.size() << endl;
            ret.insert(ret.end(), tmp.begin(), tmp.end());
        }

        cerr << "ret.size(): " << ret.size() << endl;;

        return ret;
    }


    tuple<hash_t*, int> allhash_unsorted_64_fast(const char* seq, int& seqlen, vector<int>& k_sizes){
        hash_t* ret;
        int ret_size = 0;
        for (auto k : k_sizes){
            ret_size += seqlen - k;
        }
        ret = new hash_t[ret_size];

        int track = 0;
        for (auto k : k_sizes){
            int i = 0;
            for (int i = 0; i < seqlen - k; i++){
                char khash[16];
                char rev_rev_khash[16];
                const char* start = seq + i;
                char* rev_rev_s = new char[k];
                reverse_reverse_complement(start, rev_rev_s, k);
                if (!canonical(start, k)){
                    ret[track + i] = 0;
                }
                else{
                    MurmurHash3_x64_128(start, k, 42, khash);
                    MurmurHash3_x64_128((const char *) rev_rev_s, k, 42, rev_rev_khash);

                    hash_t tmp_rev = *((hash_t *) rev_rev_khash);
                    hash_t tmp_for = *((hash_t *) khash);
                    delete [] rev_rev_s;
                    ret[ track + i ] =  tmp_for < tmp_rev ? tmp_for : tmp_rev;
                }
            }
            track += i;
        }

        return std::make_tuple(ret, ret_size);

    }


    vector<hash_t> minhashes(hash_t* hashes, int num_hashes, int sketch_size, bool useBottom){
       vector<hash_t> x = vector<hash_t>(hashes, hashes + num_hashes); 
       std::sort(x.begin(), x.end());

       int valid_ind = 0;
       while (x[valid_ind] == 0){
            valid_ind++;
       }

       /*for (auto xx : x){
           if (xx != 0){
            cerr << xx << endl;
           }
       }*/

       int hashmax = valid_ind + sketch_size < num_hashes ? valid_ind + sketch_size : num_hashes - 1;
       return std::vector<hash_t>(x.begin() + valid_ind, x.begin() + hashmax);
    }

    vector<hash_t> minhash_64(string& seq, vector<int>& k, int hashSize, bool useBottom){
        vector<hash_t> ret;
        ret.reserve(k.size() * seq.size());

        for (auto km_sz : k){
            vector<hash_t> tmp = calc_hashes(seq, km_sz);
            ret.insert(ret.end(), tmp.begin(), tmp.end());
        }
        std::sort(ret.begin(), ret.end());


        int hashmax = hashSize < ret.size() ? hashSize : ret.size() - 1 ;

        return useBottom ?
            vector<hash_t> (ret.begin(), ret.begin() + hashmax) :
            vector<hash_t> (ret.rbegin(), ret.rbegin() + hashmax);

    }

    vector<hash_t> minhash_64_depth_filter(vector<hash_t>& hashes, int hashSize, bool useBottom,
            int min_depth, unordered_map<hash_t, int>& hash_to_depth){
        vector<hash_t> ret;
        ret.reserve(hashes.size() / 2);
        for (int i = 0; i < hashes.size(); i++){
            if (hash_to_depth[hashes[i]] > min_depth){
#pragma omp critical
                ret.push_back(hashes[i]);
            }
        }

        /**
         * Special case if no hashes pass the depth filter.
         */
        if (ret.size() == 0){
            return ret;
        }

        std::sort(ret.begin(), ret.end());

        int hashmax = hashSize < ret.size() ? hashSize : ret.size() - 1 ;
        return useBottom ?
            vector<hash_t> (ret.begin(), ret.begin() + hashmax) :
            vector<hash_t> (ret.rbegin(),ret.rbegin() + hashmax);
    }

    vector<hash_t> minhash_64_depth_filter(string& seq, vector<int>& k,
            int hashSize, bool useBottom, int minDepth,
            unordered_map<hash_t, int>& hash_to_depth){

        vector<string> kmers = multi_kmerize(seq, k);

        vector<hash_t> ret;
        ret.reserve(kmers.size());

        //#pragma omp parallel for
        for (int i = 0; i < kmers.size(); i++){
            uint32_t khash[4];
            uint32_t rev_rev_khash[4];
            int str_length = kmers[i].size();
            const char* forward = kmers[i].c_str();

            string rrf = reverse(reverse_complement(kmers[i]));
            const char* rev_rev_forward = rrf.c_str();

            MurmurHash3_x64_128(forward, str_length, 42, khash);
            MurmurHash3_x64_128(rev_rev_forward, str_length, 42, rev_rev_khash);

            hash_t tmp_for = hash_t(khash[2]) << 32 | hash_t(khash[1]);
            hash_t tmp_rev = hash_t(rev_rev_khash[2]) << 32 | hash_t(rev_rev_khash[1]);

            hash_t r_hash = tmp_for < tmp_rev ? tmp_for : tmp_rev; //ret.push_back(r_hash);
            if (hash_to_depth[r_hash] > minDepth){
#pragma omp critical
                ret.push_back(r_hash);
            }
        }

        if (ret.size() == 0){
            return ret;
        }
        std::sort(ret.begin(), ret.end());

        int hashmax = hashSize < ret.size() ? hashSize : ret.size() - 1 ;

        return useBottom ?
            vector<hash_t> (ret.begin(), ret.begin() + hashmax) :
            vector<hash_t> (ret.rbegin(),ret.rbegin() + hashmax);
    }

    vector<hash_t> minhash_64(string seq, int k, int hashSize, bool useBottom){
        vector<string> kmers = kmerize(seq, k);
        vector<hash_t> ret(kmers.size(), 0);

        for (int i = 0; i < kmers.size(); i++){

            uint32_t seed = 42;
            uint32_t khash[4];
            uint32_t rev_rev_khash[4];
            int str_length = kmers[i].size();
            const char*
                forward = kmers[i].c_str();
            string rrf = reverse(reverse_complement(kmers[i]));
            const char*
                rev_rev_forward = rrf.c_str();

            MurmurHash3_x64_128(forward, str_length, seed, khash);
            MurmurHash3_x64_128(rev_rev_forward, str_length, seed, rev_rev_khash);

            hash_t tmp_rev = *((hash_t *) rev_rev_khash);
            hash_t tmp_for = *((hash_t *) khash);

            ret[i] = tmp_for < tmp_rev ? tmp_for : tmp_rev;
        }

        std::sort(ret.begin(), ret.end());


        if (useBottom){
            return vector<hash_t> (ret.begin(), ret.begin() + hashSize);
        }
        else{
            return vector<hash_t> (ret.rbegin(), ret.rbegin() + hashSize);
        }
    }

    vector<hash_t> top_minhash_64(string seq, int k, int hashSize){
        return minhash_64(seq, k, hashSize, false);
    }

    vector<hash_t> bottom_minhash_64(string seq, int k, int hashSize){
        return minhash_64(seq, k, hashSize, true);
    }

    vector<hash_t> hash_intersection(vector<hash_t> alpha, vector<hash_t> beta){
        vector<hash_t> ret;
        ret.reserve(alpha.size());
        int i = 0;
        int j = 0;
        while (i < alpha.size() && j < beta.size()){
            if (alpha[i] == beta[j]){
                ret.push_back(alpha[i]);
                i++;
                j++;
            }
            else if (alpha[i] > beta[j]){
                j++;
            }
            else{
                i++;
            }
        }

        return ret;
    }

    /**
     * If beta_len < alpha_len
     * j = 0
     * for (i in beta):
     *    if i > alpha[j]:
     *      j++
     *    elif i == alpha[j]:
     *      matches.add(alpha[j])
     *      i++
     *      j++
     */

    std::tuple<hash_t*, int> hash_intersection(hash_t alpha [ ], int alpha_start, int alpha_len,
                                                            hash_t beta[ ], int beta_start, int beta_len,
                                                            int sketch_size){
        hash_t* ret = new hash_t[ sketch_size ];
        int ret_len = 0;

        int i = alpha_start;
        int j = beta_start;
        while(i < alpha_len && j < beta_len){
            if (alpha[i] == beta[j]){
                ret[ret_len] = alpha[i];
                ++i;
                ++j;
                ++ret_len;
            }
            else if (alpha[i] > beta[j]){
                ++j;
            }
            else{
                ++i;
            }

        }
        
        //vector<hash_t> ret;
        //set_intersection(alpha + alpha_start, alpha + alpha_len,
        //                beta + beta_start, beta + beta_len,
        //                ret.begin());


        //return std::make_tuple(&(*ret.begin()), ret.size());
        return std::make_tuple(ret, ret_len);
    }

    vector<string> kmer_intersection(vector<string> alpha, vector<string> beta){
        vector<string> ret;
        ret.reserve(alpha.size());
        int i = 0;
        int j = 0;
        while(i < alpha.size() && j < beta.size()){
            if (alpha[i] == beta[j]){
                ret.push_back(alpha[i]);
                i++;
                j++;
            }
            else if (alpha[i] > beta[j]){
                j++;
            }
            else{
                i++;
            }
        }

        return ret;
    }

    vector<hash_t> hash_union(vector<hash_t> alpha, vector<hash_t> beta){
        vector<hash_t> ret;


        ret.reserve(alpha.size() + beta.size());
        ret = vector<hash_t> (alpha.begin(), alpha.end());
        ret.insert(ret.end(), beta.begin(), beta.end());
        return ret;
    }

    vector<hash_t> hash_set_intersection(vector<hash_t> alpha, vector<hash_t> beta){
        return hash_intersection(v_set(alpha), v_set(beta));

    }

    vector<hash_t> hash_set_union(vector<hash_t> alpha, vector<hash_t> beta){
        return v_set(hash_union(alpha, beta));
    }

    priority_queue<string> kmer_heap_intersection(priority_queue<string> alpha, priority_queue<string> beta){ 
        vector<string> base;
        base.reserve(alpha.size());
        while (alpha.size() != 0 && beta.size() != 0){
            if (alpha.top() == beta.top()){
                base.push_back(alpha.top());
                alpha.pop();
                beta.pop();
            }
            else if (alpha.top() > beta.top()){
                alpha.pop();
            }
            else if (alpha.top() < beta.top()){
                beta.pop();
            }
        }

        priority_queue<string> ret(base.begin(), base.end());
        return ret;
    }

}
