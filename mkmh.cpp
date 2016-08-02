#include "mkmh.hpp"
#include <string.h>

namespace mkmh{

    using namespace std;

    string reverse_complement(string seq){
        stringstream ret;

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
                    /* Handle X, N, Y, all that stuff. */
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

    const char* to_upper(char* seq, int length){

        std::locale loc;
        stringstream ret;
        for (int i = 0; i < length; i++){
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

        return ret.str().c_str();
    }

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

    int64_t kmer_to_integer(string kmer){
        //00 - 'G', 01 - 'A', 10 - 'T', and 11 - 'C'
        int8_t ret[8];
        int64_t conv;
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

    // vector<int64_t> preserve_kmer_mh64(string seq, vector<int> kSizes, int hashSize);
    vector<int64_t> allhash_unsorted_64(string& seq, vector<int>& k){
        vector<string> kmers = multi_kmerize(seq, k);
        //ret.reserve(kmers.size());

        vector<int64_t> ret(kmers.size(), 0);
        //#pragma omp parallel for
        for (int i = 0; i < kmers.size(); i++){
            uint32_t khash[4];
            uint32_t rev_rev_khash[4];
            //uint32_t seed = 101;
            int str_length = kmers[i].size();
            const char* forward = kmers[i].c_str();

            string rrf = reverse(reverse_complement(kmers[i]));
            const char* rev_rev_forward = rrf.c_str();

            MurmurHash3_x64_128(forward, str_length, 42, khash);
            MurmurHash3_x64_128(rev_rev_forward, str_length, 42, rev_rev_khash);

            int64_t tmp_for = int64_t(khash[2]) << 32 | int64_t(khash[1]);
            int64_t tmp_rev = int64_t(rev_rev_khash[2]) << 32 | int64_t(rev_rev_khash[1]);

            ret[i] = tmp_for < tmp_rev ? tmp_for : tmp_rev; //ret.push_back(r_hash);
        }

        return ret;
    }

    vector<int64_t> minhash_64_fast(string seq, vector<int> kmer, int hashSize, bool useBottom){
        vector<int64_t> ret;
        ret.reserve(seq.length() * kmer.size());
        for (auto k : kmer){
            vector<int64_t> tmp = calc_hashes(seq, k);
            ret.insert(ret.end(), tmp.begin(), tmp.end());
        }
        std::sort(ret.begin(), ret.end());
        
        int hashmax = hashSize < ret.size() ? hashSize : ret.size() - 1 ;

        return useBottom ?
            vector<int64_t> (ret.begin(), ret.begin() + hashmax) :
            vector<int64_t> (ret.rbegin(), ret.rbegin() + hashmax);
    }

    vector<int64_t> calc_hashes(string seq, int k){
        const char* x = seq.c_str();
        return calc_hashes(x, int(seq.length()), k);
    }

    vector<int64_t> calc_hashes(const char* x, int seq_length, int k){
        vector<int64_t> ret(seq_length - k, 0);        
        //#pragma omp parallel for
        for (int i = 0; i < seq_length - k; i++){
            uint32_t khash[4];
            uint32_t rev_rev_khash[4];
            const char* start = x + i;
            string rr_string = reverse(reverse_complement(string(start, k)));
            const char* rev_rev_s = rr_string.c_str();
            auto canonical = [](string x){
                bool allATGC = true;
                for (int i = 0; i < x.length(); i++){
                    char c = x[i];
                    switch (c){
                        case 'A':
                        case 'a':
                        case 'T':
                        case 't':
                        case 'C':
                        case 'c':
                        case 'G':
                        case 'g':
                            continue;
                            break;
                        default:
                            allATGC = false;
                            break;
                    }
                }
                return allATGC;
            };
            if (!canonical(rr_string)){
                //cerr << "Noncanonical bases found; exluding... " << rr_string << endl;
                continue;     
            }
            // need to handle reverse of char*
            MurmurHash3_x64_128(start, k, 42, khash);
            MurmurHash3_x64_128(rev_rev_s, k, 42, rev_rev_khash);
            int64_t tmp_for = int64_t(khash[2]) << 32 | int64_t(khash[1]);
            int64_t tmp_rev = int64_t(rev_rev_khash[2]) << 32 | int64_t(rev_rev_khash[1]);
            #pragma omp critical
            ret.push_back( tmp_for < tmp_rev ? tmp_for : tmp_rev );
        }
        return ret;
    }
    /*vector<int64_t> calc_hashes(const char* x, int seq_length, int k){
        vector<int64_t> ret(seq_length - k, 0);        
        //#pragma omp parallel for
        for (int i = 0; i < seq_length - k; i++){
            uint32_t khash[4];
            uint32_t rev_rev_khash[4];
            const char* start = x + i;
            string rr_string = reverse(reverse_complement(string(start, k)));
            const char* rev_rev_s = rr_string.c_str();
            std::function<bool(string)> canonical = [](string x){
                bool allATGC = true;
                for (int i = 0; i < x.length(); i++){
                    switch (i){
                        case 'A':
                        case 'a':
                        case 'T':
                        case 't':
                        case 'C':
                        case 'c':
                        case 'G':
                        case 'g':
                            continue;
                            break;
                        default:
                            allATGC = false;
                            break;
                    }
                }
                return allATGC;
            }
            if (!canonical(rr_string)){
                
            }
            // need to handle reverse of char*
            MurmurHash3_x64_128(start, k, 42, khash);
            MurmurHash3_x64_128(rev_rev_s, k, 42, rev_rev_khash);
            int64_t tmp_for = int64_t(khash[2]) << 32 | int64_t(khash[1]);
            int64_t tmp_rev = int64_t(rev_rev_khash[2]) << 32 | int64_t(rev_rev_khash[1]);
            ret[i] = tmp_for < tmp_rev ? tmp_for : tmp_rev;
        }
        return ret;
    }
*/
    vector<int64_t> allhash_unsorted_64_fast(const char* seq, vector<int>& k_sizes){
        //vector<string> kmers = multi_kmerize(seq, k);
        //ret size: ((len(seq) - k) for k in k_sizes)
        int ret_size = 0;
        vector<int> k_lens(k_sizes.size(), 0);
        for (int i = 0; i < k_sizes.size(); i++){
            ret_size += strlen(seq) - k_sizes[i];
            k_lens[i] = strlen(seq) - k_sizes[i];
        }
        int track_i = 0;
        vector<int64_t> ret;
        ret.reserve(ret_size);
        for (int k_track = 0; k_track < k_sizes.size(); k_track++){
            vector<int64_t> i_ret(k_lens[k_track], 0);
            //should replace loop end with klens[i]
            #pragma omp parallel for
            for (int i = 0; i < strlen(seq) - k_sizes[k_track]; i++){
                uint32_t khash[4];
                uint32_t rev_rev_khash[4];
                const char* start = seq + i;
                string rrs = reverse(reverse_complement(string(start, k_sizes[k_track]))).c_str();
                //cerr << string(start, k_sizes[k_track]) << endl;
                //cerr << rrs << endl;
                const char* rev_rev_s = rrs.c_str();                // need to handle reverse of char*
                MurmurHash3_x64_128(start, k_sizes[k_track], 42, khash);
                MurmurHash3_x64_128(rev_rev_s, k_sizes[k_track], 42, rev_rev_khash);
                int64_t tmp_for = int64_t(khash[2]) << 32 | int64_t(khash[1]);
                int64_t tmp_rev = int64_t(rev_rev_khash[2]) << 32 | int64_t(rev_rev_khash[1]);
                i_ret[i] = tmp_for < tmp_rev ? tmp_for : tmp_rev;
            }

            ret.insert(ret.end(), i_ret.begin(), i_ret.end());
        }

        return ret;
    }

    vector<int64_t> minhash_64(string& seq, vector<int>& k, int hashSize, bool useBottom){
        vector<int64_t> ret;
        ret.reserve(k.size() * seq.size());

        for (auto km_sz : k){
           vector<int64_t> tmp = calc_hashes(seq, km_sz);
           ret.insert(ret.end(), tmp.begin(), tmp.end());
        }
        std::sort(ret.begin(), ret.end());

        int hashmax = hashSize < ret.size() ? hashSize : ret.size() - 1 ;

        return useBottom ?
            vector<int64_t> (ret.begin(), ret.begin() + hashmax) :
            vector<int64_t> (ret.rbegin(), ret.rbegin() + hashmax);

    }

    vector<int64_t> minhash_64_depth_filter(vector<int64_t>& hashes, int hashSize, bool useBottom,
            int min_depth, unordered_map<int64_t, int>& hash_to_depth){
        vector<int64_t> ret;
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
            vector<int64_t> (ret.begin(), ret.begin() + hashmax) :
            vector<int64_t> (ret.rbegin(),ret.rbegin() + hashmax);
    }

    vector<int64_t> minhash_64_depth_filter(string& seq, vector<int>& k,
            int hashSize, bool useBottom, int minDepth,
            unordered_map<int64_t, int>& hash_to_depth){

        vector<string> kmers = multi_kmerize(seq, k);

        vector<int64_t> ret;
        ret.reserve(kmers.size());
        //#pragma omp parallel for
        for (int i = 0; i < kmers.size(); i++){
            uint32_t khash[4];
            uint32_t rev_rev_khash[4];
            //uint32_t seed = 101;
            int str_length = kmers[i].size();
            const char* forward = kmers[i].c_str();

            string rrf = reverse(reverse_complement(kmers[i]));
            const char* rev_rev_forward = rrf.c_str();

            MurmurHash3_x64_128(forward, str_length, 42, khash);
            MurmurHash3_x64_128(rev_rev_forward, str_length, 42, rev_rev_khash);

            int64_t tmp_for = int64_t(khash[2]) << 32 | int64_t(khash[1]);
            int64_t tmp_rev = int64_t(rev_rev_khash[2]) << 32 | int64_t(rev_rev_khash[1]);

            //int64_t r_hash = khash[2] < rev_rev_khash[2] ? int64_t(khash[2]) << 32 | int64_t(khash[1]) : int64_t(rev_rev_khash[2]) << 32 | int64_t(rev_rev_khash[1]);
            //ret[i] = r_hash;
            int64_t r_hash = tmp_for < tmp_rev ? tmp_for : tmp_rev; //ret.push_back(r_hash);
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
            vector<int64_t> (ret.begin(), ret.begin() + hashmax) :
            vector<int64_t> (ret.rbegin(),ret.rbegin() + hashmax);
    }

    vector<int64_t> minhash_64(string seq, int k, int hashSize, bool useBottom){
        vector<string> kmers = kmerize(seq, k);
        vector<int64_t> ret(kmers.size(), 0);

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

            int64_t tmp_for = int64_t(khash[2]) << 32 | int64_t(khash[1]);
            int64_t tmp_rev = int64_t(rev_rev_khash[2]) << 32 | int64_t(rev_rev_khash[1]);
            ;

            ret[i] = tmp_for < tmp_rev ? tmp_for : tmp_rev;
            //cerr << "kmer: " << kmers[i] << tmp_for << " " << tmp_rev << endl;
            //ret.push_back(tmp_for < tmp_rev ? tmp_for : tmp_rev);
        }

        std::sort(ret.begin(), ret.end());


        if (useBottom){
            return vector<int64_t> (ret.begin(), ret.begin() + hashSize);
        }
        else{
            return vector<int64_t> (ret.rbegin(), ret.rbegin() + hashSize);
        }
    }

    vector<int64_t> top_minhash_64(string seq, int k, int hashSize){
        return minhash_64(seq, k, hashSize, false);
    }

    vector<int64_t> bottom_minhash_64(string seq, int k, int hashSize){
        return minhash_64(seq, k, hashSize, true);
    }

    vector<int64_t> hash_intersection(vector<int64_t> alpha, vector<int64_t> beta){
        vector<int64_t> ret;
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

    vector<int64_t> hash_union(vector<int64_t> alpha, vector<int64_t> beta){
        vector<int64_t> ret;


        ret.reserve(alpha.size() + beta.size());
        ret = vector<int64_t> (alpha.begin(), alpha.end());
        ret.insert(ret.end(), beta.begin(), beta.end());
        return ret;
    }

    vector<int64_t> hash_set_intersection(vector<int64_t> alpha, vector<int64_t> beta){
        return hash_intersection(v_set(alpha), v_set(beta));

    }

    vector<int64_t> hash_set_union(vector<int64_t> alpha, vector<int64_t> beta){
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
