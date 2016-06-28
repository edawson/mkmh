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

    /* Returns the forward and reverse-reverse complement kmers of a sequence */
    vector<string> kmerize(string seq, int k){
        vector<string> ret(seq.length() - k, "");
        
        #pragma omp parallel for if (seq.length() > 1000)
        for (int i = 0; i < seq.length() - k; i++){
            string s = seq.substr(i, k);
            //#pragma omp atomic read
            ret[i] = s;
            //ret.push_back(s);
            //ret.push_back(reverse(reverse_complement(s)));
        }
        return ret;
    }
/*
    vector<string> kmerize(string seq, int k){
        int i = 0;
        int len_kmers = 2 * seq.length() - k;
        vector<string> ret(0, len_kmers);
        #pragma omp parallel for
        for (i = 0; i + k < seq.length(); i++){
            string s = seq.substr(i, k);
            ret[i] = s;
            ret[i + len_kmers] = reverse(reverse_complement(s));
        }
        return ret;
    }
*/

    vector<string> multi_kmerize(string seq, vector<int> kSizes){
        int i = 0;
        vector<string> ret;
        ret.reserve(kSizes.size() * 1000);
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
    

    vector<int64_t> minhash_64(string seq, vector<int> k, int hashSize, bool useBottom){
        vector<string> kmers = multi_kmerize(seq, k);
        //ret.reserve(kmers.size());

        vector<int64_t> ret(kmers.size(), 0);
        //const char* forward;
        //const char* rev_rev_forward;
        #pragma omp parallel for
        for (int i = 0; i < kmers.size(); i++){
            uint32_t khash[4];
            uint32_t rev_rev_khash[4];
            //uint32_t seed = 101;
            const char* forward = kmers[i].c_str();

            string rrf = reverse(reverse_complement(kmers[i]));
            const char* rev_rev_forward = rrf.c_str();

            MurmurHash3_x64_128(forward, strlen(forward), 101, khash);
            MurmurHash3_x64_128(rev_rev_forward, strlen(rev_rev_forward), 101, rev_rev_khash);
            
            int64_t tmp_for = int64_t(khash[2]) << 32 | int64_t(khash[1]);
            int64_t tmp_rev = int64_t(rev_rev_khash[2]) << 32 | int64_t(rev_rev_khash[1]);
            
            //int64_t r_hash = khash[2] < rev_rev_khash[2] ? int64_t(khash[2]) << 32 | int64_t(khash[1]) : int64_t(rev_rev_khash[2]) << 32 | int64_t(rev_rev_khash[1]);
            //ret[i] = r_hash;
            ret[i] = tmp_for < tmp_rev ? tmp_for : tmp_rev; //ret.push_back(r_hash);
        }

        std::sort(ret.begin(), ret.end());

        return useBottom ?
            vector<int64_t> (ret.begin(), ret.begin() + hashSize) :
            vector<int64_t> (ret.rbegin(),ret.rbegin() + hashSize);

    }

    vector<int64_t> minhash_64(string seq, int k, int hashSize, bool useBottom){
        vector<string> kmers = kmerize(seq, k);
        vector<int64_t> ret(kmers.size(), 0);

        //ret.reserve(kmers.size());
        //vector<string>::iterator it;
        //const char* forward;
        //const char* rev_rev_forward;
        //for (it = kmers.begin(); it != kmers.end(); it++){
        for (int i = 0; i < kmers.size(); i++){

            uint32_t seed = 101;
            uint32_t khash[4];
            uint32_t rev_rev_khash[4];
            const char*
            forward = kmers[i].c_str();
            string rrf = reverse(reverse_complement(kmers[i]));
            const char*
            rev_rev_forward = rrf.c_str();
            
            MurmurHash3_x64_128(forward, strlen(forward), seed, khash);
            MurmurHash3_x64_128(rev_rev_forward, strlen(rev_rev_forward), seed, rev_rev_khash);
            
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
        //while (i < alpha.size() && j < beta.size()){
        //#pragma omp parallel for private(i, j)
        for (i = 0, j = 0; i < alpha.size(), j < beta.size();){
            if (alpha[i] == beta[j]){
                //#pragma omp critical
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

}
