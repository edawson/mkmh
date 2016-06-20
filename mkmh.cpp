#include "mkmh.hpp"

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
                    ret << "c";
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

    /* Returns the forward and reverse-reverse complement kmers of a sequence */
    vector<string> kmerize(string seq, int k){
        int i = 0;
        vector<string> ret;
        for (i = 0; i + k < seq.length(); i++){
            ret.push_back(seq.substr(i, i+k));
            ret.push_back(reverse(reverse_complement(seq.substr(i, i+k))));
        }
        return ret;
    }

    vector<string> multi_kmerize(string seq, vector<int> kSizes){
        int i = 0;
        vector<string> ret;
        for (auto k : kSizes){
            ret.extend(kmerize(seq, k));
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
            ret.push_back(seq.substr(i, i+k));
        }
        return ret;
    }

    vector<string> multi_shingle(string seq, vector<int> kSizes){
        int i = 0;
        vector<string> ret;
        for (auto k : kSizes){
            for (i = 0; i + k < seq.length(); i++){
                ret.push_back(seq.substr(i, i+k));
            }
        }
        return ret;
    }

    // vector<int64_t> preserve_kmer_mh64(string seq, vector<int> kSizes, int hashSize);

    vector<int64_t> minhash_64(string seq, int k, int hashSize, bool useBottom){
        vector<int64_t> ret;

        return ret;
    }

    vector<int64_t> top_minhash_64(string seq, int k, int hashSize){
        return minhash_64(seq, k, hashSize, false);
    }

    vector<int64_t> bottom_minhash_64(string seq, int k, int hashSize){
        return minhash_64(seq, k, hashSize, true);
    }

    vector<int64_t> hash_union(vector<int64_t> alpha, vector<int64_t> beta){
        vector<int64_t> ret;

        return ret;
    }

    vector<int64_t> hash_intersection(vector<int64_t> alpha, vector<int64_t> beta){
        vector<int64_t> ret;

        return ret;
    }
}
