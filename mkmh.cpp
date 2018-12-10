#include "mkmh.hpp"
#include <string.h>

namespace mkmh{

    using namespace std;





    // void calc_hashes(const char* seq, const int& len,
    //         const int& k, hash_t*& hashes, int& numhashes, unordered_map<hash_t, int> counts){
    //     char* reverse = new char[k];
    //     uint32_t rhash[4];
    //     uint32_t fhash[4];
    //     //hash_t tmp_fwd;
    //     //hash_t tmp_rev;
    //     numhashes = len - k;
    //     hashes = new hash_t[numhashes];
    //     for (int i = 0; i < numhashes; ++i){
    //         if (canonical(seq + i, k)){
    //             reverse_complement(seq + i, reverse, k);
    //             MurmurHash3_x64_128(seq + i, k, 42, fhash);
    //             MurmurHash3_x64_128(reverse, k, 42, rhash);
    //             hash_t tmp_fwd = *((hash_t*) fhash);
    //             hash_t tmp_rev = *((hash_t*) rhash);
    //             hashes[i] = (tmp_fwd < tmp_rev ? tmp_fwd : tmp_rev);
    //             counts[hashes[i]] += 1;
    //         }
    //         else{
    //             hashes[i] = 0;
    //         }

    //     }
    //     delete reverse;
    // }


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

    tuple<vector<string>, vector<double>> sort_by_similarity(vector<hash_t> alpha, vector<vector<hash_t>> comps, vector<string> comp_names){
        vector<string> ret_names;
        vector<double> ret_sims;
        vector<pair<int, double>> helper_vec;

        for (int i = 0; i < comps.size(); ++i){
            int divisor = alpha.size();
            int shared = hash_intersection(alpha, comps[i]).size();
            double pct_id = (double) shared / (double) divisor;
            helper_vec.push_back(make_pair(i, pct_id));
        }
        sort( helper_vec.begin( ), helper_vec.end( ), [ ]( const pair<int, double>& lhs, const pair<int, double>& rhs )
        {
            return lhs.second > rhs.second;
        });
        
        for (int i = 0; i < helper_vec.size(); i++){
            ret_names.push_back(comp_names[helper_vec[i].first]);
            ret_sims.push_back(helper_vec[i].second);
        }

        return make_tuple(ret_names, ret_sims);
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
