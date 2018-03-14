#include "mkmh.hpp"
#include <string.h>

namespace mkmh{

    using namespace std;

    vector<string> kmer_set(vector<string> kmers){
        set<string> uniqs = set<string> (kmers.begin(), kmers.end());
        vector<string> ret = vector<string> (uniqs.begin(), uniqs.end());
        return ret;
    }

    priority_queue<string> kmer_heap(string seq, vector<int> kmer){

        vector<string> base;

        //priority_queue<string> ret(base.begin(), base.end());
        for (auto k : kmer){
            vector<string> outmers(seq.length() - k, "");
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

    /* Returns the forward kmers of a sequence */
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

    mkmh_kmer_list_t kmerize(char* seq, int seq_len, int k){
        mkmh_kmer_list_t ret;
        ret.kmers = new char*[seq_len - k];
        ret.k = k;
        ret.length = seq_len - k;

        for (int i = 0; i < ret.length; ++i){
            char* km = new char[k + 1];
            memcpy(km, seq + i, k);
            ret.kmers[i] = new char[k + 1];
            ret.kmers[i] = km;
            ret.kmers[i][k] = '\0';
        }
        return ret;
    }

    void kmerize(char* seq, const int& seq_len, const int& k, char** kmers, int& kmer_num){
        char** ret = new char*[seq_len - k];
        kmer_num = seq_len - k;
        for (int i = 0; i < kmer_num; ++i){
            ret[i] = new char[ k + 1 ];
            memcpy(ret[i], seq + i, k);
        }
    }

    void calc_hashes(const char* seq, const int& len,
            const int& k, hash_t*& hashes, int& numhashes, unordered_map<hash_t, int> counts){
        char* reverse = new char[k];
        uint32_t rhash[4];
        uint32_t fhash[4];
        //hash_t tmp_fwd;
        //hash_t tmp_rev;
        numhashes = len - k;
        hashes = new hash_t[numhashes];
        for (int i = 0; i < numhashes; ++i){
            if (canonical(seq + i, k)){
                reverse_complement(seq + i, reverse, k);
                MurmurHash3_x64_128(seq + i, k, 42, fhash);
                MurmurHash3_x64_128(reverse, k, 42, rhash);
                hash_t tmp_fwd = *((hash_t*) fhash);
                hash_t tmp_rev = *((hash_t*) rhash);
                hashes[i] = (tmp_fwd < tmp_rev ? tmp_fwd : tmp_rev);
                counts[hashes[i]] += 1;
            }
            else{
                hashes[i] = 0;
            }

        }
        delete reverse;
    }

    void print_kmers(char* seq, const int& len, int k){
        int kmerized_length = len - k;
        stringstream st;
        for (int i = 0; i < kmerized_length - 1; ++i){
            int j = 0;
            while (j < k){
                st << seq[i + j];
                ++j;
            }
            st << "\t";
        }
        int j = 0;
        while(j < k){
            st << seq[kmerized_length - 1 + j];
            ++j;
        }
        st << endl;
        cout << st.str();
        st.str("");
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

    /* Returns only the forward shingles size k of a sequence */
    vector<string> shingle(string seq, int k){
        int i = 0;
        vector<string> ret;
        for (i = 0; i < seq.length() - k; i++){
            ret.push_back(seq.substr(i, k));
        }
        return ret;
    }

    vector<mkmh_minimizer> kmer_tuples(string seq, int k){
        vector<string> kmers = kmerize(seq, k);
        vector<mkmh_minimizer> tups (kmers.size());
        for (int i = 0; i < kmers.size(); i++){
            mkmh_minimizer mm;
            mm.seq = kmers[i];
            mm.pos = i;
            mm.length = k;
            tups[i] = mm;
        }

        return tups;
    }

    vector<hash_t> allhash_64_linkmer(string seq, int k, int skip){
        vector<hash_t> ret(seq.size() - (k*2));
        int last_kmer_ind = seq.size() - (skip + (2 * k) );

        #pragma omp parallel for
        for (int i = 0; i <= last_kmer_ind; ++i){
            char * linkmer = new char [k * 2];
            for (int j = 0; j < k; ++j){
                linkmer[j] = seq[i + j];
                linkmer[k + j] = seq[i + skip + k + j];
            }
            hash_t c_hash = calc_hash(linkmer, 2*k);
            ret[i] = c_hash;
            delete [] linkmer;
        }
        
        return ret;
    }

    vector<mkmh_minimizer> minimizers(string seq, int k, int w){
        vector<mkmh_minimizer> ret;
        vector<mkmh_minimizer> kmert = kmer_tuples(seq, k);
        int i = 0;
        for (i = 0; i + w < kmert.size(); ++i){
            // get and sort kmers in window (i, i + w)
            vector<mkmh_minimizer> window_kmers(kmert.begin() + i, kmert.begin() + i + w);
            std::sort(window_kmers.begin(), window_kmers.end());
            // TODO filter minimizers if needed, e.g. to remove poly-As
            // create an mkmh_minimizer struct
            // tuck minimizer in ret
            ret.push_back(*(window_kmers.begin()));
        }
        return v_set(ret);
    }



    vector<mkmh_minimizer> unreduced_minimizers(string seq, int k, int w){
        vector<mkmh_minimizer> ret;
        vector<mkmh_minimizer> kmert = kmer_tuples(seq, k);
        int i = 0;
        for (i = 0; i + w < kmert.size(); ++i){
            // get and sort kmers in window (i, i + w)
            vector<mkmh_minimizer> window_kmers(kmert.begin() + i, kmert.begin() + i + w);
            std::sort(window_kmers.begin(), window_kmers.end());
            // TODO filter minimizers if needed, e.g. to remove poly-As
            // create an mkmh_minimizer struct
            // tuck miimizer in ret
            ret.push_back(*(window_kmers.begin()));
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
        //ret.reserve(k.size() * seq.size());

        for (auto km_sz : k){
            vector<hash_t> tmp = calc_hashes(seq, km_sz);
            ret.insert(ret.end(), tmp.begin(), tmp.end());
        }
        std::sort(ret.begin(), ret.end());
        int nonzero_ind = 0;
        while (ret[nonzero_ind] == 0){
            nonzero_ind++;
        }

        hashSize += nonzero_ind;
        int hashmax = hashSize < ret.size() ? hashSize : ret.size() - 1 ;



        return useBottom ?
            vector<hash_t> (ret.begin() + nonzero_ind, ret.begin() + hashmax) :
            vector<hash_t> (ret.rbegin() + nonzero_ind, ret.rbegin() + hashmax);

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
        while (alpha[i] == 0){
            i++;
        }
        while(beta[j] == 0){
            j++;
        }
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
