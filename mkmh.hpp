#ifndef MKMH_D
#define MKMH_D

#include <vector>
#include <queue>
#include <set>
#include <unordered_map>
#include <string>
#include <sstream>
#include <locale>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <omp.h>
#include <assert.h>

#include "murmur3.hpp"
#include "HASHTCounter.hpp"

#define DBGG

namespace mkmh{
    using namespace std;
    using namespace mkmh;

    typedef uint64_t hash_t;

    struct mkmh_kmer_list_t{
        char** kmers;
        int length;
        int k;

        ~mkmh_kmer_list_t(){
            for (int i = 0; i < length; ++i){
                delete [] kmers[i];
            }
            delete [] kmers;
        }
    };

    struct mkmh_minimizer {
        uint64_t pos;
        uint32_t length;
        string seq;
        bool operator<(const mkmh_minimizer& rhs) const {return seq < rhs.seq;};
    };

    // Crazy hack char table to test for canonical bases
    static const int valid_dna[127] = {
        1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,0,1,0,1,1,1,
        0,1,1,1,1,1,1,1,1,1,
        1,1,1,0,1,1,1,1,1,1,
        1,1,1,1,1,1,0,1,0,1,
        1,1,0,1,1,1,1,1,1,1,
        1,1,1,1,1,0,1,1,1,1,
        1,1,1,1,1,1 
    };

    // Reverse complement lookup table
    static char rev_arr [26] = {
        84, 66, 71, 68, 69,
        70, 67, 72, 73, 74,
        75, 76, 77, 78, 79,
        80, 81, 82, 83, 65,
        85, 86, 87, 88, 89, 90
    };


    // Check a string (as a char*) for non-canonical DNA bases
    inline bool canonical(const char* x, int len){
        bool trip = false;
        for (int i = 0; i < len; ++i){
            trip |= valid_dna[x[i]];   
        }
        return !trip;
    };

    inline bool canonical(string seq){
        const char* x = seq.c_str();
        int len = seq.length();
        return canonical(x, len);
    };

    /* Reverse complement the string seq
     * (assumes seq is DNA, and returns non-ACTG letters as-is*/

    /* Reverse complement a C string
     * NB: does not check safety of string lengths.
     * NB: ret is modified to hold the reverse complement of seq.
     * */

    inline void reverse_complement(const char* seq, char* ret, int len){
        
        //assert(seq != ret);
        
        for (int i = len - 1; i >=0; i--){
            ret[ len - 1 - i ] = (char) rev_arr[ (int) seq[i] - 65];
        }
        ret[len] = '\0';
    };

    /* Reverse complement a string */
    inline string reverse_complement(string& seq){
        const char* s = seq.c_str();
        int seqlen = seq.length();
        char* ret = new char[seqlen];
        reverse_complement(s, ret, seqlen);
        string s_revc(ret);
        delete ret;

        return s_revc;
    };

    /* Capitalize all characters in a string */
    /* Capitalize a C string */
    inline void to_upper(char* seq, int length){
        for (int i = 0; i < length; i++){
            char c = seq[i];
            seq[i] = ( (c - 91) > 0 ? c - 32 : c);
        }
    };

    /* Capitalize a string */
    inline string to_upper(string& seq){
        for (int i = 0; i < seq.length(); i++){
            char c = seq[i];
            seq[i] =  ((c - 91) > 0 ? c - 32 : c);
        }
        return seq;
    };

    /* Reverse a string */
    inline string reverse(string seq){
        string copy = string(seq);
        std::reverse(copy.begin(), copy.end());
        return copy;
    };

    /* Reverse a C string*/
    inline void reverse(char* seq, const int& len){
        char tmp;
        for (int i = 0; i < len; ++i){
            tmp = seq[len - i];
            seq[ len - i ] = seq[i];
            seq[i] = tmp;
        }
    };

    inline void sort(hash_t* hashes, int len, bool descending = false){
        if (!descending){
           std::sort(hashes, hashes + len);
        }else{
            std::sort(hashes, hashes + len, std::less<uint64_t>());

        }
    };

    inline void sort(vector<hash_t>& hashes){
        std::sort(hashes.begin(), hashes.end());
    };



    /* Returns the forward and reverse-reverse complement kmers of a sequence */
    vector<string> kmerize(string seq, int k);

    mkmh_kmer_list_t kmerize(char* seq, int seq_len, int k);

    void kmerize(char* seq, const int& seq_len, const int& k, char** kmers, int& kmer_num);

    /* Print the kmers of a string, tab separated, to cout 
     *   avoids allocating any new memory. */
    void print_kmers(char* seq, const int& seq_len, int k);

    /* Returns the forward and reverse-reverse complement kmers for all kmer sizes in k */
    vector<string> multi_kmerize(string seq, vector<int> k);

    /* Returns a deduplicated set of string kmers */
    vector<string> kmer_set(vector<string> kmers);

    /* Returns a heap (priority queue) of the kmers of a read converted to ints. */

    /* Returns a heap (priority queue) of the kmers of the read */
    priority_queue<string> kmer_heap(string seq, vector<int> k);

    /* Converts a string kmer to an integer representation */
    //hash_t kmer_to_integer(string kmer);

    /* Returns a deduplicated set of kmers or hashes as a vector<T> */
    template<typename T>
        inline vector<T> v_set(vector<T> kmers){
            set<T> s = set<T>(kmers.begin(), kmers.end());
            vector<T> ret = vector<T>(s.begin(), s.end());
            return ret;
        }

    /* Returns the forward shingles size k of a sequence */
    vector<string> shingle(string seq, int k);

    /** Returns an mkmh_minimizer struct, equivalent to a tuple(kmer, position, kmer length), for every position in the genome **/
    vector<mkmh_minimizer> kmer_tuples(string seq, int k);

    /** Finds the (w, k) minimizers of a string **/
    vector<mkmh_minimizer> minimizers(string seq, int k, int w); 

    /** Finds the (w,k) minimizers and reports all of them (including duplicates) **/
    vector<mkmh_minimizer> unreduced_minimizers(string seq, int k, int w); 


    // Calculate a single hash of a sequence (usually a kmer).
    // Takes in:
    //      seq: a char* (not affected)
    //      len: the length of seq (not affected)
    //      reverse: a char* of length(seq); modified to hold reverse_complement of (seq)
    //      forhash: a 128-bit int (e.g. uint32_t[4]) for holding forward hash.
    //      revhash: a 128-bit int (e.g. uint32_t[4]) for holding reverse hash.
    //      finhash: a 64-bit int (e.g. uint64_t) which holds the lesser of (forhash, revhash);
    inline void calc_hash(const char* seq, const int& len,
            char*& reverse,
            hash_t*& forhash, hash_t*& revhash,
            hash_t*& fin_hash){


        if (canonical(seq, len)){
            reverse_complement(seq, reverse, len);
            MurmurHash3_x64_128(seq, len, 42, forhash);
            MurmurHash3_x64_128(reverse, len, 42, revhash);
            *fin_hash = (*((hash_t*)forhash) < *((hash_t*) revhash) ? *((hash_t*)forhash) : *((hash_t*)revhash));
        }
        else{
            *forhash = 0;
            *revhash = 0;
            *fin_hash = 0;
        }
    };



    /* Calculate the 64-bit hash for a string defined by seq and the length of seq */
    inline hash_t calc_hash(const char* seq, int seqlen){
        char* reverse = new char[seqlen];
        hash_t* fhash = new hash_t [4];
        hash_t* rhash = new hash_t [4];
        hash_t* fin_hash = new hash_t [1];
        calc_hash(seq, seqlen, reverse, fhash, rhash, fin_hash);
        hash_t ret = *(fin_hash);
        delete [] fhash;
        delete [] rhash;
        delete [] fin_hash;
        return ret;
    }

    /* Calculate the hash for a string seq */
    inline hash_t calc_hash(string seq){
        int k = seq.length();
        const char* x = seq.c_str();
        return calc_hash(x, k);
    };



    /** Primary calc_hashes function **/
    /** Takes the string to be hashed, its length,
     *  a single kmer size, a pointer to hold the hashes,
     *  and an integer to hold the number of hashes.
     *  
     *  Possibly thread safe:
     *      seq, len and k are not modified
     *      new [] operator is known threadsafe
     *      User must handle hashes and numhashes properly in calling function.
     **/
    inline void calc_hashes(const char* seq, const int& len,
            const int& k, hash_t*& hashes, int& numhashes){
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
            }
            else{
                hashes[i] = 0;
            }

        }
        delete reverse;
    };

    /** calc_hashes for multiple kmers sizes **/
    /** returns an array, with hashes for each kmer size concatenated
     * to those of the previous kmer size
     **/
    inline void calc_hashes(const char* seq, int seq_length,
     vector<int> kmer_sizes,
     hash_t*& hashes, int& numhashes){
        numhashes = 0;

        // This holds the number of hashes preceeding the
        // kmer size currently being hashed.
        vector<int> offsets;
        for (auto k : kmer_sizes){
            offsets.push_back(numhashes);
            numhashes += seq_length - k;
        }
        hashes = new hash_t [numhashes];

        for (int i = 0; i < kmer_sizes.size(); ++i){
            int k = kmer_sizes[i];
            int local_numhash;
            hash_t* l_start = hashes + offsets[i];
            calc_hashes(seq, seq_length, k, l_start, local_numhash);
        }
    }


    /* Calculate all the hashes of the kmers length k of seq */
    inline vector<hash_t> calc_hashes(const char* seq, int seq_length, int k){
        int numhashes = 0;
        hash_t* hashes;
        calc_hashes(seq, seq_length, k, hashes, numhashes);
        vector<hash_t> ret(numhashes);
        for (int i = 0; i < numhashes; i++){
            ret[i] = *(hashes + i);
        }

        delete [] hashes;

        return ret;
    };

    /* Calculate all the hashes of the kmers length k of seq */
    inline vector<hash_t> calc_hashes(string seq, int k){
        const char* x = seq.c_str();
        int l = seq.length();
        return calc_hashes(x, l, k);
    }
    
    inline vector<hash_t> calc_hashes(const char* seq, const int& len, const vector<int>& k_sizes){
        vector<hash_t> ret;
        for (auto k : k_sizes){
            vector<hash_t> t = calc_hashes(seq, len, k);
            ret.insert(ret.end(), t.begin(), t.end());
        }
        return ret;
    };

    /** Calculate the hashes of seq
     *  and fill in a HASHTCounter htc so that
     *  hashes can be kept or removed based on the number of times
     *  they occur in seq.
     **/
    void calc_hashes(const char* seq, const int& len,
            const int& k, hash_t*& hashes, int& numhashes, HASHTCounter htc);
    
    /** Calculate the hashes for kmers of multiple lengths in <kmer>
    */
    inline vector<hash_t> calc_hashes(string seq, const vector<int>& k_sizes){
        const char* x = seq.c_str();
        int l = seq.length();
        return calc_hashes(x, l, k_sizes);
    };

    inline void hash_intersection_size(const hash_t* alpha, const int& alpha_size, const hash_t* beta, const int& beta_size, int& ret){
        int a_ind = 0;
        int b_ind = 0;
        ret = 0;
        while (a_ind < alpha_size && b_ind < beta_size){
            if (alpha[a_ind] == beta[b_ind]){
                ++ret;
                ++a_ind;
                ++b_ind;
            }
            else if (alpha[a_ind] < beta[b_ind]){
                ++a_ind;
            }
            else{
                ++b_ind;
            }
        }
        
    };

    /** Calculate a MinHash sketch for kmers length (2 * k) with skip bases in between the two k-length halves **/
    vector<hash_t> allhash_64_linkmer(string seq, int k, int skip = 0);

    /* Returns the forward shingles of all k sizes of a sequence */
    /* Shingles are just forward-only kmers */
    vector<string> multi_shingle(string seq, vector<int> k);

    /** MinHash - given an array of hashes, modify the mins array to hold 
     * the lowest/highest N (excluding zeros) **/
    inline void minhashes(hash_t* hashes, int num_hashes,
             int sketch_size,
             hash_t*& ret,
             int& retsize,
            bool use_bottom=true){
        
        ret = new hash_t[sketch_size];
        mkmh::sort(hashes, num_hashes, !use_bottom);

        int maxlen = min(num_hashes, sketch_size);
        int start = 0;
        while (retsize < sketch_size && start < num_hashes){
            if (hashes[start] != 0){
                ret[retsize] = hashes[start];
                ++retsize;
            }
            ++start;
        }
    };

    /** Base MinHash function - return the lowest n = min(num_hashes, sketch_size) hashes. **/
    vector<hash_t> minhashes(hash_t* hashes, int num_hashes, int sketch_size, bool useBottom=true);    

    /* Returns the lowest hashSize hashes of the kmers (length k...k` in k) of seq */
    vector<hash_t> minhash_64(string& seq, vector<int>& k, int hashSize, bool useBottom=true);

    /* Returns the bottom/top hashSize hashes of kmers size k in seq */ 
    vector<hash_t> minhash_64(string seq, int k, int hashSize, bool useBottom=true);

    /* helper function: returns the top hashSize hashes of the kmers size k in seq */
    vector<hash_t> top_minhash_64(string seq, int k, int hashSize);

    /* helper function: returns the bottom hashSize hashes of the kmers size k in seq */
    vector<hash_t> bottom_minhash_64(string seq, int k, int hashSize);

    /* Returns the bottom/top hashSize hashes of kmers size k which 
     * occur more than minDepth times, based on the depth in hash_to_depth */
    vector<hash_t> minhash_64_depth_filter(string& seq, vector<int>& k,
            int hashSize, bool useBottom, int minDepth,
            unordered_map<hash_t, int>& hash_to_depth);

    /* Takes in a list of pre-computed hashes and returns the MinHash (size hashSize)
     * of the hashes that pass the depth filter */
    vector<hash_t> minhash_64_depth_filter(vector<hash_t>& hashes, int hashSize, bool useBottom,
            int min_depth, unordered_map<hash_t, int>& hash_to_depth);

    vector<hash_t> minhash_64_fast(string seq, vector<int> kmer, int sketchSize, bool isBottom=true);

    /* Returns the union of the hashes in alpha and beta, including duplicates */
    vector<hash_t> hash_union(vector<hash_t> alpha, vector<hash_t> beta);

    /* Returns the intersection of alpha and beta, including duplicates */   
    std::tuple<hash_t*, int> hash_intersection(hash_t* alpha, int alpha_start, int alpha_len,
            hash_t* beta, int beta_start, int beta_len,
            int sketch_size);

    /* Returns the intersection of alpha and beta, including duplicates the number of
       times they appear in both vectors */
    vector<hash_t> hash_intersection(vector<hash_t> alpha, vector<hash_t> beta);

    /* Returns the union of the two sets after deduplicating all duplicates */
    vector<hash_t> hash_set_union(vector<hash_t> alpha, vector<hash_t> beta);

    /* Returns the intersection of both sets. Duplicates are included only once */
    vector<hash_t> hash_set_intersection(vector<hash_t> alpha, vector<hash_t> beta);

    /* Returns two vectors, one of sequence names and one of percent similarity, sorted by percent similarity to alpha.
     * NB: input vectors should be sorted. */
    tuple<vector<string>, vector<double>> sort_by_similarity(vector<hash_t> alpha, vector<vector<hash_t>> comps, vector<string> comp_names);

    /** Return the intersection of two kmer heaps **/
    priority_queue<string> kmer_heap_intersection(priority_queue<string> alpha, priority_queue<string> beta);

    /** Return the intersection of two lists of kmers **/
    vector<string> kmer_intersection(vector<string> alpha, vector<string> beta);
}

#endif
