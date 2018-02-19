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


namespace mkmh{
    using namespace std;

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
            cout << x[i] << " " << (int) x[i] << " " << valid_dna[ (int) x[i] ] << endl;
            trip |= valid_dna[x[i]];   
       }
        return !trip;
    };

    inline bool canonical(string x){
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

    /* Reverse complement the string seq
    * (assumes seq is DNA, and returns non-ACTG letters as-is*/

    /* Reverse complement a C string */
    inline void reverse_complement(const char* seq, char* ret, int len){
        assert(seq != ret);
        for (int i = len - 1; i >=0; i--){
            ret[ len - 1 - i ] = (char) rev_arr[ (int) seq[i] - 65];
        }
    };

    /* Reverse complement a string */
    inline string reverse_complement(string& seq){
        const char* s = seq.c_str();
        int seqlen = seq.length();
        char* ret = new char[seqlen];
        reverse_complement(s, ret, seqlen);
        string s_revc(ret);

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
    hash_t kmer_to_integer(string kmer);

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

    /** Calculate a MinHash sketch for kmers length (2 * k) with skip bases in between the two k-length halves **/
    vector<hash_t> allhash_64_linkmer(string seq, int k, int skip = 0);

    /* Returns the forward shingles of all k sizes of a sequence */
    /* Shingles are just forward-only kmers */
    vector<string> multi_shingle(string seq, vector<int> k);

    /* Return all hashes, unsorted*/
    vector<hash_t> allhash_unsorted_64(string& seq, vector<int>& k);

    vector<hash_t> allhash_unsorted_64_fast(const char* seq, vector<int>& k_sizes);
    
    
    tuple<hash_t*, int> allhash_unsorted_64_fast(const char* seq, int& seqlen, vector<int>& k_sizes);
   

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
    
    /* Calculate the hash for a string seq */
    hash_t calc_hash(string seq);

    /* Calculate the 64-bit hash for a string defined by seq and the length of seq */
    hash_t calc_hash(char* seq, int seqlen);

    inline void calc_hash(const char* seq, const int& len,
        char* reverse,
        hash_t* forhash, hash_t* revhash,
        hash_t* fin_hash){

        
        if (canonical(seq, len)){
            reverse_complement(seq, reverse, len);
            MurmurHash3_x64_128(seq, len, 42, forhash);
            MurmurHash3_x64_128(reverse, len, 42, revhash);
            *fin_hash = *forhash < *revhash ? *forhash : *revhash;
        }
        else{
            *forhash = 0;
            *revhash = 0;
            *fin_hash = 0;
        }
    };

    /* Calculate all the hashes of the kmers length k of seq */
    vector<hash_t> calc_hashes(string seq, int k);

    /* Calculate all the hashes of the kmers length k of seq */
    vector<hash_t> calc_hashes(const char* seq, int seq_length, int k);

    /** Calculate the hashes of the kmers length k of seq **/
    void calc_hashes(const char* seq, const int& len,
            const int& k, hash_t* hashes, int& numhashes);

    /** Calculate the hashes for kmers of multiple lengths in <kmer>
     */
    inline void calc_hashes(const char* seq, const int& len,
                vector<int>& kmer, hash_t* hashes, int& num_hashes){
        
        assert(num_hashes == 0);
        hash_t* tmp;
        for (int i = 0; i < kmer.size(); ++i){
            int k = kmer[i];
            hash_t* khashes = new hash_t[len - k];
            int n = len - k;
            calc_hashes(seq, len, k, khashes, n);
            tmp = new hash_t[num_hashes + n];
            copy(hashes, hashes + num_hashes, tmp);
            copy(khashes, khashes + n, tmp + num_hashes);
            delete [] hashes;
            delete [] khashes;
            hashes = tmp;
            num_hashes += n;
        }
    };

    inline void calc_hashes(const char* seq, const int& len,
            const int& k, hash_t* hashes, int& numhashes){
        char* reverse = new char[k];
        hash_t* rhash = new hash_t[k];
        numhashes = len - k;
        hashes = new hash_t[numhashes];
        if (canonical(seq, len)){
            for (int i = 0; i < numhashes; ++i){
                reverse_complement(seq + i, reverse, k);
                MurmurHash3_x64_128(seq + i, k, 42, hashes + i);
                MurmurHash3_x64_128(reverse, k, 42, rhash);
                *(hashes + i) = *(hashes + i) < *(rhash) ? *(hashes + i) : *(rhash);
            }
        }
        delete reverse;
        delete rhash;
    };

   
}

#endif
