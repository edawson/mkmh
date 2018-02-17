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
#include <omp.h>
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

    inline bool canonical(const char* x, int len){
        bool allATGC = true;
        for (int i = 0; i < len; ++i){
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

    /* Reverse the string seq */
    inline string reverse(string seq){
        string copy = string(seq);
        std::reverse(copy.begin(), copy.end());
        return copy;
    }

    /* Reverse complement the string seq (assumes seq is DNA, and returns non-ACTG letters as-is*/
    string reverse_complement(string& seq);
    
    void reverse_complement(const char* seq, char* ret, int len);

    /* Capitalize all characters in a string */
    string to_upper(string seq);
    
    void to_upper(char* seq, int length);

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


    //TODO void calc_hash(char* seq, int& len, hash_t* h);
    //

    /* Calculate all the hashes of the kmers length k of seq */
    vector<hash_t> calc_hashes(string seq, int k);

    /* Calculate all the hashes of the kmers length k of seq */
    vector<hash_t> calc_hashes(const char* seq, int seq_length, int k);


}

#endif
