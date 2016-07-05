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

    /* Reverse the string seq */
    string reverse(string seq);
    
    /* Reverse complement the string seq (assumes seq is DNA, and returns non-ACTG letters as-is*/
    string reverse_complement(string seq);

    /* Capitalize all characters in a string */
    string to_upper(string seq);

    /* Returns the forward and reverse-reverse complement kmers of a sequence */
    vector<string> kmerize(string seq, int k);

    /* Returns the forward and reverse-reverse complement kmers for all kmer sizes in k */
    vector<string> multi_kmerize(string seq, vector<int> k);

    /* Returns a deduplicated set of string kmers */
    vector<string> kmer_set(vector<string> kmers);

    /* Returns a heap (priority queue) of the kmers of a read converted to ints. */
    
    /* Returns a heap (priority queue) of the kmers of the read */
    priority_queue<string> kmer_heap(string seq, vector<int> k);

    /* Converts a string kmer to an integer representation */
    int64_t kmer_to_integer(string kmer);


    /* Returns a deduplicated set of kmers or hashes as a vector<T> */
    template<typename T>
    inline vector<T> v_set(vector<T> kmers){
        set<T> s = set<T>(kmers.begin(), kmers.end());
        vector<T> ret = vector<T>(s.begin(), s.end());
        return ret;
    }

    /* Returns the forward shingles size k of a sequence */
    vector<string> shingle(string seq, int k);

    /* Returns the forward shingles of all k sizes of a sequence */
    vector<string> multi_shingle(string seq, vector<int> k);

    /* Return all hashes, unsorted, as fast as possible */
    vector<int64_t> allhash_unsorted_64(string& seq, vector<int>& k);

    /* Returns the lowest hashSize hashes of the kmers (length k...k` in k) of seq */
    vector<int64_t> minhash_64(string& seq, vector<int>& k, int hashSize, bool useBottom=true);

    /* Returns the bottom/top hashSize hashes of kmers size k in seq */ 
    vector<int64_t> minhash_64(string seq, int k, int hashSize, bool useBottom=true);

    /* Returns the bottom/top hashSize hashes of kmers size k which 
     * occur more than minDepth times, based on the depth in hash_to_depth */
    vector<int64_t> minhash_64_depth_filter(string& seq, vector<int>& k,
                            int hashSize, bool useBottom, int minDepth,
                            unordered_map<int64_t, int>& hash_to_depth);
    /* Takes in a list of pre-computed hashes and returns the MinHash (size hashSize)
     * of the hashes that pass the depth filter */
    vector<int64_t> minhash_64_depth_filter(vector<int64_t>& hashes, int hashSize, bool useBottom,
                                    int min_depth, unordered_map<int64_t, int>& hash_to_depth);
    /* helper function: returns the top hashSize hashes of the kmers size k in seq */
    vector<int64_t> top_minhash_64(string seq, int k, int hashSize);

    /* helper function: returns the bottom hashSize hashes of the kmers size k in seq */
    vector<int64_t> bottom_minhash_64(string seq, int k, int hashSize);

    /* Returns the union of the hashes in alpha and beta, including duplicates */
    vector<int64_t> hash_union(vector<int64_t> alpha, vector<int64_t> beta);

    /* Returns the intersection of alpha and beta, including duplicates the number of times they appear in both vectors */
    vector<int64_t> hash_intersection(vector<int64_t> alpha, vector<int64_t> beta);

    vector<string> kmer_intersection(vector<string> alpha, vector<string> beta);

    /* Returns the union of the two sets after deduplicating all duplicates */
    vector<int64_t> hash_set_union(vector<int64_t> alpha, vector<int64_t> beta);

    /* Returns the intersection of both sets. Duplicates are included only once */
    vector<int64_t> hash_set_intersection(vector<int64_t> alpha, vector<int64_t> beta);

    priority_queue<string> kmer_heap_intersection(priority_queue<string> alpha, priority_queue<string> beta);
}

#endif
