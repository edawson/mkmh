#include <vector>
#include <set>
#include <string>
#include <sstream>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include "murmur3.hpp"


namespace mkmh{
    using namespace std;

    string reverse(string seq);

    string reverse_complement(string seq);

    /* Returns the forward and reverse-reverse complement kmers of a sequence */
    vector<string> kmerize(string seq, int k);

    vector<string> multi_kmerize(string seq, vector<int> k);

    vector<string> kmer_set(vector<string> kmers);

    /* Returns the forward shingles size k of a sequence */
    vector<string> shingle(string seq, int k);

    vector<string> multi_shingle(string seq, vector<int> k);

    vector<int64_t> minhash_64(string seq, vector<int> k, int hashSize, bool useBottom=true);

    vector<int64_t> minhash_64(string seq, int k, int hashSize, bool useBottom=true);

    vector<int64_t> top_minhash_64(string seq, int k, int hashSize);

    vector<int64_t> bottom_minhash_64(string seq, int k, int hashSize);

    vector<int64_t> hash_union(vector<int64_t> alpha, vector<int64_t> beta);

    vector<int64_t> hash_intersection(vector<int64_t> alpha, vector<int64_t> beta);
}
