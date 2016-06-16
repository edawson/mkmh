#include <vector>

/* Returns the forward and reverse-reverse complement kmers of a sequence */
vector<string> kmerize(string seq, int k);

vector<string> multi_kmerize(string seq, vector<int> k);

/* Returns the forward shingles size k of a sequence */
vector<string> shingle(string seq, int k);

vector<string> multi_shingle(string seq, vector<int> k);

vector<int64_t> top_minhash_64(string seq, int k, int hashSize);

vector<int64_t> bottom_minhash_64(string seq, int k, int hashSize);

vector<int64_t> multi_bottom_minhash_64(string seq, vector<int> kSizes, int hashSize);
