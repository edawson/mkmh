#include "mkmh.hpp"
#include <string>
#include <iostream>
#include <queue>

using namespace std;
using namespace mkmh;

bool testify(int t_num, string test, string obs, string exp){
    if (obs == exp){
        cout << "PASS: " << t_num << " " << test << endl;
        return true;
    }
    else{
        cout << "FAIL: " << t_num << " " << test << endl;
        return true;
    }
}

bool testify(int t_num, string test, bool x){
    if (x){
        testify(t_num, test, "", "");
        return true;
    }
    else{
        testify(t_num, test, "1", "2");
        return false;
    }   
}

bool same(vector<int64_t> x, vector<int64_t> y){
    if (x.size() != y.size()){
        return false;
    }
    
    std::sort(x.begin(), x.end());
    std::sort(y.begin(), y.end());
    for (int i = 0; i < x.size(); i++){
        if (x[i] != y[i]){
            return false;
        }
    }
    return true;

}

int main(){
    int t_num = 0;

    string seq = "ATGCATGCATGCATGCATGC";
    vector<string> kmers = kmerize(seq, 5);
    bool x = true;
    for (auto e : kmers){
        if (e.length() != 5){
            x = false;
        }
    }

    /* 
     *
     * Test kmerize, shingle, and minhash64
     * 
     * */
    testify(t_num++, "The kmers produced by kmerize are the right size", x);

    vector<string> shingles = shingle(seq, 5);
    x = true;
    for (auto e : shingles){
        if (e.length() != 5){
            x = false;
        }
    }
    testify(t_num++, "The shingles produced by kmerize are the right size", x);

    vector<int64_t> ret = minhash_64(seq, 5, 5);
    testify(t_num++, "minhash_64 produces the right number of hashes", ret.size() == 5);
    vector<int64_t> o_ret = minhash_64(seq, 5, 5);
    testify(t_num++, "minhash_64 produces the same top and bottom values in the hash", (ret[0] == o_ret[0] & ret[4] == o_ret[4]));
    
    vector<int> ret_test = {5};
    ret = minhash_64(seq, ret_test, 5);
    testify(t_num++, "multi_kmer minhash_64 produces the right number of hashes", ret.size() == 5);
    for (int iii = 1; iii < 1000; iii++){
        continue;
    }
    o_ret = minhash_64(seq, ret_test, 5);
    testify(t_num++, "multi_kmer minhash_64 produces the same top and bottom values in the hash", (ret[0] == o_ret[0] & ret[4] == o_ret[4]));


    vector<string> k_set = kmer_set(kmers);
    testify(t_num++, "Kmer set removes duplicate kmers", k_set.size() < kmers.size());

    /* Test top_minhash_64 and bottom_minhash_64 */
    x = true;
    vector<int64_t> tops = top_minhash_64(seq, 4, 4);
    vector<int64_t> bottoms = bottom_minhash_64(seq, 4, 4);
    for (auto e : tops){
        for (auto f : bottoms){
            if (e < f){
                x = false;
            }
        }
    }
    for (int i = 0; i < bottoms.size(); i++){
        //cerr << kmerize(seq, 4)[i] << endl;
        //cerr << bottoms[i] << " " << tops[i] << endl;
    }
    testify(t_num++, "top_minhash64 produces the bigger values than bottom_minhash_64", x);


    string seq2 = "ACTGaaatttt";
    vector<int> ks;
    ks.push_back(4);
    ks.push_back(4);
    kmers = kmerize(seq2, 4);
    vector<string> m_kmers = multi_kmerize(seq2, ks);

    /* Test multi_kmerize */
    testify(t_num++, "multi_kmerize produces twice as many kmers with two kmer sizes", kmers.size() * 2 == m_kmers.size());

    /* Test hash union / intersection */
    vector<int64_t> t1 = {1, 2, 3, 4, 5, 6};
    vector<int64_t> t2 = {4, 5, 6, 7, 8, 9};
    testify(t_num++, "Hash intersection of two sets is the expected size.", hash_intersection(t1, t2).size() == 3);

    t1 = {1, 2, 3, 4, 4, 4, 5, 6};
    t2 = {4, 4, 5, 6, 7, 8, 9};
    testify(t_num++, "Hash intersection counts duplicate values", hash_intersection(t1, t2).size() == 4);

    testify(t_num++, "Hash union yields the correct size with duplicate values", hash_union(t1, t1).size() ==   2 * t1.size());
    //testify(t_num++ "Union / Intersection produces expected value", hash_intersection(t1, t1).size() / hash_intersection(
    testify(t_num++, "Hash union yields the correct size when duplicates are removed", v_set(hash_union(t1, t1)).size() == 6);

    testify(t_num++, "Hash set union produces the expected number of values", hash_set_union(t1, t2).size() == 9);
    testify(t_num++, "Hash set intersection produces the expected number of values", hash_set_intersection(t1, t2).size() == 3);

    string a = "AAATGCTTTTGCA";
    vector<int> three;
    three.push_back(3);
    priority_queue<string> a_heap = kmer_heap(a, three);
    //    cerr << a_heap.top() << endl;
    testify(t_num++, "Kmer heap produces expected lowest kmer", (a_heap.top() == "GCA"));
    while (a_heap.size() > 1){
        //cerr << a_heap.top() << endl;
        a_heap.pop();
    }
    testify(t_num++, "Kmer heap produces expected highest kmer", (a_heap.top() == "AAA"));

    vector<int64_t> a_allhash_fast = allhash_unsorted_64_fast(a.c_str(), three);
    testify(t_num++, "allhash_unsorted_64_fast produces a hash vector of the proper length", (a_allhash_fast.size() == 10));
    for (auto ll : a_allhash_fast){
        cerr << ll << endl;
    }
    
    vector<int> four;
    four.push_back(4);

    vector<int64_t> bottoms_fast = minhash_64_fast(seq, four, 4, true);
    testify(t_num++, "minhash_64 and minhash_64_fast produce the same hashes", same(bottoms, bottoms_fast));

    return 0;
}
