#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "mkmh.hpp"
#include <fstream>
#include <string>


using namespace std;
using namespace mkmh;
TEST_CASE("Reverse complement function works", "[reverse_complement]"){
    string t = "ACTGGCC";
    string rev = reverse_complement(t);
    //REQUIRE(rev == "GGCCAGT");

    char k [6] = "AGGTC";
    char* ret = new char[5];
    char* retret = new char[5];
    reverse_complement(k, ret, 5);
    REQUIRE(strcmp(ret, "GACCT") == 0);
    reverse_complement(ret, retret, 5);
    REQUIRE(ret == ret);
    REQUIRE(*retret == *k);


}

TEST_CASE("Canonical function catches non-valid chars", "[canonical]"){
    string t = "ACTGGCNNNN";
    REQUIRE(canonical(t) == false);

    char k [27] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    REQUIRE(canonical(k, 26) == false);

    char o [8] = "ACCCCTG";
    REQUIRE(canonical(o, 7) == true);


    char low[8] = "acccctg";
    REQUIRE(canonical(low, 7) == true);
}

TEST_CASE("Upper works for strings and chars"){

}

TEST_CASE("v_set removes duplicates and returns a vector", "[v_set]"){
    
}

TEST_CASE("kmerize works as expected for strings", "[kmerize(string, ...)]"){

}

TEST_CASE("kmerize functions for char* work as expected", "[kmerize(char*, ...)]"){

}

TEST_CASE("minimizers behave as expected", "[minimizers]"){

}

TEST_CASE("Calc_hashes functions produce the right hashes", "[calc_hashes]"){
    char  o [8] = "ACCCCTG";

    SECTION("Hashes for char* are consistent"){
    hash_t* h;
    int numhashes;
    calc_hashes((const char*) o, 7, 4, h, numhashes);
    
    bool trip = false;
    for (int i = 0; i < 7 - 4; i++){
        trip = *(h + i) != calc_hash(o + i, 4);
    }
    REQUIRE(trip == false);
    
    }
    
    
    vector<hash_t> x = calc_hashes((const char*) o, 7, 4);
    for (auto y : x){
        cout << y << endl;
    }
}

TEST_CASE("Calc_hash family of functions produce the right hash", "[calc_hash()]"){

}
