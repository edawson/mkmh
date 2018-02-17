# mkmh
Make kmers, minimizers, hashes, and MinHash sketches (with multiple k), and compare them. 

## Usage
To use mkmh functions in your code:  
1. Include the header file in your code  
    ```#include "mkmh.hpp"```      
2. Compile the library:  
    `` cd mkmh && make lib``  
3. Make sure the lib and header are on the LD include/lib paths (e.g. in your makefile):  
    `` gcc -o my_code my_code.cpp -L/path/to/mkmh -I/path/to/mkmh -lmkmh  
4. That's it!

## Available functions
Convenience funtions:

Get the shingles / kmers / minimizers / hashes of a string:

Compare sets of shingles / kmers / minimizers / hashes:

Fun extras:

## Getting help
Please reach out through [github](https://github.com/edawson/mkmh) to post an issue,
shoot me an email, or file a bug report.
