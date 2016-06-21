# mkmh
Generate and compare MinHash signatures with multiple kmer sizes.

## Usage
To use mkmh functions in your code:  
    1. Include the header file in your code  
    ```#include "mkmh.hpp"```  
    
    2. Compile the library:  
    `` cd mkmh && make lib``

    3. Make sure the lib and header are on the LD include/lib paths (e.g. in your makefile):  
    `` gcc -o my_code my_code.cpp -L/path/to/mkmh -I/path/to/mkmh -lmkmh

