I upload a prototype of keccak256 in C here, which fully reproduce what Binius does.

Actually, each algorithm used for ZKP-Binius64 coresponds to a prototype, or it is called "circuit" in Binius. To compile and run, commend as 

    
    gcc -O0 -g keccak256_witness.c -o state

