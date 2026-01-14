Here is the basic idea of migerating the Rust project to hardware platform. 

# step 1: (Done)

I compile and execute the rust project of keccak and print out the witness as a golden reference. The main src file I wrok on is in dir "crates/circuits/src/keccak/fixed_length.rs", 
remark here so that anyone can reproduce the for the next. Also remember to install extensions "rust-analyzer" for a quick benchmark.


It supports to compress arbitrary length of message, recorded in "Keccak_witness_dump_xxbyte.txt".

# step 2: (under prosessing)

I try to migrate the whole rust project to the edge platform, and one of the challenge is to reproduce the whole project in C++. 

Based on the golden reference, I should first build keccak, and then find the key spot where to output witness.

# step 3

Build the hardware accelerator and implement an end-to-end Binius-64 implementation.
