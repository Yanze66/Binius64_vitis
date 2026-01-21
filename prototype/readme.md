An algorithm corresponds to a specific circuit used for the zk proof. The question is, it's a protocol efficient for hardware, but Binius64 still compile and execute in software (Rust), which doesn't take the advantages of hardware properties. As a result, I build really hardware circuit to perform Binius64.  


The logic is, any given algorithm ->  circuit of the alforithm (in Rust) -> prototype of the circuit (in C) -> build really hardware circuit.

This direction stalls the prototype of the circuit in C. I build a full implemented SHA3 primitive -- Keccak256 as an aexample.
