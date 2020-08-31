# Circuit Amortization Friendly Encodings (CAFEs)

To run benchmarks, first compile using the Makefile using

    make shamir-ds-bm

This will compile benchmarking for SIMD double shares, FLEX double shares, and "regular" double shares (referred to as "naive" in the paper).

The script benchmark-ds.sh can be used to run the benchmarks with variable number of parties

    ./benchmark-ds.sh [program] [first ID] [last ID] [numparties] [network file]

Example network files (e.g. network_info5.txt, network_info7.txt, ...) are given.
