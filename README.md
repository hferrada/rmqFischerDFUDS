# rmqFischerDFUDS
This code creates a data structure to compute, in constant time, the RMQ(i,j) on an array of long integers.  The size for array of considerable length is around of 2.3n bits.

Authors: Hector Ferrada and Gonzalo Navarro 
hferrada@dcc.uchile.cl, gnavarro@dcc.uchile.cl

Description:<br />
This is an RMQ compressed data structure. The implementation is based on the method of Fischer and Heun [1]. 
Their method uses a tree representation with DFUDS. We use Range Min-Max Tree of Sadakane and Navarro [2]. 
In order to reduce the size, our rmq uses a simplified version of the Range min-max tree of Navarro and Sadakane [2], 
we used the forward minimum array to store the ranges (only for forward) and maximum backward array (only backward). 
We included an small example ("rmqFischerDFUDS.cpp") to show how to use it and included an small experiment for the time. 
This works only for 64 bits and supporting long sequences.

Make: 
To make the library just give the command 'make', this will create the lib: 'rmqFischer.a'.

Compile: 
To use the library you must compile your program linking 'rmqFischer.a' and include the file "DFUDSrmq.h" in your sorce code. 
For example to compiling the file rmqFischerDFUDS.cpp (included here) we will run: 
g++ rmqFischerDFUDS.cpp -o rmqFischerTest -O3 rmqFischer.a 
or simply run the command 'make test'. it will create the binary 'rmqFischerTest'. 
This binary have to recieve 7 parameter: 
1.- n: the length of random sequence 
2.- 0/1: to load from a file(0) or create(1) the complete structure
3.- saveLoadFile: the file (the path will be included) to store the data structure in order that you can load this later. 
4.- repetitions: number of repetitions for experiments
5.- pseudoSorted: 1 indicated that the input array will be encrease pseudo-sorted
6.- RandomWeight: if the file is in mode pseudo-sorted, then A[i] will be a arondom value in the range [i-weight, i+weight]
7.- resultsFile: the file to store the esperiments's results

For example, this line execute the code for n=10^4, stores the data in './rmqF-data.rmq' and the results in rmq-fischerRandom.txt: 
./rmqFischerTest 10000 1 rmqF-data.rmq 2000000 0 10000 rmq-fischerRandom.txt

References: Please, if you want to include this tool as part of your experiments, in your references, please to include the two papers above. 
Later, it will appear another publication to replace these ones.

[1]. J. Fischer and V. Heun. Space-efficient preprocessing schemes for range minimum queries on static arrays. 
SIAM Journal on Computing, 40(2):465–492, 2011.

[2]. K. Sadakane and G. Navarro. Fully-Functional Static and Dynamic Succinct Trees. 
ACM Transactions on Algorithms 10(3):article 16, 2014
