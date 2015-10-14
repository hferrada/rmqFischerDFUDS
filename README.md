# rmqFischerDFUDS
This code creates a data structure to compute, in constant time, the RMQ(i,j) on an array of long integers.  The size for array of considerable length is around of 2.3n bits.

Authors: Hector Ferrada and Gonzalo Navarro 
hferrada@dcc.uchile.cl, gnavarro@dcc.uchile.cl

Description:<br />
This is an RMQ compressed data structure. The implementation is based on the method of Fischer and Heun [1].<br />
Their method uses a tree representation with DFUDS. We use Range Min-Max Tree of Sadakane and Navarro [2].<br />
In order to reduce the size, our rmq uses a simplified version of the Range min-max tree of Navarro and Sadakane [2], we only used the forward minimum array to store the ranges (only for forward).<br />
We included an small example ("rmqFischerDFUDS.cpp") to show how to use it and included an small experiment for the time. This works only for 64 bits and supporting long sequences.

Make:<br />
To make the library just give the command 'make', this will create the lib: 'rmqFischer.a'.

Compile:<br />
To use the library you must compile your program linking 'rmqFischer.a' and include the file "DFUDSrmq.h" in your sorce code.<br />
For example to compiling the file rmqFischerDFUDS.cpp (included here) we will run:<br />
g++ rmqFischerDFUDS.cpp -o rmqFischerTest -O3 rmqFischer.a or simply run the command 'make test'. it will create the binary 'rmqFischerTest'. This binary have to recieve 7 parameter:<br />
1.- n: the length of random sequence.<br />
2.- 0/1: to load from a file(0) or create(1) the complete structure.<br />
3.- saveLoadFile: the file (the path will be included) to store the data structure in order that you can load this later.<br />
4.- repetitions: number of repetitions for experiments.<br />
5.- pseudoSorted: 1 indicated that the input array will be encrease pseudo-sorted.<br />
6.- RandomWeight: if the file is in mode pseudo-sorted, then A[i] will be a arondom value in the range [i-weight, i+weight].<br />
7.- resultsFile: the file to store the esperiments's results.<br />

For example, this line execute the code for n=10^4, stores the data in './rmqF-data.rmq' and the results in rmq-fischerRandom.txt:<br />
./rmqFischerTest 10000 1 rmqF-data.rmq 2000000 0 10000 rmq-fischerRandom.txt

References:<br />
Please, if you want to include this tool as part of your experiments, in your references, please to include the two papers above. Later, it will appear another publication to replace these ones.<br />

[1]. J. Fischer and V. Heun. Space-efficient preprocessing schemes for range minimum queries on static arrays. 
SIAM Journal on Computing, 40(2):465–492, 2011.<br />
[2]. K. Sadakane and G. Navarro. Fully-Functional Static and Dynamic Succinct Trees. 
ACM Transactions on Algorithms 10(3):article 16, 2014.
