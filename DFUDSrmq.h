/*
 * DFUDSrmq.h
 *
 *  Created on: 29-08-2015
 *      Author: hector
 */

#ifndef DFUDSRMQ_H_
#define DFUDSRMQ_H_z

#include "Basicrmq.h"

using namespace std;
using namespace dfudsrmq;


#define Srmq 256	// size of blocks (s bits each one), (power of 2 >= W)
#define PotSrmq 8	// power for block = log(Srmq)
#define SrmqD 512	// 2*Srmq
#define SrmqM 128	// Srmq/2;
#define N8Srmq 32 	// Srmq/8;
#define SaZe 512	// Sampling size for count zeros, it is the number of super blocks

//test values with small blocks ...
/*#define Srmq 64		// size of blocks (s bits each one), (power of 2 >= W)
#define PotSrmq 6	// power for block = log(Srmq)
#define SrmqD 128	// 2*Srmq
#define SrmqM 32	// Srmq/2;
#define N8Srmq 8 	// Srmq/8;*/

// FIXED VALUES:
#define SuBrmq 2	// number of leaves for each super block (power of 2 > RB)
#define BSrmq 8		// bits for each little block
#define BrmqMOne 7	// BSrmq minus one
const ulong RMMMasks[] = {0xFF00000000000000, 0x00FF000000000000, 0x0000FF0000000000, 0x000000FF00000000,
							 0x00000000FF000000, 0x0000000000FF0000, 0x000000000000FF00, 0x00000000000000FF,};

// -8 <= sum <= 8; // 2 bytes per cell --> 512 bytes
const short int T_SUM_BLOCK[] = {
		-8,-6,-6,-4,-6,-4,-4,-2,-6,-4,-4,-2,-4,-2,-2,0,
		-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,
		-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,
		-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
		-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,
		-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
		-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
		-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
		-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,
		-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
		-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
		-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
		-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
		-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
		-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
		0,2,2,4,2,4,4,6,2,4,4,6,4,6,6,8,
};

// 1 Bytes per cell --> 256 bytes
const uchar T_MIN_FWDI[] = {
		8,7,6,6,6,5,5,5,6,5,4,4,4,4,4,4,
		6,5,4,4,4,3,3,3,4,3,3,3,3,3,3,3,
		6,5,4,4,4,3,3,3,4,3,2,2,2,2,2,2,
		4,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
		6,5,4,4,4,3,3,3,4,3,2,2,2,2,2,2,
		4,3,2,2,2,1,1,1,2,1,1,1,1,1,1,1,
		4,3,2,2,2,1,1,1,2,1,1,1,1,1,1,1,
		2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
		6,5,4,4,4,3,3,3,4,3,2,2,2,2,2,2,
		4,3,2,2,2,1,1,1,2,1,1,1,1,1,1,1,
		4,3,2,2,2,1,1,1,2,1,0,0,0,0,0,0,
		2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		4,3,2,2,2,1,1,1,2,1,0,0,0,0,0,0,
		2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};

// 512 bytes
const short int T_MAX_BCKDI[] = {
		-8,1,0,2,-1,1,1,3,-2,1,0,2,0,2,2,4,
		-3,1,0,2,-1,1,1,3,-1,1,1,3,1,3,3,5,
		-4,1,0,2,-1,1,1,3,-2,1,0,2,0,2,2,4,
		-2,1,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
		-5,1,0,2,-1,1,1,3,-2,1,0,2,0,2,2,4,
		-3,1,0,2,-1,1,1,3,-1,1,1,3,1,3,3,5,
		-3,1,0,2,-1,1,1,3,-1,1,1,3,1,3,3,5,
		-1,1,1,3,1,3,3,5,1,3,3,5,3,5,5,7,
		-6,1,0,2,-1,1,1,3,-2,1,0,2,0,2,2,4,
		-3,1,0,2,-1,1,1,3,-1,1,1,3,1,3,3,5,
		-4,1,0,2,-1,1,1,3,-2,1,0,2,0,2,2,4,
		-2,1,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
		-4,1,0,2,-1,1,1,3,-2,1,0,2,0,2,2,4,
		-2,1,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
		-2,1,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
		0,2,2,4,2,4,4,6,2,4,4,6,4,6,6,8,
};

// 1 Bytes per cell --> 256 bytes
const uchar PT_MIN_FWDI[] = {
		7,6,5,5,7,4,4,4,7,6,3,3,3,3,3,3,
		7,6,5,5,7,2,2,2,7,2,2,2,2,2,2,2,
		7,6,5,5,7,4,4,4,7,6,1,1,1,1,1,1,
		7,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
		7,6,5,5,7,4,4,4,7,6,3,3,3,3,3,3,
		7,6,5,5,7,0,0,0,7,0,0,0,0,0,0,0,
		7,6,5,5,7,0,0,0,7,0,0,0,0,0,0,0,
		7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		7,6,5,5,7,4,4,4,7,6,3,3,3,3,3,3,
		7,6,5,5,7,2,2,2,7,2,2,2,2,2,2,2,
		7,6,5,5,7,4,4,4,7,6,1,1,1,1,1,1,
		7,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
		7,6,5,5,7,4,4,4,7,6,3,3,3,3,3,3,
		7,6,5,5,7,0,0,0,7,0,0,0,0,0,0,0,
		7,6,5,5,7,0,0,0,7,0,0,0,0,0,0,0,
		7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};

// 1 Bytes per cell --> 256 bytes
const uchar POPC0[] = {
		8,7,7,6,7,6,6,5,7,6,6,5,6,5,5,4,
		7,6,6,5,6,5,5,4,6,5,5,4,5,4,4,3,
		7,6,6,5,6,5,5,4,6,5,5,4,5,4,4,3,
		6,5,5,4,5,4,4,3,5,4,4,3,4,3,3,2,
		7,6,6,5,6,5,5,4,6,5,5,4,5,4,4,3,
		6,5,5,4,5,4,4,3,5,4,4,3,4,3,3,2,
		6,5,5,4,5,4,4,3,5,4,4,3,4,3,3,2,
		5,4,4,3,4,3,3,2,4,3,3,2,3,2,2,1,
		7,6,6,5,6,5,5,4,6,5,5,4,5,4,4,3,
		6,5,5,4,5,4,4,3,5,4,4,3,4,3,3,2,
		6,5,5,4,5,4,4,3,5,4,4,3,4,3,3,2,
		5,4,4,3,4,3,3,2,4,3,3,2,3,2,2,1,
		6,5,5,4,5,4,4,3,5,4,4,3,4,3,3,2,
		5,4,4,3,4,3,3,2,4,3,3,2,3,2,2,1,
		5,4,4,3,4,3,3,2,4,3,3,2,3,2,2,1,
		4,3,3,2,3,2,2,1,3,2,2,1,2,1,1,0,
};
class DFUDSrmq {
	typedef struct StackDFUDS{
		long int val;
		StackDFUDS *next;
	} StackDFUDS;
private:
	ulong nW;				// number of words to store P
	ulong *P;				// sequence of 2n' balanced parentheses, which represents a tree with n/2 nodes.

	ulong rank1_Bin;		// number of 1's from 0 to nBin
	ulong nBin;				// 0 < nBin <= n. This determine the partition of the interval for the Binary Tree and the last block
							// 'n - nBin' is the number of bits that must be considered as a single last block of length 'n - nBin'

	uint lenLB;				// length of lastBlock (number of bits)
	uint bitsSuB, bitsRB;	// by default bitsSuB=256 and bitsRB=128
	uint h;					// nim-max tree's height

	ulong cantN;			// total nodes = leaves + cantIN
	ulong cantIN;			// number of internal nodes of the min-max tree (nodes = cantIN + leaves)
	ulong leaves;			// number of leaves of min-max tree
	ulong leavesBottom;		// number of leaves in the last level h (perhaps there are leaves in the level h-1 too)
	ulong firstLeaf;		// position of the first leaf (the left-most leaf)

	ulong *Fwd_MinIN;		// the minimum excess value for internal nodes (stored as positive number)
	int MIN_Fwd;			// the lowest value for forward intervals
	uint lgMIN_Fwd;

	ulong *TSBlock;			// Table of global excess for each superblock in P (groups of k blocks). Size = (n/ks)*MAXEXC
	ulong lenSB;			// number of super/relative blocks
	ulong zeSB;				// number of 0 in all SB

	ulong *SZ;				// Table of Sampling count for zeros
	uint MAX_SupB;
	uint lgMAX_SupB;
	char *TRBlock;			// Table of excess relative only for the first block in each super block.
							// TRBlock[i] = sum_i/2. where, sim_i is the relative sum for the first block of the superblock i,
							// also -255 <= sum_i <= 255 and sum_i is an even number. We need 8 bits for each value.
	ulong *Bfull;			// this indicates the leaves which contain in TRBlock the numbers -256 or 256
	uchar *TPMinB;			// Table of minimum positions for each leaf. values between 0 and 255
	ulong *TMinB;			// Table of leaf minimum

	uint MAX_B;				// the greater global excess for all block.
	uint lgMAX_B;

	uint sizeDS;			// in bytes

public:
	ulong nP;				// Length of sequence P (n parentheses and n/2 nodes)

	static bool TRACE;		// true: print all details for console
	static bool RUNTEST;
	static uint TEST;

	DFUDSrmq(long int *A, ulong len);
	DFUDSrmq(char *fileName);

	void createMinMaxTree();
	void createTables();
	void printTree();

	ulong binRank_1(ulong i);
	ulong rank_1(ulong i);

	ulong binSelect_1(ulong i);
	ulong select_1(ulong i);

	ulong binSelect_0(ulong i);
	ulong select_0(ulong i);
	ulong select_0_old(ulong i);

	// return the position of the "(" for the ")" at position i
	ulong open_0(ulong i);

	// give the excess from 0 to pos
	long int sumAtPos(ulong pos);

	// give the excess of the internal node 'node=preorder+1' that has a distance 'dist' to the tree's depth
	long int computeSumOfNode(ulong node, ulong dist);

	ulong computeLeavesOfNode(ulong node, ulong dist);
	ulong cantLeavesOfNode(ulong node);

	// give the excess of the internal node 'node=preorder+1' that has a distance 'dist' to the tree's depth
	void search_min_block(ulong x1, ulong x2, long int *min, long int *curSum, ulong *position);

	// return the position of the open parenthesis closet to the root between i and j
	ulong rmqi(ulong i, ulong j);
	ulong rmqi_rmm(ulong x1, ulong x2, long int *min, long int *currSum, ulong posMin);

	bool backward_search_block(ulong x1, ulong x2, long int exc, long int *sumPos, ulong *pos);
	ulong backward_search_rmm(long int exc, ulong j);
	ulong backward_search(long int exc, ulong j);
	void test_backward_search_block();
	void test_backward_search();

	// query for RMQ
	ulong queryRMQ(ulong i, ulong j);

	uint getSize();

	// save the Data Structure in file 'fileName'
	void saveDS(char *fileName);

	// load the Data Structure from the file 'fileName'
	void loadDS(char *fileName);

	virtual ~DFUDSrmq();

	void test_rank_1();
	void test_select_1();
	void test_select_0();
	void test_sumAtPos();
	void test_rmqi();
	void test_search_min_block();
	void test_findclose();
};

#endif /* DFUDSRMQ_H_ */
