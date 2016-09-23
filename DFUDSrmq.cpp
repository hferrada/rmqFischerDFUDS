/*
 * DFUDSrmq.cpp
 *
 *  Created on: 29-08-2015
 *      Author: hector
 *
 *
 *      The 2n bits OF TOPOLOGY ARE INCLUDED  !!
 *
 */

#include "DFUDSrmq.h"

bool DFUDSrmq::TRACE = false;
bool DFUDSrmq::RUNTEST = false;
uint DFUDSrmq::TEST = 10000;

DFUDSrmq::DFUDSrmq(char *fileName){
	loadDS(fileName);

	if (RUNTEST){
		test_sumAtPos();
		test_rank_1();
		test_select_1();
		test_select_0();
		test_search_min_block();
		test_backward_search_block();
		test_backward_search();
		if (leaves>0) test_rmqi();
	}
}

DFUDSrmq::DFUDSrmq(long int *A, ulong len) {
	if (TRACE){
		cout << "*** RMQ structure for A[0.." << len-1 << "] and leaf length S = " << Srmq << endl;
		cout << "Create BP sequence for 2D Min Heap..." << endl;
	}

	sizeDS = 2*512 + 3*256;										// size for T_SUM_BLOCK[] + T_MIN_FWDI[] + T_MAX_BCKDI[] + PT_MIN_FWDI[] + POPC0[]
	sizeDS += 11*sizeof(ulong) + 10*sizeof(uint) + sizeof(int);	// size for variables

	nP = (len+1)<<1;
	ulong lenP = nP >> BW64;
	if (nP % W64)
		lenP++;
	P = new ulong[lenP];
	sizeDS += lenP*sizeof(ulong);		// 2n bits OF TOPOLOGY ARE INCLUDED  !!
	if (true || TRACE) cout << " ** size of topology " << lenP*sizeof(ulong) << " Bytes" << endl;

	cout << "Creating DFUDS sequence in P[1..2*len]" << endl;
	ulong k, pos, items, childRoot=1;
	long int i;

	pos=nP-1;
	items = 0;
	StackDFUDS *Q = new StackDFUDS();
	StackDFUDS *R;

	for(i=len-1; i>=0; i--){
		cleanBit64(P, pos);	 	// close parenthesis for A[i]
		pos--;

		// put a opening parenthesis for each value stored in Q, such that A[i]>Q->val, deleting these values from Q.
		while(items && A[Q->val] > A[i]){
			setBit64(P, pos);
			pos--;
			R = Q;
			Q = Q->next;
			delete R;
			items--;
		}

		// insert A[i] into Q
		R = new StackDFUDS();
		R->val = i;
		R->next = Q;
		Q = R;
		items++;
	}
	cleanBit64(P, pos);	 	// close parenthesis for -infinity
	pos--;
	// put a opening parenthesis for each value stored in Q, such that A[i]>Q->val, deleting these values from Q.
	while(items){
		setBit64(P, pos);
		pos--;
		R = Q;
		Q = Q->next;
		delete R;
		items--;
	}
	setBit64(P, pos);		// additional open parenthesis for -infinity

	cout << "DFUDS TREE REPRESENTATION OK !!" << endl;

	if(RUNTEST){ // test for balanced sequence
		long int sum = 0;
		ulong r=0;

		for (; r<nP; r++){
			if(readBit64(P, r))
				sum++;
			else sum--;
		}
		if(sum != 0){
			cout << " ERROR. DFUDS-P[] is not a balanced sequence of parentheses !! " << endl;
			exit(1);
		}else
			cout << " DFUDS-P[] is a well balanced sequence of parentheses !! " << endl;
	}

	if (TRACE){
		cout << " DEFUDS sequence, P[0.." << nP-1 << "]" << endl;
		for (pos=0; pos<nP; pos++){
			if(readBit64(P, pos))
				cout << "(";
			else
				cout << ")";
		}
		cout << endl;
	}

	ulong nb = nP / Srmq;
	if (nb%2 && nb > 1)
		nb--;
	nBin = nb*Srmq;
	nW = nP/W64;
	if (nP%W64)
		nW++;

	if (nBin >= nP){
		nBin = nP;
		lenLB = 0;
	}else
		lenLB = nP - nBin;

	if (TRACE)
		cout << "Create RangeMinMaxTree_Bin with length N = 2n = " << nP << ", nBin: " << nBin << ", nW: " << nW << endl;

	this->leaves = nBin/Srmq;
	this->h = ceilingLog64(this->leaves, 2);
	this->firstLeaf = (ulong)pow(2, (double)h) - 1;	// the left-most leaf

	if (TRACE){
		cout << "___________________________________________________________" << endl;
		cout << "P_bin :";
		for (i=0; i<nBin; i++){
			cout << readBit64(P, i);
			if ((i+1)%Srmq == 0)
				cout << "-";
		}
		cout << endl << "Last Block :" << endl;
		for (; i<nP; i++)
			cout << readBit64(P, i);
		cout << endl;
	}

	createMinMaxTree();
	if (TRACE){
		if (TRACE) cout << " *** TOTAL DFUDSrmq SIZE: " << sizeDS << " Bytes = " << (float)sizeDS/(1024.0*1024.0) << " MB." << endl;
		if (leaves>0) printTree();
		cout << "# blocks " << nb << endl;
		cout << "# blocks of the binary tree " << nb << endl;
		cout << "rank1_Bin (1's) " << rank1_Bin << endl;
		cout << "last bit of binary tree: " << nBin << endl;
		cout << "Length of the last block " << lenLB << endl;
		cout << "Sum relative to binary tree: " << sumAtPos(nBin) << endl;
		cout << "======================"<< endl;
	}
	if (RUNTEST){
		test_sumAtPos();
		test_rank_1();
		test_select_1();
		test_select_0();
		test_search_min_block();
		test_backward_search_block();
		test_backward_search();
		if (leaves>0) test_rmqi();
	}
}

void DFUDSrmq::createMinMaxTree(){
	ulong groups, leavesUp, sizeAux;
	ulong i, j, position, father, node, cont, child, rb;
	ulong *auxL, *auxR;
	int *Aux_Fwd_MinIN;

	if (h>0){
		groups = (ulong)pow(2, (double)h-1.0);	// number of nodes at level (h-1)
		leavesBottom = (leaves - groups)<<1;
		leavesUp = leaves - leavesBottom;
	}else{
		leavesBottom = leaves;
		leavesUp = 0;
	}

	firstLeaf = (ulong)pow(2, (double)h) - 1;
	cantIN = firstLeaf - leavesUp;
	cantN = cantIN + leaves;

	if (TRACE && leaves>0){
		cout << "leaves: " << leaves << endl;
		cout << "leaves Up: " << leavesUp << endl;
		cout << "leaves Bottom: " << leavesBottom << endl;
		cout << "internal nodes: " << cantIN << endl;
		cout << "first leaf: " << firstLeaf << endl;
		cout << "total nodes: " << cantN << endl;
	}

	if (leaves>0){
		auxL = new ulong[cantN];	// these are auxiliary vectors to facility the compute of intervals in each internal nodes. These will be delete later.
		auxR = new ulong[cantN];

		// Step 1: process the n bits in groups of size S, and create the leaves in each array:	intervals of relatives excess per block, in MinLeaf and MaxLeaf,
		//         and the relative sum of excess per block, in ELeaf
		// 'i' is for bits in P, 'cont' is for count the size of each block, and 'node' is the current leaf position to set min-max values...
		//cout << "Step 1: process the n bits in groups of size s, and create the leaves in each array..." << endl;
		j = firstLeaf; // this is for auxiliary vectors
		ulong bitBottom = leavesBottom*Srmq;
		for (i=node=cont=0; i<nBin; i++){
			if (i == bitBottom){
				node = leavesBottom;	// we move up one level
				j = cantIN;
			}
			if(cont==0){
				auxL[j] = i;
				position = j;
				// set left boundaries, when the index j is the first child.
				while(position && (position+1)%2==0){
					father = (position+1)>>1;
					if ((position+1)%2 > 1)
						father++;
					father--;
					auxL[father] = i;
					position = father;
				}
				cont=1;
			}else{
				cont++;
				if(cont == Srmq){
					auxR[j] = i;
					position = j;
					// set right boundaries, when the index j is the last child. The last child always has rest == 1.
					while(position && (position+1)%2==1){
						father = (position+1)>>1;
						if ((position+1)%2 > 1)
							father++;
						father--;
						auxR[father] = i;
						position = father;
					}
					node++;
					j++;
					cont = 0;
				}
			}
		}
		if (cont)
			auxR[node] = i;
	}

	// Step 2: create arrays of super/relative blocks...
	//cout << "Step 2: create arrays of Super/Relative blocks..." << endl;
	createTables();
	if (TRACE && leaves>0){
		cout << endl << "min-max Intervals..." << endl;
		for (i=0; i<cantN; i++)
			cout << "[" << auxL[i] << "," << auxR[i] << "] ";
		cout << endl;
	}

	if (leaves>0){
		int miniFwd;
		ulong segment;
		long int currSumBack, currSumFwd;
		Aux_Fwd_MinIN = new int[cantIN];
		//Aux_Back_MaxIN = new int[cantIN];
		MAX_B = MIN_Fwd = 0;
		// Step 3: set the relative min values for each internal node.
		//cout << "Step 3: set the relative min values for each internal node..." << endl;
		for(i=cantIN; i>0; i--){
			child = i<<1;		// child is the second child of node i
			if (child > cantIN){
				if (child > firstLeaf)
					node = child - firstLeaf;
				else
					node = child - cantIN + leavesBottom;
				currSumFwd = currSumBack = 0;
				miniFwd = pow(2,31)-1;

				// Forward...
				for (j=rb=0; j<N8Srmq; j++){
					//segment = (P[(node*Srmq-1-BSrmq*j)>>BW64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
					segment = (P[((node-1)*Srmq+BSrmq*j)>>BW64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
					if (currSumFwd - T_MIN_FWDI[segment] < miniFwd)
						miniFwd = currSumFwd - T_MIN_FWDI[segment];

					currSumFwd += T_SUM_BLOCK[segment];
					if (rb == N8W64-1) rb=0;
					else rb++;
				}
				for (j=0; j<N8Srmq; j++){
					//segment = (P[((node+1)*Srmq-1-BSrmq*j)>>BW64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
					segment = (P[(node*Srmq+BSrmq*j)>>BW64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
					if (currSumFwd - T_MIN_FWDI[segment] < miniFwd)
						miniFwd = currSumFwd - T_MIN_FWDI[segment];

					currSumFwd += T_SUM_BLOCK[segment];
					if (rb == N8W64-1) rb=0;
					else rb++;
				}

			}else{
				miniFwd = Aux_Fwd_MinIN[child-1];

				currSumFwd = sumAtPos(auxR[child-1]);
				if (auxL[child-1])
					currSumFwd -= sumAtPos(auxL[child-1]-1);
				if (currSumFwd + Aux_Fwd_MinIN[child] < miniFwd)
					miniFwd = currSumFwd + Aux_Fwd_MinIN[child];
			}
			Aux_Fwd_MinIN[i-1] = miniFwd;
			if(miniFwd<0){
				if (-1*miniFwd > MIN_Fwd)
					MIN_Fwd = -1*miniFwd;
			}else{
				if (miniFwd > MIN_Fwd)
					MIN_Fwd = miniFwd;
			}
		}
		delete [] auxL;
		delete [] auxR;
	}

	lgMIN_Fwd = ceilingLog64(MIN_Fwd+1, 2);
	if(lgMIN_Fwd==0)
		lgMIN_Fwd=1;
	if (TRACE){
		cout << "MIN_Fwd (*-1) = " << MIN_Fwd << ", lgMIN_Fwd = " << lgMIN_Fwd << endl;
		//cout << "MAX_Fwd (*-1) = " << MAX_Bck << ", lgMAX_Fwd = " << lgMAX_Bck << endl;
	}

	sizeAux = cantIN*lgMIN_Fwd/W64;
	if ((cantIN*lgMIN_Fwd)%W64)
		sizeAux++;
	Fwd_MinIN = new ulong[sizeAux];
	sizeAux *= sizeof(ulong);
	sizeDS += sizeAux;
	if (true || TRACE) cout << " ** size of Fwd_MinIN[] " << sizeAux << " Bytes" << endl;
	cont = 0;
	for(i=0; i<cantIN; i++){
		setNum64(Fwd_MinIN, cont, lgMIN_Fwd, -1*Aux_Fwd_MinIN[i]);
		cont += lgMIN_Fwd;
	}
	if (cantIN)
		delete [] Aux_Fwd_MinIN;
}

void DFUDSrmq::createTables(){
	ulong sizeAux;
	ulong i, j, jS, jR, cont, rb, segment, posMin;
	long int sum, sumBlock;
	int Min;
	lenSB = (nP/Srmq)/SuBrmq + 1;
	bitsSuB = Srmq*SuBrmq;
	ulong *AuxTSBlock = new ulong[lenSB];



	TRBlock = new char[lenSB];				// ???
	sizeAux = (lenSB-1)*sizeof(char);
	sizeDS += sizeAux;
	if (true || TRACE) cout << " ** size of TRBlock[] " << sizeAux << " Bytes" << endl;

	ulong lenWRel = (lenSB-1)/W64;
	if ((lenSB-1)%W64)
		lenWRel++;
	Bfull = new ulong[lenWRel];
	sizeAux = lenWRel*sizeof(ulong);
	sizeDS += sizeAux;
	if (true || TRACE) cout << " ** size of Bfull[] " << sizeAux << " Bytes" << endl;

	j = leaves/W64;
	if (j%W64)
		j++;
	TPMinB = new uchar[leaves];		// positions explicitly  ????
	sizeAux = leaves*sizeof(uchar);
	sizeDS += sizeAux;
	if (true || TRACE) cout << " ** size of TPMinB[] " << sizeAux << " Bytes" << endl;

	int *AuxTMinB = new int[leaves];
	cont = sum = 0;
	jR = 0;
	jS = 1;
	AuxTSBlock[0] = MAX_SupB = 0;
	// here.... always sum >= 0, because is a BP sequence
	rb = 0;
	for (i=0; jS<lenSB; i++){
		sumBlock = 0;
		for (j=0; j<N8Srmq; j++){
			segment = (P[(i*Srmq+BSrmq*j)/W64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
			//printBitsUlong(segment);cout<<endl;
			sumBlock += T_SUM_BLOCK[segment];
			if (rb == N8W64-1) rb=0;
			else rb++;
		}

		sum += sumBlock;
		if(cont==0){
			if (sumBlock==-Srmq || sumBlock==Srmq){
				setBit64(Bfull,jR);
				if(sumBlock==-Srmq)
					TRBlock[jR] = 0;
				else
					TRBlock[jR] = 1;
			}else{
				cleanBit64(Bfull,jR);
				TRBlock[jR] = (char)(sumBlock>>1);
			}
			jR++;
		}
		cont++;
		if (cont == SuBrmq){
			AuxTSBlock[jS] = sum;
			if (sum > (int)MAX_SupB)
				MAX_SupB = sum;
			jS++;
			cont = 0;
		}
	}
	rb = 0;
	if((Srmq+(lenSB-1)*bitsSuB-1)<nP){
		sumBlock = 0;
		for (j=0; j<N8Srmq; j++){
			segment = (P[(jR*bitsSuB+BSrmq*j)/W64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
			sumBlock += T_SUM_BLOCK[segment];
			if (rb == N8W64-1) rb=0;
			else rb++;
		}
		if (sumBlock==-Srmq || sumBlock==Srmq){
			setBit64(Bfull,jR);
			if(sumBlock==-Srmq)
				TRBlock[jR] = 0;
			else
				TRBlock[jR] = 1;
		}else{
			cleanBit64(Bfull,jR);
			TRBlock[jR] = (char)(sumBlock>>1);
		}
	}

	if (lenSB == 1){
		rank1_Bin = 0;
		for (j=0; j<nBin; j++){
			if(readBit64(P, j))
				rank1_Bin++;
		}
	}else
		rank1_Bin = ((nBin-sum)>>1) + sum;
	sum = 0;
	for (i=0; i<h; i++){
		cont = (ulong)pow(2, (double)(i));
		sum += cont;
	}

	MAX_B = 0;
	if (leaves){
		bool isMin;
		for (i=0; i<leaves; i++){
			sumBlock = posMin = 0;
			Min = pow(2,31)-1;
			isMin = false;

			cont = i*Srmq;
			for (j=0; j<Srmq; j++, cont++){
				if (readBit64(P, cont))
					sumBlock++;
				else{
					sumBlock--;
					if(sumBlock < Min){
						isMin = true;
						Min = sumBlock;
						posMin = cont%Srmq;
					}
				}
			}
			TPMinB[i] = posMin;
			if(isMin == false)
				AuxTMinB[i] = 0;
			else
				AuxTMinB[i] = Min;

			if (Min<0 && MAX_B<(uint)(-1*Min))
				MAX_B = -1*Min;
		}
		lgMAX_SupB = ceilingLog64(MAX_SupB+1, 2);
	}else
		lgMAX_SupB = 1;

	lgMAX_B = ceilingLog64(MAX_B+1, 2);

	cont = leaves*lgMAX_B/W64;
	if ((leaves*lgMAX_B)%W64)
		cont++;
	TMinB = new ulong[cont];
	sizeAux = cont*sizeof(ulong);
	sizeDS += sizeAux;
	if (true || TRACE) cout << " ** size of TMinB[] " <<  sizeAux << " Bytes" << endl;

	for (i=cont=0; i<leaves; i++){
		if (AuxTMinB[i] < 0)
			setNum64(TMinB, cont, lgMAX_B, (ulong)(-1*AuxTMinB[i]));
		else
			setNum64(TMinB, cont, lgMAX_B, 0);
		cont += lgMAX_B;
	}
	if (leaves)
		delete []AuxTMinB;

	cont = lenSB*lgMAX_SupB/W64;
	if ((lenSB*lgMAX_SupB)%W64)
		cont++;
	TSBlock = new ulong[cont];
	sizeAux = cont*sizeof(ulong);
	sizeDS += sizeAux;
	if (true || TRACE) cout << " ** size of TSBlock[] " <<  sizeAux << " Bytes" << endl;

	for (i=cont=0; i<lenSB; i++, cont+=lgMAX_SupB)
		setNum64(TSBlock, cont, lgMAX_SupB, AuxTSBlock[i]);

	zeSB = ((lenSB-1)*bitsSuB - AuxTSBlock[lenSB-1])>>1;

	if (lenSB)
		delete [] AuxTSBlock;

	if (TRACE){
		cout << "MAX_SupB " << MAX_SupB << ", lgMAX_SupB " << lgMAX_SupB << ", MAX_B " << MAX_B << ", lgMAX_B " << lgMAX_B << endl;
		cout << "TSBlock[1.." <<lenSB<< "]..." << endl;
		for (i=0; i<lenSB; i++)
			cout << getNum64(TSBlock, i*lgMAX_SupB, lgMAX_SupB) << " ";
		cout << endl;
		cout << "TRBlock[1.." <<lenSB-1<< "]..." << endl;
		for (i=0; i<lenSB-1; i++)
			cout << (int)TRBlock[i] << " ";
		cout << endl;
		cout << "Bfull[1.." <<lenSB-1<< "]..." << endl;
		for (i=0; i<lenSB-1; i++)
			cout << readBit64(Bfull, i);
		cout << endl;
		cout << "TPMinB[1.." <<leaves<< "]..." << endl;
		for (i=0; i<leaves; i++)
			cout << (uint)TPMinB[i] << " ";
		cout << endl;
		cout << "TMinB[1.." <<leaves<< "]..." << endl;
		for (i=0; i<leaves; i++)
			cout << getNum64(TMinB, i*lgMAX_B, lgMAX_B) << " ";
		cout << endl;
	}
}

// *******************************************************************************
// ********************************* BASIC OPERATIONS ****************************

// for 0 <= i < n
ulong DFUDSrmq::binRank_1(ulong i){
	return (sumAtPos(i)+i+1)>>1;
}

// for 0 <= i < n
ulong DFUDSrmq::rank_1(ulong i){
	if(i < nBin-1)
		return binRank_1(i);
	else{
		if (nBin-1 == i) return rank1_Bin;
		if(i >= nP-1) return nP>>1;

		ulong blk=i/bitsSuB;
		ulong x = blk*bitsSuB;
		ulong rank = (x + getNum64(TSBlock, blk*lgMAX_SupB, lgMAX_SupB))>>1;
		if((i-x) > Srmq){
			if(readBit64(Bfull, blk)){
				if(TRBlock[blk])
					rank += Srmq;
			}else
				rank += SrmqM + TRBlock[blk];
			x+=Srmq;
		}
		ulong b, rest, q;
		while(x+BrmqMOne <= i){
			b = x>>BW64;
			rest = (x+BSrmq)%W64;
			q = (P[b] >> (W64-rest)) & 0xff;
			rank += __popcount_tab[q];
			x += BSrmq;
		}

		// check last segment (< S) bit by bit...
		while (x<=i){
			if(readBit64(P,x))
				rank++;
			x++;
		}

		return rank;
	}
}

ulong DFUDSrmq::binSelect_1(ulong i){
	ulong blk, blk2, rank_B, raAux, q, rb, curr;

	// search on super blocks...
	blk = (i/bitsSuB)<<1;  // proportional search considering that it is a sequence balanced of 1's and 0's
	if (blk){
		raAux = (blk*bitsSuB + getNum64(TSBlock, blk*lgMAX_SupB, lgMAX_SupB))>>1;
		if(raAux<i){
			rank_B = raAux;
			blk2 = blk+1;
			raAux = (blk2*bitsSuB + getNum64(TSBlock, blk2*lgMAX_SupB, lgMAX_SupB))>>1;
			while(raAux<i){
				rank_B = raAux;
				blk = blk2;
				blk2++;
				raAux = (blk2*bitsSuB + getNum64(TSBlock, blk2*lgMAX_SupB, lgMAX_SupB))>>1;
			}
		}else{
			while(raAux >= i){		// I exceeded! come back...
				blk--;
				raAux = (blk*bitsSuB + getNum64(TSBlock, blk*lgMAX_SupB, lgMAX_SupB))>>1;
			}
			rank_B = raAux;
		}

		// search in blocks...
		if(readBit64(Bfull, blk)){
			if(TRBlock[blk])
				raAux = rank_B + Srmq;
			else
				raAux = rank_B;
		}else
			raAux = rank_B + TRBlock[blk] + SrmqM;

		curr = SrmqD*blk;
		if(raAux<i){
			rank_B = raAux;
			curr += Srmq;
		}
	}else{
		if(readBit64(Bfull, 0)){
			if(TRBlock[0])
				raAux = Srmq;
			else
				raAux = 0;
		}else
			raAux = SrmqM + TRBlock[0];

		if(raAux<i){
			rank_B = raAux;
			curr=Srmq;
		}else
			rank_B = curr = 0;
	}

	// search in the leaves... go from the left... it is possible go from the right !!
	q = (P[curr>>BW64] & RMMMasks[0]) >> W64m8;
	raAux = rank_B + ((BSrmq + T_SUM_BLOCK[q])>>1);
	rb = 0;
	while (raAux < i){
		rank_B = raAux;
		curr+=BSrmq;
		rb++;
		if (rb == N8W64)
			rb=0;
		q = (P[curr>>BW64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
		raAux += (BSrmq + T_SUM_BLOCK[q])>>1;
	}
	for (; rank_B<i; curr++){
		if (readBit64(P, curr))
			rank_B++;
	}

	return curr-1;
}

ulong DFUDSrmq::select_1(ulong i){
	ulong pos = 0;

	if(i <= rank1_Bin)
		return binSelect_1(i);
	else{
		ulong rb, l, q, curr;
		ulong rank_A, rank_B = rank1_Bin;

		curr = nBin+BrmqMOne;
		q = (P[curr>>BW64] & RMMMasks[0]) >> W64m8;
		rank_A = rank_B + ((BSrmq + T_SUM_BLOCK[q])>>1);
		rb = l = 0;
		while (rank_A < i){
			rank_B = rank_A;
			rb++;
			if (rb == N8W64)
				rb=0;
			curr+=BSrmq;
			q = (P[curr>>BW64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
			rank_A += (BSrmq + T_SUM_BLOCK[q])>>1;
		}
		pos=curr-BrmqMOne;
		for (; rank_B<i; pos++){
			if (readBit64(P, pos))
				rank_B++;
		}
	}

	return pos-1;
}

ulong DFUDSrmq::binSelect_0(ulong i){
	ulong blk, blk2, rank_B, raAux, q, rb, curr;

	// search on super blocks...
	blk = (i/bitsSuB)<<1;  // proportional search considering that it is a sequence balanced of 1's and 0's
	if (blk){
		raAux = (blk*bitsSuB - getNum64(TSBlock, blk*lgMAX_SupB, lgMAX_SupB))>>1;
		if(raAux<i){
			rank_B = raAux;
			blk2 = blk+1;
			raAux = (blk2*bitsSuB - getNum64(TSBlock, blk2*lgMAX_SupB, lgMAX_SupB))>>1;
			while(raAux<i){
				rank_B = raAux;
				blk = blk2;
				blk2++;
				raAux = (blk2*bitsSuB - getNum64(TSBlock, blk2*lgMAX_SupB, lgMAX_SupB))>>1;
			}
		}else{
			while(raAux >= i){		// I exceeded! come back...
				blk--;
				raAux = (blk*bitsSuB - getNum64(TSBlock, blk*lgMAX_SupB, lgMAX_SupB))>>1;
			}
			rank_B = raAux;
		}

		// search in blocks...
		if(readBit64(Bfull, blk)){ // si esta marcado como lleno
			if(TRBlock[blk])
				raAux = rank_B;
			else
				raAux = rank_B + Srmq;
		}else
			raAux = rank_B + SrmqM - TRBlock[blk];

		curr = SrmqD*blk;
		if(raAux<i){
			rank_B = raAux;
			curr += Srmq;
		}
	}else{
		if(readBit64(Bfull, 0)){
			if(TRBlock[0])
				raAux = 0;
			else
				raAux = Srmq;
		}else
			raAux = SrmqM - TRBlock[0];

		if(raAux<i){
			rank_B = raAux;
			curr=Srmq;
		}else
			rank_B = curr = 0;
	}

	// search in the leaves... go from the left... it is possible go from the right !!
	q = (P[curr>>BW64] & RMMMasks[0]) >> W64m8;
	raAux = rank_B + ((BSrmq - T_SUM_BLOCK[q])>>1);
	rb = 0;
	while (raAux < i){
		rank_B = raAux;
		curr+=BSrmq;
		rb++;
		if (rb == N8W64)
			rb=0;
		q = (P[curr>>BW64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
		raAux += (BSrmq - T_SUM_BLOCK[q])>>1;
	}
	for (; rank_B<i; curr++){
		if (!readBit64(P, curr))
			rank_B++;
	}

	return curr-1;
}

ulong DFUDSrmq::select_0(ulong i){
	ulong z, nxt, curr, pos, l, r, m;

	if (i <= zeSB){
		l=(i>>1)/bitsSuB;
		r=lenSB;
		m=(lenSB-l)>>1;
		z = (m*bitsSuB-getNum64(TSBlock, m*lgMAX_SupB, lgMAX_SupB))>>1;
		nxt = ((m+1)*bitsSuB-getNum64(TSBlock, (m+1)*lgMAX_SupB, lgMAX_SupB))>>1;
		while ((z>=i || nxt<i) && m<lenSB){
			if(z>=i)
				r = m-1;
			else
				l = m+1;
			m = l+((r-l)>>1);
			z = (m*bitsSuB - getNum64(TSBlock, m*lgMAX_SupB, lgMAX_SupB))>>1;
			nxt = ((m+1)*bitsSuB - getNum64(TSBlock, (m+1)*lgMAX_SupB, lgMAX_SupB))>>1;
		}
		curr = m*bitsSuB + BrmqMOne;
	}else{
		z = zeSB;
		m = (lenSB-1);
		curr = m*bitsSuB + BrmqMOne;
	}

	if(readBit64(Bfull, m)){
		if(TRBlock[m])
			nxt = z;
		else
			nxt = z + Srmq;
	}else
		nxt = z + SrmqM - TRBlock[m];

	if (nxt < i){
		z = nxt;
		curr += Srmq;
	}

	l = 0;
	r = (P[curr>>BW64] & RMMMasks[0]) >> W64m8;
	nxt = z + POPC0[r];
	while (nxt < i){
		z = nxt;
		l++;
		if (l == N8W64)
			l=0;
		curr+=BSrmq;
		r = (P[curr>>BW64] & RMMMasks[l]) >> (W64m8-BSrmq*l);
		nxt += POPC0[r];
	}
	pos=curr-BrmqMOne;
	for (; z<i; pos++){
		if (!readBit64(P, pos))
			z++;
	}

	return pos-1;
}

ulong DFUDSrmq::select_0_old(ulong i){
	ulong pos = 0;

	if(i <= (nBin-rank1_Bin))
		return binSelect_0(i);
	else{
		ulong rb, l, q, curr;
		ulong rank_A, rank_B = nBin-rank1_Bin;

		curr = nBin+BrmqMOne;
		q = (P[curr>>BW64] & RMMMasks[0]) >> W64m8;
		rank_A = rank_B + ((BSrmq - T_SUM_BLOCK[q])>>1);
		rb = l = 0;
		while (rank_A < i){
			rank_B = rank_A;
			rb++;
			if (rb == N8W64)
				rb=0;
			curr+=BSrmq;
			q = (P[curr>>BW64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
			rank_A += (BSrmq - T_SUM_BLOCK[q])>>1;
		}
		pos=curr-BrmqMOne;
		for (; rank_B<i; pos++){
			if (!readBit64(P, pos))
				rank_B++;
		}
	}

	return pos-1;
}

ulong DFUDSrmq::open_0(ulong i){
	ulong pos=0;

	// not necesary!
	return pos;
}

// give the excess from 0 to pos
long int DFUDSrmq::sumAtPos(ulong pos){
	ulong rb, q, l;
	ulong blk = (pos+1)>>PotSrmq;
	ulong sup = blk/SuBrmq;
	long int sum = getNum64(TSBlock, sup*lgMAX_SupB, lgMAX_SupB);

	if(blk%2){
		if(readBit64(Bfull,sup)){
			if(TRBlock[sup])
				sum += Srmq;
			else
				sum -= Srmq;
		}else
			sum += TRBlock[sup] << 1;
	}

	l=blk<<PotSrmq;
	for (rb=0; l+BSrmq<=pos; l+=BSrmq){
		q = (P[l>>BW64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
		sum += T_SUM_BLOCK[q];
		if (rb == N8W64-1) rb=0;
		else rb++;
	}

	for (; l<=pos; l++){
		if (readBit64(P, l))
			sum++;
		else
			sum--;
	}
	return sum;
}

// *******************************************************************************
// ****************************** TREE OPERATIONS ********************************

// give the excess of the internal node 'node=preorder+1' that has a distance 'dist' to the tree's depth
long int DFUDSrmq::computeSumOfNode(ulong node, ulong dist){
	long int sum = 0;
	ulong ini=node<<dist, end=(node+1)<<dist;

	if(ini>cantN){
		ini>>=1;
		end>>=1;
	}
	if (ini >= firstLeaf)
		ini -= firstLeaf+1;
	else
		ini = ini-cantIN+leavesBottom-1;

	if (end >= cantN){
		end>>=1;

		if(end < firstLeaf)
			end = leavesBottom+end-cantIN-1;
		else
			end = leaves;
	}else{
		if (end > firstLeaf)
			end -= firstLeaf+1;
		else
			end = end-cantIN+leavesBottom-1;
	}
	if(ini > end) end = leaves;

	if (ini)
		sum += sumAtPos((end<<PotSrmq)-1) - sumAtPos((ini<<PotSrmq)-1);
	else
		sum += sumAtPos((end<<PotSrmq)-1);
	return sum;
}

// give the number of leaves of this node 'node=preorder+1' that has a distance 'dist' to the tree's depth
ulong DFUDSrmq::computeLeavesOfNode(ulong node, ulong dist){
	ulong ini=node<<dist, end=(node+1)<<dist;

	if(ini>cantN){
		ini>>=1;
		end>>=1;
	}
	if (ini >= firstLeaf)
		ini -= firstLeaf+1;
	else
		ini = ini-cantIN+leavesBottom-1;

	if (end >= cantN){
		end>>=1;

		if(end < firstLeaf)
			end = leavesBottom+end-cantIN-1;
		else
			end = leaves;
	}else{
		if (end > firstLeaf)
			end -= firstLeaf+1;
		else
			end = end-cantIN+leavesBottom-1;
	}
	if(ini > end) end = leaves;

	return end-ini;
}

// give the number of leaves of this node 'node=preorder+1' that has a distance 'dist' to the tree's depth
ulong DFUDSrmq::cantLeavesOfNode(ulong node){
	ulong ini, end;

	for(ini=node; ini<cantIN; ini=2*ini+1);
	if(ini >= firstLeaf)
		ini -= firstLeaf;
	else
		ini = ini - cantIN + leavesBottom;

	for(end=node; end<cantIN; end=2*(end+1));
	if(end >= firstLeaf)
		end -= firstLeaf;
	else
		end = end - cantIN + leavesBottom;

	if(end>leaves)
		end=leaves-1;

	return end-ini+1;
}


// search the minimum in a block (sequential search)
void DFUDSrmq::search_min_block(ulong x1, ulong x2, long int *min, long int *curSum, ulong *position){
	long int Min, sum, auxMin;
	uint rest, q;
	ulong b, len = x2-x1+1, posMin = *position;

	Min = *min;
	sum = *curSum;
	while(len > BrmqMOne){
		b = x1>>BW64;
		rest = (x1+8)%W64;
		if (b == (x1+BrmqMOne)>>BW64)				// x1 and (x1+7) are in the same word...
			q = (P[b] >> (W64-rest)) & 0xff;
		else
			q = ((P[b] << rest) | (P[b+1] >> (W64-rest))) & 0xff;
		auxMin = T_MIN_FWDI[q];
		if (sum-auxMin < Min){
			Min = sum-auxMin;
			posMin = x1+PT_MIN_FWDI[q];
		}
		sum += T_SUM_BLOCK[q];
		len -= BSrmq;
		x1 += BSrmq;
	}
	if (len){
		// check last segment (len < 8) bit by bit...
		while (len>=1){
			if(readBit64(P,x1))
				sum++;
			else{
				sum--;
				if (sum < Min){
					Min = sum;
					posMin = x1;
				}
			}
			x1++;
			len--;
		}
	}
	*curSum = sum;
	*min = Min;
	*position = posMin;
}

// it gives the position of the ")" in where occurs the minimum excess. Forward search (the first minimum)
ulong DFUDSrmq::rmqi(ulong i, ulong j){
	long int min, sum;
	min = sum = 0;

	if(j < nBin)
		return rmqi_rmm(i, j, &min, &sum, i);
	else{
		ulong posMin;

		if(i >= nBin){
			for (sum=0; i<=j && readBit64(P, i); i++, sum++);
			sum--;
			min=sum;
			posMin=i;
			search_min_block(i+1, j, &min, &sum, &posMin);
			return posMin;
		}else{
			if (nBin){
				ulong posMin = rmqi_rmm(i, nBin-1, &min, &sum, i);
				search_min_block(nBin, j, &min, &sum, &posMin);
				return posMin;
			}
		}
	}
	return 0;
}

// return the position of the close parenthesis closet to the root between i and j
ulong DFUDSrmq::rmqi_rmm(ulong x1, ulong x2, long int *min, long int *currSum, ulong posMin){
	ulong rLeaf, aLeaf, node, nodeMin;
	ulong group, distNodeMin, dist = 1;
	long int mini, sumMin;
	bool isInNode = false;
	bool isMinInLeaf = false;

	for (sumMin=0; x1<=x2 && readBit64(P, x1); x1++, sumMin++);
	if(x1>=x2) return x2;
	sumMin--;
	posMin=x1;
	x1++;
	*currSum = sumMin;
	*min=sumMin;
	node = x1>>PotSrmq;
	rLeaf = x2>>PotSrmq;

	// [1]- search the minimum inside the leftmost block...
	if (rLeaf == node){
		// here is the answer
		search_min_block(x1, x2, min, currSum, &posMin);
		return posMin;
	}else
		search_min_block(x1, ((node+1)<<PotSrmq)-1, min, currSum, &posMin);

	long int Min = *min, sum = *currSum;

	// [2]- We search in the next leaf -->
	if (node%2==0){	// this is a left leaf --> we go to the right
		node++;
		if (node == rLeaf){		// here is the answer
			search_min_block((node<<PotSrmq), x2, min, currSum, &posMin);
			return posMin;
		}
		mini = getNum64(TMinB, node*lgMAX_B, lgMAX_B);
		if (sum-mini<Min){
			Min = sum-mini;
			posMin = (node<<PotSrmq) + TPMinB[node];
			isInNode = isMinInLeaf = true;
		}
		aLeaf = node>>1;
		sum += getNum64(TSBlock, (aLeaf+1)*lgMAX_SupB, lgMAX_SupB) - getNum64(TSBlock, aLeaf*lgMAX_SupB, lgMAX_SupB);
		if(readBit64(Bfull, (node>>1))){
			if(TRBlock[node>>1])
				sum += Srmq;
			else
				sum -= Srmq;
		}else
			sum -= TRBlock[aLeaf]<<1;
	}
	if (node == leaves){
		*min = Min;
		return posMin;
	}

	// [3]- We climb recomputing the Min until the position j.
	aLeaf = node;
	if (node < leavesBottom){
		if (node == leavesBottom-1)
			node = cantIN-1;
		else
			node += firstLeaf;
	}else
		node += cantIN-leavesBottom;

	node >>= 1;
	group = cantLeavesOfNode(node);

	while (rLeaf > aLeaf+group){
		aLeaf += group;
		mini = sum-getNum64(Fwd_MinIN, node*lgMIN_Fwd, lgMIN_Fwd);
		if (mini < Min){
			Min = mini;
			isInNode = true;
			isMinInLeaf = false;
			nodeMin = node;
			sumMin = sum;
			distNodeMin = dist;
		}
		sum += computeSumOfNode(node+1, dist);
		if (node%2)
			node++;
		else{
			node>>=1;		// go to my uncle
			dist++;
		}

		group = cantLeavesOfNode(node);
	}

	// [4]- We move down recomputing the Min and reach to the left boundaries leaf.
	node<<= 1;
	node++;
	dist--;
	while (node < cantIN){
		group = cantLeavesOfNode(node);
		if(rLeaf > aLeaf+group){
			aLeaf += group;
			mini = sum-getNum64(Fwd_MinIN, node*lgMIN_Fwd, lgMIN_Fwd);
			if (mini<Min){
				Min = mini;
				isInNode = true;
				isMinInLeaf = false;
				nodeMin = node;
				sumMin = sum;
				distNodeMin = dist;
			}
			sum += computeSumOfNode(node+1, dist);
			node++;
		}else{
			node=2*node+1;
			dist--;
		}
	}


	// at this point, 'node' is a leaf; node >= firstLeaf
	if (node > firstLeaf)			// at this point, 'node' is a leaf; node >= firstLeaf
		node -= firstLeaf;
	else
		node += leavesBottom - cantIN;
	// [5] looking for in the leaves of the rightmost side
	if (node == rLeaf){
		if (isInNode){
			mini = Min;
			search_min_block(node<<PotSrmq, x2, &mini, &sum, &posMin);
			if (mini < Min){
				*min = mini;
				*currSum = sum;
				return posMin;
			}
		}else{
			*min = Min;
			*currSum = sum;
			search_min_block(node<<PotSrmq, x2, min, currSum, &posMin);
			return posMin;
		}
	}else{
		if (isInNode){
			aLeaf = node>>1;
			mini = getNum64(TMinB, node*lgMAX_B, lgMAX_B);
			if (sum - mini < Min){
				Min = sum - mini;
				posMin = (node<<PotSrmq) + TPMinB[node];
				if(readBit64(Bfull, aLeaf)){
					if(TRBlock[aLeaf])
						sum -= Srmq;
					else
						sum += Srmq;
				}else
					sum += (TRBlock[aLeaf]<<1);
				node++;
				*min = Min;
				*currSum = sum;
				search_min_block(node<<PotSrmq, x2, min, currSum, &posMin);
				return posMin;
			}

			if(readBit64(Bfull, aLeaf)){
				if(TRBlock[aLeaf])
					sum -= Srmq;
				else
					sum += Srmq;
			}else
				sum += (TRBlock[aLeaf]<<1);
			mini = Min;
			node++;
			search_min_block(node<<PotSrmq, x2, &mini, &sum, &posMin);
			if (mini < Min){
				*currSum = sum;
				*min = mini;
				return posMin;
			}
		}else{
			aLeaf = node>>1;
			mini = getNum64(TMinB, node*lgMAX_B, lgMAX_B);
			if (sum - mini < Min){
				Min = sum - mini;
				posMin = (node<<PotSrmq) + TPMinB[node];
			}
			if(readBit64(Bfull, aLeaf)){
				if(TRBlock[aLeaf])
					sum -= Srmq;
				else
					sum += Srmq;
			}else
				sum += 2*TRBlock[aLeaf];

			*min = Min;
			*currSum = sum;
			search_min_block((node+1)<<PotSrmq, x2, min, currSum, &posMin);
			return posMin;
		}
	}

	*currSum = sum;
	if(isMinInLeaf){
		*min = Min;
		return posMin;
	}


	// [6] Here the minimum value is in a block which descent of the internal node 'nodeMin',
	// then we must descent in order to find its position.
	nodeMin<<= 1;
	nodeMin++;
	distNodeMin--;
	sum = sumMin;
	while (nodeMin < cantIN){
		mini = getNum64(Fwd_MinIN, nodeMin*lgMIN_Fwd, lgMIN_Fwd);

		if (Min == sum-mini){   // is it ? then go to the right without updating the Min
			nodeMin<<= 1;
			nodeMin++;
			distNodeMin--;
		}else{					// else, go to the right updating the sum
			sumMin = computeSumOfNode(nodeMin+1, distNodeMin);
			nodeMin++;
			sum += sumMin;
		}
	}
	// at this point, 'nodeMin' is a leaf
	if (nodeMin > firstLeaf)
		nodeMin -= firstLeaf;
	else
		nodeMin += leavesBottom - cantIN;

	// Here, we get the left leaf...
	*min = Min;
	aLeaf = nodeMin>>1;
	if(readBit64(Bfull, aLeaf)){
		if(TRBlock[aLeaf-1])
			sumMin = Srmq;
		else
			sumMin = Srmq;
	}else
		sumMin = (TRBlock[aLeaf]<<1);

	if ((long int)getNum64(TMinB, nodeMin*lgMAX_B, lgMAX_B) >= (long int)getNum64(TMinB, (nodeMin+1)*lgMAX_B, lgMAX_B)-sumMin){
		return (nodeMin<<PotSrmq)+TPMinB[nodeMin];
	}else
		return ((nodeMin+1)<<PotSrmq)+TPMinB[nodeMin+1];
}

// exc < 0. return the position with excess = exc from x2 to x1
bool DFUDSrmq::backward_search_block(ulong x1, ulong x2, long int exc, long int *sumPos, ulong *pos){
	long int sum, aux;
	uint rest, q;
	ulong b, len = x2-x1+1;
	sum = *sumPos;

	while(len > BrmqMOne){
		b = x2>>BW64;
		rest = (x2+1)%W64;
		if (b == (x2-BrmqMOne)>>BW64)				// x2 and (x2-7) are in the same word...
			q = (P[b] >> (W64-rest)) & 0xff;
		else
			q = ((P[b-1] << rest) | (P[b] >> (W64-rest))) & 0xff;

		aux = T_MAX_BCKDI[q];
		if (sum-aux <= exc){
			while(sum > exc){
				if(readBit64(P,x2)){
					sum--;
					if(sum==exc){
						*pos = x2;
						*sumPos = sum;
						return true;
					}
				}else
					sum++;
				x2--;
			}
		}
		sum -= T_SUM_BLOCK[q];
		len -= BSrmq;
		x2 -= BSrmq;
	}

	if (len){
		// check last segment (len < 8) bit by bit...
		while (len>=1){
			if(readBit64(P,x2)){
				sum--;
				if(sum==exc){
					*pos = x2;
					*sumPos = sum;
					return true;
				}
			}else
				sum++;
			x2--;
			len--;
		}
	}

	*sumPos = sum;
	return false;
}

ulong DFUDSrmq::backward_search_rmm(long int exc, ulong j){
	ulong node, posExc, dist = 1;
	long int currSum, min, sNod;

	currSum = 0;
	node = j>>PotSrmq;

	// [1]- search the excess inside the rightmost block...
	if (!node){
		// here is the answer
		if (backward_search_block(0, j, exc, &currSum, &posExc))
			return posExc;
	}else{
		if (backward_search_block(node<<PotSrmq, j, exc, &currSum, &posExc))
			return posExc;
	}


	// [2]- We search in the next leaf <--
	if (node%2){	// this is a right leaf --> we go to the left
		if (backward_search_block((node-1)<<PotSrmq, (node<<PotSrmq)-1, exc, &currSum, &posExc))
			return posExc;
		node--;
	}

	// [3]- We climb recomputing the Min until the position i.
	if (node < leavesBottom){
		if (node == leavesBottom-1)
			node = cantIN-1;
		else
			node += firstLeaf;
	}else
		node += cantIN-leavesBottom;
	node>>= 1;
	node--;
	dist++;

	sNod = computeSumOfNode(node+1, dist);
	min = -1*sNod-getNum64(Fwd_MinIN, node*lgMIN_Fwd, lgMIN_Fwd);
	while (min > exc-currSum){
		currSum -= sNod;
		if (node%2==0)			// go to my left sibling
			node--;
		else{
			node=node/2-1;		// go to my uncle
			dist++;
		}
		sNod = computeSumOfNode(node+1, dist);
		min = -1*sNod-getNum64(Fwd_MinIN, node*lgMIN_Fwd, lgMIN_Fwd);
	}

	// [4]- We move down recomputing the sum until to reach to the left that store the excess
	node++;
	node<<= 1;
	dist--;
	while (node < cantIN){
		sNod = computeSumOfNode(node+1, dist);
		min = -1*sNod-getNum64(Fwd_MinIN, node*lgMIN_Fwd, lgMIN_Fwd);
		if (min <= exc-currSum){
			node++;			// go to my nephew
			node<<= 1;
			dist--;
		}else{
			currSum -= sNod;
			node--;			// go to my left sibling
		}
	}
	// at this point, 'node' is a the leaf that contains the excess searched
	if (node > firstLeaf)
		node -= firstLeaf;
	else
		node += leavesBottom - cantIN;

	// this covers two leaves:node-1 and node
	backward_search_block((node-1)<<PotSrmq, ((node+1)<<PotSrmq)-1, exc, &currSum, &posExc);
	return posExc;

}

ulong DFUDSrmq::backward_search(long int exc, ulong j){
	if(j < nBin)
		return backward_search_rmm(exc, j);
	else{
		// 1. si esta excess en el ultimo bloque--> retornamos su posicion
		ulong posExc;
		long int sum = 0;
		if (backward_search_block(nBin, j, exc, &sum, &posExc))
			return posExc;

		if (nBin)
			// 2. sino esta--> buscamos en el arbol
			return backward_search_rmm(exc-sum, nBin-1);
	}
	return j;
}

ulong DFUDSrmq::queryRMQ(ulong i, ulong j){
	if (i>=j)
		return i;

	ulong x = select_0(i+2);
	ulong y = select_0(j+1);
	ulong w = rmqi(x,y);
	ulong openW = backward_search(-1,w-1);

	if(openW-rank_1(openW) == i)
		return i;

	return w-rank_1(w);
}

uint DFUDSrmq::getSize(){
	return sizeDS;
}

void DFUDSrmq::saveDS(char *fileName){
	cout << "Save data structure in " << fileName << endl;
	ofstream os (fileName, ios::binary);
	cout << "   Data structure size: " << sizeDS << endl;

	if(TRACE){
		cout << "Variables load: " << endl;
		cout << "nP " << nP << endl;
		cout << "nW " << nW << endl;
		cout << "rank1_Bin " << rank1_Bin << endl;
		cout << "nBin " << nBin << endl;
		cout << "cantN " << cantN << endl;
		cout << "cantIN " << cantIN << endl;
		cout << "leaves " << leaves << endl;
		cout << "leavesBottom " << leavesBottom << endl;
		cout << "firstLeaf " << firstLeaf << endl;
		cout << "lenSB " << lenSB << endl;
		cout << "zeSB " << zeSB << endl;
		cout << "lenLB " << lenLB << endl;
		cout << "bitsSuB " << bitsSuB << endl;
		cout << "bitsRB " << bitsRB << endl;
		cout << "h " << h << endl;
		cout << "MAX_B " << MAX_B << endl;
		cout << "lgMAX_B " << lgMAX_B << endl;
		cout << "MAX_SupB " << MAX_SupB << endl;
		cout << "lgMAX_SupB " << lgMAX_SupB << endl;
		cout << "lgMIN_Fwd " << lgMIN_Fwd << endl;
		cout << "MIN_Fwd " << MIN_Fwd << endl;
	}

	os.write((const char*)&nP, sizeof(ulong));
	os.write((const char*)&nW, sizeof(ulong));
	os.write((const char*)&rank1_Bin, sizeof(ulong));
	os.write((const char*)&nBin, sizeof(ulong));
	os.write((const char*)&cantN, sizeof(ulong));
	os.write((const char*)&cantIN, sizeof(ulong));
	os.write((const char*)&leaves, sizeof(ulong));
	os.write((const char*)&leavesBottom, sizeof(ulong));
	os.write((const char*)&firstLeaf, sizeof(ulong));
	os.write((const char*)&lenSB, sizeof(ulong));
	os.write((const char*)&zeSB, sizeof(ulong));
	os.write((const char*)&lenLB, sizeof(uint));
	os.write((const char*)&bitsSuB, sizeof(uint));
	os.write((const char*)&bitsRB, sizeof(uint));
	os.write((const char*)&h, sizeof(uint));
	os.write((const char*)&MAX_B, sizeof(uint));
	os.write((const char*)&lgMAX_B, sizeof(uint));
	os.write((const char*)&MAX_SupB, sizeof(uint));
	os.write((const char*)&lgMAX_SupB, sizeof(uint));
	os.write((const char*)&lgMIN_Fwd, sizeof(uint));
	os.write((const char*)&MIN_Fwd, sizeof(int));

	TRACE = true;

	ulong sizeDT = 11*sizeof(ulong) + 10*sizeof(uint) + sizeof(int);
	sizeDT +=  2*(512 + 256);				// size for T_SUM_BLOCK[] + T_MIN_FWDI[] + T_MAX_BCKDI[] + PT_MIN_FWDI[]
	if(TRACE) cout << " .- T_SUM_BLOCK[] + T_MIN_FWDI[] + T_MAX_BCKDI[] + PT_MIN_FWDI[] + Variables " << sizeDT << " Bytes" << endl;

	ulong size = nP >> BW64;
	if (nP % W64)
		size++;
	os.write((const char*)P, size*sizeof(ulong));				// save P[]
	sizeDT += size*sizeof(ulong);
	if(TRACE) cout << " .- P[] " << size*sizeof(ulong) << " Bytes" << endl;

	size = cantIN*lgMIN_Fwd/W64;
	if ((cantIN*lgMIN_Fwd)%W64)
		size++;
	os.write((const char*)Fwd_MinIN, size*sizeof(ulong));		// save Fwd_MinIN[]
	sizeDT += size*sizeof(ulong);
	if(TRACE) cout << " .- Fwd_MinIN[] " << size*sizeof(ulong) << " Bytes" << endl;

	size = lenSB*lgMAX_SupB/W64;
	if ((lenSB*lgMAX_SupB)%W64)
		size++;
	os.write((const char*)TSBlock, size*sizeof(ulong));			// save TSBlock[]
	sizeDT += size*sizeof(ulong);
	if(TRACE) cout << " .- TSBlock[] " << size*sizeof(ulong) << " Bytes" << endl;

	os.write((const char*)TRBlock, lenSB*sizeof(char));			// save TRBlock[]
	sizeDT += lenSB*sizeof(char);
	if(TRACE) cout << " .- TRBlock[] " << lenSB*sizeof(char) << " Bytes" << endl;

	size = lenSB/W64;
	if (lenSB%W64)
		size++;
	os.write((const char*)Bfull, size*sizeof(ulong));			// save Bfull[]
	sizeDT += size*sizeof(ulong);
	if(TRACE) cout << " .- Bfull[] " << size*sizeof(ulong) << " Bytes" << endl;

	os.write((const char*)TPMinB, leaves*sizeof(uchar));		// save TPMinB[]  = new uchar[leaves]
	sizeDT += leaves*sizeof(uchar);
	if(TRACE) cout << " .- TPMinB[] " << leaves*sizeof(uchar) << " Bytes" << endl;

	size = leaves*lgMAX_B/W64;
	if ((leaves*lgMAX_B)%W64)
		size++;
	os.write((const char*)TMinB, size*sizeof(ulong));							// save TMinB[]
	sizeDT += size*sizeof(ulong);
	if(TRACE) cout << " .- TMinB[] " << size*sizeof(ulong) << " Bytes" << endl;

	os.close();
	cout << "   Total bytes saved from data structure: " << sizeDT << endl;
	TRACE = false;
}

void DFUDSrmq::loadDS(char *fileName){
	cout << " Load data structure from " << fileName << endl;
	ifstream is(fileName, ios::binary);

	is.read((char*)&nP, sizeof(ulong));
	is.read((char*)&nW, sizeof(ulong));
	is.read((char*)&rank1_Bin, sizeof(ulong));
	is.read((char*)&nBin, sizeof(ulong));
	is.read((char*)&cantN, sizeof(ulong));
	is.read((char*)&cantIN, sizeof(ulong));
	is.read((char*)&leaves, sizeof(ulong));
	is.read((char*)&leavesBottom, sizeof(ulong));
	is.read((char*)&firstLeaf, sizeof(ulong));
	is.read((char*)&lenSB, sizeof(ulong));
	is.read((char*)&zeSB, sizeof(ulong));
	is.read((char*)&lenLB, sizeof(uint));
	is.read((char*)&bitsSuB, sizeof(uint));
	is.read((char*)&bitsRB, sizeof(uint));
	is.read((char*)&h, sizeof(uint));
	is.read((char*)&MAX_B, sizeof(uint));
	is.read((char*)&lgMAX_B, sizeof(uint));
	is.read((char*)&MAX_SupB, sizeof(uint));
	is.read((char*)&lgMAX_SupB, sizeof(uint));
	is.read((char*)&lgMIN_Fwd, sizeof(uint));
	is.read((char*)&MIN_Fwd, sizeof(int));

	if(TRACE){
		cout << "Variables load: " << endl;
		cout << "nP " << nP << endl;
		cout << "nW " << nW << endl;
		cout << "rank1_Bin " << rank1_Bin << endl;
		cout << "nBin " << nBin << endl;
		cout << "cantN " << cantN << endl;
		cout << "cantIN " << cantIN << endl;
		cout << "leavesBottom " << leavesBottom << endl;
		cout << "leaves " << leaves << endl;
		cout << "firstLeaf " << firstLeaf << endl;
		cout << "lenSB " << lenSB << endl;
		cout << "zeSB " << zeSB << endl;
		cout << "lenLB " << lenLB << endl;
		cout << "bitsSuB " << bitsSuB << endl;
		cout << "bitsRB " << bitsRB << endl;
		cout << "h " << h << endl;
		cout << "MAX_B " << MAX_B << endl;
		cout << "lgMAX_B " << lgMAX_B << endl;
		cout << "MAX_SupB " << MAX_SupB << endl;
		cout << "lgMAX_SupB " << lgMAX_SupB << endl;
		cout << "lgMIN_Fwd " << lgMIN_Fwd << endl;
		cout << "MIN_Fwd " << MIN_Fwd << endl;
	}

	TRACE = true;
	// size for variables
	sizeDS = 11*sizeof(ulong) + 10*sizeof(uint) + sizeof(int);
	sizeDS += 2*(512 + 256);									// size for T_SUM_BLOCK[] + T_MIN_FWDI[] + T_MAX_BCKDI[] + PT_MIN_FWDI[]
	if(TRACE) cout << " .- T_SUM_BLOCK[] + T_MIN_FWDI[] + T_MAX_BCKDI[] + PT_MIN_FWDI[] + Variables " << sizeDS << " Bytes" << endl;

	ulong sizeAux = nP >> BW64;
	if (nP % W64)
		sizeAux++;
	P = new ulong[sizeAux];
	is.read((char*)P, sizeAux*sizeof(ulong));
	sizeDS += sizeAux*sizeof(ulong);
	if(TRACE) cout << " .- P[] " << sizeAux*sizeof(ulong) << " Bytes" << endl;

	sizeAux = cantIN*lgMIN_Fwd/W64;
	if ((cantIN*lgMIN_Fwd)%W64)
		sizeAux++;
	Fwd_MinIN = new ulong[sizeAux];
	is.read((char*)Fwd_MinIN, sizeAux*sizeof(ulong));
	sizeDS += sizeAux*sizeof(ulong);
	if(TRACE) cout << " .- Fwd_MinIN[] " << sizeAux*sizeof(ulong) << " Bytes" << endl;

	sizeAux = lenSB*lgMAX_SupB/W64;
	if ((lenSB*lgMAX_SupB)%W64)
		sizeAux++;
	TSBlock = new ulong[sizeAux];
	is.read((char*)TSBlock, sizeAux*sizeof(ulong));
	sizeDS += sizeAux*sizeof(ulong);
	if(TRACE) cout << " .- TSBlock[] " << sizeAux*sizeof(ulong) << " Bytes" << endl;

	TRBlock = new char[lenSB];
	is.read((char*)TRBlock, lenSB*sizeof(char));
	sizeDS += lenSB*sizeof(char);
	if(TRACE) cout << " .- TRBlock[] " << lenSB*sizeof(char) << " Bytes" << endl;

	sizeAux = lenSB/W64;
	if (lenSB%W64)
		sizeAux++;
	Bfull = new ulong[sizeAux];
	is.read((char*)Bfull, sizeAux*sizeof(ulong));
	sizeDS += sizeAux*sizeof(ulong);
	if(TRACE) cout << " .- Bfull[] " << sizeAux*sizeof(ulong) << " Bytes" << endl;

	TPMinB = new uchar[leaves];
	is.read((char*)TPMinB, leaves*sizeof(uchar));
	sizeDS += leaves*sizeof(uchar);
	if(TRACE) cout << " .- TPMinB[] " << leaves*sizeof(uchar) << " Bytes" << endl;

	sizeAux = leaves*lgMAX_B/W64;
	if ((leaves*lgMAX_B)%W64)
		sizeAux++;
	TMinB = new ulong[sizeAux];
	is.read((char*)TMinB, sizeAux*sizeof(ulong));
	sizeDS += sizeAux*sizeof(ulong);
	if(TRACE) cout << " .- TMinB[] " << sizeAux*sizeof(sizeDS) << " Bytes" << endl;

	is.close();
	cout << " Data Structure loaded !!" << endl;
	TRACE = false;
}

// *************************** PRINT THE TREE ************************************
void DFUDSrmq::printTree(){
	ulong i, j, cant, acum, rb, segment;
	long int currSum, mini;
	uint cMIN_BCK;
	cMIN_BCK = 0;

	cout << "Min-max Binary Tree...[min bck]i(sum)" << endl;
	cant = 1;
	for (acum=j=i=0; i<cantIN; i++, j++, acum++){
		if (j==cant){
			cout << endl;
			j = acum = 0;
			cant *= 2;
		}
		if (acum == 2){
			cout << " | ";
			acum = 0;
		}

		mini = getNum64(Fwd_MinIN, cMIN_BCK, lgMIN_Fwd);
		cout << i << "_[" << mini << "]";
		cMIN_BCK += lgMIN_Fwd;
	}

	for (i=leavesBottom; i<leaves; i++, j++, acum++){
		if (j==cant){
			cout << endl;
			j = acum = 0;
			cant *= 2;
		}
		if (acum == 2){
			cout << " | ";
			acum = 0;
		}

		rb = 3;
		currSum = 0;
		mini = 0;
		for (j=0; j<N8Srmq; j++){
			segment = (P[((i+1)*Srmq-1-BSrmq*j)>>BW64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
			if (currSum + T_MIN_FWDI[segment] > mini)
				mini = currSum + T_MIN_FWDI[segment];
			currSum += T_SUM_BLOCK[segment];
			if (rb == 0) rb=N8W64-1;
			else rb--;
		}
		cout << i << "_h[" << mini << "](" << currSum << ") ";
	}
	cout << endl;
	for (i=0; i<leavesBottom; i++, j++, acum++){
		if (j==cant){
			cout << endl;
			j = acum = 0;
			cant *= 2;
		}
		if (acum == 2){
			cout << " | ";
			acum = 0;
		}
		rb = 3;
		currSum = 0;
		mini = 0;
		for (j=0; j<N8Srmq; j++){
			segment = (P[((i+1)*Srmq-1-BSrmq*j)>>BW64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
			if (currSum + T_MIN_FWDI[segment] > mini)
				mini = currSum + T_MIN_FWDI[segment];
			currSum += T_SUM_BLOCK[segment];
			if (rb == 0) rb=N8W64-1;
			else rb--;
		}
		cout << i << "_h[" << mini << "]"<< i << "(" << currSum << ") ";
	}
	cout << endl;
}

DFUDSrmq::~DFUDSrmq() {
	nP = nW = rank1_Bin = nBin = cantN = cantIN = leaves =
	leavesBottom = firstLeaf = lenSB = lenLB = bitsSuB =
	bitsRB = h = MAX_B = lgMAX_B = MAX_SupB = lgMAX_SupB =
	lgMIN_Fwd = MIN_Fwd = sizeDS = 0;

	delete [] Fwd_MinIN;
	//delete [] Bck_MaxIN;
	delete [] TSBlock;
	delete [] TRBlock;
	delete [] Bfull;
	delete [] TPMinB;
	delete [] TMinB;
	cout << " ~ DFUDSrmq destroyed !!" << endl;
}

// *******************************************************************************
// **************************** TEST TREE OPERATIONS *****************************
void DFUDSrmq::test_rank_1(){
	long int sumR, rank;
	ulong i,j,k;
	cout << "DFUDSrmq::test_rank_1..." << endl;

	/*i=224;
	for (j=sumR=0; j<=i; j++){
		if(readBit64(P, j))
			sumR++;
	}
	rank = rank_1(i);
	cout << "rank_1(" << i <<") = " << rank << ", ? sumR =" << sumR << endl;
	exit(0);*/

	if (nBin>0){
		i = nBin-1;
		for (j=sumR=0; j<=i; j++){
			if(readBit64(P, j))
				sumR++;
		}
		rank = rank_1(i);
		if (sumR != rank){
			cout << "sumR("<<i<<") = " << sumR << endl;
			cout << "rank_1("<<i<<") = " << rank << endl;
			exit(0);
		}
	}

	for (k=0; k<TEST; k++){
		i = (rand() % (nP-1));
		for (j=sumR=0; j<=i; j++){
			if(readBit64(P, j))
				sumR++;
		}
		rank = rank_1(i);
		if (sumR != rank){
			cout << "ERROR !! rank1(" << i << ") = " << rank << " != sumR = " << sumR << endl;
			exit(1);
		}
	}
	cout << "  test_rank_1 OK !!" << endl;
}

void DFUDSrmq::test_select_1(){
	ulong i, j, k, sum, pos;

	cout << "DFUDSrmq::test_select_1..." << endl;
	/*i=223;
	for (j=sum=0; sum<i; j++){
		if(readBit64(P, j))
			sum++;
	}
	j--;
	pos = select_1(i);
	cout << "select_1(" << i <<") = " << pos << ", ? j=" << j << endl;
	exit(0);*/

	for (k=0; k<TEST; k++){
		i = (rand() % ((nP>>1)-2)) + 1;
		for (j=sum=0; sum<i; j++){
			if(readBit64(P, j))
				sum++;
		}
		j--;
		pos = select_1(i);
		if (j != pos){
			cout << "ERROR !! select1(" << i << ") = " << pos << " != j = " << j << endl;
			exit(1);
		}
	}
	cout << "  test_select_1 OK !!" << endl;
}

void DFUDSrmq::test_select_0(){
	ulong sum0,sel,i,j,k;

	cout << "DFUDSrmq::test_select_0..." << endl;
	/*i=4995;
	for (j=sum0=0; sum0<i; j++){
		if(!readBit64(P, j))
			sum0++;
	}
	j--;
	cout << "brute select_0("<<i<<") = " << j << endl;
	sel = select_0(i);
	cout << "funct. select_0 = " << sel << endl;*/
	//exit(0);

	for (k=0; k<TEST; k++){
		i = (rand() % ((nP>>1)-2)) + 1;
		for (j=sum0=0; sum0<i; j++){
			if(!readBit64(P, j))
				sum0++;
		}
		j--;
		sel = select_0(i);
		if (j != sel){
			cout << "ERROR !! sel0(" << i << ") = " << sel << " != select_brute = " << j << endl;
			exit(1);
		}
	}
	cout << "  test_select_0 OK !!" << endl;
}

void DFUDSrmq::test_sumAtPos(){
	ulong i, j, k;
	long int sum, sumPos;

	cout << "DFUDSrmq::test_sumAtPos..." << endl;
	/*i=59;
	for (j=sum=0; j<=i; j++){
		if(readBit64(P, j))
			sum++;
		else
			sum--;
	}
	sumPos = sumAtPos(i);
	cout << "sumAtPos(" << i << ") = " << sumPos << " ?= sum = " << sum << endl;
	if (sum != sumPos){
		cout << "ERROR !! sumAtPos(" << i << ") = " << sumPos << " != sum = " << sum << endl;
		exit(1);
	}*/

	for (k=0; k<TEST; k++){
		i = (rand() % (nP-2));
		for (j=sum=0; j<=i; j++){
			if(readBit64(P, j))
				sum++;
			else
				sum--;
		}
		sumPos = sumAtPos(i);
		if (sum != sumPos){
			cout << "ERROR !! sumAtPos(" << i << ") = " << sumPos << " != sum = " << sum << endl;
			exit(1);
		}
	}
	cout << "  test_sumAtPos OK !!" << endl;
}

void DFUDSrmq::test_rmqi(){
	long long i, j, posMin, posRmqi, sum, Min, k;

	/*i=41846, j=50411;
	for (k=i, sum=0; k<=j && readBit64(P, k); k++, sum++);
	sum--;
	Min=sum;
	posMin=k;
	for (k++; k<=j; k++){
		if(readBit64(P, k))
			sum++;
		else{
			sum--;
			if (sum < Min){
				Min = sum;
				posMin = k;
			}
		}
	}
	cout << "MINIMO(" << i << ", " << j << ") = " << Min << "... posMin = " << posMin << endl;
	posRmqi = rmqi(i,j);
	if (posRmqi != posMin){
		cout << "ERROR !!... rmqi(" << i << ", " << j << ") = " << posRmqi << " != " << posMin << endl;
		exit(1);
	}*/

	cout << "DFUDSrmq::test_rmqi..." << endl;
	for (ulong t=0; t<TEST; t++){
		i = (rand() % (nP/2));
		j = nP/2 + (rand() % (nP/2)-2);

		for (k=i, sum=0; k<=j && readBit64(P, k); k++, sum++);
		sum--;
		Min=sum;
		posMin=k;
		if (k>j) continue;
		for (k++; k<=j; k++){
			if(readBit64(P, k))
				sum++;
			else{
				sum--;
				if (sum < Min){
					Min = sum;
					posMin = k;
				}
			}
		}

		posRmqi = rmqi(i,j);
		if (posRmqi != posMin){
			cout << "ERROR !!... rmqi(" << i << ", " << j << ") = " << posRmqi << " != " << posMin << endl;
			exit(1);
		}
	}

	cout << "  test_rmqi OK !!" << endl;
}

void DFUDSrmq::test_search_min_block(){
	long int sum, Min, k, i, j, posMin;

	/*i=24, j=65;
	for (k=i, sum=0; k<=(long int)j && readBit64(P, k); k++, sum++);
	sum--;
	Min=sum;
	posMin=k;
	for (k++; k<=(long int)j; k++){
		if(readBit64(P, k))
			sum++;
		else{
			sum--;
			if (sum < Min){
				Min = sum;
				posMin = k;
			}
		}
	}
	cout << "MINIMO(" << i << ", " << j << ") = " << Min << ", posMin = " << posMin << endl;
	long int min, curSum;
	ulong position;
	min = curSum = 0;
	position = i;
	for (k=i, curSum=0; k<=(long int)j && readBit64(P, k); k++, curSum++);
	curSum--;
	min=curSum;
	position=k;
	search_min_block(k+1,j,&min, &curSum, &position);
	if (position != posMin || min != Min){
		cout << "ERROR !!... search_min_block(" << i << ", " << j << ") = " << position << " != " << posMin << endl;
		cout << "minimum found = " << min << " =? " << Min << endl;
		exit(1);
	}*/


	cout << "DFUDSrmq::test_search_min_block()..." << endl;

	for (ulong t=0; t<TEST; t++){
		i = (rand() % (nP/2));
		j = nP/2 + (rand() % (nP/2)-2);

		for (k=i, sum=0; k<=(long int)j && readBit64(P, k); k++, sum++);
		sum--;
		Min=sum;
		posMin=k;
		if (k>j) continue;
		for (k++; k<=(long int)j; k++){
			if(readBit64(P, k))
				sum++;
			else{
				sum--;
				if (sum < Min){
					Min = sum;
					posMin = k;
				}
			}
		}

		//if (TRACE) cout << "search_min_block(" << i << ", " << j << ") " << endl;
		long int min, curSum;
		ulong position;
		min = curSum = 0;
		position = i;

		for (k=i, curSum=0; k<=(long int)j && readBit64(P, k); k++, curSum++);
		curSum--;
		min=curSum;
		position=k;
		search_min_block(k+1,j,&min, &curSum, &position);
		if (position != posMin || min != Min){
			cout << "ERROR !!... search_min_block(" << i << ", " << j << ") = " << position << " != " << posMin << endl;
			cout << "minimum found = " << min << " =? " << Min << endl;
			exit(1);
		}
	}

	cout << "  test_search_min_block() OK !!" << endl;
}




void DFUDSrmq::test_backward_search_block(){
	long int sum, exc, k;
	ulong j, posExc;

	/*j=62;
	exc = 1;
	sum = 0;
	k = j;
	while(k && sum != exc){
		k--;
		if(readBit64(P, k)){
			sum++;
		}else
			sum--;
	}
	cout << "OPEN PAR (" << j << ") = " << k << endl;

	sum = 0;
	if(backward_search_block(0, j-1, exc, &sum, &posExc)){
		if (k != posExc){
			cout << "ERROR !!... backward_search_block(0, " << j << ", " << exc << ") = " << posExc << " != " << k << endl;
			exit(1);
		}
	}*/

	cout << "DFUDSrmq::test_backward_search_block()..." << endl;
	for (ulong t=0; t<TEST; t++){
		j = nP/2 + (rand() % (nP/2)-2);
		for (; readBit64(P, j); j++);
		exc = -1;
		sum = 0;
		k = j;
		while(k && sum != exc){
			k--;
			if(readBit64(P, k)){
				sum--;
			}else
				sum++;
		}
		// k is the position !!

		sum = 0;
		if(backward_search_block(0, j-1, exc, &sum, &posExc)){
			if (k != posExc){
				cout << "ERROR !!... backward_search_block(0, " << j << ", " << exc << ") = " << posExc << " != " << k << endl;
				exit(1);
			}
		}
	}
	cout << "  test_backward_search_block() OK !!" << endl;
}

void DFUDSrmq::test_backward_search(){
	long int sum, exc, k;
	ulong j, posExc;

	/*j=558;
	if (readBit64(P, j))
		cout << "OPEN" << endl;
	for (; readBit64(P, j); j++);
	exc = -1;
	sum = 0;
	k = j;
	while(k && sum != exc){
		k--;
		if(readBit64(P, k)){
			sum--;
		}else
			sum++;
	}
	cout << "OPEN PAR (" << j << ") = " << k << endl;

	sum = 0;
	posExc = backward_search(exc, j-1);
	if (k != posExc){
		cout << "ERROR !!... backward_search(0, " << j << ", " << exc << ") = " << posExc << " != " << k << endl;
		exit(1);
	}*/

	cout << "DFUDSrmq::test_backward_search()..." << endl;
	for (ulong t=0; t<TEST; t++){
		j = nP/2 + (rand() % (nP/2)-2);

		for (; readBit64(P, j); j++);
		exc = -1;
		sum = 0;
		k = j;
		while(k && sum != exc){
			k--;
			if(readBit64(P, k)){
				sum--;
			}else
				sum++;
		}
		// k is the position !!

		sum = 0;
		posExc = backward_search(exc, j-1);
		if (k != posExc){
			cout << "ERROR !!... backward_search(0, " << j << ", " << exc << ") = " << posExc << " != " << k << endl;
			exit(1);
		}
	}
	cout << "  test_backward_search() OK !!" << endl;
}
