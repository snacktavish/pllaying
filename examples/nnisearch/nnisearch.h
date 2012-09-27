void read_msa (tree * tr, char * filename);

void makeParsimonyTree(tree *tr)
{
  allocateParsimonyDataStructures(tr);
  makeParsimonyTreeFast(tr);
  freeParsimonyDataStructures(tr);
}

typedef struct {
	tree* tr;
	nodeptr p;
	int nniType;
	double z[NUM_BRANCHES];
	double likelihood;
	double deltaLH;
} nniMove;

int cmp_nni(const void* nni1, const void* nni2) {
	nniMove* myNNI1 = (nniMove*) nni1;
	nniMove* myNNI2 = (nniMove*) nni2;
	return (int) (1000000.f*myNNI1->deltaLH - 1000000.f*myNNI2->deltaLH);
}
/*
 *  Find the best NNI move for the current branch
 *  Return NULL if no positive NNI is found
 *  Otherwise return the best positive NNI move found
 */
nniMove getBestNNIForBran(tree* tr, nodeptr p);

double doOneNNI(tree * tr, nodeptr p, int swap, int optBran);

/*
 *  Go through all 2(n-3) internal branches of the tree and
 *  evaluate all possible NNI moves
 */
void evalAllNNI(tree* tr);

/*
 *  do a full round of fast NNI
 */
double doNNISearch(tree* tr);

/*
 *  cnt: number of internal branches that have been visited
 *  cnt_nni: number of positive NNI found
 */
void evalNNIForSubtree(tree* tr, nodeptr p, nniMove* nniList, int* cnt_bran, int* cnt_nni);
/*
 *  Save the likelihood vector of p and q to the 2 pointer p_lhsave and
 *  q_lhsave.
 *  Should I use memcpy or just copy the pointer ?
 */
//void saveLHVector(nodeptr p, nodeptr q, double* p_lhsave, double* q_lhsave);
