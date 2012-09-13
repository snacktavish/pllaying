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
	double len0; /* original length of the NNI branch */
	double len1; /* optimized length of the NNI branch before the NNI move */
	double len2; /* optimized length of the NNI branch after the NNI move */
} nniMove;

typedef struct {
	tree* tree;
	nniMove** data;
} nniList;

/*
 *  Find the best NNI move for the current branch
 *  Return NULL if no positive NNI is found
 *  Otherwise return the best positive NNI move found
 */
nniMove* evalNNIForBranch(tree* tr, nodeptr p);

double doNNI(tree * tr, nodeptr p, int swap, int optBran);

/*
 *  Go through all 2(n-3) internal branches of the tree and
 *  evalute all possible NNI moves
 */
void evalAllNNI(tree* tr);

void evalNNIForSubtree(tree* tr, nodeptr p, nniMove** nniList, int* cnt);
/*
 *  Save the likelihood vector of p and q to the 2 pointer p_lhsave and
 *  q_lhsave.
 *  Should I use memcpy or just copy the pointer ?
 */
//void saveLHVector(nodeptr p, nodeptr q, double* p_lhsave, double* q_lhsave);
