#include <stdio.h>
#include <stdlib.h>
#define GLOBAL_VARIABLES_DEFINITION
#include "axml.h"
#include "globalVariables.h"
#include "nnisearch.h"


/*
 *  Traverse the tree to find a positive NNI.
 *  If such an NNI is found, apply it right away and return.
 */
void doNNISimple(tree* tr) {
	nodeptr p,q;
	p = tr->start;
	if ( !isTip(p->number, tr->mxtips) ) {
		q = p->next;
		while (q != p) {
			nniMove* posNNI = evalNNIForBranch(tr, q);
			q = q->next;
		}
	}
}

double doNNI(tree * tr, nodeptr p, int swap, int optBran) {
	nodeptr q;
	nodeptr tmp;

	q = p->back;
	assert(!isTip(q->number, tr->mxtips));
	assert(!isTip(p->number, tr->mxtips));

	if (swap == 1) {
		tmp = p->next->back;
		hookup(p->next, q->next->back, q->next->z, tr->numBranches);
		hookup(q->next, tmp, tmp->z, tr->numBranches);
	} else {
		tmp = p->next->next->back;
		hookup(p->next->next, q->next->back, q->next->z, tr->numBranches);
		hookup(q->next, tmp, tmp->z, tr->numBranches);
	}
	if (optBran) {
		newviewGeneric(tr, p, FALSE);
		newviewGeneric(tr, q, FALSE);
		update(tr, p);
		//printf("New branch length %f \n", getBranchLength(tr, 0, p) );
		evaluateGeneric(tr,p, FALSE);
		return tr->likelihood;
	} else {
		newviewGeneric(tr, p, FALSE);
		newviewGeneric(tr, q, FALSE);
		evaluateGeneric(tr, p, FALSE);
		return tr->likelihood;
	}
}


saveLHVector(p,q, p_lhsave, q_lhsave) {

}


nniMove* evalNNIForBranch(tree* tr, nodeptr p) {

	if ( isTip(p->number, tr->mxtips) ) {
		return NULL;
	}

	Tree2String(tr->tree_string, tr, tr->start->back, TRUE, FALSE, 0, 0, 0, SUMMARIZE_LH, 0,0);
	fprintf(stderr, "%s\n", tr->tree_string);
	/* Backup the current branch length */
	double zOld = p->z[0];
	//printf("zOld: %f \n", getBranchLength(tr, 0, p));
	double lhOld = tr->likelihood;
	printf("lhOld: %f \n", lhOld);

	/* Optimize the NNI branch */
	update(tr, p);

	double zNew = p->z[0];
	evaluateGeneric(tr,tr->start, TRUE);
	double lhNew = tr->likelihood;
	//printf("zNew: %f \n", getBranchLength(tr, 0, p));
	printf("lhNew: %f \n", lhNew);

	/* Save the likelihood vector at node p and q */
	//saveLHVector(p, q, p_lhsave, q_lhsave);
	/* Save the scaling factor */

	/******************** NNI part *******************************/
	/* Now try to do an NNI move of type 1 */
	double lh_nni1 = doNNI(tr, p, 1, TRUE);
	printf("likelihood of the 1.NNI move: %f\n", lh_nni1);
	//printTopology(tr, TRUE);

	/* Restore previous NNI move */
	doNNI(tr, p, 1, FALSE);
	/* Restore the old branch length */
	p->z[0] = zOld;
	p->back->z[0] = zOld;
	printf("Restore topology\n");
	Tree2String(tr->tree_string, tr, tr->start->back, TRUE, FALSE, 0, 0, 0, SUMMARIZE_LH, 0,0);
	fprintf(stderr, "%s\n", tr->tree_string);
	evaluateGeneric(tr, tr->start, TRUE);
	printf("Likelihood after restoring from NNI 1: %f\n", tr->likelihood);

	/* Try to do an NNI move of type 2 */
	double lh_nni2 = doNNI(tr, p, 0, TRUE);
	printf("likelihood of the 2.NNI move: %f\n", lh_nni2);
	//printTopology(tr, TRUE);

	/* Restore previous NNI move */
	doNNI(tr, p, 0, FALSE);
	/* Restore the old branch length */
	p->z[0] = zOld;
	p->back->z[0] = zOld;
	printf("Restore topology\n");
	Tree2String(tr->tree_string, tr, tr->start->back, TRUE, FALSE, 0, 0, 0, SUMMARIZE_LH, 0,0);
	fprintf(stderr, "%s\n", tr->tree_string);

	evaluateGeneric(tr, tr->start, TRUE);
	printf("Likelihood after restoring from NNI 2: %f\n", tr->likelihood);

	int bestType = -1;
	double best = lhNew;

	if (lh_nni1 > best || lh_nni2 > best) {
		if (lh_nni1 > lh_nni2) {
			best = lh_nni1;
			bestType = 0;
		} else {
			best = lh_nni2;
			bestType = 1;
		}
	} else {
		return NULL;
	}

	nniMove *bestMove = (nniMove *)malloc(sizeof(nniMove));
	bestMove->p = p;
	bestMove->nniType = bestType;

	/******************** NNI part *******************************/

	/* Restore the likelihood vector */
	//restoreLHVector(p,q, p_lhsave, q_lhsave);

	return bestMove;
}

int main(int argc, char * argv[])
{
	tree        * tr;
	if (argc != 2)
	{
		fprintf (stderr, "syntax: %s [binary-alignment-file]\n", argv[0]);
		return (1);
	}
	tr = (tree *)malloc(sizeof(tree));

	/* read the binary input, setup tree, initialize model with alignment */
	read_msa(tr,argv[1]);
	tr->randomNumberSeed = 665;

	/* Create random tree */

	makeRandomTree(tr);
	printf("RANDOM TREE: Number of taxa: %d\n", tr->mxtips);
	printf("RANDOM TREE: Number of partitions: %d\n", tr->NumberOfModels);

	/* compute the LH of the full tree */
	printf ("Virtual root: %d\n", tr->start->number);
	evaluateGeneric(tr, tr->start, TRUE);
	int printBranchLengths=TRUE;
	//Tree2String(tr->tree_string, tr, tr->start->back, printBranchLengths, FALSE, 0, 0, 0, SUMMARIZE_LH, 0,0);
	//fprintf(stderr, "%s\n", tr->tree_string);
	//printf("Likelihood: %f\n", tr->likelihood);

	/* Model optimization */
	modOpt(tr, 0.1);
	evaluateGeneric(tr, tr->start, TRUE);
	printf("Likelihood after model optimization: %f\n", tr->likelihood);

	/* 8 rounds of branch length optimization */
	//smoothTree(tr, 32);
	//evaluateGeneric(tr, tr->start, TRUE);
	//printf("Likelihood after branch length optimization: %f\n", tr->likelihood);
	printBranchLengths=TRUE;
	//Tree2String(tr->tree_string, tr, tr->start->back, printBranchLengths, FALSE, 0, 0, 0, SUMMARIZE_LH, 0,0);
	//fprintf(stderr, "%s\n", tr->tree_string);

	/* Initialize the NNI list */

	nniMove** nniList = malloc( (tr->ntips - 3) * sizeof(int));
	int i;
	for (i = 0; i < tr->ntips-2; i++) {
		nniList[i] = NULL;
	}

	/* fill up the NNI list */

	nodeptr p = tr->start->back;
	nodeptr q = p->next;
	//evalNNIForBranch(tr, q->back);
	int cnt = 1;
	while ( q != p ) {
		evalNNIForSubtree(tr, q->back, nniList, &cnt);
		q = q->next;
	}
	evaluateGeneric(tr, tr->start, TRUE);
	printf("Likelihood before evaluating all NNIs: %f\n", tr->likelihood);

	/*
	nodeptr p;
	do {
		p = pickRandomSubtree(tr);
	} while (isTip(p->back->number, tr->mxtips));
	printf("Trying to do NNI from branch %d - %d\n", p->number, p->back->number);
	printTopology(tr, TRUE);
	nniMove* bestMove = getBestNNI(tr, p);
	if (bestMove == NULL)
		printf("No improving NNI move found\n");
	else
		printf("Found improving NNI move of type: %d\n", bestMove->nniType);
	*/
	printTopology(tr, TRUE);
	return (0);
}

void evalNNIForSubtree(tree* tr, nodeptr p, nniMove** nniList, int* cnt) {
	printf("At node %d\n", p->number);
	if ( ! isTip(p->number, tr->mxtips) ) {
		newviewGeneric(tr, p, FALSE);
		newviewGeneric(tr, p->back, FALSE);
		nniMove* bestnni = evalNNIForBranch(tr, p);
		printf("cnt = %d\n", *cnt);
		*cnt = *cnt + 1;
		nniList[*cnt] = bestnni;
		nodeptr q = p->next;
		while ( q != p ) {
			evalNNIForSubtree(tr, q->back, nniList, cnt);
			q = q->next;
		}
	}
}




