/* from utils.c */
void getDataTypeString(pllInstance *tr, pInfo *partitionInfo, char typeOfData[1024])
{
  switch(partitionInfo->dataType)
  {
    case AA_DATA:
      strcpy(typeOfData,"AA");
      break;
    case DNA_DATA:
      strcpy(typeOfData,"DNA");
      break;
    case BINARY_DATA:
      strcpy(typeOfData,"BINARY/MORPHOLOGICAL");
      break;
    case SECONDARY_DATA:
      strcpy(typeOfData,"SECONDARY 16 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case SECONDARY_DATA_6:
      strcpy(typeOfData,"SECONDARY 6 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case SECONDARY_DATA_7:
      strcpy(typeOfData,"SECONDARY 7 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case GENERIC_32:
      strcpy(typeOfData,"Multi-State");
      break;
    case GENERIC_64:
      strcpy(typeOfData,"Codon"); 
      break;
    default:
      assert(0);
  }
}


/* removed the static keyword for using this function in the examples */
static boolean setupTree (pllInstance *tr, boolean doInit, partitionList *partitions)
{
  nodeptr  p0, p, q;
  int
    i,
    j;

  int
    tips,
    inter; 

  if(doInit)
    init_default(tr);

  tr->bigCutoff = PLL_FALSE;

  tr->maxCategories = MAX(4, tr->categories);

  tips  = (size_t)tr->mxtips;
  inter = (size_t)(tr->mxtips - 1);

  tr->treeStringLength = tr->mxtips * (PLL_NMLNGTH + 128) + 256 + tr->mxtips * 2;

  tr->tree_string  = (char*)rax_calloc((size_t)tr->treeStringLength, sizeof(char)); 
  tr->tree0 = (char*)rax_calloc((size_t)tr->treeStringLength, sizeof(char));
  tr->tree1 = (char*)rax_calloc((size_t)tr->treeStringLength, sizeof(char));


  /*TODO, must that be so long ?*/

  tr->td[0].count = 0;
  tr->td[0].ti    = (traversalInfo *)rax_malloc(sizeof(traversalInfo) * (size_t)tr->mxtips);
  tr->td[0].executeModel = (boolean *)rax_malloc(sizeof(boolean) * (size_t)NUM_BRANCHES);
  tr->td[0].parameterValues = (double *)rax_malloc(sizeof(double) * (size_t)NUM_BRANCHES);

  tr->fracchange = -1.0;

  tr->constraintVector = (int *)rax_malloc((2 * (size_t)tr->mxtips) * sizeof(int));

  tr->nameList = (char **)rax_malloc(sizeof(char *) * (tips + 1));


  p0 = (nodeptr)rax_malloc((tips + 3 * inter) * sizeof(node));
  assert(p0);

  tr->nodeBaseAddress = p0;


  tr->nodep = (nodeptr *) rax_malloc((2* (size_t)tr->mxtips) * sizeof(nodeptr));
  assert(tr->nodep);    

  tr->nodep[0] = (node *) NULL;    /* Use as 1-based array */

  for (i = 1; i <= tips; i++)
  {
    p = p0++;

    p->hash   =  KISS32(); /* hast table stuff */
    p->x      =  0;
    p->xBips  = 0;
    p->number =  i;
    p->next   =  p;
    p->back   = (node *)NULL;
    p->bInf   = (branchInfo *)NULL;            
    tr->nodep[i] = p;
  }

  for (i = tips + 1; i <= tips + inter; i++)
  {
    q = (node *) NULL;
    for (j = 1; j <= 3; j++)
    {	 
      p = p0++;
      if(j == 1)
      {
        p->xBips = 1;
        p->x = 1;
      }
      else
      {
        p->xBips = 0;
        p->x =  0;
      }
      p->number = i;
      p->next   = q;
      p->bInf   = (branchInfo *)NULL;
      p->back   = (node *) NULL;
      p->hash   = 0;       
      q = p;
    }
    p->next->next->next = p;
    tr->nodep[i] = p;
  }

  tr->likelihood  = PLL_UNLIKELY;
  tr->start       = (node *) NULL;  

  tr->ntips       = 0;
  tr->nextnode    = 0;

  for(i = 0; i < NUM_BRANCHES; i++)
    tr->partitionSmoothed[i] = PLL_FALSE;

  tr->bitVectors = (unsigned int **)NULL;

  tr->vLength = 0;

  tr->h = (hashtable*)NULL;

  tr->nameHash = initStringHashTable(10 * tr->mxtips);

  for (i = 0; i < partitions->numberOfPartitions; i++) {
	partitions->partitionData[i] = (pInfo*)rax_malloc (sizeof(pInfo));
	partitions->partitionData[i]->partitionContribution = -1.0;
	partitions->partitionData[i]->partitionLH = 0.0;
	partitions->partitionData[i]->fracchange = 1.0;
  }

  return PLL_TRUE;
}


boolean modelExists(char *model, pllInstance *tr)
{
  /********** BINARY ********************/

  if(strcmp(model, "PSR") == 0)
  {
    tr->rateHetModel = CAT;
    return PLL_TRUE;
  }

  if(strcmp(model, "GAMMA") == 0)
  {
    tr->rateHetModel = GAMMA;
    return PLL_TRUE;
  }


  return PLL_FALSE;
}




