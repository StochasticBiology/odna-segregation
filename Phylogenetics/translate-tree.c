// phrase a text tree in old NCBI format (not Newick!) as edges and nodes

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define LEN 200
#define MAXL 1000

void ReadTreeNCBI(char *names, int *n, int *src, int *dest, int *nc, char *fstr)
{
  FILE *fp;
  int ref;
  int lastlevel[10000];
  char str[1000];
  int level;
  int i, j;
  FILE *fpout1, *fpout2;
  
  fp = fopen(fstr, "r");
  if(fp == NULL)
    {
      printf("Couldn't read common tree file %s\n", fstr);
      exit(0);
    }
  ref = 0; lastlevel[0] = 0;

  sprintf(str, "%s-labels", fstr);
  fpout1 = fopen(str, "w");
  sprintf(str, "%s-edges", fstr);
  fpout2 = fopen(str, "w");
  
  // level stores level in tree
  // assign the previous node at one higher level to be this node's parent

  while(!feof(fp))
    {
      fgets(str, 1000, fp);
      if(feof(fp)) break;
      ref++;
      level = 0;
      for(i = 0; str[i] == '+' || str[i] == ' '; i++) level += (str[i] == '+');
      j = 0; for(; str[i] != '\n'; i++) names[ref*LEN+(j++)] = str[i];
      names[ref*LEN+j-1] = '\0';
      //  fprintf(fp2, "%i %i %s", ref, level, &str[i]);
      lastlevel[level] = ref;
      if(level > 0)
	{
	  dest[*nc] = ref;
	  src[*nc] = lastlevel[level-1];
	  (*nc)++;
	  //  fprintf(fp1, "%i %i\n", ref, lastlevel[level-1]);
	}
    }
  *n = ref;
  fclose(fp);
  //  fclose(fp1);
  //fclose(fp2);
  // output
  printf("\n{");
  for(i = 1; i < *n; i++)
    fprintf(fpout1, "%s\n", &names[LEN*i]);

  for(i = 0; i < *nc; i++)
    fprintf(fpout2, "%i %i\n", src[i], dest[i]);

  printf("%i nodes overall\n", *n);

  fclose(fpout1);
  fclose(fpout2);
  
}

int main(int argc, char *argv[])
{
  FILE *fp, *fp1;
  int change;
  int refs[10000];
  int line;
  char str[10000];
  int *matrix;
  int *traits;
  int foundtraits[10000];
  int cando;
  int *src, *dest;
  char *names;
  int n, nc;
  int i, j, k;
  int ref[1000], ref1, ref2, nref;
  int leaf;
  char fstr[200];
  int L;
  
  if(argc != 2)
    {
      printf("Usage: ./tree.ce [common tree]\n");
      return 0;
    }

  matrix = (int*)malloc(sizeof(int)*100000*MAXL);
  traits = (int*)malloc(sizeof(int)*100000*MAXL);
  names = (char*)malloc(sizeof(char)*100000*LEN);
  src = (int*)malloc(sizeof(int)*1000000);
  dest = (int*)malloc(sizeof(int)*1000000);

  ReadTreeNCBI(names, &n, src, dest, &nc, argv[1]);

  return 0;
}
