# NCBI returns text trees where some sub-species layers awkwardly lie far beneath the most common species layer
# given that we are just visualising a taxonomy, this code pads edge lengths to "catch up" with these deep layers, so all leaves appear on the same level
# once more, all of this is just for visualisation and it is recommended to use a neater pipeline with e.g. phytools or ETE3

awk -v maxl=$2 'BEGIN{currentdepth = 0;ghostnum = 0;}{
tmp = $0;
newdepth = gsub("+ ", "+ ", tmp);
if(newdepth < currentdepth)
{
  /* first print ghost ancestors */
  for(j = currentdepth; j < maxl; j++)
  {
    for(k = 0; k < j; k++)
    {
      printf("+ ");
    }
    printf("ghost%i\n", ghostnum++);
  }
  j = maxl;
  for(i = 0; i < ntoprint; i++)
  {
    for(k = 0; k < j-currentdepth; k++)
    {
      printf("+ ");
    }
    printf("%s\n", toprint[i]);
  }
  ntoprint = 1;
  toprint[0] = $0;
}
else if(newdepth == currentdepth && newdepth != 0)
{
  toprint[ntoprint] = $0;
  ntoprint++;
}
else 
{
  for(i = 0; i < ntoprint; i++)
  {
    printf("%s\n", toprint[i]);
  }
  ntoprint = 1;
  toprint[0] = $0;
}
currentdepth = newdepth;
}' $1 > $1.padded
  
