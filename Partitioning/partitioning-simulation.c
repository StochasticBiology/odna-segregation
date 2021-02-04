// simulation of genetics of oDNA inheritance given physical arrangement of oDNAs pre-division

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NMAX 2000    // max number of mitos for memory allocations
#define NIT 1000     // number of iterations
#define MAXT 200     // time for position equilibration

#define RND drand48()
#define SRND (RND-0.5)

int main(void)
{
  double N;
  double x[NMAX], y[NMAX], z[NMAX];
  double count[NIT], het[NIT];
  int i, j, k, t;
  double dx, dy, dz;
  double mean, sd;
  double meanh, sdh;
  double scale = 0.2, pairscale = 0.1;
  double d, norm;
  double closest;
  int closestref;
  int expt;
  char str[100];
  FILE *fp;

  // wipe output file
  fp = fopen("partitioning-stats-out.txt", "w"); fclose(fp);
  
  // loop through numbers of mitos
  for(N = 2; N < 1500; N*= 1.5)
    {
      // loop through interaction protocol
      for(expt = 0; expt <= 1; expt++)
	{
	  // stochastic iterations
	  for(i = 0; i < NIT; i++)
	    {
	      // initialise position of each mito
	      for(j = 0; j < N; j++)
		{
		  x[j] = y[j] = z[j] = 0;
		}
	      
	      // loop through simulation time
	      for(t = 0; t < MAXT; t++)
		{
		  // loop over mitos
		  for(j = 0; j < N; j++)
		    {
		      // apply a random position shift as long as it doesn't take us outside of the unit sphere
		      dx = SRND*scale; dy = SRND*scale; dz = SRND*scale;
		      if((x[j]+dx)*(x[j]+dx)+(y[j]+dy)*(y[j]+dy)+(z[j]+dz)*(z[j]+dz) < 1)
			{
			  x[j] += dx; y[j] += dy; z[j] += dz;
			}

		      // if we're interacting
		      if(expt == 1)
			{
			  // initialise overall shift
			  dx = dy = dz = 0;

			  // loop over other mitos
			  for(k = 0; k < N; k++)
			    {
			      if(k != j)
				{
				  // calculate distance, and increment overall displacement vector corresponding scaled by 1/d^2
				  d = (x[j]-x[k])*(x[j]-x[k]) + (y[j]-y[k])*(y[j]-y[k]) + (z[j]-z[k])*(z[j]-z[k]);
				  dx += (x[j]-x[k])/(d*d); dy += (y[j]-y[k])/(d*d); dz += (z[j]-z[k])/(d*d);
				}
			    }
			  
			  // normalise overall displacement vector
			  norm = sqrt(dx*dx + dy*dy + dz*dz);
			  dx *= pairscale/norm; dy *= pairscale/norm; dz *= pairscale/norm;

			  // apply this shift, if it keeps us in the unit sphere
			  if((x[j]+dx)*(x[j]+dx)+(y[j]+dy)*(y[j]+dy)+(z[j]+dz)*(z[j]+dz) < 1)
			    {
			      x[j] += dx; y[j] += dy; z[j] += dz;
			    }
			}
		    }
		}

	      // model partitioning with the x=0 division plane
	      // we're picturing mitos with ref < N/2 as mutants
	      count[i] = het[i] = 0;
	      for(j = 0; j < N; j++)
		{
		  if(x[j] > 0) count[i]++;
		  if(x[j] > 0 && j < N/2) het[i]++;
		}
	      het[i] /= count[i];
	      count[i] /= N;

	      // output co-ordinates as an example for the first simulation each time
	      if(i == 0)
		{
		  sprintf(str, "ex-%.0f-%i.txt", N, expt);
		  fp = fopen(str, "w");
		  for(j = 0; j < N; j++)
		    fprintf(fp, "%f %f %f\n", x[j], y[j], z[j]);
		  fclose(fp);
		}
	    }

	  // compute statistics
	  mean = sd = meanh = sdh = 0;
	  for(i = 0; i < NIT; i++)
	    {
	      mean += count[i];
	      meanh += het[i];
	    }
	  mean /= NIT;
	  meanh /= NIT;
	  for(i = 0; i < NIT; i++)
	    {
	      sd += (count[i]-mean)*(count[i]-mean);
	      sdh += (het[i]-meanh)*(het[i]-meanh);
	    }
	  sd = sqrt(sd/(NIT-1));
	  sdh = sqrt(sdh/(NIT-1));

	  fp = fopen("partitioning-stats-out.txt", "a"); 
	  fprintf(fp, "%.0f %i %f %f %f %f\n", N, expt, mean, sd, meanh, sdh);
	  fclose(fp);
	}
    }
  
  return 0;
}
  
  
  
