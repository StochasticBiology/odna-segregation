// stochastic simulation to characterise the influence of different mechanisms on heteroplasmy variance

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND drand48()

#define MAXT 30       // timescale of simulation (used for memory allocation)
#define NIT 5000      // number of instances of stochastic simulation used to build up statistics
#define NUMR 6        // number of individual processes (reactions) affecting mtDNA
#define MAXPHASE 100  // number of possible dynamic phases

// structure storing a population state (w wildtype, m mutant mtDNAs)
typedef struct 
{
  int w, m;
} Population;

// structure storing parameters for a particular dynamic model
typedef struct 
{
  int index;              // 0: no feedback (no use for nstar); 1: feedback control using nstar
  double lambda, nu;      // replication and degradation rates
  double kappa;           // gene conversion rate
  double nstar;           // target copy number for feedback control
  double tau;             // timescale for this phase
} Model;

  
// make a new population "dest" by taking a sample of size "endsize" from population "src"
// "random": 0 for deterministic; 1 for sampling without replacement (hypergeometric); 2 for sampling with replacement (binomial)
void Subsample(Population src, Population *dest, int endsize, int random, int cluster)
{
  int i;
  double h;
  int nelements, effectiveendsize;
  
  // initialise empty new population
  dest->w = dest->m = 0;

  // "random" controls partitioning regime:
  // 0: deterministic partitioning; 1: hypergeometric; 2: binomial; 3: hypergeometric with binomial choice of daughter copy number; 4: binomial with binomial choice of daughter copy number
  if(random == 0)
    {
      // get current heteroplasmy 
      h = (double)src.m/(src.w+src.m);

      // deterministic sampling
      dest->m = endsize*h;
      dest->w = endsize*(1-h);
      return;
    }

  // lump into clusters
  src.w /= cluster; src.m /= cluster;

  // deal with binomial copy number choice
  if(random >= 3)
    {
      nelements = endsize/cluster*2;
      effectiveendsize = 0;
      for(i = 0; i < nelements; i++)
	{
	  if(RND < 0.5)
	    effectiveendsize++;
	}
      endsize = effectiveendsize*cluster;
    }
  
  // build up new population one at a time
  for(i = 0; i < endsize/cluster; i++)
    {
      // get current heteroplasmy (this will change if we're sampling without replacement)
      h = (double)src.m/(src.w+src.m);
      
      if(RND < h)
	{
	  // we've sampled a mutant
	  dest->m++;
	  if(random%2 == 1)
	    src.m--;
	}
      else
	{
	  // we've sampled a wildtype
	  dest->w++;
	  if(random%2 == 1)
	    src.w--;
	}
    }

  // build back from clusters
  dest->w *= cluster; dest->m *= cluster;
}

// make a new population "dest" by running a Polya urn from original population "src" until we get to size "endsize"
void Amplify(Population src, Population *dest, int endsize, int random)
{
  int i;
  double h;

  // initialise new population identical to old
  dest->w = src.w; dest->m = src.m;

  if(random == 0)
    {
      // deterministic amplification
      dest->w = src.w*((double)endsize/(src.w+src.m));
      dest->m = src.m*((double)endsize/(src.w+src.m));
      return;
    }
      
  // build up population one at a time
  for(; dest->w+dest->m < endsize;)
    {
      // get current heteroplasmy (changes as the urn fills)
      h = (double)dest->m/(dest->w+dest->m);
      
      if(RND < h)
        dest->m++;  // sampled a mutant
      else
	dest->w++;  // sampled a wildtype
    }
}  

// return replication rate given model "M" and population state "P"
double Lambda(Population P, Model M)
{
  if(M.index == 0) return M.lambda;                 // non-feedback case
  else return M.lambda*(1. - (P.w+P.m)/M.nstar);    // feedback control based on population size
}

// return degradation rate given model "M" and population state "P"
double Nu(Population P, Model M)
{
  return M.nu;
}

// return gene conversion rate given model "M" and population state "P"
double Kappa(Population P, Model M)
{
  return M.kappa;
}

// run Gillespie algorithm to apply individual mtDNA population dynamics to a population
// start with state "src" and end at state "dest" after time "M.tau"
// model "M" describes model structure and parameters
void Gillespie(Population src, Population *dest, Model M)
{
  double t, dt;                       // timing variables
  double r;                           // roulette ball for choosing a reaction
  double rates[NUMR], cumsum[NUMR];   // rates and cumulative sum of rates
  Population P;                       // evolving population
  int i;                              // counter

  // initialise population "P" identical to "src"
  P.w = src.w; P.m = src.m;

  // start timer and run until "M.tau"
  for(t = 0; t < M.tau; )
    {
      // assign rates for each process
      rates[0] = Lambda(P, M)*P.w;          // W -> W+W
      rates[1] = Nu(P, M)*P.w;              // W -> 0
      rates[2] = Lambda(P, M)*P.m;          // M -> M+M
      rates[3] = Nu(P, M)*P.m;              // M -> 0
      rates[4] = Kappa(P, M)*P.w*P.m;       // W + M -> M + M
      rates[5] = Kappa(P, M)*P.w*P.m;       // W + M -> W + W

      // build cumulative sum of rates
      cumsum[0] = rates[0];
      for(i = 1; i < NUMR; i++)
	cumsum[i] = cumsum[i-1]+rates[i];

      // roll roulette ball
      r = RND*cumsum[NUMR-1];

      // choose reaction corresponding to roulette roll
      for(i = 0; cumsum[i] < r; i++);

      // apply population changes based on which reaction was chosen
      switch(i)
	{
	case 0: P.w++; break;
	case 1: P.w--; break;
	case 2: P.m++; break;
	case 3: P.m--; break;
	case 4: P.w--; P.m++; break;
	case 5: P.w++; P.m--; break;
	}

      // update timer
      dt = -log(RND)/cumsum[NUMR-1];
      t += dt;
    }

  // return final population
  dest->w = P.w; dest->m = P.m;
}

int main(void)
{
  Population P, newP;
  int m0 = 500, w0 = 500;
  Model M;
  double *h;
  int i, t;
  double mean, var;
  FILE *fp;
  int gillespie, divisions, rsample, ramplify;
  char str[100];
  int cluster;
  int lowsize;
  
  // allocate memory to store heteroplasmy values
  h = (double*)malloc(sizeof(double)*NIT*MAXT);

  // parameters for simulation
  // "gillespie" [0-6]: different model structures for cellular turnover. 0: no turnover; [1-6]: see switch below. M.index determines feedback (0 none; 1 relaxed replication); other M.* parameters follow accordingly
  // "divisions" [0-14]: divisions % 5 gives cluster sizes, divisions / 5 gives population size of daughter cell (see switches below)
  // "rsample" [0-4]: 0: deterministic partitioning; 1: hypergeometric; 2: binomial; 3: hypergeometric with binomial choice of daughter copy number; 4: binomial with binomial choice of daughter copy number
  // "ramplify" [0-1]: 0: deterministic amplification; 1: random (Polya) reamplification

  // we will only simulate a subset of possibilities here for the illustration plot; others can be (and have been) investigated too!

  for(gillespie = 0; gillespie <= 0; gillespie++)
    {
      for(divisions = 0; divisions <= 14; divisions++)
	{
	  for(rsample = 0; rsample <= 2; rsample++)
	    {
	      for(ramplify = 0; ramplify <= 1; ramplify++)
		{
		  // interpret simulation parameters as model parameters
		  M.tau = 1;
		  switch(gillespie)
		    {
		    case 1: M.index = 0; M.lambda = 0.5; M.nu = 0.5; M.kappa = 0; break;
		    case 2: M.index = 0; M.lambda = 0; M.nu = 0; M.kappa = 0.001; break;
		    case 3: M.index = 0; M.lambda = 0.5; M.nu = 0.5; M.kappa = 0.001; break;
		    case 4: M.index = 1; M.nstar = 2000; M.lambda = 1; M.nu = 0.5; M.kappa = 0; break;
		    case 5: M.index = 1; M.nstar = 2000; M.lambda = 0; M.nu = 0; M.kappa = 0.001; break;
		    case 6: M.index = 1; M.nstar = 2000; M.lambda = 1; M.nu = 0.5; M.kappa = 0.001; break;
		    }

		  switch(divisions % 5)
		    {
		    case 1: cluster = 1; break;
		    case 2: cluster = 2; break;
		    case 3: cluster = 4; break;
		    case 4: cluster = 16; break;
		    }

		  switch(divisions/5)
		    {
		    case 0: lowsize = 500; break;
		    case 1: lowsize = 800; break;
		    case 2: lowsize = 200; break;
		    }
		  
		  printf("%i %i %i %i\n", gillespie, divisions, rsample, ramplify);
	      
		  // loop through "NIT" instances of the stochastic simulation algorithm
		  for(i = 0; i < NIT; i++)
		    {
		      // initialise our working population "P"
		      P.m = m0; P.w = w0;
		      // loop through time
		      for(t = 0; t < MAXT; t++)
			{
			  // record heteroplasmy
			  h[i*MAXT+t] = (double)P.m/(P.w+P.m);
			  if(gillespie)
			    {
			      Gillespie(P, &newP, M);
			      P = newP;
			    }
			  if(divisions)
			    {
			      Subsample(P, &newP, lowsize, rsample, cluster);
			      P = newP;
			      Amplify(P, &newP, 1000, ramplify);
			      P = newP;
			    }
			}
		    }

		  // open file for output
		  sprintf(str, "simple-output-%i-%i-%i-%i.txt", gillespie, divisions, rsample, ramplify);
		  fp = fopen(str, "w");

		  // loop through simulation timepoints
		  for(t = 0; t < MAXT; t++)
		    {
		      // initialise statistics
		      mean = var = 0;

		      // build up mean heteroplasmy, looping through stochastic instances
		      for(i = 0; i < NIT; i++)
			mean += h[i*MAXT+t];
		      mean /= NIT;

		      // build up heteroplasmy variance, looping through stochastic instances
		      for(i = 0; i < NIT; i++)
			var += (mean-h[i*MAXT+t])*(mean-h[i*MAXT+t]);
		      var /= (NIT-1);

		      // output to file
		      fprintf(fp, "%i %f %f %f\n", t, mean, var, var/(mean*(1.-mean)));
		    }
		  fclose(fp);
		}
	    }
	}
    }
  
  return 0;
}

  
