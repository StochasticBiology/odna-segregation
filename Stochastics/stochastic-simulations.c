// stochastic simulation to characterise the influence of different mechanisms on heteroplasmy variance

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND drand48()

#define MAXT 10       // timescale of simulation (used for memory allocation)
#define NIT 5000      // number of instances of stochastic simulation used to build up statistics
#define NUMR 24       // number of individual processes (reactions) affecting mtDNA
#define MAXPHASE 100  // number of possible dynamic phases

// structure storing a population state (w wildtype, m mutant mtDNAs)
typedef struct 
{
  int ws, wf, ms, mf;
} Population;

// structure storing parameters for a particular dynamic model
typedef struct 
{
  int divisions, gillespie;
  int cluster;
  int lowsize;
  int rsample, ramplify;
  
  double alphafuse, alphafrag;
  double alpha;
  double kappa;
  double beta1, beta2;
  double delta;
  double epsilon;
  int nd, ndf;
  double lambda, nu;      // replication and degradation rates
  double tau;             // timescale for this phase
  double nuf;
} Model;

// make a new population "dest" by taking a sample of size "endsize" from population "src"
// "random": 0 for deterministic; 1 for sampling without replacement (hypergeometric); 2 for sampling with replacement (binomial), 3 hypergeom with binomial end number; 4 binomial with binomial end number
void Subsample(Population src, Population *dest, int endsize, int random, int cluster)
{
  int i;
  double h;
  int nelements, effectiveendsize;
  double f;

  f = (double)(src.ws+src.ms)/(src.wf+src.mf+src.ws+src.ms);
  h = (double)(src.ms+src.mf)/(src.wf+src.mf+src.ws+src.ms);

  src.wf += src.ws; src.ws = 0;
  src.mf += src.ms; src.ms = 0;
  
  // initialise empty new population
  dest->wf = dest->mf = dest->ws = dest->ms = 0;

  if(random == 0)
    {
      // deterministic sampling
      dest->mf = endsize*h*(1-f);
      dest->ms = endsize*h*f;
      dest->wf = endsize*(1-h)*(1-f);
      dest->ws = endsize*(1-h)*f;
      return;
    }

  // lump into clusters
  src.wf /= cluster; src.mf /= cluster;

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
      h = (double)src.mf/(src.wf+src.mf);
      
      if(RND < h)
	{
	  // we've sampled a mutant
	  dest->mf++;
	  if(random%2 == 1)
	    src.mf--;
	}
      else
	{
	  // we've sampled a wildtype
	  dest->wf++;
	  if(random%2 == 1)
	    src.wf--;
	}
    }

  // build back from clusters
  dest->wf *= cluster; dest->mf *= cluster;

  dest->ws = dest->wf*f; dest->wf *= (1-f);
  dest->ms = dest->mf*f; dest->mf *= (1-f);
}

// make a new population "dest" by running a Polya urn from original population "src" until we get to size "endsize"
void Amplify(Population src, Population *dest, int endsize, int random)
{
  int i;
  double h;
  double n;
  double f;

  n = (src.wf+src.mf+src.ws+src.ms);
  
  // initialise new population identical to old
  dest->wf = src.wf; dest->mf = src.mf;  dest->ws = src.ws; dest->ms = src.ms;

  if(random == 0)
    {
      // deterministic amplification
      dest->wf = src.wf*((double)endsize/n);
      dest->mf = src.mf*((double)endsize/n);
      dest->ws = src.ws*((double)endsize/n);
      dest->ms = src.ms*((double)endsize/n);
      return;
    }
      
  // build up population one at a time
  for(; dest->wf+dest->mf+dest->ws+dest->ms < endsize;)
    {
      // get current heteroplasmy (changes as the urn fills)
      h = (double)(dest->mf+dest->ms)/(dest->wf+dest->mf+dest->ws+dest->ms);
      f = (double)(dest->ws+dest->ms)/(dest->wf+dest->mf+dest->ws+dest->ms);
      
      if(RND < h)
	{
	  if(RND < f)
	    dest->ms++;  // sampled a mutant
	  else
	    dest->mf++;
	}
      else
	{
	  if(RND < f)
	    dest->ws++;  // sampled a wildtype
	  else
	    dest->wf++;
	}
    }
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

  P.ws = src.ws; P.wf = src.wf; P.ms = src.ms; P.mf = src.mf;
   
  // start timer and run until "M.tau"
  for(t = 0; t < M.tau; )
    {
      // assign rates for each process
      rates[0] = P.wf*(M.lambda+M.delta)*(1.-M.alpha*(P.mf+P.ms+P.wf+P.ws));
      rates[1] = P.mf*M.lambda*(1.-M.alpha*(P.mf+P.ms+P.wf+P.ws));
      rates[2] = P.wf*M.nuf/M.ndf;
      rates[3] = P.mf*M.nuf/M.ndf;
      rates[4] = P.ws*M.beta1*(1.-M.alpha*(P.mf+P.ms+P.wf+P.ws))*(M.lambda+M.delta);
      rates[5] = P.ms*M.beta1*(1.-M.alpha*(P.mf+P.ms+P.wf+P.ws))*M.lambda;
      rates[6] = P.ws*(1.-M.beta1)*(1.-M.alpha*(P.mf+P.ms+P.wf+P.ws))*(M.lambda+M.delta);
      rates[7] = P.ms*(1.-M.beta1)*(1.-M.alpha*(P.mf+P.ms+P.wf+P.ws))*M.lambda;
      rates[8] = P.ws*M.nu/M.nd;
      rates[9] = P.ms*M.nu/M.nd;

      rates[10] = P.ws*(P.ws-1)/2.*M.alphafuse;
      rates[11] = P.wf*P.ws*M.alphafuse;
      rates[12] = P.ms*(P.ms-1)/2.*M.alphafuse;
      rates[13] = P.mf*P.ms*M.alphafuse;
      rates[14] = P.wf*M.beta2*M.alphafrag;
      rates[15] = P.mf*M.beta2*M.alphafrag;
      rates[16] = P.wf*(P.wf-1)/2.*(1.-M.beta2)*M.alphafrag;
      rates[17] = P.mf*(P.mf-1)/2.*(1.-M.beta2)*M.alphafrag;
      rates[18] = P.wf*P.mf*(1.-M.beta2)*M.alphafrag;
      rates[19] = P.wf*P.ms*M.alphafuse;
      rates[20] = P.ws*P.mf*M.alphafuse;
      rates[21] = P.ws*P.ms*M.alphafuse;
      
      rates[22] = P.wf*P.mf*M.kappa;
      rates[23] = P.wf*P.mf*(M.kappa+M.epsilon);

      for(i = 0; i < NUMR; i++)
	{
	  if(rates[i] < 0)
	    rates[i] = 0;
	}
      
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
	case 0: P.wf+=1; break;
	case 1: P.mf+=1; break;
	case 2: P.wf-=M.ndf; break;
	case 3: P.mf-=M.ndf; break;
	case 4: P.ws+=1; break;
	case 5: P.ms+=1; break;
	case 6: P.ws-=1; P.wf+=2; break;
	case 7: P.ms-=1; P.mf+=2; break;
	case 8: P.ws-=M.nd; break;
	case 9: P.ms-=M.nd; break;
	  
	case 10:P.ws-=2; P.wf+=2; break;
	case 11:P.ws-=1; P.wf+=1; break;
	case 12:P.ms-=2; P.mf+=2; break;
	case 13:P.ms-=1; P.mf+=1; break;
	case 14:P.wf-=1; P.ws+=1; break;
	case 15:P.mf-=1; P.ms+=1; break;
	case 16:P.wf-=2; P.ws+=2; break;
	case 17:P.mf-=2; P.ms+=2; break;
	case 18:P.wf-=1; P.mf-=1; P.ws+=1; P.ms+=1; break;	  
	case 19:P.ms-=1; P.mf+=1; break;
	case 20:P.ws-=1; P.wf+=1; break;
	case 21:P.ws-=1; P.ms-=1; P.wf+=1; P.mf+=1; break;	  

	case 22:P.wf-=1; P.mf+=1; break;
	case 23:P.wf+=1; P.mf-=1; break;
	}

      // update timer
      dt = -log(RND)/cumsum[NUMR-1];
      t += dt;
    }

  // return final population
  dest->ws = P.ws; dest->wf = P.wf; dest->ms = P.ms; dest->mf = P.mf;
}

double Predicth(Model M, double n, double f, double t)
{
  double gamma1, gamma2;
  double predict;
  double h0 = 0.5;
  
  gamma1 = (n-1)*(M.alpha*n - 1)*M.delta;
  gamma2 = (1.-f)*(1.-f)*n*n*M.epsilon;

  predict = 1./(1 + (1./h0-1)*exp((-gamma1/n + gamma2/n)*t));
}
  
double PredictNeutral(Model M, int popn, int compartments, double f, double t)
{
  double predict = 0;
  
  if(M.divisions)
    {
      if(M.rsample == 1)
	{
	  predict += M.cluster*(1./M.lowsize - 1./popn);
	}
      if(M.ramplify == 1)
	{
	  predict += (1./M.lowsize - 1./popn);
	}
    }
  if(M.gillespie)
    {
      if(compartments)
	{
	  predict += t*(f*M.nu*(1+M.nd)/popn + 2.*(1.-f)*(1.-f)*M.kappa);
	}
      else
	{
	  predict += t*((1+M.nd)*M.nu/popn + 2.*M.kappa);
	}
    }

  return predict;

}

double PredictSelection(Model M, double n, double f, double t)
{
  double rho1, rho2, b;
  double gamma1, gamma2;
  double predict;
  
  gamma1 = (n-1)*(M.alpha*n - 1)*M.delta;
  gamma2 = (1.-f)*(1.-f)*n*n*M.epsilon;

  b = (-gamma1+gamma2)/n;
  
  rho1 = (-1+n*M.alpha)*M.delta - (-1+f)*(-1+f)*n*(M.epsilon+2*M.kappa) - M.lambda + n*M.alpha*M.lambda - f*M.nd*M.nu;
  rho2 = ((-1+f)*(-1+f)*n*(M.epsilon+2*M.kappa) + M.lambda - n*M.alpha*M.lambda + f*M.nd*M.nu);

  predict = 2*exp(2./(1+exp(b*t)))*rho1 + 2*exp(2./(1+exp(b*t)) + b*t)*rho1 + exp(1+b*t)*(M.delta - n*M.alpha*M.delta + 2*rho2) + exp(1)*((3 - 3*n*M.alpha)*M.delta + 2*rho2);
  predict /= -4*exp(1)*(1+exp(b*t))*(gamma1-gamma2);
  
  return predict;

}

int main(void)
{
  Population P, newP, initP;
  int m0 = 500, w0 = 500;
  Model M;
  double *h, *n, *f;
  int i, t;
  double meanh, varh, meann, varn, meanf, varf;
  FILE *fp;
  char str[100];
  int expt;
  double predict, predicth, predicts;
  int popn;
  double popf, tpopn;
  double nstart;
  int equilibrate;
  double TEQ = 2;
  
  // allocate memory to store heteroplasmy values
  h = (double*)malloc(sizeof(double)*NIT*MAXT);
  n = (double*)malloc(sizeof(double)*NIT*MAXT);
  f = (double*)malloc(sizeof(double)*NIT*MAXT);

  // loop through "NIT" instances of the stochastic simulation algorithm

  sprintf(str, "stochastic-simulation-output.txt");
  printf("%s\n", str);
  fp = fopen(str, "w");

  for(expt = 0; expt <= 12; expt++)
    {

      for(popn = 500; popn <= 2000; popn += 500)
	{
	  // default params
	  M.beta1 = 0; M.beta2 = 0; M.delta = 0; M.alphafuse = 0; M.alphafrag = 0;
	  M.lambda = 2; M.nu = 1; M.nuf = 1; M.alpha = 1./1000; M.kappa = 0.002;  M.tau = 1;
	  M.ndf = M.nd = 1;
	  M.epsilon = 0;
  
	  M.gillespie = M.divisions = 0;
	  M.cluster = 1; M.rsample = 1; M.ramplify = 1;

	  equilibrate = 0;
      
	  M.alpha = 1./popn;
	  //nstart = (M.lambda-M.nu)/(M.lambda*M.alpha);
	  switch(expt)
	    {
	    case 0: break;
	    case 1: M.divisions = 1; M.ramplify = 0; M.lowsize = popn/2.; M.cluster = 1; break;
	    case 2: M.divisions = 1; M.lowsize = popn/2.; M.cluster = 1; break;
	    case 3: M.divisions = 1; M.lowsize = popn/2.; M.cluster = 5; break;
	    case 4: M.divisions = 1; M.lowsize = popn/4.; M.cluster = 1; break;
	    case 5: M.gillespie = 1; M.nuf = 0; M.alphafuse = 1; break;
	    case 6: M.gillespie = 1; M.nuf = 0; M.alphafrag = 1; break;
	    case 7: M.gillespie = 1; M.nuf = 0; M.alphafuse = 0.005; M.alphafrag = 0.01; break;
	    case 8: M.gillespie = 1; equilibrate = 1; M.nuf = 0; M.alphafuse = 0.005; M.alphafrag = 0.01; M.delta = -0.1; break;
	    case 9: M.gillespie = 1; equilibrate = 1; M.nuf = 0; M.alphafuse = 0.005; M.alphafrag = 0.01; M.epsilon = -0.0001; break;
	    case 10: M.gillespie = 1; equilibrate = 1; M.nuf = 0; M.alphafuse = 0.005; M.alphafrag = 0.01; M.delta = -0.1; M.epsilon = -0.0001; break;
	    case 11: M.gillespie = 1; M.nuf = 0; M.alphafuse = 0.005; M.alphafrag = 0.01; M.nd = 2; break;
	    case 12: M.gillespie = 1; M.nuf = 0; M.alphafuse = 0.005; M.alphafrag = 0.01; M.nd = 4; break;
	    }
    
	  printf("%i %i\n", expt, popn);
	  // loop through "NIT" instances of the stochastic simulation algorithm
	  for(i = 0; i < NIT; i++)
	    {
	      P.wf = popn/2.; P.mf = popn/2.; P.ws = 0; P.ms = 0;

	      if(equilibrate)
		{
		  for(t = 0; t <= TEQ; t++)
		    {
		      if(M.gillespie)
			{
			  Gillespie(P, &newP, M);
			  P = newP;
			}
		    }
		  tpopn = P.wf+P.mf+P.ws+P.ms;
		  popf = P.ms+P.ws;
		  popf /= tpopn;

		  P.wf = (1-popf)*tpopn/2; P.mf = (1-popf)*tpopn/2; P.ws = popf*tpopn/2; P.ms = popf*tpopn/2;
		}
	      
	      // loop through time
	      for(t = 0; t <= MAXT; t++)
		{
  		  // record heteroplasmy
		  h[i*MAXT+t] = (double)(P.mf+P.ms)/(P.wf+P.ws+P.mf+P.ms);
		  n[i*MAXT+t] = (double)(P.wf+P.ws+P.mf+P.ms);
		  f[i*MAXT+t] = (double)(P.ms+P.ws)/(P.wf+P.ws+P.mf+P.ms);

		  if(M.gillespie)
		    {
		      Gillespie(P, &newP, M);
		      P = newP;
		    }
		  if(M.divisions)
		    {
		      Subsample(P, &newP, M.lowsize, M.rsample, M.cluster);
		      P = newP;
		      Amplify(P, &newP, popn, M.ramplify);
		      P = newP;
		    }
		}
	    }

	  // loop through simulation timepoints
	  t = 0;
	  for(t = 0; t < MAXT; t++)
	    {
	      // initialise statistics
	      meanh = varh = meann = varn = meanf = varf = 0;

	      // build up mean heteroplasmy, looping through stochastic instances
	      for(i = 0; i < NIT; i++)
		{
		  meanh += h[i*MAXT+t];
		  meann += n[i*MAXT+t];
		  meanf += f[i*MAXT+t];
		}
	      meanh /= NIT;
	      meann /= NIT;
	      meanf /= NIT;

	      // build up heteroplasmy variance, looping through stochastic instances
	      for(i = 0; i < NIT; i++)
		{
		  varh += (meanh-h[i*MAXT+t])*(meanh-h[i*MAXT+t]);
		  varn += (meann-n[i*MAXT+t])*(meann-n[i*MAXT+t]);
		  varf += (meanf-f[i*MAXT+t])*(meanf-f[i*MAXT+t]);
		}
	      varh /= (NIT-1);
	      varn /= (NIT-1);
	      varf /= (NIT-1);

	      // output to file


	      if(M.delta != 0 || M.epsilon != 0)
		predict = PredictSelection(M, meann, meanf, t);
	      else
		predict = PredictNeutral(M, popn, (M.alphafuse != 0 || M.alphafrag != 0), meanf, t)*meanh*(1.-meanh);
	      predicth = Predicth(M, meann, meanf, t);

	      fprintf(fp, "%i %i %i %f %f %f %f %f %f %f %f %f %f\n", expt, popn, t, meanh, varh, varh/(meanh*(1.-meanh)), meann, varn, meanf, varf, predicth, predict, predict/(meanh*(1.-meanh)));
	    }
	}
    }
  fclose(fp);

  return 0;
}

  
