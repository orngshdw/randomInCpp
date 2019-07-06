#include "stdio.h"
#include "stdlib.h"
#include "math.h"

void   RandomInitialise(int,int);
double RandomUniform(void);
double RandomGaussian(double,double);
int    RandomInt(int,int);
double RandomDouble(double,double);

int main(int argc,char **argv)
{
    int tried2move, z, remadegrid, failure, porocity, i, j, k, p, v=1 /*particle velocity*/;
	int repeat, randpart, clogcheck, loc, moves;
	int percentclog; int clognum, times;
	int clogcheck_1, inte;
    int B[27][77], Ptcle_j[625], Ptcle_k[625];
    double e, accept;
    double r,rmin=1e32,rmax=-1e32;
    double log(double root);
    double pow (double x, double y);
    FILE *fp;

	fp = fopen("Clog Percentage 0 to 50%.txt", "w+");
	fprintf(fp, "porocity percentclog\n");

    /* Generate 20000 random numbers */
    RandomInitialise(1802,9373);
    for (i=0;i<20000;i++) {
       r = RandomUniform();
       if (r < rmin) rmin = r;
       if (r > rmax) rmax = r;
    }
    fprintf(stderr,"Numbers range from %g to %g\n",rmin,rmax);

	//changes the porocity
	for(porocity=0;porocity<=50;porocity++)
	{
		inte=0;
		failure=0;   	
		clognum=0;
		percentclog=0;
		for(times=1;times<=100;times++)//100 times for each porocity
		{
			//Sets grid to 0
			//j=1-25,k=1-25 is where the air particles are originally
			//j=1-25,k=26-50 is the porous filter
			//j=1-25,k=51-75 is the other of side of the filter
			for(j=1;j<=25;j++)
			{
				for(k=1;k<=75;k++) 
				{
					B[j][k]=0;
				}
			}

			//Creates a border numbered 95
			for(k=0;k<=76;k++)
			{
				B[0][k]=95;
				B[26][k]=95;
			}
			for(j=0;j<=26;j++)
			{
 				B[j][0]=95;
				B[j][76]=95;
			}

			//creating the particles; j=1-25,k=1-25
			p=0;
			for(j=1;j<=25;j++)
			{
				for(k=1;k<=25;k++)
				{
					r=(int)(101*RandomUniform());

					if(r<50)
					{
						B[j][k]=1;
						//stores coordinates of the particle
						Ptcle_j[p]=j;
						Ptcle_k[p]=k;
						p++;
						//total particles = p
					}
				}
			}
		    
			//creating a non-clogged porous grid; j=1-25,k=26-50
			remadegrid=0;
			failure=0;
			repeat=10;
			moves=1;
			while (repeat>0)
			{
				for(j=1;j<=25;j++)
				{
					for(k=26;k<=50;k++)
					{
						r=(int)(101*RandomUniform());

						if(r>=porocity)
						{
							B[j][k]=99;
						}
					}
				}
				
				//checks to see if the filter is clogged				
				for(k=26;k<=50;k++)//j=1-25,k=26-50 is the porous filter
				{
					clogcheck_1=0;
					for(j=1;j<=25;j++)
					{
						if(B[j][k]==0 && B[j][k+1]==0)
						{
							clogcheck_1++;
						}
					}
					if(clogcheck_1>0)//successful nonclog grid stops "while" loop
					{
						repeat=0;
					}
				}
				remadegrid++;
				if(remadegrid==100)//halts attempt to create nonclogged grid due to many failures
				{
					repeat=0;
					failure=1;
				}
				//printf("times=%d failure=%d\n", times, failure);
			}
		


			//moves the partricles, if creating the nonporous grid was sucessful
			moves=0;
			repeat=10;
			while(repeat>0 && failure==0)
			{
				//randomly selects a particle
				randpart=(int)((p+1)*RandomUniform());

				//denotes location of particle
				//particle has not entered filter=1;particle in filter=2;particle has left filter=3
				if(Ptcle_k[randpart]<=25)//j=1-25,k=1-25 is where the air particles are originally
				{
					loc=1;
				}
				if(Ptcle_k[randpart]>=26 && Ptcle_k[randpart]<=50)//j=1-25,k=26-50 is the porous filter
				{
					loc=2;
				}
				if(Ptcle_k[randpart]>=51)//j=1-25,k=51-75 is the other of side of the filter
				{
					loc=3;
				}

				//randomly selects a direction to move the particle. Particle must move.
				//0=left ; 1=up ; 2=forward, always accepted ; 3=down
				tried2move=0;
				while (repeat>0)
				{
					repeat=10;
					r=(int)(4*RandomUniform());
					accept=RandomUniform();
					e=pow(2.718281828459, -v);
					tried2move++;
					if(r==0 && B[Ptcle_j[randpart]][Ptcle_k[randpart]-1]==0 && accept>e)//left
					{
						//move particle to unoccupied space, record the postion, if it is within filter check dead end
						B[Ptcle_j[randpart]][Ptcle_k[randpart]]--;
						B[Ptcle_j[randpart]][Ptcle_k[randpart]-1]++;
						Ptcle_k[randpart]=Ptcle_k[randpart]-1;
						repeat=0;
					}
					if(r==1 && B[Ptcle_j[randpart]-1][Ptcle_k[randpart]]==0 && accept>e)//up
					{
						B[Ptcle_j[randpart]][Ptcle_k[randpart]]--;
						B[Ptcle_j[randpart]-1][Ptcle_k[randpart]]++;
						Ptcle_j[randpart]=Ptcle_j[randpart]-1;
						repeat=0;
					}
					if(r==2 && B[Ptcle_j[randpart]][Ptcle_k[randpart]+1]==0)//right
					{
						B[Ptcle_j[randpart]][Ptcle_k[randpart]]--;
						B[Ptcle_j[randpart]][Ptcle_k[randpart]+1]++;
						Ptcle_k[randpart]=Ptcle_k[randpart]+1;
						repeat=0;
					}
					if(r==3 && B[Ptcle_j[randpart]+1][Ptcle_k[randpart]]==0 && accept>e)//down
					{
						B[Ptcle_j[randpart]][Ptcle_k[randpart]]--;
						B[Ptcle_j[randpart]+1][Ptcle_k[randpart]]++;
						Ptcle_j[randpart]=Ptcle_j[randpart]+1;
						repeat=0;
					}
					if(loc==2 && B[Ptcle_j[randpart]+1][Ptcle_k[randpart]]>0
					&& B[Ptcle_j[randpart]][Ptcle_k[randpart]+1]>0
					&& B[Ptcle_j[randpart]][Ptcle_k[randpart]-1]>0
					&& B[Ptcle_j[randpart]-1][Ptcle_k[randpart]]>0)//ends if the particle is trapped
					{
						repeat=0;
					}
					if(loc==3 && Ptcle_k[randpart]>=71)//particle is already out of filter
					{
						repeat=0;
					}
					if(tried2move==10000)
					{
						repeat=0;
					}
				}//repeat=0 if loop is done

				//checks to see if there are any particles on the left side 
				repeat=0;
				for(j=1;j<=25;j++)
				{
					for(k=1;k<=25;k++)
					{
						if(B[j][k]==1)
						{
							repeat++;
						}
					}
				}//repeat=0 if all clear. else repeat>0

				//total number of times attempting to move particles across porous grid
				if(moves==1000)
				{
					repeat=0;
					clogcheck=0;
					for(z=0;z<p;z++)
					{
						//j=1-25,k=51-75 is the other of side of the filter
						if(Ptcle_k[z]<25)
						{
							clogcheck++;
						}
					}
					if(clogcheck>0)
					{
						clognum++;
					}
					//printf("%d ",clognum);
				}
				moves++;
				/*
				//checks to see if the filter is clogged
				for(k=26;k<=50;k++)//j=1-25,k=26-50 is the porous filter
				{
					clogcheck=0;
					for(j=1;j<=25;j++)
					{
						if(B[j][k]>0)
						{
							clogcheck++;
						}
					}
					if(clogcheck>=25)
					{
						repeat=0;
					}
				}//if the filter is clogged, repeat=0, and moving particle loop stops.
				//if filter is not clogged, repeat>0 unless there is no particles to enter the filter
				*/

                if(failure==0)
				{
					inte++;
					percentclog=clognum;
				}
			}
		}
		if(inte==0)
		{
			percentclog=100;
		}
		printf("failure=%d ", failure);
		printf("moves=%d ", moves-1);
		printf("porocity=%d percentclog=%d\n", porocity, percentclog);
		fprintf(fp, "porocity=%d percentclog=%d\n", porocity, percentclog);
		//printf("failure=%d", failure);
/*		if(clognum==0)
		{
			percentclog=100;
			moves=1;
		}*/

	}
	printf("p=%d\n",p);//total number of particles

	fclose(fp);
        
    //e=pow(2.718281828459,);
}
#define FALSE 0
#define TRUE 1
/* Globals */
double u[97],c,cd,cm;
int i97,j97;
int test = FALSE;

/*
   This is the initialization routine for the random number generator.
   NOTE: The seed variables can have values between:    0 <= IJ <=
 31328
                                                        0 <= KL <=
 30081
   The random number sequences created by these two seeds are of
 sufficient
   length to complete an entire calculation with. For example, if
 sveral
   different groups are working on different parts of the same
 calculation,
   each group could be assigned its own IJ seed. This would leave each
 group
   with 30000 choices for the second seed. That is to say, this random
   number generator can create 900 million different subsequences --
 with
   each subsequence having a length of approximately 10^30.
*/
void RandomInitialise(int ij,int kl)
{
   double s,t;
   int ii,i,j,k,l,jj,m;

   /*
      Handle the seed range errors
         First random number seed must be between 0 and 31328
         Second seed must have a value between 0 and 30081
   */
   if (ij < 0 || ij > 31328 || kl < 0 || kl > 30081) {
		ij = 1802;
		kl = 9373;
   }

   i = (ij / 177) % 177 + 2;
   j = (ij % 177)       + 2;
   k = (kl / 169) % 178 + 1;
   l = (kl % 169);

   for (ii=0; ii<97; ii++) {
      s = 0.0;
      t = 0.5;
      for (jj=0; jj<24; jj++) {
         m = (((i * j) % 179) * k) % 179;
         i = j;
         j = k;
         k = m;
         l = (53 * l + 1) % 169;
         if (((l * m % 64)) >= 32)
            s += t;
         t *= 0.5;
      }
      u[ii] = s;
   }

   c    = 362436.0 / 16777216.0;
   cd   = 7654321.0 / 16777216.0;
   cm   = 16777213.0 / 16777216.0;
   i97  = 97;
   j97  = 33;
   test = TRUE;
}
double RandomUniform(void)
{
   double uni;

   /* Make sure the initialisation routine has been called */
   if (!test) 
   	RandomInitialise(1802,9373);

   uni = u[i97-1] - u[j97-1];
   if (uni <= 0.0)
      uni++;
   u[i97-1] = uni;
   i97--;
   if (i97 == 0)
      i97 = 97;
   j97--;
   if (j97 == 0)
      j97 = 97;
   c -= cd;
   if (c < 0.0)
      c += cm;
   uni -= c;
   if (uni < 0.0)
      uni++;

   return(uni);
}

/*
  ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
  THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
  VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
  The function returns a normally distributed pseudo-random number
  with a given mean and standard devaiation.  Calls are made to a
  function subprogram which must return independent random
  numbers uniform in the interval (0,1).
  The algorithm uses the ratio of uniforms method of A.J. Kinderman
  and J.F. Monahan augmented with quadratic bounding curves.
*/
double RandomGaussian(double mean,double stddev)
{
   double  q,u,v,x,y;

	/*  
		Generate P = (u,v) uniform in rect. enclosing acceptance region 
      Make sure that any random numbers <= 0 are rejected, since
      gaussian() requires uniforms > 0, but RandomUniform() delivers >=
 0.
	*/
   do {
      u = RandomUniform();
      v = RandomUniform();
   	if (u <= 0.0 || v <= 0.0) {
       	u = 1.0;
       	v = 1.0;
   	}
      v = 1.7156 * (v - 0.5);

      /*  Evaluate the quadratic form */
      x = u - 0.449871;
   	y = fabs(v) + 0.386595;
      q = x * x + y * (0.19600 * y - 0.25472 * x);

      /* Accept P if inside inner ellipse */
      if (q < 0.27597)
			break;

      /*  Reject P if outside outer ellipse, or outside acceptance
 region */
    } while ((q > 0.27846) || (v * v > -4.0 * log(u) * u * u));

    /*  Return ratio of P's coordinates as the normal deviate */
    return (mean + stddev * v / u);
}

/*
   Return random integer within a range, lower -> upper INCLUSIVE
*/
int RandomInt(int lower,int upper)
{
   return((int)(RandomUniform() * (upper - lower + 1)) + lower);
}

/*
   Return random float within a range, lower -> upper
*/
double RandomDouble(double lower,double upper)
{
   return((upper - lower) * RandomUniform() + lower);
}