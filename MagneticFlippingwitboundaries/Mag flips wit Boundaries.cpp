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
   int i, j, k, spin, h, times, mag;
   int B[32][32];
   double e;
   double r,rmin=1e32,rmax=-1e32;
   double log(double root);
   double pow (double x, double y);
   FILE *fp;

   /* Generate 20000 random numbers */
   RandomInitialise(1802,9373);
   for (i=0;i<20000;i++) {
      r = RandomUniform();
      if (r < rmin) rmin = r;
      if (r > rmax) rmax = r;
   }
   fprintf(stderr,"Numbers range from %g to %g\n",rmin,rmax);

   //32x32 array with all values 0 DONE
   
   for(j=0;j<32;j++)
	{
		for(k=0;k<32;k++)
		{
			B[j][k]=0;
		}
	}


	//creates a 32x32 board of random spin values for each array box DONE
	
   for(j=0;j<=31;j++)
	{
		for(k=0;k<=31;k++)
		{
			r=RandomUniform();

			if(r<0.5)
			{
				spin=-1;
			}
			else
			{
				spin=1;
			}
			B[j][k]=spin;

			spin=0;
		}
	}

	fp=fopen("Magnetism Values_5000_boundaries.txt", "w+");
	fprintf(fp, "Magnetism Number of times done\n");
	
	for(times=0;times<=5000;times++)
	{
		mag=0;
        for(i=1;i<=32;i++)//Randomly selects 32 cells
		{
            //picks a random cell in the array. Note: 1156 cells in a 32x32 array
            r=RandomUniform();
			j=(int)(r*32);//max r value ~.99 Need 0-31 as range
			r=RandomUniform();
			k=(int)(r*32);

			B[j][k]=-1*B[j][k];//reverses the flip
			
			//finds delta H
			//h=(B[j][k])*(B[( ((j+1) + (5*32))%32 )][k]+B[j-1][k]+B[j][k+1]+B[j][k-1]);
			
			//finds delta H
			h=(B[j][k])*(B[( ((j+1) + (5*32))%32 )][k]
						+B[( ((j-1) + (5*32))%32 )][k]
						+B[j][( ((k+1) + (5*32))%32 )]
						+B[j][( ((k-1) + (5*32))%32 )]);

			if(h>0)
			{
                e=pow(2.718281828459,-h);

                r=RandomUniform();

				if(r<e)
				{
					B[j][k]=-1*B[j][k];//rejects
				}
			}
            
		}

		for(j=0;j<=31;j++)//finds the resulting magnetism of the whole array DONE
		{
            for(k=0;k<=31;k++)
			{
				mag=mag+B[j][k];
			}
		}

		//printf("t=%d m=%d\n", times, mag);
		fprintf(fp, "%d %d\n", times, mag);
	}

	fclose(fp);

}

#define FALSE 0
#define TRUE 1
/* Globals */
double u[97],c,cd,cm;
int i97,j97;
int test = FALSE;

void RandomInitialise(int ij,int kl)
{
   double s,t;
   int ii,i,j,k,l,jj,m;

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

double RandomGaussian(double mean,double stddev)
{
   double  q,u,v,x,y;

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