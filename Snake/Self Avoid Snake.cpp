#include "stdio.h"
#include "stdlib.h"
#include "math.h"
void   RandomInitialise(int,int);
double RandomUniform(void);
double RandomGaussian(double,double);
int    RandomInt(int,int);
double RandomDouble(double,double);

int plus ( int num1 );
int minus ( int num2 );


int main(int argc,char **argv)
{

	int i, j, k, z, repeat, occ, endpt, times, testvar;
	double r,rmin=1e32,rmax=-1e32;
	int B[40][40], x0_place[1600], y0_place[1600], x_place[1600], y_place[1600];
	FILE *fp;

	fp=fopen("Self Avoid Snake 01.txt", "w+");

	/* Generate 20000 random numbers */
	RandomInitialise(1802,9373);
	for (i=0;i<20000;i++) 
	{
		r = RandomUniform();
		if (r < rmin) rmin = r;
		if (r > rmax) rmax = r;
	}
	fprintf(stderr,"Numbers range from %g to %g\n",rmin,rmax);

	//pre-sets everything to 0
	for(j=0;j<=39;j++)
	{
		for(k=0;k<=39;k++)
		{
			B[j][k]=0;
		}
	}
	
	//starting point will be (arbitrarily) at B[20][20].
	//chain will be 16 units long
	B[20][20]=1;
	j=20;
	k=20;
	y0_place[0]=j;
	x0_place[0]=k;

	for(i=1;i<=15;i++)
	{
		r=(int)(RandomUniform()*4);

		if(r==0)
		{
			j=((j+1) + (5*40))%40;
			B[j][k]=B[j][k]++;
			y0_place[i]=j;
			x0_place[i]=k;
		}
		if(r==1)
		{
			j=((j-1) + (5*40))%40;
			B[j][k]=B[j][k]++;
			y0_place[i]=j;
			x0_place[i]=k;
		}
		if(r==2)
		{
			k=((k+1) + (5*40))%40;
			B[j][k]=B[j][k]++;
			y0_place[i]=j;
			x0_place[i]=k;
		}
		if(r==3)
		{
			k=((k-1) + (5*40))%40;
			B[j][k]=B[j][k]++;
			y0_place[i]=j;
			x0_place[i]=k;
		}
	}

	//The deleting and adding  part:
	repeat=1;
	times=1;
	while(repeat==1)
	{
		endpt=(int)(RandomUniform()*2);

		printf("times=%d\n",times);times++;
		if(endpt==0)//deletes point at the end -> adds to beginning
		{

			for(z=1; z<=15; z++)
			{
				y_place[z]=y0_place[z-1];
				x_place[z]=x0_place[z-1];
			}
			//removes end point
			y0_place[15]--;
			x0_place[15]--;

			//determines which spaces are occupied at beginning point
			if(B[y0_place[0]][minus(x0_place[0])]!=0)
			{
				occ=0;
			}
			if(B[minus(y0_place[0])][x0_place[0]]!=0)
			{
				occ=1;
			}
			if(B[y0_place[0]][plus(x0_place[0])]!=0)
			{
				occ=2;
			}
			if(B[plus(y0_place[0])][x0_place[0]]!=0)
			{
				occ=3;
			}

			while(repeat==1)
			{
				r=(int)(RandomUniform()*4);

				if(r==0 && occ!=0)
				{
					y_place[0]=y0_place[0];
					x_place[0]=minus(x0_place[0]);
					repeat=0;
				}
				if(r==1 && occ!=1)
				{
					y_place[0]=minus(y0_place[0]);
					x_place[0]=x0_place[0];
					repeat=0;
				}
				if(r==2 && occ!=2)
				{
					y_place[0]=y0_place[0];
					x_place[0]=plus(x0_place[0]);
					repeat=0;
				}
				if(r==3 && occ!=3)
				{
					y_place[0]=plus(y0_place[0]);
					x_place[0]=x0_place[0];
					repeat=0;
				}
			}

		}
		if(endpt==1)//deletes point at the beginning -> adds to end
		{
			for(z=0; z<=14; z++)
			{
				y_place[z]=y0_place[z+1];
				x_place[z]=x0_place[z+1];
			}
			//removes initial point
			y0_place[0]--;
			x0_place[0]--;

			//determines which spaces are occupied at last point
			if(B[y0_place[15]][minus(x0_place[15])]!=0)
			{
				occ=0;
			}
			if(B[minus(y0_place[15])][x0_place[15]]!=0)
			{
				occ=1;
			}
			if(B[y0_place[15]][plus(x0_place[15])]!=0)
			{
				occ=2;
			}
			if(B[plus(y0_place[15])][x0_place[15]]!=0)
			{
				occ=3;
			}

			while(repeat==1)
			{
				r=(int)(RandomUniform()*4);

				if(r==0 && occ!=0)
				{
					y_place[15]=y0_place[15];
					x_place[15]=minus(x0_place[15]);
					repeat=0;
				}
				if(r==1 && occ!=1)
				{
					y_place[15]=minus(y0_place[15]);
					x_place[15]=x0_place[15];
					repeat=0;
				}
				if(r==2 && occ!=2)
				{
					y_place[15]=y0_place[15];
					x_place[15]=plus(x0_place[15]);
					repeat=0;
				}
				if(r==3 && occ!=3)
				{
					y_place[15]=plus(y0_place[15]);
					x_place[15]=x0_place[15];
					repeat=0;
				}
			}
		}

		for(j=0;j<=39;j++)
		{
			for(k=0;k<=39;k++)
			{
				if(B[j][k]>1)
				{
					repeat=1;
				}
				else
				{
					repeat=0;
				}
			}
		}
	}
	for(i=0;i<=15;i++)
	{
		fprintf(fp, "%d %d\n", x_place[i], y_place[i]);
	}
	fclose(fp);
}

int plus ( int num1 )
{
	int movement1;
	movement1=((num1+1) + (5*40))%40;

	return movement1;
}

int minus ( int num2 )
{
	int movement2;
	movement2=((num2-1) + (5*40))%40;

	return movement2;
}

int minus1( int num2 );
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