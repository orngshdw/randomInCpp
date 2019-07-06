#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "randomlib.h"

int main(int argc,char **argv)
{
	FILE *fp;
   long int i, k, ident, j=0, med[100][100], col, row, perc, num, clus[10000], clusnum;
   int dum1, dum2;//dummies
   double r,rmin=1e32,rmax=-1e32;

   /* Generate 20000 random numbers */
   RandomInitialise(1802,9373);
   for (i=0;i<20000;i++) {
      r = RandomUniform();
      if (r < rmin) rmin = r;
      if (r > rmax) rmax = r;
   }
   fprintf(stderr,"Numbers range from %g to %g\n",rmin,rmax);

   fp = fopen("Number of clusters per percent pores.txt", "w+");
   fprintf(fp, "percentpores cluster\n");

   for(perc=1; perc<=99; perc++)
   {

   //perc=90;
   //pre-sets all the values in the medium/grid to be 3000
       for(row=0; row<=51; row++)
       {
           for(col=0; col<=51; col++)
		   {
               med[row][col]=99999;
           }
       }

//////////creates the pores!
	   num=1;//so that each new grid will have the 1st pore labelled 1
	   for(row=1; row<=50; row++)
       {
           for(col=1; col<=50; col++)
		   {
			   r=(int)(RandomUniform()*100);

               if(r<=perc)
			   {
				   med[row][col]=num;
				   num++;
			   }
			   else
			   {
				   med[row][col]=99999;
			   }
           }
       }


//checks for clusters from left to right starting at the top, cross

	   for(row=1; row<=50; row++)
       {
           for(col=1; col<=50; col++)
		   {
			   
			   if(med[row][col]<3000)
			   {
				   			 
				   if(med[row+1][col]<3000)
				   {
					   med[row+1][col]=med[row][col];
					   if(med[row+1][col-1]<3000)
					   {
						   med[row+1][col-1]=med[row][col];
					   }
					   if(med[row+1][col+1]<3000)
					   {
						   med[row+1][col+1]=med[row][col];
					   }
					   if(med[row+2][col]<3000)
					   {
						   med[row+2][col]=med[row][col];
					   }
				
				   }
				   			  
				   if(med[row-1][col]<3000)
				   {
					   med[row-1][col]=med[row][col];

					   if(med[row-1][col-1]<3000)
					   {
						   med[row-1][col-1]=med[row][col];
					   }
					   if(med[row-1][col+1]<3000)
					   {
						   med[row-1][col+1]=med[row][col];
					   }
					   if(med[row-2][col]<3000)
					   {
						   med[row-2][col]=med[row][col];
					   }
					   			   
				   }
				   	
				   if(med[row][col+1]<3000)
				   {
					   med[row][col+1]=med[row][col];

					   if(med[row+1][col+1]<3000)
					   {
						   med[row+1][col+1]=med[row][col];
					   }
					   if(med[row-1][col+1]<3000)
					   {
						   med[row-1][col+1]=med[row][col];
					   }
					   if(med[row][col+2]<3000)
					   {
						   med[row][col+2]=med[row][col];
					   }
					   			   
				   }
				   			  
				   if(med[row][col-1]<3000)
				   {
					   med[row][col-1]=med[row][col];
					   
					   if(med[row+1][col-1]<3000)
					   {
						   med[row+1][col-1]=med[row][col];
					   }
					   if(med[row-1][col-1]<3000)
					   {
						   med[row-1][col-1]=med[row][col];
					   }
					   if(med[row][col-2]<3000)
					   {
						   med[row][col-2]=med[row][col];
					   }
					   			   
				   }

			   }
		   }
	   }


//checks for clusters from top to bottom starting from the left, cross 
	   for(col=1; col<=50; col++)
       {
           for(row=1; row<=50; row++)
		   {
			   
			   if(med[row][col]<3000)
			   {
				   			 
				   if(med[row+1][col]<3000)
				   {
					   med[row+1][col]=med[row][col];
					   if(med[row+1][col-1]<3000)
					   {
						   med[row+1][col-1]=med[row][col];
					   }
					   if(med[row+1][col+1]<3000)
					   {
						   med[row+1][col+1]=med[row][col];
					   }
					   if(med[row+2][col]<3000)
					   {
						   med[row+2][col]=med[row][col];
					   }
				
				   }
				   			  
				   if(med[row-1][col]<3000)
				   {
					   med[row-1][col]=med[row][col];

					   if(med[row-1][col-1]<3000)
					   {
						   med[row-1][col-1]=med[row][col];
					   }
					   if(med[row-1][col+1]<3000)
					   {
						   med[row-1][col+1]=med[row][col];
					   }
					   if(med[row-2][col]<3000)
					   {
						   med[row-2][col]=med[row][col];
					   }
					   			   
				   }
				   	
				   if(med[row][col+1]<3000)
				   {
					   med[row][col+1]=med[row][col];

					   if(med[row+1][col+1]<3000)
					   {
						   med[row+1][col+1]=med[row][col];
					   }
					   if(med[row-1][col+1]<3000)
					   {
						   med[row-1][col+1]=med[row][col];
					   }
					   if(med[row][col+2]<3000)
					   {
						   med[row][col+2]=med[row][col];
					   }
					   			   
				   }
				   			  
				   if(med[row][col-1]<3000)
				   {
					   med[row][col-1]=med[row][col];
					   
					   if(med[row+1][col-1]<3000)
					   {
						   med[row+1][col-1]=med[row][col];
					   }
					   if(med[row-1][col-1]<3000)
					   {
						   med[row-1][col-1]=med[row][col];
					   }
					   if(med[row][col-2]<3000)
					   {
						   med[row][col-2]=med[row][col];
					   }
					   			   
				   }

			   }
		   }
	   }


/////prints the grid for a certain percentage
	   /*if(perc==40)
	   {

		   printf("Grid for %d percentage\n", perc);

           for(row=0; row<=51; row++)
           {
               for(col=0; col<=51; col++)
               {
                   printf("%d ",med[row][col]);
               }

               printf("\n\n");
           }
	   }*/

////counts the number of clusters

	   ident=0;
	   clus[0]=0;
	   for(row=0; row<=51; row++)
       {
           for(col=0; col<=51; col++)
		   {
               if(med[row][col]<3000)
			   {
				   clus[ident]=med[row][col];
				   ident++;
			   }
           }
       }

	   	   if(perc==40)
	   {
           for(i=0;i<=ident-1;i++)
		   {
				printf("%d ", clus[ident-1]);
		   }
	   }


////rearranges the list of cluster numbers
	   for(i=0;i<=ident-1;i++)
	   {
		   for(k=i+1;k<=ident-1;k++)
		   {
			   if(clus[i]>clus[k])
			   {
				   dum1=clus[i];
				   clus[i]=clus[k];
				   clus[k]=dum1;
			   }
		   }
	   }



		   

	   

////will count the number of clusters
	   dum2=2;
	   for(i=0;i<=ident-1;i++)
	   {
		   for(k=i+1;k<=ident-1;k++)
		   {
			   if(clus[i]!=clus[k])
			   {
				   dum2++;
				   i=k;
			   }
		   }
	   }
	   clusnum=dum2-1;

////Writes in the data down
	  // fprintf(fp, "%d %d\n", perc, clusnum);
	}
	fclose(fp);
}

#define FALSE 0
#define TRUE 1

/*
   This Random Number Generator is based on the algorithm in a FORTRAN
   version published by George Marsaglia and Arif Zaman, Florida State
   University; ref.: see original comments below.
   At the fhw (Fachhochschule Wiesbaden, W.Germany), Dept. of Computer
   Science, we have written sources in further languages (C, Modula-2
   Turbo-Pascal(3.0, 5.0), Basic and Ada) to get exactly the same test
   results compared with the original FORTRAN version.
   April 1989
   Karl-L. Noell <NOELL@DWIFH1.BITNET>
      and  Helmut  Weber <WEBER@DWIFH1.BITNET>

   This random number generator originally appeared in "Toward a Universal
   Random Number Generator" by George Marsaglia and Arif Zaman.
   Florida State University Report: FSU-SCRI-87-50 (1987)
   It was later modified by F. James and published in "A Review of Pseudo-
   random Number Generators"
   THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
   (However, a newly discovered technique can yield
   a period of 10^600. But that is still in the development stage.)
   It passes ALL of the tests for random number generators and has a period
   of 2^144, is completely portable (gives bit identical results on all
   machines with at least 24-bit mantissas in the floating point
   representation).
   The algorithm is a combination of a Fibonacci sequence (with lags of 97
   and 33, and operation "subtraction plus one, modulo one") and an
   "arithmetic sequence" (using subtraction).

   Use IJ = 1802 & KL = 9373 to test the random number generator. The
   subroutine RANMAR should be used to generate 20000 random numbers.
   Then display the next six random numbers generated multiplied by 4096*4096
   If the random number generator is working properly, the random numbers
   should be:
           6533892.0  14220222.0  7275067.0
           6172232.0  8354498.0   10633180.0
*/

/* Globals */
double u[97],c,cd,cm;
int i97,j97;
int test = FALSE;

/*
   This is the initialization routine for the random number generator.
   NOTE: The seed variables can have values between:    0 <= IJ <= 31328
                                                        0 <= KL <= 30081
   The random number sequences created by these two seeds are of sufficient
   length to complete an entire calculation with. For example, if sveral
   different groups are working on different parts of the same calculation,
   each group could be assigned its own IJ seed. This would leave each group
   with 30000 choices for the second seed. That is to say, this random
   number generator can create 900 million different subsequences -- with
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

/* 
   This is the random number generator proposed by George Marsaglia in
   Florida State University Report: FSU-SCRI-87-50
*/
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
      gaussian() requires uniforms > 0, but RandomUniform() delivers >= 0.
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

      /*  Reject P if outside outer ellipse, or outside acceptance region */
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

