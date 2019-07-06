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

	int i, j, k, p, occ, repeat, check, line;
	double r,rmin=1e32,rmax=-1e32;
	int B[40][40], x_place[1600], y_place[1600];
	FILE *fp;

	/* Generate 20000 random numbers */
	RandomInitialise(1802,9373);
	for (i=0;i<20000;i++) 
	{
		r = RandomUniform();
		if (r < rmin) rmin = r;
		if (r > rmax) rmax = r;
	}
	fprintf(stderr,"Numbers range from %g to %g\n",rmin,rmax);


	
	//pre-set everything to 0
	for(j=0;j<=39;j++)
	{
		for(k=0;k<=39;k++)
		{
			B[j][k]=0;
		}
	}
	
//starting point will be (arbitrarily) at B[20][20].
	B[20][20]=1;
	j=20;
	k=20;
	y_place[0]=j;
	x_place[0]=k;

	for(i=1;i<=15;i++)//creates the random walk of length 16
	{
		r=(int)(RandomUniform()*4);

		if(r==0)//next point is one coordinate up
		{
			j=plus(j);
			B[j][k]=B[j][k]++;//adding the point on the grid
			
			//updating the coordinates:
			y_place[i]=j;
			x_place[i]=k;
		}
		if(r==1)//next point is one coordinate down
		{
			j=minus(j);
			B[j][k]=B[j][k]++;
			y_place[i]=j;
			x_place[i]=k;
		}
		if(r==2)//next point is one coordinate to the right
		{
			k=plus(k);
			B[j][k]=B[j][k]++;
			y_place[i]=j;
			x_place[i]=k;
		}
		if(r==3)//next point is one coordinate left
		{
			k=minus(k);
			B[j][k]=B[j][k]++;
			y_place[i]=j;
			x_place[i]=k;
		}
	}
	for(j=0;j<=15;j++)
	{
		printf("%d %d  ", y_place[j],x_place[j]);

	}
		printf("\n");
	//the flipping of points in the walk
	repeat=1;
	while(repeat==1)
	{
		//picks a random point on walk
		i=(int)(RandomUniform()*16);
		printf("%d ", i);
		/*j=y_place[i];
		k=x_place[i];

	
		//Moves end pt.
		if( i == 0 || i == 15)
		{
			B[y_place[i]][x_place[i]]--;//subtract a point from selected location
			//checks which boxes are open
			if(B[ j, minus(k) ]!=0)//left side is occupied
			{
				occ=1;
			}
			if(B[ minus(j), k ]!=0)//top side is occupied
			{
				occ=2;
			}
			if(B[ j, plus(k) ]!=0)//right side is occupied
			{
				occ=3;
			}
			if(B[ plus(j),k ]!=0)//bottom side is occupied
			{
				occ=4;
			}			

			p=(int)(RandomUniform()*3);//picks the direction to add the point

			//Randomly moves end point
			if(occ==1)
			{
				if(p==0)
				{

					y_place[i]=plus(j);
					x_place[i]=k;
				}	
				if(p==1)
				{
					y_place[i]=minus(j);
					x_place[i]=k;
				}	
				if(p==2)
				{
					y_place[i]=j;
					x_place[i]=plus(k);
				}	
			}
			if(occ==2)
			{
				if(p==0)
				{
					y_place[i]=plus(j);
					x_place[i]=k;
				}	
				if(p==1)
				{
					y_place[i]=j;
					x_place[i]=minus(k);
				}	
				if(p==2)
				{
					y_place[i]=j;
					x_place[i]=plus(k);
				}	
			}
			if(occ==3)
			{
				if(p==0)
				{
					y_place[i]=plus(j);
					x_place[i]=k;
				}	
				if(p==1)
				{
					y_place[i]=minus(j);
					x_place[i]=k;
				}	
				if(p==2)
				{
					y_place[i]=j;
					x_place[i]=minus(k);
				}
			}
			if(occ==4)
			{
				if(p==0)
				{
					y_place[i]=j;
					x_place[i]=minus(k);
				}	
				if(p==1)
				{
					y_place[i]=minus(j);
					x_place[i]=k;
				}	
				if(p==2)
				{
					y_place[i]=j;
					x_place[i]=plus(k);
				}	
			}
			B[y_place[i]][x_place[i]]++;//add a point from selected location
		}
*/
		//Moving Corners
		if(i!=0 && i!=15)
		{
			line=0;
			if(y_place[i]==y_place[i-1] && y_place[i]==y_place[i+1])//point is part of vertical line
			{
				line++;
		
			}
			if(x_place[i]==x_place[i-1] && x_place[i]==x_place[i+1])//point is part of horizon. line
			{
				line++;
		
			}
			if(line==0)
			{

				if(x_place[i]!=x_place[i-1] && y_place[i]!=y_place[i-1] && B[x_place[i-1]][y_place[i-1]]==0)
				{
					//changing the grid:
					B[y_place[i]][x_place[i]]--;//subtract a point from selected location
                    B[x_place[i-1]][y_place[i-1]]++;//adding a new corner
					
					//update the array that stores the random walk coordinates:
					x_place[i]=x_place[i-1];
					y_place[i]=y_place[i-1];
				}
				if(x_place[i]!=x_place[i+1] && y_place[i]!=y_place[i-1] && B[x_place[i+1]][y_place[i-1]]==0)
				{
					B[y_place[i]][x_place[i]]--;
					B[x_place[i+1]][y_place[i-1]]++;
					x_place[i]=x_place[i+1];
					y_place[i]=y_place[i-1];
				}
				if(x_place[i]!=x_place[i-1] && y_place[i]!=y_place[i+1] && B[x_place[i-1]][y_place[i+1]]==0)
				{
					B[y_place[i]][x_place[i]]--;
					B[x_place[i-1]][y_place[i+1]]++;
					x_place[i]=x_place[i-1];
					y_place[i]=y_place[i+1];
				}
				if(x_place[i]!=x_place[i+1] && y_place[i]!=y_place[i+1] && B[x_place[i+1]][y_place[i+1]]==0)
				{
					B[y_place[i]][x_place[i]]--;
					B[x_place[i+1]][y_place[i+1]]++;
					x_place[i]=x_place[i+1];
					y_place[i]=y_place[i+1];
				}
			}
		}

		
		//checks for self avoiding.
		repeat=1;//ensures the while loop will repeat, will be negated if walk is self avoiding
		check=0;//will remain 0 if walk is self avoiding
		
		for(j=0;j<=39;j++)
		{
			for(k=0;k<=39;k++)
			{
				if(B[j][k]>1)//i.e. walk is not self avoiding
				{
					check=check++;//check is changed
				}
			}
		}
		if(check==0)//negates the repeat
		{
			repeat=0;
		}
	}
	fp=fopen("Self Avoid Flip #4.txt", "w");
	fprintf(fp, "x y\n");
	for(i=0;i<=15;i++)
	{
		fprintf(fp, "%d %d\n", x_place[i], y_place[i]);
	}
	fclose(fp);
	printf("\n");
	for(j=0;j<=39;j++)
	{
		for(k=0;k<=39;k++)
		{
			printf("%d ",B[j][k]);
		}
		printf("\n");
	}
	printf("done\n");
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
