#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <conio.h>
#include <ctype.h>
#include <sys/types.h>

#define max 95
#define s_num 198*3

int main ()
{
static const char filename[] = "data.txt";
FILE *file = fopen ( filename, "r" );
FILE *f1 = fopen("output.txt","wt");
int i, j,*x,*y,count=0,gap_count,k=0,l=0;
char b_gap[s_num][max];

//emission matices
int **emission;
double **e_matrix;

//transition matrices
double **transition;

//amino acids
char amino_acids[4]={'A','T','C','G'};

int amino_count[48][4];

int seq[max];

char line[max]; /* or other suitable maximum line size */ 
int len=0,length;

//initializing arrays with null values
for(i=0;i<s_num;i++)
{
	for(j=-1;j<max;j++)
	{
		b_gap[i][j]='\0';
	}
}

//initialization to x&&y
x=(int *) malloc(6*sizeof(int));
y=(int *) malloc(6*sizeof(int));

//Calculating emission matrices
//initialization && delete states
for(i=0;i<6;i++)
{
	y[i]=0;
}

for(i=0;i<6;i++)
{
	x[i]=0;
}
//initializing character array with null values
for(i=0;i<max;i++)
{
	line[i] = '\0';
	seq[i] = 0;
}

emission=(int **) malloc(48*sizeof(int*));

//Calculating emission matrices
//initialization && delete states
for(i=0;i<48;i++)
{
	emission[i] = (int* )malloc( 4*sizeof(int) );
	for(j=0;j<4;j++)
	{
		emission[i][j]=0;
	}
}

//checking the file 
if ( file != NULL )
{ 
	j=0;
	while ( fgets ( line, sizeof line, file ) != NULL )
	{
		if (line[0] != '\n')
		{

			len=strlen(line);
			if((j%2)==0) 
				{
					seq[k]=len;
					k++;
				}

			//printing all initial sequences
			for(i=0;i<len;i++)
			{
				b_gap[j][i]=line[i];
			//	printf("%c",b_gap[j][i]);
			}
			//	printf("\n"); 
			j++;
		}
	}
}

else
{
perror ( filename ); /* why didn't the file open? */
}

//counting frequencies of A/T/C/G
for(i=0;i<198*2;i=i+2)
{
	count=0;
	for(j=0;j<48;j++)
	{
	for(l=count;l<count+1;l++)
	{
		if(((b_gap[i+1][l-1]=='_') && (b_gap[i+1][l]=='_'))||
		   ((b_gap[i+1][l-1]=='d') && (b_gap[i+1][l]=='d'))||
		   ((b_gap[i+1][l-1]=='a') && (b_gap[i+1][l]=='a'))||
		   ((b_gap[i+1][l-1]=='v') && (b_gap[i+1][l]=='v'))||
		   ((b_gap[i+1][l-1]=='t') && (b_gap[i+1][l]=='t')))
			{
				for(k=0;k<4;k++)
				{
					if(b_gap[i][l]==amino_acids[k])
					{
						emission[j-1][k]= (emission[j-1][k])+1;						
						j--;
					}
				}					
			}
			else 
				{					
					for(k=0;k<4;k++)
				{
					if(b_gap[i][l]==amino_acids[k])
					{
						emission[j][k]= (emission[j][k])+1;										
					}
				}
			    }		
	}
	count=count+1;
	}
}

//printf("\n%d\n",j);
printf("\n");
fprintf(f1,"Emission Matrix:\n");
fprintf(f1,"----------------\n");

//e_matrix initialization
e_matrix=(double **) malloc(48*sizeof(double*));
for(i=0;i<48;i++)
{
	e_matrix[i] = (double* )malloc( 4*sizeof(double) );
	for(j=0;j<4;j++)
	{
		e_matrix[i][j]=0.0;
	}
}

for(i=0;i<48;i++)
{
	for(j=0;j<4;j++)
	{
		fprintf(f1,"[%d][%d]: %0.3lf\t",i,j,(double)(emission[i][j]/198.00));
	}
	fprintf(f1,"\n");
}

//initializing transition matrices
transition=(double **) malloc(48*sizeof(double*));

//initializing 0s to transition matrices
for(i=0;i<48;i++)
{
	transition[i] = (double* )malloc(48*sizeof(double) );
	for(j=0;j<48;j++)
	{
		transition[i][j]=0.0;
	}
}

for(i=0;i<47;i++)
{
	transition[i][i+1]=1.0;
}

//counting transitions of A/T/C/G
for(i=0;i<198*2;i=i+2)
{
	count=0;
	for(j=0;j<48;j++)
	{
	for(l=count;l<count+1;l++)
	{
		if ((b_gap[i+1][l]=='_') && (j==7))
				{
					x[0]=x[0]+1;						
					if (b_gap[i+1][l+1]=='_') j=j-1;
				}
		if ((b_gap[i+1][l]=='d'))
				{
					x[1]=x[1]+1;
					if (b_gap[i+1][l+1]=='d') j=j-1;
		}
		if ((b_gap[i+1][l]=='_') && (j==17))
		{
			x[2]=x[2]+1;
			if (b_gap[i+1][l+1]=='_') j=j-1;
		}
		if ((b_gap[i+1][l]=='a'))
		{
			x[3]=x[3]+1;
			if (b_gap[i+1][l+1]=='a') j=j-1;
		}
		if ((b_gap[i+1][l]=='v'))
			{
				x[4]=x[4]+1;
				if (b_gap[i+1][l+1]=='v') j=j-1;
		   }
		if ((b_gap[i+1][l]=='t'))
				{
					x[5]=x[5]+1;
					if (b_gap[i+1][l+1]=='t') j=j-1;
		  }		
	}
	count=count+1;
	}
}

k=0;
for(i=0;i<48;i++)
{
	if(i==7) 
	{
		transition[i][i]=(double)(x[k]+1.00)/((x[k]+198.00)+2.00);		
		transition[i][i+1]=(double)(198.00+1.00)/((x[k]+198.00)+2.00);;		
		k++;
	}
	
	if(i==12) 
	{
		transition[i][i]=(double)(x[k]+1.00)/((x[k]+198.00)+2.00);;	
		transition[i][i+1]=(double)(198.00+1.00)/((x[k]+198.00)+2.00);;
		k++;
	}

	if(i==17) 
	{
		transition[i][i]=(double)(x[k]+1.00)/((x[k]+198.00)+2.00);;
		transition[i][i+1]=(double)(198.00+1.00)/((x[k]+198.00)+2.00);;
		k++;
	}

	if(i==23) 
	{
		transition[i][i]=(double)(x[k]+1.00)/((x[k]+198.00)+2.00);;
		transition[i][i+1]=(double)(198.00+1.00)/((x[k]+198.00)+2.00);;
		k++;
	}

	if(i==29) 
	{
		transition[i][i]=(double)(x[k]+1.00)/((x[k]+198.00)+2.00);;
		transition[i][i+1]=(double)(198.00+1.00)/((x[k]+198.00)+2.00);;
		k++;
	}

	if(i==35) 
	{
		transition[i][i]=(double)(x[k]+1.00)/((x[k]+198.00)+2.00);;
		transition[i][i+1]=(double)(198.00+1.00)/((x[k]+198.00)+2.00);;
		k++;
	}
}

fprintf(f1,"\nTransition Matrix:\n");
fprintf(f1,"------------------\n");

for(i=0;i<48;i++)
{
	for(j=0;j<48;j++)
	{
			fprintf(f1,"[%d][%d]:%0.3lf\t",i,j,transition[i][j]);
	}
	fprintf(f1,"\n");
}
//free memory
/*for(i = 0; i < 48; i++)
	 {
          free(emission[i]);
     }
free(emission);
free(transition);*/
fclose (file);
fclose (f1);
getch();
return 0;
}