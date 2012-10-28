#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <conio.h>
#include <ctype.h>
#include <sys/types.h>
#include <math.h>

#define max 128
#define s_num 4

int main ()
{
static const char filename[] = "homologene.txt";
FILE *file = fopen ( filename, "r" );
FILE *f1 =fopen ("output_2.txt", "wt");
int i, j,count=0,gap_count,k=0,l=0;
char b_gap[s_num][max];
char w_gap[s_num][max];

//emission matices
double **emission;

//transition matrices
double **transition;

//amino acids
char amino_acids[20]={'A', 'Q', 'L', 'S', 
					  'R', 'E',	'K', 'T',
					  'N', 'G', 'M', 'W', 
					  'D', 'H', 'F', 'Y', 
					  'C', 'I',	'P', 'V'};
int amino_count[max][20];
int gap[max];
//background probabilities
float b_prob[20]={8.27,3.93, 9.67,	6.52,
				  5.53, 6.76, 5.85, 5.33,
				  4.05, 7.09, 2.42, 1.08,
				  5.45, 2.27, 3.86, 2.92,
				  1.36, 5.98, 4.69, 6.87};

char line[max]; /* or other suitable maximum line size */ 
int len,length;

//for counting match-match,match-delete,mathc-insertion state
int mi[max],mm[max],md[max];
//for counting transition matrix
int seq[max];
char ch;

//initializing all above
for(i=0;i<max;i++)
{
	mi[i] =0;
	mm[i] =0;
	md[i] =0;
	seq[i]=s_num;
}

//initializing arrays with null values
for(i=0;i<s_num;i++)
{
	for(j=0;j<max;j++)
	{
		b_gap[i][j]='\0';
		w_gap[i][j]='\0';
	}
}

//initializing frequecy array
for(i=0;i<max;i++)
{
	for(j=0;j<20;j++)
	{
		amino_count[i][j]=0;
	}
}
//initializing character array with null values
for(i=0;i<max;i++)
{
	line[i] = '\0';
}

//checking the file 
if ( file != NULL )
{ 
	j=0;
	while ( fgets ( line, sizeof line, file ) != NULL )
	{
		for(i=0;line[i]!=0;i++)//read untill null terminator is reached
		{	
			ch=line[i];
		}

		len=strlen(line);

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
else
{
perror ( filename ); /* why didn't the file open? */
}

printf("\nInitial Length:%d\n\n",len);    

/*Discarding gaps when (number of gaps)>= 50% of the (total number of seqs)*/
for(i=0;i<len;i++)
{
gap_count=0;
	for(j=0;j<s_num;j++)
	{
		if(b_gap[j][i]=='-')
		gap_count++;
	}
	if(gap_count<(0.5*j)) 
	{	

		for(k=0;k<s_num;k++)
		{
			w_gap[k][l]=b_gap[k][i];
		}
		l++;
	}
}

printf("\nLength after discarding gaps:%d\n\n",l); 

emission=(double **) malloc(((l*3)+1)*sizeof(double*));

//Calculating emission matrices
//initialization && delete states
for(i=0;i<(l*3)+1;i++)
{
	emission[i] = (double* )malloc( 20*sizeof(double) );
	for(j=0;j<20;j++)
	{
		emission[i][j]=0.0;
	}
}

//handling insertion states
for(i=1;i<(l*3)+1;i=i+3)
{
	for(j=0;j<20;j++)
	{
		emission[i][j]=(b_prob[j]/100.00);
		//printf("%f",emission[i][j]);
	}
}

	for(j=0;j<20;j++)
	{
		emission[l*3][j]=(b_prob[j]/100.00);
	}

//initializing gap
	for(i=0;i<max;i++)
	{
		gap[i]=0;
	}

//counting the gaps for emission matrix of matching states
//counting frequencies of amino acids
for(i=0;i<l;i++)
{
	for(j=0;j<4;j++)
	{
		for(k=0;k<20;k++)
		{
			if(w_gap[j][i]==amino_acids[k])
				amino_count[i][k]= (amino_count[i][k])+1;
		}
		if(w_gap[j][i]=='-')
		gap[i]=(gap[i])+1;
	}
}

//test
/*
for(i=0;i<2;i++)
{
	for(k=0;k<20;k++)
	{	
		printf("%c",amino_acids[k]);
		printf("%d\t",amino_count[i][k]);
	}
	printf("\n");
}
*/
//counting the matching states[Emission matrix]
k=0;
for(i=2;i<(l*3)+1;i=i+3)
{
	for(j=0;j<20;j++)
	{
		emission[i][j]=(double)(((amino_count[k][j])+1.00)/(s_num-gap[k]+20.00));
	//	printf("[%d][%d]: %lf\t",i,j,emission[i][j]);
	}
	k++;
}

fprintf(f1,"Emission Matrix\n");
fprintf(f1,"------------------\n");

for(i=0;i<(l*3)+1;i++)
{
	for(j=0;j<20;j++)
	{
		
		fprintf(f1,"[%d][%d]: %0.3lf\t",i,j,log(emission[i][j]));
	}
	fprintf(f1,"\n");
}

printf("\n------------------------------------------------------------");

//initializing transition matrices
transition=(double **) malloc(((l*3)+1)*sizeof(double*));

//initializing 0s to transition matrices
for(i=0;i<(l*3)+1;i++)
{
	transition[i] = (double* )malloc(((l*3)+1)*sizeof(double) );
	for(j=0;j<(l*3)+1;j++)
	{
		transition[i][j]=0.0;
	}
}

//Comparing complete sequences for matching states
for(i=0;i<=l;i++)
{
	for(j=0;j<4;j++)
	{
		if(w_gap[j][i]!='-' && w_gap[j][i+1]!='-')
			mm[i]=mm[i]+1;

		if(w_gap[j][i]=='-' && w_gap[j][i+1]!='-')
			md[i]=md[i]+1;

		if(w_gap[j][i]!='-' && w_gap[j][i+1]=='-')
			mi[i]=mi[i]+1;
		
		if(w_gap[j][i]=='-' && w_gap[j][i+1]=='-')
			seq[i]=seq[i]+1;
	}
}

//calculating transition matrix for matching states
j=0;
for(i=2;i<(l*3)-3;i=i+3)
{
	transition[i][i+1]= (double)((md[j]+1.0)/(seq[j]+3.0));
//	printf("%lf\t",transition[i][i+1]);
	transition[i][i+2]= (double)((mi[j]+1.0)/(seq[j]+3.0));
//	printf("%lf\t",transition[i][i+2]);
	transition[i][i+3]= (double)((mm[j]+1.0)/(seq[j]+3.0));
//	printf("%lf\t",transition[i][i+3]);
	j++;
}

//calculating last match state
transition[l*3-1][l*3]=(double)((mi[j]+1.0)/(seq[j]+3.0));

//calculating insertion state
for(i=1;i<(l*3);i=i+3)
	{
		transition[i][i]=0.135;
		transition[i][i+1]=0.865;		
	}

//handling last insertion state
transition[l*3][l*3]=0.135;
//printf("%lf\t",transition[l*3][l*3]);

//calculating delete states
for(i=0;i<(l*3)-3;i=i+3)
	{
		transition[i][i+3]=0.135;
		transition[i][i+5]=0.865;		
	}

fprintf(f1,"Transition Matrix\n");
fprintf(f1,"------------------\n");
//printing transition matrix
for(i=0;i<(l*3)+1;i++)
{
	for(j=0;j<(l*3)+1;j++)
	{
		fprintf(f1,"[%d][%d]: %0.3lf\t",i,j,log(transition[i][j]));
		if (transition[i][j]>0.0) count++;
	}
	fprintf(f1,"\n");			
}
printf("\nNumber of non-zero transitions[without start and end states]: %d\n",count);
printf("\nTotal number of non-zero transitions: %d\n",count+6);	//(For start states: +3 && for end end state: +3)

//free memory
/*     for(i = 0; i < (l*3); i++)
	 {
          free(emission[i]);
     }
free(emission);
free(transition);*/
fclose (file);
fclose (f1);
//getch();
return 0;
}