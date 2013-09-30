/**************************************************************************************/
/*    gcc stadio5.c -o stadio5 -I CSparse/Include CSparse/Lib/libcsparse.a -lm      */
/***********************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include "strlib.h"
#include <string.h>
#include <ctype.h>
#include "genlib.h"
#include <math.h>
#include "uthash.h"
//#include "csparse.h"
//#include "csparse.cpp"
//#include "cs_lu.c"
//#include "cs_sqr.c"
#include "CSparse/Include/cs.h"
void dianysma_agnwstwn();
void Bi_CG();
void alg_CG();
void anastrofos();
void mna_tran_C();
void backward_sustitution();
void forward_substitution();
void paragontopoihsh_LU();
void  Sunarthsh_sparse();
void findNonZero();
double* solve_A_Tp( double *p);
int nz=0;
double* solve_Mz(double *M,double *r);
void Chollesky();
int svd=0;
double **table;
double **tranC;
double *tableB;
double *tranB;
double *y;
int *P;
double *ch;
double *x;
char **tableX;
int flag3=0;
char PLtypos[5];

/****Variables for CG method ****/

double *r;
double *z;
double rho_previous;
double itol=1e-3;
double normr;		//norma tou r
double normb;		//norma tou b
double *p;
double *q;
double *M;			//M=diag(A) where A=table
double sumr;
double sumb;
double beta;

/****Variables for Bi_CG method ****/

double **ana_table;
double omega;
double alpha;
double *pp;
double *qp;
double *rp;
double *zp;
double EPS=1e-14;
int flag_CG=0;
int SPARSE=0;		//SPARSE flag
cs *A;
cs *C;
char *sparsefile1;
css *S;
double *MS;		//M=diag(A) where A=table
csn *N;
void Sparse_LU();
void Sparse_Cholesky();
void sparse_alg_CG();
void sparse_Bi_CG();
void plot_void(char linePWL[]);

/****Transient Analysis ****/
double time_step;
double fin_time;
int method=0; 
int tran=0;
void mna_b_tran();
void trapezoidal();
void backeuler();
double exp_tran(double i1,double i2,double td1,double tc1,double td2,double tc2);
double sin_tran(double i1,double ia,double fr,double td,double df,double ph);
void pwl_sunarthsh(char line[]);
double pulse_tran(double i1,double i2,double td,double tr,double pw,double tf,double per);
double *xpre;
double *tranBpre;
double t=0;
double *t_pwl;
double *i_pwl;
int *plot_array;
double pwl_tran();
int count_pwl=0;
int countPlot=0;
int count_plot=0;

struct list{
	char type;
	char *counter;
	int first;
	int second;
	int third;
	int fourth;
	double value1;
	double value2;
	char modelname[30];
	char area[30];
	char input_var[5];
	double i1,i2,td1, tc1,td2, tc2,ia,fr,td, df,ph,tr,tf,pw,per;

	char trans_char;
	struct list *nxt;
	struct list *prv;

};

char *substring(char *string, int position, int length) 
{
	char *pointer;
	int c;

	pointer = malloc(length+1);

	if( pointer == NULL )
	{
		return;
	}

	for( c = 0 ; c < position -1 ; c++ ) 
		string++; 

	for( c = 0 ; c < length ; c++ ) 
	{
		*(pointer+c) = *string;      
		string++;   
	}

	*(pointer+c) = '\0'; 
	return pointer;
}

struct list *root;
int group1=0;
int group2=0;
int n=0;			//number of nodes
int i;


struct my_hashtable {
	int id;                    /* key */
	int node;
	UT_hash_handle hh;         /* makes this structure hashable */
};

struct my_hashtable *new_nodes = NULL;

void add_node(int user_id, int mnode) {
	struct my_hashtable *s;

	s = malloc(sizeof(struct my_hashtable));
	s->id = user_id;
	s->node=mnode;
	HASH_ADD_INT( new_nodes, id, s );  /* id: name of key field */
}

struct my_hashtable *find_user(int user_id) {
	struct my_hashtable *s;

	HASH_FIND_INT( new_nodes, &user_id, s );  /* s: output pointer */
	return s;
}



void print_users() {
	struct my_hashtable *s;

	for(s=new_nodes; s != NULL; s=s->hh.next) {
		printf("user id %d: node %d\n", s->id, s->node);
	}
}

void delete_user(struct my_hashtable *user) {
	HASH_DEL( new_nodes, user);  
}


/****Output File of the Results ****/
void paragontopoihsh_output(){

	FILE *fp2;
	fp2=fopen("paragontopoihsh.txt","w");
	int i,j,k;
	fputs("O PINAKAS MNA META THN PARAGONTOPOIHSH\n\n",fp2);

	for(i=0;i<n-1+group2;i++){
		for(k=0;k<n-1+group2;k++){
			char c[10];

			sprintf(c, "%f", table[i][k]);

			for(j=0;j<strlen(c);j++){
				putc(c[j],fp2);	
			}
			putc('\t',fp2);
		}
		putc('\n',fp2);
	}

	fputs("\n\n\n",fp2);


	//fputs("DIANYSMA y (Ly=b)\n\n",fp2);
	//for(i=0;i<n-1+group2;i++){
	//	char c[10];
	//	sprintf(c, "%f", y[i]);
	//	for(j=0;j<strlen(c);j++){
	//		putc(c[j],fp2);
	//	}
	//	putc('\t',fp2);
	//}
	//fputs("\n\n\n",fp2);
	//fputs("DIANYSMA x (Ux=y || L*x=y)\n\n",fp2);


	fputs("VECTOR x \n\n",fp2);
	for(i=0;i<n-1+group2;i++){
		char c[10];
		sprintf(c, "%f", x[i]);
		for(j=0;j<strlen(c);j++){
			putc(c[j],fp2);	
		}
		putc('\t',fp2);

	}
	putc('\n',fp2);
	fclose(fp2);

}


void read_list(){
	int id=1;
	struct my_hashtable *s;
	struct my_hashtable *stemp;
	struct my_hashtable *swaptemp;	
	int stoixeio;
	int j;
	struct list *l2;
	struct list *l3;
	int flag=0;
	int flag1=0;
	int flag2=0;
	int flag3=0;
	int flag4=0;
	int *komvoi= (int*)malloc(sizeof(int));

	for(l2=root->nxt; l2!=root; l2=l2->nxt){

		if((l2->type=='R')||(l2->type=='C')||(l2->type=='I')){
			group1++;
		}
		if((l2->type=='V')||(l2->type=='L')||(l2->type=='Q')){
			group2++;
		}
		/**** Hashtable ****/

		if(flag==0){
			flag=1;
			if(l2->first!=0){
				add_node(id++,l2->first);
			}
			if(l2->second!=0){
				add_node(id++,l2->second);
			}
			if(l2->third>0){
				add_node(id++,l2->third);
			}
			if(l2->fourth>0){
				add_node(id++,l2->fourth);
			}

		} 
		else{
			for(i=1;i<id;i++){
				s=find_user(i);
				if((s->node==l2->first)||(l2->first==0)){
					flag1=1;

				}
				if((s->node==l2->second)||(l2->second==0)){
					flag2=1;

				}
				if((s->node==l2->third)||(l2->third==0)){
					flag3=1;

				}
				if((s->node==l2->fourth)||(l2->fourth==0)){
					flag4=1;

				}
			}


			if(flag1==0){
				add_node(id++,l2->first);
			}
			if(flag2==0){
				add_node(id++,l2->second);
			}
			if((l2->third>=0)&&(flag3==0)){
				add_node(id++,l2->third);
			}
			if((l2->fourth>=0)&&(flag4==0)){
				add_node(id++,l2->fourth);
			}

		}
		flag1=0;
		flag2=0;
		flag3=0;
		flag4=0;     

	}

	printf("nodes before sorting\n");
	print_users();
	printf("nodes after sorting\n");
	/**** SORTING THE HASH TABLE ****/
	for(i=1;i<id-1;i++){		
		for (j=i+1; j<id; j++){
			stemp=find_user(j);  
			s=find_user(i);
			if (s->node > stemp->node){
				swaptemp->node=s->node;
				delete_user(s);
				add_node(i,stemp->node);
				delete_user(stemp);				
				add_node(j,swaptemp->node);				
			}
		}
	}

	print_users();

	for(i=1;i<id;i++){
		s=find_user(i);
		for(l2=root->nxt; l2!=root; l2=l2->nxt){
			if(s->node==l2->first){
				l2->first=s->id;
			}
			else if(s->node==l2->second){
				l2->second=s->id;
			}
			else if(s->node==l2->third){
				l2->third=s->id;
			}
			else if(s->node==l2->fourth){
				l2->fourth=s->id;
			}
		}
	}

	for(l2=root->prv; l2!=root; l2=l2->prv){
		printf("%s\t %c\t%d\t%d\t%d\t%d\t%f\t%f\t%s\t%s",l2->input_var,l2->type,l2->first,l2->second,l2->third,l2->fourth,l2->value1,l2->value2,l2->modelname,l2->area);
		//i1,i2,td1, tc1,td2, tc2,ia,fr,td, df,ph,tr,tf,pw,per;
		printf("\t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f\t%f \t%f \t%f \t%f\n",l2->i1,l2->i2,l2->td1,l2->tc1,l2->td2,l2->tc2,l2->ia,l2->fr,l2->td,l2->df,l2->ph,l2->tr,l2->tf,l2->pw,l2->per);

	}	
	
	printf("\nThe elements of GROUP1 %d\n The elements of GROUP2 %d\n Number of Nodes %d\n",group1,group2,id);
	n=id;

}

/**** Create the table MNA_B ****/
void mna_b(){

	int i,j;
	tableB=(double*)calloc(n-1+group2,sizeof(double));

	for(i=0;i<n-1+group2;i++){
		tableB[i]=0;}

	struct list *l2;
	int k=0;
	for(l2=root->prv; l2!=root; l2=l2->prv){

		switch ( l2->type ) {
		case 'R':

			break;

		case 'V':
			k++;
			tableB[n-2+k]+=l2->value1;
			break;

		case 'I':
			if((l2->first!=0)&&(l2->second!=0)){

				tableB[l2->first-1]+=-l2->value1;    
				tableB[l2->second-1]+=l2->value1; 
			}
			else if((l2->first==0)&&(l2->second!=0)){

				tableB[l2->second-1]+=l2->value1; 
			}
			else if((l2->first!=0)&&(l2->second==0)){

				tableB[l2->first-1]+=-l2->value1;    

			}
			else{}
			break;

		case 'C':
			break;

		case 'L':
			k++;
			tableB[n-2+k]+=0;
			break;

		case 'W':
			break;

		default:
			break;
		} 

	}
}


/**** Create the table MNA_A ****/
void mna(){
	int i,j;

	table = (double**)calloc(n-1+group2,sizeof(double*));
	for(i=0;i<n-1+group2;i++){
		table[i]= (double*)calloc(n-1+group2,sizeof(double));
	}

	for(i=0;i<n-1+group2;i++){
		for(j=0;j<n-1+group2;j++){
			table[i][j]=0;
		}
	}

	struct list *l2;
	int k=0;
	for(l2=root->prv; l2!=root; l2=l2->prv){

		switch ( l2->type ) {
		case 'R':

			if((l2->value1!=0)&&(l2->first!=0)&&(l2->second!=0)){

				table[l2->first-1][l2->first-1]+=1/l2->value1;
				table[l2->first-1][l2->second-1]+=-1/l2->value1;
				table[l2->second-1][l2->first-1]+=-1/l2->value1;
				table[l2->second-1][l2->second-1]+=1/l2->value1;
			}
			else if((l2->value1!=0)&&(l2->first==0)&&(l2->second!=0)){

				table[l2->second-1][l2->second-1]+=1/l2->value1;  
			}
			else if((l2->value1!=0)&&(l2->first!=0)&&(l2->second==0)){

				table[l2->first-1][l2->first-1]+=1/l2->value1;
			}
			else{}
			break;

		case 'V':
			k++;
			if((l2->first!=0)&&(l2->second!=0)){
				table[l2->first-1][n+k-2]+=1.00;
				table[l2->second-1][n+k-2]+=-1.00;
				table[n+k-2][l2->first-1]+=1.00;
				table[n+k-2][l2->second-1]+=-1.00;
			}
			else if((l2->first==0)&&(l2->second!=0)){

				table[l2->second-1][n+k-2]+=-1.00;
				table[n+k-2][l2->second-1]+=-1.00;
			}
			else if((l2->first!=0)&&(l2->second==0)){
				table[l2->first-1][n+k-2]+=1.00;
				table[n+k-2][l2->first-1]+=1.00;
			}
			else{}
			break;

		case 'I':

			break;

		case 'C':

			break;

		case 'L':
			k++;
			if((l2->first!=0)&&(l2->second!=0)){

				table[l2->first-1][n+k-2]+=1.00;
				table[l2->second-1][n+k-2]+=-1.00;
				table[n+k-2][l2->first-1]+=1.00;
				table[n+k-2][l2->second-1]+=-1.00;
			}
			else if((l2->first==0)&&(l2->second!=0)){

				table[l2->second-1][n+k-2]+=-1.00;
				table[n+k-2][l2->second-1]+=-1.00;
			}
			else if((l2->first!=0)&&(l2->second==0)){

				table[l2->first-1][n+k-2]+=1.00;
				table[n+k-2][l2->first-1]+=1.00;
			}
			else{}
			break;
		case 'W':

			break;

		default:

			break;
		}


	}	

}


void mna_tran_C(){

	int i,j;

	tranC = (double**)calloc(n-1+group2,sizeof(double*));
	for(i=0;i<n-1+group2;i++){
		tranC [i]= (double*)calloc(n-1+group2,sizeof(double));
	}

	for(i=0;i<n-1+group2;i++){
		for(j=0;j<n-1+group2;j++){
			tranC[i][j]=0;
		}
	}

	struct list *l2;
	int k=0;
	for(l2=root->prv; l2!=root; l2=l2->prv){

		switch ( l2->type ) {
		case 'C':

			if((l2->value1!=0)&&(l2->first!=0)&&(l2->second!=0)){

				tranC [l2->first-1][l2->first-1]+=l2->value1;
				tranC [l2->first-1][l2->second-1]+=-l2->value1;
				tranC [l2->second-1][l2->first-1]+=-l2->value1;
				tranC [l2->second-1][l2->second-1]+=l2->value1;
			}
			else if((l2->value1!=0)&&(l2->first==0)&&(l2->second!=0)){

				tranC [l2->second-1][l2->second-1]+=l2->value1;  
			}
			else if((l2->value1!=0)&&(l2->first!=0)&&(l2->second==0)){

				tranC [l2->first-1][l2->first-1]+=l2->value1;
			}
			else{}
			break;

		case 'V':

			break;

		case 'I':

			break;

		case 'R':

			break;

		case 'L':
			k++;

			tranC [n+k-2][n+k-2]+=-l2->value1;

			break;
		case 'W':

			break;

		default:

			break;
		}


	}

}


////UNKNOWN VECTOR/////////
void dianysma_agnwstwn()
{
	int i,j,k;
	struct my_hashtable *s;
	tableX = (char**)calloc(n-1+group2,sizeof(char*));
	char **swap_x;
	swap_x=(char**)calloc(n-1+group2,sizeof(char));
	for(i=0;i<n-1+group2;i++){
		swap_x[i]= (char*)calloc(5,sizeof(char));
	}
	for(i=0;i<n-1+group2;i++){
		tableX[i]= (char*)calloc(5,sizeof(char));
	}
	char c[5];
	for(i=0;i<n-1;i++){
		tableX[i][0]='V';
	}
	
	for(i=0;i<n-1;i++){
		s=find_user(i+1);
		sprintf(c, "%d",s->node);
		for(j=1;j<5;j++){
			tableX[i][j]=c[j-1];
		}
	}

	for(i=n-1;i<n-1+group2;i++){
		tableX[i][0]= 'I';

	}
	for(i=n-1;i<n-1+group2;i++){
		for(j=1;j<5;j++){
			sprintf(c, "%d",i+1);
			tableX[i][j]=c[j-1];

		}


	}
}


FILE *eksodos;
void outputfile()
{
	int i,j,temp1;
	int k=0;
	int l;
	char temp[5];
	struct my_hashtable *s;
	i=0;
	char c[11];
	char time[11];
	if(tran==0){       
		for(j=1;j<strlen(PLtypos);j++){
			temp[i]=PLtypos[j];
			i++;
		}
		temp1=atoi(temp);
		for(i=1;i<n;i++)
		{
			s=find_user(i);
			if(s->node==temp1){				
				break;
			}

		}
		if(SPARSE==0){
			sprintf(c, "%f", x[(s->id)-1]);
			for(j=0;j<strlen(c);j++){
				putc(c[j],eksodos);
			}
		}
		else{
			if(flag_CG==0){
				sprintf(c, "%f", tableB[(s->id)-1]);
				for(j=0;j<strlen(c);j++){
					putc(c[j],eksodos);
				}
			}
			else{
				sprintf(c, "%f", x[(s->id)-1]);
				for(j=0;j<strlen(c);j++){
					putc(c[j],eksodos);
				}
			}
		}
	}
	else{
		for(l=0; l<count_plot;l++){
			fputs("\n",eksodos);
			for(i=1;i<n;i++){
				s=find_user(i);
				if(s->node==plot_array[l]){					
					break;
				}
			}			
			fputs("NODE: V",eksodos);
			sprintf(c,"%d", plot_array[l]);
			for(j=0;j<strlen(c);j++){
				putc(c[j],eksodos);
			}
			fputs("\n",eksodos);

			if(SPARSE==0){
				sprintf(c,"%f", x[(s->id)-1]);
				sprintf(time,"%f",t);
				for(j=0;j<strlen(c);j++){
					putc(time[j],eksodos);
				}
				fputs("\t",eksodos);
				for(j=0;j<strlen(c);j++){
					putc(c[j],eksodos);
				}
			}
			else{
				if(flag_CG==0){
					sprintf(c, "%f", tableB[(s->id)-1]);
					for(j=0;j<strlen(c);j++){
						putc(c[j],eksodos);
					}
				}
				else{
					sprintf(c, "%f", x[(s->id)-1]);
					for(j=0;j<strlen(c);j++){
						putc(c[j],eksodos);
					}
				}
			}
		}
	}
	putc('\n',eksodos);
}


void DC_sunarthsh(char DCtypos[5],double startvalue,double endvalue,double increment){

	struct list *l2;
	double i;
	int j;
	for(l2=root->nxt; l2!=root; l2=l2->nxt){
		if(!strcmp(l2->input_var,DCtypos)){
			break;
		}
	}

	eksodos=fopen("APOTELESMATA_PLOT.txt","w");
	fputs("Ta zhtoumena apotelesmata:\n",eksodos);
	x=(double*)calloc(n-1+group2,sizeof(double));
	/**** x : initial guess x[0]=0 ****/
	for(j=0;j<n-1+group2;j++){
		x[j]=0.0;
	}
	if(SPARSE==1){
		Sunarthsh_sparse();
	}
	for(i=startvalue;i<=endvalue;i=i+increment){
		l2->value1=i;
		if(flag_CG==0){
			if(SPARSE==0){
				if(svd==0){
					mna_b();
					forward_substitution();
					backward_sustitution();
				}
				else{
					mna_b();
					Chollesky();
				}
			}
			else{
				if(svd==0){
					Sparse_LU();
				}
				else{
					Sparse_Cholesky();
				}
			}
		}
		else{
			if(SPARSE==0){
				if(svd==0){
					mna_b();
					Bi_CG();
				}
				else{
					mna_b();
					alg_CG();
				}
			}
			else{
				if(svd==0){
					x=(double*)calloc(n-1+group2,sizeof(double));
					/**** x : initial guess x[0]=0 ****/
					for(j=0;j<n-1+group2;j++){
						x[j]=0;
					}
					sparse_Bi_CG();

				}
				else{
					x=(double*)calloc(n-1+group2,sizeof(double));
					/**** x   x : initial guess x[0]=0 ****/
					for(j=0;j<n-1+group2;j++){
						x[j]=0;
					}
					sparse_alg_CG();
				}
			}
		}
		if(flag3==1){
			outputfile();
		}
	}
	fclose(eksodos);
}


int main ( void ){

	root=(struct list*)malloc(sizeof(struct list));
	root->nxt=root;
	root->prv=root;
	char *str1;
	int flag2=0;

	t_pwl=(double*)calloc(50,sizeof(double));
	i_pwl=(double*)calloc(50,sizeof(double));
	plot_array=(int*)calloc(50,sizeof(int));
	char DCtypos[5];

	double startvalue;
	double endvalue;
	double increment;
	static const char filename[] = "netlist.txt";
	FILE *file = fopen ( filename, "r" );
	if ( file != NULL ){	
		char line [128]; /* or other suitable maximum line size */
		char linePWL [128];
		while ( fgets ( line, sizeof line, file ) != NULL ) /* read a line */{
			strcpy(linePWL,line);			
			int flag=0;
			struct list *l;
			l=(struct list*)malloc(sizeof(struct list));

			//initialization -1
			l->third=-1;
			l->fourth=-1;			
			str1 = strtok(line, " ");
			char c=*str1;

			if((c!='*')&&(c!='.')){
				//non-case sensitive
				c=toupper(c);
				l->type=c;

				for(i=0; i<strlen(str1); i++){l->input_var[i]=str1[i];}
				l->input_var[i]='\0';
			}
			int count=1;
			int i,j;
			while (1){
				if((c=='*')){
					free(l);
					flag=1;
					break;
				}
				else if(c=='.'){
					if(strcmp( str1, ".OPTIONS" )==0){						
						while (1){
							str1 = strtok(NULL, " ");							
							if((str1!=NULL)&&(str1[strlen(str1)-1]=='\n')){	
								str1=substring(str1, 0, strlen(str1)-2);
							}							
							if (str1 == NULL){  
								break;
							}
							else if(strcmp( str1, "SPARSE" )==0){								
								SPARSE=1;  
							}
							else if(strcmp( str1, "SPD" )==0){								
								svd=1;  
							}
							else if(strcmp(str1, "ITER" )==0){								
								flag_CG=1;  
							}
							else if(strcmp( str1, "ITOL" )==0){
								str1 = strtok(NULL, " ");
								if (str1 == NULL){
									break;
								}
								else{
									itol=atof(str1);

								}
							}
							else if(!strcmp(substring(str1, 0, strlen(str1)-2),"METHOD=")){								
								if(!strcmp(substring(str1, strlen(str1)-1, strlen(str1)),"BE")){
									method=1;									
								}
							}
						}
					}
					else if(!strcmp( str1, ".DC" )){
						while (1){
							str1 = strtok(NULL, " ");

							if (str1 == NULL){  

								break;
							}
							else{
								if(count==1){

									for(i=0; i<strlen(str1); i++){DCtypos[i]=str1[i];}
									DCtypos[i]='\0';

								}
								else if (count==2){
									startvalue=atof(str1);

								}
								else if (count==3){
									endvalue=atof(str1);

								} 
								else if (count==4){
									increment=atof(str1);

								}
								count++;
							}
						}
						flag2=1;

					}
					else if(!strcmp( str1, ".PLOT" )){
						if(tran==0){							
							while (1){
								str1 = strtok(NULL, " ");

								if (str1 == NULL){ 

									break;
								}
								else{
									if(count==1){										
										j=0;
										for(i=0; i<strlen(str1); i++){
											if((str1[i]=='(')||(str1[i]==')')){}
											else{PLtypos[j]=str1[i];j++;}

										}
										PLtypos[i]='\0';
										count++;
										flag3=1;
									}
								}
							}
						}
						else{
							plot_void(linePWL);

						}

					}
					else if(!strcmp( str1, ".TRAN" )){						
						tran=1;
						while (1){
							str1 = strtok(NULL, " ");
							if (str1 == NULL){ 
								break;
							}
							else{
								if(count==1){

									time_step=atof(str1);									
									count++;
								}
								else if (count==2){
									fin_time=atof(str1);									
								}
							}
						}
					}
					free(l);
					flag=1;
					break;

				}
				/* extract string from string sequence */
				str1 = strtok(NULL, " ");

				/* check if there is nothing else to extract */
				if (str1 == NULL){
					break;
				}
				
				/* print string after tokenized */
				if((c=='V') || (c=='I') ||  (c=='R') || (c=='C') || (c=='L')){

					if(count==1){
						l->first=atoi(str1);

					}
					else if (count==2){
						l->second=atoi(str1);
					}
					else if (count==3){
						l->value1=atof(str1);
					}
					else if (count==4){

						if(!strcmp(str1,"EXP")){							
							l->trans_char=*str1;
						}
						else if(!strcmp(str1,"PULSE")){							
							l->trans_char=*str1;
						}
						else if(!strcmp(str1,"SIN")){	
							l->trans_char=*str1;
						}
						else if(!strcmp(str1,"PWL")){
							l->trans_char='W';
						}						
					}
					else if (count==5){
						str1=substring(str1, 2, strlen(str1));
						if(l->trans_char=='E'){
							l->i1=atof(str1);							
						}
						else if (l->trans_char=='S'){
							l->i1=atof(str1);
						}
						else if (l->trans_char=='P'){
							l->i1=atof(str1);
						}
						else if (l->trans_char=='W'){
							pwl_sunarthsh(linePWL);
						}
					}
					else if (count==6){
						if(l->trans_char=='E'){
							l->i2=atof(str1);
						}
						else if (l->trans_char=='S'){
							l->ia=atof(str1);
						}
						else if (l->trans_char=='P'){
							l->i2=atof(str1);
						}
					}
					else if (count==7){
						if(l->trans_char=='E'){
							l->td1=atof(str1);
						}
						else if (l->trans_char=='S'){
							l->fr=atof(str1);
						}
						else if (l->trans_char=='P'){
							l->td=atof(str1);
						}
					}
					else if (count==8){
						if(l->trans_char=='E'){
							l->tc1=atof(str1);
						}
						else if (l->trans_char=='S'){
							l->td=atof(str1);
						}
						else if (l->trans_char=='P'){
							l->tr=atof(str1);
						}
					}
					else if (count==9){
						if(l->trans_char=='E'){
							l->td2=atof(str1);
						}
						else if (l->trans_char=='S'){
							l->df=atof(str1);
						}
						else if (l->trans_char=='P'){
							l->tf=atof(str1);
						}
					}
					else if (count==10){
						if(l->trans_char=='E'){
							str1=substring(str1, 0, strlen(str1)-1);
							l->tc2=atof(str1);
						}
						else if (l->trans_char=='S'){
							str1=substring(str1, 0, strlen(str1)-1);
							l->ph=atof(str1);
						}
						else if (l->trans_char=='P'){
							l->pw=atof(str1);
						}
					}
					else if (count==11){
						if (l->trans_char=='P'){
							str1=substring(str1, 0, strlen(str1)-1);
							l->per=atof(str1);
						}
					}
				}

				if(c=='D'){
					if(count==1){
						l->first=atoi(str1);
					}
					else if (count==2){
						l->second=atoi(str1);
					}
					else if(count==3){
						for(i=0; i<strlen(str1); i++){l->modelname[i]=str1[i];}
						l->modelname[i]='\0';
					}
					else if(count==4){
						for(i=0; i<strlen(str1); i++){l->area[i]=str1[i];}
						l->area[i]='\0';
					}
					else{}
				}


				if(c=='M'){

					if(count==1){
						l->first=atoi(str1);

					}
					else if (count==2){
						l->second=atoi(str1);
					}
					else if (count==3){
						l->third=atoi(str1);
					}
					else if (count==4){
						l->fourth=atoi(str1);
					}
					else if(count==5){
						for(i=0; i<strlen(str1); i++){l->modelname[i]=str1[i];}
						l->modelname[i]='\0';
					}
					else if (count==6){
						l->value1=atof(str1);
					}
					else if (count==7){
						l->value2=atof(str1);
					}
					else{}
				}

				if(c=='Q'){
					if(count==1){
						l->first=atoi(str1);

					}
					else if (count==2){
						l->second=atoi(str1);
					}
					else if (count==3){
						l->third=atoi(str1);
					}
					else if(count==4){
						for(i=0; i<strlen(str1); i++){l->modelname[i]=str1[i];}
						l->modelname[i]='\0';
					}
					else if(count==5){
						for(i=0; i<strlen(str1); i++){l->area[i]=str1[i];}
						l->area[i]='\0';
					}
					else{}
				}
				count++;
			}
			if(flag==0){
				l->nxt=root->nxt;
				l->prv=root;
				l->nxt->prv=l;
				l->prv->nxt=l;

			}
		}
		fclose ( file );
		read_list();
		mna();
		if(tran==1){
			mna_tran_C();
			xpre=(double*)calloc(n-1+group2,sizeof(double));
			tranBpre=(double*)calloc(n-1+group2,sizeof(double));
			tableB=(double*)calloc(n-1+group2,sizeof(double));			
			x=(double*)calloc(n-1+group2,sizeof(double));
			/**** x : initial guess x[0]=0 ****/
			int j;
			for(j=0;j<n-1+group2;j++){
				x[j]=0.0;
			}			
			eksodos=fopen("APOTELESMATA_PLOT.txt","w");
			fputs("Results according to options:\n",eksodos);
			for(t=0; t<=fin_time; t=t+time_step){
				mna_b_tran();
				if(method==0){
					for(i=0;i<n-1+group2;i++){
						xpre[i]= x[i];
						tranBpre[i]=tranB[i];
					}
					trapezoidal();
				}
				else{
					for(i=0;i<n-1+group2;i++){
						xpre[i]= x[i];
					}
					backeuler();
				}
				M=(double*)calloc(n-1+group2,sizeof(double));
				if(SPARSE==0){
					if((svd==1)&&(flag_CG==0)){
						printf("CALL Chollesky\n");
						Chollesky();
						dianysma_agnwstwn();
					}
					else if((svd==1)&&(flag_CG==1)){
						printf("CALL CG\n");					
						for(i=0;i<n-1+group2;i++)
						{
							
							if(table[i][i]==0){
								M[i]=1;
							}
							else{
								M[i]=1/table[i][i];		// M=1/diag(A) 
							}
						}
					}
					else if((svd==0)&&(flag_CG==1)){
						printf("CALL BI-CG\n");			
						for(i=0;i<n-1+group2;i++)
						{							
							if(table[i][i]==0){
								M[i]=1;
							}
							else{
								M[i]=1/table[i][i];		//M=1/diag(A)
							}	
						}
						anastrofos(); 

					}
					else{
						printf("CALL LU\n");
						paragontopoihsh_LU();
						dianysma_agnwstwn();
					}

					if(flag2==1){
						printf("CALL dc\n");
						DC_sunarthsh(DCtypos,startvalue,endvalue,increment);
					}
					else if((flag2==0)&&(flag_CG==0)){

						mna_b();
						forward_substitution();
						backward_sustitution();
					}
					else if((flag2==0)&&(flag_CG==1)){

						mna_b();
						if(svd==0){
							x=(double*)calloc(n-1+group2,sizeof(double));
							for(i=0;i<n-1+group2;i++){
								x[i]=0;
							}

							Bi_CG();
						}
						else{
							x=(double*)calloc(n-1+group2,sizeof(double));
							for(i=0;i<n-1+group2;i++){
								x[i]=0;
							}
							alg_CG();
						}
					}
					paragontopoihsh_output();

				}
				else{
					Sunarthsh_sparse();
					if((svd==1)&&(flag_CG==0)){
						printf("CALL Sparse Cholesky\n");
						Sparse_Cholesky();
					}
					else if((svd==1)&&(flag_CG==1)){
						printf("CALL sparse CG\n");						
						x=(double*)calloc(n-1+group2,sizeof(double));						
						for(i=0;i<n-1+group2;i++){
							x[i]=0;
						}
						sparse_alg_CG();
					}
					else if((svd==0)&&(flag_CG==1)){

						printf("CALL sparse Bi CG\n");						
						x=(double*)calloc(n-1+group2,sizeof(double));
						/**** x : initial guess x[0]=0 ****/
						for(i=0;i<n-1+group2;i++){
							x[i]=0;
						}
						sparse_Bi_CG();						

					}
					else{
						printf("CALL LU\n");
						Sparse_LU();
					}

					if(flag2==1){
						printf("CALL dc\n");
						DC_sunarthsh(DCtypos,startvalue,endvalue,increment);
					}
				}
				if(flag3==1){
					outputfile();
				}
			}
			fclose(eksodos);
		}		

		if(tran==0){
			M=(double*)calloc(n-1+group2,sizeof(double));
			if(SPARSE==0){
				if((svd==1)&&(flag_CG==0)){
					printf("CALL Chollesky\n");
					Chollesky();
					dianysma_agnwstwn();
				}
				else if((svd==1)&&(flag_CG==1)){
					printf("CALL CG\n");					
					/**** Compute M=diag(table) for CG ****/
					for(i=0;i<n-1+group2;i++)
					{
						if(table[i][i]==0){
							M[i]=1;
						}
						else{
							M[i]=1/table[i][i];		// M=1/diag(A)
						}
					}

				}
				else if((svd==0)&&(flag_CG==1)){
					printf("CALL BI-CG\n");					
					/**** Compute M=diag(table) for BI-CG ****/

					for(i=0;i<n-1+group2;i++)
					{
						
						if(table[i][i]==0){
							M[i]=1;
						}
						else{
							M[i]=1/table[i][i];		// M=1/diag(A)
						}	
					}
					anastrofos(); 
				}
				else{
					printf("CALL LU\n");
					paragontopoihsh_LU();
					dianysma_agnwstwn();
				}

				if(flag2==1){
					printf("CALL dc\n");
					DC_sunarthsh(DCtypos,startvalue,endvalue,increment);
				}
				else if((flag2==0)&&(flag_CG==0)){
					mna_b();
					forward_substitution();
					backward_sustitution();
				}
				else if((flag2==0)&&(flag_CG==1)){
					mna_b();
					if(svd==0){
						x=(double*)calloc(n-1+group2,sizeof(double));
						/**** x : initial guess x[0]=0 ****/
						for(i=0;i<n-1+group2;i++){
							x[i]=0;
						}
						Bi_CG();
					}
					else{
						x=(double*)calloc(n-1+group2,sizeof(double));
						/**** x : initial guess x[0]=0 ****/
						for(i=0;i<n-1+group2;i++){
							x[i]=0;
						}
						alg_CG();
					}
				}
				/**** OUTPU RESULT ****/
				paragontopoihsh_output();

			}
			else{
				Sunarthsh_sparse();
				if((svd==1)&&(flag_CG==0)){
					printf("CALL Sparse Cholesky\n");
					Sparse_Cholesky();
				}
				else if((svd==1)&&(flag_CG==1)){
					printf("CALL sparse CG\n");					
					x=(double*)calloc(n-1+group2,sizeof(double));
					/****  x : initial guess x[0]=0 ****/
					for(i=0;i<n-1+group2;i++){
						x[i]=0;
					}
					sparse_alg_CG();
				}
				else if((svd==0)&&(flag_CG==1)){
					printf("CALL sparse Bi CG\n");
					
					x=(double*)calloc(n-1+group2,sizeof(double));
					/**** x : initial guess x[0]=0 ****/
					for(i=0;i<n-1+group2;i++){
						x[i]=0;
					}
					sparse_Bi_CG();					

				}
				else{
					printf("CALL LU\n");
					Sparse_LU();
				}

				if(flag2==1){
					printf("CALL dc\n");
					DC_sunarthsh(DCtypos,startvalue,endvalue,increment);
				}
			}
		}
	}
	else{
		perror ( filename );	// why didn't the file open? 
	}
	
	return 0;
}


void findNonZero(){

	struct list *l2;
	for(l2=root->nxt; l2!=root; l2=l2->nxt){
		if (l2->type=='R'){
			if((l2->first!=0)&&(l2->second!=0)){
				nz+=4;
			}
			else{
				nz++;
			}
		}   
		if(l2->type=='V'){
			if((l2->first!=0)&&(l2->second!=0)){
				nz+=4;
			}
			else{
				nz+=2;
			}
		}  
	}
	printf("Non zero %d\n",nz);
} 



void paragontopoihsh_LU()
{

	int k,i,l,m,e,j,tmp2;

	double x;
	double tmp;
	P=(int*)calloc(n-1+group2,sizeof(int));
	for(j=0;j<n-1+group2;j++){
		P[j]=j;
	}
	for(k=0;k<n-1+group2;k++){
		x=fabs(table[k][k]);
		for(i=k;i<n-1+group2;i++){
			if(fabs(table[i][k])>=x){
				x=fabs(table[i][k]);
				m=i;

				tmp2=P[m];
				P[m]=P[k];
				P[k]=tmp2;
				for(l=0;l<n-1+group2;l++){
					tmp=table[k][l];
					table[k][l]=table[m][l];
					table[m][l]=tmp;}
			} 

		}
		for(i=k+1;i<n-1+group2;i++){
			if(table[k][k]==0){	
				table[i][k]=table[i][k];
			}
			else{
				table[i][k]=table[i][k]/table[k][k];
			}
		}
		for(i=k+1;i<n-1+group2;i++){
			for(j=k+1;j<n-1+group2;j++){
				table[i][j]=table[i][j]-table[i][k]*table[k][j];
			}
		} 

	}
	svd=0;
	return ;
}


void Chollesky()
{
	int k,j,i;
	double w;
	double sum2,sum; 

	ch=(double*)calloc(n-1+group2,sizeof(double));
	for(i=0;i<n-1+group2;i++){
		ch[i]=0;
	}

	for(k=0; k<=n-1+group2-1; k++){
		sum=0;
		sum2=0;

		for(j=0; j<=k-1; j++){
			sum=sum+((table[k][j])*(table[k][j]));			
		}
		w=table[k][k]-sum;
		if (w>0){
			ch[k]=sqrt(w);			
		}
		else{			
			printf("\n !!Arnhtiko uporrizo!!\n");
			//exit(0);
		}
		for (i=k+1; i<=n-1+group2-1; i++){
			for(j=0; j<=k-1; j++){
				sum2=sum2+(table[k][j]*table[i][j]);
			}
			table[i][k]=(1/ch[k])*(table[i][k]-sum2);
		}
		sum=0;
		sum2=0;
	}	
	printf("\n");
	for(k=0;k<n-1+group2;k++){
		printf("\t%f",ch[k]);
	}
	printf("\n");
	svd=1;
	
}

void forward_substitution()
{
	int k,j,e;
	int grammh,tmp;
	double *swap_b;
	swap_b=(double*)calloc(n-1+group2,sizeof(double));
	y=(double*)calloc(n-1+group2,sizeof(double));
	for(i=0;i<n-1+group2;i++){
		y[i]=0;
	}


	if(svd==0){
		for(k=0;k<n-1+group2;k++){
			swap_b[k]=tableB[P[k]];
		}
		for(k=0;k<n-1+group2;k++){
			tableB[k]=swap_b[k];
		}

	}

	for(k=0;k<n-1+group2;k++){
		for(j=0;j<k;j++){
			tableB[k]=tableB[k]-table[k][j]*y[j];      

		}

		if(svd==0){
			y[k]=tableB[k];
		}
		else{
			y[k]=tableB[k]/ch[k];
		}
	}
	
	return ;

}


void backward_sustitution()
{

	int k,j;
	int grammh,tmp;

	x=(double*)calloc(n-1+group2,sizeof(double));
	for(i=0;i<n-1+group2;i++){
		x[i]=0;
	}

	for(k=n-1+group2-1;k>=0;k--){
		for(j=k+1;j<n-1+group2;j++){
			y[k]=y[k]-table[k][j]*x[j];   

		}
		if(svd==0){
			if(table[k][k]==0){
				x[k]=y[k];
			}
			else{
				x[k]=y[k]/table[k][k];
			}
		}
		else{
			x[k]=y[k]/ch[k];;
		}
	}

	return;
}



void alg_CG(){

	int i,j;
	z=(double*)calloc(n-1+group2,sizeof(double));
	int iter=0; //arithmos epanalhpsewn ths methodou
	r=(double*)calloc(n-1+group2,sizeof(double));//initial residual
	p=(double*)calloc(n-1+group2,sizeof(double));
	q=(double*)calloc(n-1+group2,sizeof(double));
	double rho=0;
	// r=b-Ax
	double sum=0;
	for(i=0;i<n-1+group2;i++){

		// ypologismos A*x
		for(j=0;j<n-1+group2;j++){
			sum=sum+table[i][j]*x[j];
		}

		r[i]=tableB[i]-sum;
		sum=0;
	}


	printf("b\n");
	for(i=0;i<n-1+group2;i++){
		printf("\t%f",tableB[i]);
	}

	printf("r=b-Ax\n");
	for(i=0;i<n-1+group2;i++){
		printf("\t%f",r[i]);
	}	

	sumr=0;
	sumb=0;
	normr=0;	//Norm(r)
	normb=0;	//Norm(b)

	//Norm(r) and Norm(b)
	for(i=0;i<n-1+group2;i++){

		sumr=sumr+r[i]*r[i];
		normr=sqrt(sumr);
		sumb=sumb+tableB[i]*tableB[i];
		normb=sqrt(sumb);	
	}

	if (normb==0)
	{
		normb=1;
	}

	while(((normr/normb) > itol)&&(iter<n-1+group2))
	{
		iter++;
		z=solve_Mz(M,r);
		rho=0;
		for(i=0;i<n-1+group2;i++)
		{
			rho=rho+r[i]*z[i];
		}

		if(iter==1){
			for(i=0;i<n-1+group2;i++){
				p[i]=z[i];
			}
		}
		else{
			beta=rho/rho_previous;
			for(i=0;i<n-1+group2;i++){
				p[i]=z[i]+beta*p[i];
			}
		}

		rho_previous=rho;
		sum=0;
		for(i=0;i<n-1+group2;i++){
			/**** A*p ****/			
			for(j=0;j<n-1+group2;j++){
				sum=sum+table[i][j]*p[j];
			}

			q[i]=sum;
			sum=0;
		}
		/**** p*q ****/		
		for(i=0;i<n-1+group2;i++){

			sum=sum+p[i]*q[i];
		}
		alpha=rho/sum;		
		for(i=0;i<n-1+group2;i++){
			x[i]=x[i]+alpha*p[i];

			r[i]=r[i]-alpha*q[i];
		}

	}
	printf("VECTOR x\n");	
	for(i=0;i<n-1+group2;i++){
		printf("\t%f",x[i]);
	}
	printf("\n");
	//free(z);
	//free(q);
	//free(r);
	//free(p);

}

void anastrofos()
{
	int i,j;
	ana_table = (double**)calloc(n-1+group2,sizeof(double*));
	for(i=0;i<n-1+group2;i++){
		ana_table[i]= (double*)calloc(n-1+group2,sizeof(double));
	}
	for(i=0;i<n-1+group2;i++){
		for(j=0;j<n-1+group2;j++){
			ana_table[i][j]=table[j][i];
		}
	}
}

/**** Compute A*p=q and return p ****/
double* solve_Ap(double **table,double *p)
{
	int i,j;
	double sum;
	double *qqp;
	qqp=(double*)calloc(n-1+group2,sizeof(double));
	sum=0;

	/**** A*p ****/
	for(i=0;i<n-1+group2;i++){

		for(j=0;j<n-1+group2;j++){
			sum=sum+table[i][j]*p[j];
		}

		qqp[i]=sum;
		sum=0;
	}
	return(qqp);	
}

/* Ypologismos toy M*z=r */
double* solve_Mz(double *M,double *r)
{
	double *zzp;
	zzp=(double*)calloc(n-1+group2,sizeof(double));
	int i;

	for(i=0;i<n-1+group2;i++)
	{
		zzp[i]=r[i]*M[i];
	}
	return(zzp);	
}

/******* METHODOS Bi_CG ********/ 
void Bi_CG()
{
	int i;
	/**** r=b-Ax ****/
	double sum=0;
	int j,iter;
	double rho;
	r=(double*)calloc(n-1+group2,sizeof(double));//initial residual
	rp=(double*)calloc(n-1+group2,sizeof(double));
	p=(double*)calloc(n-1+group2,sizeof(double));//initial residual
	pp=(double*)calloc(n-1+group2,sizeof(double));
	q=(double*)calloc(n-1+group2,sizeof(double));//initial residual
	qp=(double*)calloc(n-1+group2,sizeof(double));

	for(i=0;i<n-1+group2;i++){

		/**** A*x ****/
		for(j=0;j<n-1+group2;j++){
			sum=sum+table[i][j]*x[j];
		}

		r[i]=tableB[i]-sum;
		rp[i]=tableB[i]-sum;
		sum=0;
	}

	iter=0;
	sumr=0;
	sumb=0;
	normr=0;	//Norm(r)
	normb=0;	//Norm(b)

	//Norm(r) and Norm(b)
	for(i=0;i<n-1+group2;i++){

		sumr=sumr+r[i]*r[i];
		sumb=sumb+tableB[i]*tableB[i];	
	}

	normr=sqrt(sumr);
	normb=sqrt(sumb);
	if (normb==0){normb=1;}

	while( ((normr/normb) > itol) && (iter<n-1+group2))
	{
		omega=0;
		iter++;
		rho=0;
		z=solve_Mz(M,r);
		zp=solve_Mz(M,rp);

		for(i=0;i<n-1+group2;i++)
		{
			rho=rho+rp[i]*z[i];
		}

		if(fabs(rho)<EPS){
			return;
		}

		if(iter==1){
			for(i=0;i<n-1+group2;i++)
			{
				p[i]=z[i];
				pp[i]=zp[i];
			}
		}
		else{
			beta=rho/rho_previous;

			for(i=0;i<n-1+group2;i++){
				p[i]=z[i]+beta*p[i];
				pp[i]=zp[i]+beta*pp[i];
			}

		}
		rho_previous=rho;

		q=solve_Ap(table,p);
		qp=solve_Ap(ana_table,pp);

		for(i=0;i<n-1+group2;i++){
			omega=omega+pp[i]*q[i];
		}

		if(fabs(omega) <= EPS) { printf(" RETURN ");return; }

		alpha=rho/omega;

		for(i=0;i<n-1+group2;i++){

			x[i]=x[i]+alpha*p[i];
			r[i]=r[i]-alpha*q[i];
			rp[i]=rp[i]-alpha*qp[i];
		}

		//Norm(r) and Norm(b)
		sumr=0;

	}

	printf("VECTOR x");
	for(j=0;j<n-1+group2;j++){
		printf("\t %f",x[j]);
	}
	printf("\n");
}





int cs_printf (const cs *A, int brief){
	int p, j, m, n, nzmax, nz, *Ap, *Ai ;
	double *Ax ;
	if (!A) { printf ("(null)\n") ; return (0) ; }
	m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
	nzmax = A->nzmax ; nz = A->nz ;

	if (nz < 0)
	{
		printf ("%d-by-%d, nzmax: %d nnz: %d, \n", m, n, nzmax,
			Ap [n]) ;

		for (j = 0 ; j < n ; j++)
		{
			printf ("    col %d : locations %d to %d\n", j, Ap [j], Ap [j+1]-1);
			for (p = Ap [j] ; p < Ap [j+1] ; p++)
			{
				printf ("      %d : %g\n", Ai [p], Ax ? Ax [p] : 1) ;
				if (brief && p > 20) { printf ("  ...\n") ; return (1) ; }
			}
		}
	}
	else
	{
		printf ("triplet: %d-by-%d, nzmax: %d nnz: %d\n", m, n, nzmax, nz) ;
		for (p = 0 ; p < nz ; p++)
		{
			printf ("    %d %d : %g\n", Ai [p], Ap [p], Ax ? Ax [p] : 1) ;
			if (brief && p > 20) { printf ("  ...\n") ; return (1) ; }
		}
	}
	return (1) ;
}


void Sunarthsh_sparse(){
	findNonZero();
	A=cs_spalloc(n-1+group2,n-1+group2,nz,1,1);
	struct list *l2;
	int k=0;
	int j=0;
	int *ir;
	int *pc;
	double *xz;

	ir=(int*)calloc(nz,sizeof(int));
	for(j=0;j<nz;j++){
		ir[j]=0;
	} 
	pc=(int*)calloc(nz,sizeof(int));
	for(j=0;j<nz;j++){
		pc[j]=0;
	}
	xz=(double*)calloc(nz,sizeof(double));
	for(j=0;j<nz;j++){
		xz[j]=0.0;
	} 

	j=0;

	for(l2=root->prv; l2!=root; l2=l2->prv){
		switch ( l2->type ) {
		case 'R':
			if((l2->value1!=0)&&(l2->first!=0)&&(l2->second!=0)){
				ir[j]=l2->first-1;
				pc[j]=l2->first-1;
				xz[j]=1/l2->value1;
				j++;
				ir[j]=l2->first-1;
				pc[j]=l2->second-1;
				xz[j]=-1/l2->value1;
				j++;
				ir[j]=l2->second-1;
				pc[j]=l2->first-1;
				xz[j]=-1/l2->value1;
				j++;
				ir[j]=l2->second-1;
				pc[j]=l2->second-1;
				xz[j]=1/l2->value1;
				j++;

			}
			else if((l2->value1!=0)&&(l2->first==0)&&(l2->second!=0)){
				ir[j]=l2->second-1;
				pc[j]=l2->second-1;
				xz[j]=1/l2->value1;  
				j++;
			}
			else if((l2->value1!=0)&&(l2->first!=0)&&(l2->second==0)){
				ir[j]=l2->first-1;
				pc[j]=l2->first-1;
				xz[j]=1/l2->value1;
				j++;
			}
			else{}
			break;

		case 'V':
			k++;
			if((l2->first!=0)&&(l2->second!=0)){
				ir[j]=l2->first-1;
				pc[j]=n+k-2;
				xz[j]=1.00;
				j++;
				ir[j]=l2->second-1;
				pc[j]=n+k-2;
				xz[j]=-1.00;  
				j++;
				ir[j]=n+k-2;
				pc[j]=l2->first-1;
				xz[j]=1.00;
				j++;
				ir[j]=n+k-2;
				pc[j]=l2->second-1;
				xz[j]=-1.00; 
				j++;
			}
			else if((l2->first==0)&&(l2->second!=0)){
				ir[j]=l2->second-1;
				pc[j]=n+k-2;
				xz[j]=-1.00;  
				j++;
				ir[j]=n+k-2;
				pc[j]=l2->second-1;
				xz[j]=-1.00;
				j++;
			}
			else if((l2->first!=0)&&(l2->second==0)){
				ir[j]=l2->first-1;
				pc[j]=n+k-2;
				xz[j]=1.00;
				j++;
				ir[j]=n+k-2;
				pc[j]=l2->first-1;
				xz[j]=1.00;
				j++;
			}
			else{}
			break;

		case 'I':
			break;

		case 'C':
			break;

		case 'L':
			break;

		case 'W':
			break;

		default:
			break;
		}


	}	
	
	for(j=0; j<nz; j++){
		cs_entry(A,ir[j],pc[j],xz[j]);		
	}

	A->nz=nz;
	printf("****AFTER THN ENTRY****\n");
	cs_printf(A,0) ;

	C=cs_compress(A);
	printf("****AFTER THN COMPRESS A****\n");
	cs_printf(A,0) ;

	cs_dupl(C);
	printf("****AFTER DUPL(C)****\n");
	cs_printf(C,0) ;

}

void Sparse_LU()
{
	int j=0;
	int i;

	/**** First step LU ****/
	S=cs_sqr(2,C,0);
	N=cs_lu(C,S,1);
	cs_spfree(C);

	mna_b();	

	/**** Second step LU ****/
	double *x;
	x=cs_malloc(n-1+group2,sizeof(double));
	for(j=0;j<n-1+group2;j++){
		x[j]=0.0;
	} 

	cs_ipvec(N->pinv,tableB,x,n-1+group2);
	cs_lsolve(N->L,x);
	cs_usolve(N->U,x);
	cs_pvec(S->q,x,tableB,n-1+group2);  

	cs_free(x);
	cs_sfree(S);
	cs_nfree(N);

}

void Sparse_Cholesky(){
	int j;

	/**** First step for Cholesky ****/
	S=cs_schol(1,C);
	N=cs_chol(C,S);
	//cs_spfree(C);


	/**** Second step for Cholesky ****/
	mna_b();
	double *x;
	x=cs_malloc(n-1+group2,sizeof(double));
	for(j=0;j<n-1+group2;j++){
		x[j]=0.0;
	} 
	cs_ipvec(S->pinv,tableB,x,n-1+group2);
	cs_lsolve(N->L,x);
	cs_ltsolve(N->L,x);
	cs_pvec(S->pinv,x,tableB,n-1+group2);

	cs_free(x);
	cs_sfree(S);
	cs_nfree(N);	

}


void sparse_alg_CG()
{
	
	int i,j;	
	double *y;	// store the C*x in y
	y=(double*)calloc(n-1+group2,sizeof(double));
	z=(double*)calloc(n-1+group2,sizeof(double));

	int iter=0; // number of iterations of method

	r=(double*)calloc(n-1+group2,sizeof(double));	// initial residual
	p=(double*)calloc(n-1+group2,sizeof(double));
	q=(double*)calloc(n-1+group2,sizeof(double));

	double rho=0;
	double sum=0;

	/**** Compute r=b-Ax where A is C ****/
	for(i=0;i<n-1+group2;i++)
	{
		y[i]=0;
	}

	/**** Compute C*x and store in y ****/
	cs_gaxpy (C,x,y);
	mna_b();
	for(i=0;i<n-1+group2;i++){
		r[i]=tableB[i]-y[i];		
	}	

	sumr=0;
	sumb=0;
	normr=0;		//Norm(r)
	normb=0;		//Norm(b)

	/**** Compute MS (1 / diag(C) )  ****/
	int h;
	MS=(double*)calloc(n-1+group2,sizeof(double));    
	int k=0;
	for(j=0;j<n-1+group2;j++)
	{
		printf("mpainei  sto j=%d kai C->p[j] einai %d\n",j,C->p[j]);
		for(h=C->p[j];h<C->p[j+1];h++){

			if( j==C->i[h]){
				MS[k]=1/C->x[h];
				k++;
			}
		}
	}
	
	//Norm(r) and Norm(b)
	for(i=0;i<n-1+group2;i++){

		sumr=sumr+r[i]*r[i];
		normr=sqrt(sumr);

		sumb=sumb+tableB[i]*tableB[i];
		normb=sqrt(sumb);	
	}

	if (normb==0)
	{
		normb=1;
	}

	while(((normr/normb) > itol)&&(iter<n-1+group2))
	{

		iter++;

		z=solve_Mz(MS,r);

		rho=0;
		for(i=0;i<n-1+group2;i++)
		{

			rho=rho+r[i]*z[i];
		}

		if(iter==1){
			for(i=0;i<n-1+group2;i++){
				p[i]=z[i];
			}
		}
		else{
			beta=rho/rho_previous;

			for(i=0;i<n-1+group2;i++){
				p[i]=z[i]+beta*p[i];
			}
		}

		rho_previous=rho;
		sum=0;

		/**** Compute A*p where A is C ****/


		for(i=0;i<n-1+group2;i++){
			y[i]=0;
		}
		cs_gaxpy (C,p,y);
		for(i=0;i<n-1+group2;i++){
			q[i]=y[i];

		}

		//ypologismos esvterikou ginomenou p*q
		for(i=0;i<n-1+group2;i++){

			sum=sum+p[i]*q[i];
		}

		alpha=rho/sum;

		for(i=0;i<n-1+group2;i++){
			x[i]=x[i]+alpha*p[i];

			r[i]=r[i]-alpha*q[i];
		}

	}
	printf("Vector x");
	printf("\n");
	for(i=0;i<n-1+group2;i++){
		printf("\t%f",x[i]);
	}
	printf("\n");
	//free(z);
	//free(q);
	//free(r);
	//free(p);\\\\\

}

/**** Ypologizei to ginomeno A^T (anastrofos tou A)*x  opou x ena dianysma  ****/
double* solve_A_Tp( double *p){      

	double *y;
	double yi;

	y=(double*)calloc(n-1+group2,sizeof(double));
	int i,j;

	for(i = 0;i <n-1+group2;i++){
		yi = 0.0;
		for(j = C->p[i];j < C->p[i+1];j++){
			yi += C->x[j] * p[C->i[j]]; 

		}
		y[i] = yi; 		
	}

	return(y);
}




void sparse_Bi_CG(){
	int i;
	double *y;
	y=(double*)calloc(n-1+group2,sizeof(double));

	//r=b-Ax
	double sum=0;
	int j,iter;
	double rho;
	r=(double*)calloc(n-1+group2,sizeof(double));//initial residual
	rp=(double*)calloc(n-1+group2,sizeof(double));
	p=(double*)calloc(n-1+group2,sizeof(double));//initial residual
	pp=(double*)calloc(n-1+group2,sizeof(double));
	q=(double*)calloc(n-1+group2,sizeof(double));//initial residual
	qp=(double*)calloc(n-1+group2,sizeof(double));

	for(i=0;i<n-1+group2;i++)
	{
		y[i]=0;
	}

	cs_gaxpy (C,x,y);
	mna_b();

	for(i=0;i<n-1+group2;i++){
		r[i]=tableB[i]-y[i];
		rp[i]=tableB[i]-y[i];

	}

	iter=0;
	sumr=0;
	sumb=0;
	normr=0;		//Norm(r)
	normb=0;		//Norm(b)

	//Norm(r) and Norm(b)
	for(i=0;i<n-1+group2;i++){
		sumr=sumr+r[i]*r[i];
		sumb=sumb+tableB[i]*tableB[i];	
	}

	/**** Compute MS (1 / diag(C) ) ****/
	int h;
	MS=(double*)calloc(n-1+group2,sizeof(double));   
	int k=0;
	for(j=0;j<n-1+group2;j++)
	{

		printf("enter in j=%d kai C->p[j] is %d\n",j,C->p[j]);
		for(h=C->p[j];h<C->p[j+1];h++){

			if( j==C->i[h]){
				MS[k]=1/C->x[h];
				k++;
			}

		}
	}



	normr=sqrt(sumr);
	normb=sqrt(sumb);
	if (normb==0){normb=1;}

	while( ((normr/normb) > itol) && (iter<n-1+group2))
	{
		omega=0;
		iter++;
		rho=0;
		z=solve_Mz(MS,r);
		zp=solve_Mz(MS,rp);

		for(i=0;i<n-1+group2;i++)
		{
			rho=rho+rp[i]*z[i];
		}

		if(fabs(rho)<EPS){
			return;
		}

		if(iter==1){
			for(i=0;i<n-1+group2;i++)
			{
				p[i]=z[i];
				pp[i]=zp[i];
			}
		}
		else{
			beta=rho/rho_previous;

			for(i=0;i<n-1+group2;i++){
				p[i]=z[i]+beta*p[i];
				pp[i]=zp[i]+beta*pp[i];
			}

		}
		rho_previous=rho;

		/**** Solve C*p and save the result in q ****/
		for(i=0;i<n-1+group2;i++)
		{
			q[i]=0;
		}
		cs_gaxpy (C,p,q);
		qp=solve_A_Tp(pp);

		for(i=0;i<n-1+group2;i++){
			omega=omega+pp[i]*q[i];
		}

		if(fabs(omega) <= EPS) { printf(" RETURN ");return; }


		alpha=rho/omega;

		for(i=0;i<n-1+group2;i++){

			x[i]=x[i]+alpha*p[i];
			r[i]=r[i]-alpha*q[i];
			rp[i]=rp[i]-alpha*qp[i];
		}

		//Norm(r) and Norm(b)
		sumr=0;

	}

	printf("Vector x\n");
	for(j=0;j<n-1+group2;j++){
		printf("%f\t",x[j]);
	}
	printf("\n");
}




double exp_tran(double i1,double i2,double td1,double tc1,double td2,double tc2){


	if((t>=0)&&(t<=td1)){
		return i1;
	}
	else if((t>=td1)&&(t<=td2)){
		return i1+(i2-i1)*(1-exp(-(t-td1)/tc1));
	}
	else{
		return  i1+(i2-i1)*((1-exp(-(t-td1)/tc1))-(1-exp(-(t-td2)/tc2)));
	}

}
double pulse_tran(double i1,double i2,double td,double tr,double pw,double tf,double per){
	double b;
	double a,i;    

	if((t>=0)&&(t<=td)){
		return i1;
	}
	else if((t>=td)&&(t<=td+tr))
	{
		a=(i2-i1)/tr;
		b=i1-(i2-i1)*td/tr;
		i=a*t+b;
		return i;
	}
	else if((t>=td+tr)&&(t<=td+tr+pw))
	{
		return i2;
	}
	else if((t>=td+tr+pw)&&(t<=td+tr+pw+tf))
	{
		return i1;
	}
	else if((t>=td+tr+pw+tf)&&(t<=td+per))
	{
		a=(i2-i1)/tf;
		b=i2-(i2-i1)*(td+tr+pw)/tf;
		i=a*t+b;
		return i;
	}

}


double sin_tran(double i1,double ia,double fr,double td,double df,double ph){
	double i;
	if((t>=0)&&(t<=td)){
		i=i1+ia*sin(2*3.14*ph/360);
		return i;
	}
	else
	{
		i=i1+ia*sin(2*3.14*fr*(t-td)+2*3.14*ph/360)*exp((-1)*(t-td)*df);
		return i;
	}

}

double pwl_tran(){	
	int i;
	for(i=0;i<count_pwl;i++){
		if((t>=t_pwl[i])&&(t<t_pwl[i+1])){
			return i_pwl[i];
			break;
		}    
	}
}

void mna_b_tran(){
	int i,j;
	double value;
	tranB=(double*)calloc(n-1+group2,sizeof(double));

	for(i=0;i<n-1+group2;i++){
		tranB[i]=0;}

	struct list *l2;
	int k=0;
	for(l2=root->prv; l2!=root; l2=l2->prv){

		switch ( l2->type ) {
		case 'R':
			break;

		case 'V':
			k++;
			value=l2->value1;
			if(l2->trans_char=='E'){
				value=exp_tran(l2->i1,l2->i2,l2->td1,l2->tc1,l2->td2,l2->tc2);
			}
			else if (l2->trans_char=='S'){
				value=sin_tran(l2->i1,l2->ia,l2->fr,l2->td,l2->df,l2->ph);
			}
			else if (l2->trans_char=='P'){
				value=pulse_tran(l2->i1,l2->i2,l2->td,l2->tr,l2->pw,l2->tf,l2->per);
			}
			else if (l2->trans_char=='W'){
				value=pwl_tran();
			}   
			tranB[n-2+k]+=value;
			break;

		case 'I':
			value=l2->value1;
			if(l2->trans_char=='E'){
				value=exp_tran(l2->i1,l2->i2,l2->td1,l2->tc1,l2->td2,l2->tc2);
			}
			else if (l2->trans_char=='S'){
				value=sin_tran(l2->i1,l2->ia,l2->fr,l2->td,l2->df,l2->ph);
			}
			else if (l2->trans_char=='P'){
				value=pulse_tran(l2->i1,l2->i2,l2->td,l2->tr,l2->pw,l2->tf,l2->per);
			}
			else if (l2->trans_char=='W'){
				value=pwl_tran();
			}
			if((l2->first!=0)&&(l2->second!=0)){

				tranB[l2->first-1]+=-value;    
				tranB[l2->second-1]+=value; 
			}
			else if((l2->first==0)&&(l2->second!=0)){

				tranB[l2->second-1]+=value; 
			}
			else if((l2->first!=0)&&(l2->second==0)){

				tranB[l2->first-1]+=-value;    

			}
			else{}
			break;

		case 'C':

			break;

		case 'L':
			k++;
			tranB[n-2+k]+=0;
			break;

		case 'W':

			break;

		default:

			break;
		}


	}
}

void trapezoidal(){
	int i,j;
	double h;
	double **tablemion;
	double *dianusmamion;

	dianusmamion=(double*)calloc(n-1+group2,sizeof(double));
	tablemion = (double**)calloc(n-1+group2,sizeof(double*));
	for(i=0;i<n-1+group2;i++){
		tablemion[i]= (double*)calloc(n-1+group2,sizeof(double));
	}

	h=2/time_step;
	for(i=0;i<n-1+group2;i++){
		for(j=0;j<n-1+group2;j++){
			table[i][j]= table[i][j]+h*tranC[i][j];
			tablemion[i][j]=  table[i][j]-h*tranC[i][j];
		}
	}

	for(i=0;i<n-1+group2;i++){
		for(j=0;j<n-1+group2;j++){
			dianusmamion[i]=tablemion[i][j] * xpre[j];
		}
	}

	for(i=0;i<n-1+group2;i++){
		tableB[i]=tranB[i]+tranBpre[i]- dianusmamion[i] ;
	}



}
void backeuler(){

	int i,j;
	double h;
	h=1/time_step;
	double *dianusmamion;

	dianusmamion=(double*)calloc(n-1+group2,sizeof(double));
	printf("BACK_EULER\n");
	for(i=0;i<n-1+group2;i++){
		for(j=0;j<n-1+group2;j++){
			table[i][j]= table[i][j]+h*tranC[i][j];
		}
	}

	for(i=0;i<n-1+group2;i++){
		for(j=0;j<n-1+group2;j++){
			dianusmamion[i]=tranC[i][j] * xpre[j];
		}
	}

	for(i=0;i<n-1+group2;i++){
		tableB[i]=tranB[i]+h*dianusmamion[i] ;
	}

}



void pwl_sunarthsh(char line[]){
	char *str1;
	int i;
	int flag=0;
	str1 = strtok(line, " ");
	while (1){
		str1 = strtok(NULL, " ");		
		if((str1!=NULL)&&(str1[strlen(str1)-1]=='\n')){			
			str1=substring(str1, 0, strlen(str1)-2);
		}
		if (str1 == NULL){ 
			break;
		}
		if(flag==1){
			char c=*str1;
			if(c=='('){
				str1=str1=substring(str1, 2, strlen(str1));
				t_pwl[count_pwl]=atof(str1);
			}
			else{				
				str1=substring(str1, 0, strlen(str1)-1);
				i_pwl[count_pwl]=atof(str1);
				count_pwl++;
			}

		}
		if(!strcmp(str1,"PWL")){
			flag=1;
		}
	}
}

void plot_void(char linePWL[]){
	char *str1;
	int i;
	int flag=0; 

	str1 = strtok(linePWL, " ");
	while (1){
		str1 = strtok(NULL, " ");		
		if((str1!=NULL)&&(str1[strlen(str1)-1]=='\n')){			
			str1=substring(str1, 0, strlen(str1)-2);
		}
		if (str1 == NULL){ 
			break;
		}
		char c=*str1;	     
		str1=substring(str1, 3, strlen(str1)-3);		
		plot_array[count_plot]=atoi(str1);		
		count_plot++;  
	}
	flag3=1;

}
