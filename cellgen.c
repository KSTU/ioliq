#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KBLU  "\x1B[34m"
#define MAXSTR  32
#define MAXATOM  500
#define MAXMOL  10
#define SIGMA  0.3

typedef struct molecula{
	char **name;
	char **atom;
	int anum;	//atom numbers
	float *x;
	float *y;
	float *z;
	float *vx;
	float *vy;
	float *vz;
}molecula;

typedef struct systemcell{
	int mnum;	//	number of molecules in cell 
	int *atype;	//	number of atoms in moltype 
	float **x;
	float **y;
	float **z;
	
}systemcell;

int readinitial(const char* grofile,struct molecula *in);
int writegro (const char* grofile,const struct molecula *smol,const struct systemcell sce,const float a, const float b, const float c);

int main(int argc, char *argv[]){
	struct molecula moltype[MAXMOL];
	struct systemcell ce;
	int i,j,k,l;
	float dens;
	int nmol;
	int Lcell;
	int Hcell;
	float Dcell;
	float Vcell;
	float dlcell;
	int id;
	int tempint;
	float rbez;
	
	if(argc<3){
		printf("usege: density filename.gro \n");
		return 0;
	}
	dens=strtof(argv[1],NULL);
	if(dens>0.0){
		printf("number density %f \n",dens);
	}
	
	if(strcmp(argv[2],"il-ph")==0){
	
	for(i=0;i<argc-3;i++){	//read initial atom configuration
		moltype[i].x=(float*)malloc(MAXATOM*sizeof(float));
		moltype[i].y=(float*)malloc(MAXATOM*sizeof(float));
		moltype[i].z=(float*)malloc(MAXATOM*sizeof(float));
		moltype[i].vx=(float*)malloc(MAXATOM*sizeof(float));
		moltype[i].vy=(float*)malloc(MAXATOM*sizeof(float));
		moltype[i].vz=(float*)malloc(MAXATOM*sizeof(float));
		moltype[i].name=(char**)malloc(MAXATOM*sizeof(char*));
		moltype[i].atom=(char**)malloc(MAXATOM*sizeof(char*));
		for(j=0;j<MAXATOM;j++){
			moltype[i].name[j]=(char*)malloc(MAXSTR*sizeof(char));
			moltype[i].atom[j]=(char*)malloc(MAXSTR*sizeof(char));
		}
		tempint=readinitial(argv[i+3],&moltype[i]);
		}
		printf("%s generation for ionic liquid phase equilibrium %s \n",KBLU,KNRM);
		Lcell=3;
		rbez=0.3;
		dens=dens/rbez/rbez/rbez;
		nmol=(int)1000.0/dens;
		Hcell=(int)pow(nmol,1.0/3.0)+1; //printf(" %f \n",pow(nmol,2));
		printf("%i %i \n",nmol,Hcell);
		nmol=Hcell*Hcell*Hcell;
		Vcell=nmol/dens;
		dlcell=pow(Vcell,1.0/3.0)/Hcell;
		ce.atype=(int*)malloc(nmol*Lcell*MAXATOM*sizeof(int));
		ce.x=(float**)malloc(nmol*Lcell*sizeof(float*));
		for(i=0;i<nmol*Lcell;i++){
			ce.x[i]=(float*)malloc(MAXATOM*sizeof(float*));
		}
		ce.y=(float**)malloc(nmol*Lcell*sizeof(float*));
		for(i=0;i<nmol*Lcell;i++){
			ce.y[i]=(float*)malloc(MAXATOM*sizeof(float));
		}
		ce.z=(float**)malloc(nmol*Lcell*sizeof(float*));
		for(i=0;i<nmol*Lcell;i++){
			ce.z[i]=(float*)malloc(MAXATOM*sizeof(float));
		}
		id=0;
		for(i=0;i<Hcell*Lcell;i++){ //z
			for(j=0;j<Hcell;j++){	//x
				for(k=0;k<Hcell;k++){	//y
					if((id+1)%2==0){
						ce.atype[id]=0;
					}
					else{
						ce.atype[id]=1;
					}
					//
					//printf("%i %i \n",id,moltype[ce.atype[id]].anum);
					for(l=0;l<moltype[ce.atype[id]].anum;l++){
						ce.x[id][l]=moltype[ce.atype[id]].x[l]+(j+0.5)*dlcell;
						ce.y[id][l]=moltype[ce.atype[id]].y[l]+(k+0.5)*dlcell;
						ce.z[id][l]=moltype[ce.atype[id]].z[l]+(i+0.5)*dlcell;
						//printf("name %s %d l %d x %f y %f z %f \n",moltype[ce.atype[id]].name[l], ce.atype[id],l,ce.x[id][l],ce.y[id][l],ce.z[id][l]);
					}
					id++;
				}
			}
			
		}
		ce.mnum=id--;
		printf("%s DONE %s \n",KBLU,KNRM);
		writegro("out.gro",moltype,ce,dlcell*Hcell,dlcell*Hcell,dlcell*Hcell*Lcell);
		//
		free(ce.atype);
		free(ce.x);
		free(ce.y);
		free(ce.z);
		//
		return 0;
	}
	else if (strcmp(argv[2],"il")==0){
		printf("%s generation for ionic liquid cube cell %s \n",KBLU,KNRM);
		return 0;
	}
	else {
		printf("Unknown type %s \n",argv[2]);
		printf("%sGeneration type: il - ionic licuid cube, il-ph - ionic liquid %s \n",KRED,KNRM);
		return 0;
	}
}

int readinitial(const char* grofile,struct molecula *in){
	FILE *infile;
	int arr;
	int tempint;
	int i;
	//
	infile=fopen(grofile,"r");
	if (infile==0){
		printf("%s  %s file open error %s \n", grofile,KRED,KNRM);
		return 1;
	}
	fscanf(infile,"%d",&in->anum);
	//in.anum=3; //&tempint;
	printf(" %d atoms in %s \n",in->anum,grofile);
	//
	for(i=0;i<in->anum;i++){
		fscanf(infile,"%5d%5s%5s%5d%f%f%f%f%f%f",
		&tempint,in->name[i],in->atom[i],&tempint,&in->x[i],&in->y[i],&in->z[i],&in->vx[i],&in->vy[i],&in->vz[i]);
		printf(" \t %d %s %s x %f y %f z %f \n", tempint,in->name[i],in->atom[i],in->x[i],in->y[i],in->z[i]);
	}
	//"%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
	fclose(infile);
	return 0;
}

int writegro (const char* grofile,const struct molecula *smol,const struct systemcell sce,
	const float a, const float b, const float c){
	FILE *outfile;
	int i,j;
	int TotMol;
	//
	TotMol=0;
	for(i=0;i<sce.mnum;i++){
		for(j=0;j<smol[sce.atype[i]].anum;j++){
			TotMol++;
		}
	}
	outfile=fopen(grofile,"w");
	fprintf(outfile," generated by iolic \n");
	fprintf(outfile," %d \n",TotMol);
	for(i=0;i<sce.mnum;i++){
		for(j=0;j<smol[sce.atype[i]].anum;j++){
			//printf("%d %s %s \n",sce.atype[i],smol[sce.atype[i]].name[j],smol[sce.atype[i]].atom[j]);
			fprintf(outfile,"%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",i+1,smol[sce.atype[i]].name[j],
				smol[sce.atype[i]].atom[j],
				j+1,sce.x[i][j],sce.y[i][j],sce.z[i][j],0.0,0.0,0.0
				);
			//printf("%d \n", sce.atype[i]);
		}
	}
	fprintf(outfile,"    %9.5f    %9.5f    %9.5f  \n", a,b,c);
	fclose(outfile);
}

