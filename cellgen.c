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

typedef struct tempcoord{
	float x;
	float y;
	float z;
	int av;	//available
}tempcoord;

int readinitial(const char* grofile,struct molecula *in);
int writegro (const char* grofile,const struct molecula *smol,const struct systemcell sce,const float a, const float b, const float c);

int main(int argc, char *argv[]){
	struct molecula moltype[MAXMOL];
	struct systemcell ce;
	struct tempcoord *TempP;
	int i,j,k,l;
	int ii,jj,kk,ll;
	float dens;
	int nmol;
	int Lcell;
	int Hcell;
	float Dcell;
	float Vcell;
	float dlcell;
	int id;
	int tempint;
	float MemDelta,MemLat;
	
	int NMem,NMol;
	int NSub;
	int Nset;
	int testm;
	int* XMol;
	int Inserted;
	char* outfile;
	
	srand (time(NULL));
	
	if(argc<3){
		printf("usege: density filename.gro \n"); //print help
		return 0;
	}
	if(strcmp(argv[1],"il-ph")==0){
	dens=strtof(argv[2],NULL);
	if(dens>0.0){
		printf("number density %f \n",dens);
	}
	else{
		printf("%s error: number density <0 \n %s",KRED,KNRM);
	}
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
		dens=dens/SIGMA/SIGMA/SIGMA;
		nmol=(int)1000.0/dens*40.0;
		Hcell=(int)pow(nmol,1.0/3.0)+2;	//printf(" %f \n",pow(nmol,2));
		printf("%i %i \n",nmol,Hcell);
		nmol=Hcell*Hcell*Hcell;
		Vcell=nmol/dens;
		dlcell=pow(Vcell,1.0/3.0)/Hcell;
		ce.atype=(int*)malloc(nmol*Lcell*MAXATOM*sizeof(int));
		ce.x=(float**)malloc(nmol*Lcell*sizeof(float*));
		ce.y=(float**)malloc(nmol*Lcell*sizeof(float*));
		ce.z=(float**)malloc(nmol*Lcell*sizeof(float*));
		for(i=0;i<nmol*Lcell;i++){
			ce.x[i]=(float*)malloc(MAXATOM*sizeof(float));
			ce.y[i]=(float*)malloc(MAXATOM*sizeof(float));
			ce.z[i]=(float*)malloc(MAXATOM*sizeof(float));
		}
		id=0;
		for(ii=0;ii<Hcell*Lcell;ii++){		//z
			for(jj=0;jj<Hcell;jj++){		//x
				for(kk=0;kk<Hcell;kk++){	//y
					if((id+1)%2==0){
						ce.atype[id]=0;
					}
					else{
						ce.atype[id]=1;
					}
					for(ll=0;ll<moltype[ce.atype[id]].anum;ll++){
						ce.x[id][ll]=moltype[ce.atype[id]].x[ll]+(j+0.5)*dlcell;
						ce.y[id][ll]=moltype[ce.atype[id]].y[ll]+(k+0.5)*dlcell;
						ce.z[id][ll]=moltype[ce.atype[id]].z[ll]+(i+0.5)*dlcell;
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
	else if (strcmp(argv[1],"il")==0){
		printf("%s generation for ionic liquid cube cell %s \n",KBLU,KNRM);
		return 0;
	}
	else if(strcmp(argv[1],"mem-dif")==0){
		if(argc<5){
			printf("%s error usege: cellgen mem-dif Deltamem MemLat dens file.gro %% %s \n",KRED,KNRM);
		}
		printf("%s generation of membranes for diffusion %s \n",KBLU,KNRM);
		//calculate voloume of cell
		MemDelta=strtof(argv[2],NULL);
		MemDelta=MemDelta*SIGMA;	//set the dimensionless form
		MemLat=strtol(argv[3],NULL,10);
		dens=strtof(argv[4],NULL);
		dlcell=MemDelta*MemLat;
		Vcell=pow(MemDelta*MemLat,3);
		NMem=MemLat*MemLat*MemLat;	//number of membrane molecules
		NMol=(int)Vcell*dens/(SIGMA*SIGMA*SIGMA);	//number of molecules
		printf("Cell voloume: %f \n", Vcell);
		//numbers of substances
		NSub=(int)(argc-5)/2;
		printf("Number of substances %d, total molecules %d, membrane atoms %d \n",NSub,NMol,NMem);
		//set coordintes
		Nset=(int)sqrt(NMol/MemLat)+1;	//(int)pow(NMol,1.0/3.0)+1;
		printf("Nset %d \n", Nset);
		TempP=(tempcoord*)malloc(NMol*sizeof(tempcoord));
		id=0;
		for(ii=0;ii<MemLat;ii++){	//set initial coordintes of atoms
			for(jj=0;jj<Nset;jj++){
				for(kk=0;kk<Nset;kk++){
					//
					if(id<NMol){
						TempP[id].x=(ii+0.5)*dlcell/MemLat;
						TempP[id].y=(jj)*dlcell/Nset;
						TempP[id].z=(kk)*dlcell/Nset;
						TempP[id].av=0;
						id++;
					}
				}
			}
		}
		XMol=(int*)malloc(NSub*sizeof(int));
		//allocate atoms
		nmol=NMem+NMol;	//molecules with membrane
		ce.atype=(int*)malloc(nmol*MAXATOM*sizeof(int));
		ce.x=(float**)malloc(nmol*sizeof(float*));
		ce.y=(float**)malloc(nmol*sizeof(float*));
		ce.z=(float**)malloc(nmol*sizeof(float*));
		for(i=0;i<nmol;i++){
			ce.x[i]=(float*)malloc(MAXATOM*sizeof(float));
			ce.y[i]=(float*)malloc(MAXATOM*sizeof(float));
			ce.z[i]=(float*)malloc(MAXATOM*sizeof(float));
		}
		printf("test \n");
		//
		id=0;	//number of inserted molecule
		for(i=0;i<NSub;i++){
			//read initial
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
			tempint=readinitial(argv[i*2+5],&moltype[i]);	//read gro file
			XMol[i]=strtof(argv[i*2+6],NULL);
			
			//
			XMol[i]=XMol[i]*NMol/100;
			printf("substance %d number of molecules %d \n",i,XMol[i]);
			Inserted=0;
			while(Inserted<XMol[i]){
				testm=rand()%NMol;
				//printf("%d \n",testm);
				if(TempP[testm].av==0){
					//insert molecules
					TempP[testm].av==1;
					ce.atype[id]=i;
					for(ii=0;ii<moltype[i].anum;ii++){
						ce.x[id][ii]=moltype[i].x[ii]+TempP[testm].x;
						ce.y[id][ii]=moltype[i].y[ii]+TempP[testm].y;
						ce.z[id][ii]=moltype[i].z[ii]+TempP[testm].z;
					}
					id++;
					Inserted++;
					//printf("inserted %d XMol %d \n",Inserted,XMol[i]);
				}
			//printf("id %d %d\n",id,i);
			}
			
		}
		moltype[NSub].x=(float*)malloc(MAXATOM*sizeof(float));
		moltype[NSub].y=(float*)malloc(MAXATOM*sizeof(float));
		moltype[NSub].z=(float*)malloc(MAXATOM*sizeof(float));
		moltype[NSub].vx=(float*)malloc(MAXATOM*sizeof(float));
		moltype[NSub].vy=(float*)malloc(MAXATOM*sizeof(float));
		moltype[NSub].vz=(float*)malloc(MAXATOM*sizeof(float));
		moltype[NSub].name=(char**)malloc(MAXATOM*sizeof(char*));
		moltype[NSub].atom=(char**)malloc(MAXATOM*sizeof(char*));
		for(j=0;j<MAXATOM;j++){
			moltype[NSub].name[j]=(char*)malloc(MAXSTR*sizeof(char));
			moltype[NSub].atom[j]=(char*)malloc(MAXSTR*sizeof(char));
		}
		//insert membrane atoms
		//Membrane moleculas
		
		moltype[NSub].anum=1;
		moltype[NSub].x[0]=0.0;
		moltype[NSub].y[0]=0.0;
		moltype[NSub].z[0]=0.0;
		moltype[NSub].vx[0]=0.0;
		moltype[NSub].vy[0]=0.0;
		moltype[NSub].vz[0]=0.0;
		moltype[NSub].name[0]="Mem";
		moltype[NSub].atom[0]="Mem";
		for(ii=0;ii<MemLat;ii++){
			for(jj=0;jj<MemLat;jj++){
				for(kk=0;kk<MemLat;kk++){
					ce.atype[id]=NSub;
					ce.x[id][0]=ii*MemDelta;
					ce.y[id][0]=jj*MemDelta;
					ce.z[id][0]=kk*MemDelta;
					id++;
				}
			}
		}
		//
		ce.mnum=NMol+NMem;
		if(argv[argc-1]==NULL){
			outfile="out.gro";
		}
		else{
			outfile=argv[argc-1];
		}
		printf("%s %s DONE %s \n",KBLU,outfile,KNRM);
		writegro(outfile,moltype,ce,dlcell,dlcell,dlcell);

		free(ce.atype);
		free(ce.x);
		free(ce.y);
		free(ce.z);
		return 0;
	}
	else {
		printf("Unknown type %s \n",argv[2]);
		printf("%sGeneration type: il - ionic licuid cube \n il-ph - ionic liquid \n mem-dif - for membrane diffusion %s \n",KRED,KNRM);
		return 0;
	}
}

int readinitial(const char* grofile,struct molecula *in){
	FILE *infile;
	int arr;
	int tempint;
	int i;
	char temps[30];
	//
	infile=fopen(grofile,"r");
	if (infile==0){
		printf("%s  %s file open error %s \n", grofile,KRED,KNRM);
		return 1;
	}
	fscanf(infile,"%s",temps);
	fscanf(infile,"%d",&in->anum);
	printf(" %d atoms in %s \n",in->anum,grofile);
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

