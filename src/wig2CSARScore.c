#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<R_ext/Print.h>        /* Rprintf */
#include<R_ext/Memory.h>       /* R_alloc */
#define NB_Chr_MAX 33


int wig2CSARScore(const char* file_Name,int* nbChr, int* chrL,char** filenames,char** chr){

	FILE* file;
	FILE* new_file;
	int i,j,k,m,n,counter,ret;
	int tab[NB_Chr_MAX];
	int tab_type[NB_Chr_MAX];
	int x,z;
	float y;
	int start_c,step_c,span_c;
	char str[1000];
	char chr_name[NB_Chr_MAX][200];
	char start[NB_Chr_MAX][250];
	char step[NB_Chr_MAX][250];
	char span[NB_Chr_MAX][200];

	if (*nbChr>NB_Chr_MAX)
	{
		Rprintf("ERROR [nbChr=%d]: FUNCTION WORKS FOR NB CHR<=33 YOU CAN CHANGE WITH THE SOURCE\n",*nbChr);
		return(-1);
	}

	file = fopen(file_Name,"r");
	if (file==NULL)
	{
		Rprintf("ERROR : OPEN FILE\n");
		return(-2);
	}

	j=0;
	i=0;
	m=0;

	Rprintf("READING [ %s ] : ", file_Name);
	while (feof(file)!=1)
	{
		i++;
		ret = fscanf(file," %s ",str);
		if (ret == -1) { break; }

		if (!strcmp(str,"variableStep"))
		{
			tab[j]=i;
			tab_type[j]=0;
			i=i+2;
			ret = fscanf(file," %s ",str);
			if (ret == -1) { break; }
			k=0;
			while(str[k+5]!='\0')
			{
				chr_name[j][k]=str[k+6];
				k++;
			}
			ret = fscanf(file," %s ",str);
			if (ret == -1) { break; }
			if((str[0]=='s')&&(str[1]=='p')&&(str[2]=='a')&&(str[3]=='n')){
				k=0;
				while(str[k+4]!='\0')
				{
					span[j][k]=str[k+5];
					k++;
				}
			}
			else{
				span[j][0]='N';
				span[j][1]='.';
				span[j][2]='D';
				span[j][3]='.';
				span[j][4]='\0';
			}
			j++;
		}

		if (!strcmp(str,"fixedStep"))
		{
			tab[j]=i;
			tab_type[j]=1;
			i=i+4;
			ret = fscanf(file," %s ",str);
			if (ret == -1) { break; }
			k=0;
			while(str[k+5]!='\0')
			{
				chr_name[j][k]=str[k+6];
				k++;
			}
			ret = fscanf(file," %s ",str);
			if (ret == -1) { break; }
			k=0;
			while(str[k+5]!='\0')
			{
				start[j][k]=str[k+6];
				k++;
			}
			ret = fscanf(file," %s ",str);
			if (ret == -1) { break; }
			k=0;
			while(str[k+4]!='\0')
			{
				step[j][k]=str[k+5];
				k++;
			}
			ret = fscanf(file," %s ",str);
			if (ret == -1) { break; }
			if((str[0]=='s')&&(str[1]=='p')&&(str[2]=='a')&&(str[3]=='n')){
				k=0;
				while(str[k+4]!='\0')
				{
					span[j][k]=str[k+5];
					k++;
				}
			}
			else{
				span[j][0]='N';
				span[j][1]='.';
				span[j][2]='D';
				span[j][3]='.';
				span[j][4]='\0';
			}
			j++;
		}
	}
	if ((j<*nbChr)||(j>*nbChr))
	{
		Rprintf("ERROR\n  -NB_Chr = %d AND NB_Chr_Find = %d\n  -Summary :",*nbChr,j);
		for (i=0;i<j;i++)
		{
			if(tab_type[i]==0){
				Rprintf("\n\t | %s | %s | %d |",chr_name[i],span[i],chrL[i]);
			}else{
				Rprintf("\n\t | %s | %s | %s | %s | %d |",chr_name[i],start[i],step[i],span[i],chrL[i]);
			}
		}
		return(-3);
	}
	Rprintf("done\n  -NB_Chr = %d\n  -Summary :",j);
	for (i=0;i<j;i++)
	{
        chr[i] = (char *) R_alloc(strlen(chr_name[i] + 1), sizeof(char));
		strcpy(chr[i],chr_name[i]);
		if(tab_type[i]==0){
			Rprintf("\n\t | %s | %s | %d |",chr_name[i],span[i],chrL[i]);
		}else{
			Rprintf("\n\t | %s | %s | %s | %s | %d |",chr_name[i],start[i],step[i],span[i],chrL[i]);
		}
	}
	Rprintf("\n");
	fclose(file);

	Rprintf("CREATING BINARY FILES [CSAR Bioconductor pkg format] :\n");

	file = fopen(file_Name,"r");
	if (file==NULL)
	{
		Rprintf("ERROR : OPEN FILE\n");
		return(-1);
	}
	n=1;

	for (m=0;m<j;m++)
	{
		if(tab_type[m]==0){
			sprintf(str,"%s_ChIPseq.CSARScore",chr_name[m]);
			filenames[m] = (char *) R_alloc(strlen(str) + 1, sizeof(char));
			strcpy(filenames[m], str);
			new_file=fopen(str,"wb");
			if (new_file==NULL)
			{
				Rprintf("ERROR : Unable to open file!\n");
				return(-4);
			}

			while (n<=tab[m]+1)
			{
				ret = fscanf(file," %s ",str);
				if (ret == -1) { break; }
				n++;
			}
			if (strcmp(span[m],"N.D.")){
				ret = fscanf(file," %s ",str);
				if (ret == -1) { break; }
				n++;
			}
			else{
				span[m][0]='1';
				span[m][1]='\0';
			}
			Rprintf("  - %s :",chr_name[m]);
			strcpy(str,"CSARScore v.1 Minimum");
			counter=0;
			while (str[counter]!='\0')
			{
				if (str[counter]==' ')
				{
					str[counter]='\0';
				}
				fwrite(&str[counter], sizeof(char), 1, new_file);
				counter++;
			}
			fwrite(&str[counter], sizeof(char), 1, new_file);
			x=300;
			fwrite(&x, sizeof(int), 1, new_file);
			x=1;
			fwrite(&x, sizeof(int), 1, new_file);
			x=0;
			fwrite(&x, sizeof(int), 1, new_file);
			counter=0;
			while (chr_name[m][counter]!='\0')
			{
				fwrite(&chr_name[m][counter], sizeof(char), 1, new_file);
				counter++;
			}
			fwrite(&chr_name[m][counter], sizeof(char), 1, new_file);
			fwrite(&chrL[m], sizeof(int), 1, new_file);
			span_c=atoi(span[m]);
			if(m+1<j)
			{
				i=1;
				x=0;
				z=0;
				y=0;
				while (n<tab[m+1])
				{
					ret = fscanf(file,"%d\t%f",&x,&y);
					z=0;
					while(i<x)
					{
						i++;
						fwrite(&z, sizeof(int), 1, new_file);
					}
					z=(int)(y+0.5);
					for(i=i;i<x+span_c;i++){
						fwrite(&z, sizeof(int), 1, new_file);
					}
					n=n+2;
				}
				z=0;
				for(i=i;i<=chrL[m];i++)
				{
					fwrite(&z, sizeof(int), 1, new_file);
				}
			}
			else
			{
				i=1;
				x=0;
				y=0;
				z=0;
				while (!feof(file))
				{
					ret = fscanf(file,"%d\t%f",&x,&y);
					z=0;
					while(i<x)
					{
						i++;
						fwrite(&z, sizeof(int), 1, new_file);
					}
					z=(int)(y+0.5);
					for(i=i;i<x+span_c;i++){
						fwrite(&z, sizeof(int), 1, new_file);
					}
					n=n+2;
				}
				z=0;
				for(i=i;i<=chrL[m];i++)
				{
					fwrite(&z, sizeof(int), 1, new_file);
				}
			}
		fclose(new_file);
		Rprintf(" done\n");
		}

		if(tab_type[m]==1){
			Rprintf(str,"%s_ChIPseq.CSARScore",chr_name[m]);
			new_file=fopen(str,"wb");
			if (new_file==NULL)
			{
				Rprintf("ERROR : Unable to open file!\n");
				return(-42);
			}

			while (n<=tab[m]+3)
			{
				ret = fscanf(file," %s ",str);
				if (ret == -1) { break; }
				n++;
			}
			if (strcmp(span[m],"N.D.")){
				ret = fscanf(file," %s ",str);
				if (ret == -1) { break; }
				n++;
			}
			else{
				span[m][0]='1';
				span[m][1]='\0';
			}
			Rprintf("  - %s :",chr_name[m]);
			strcpy(str,"CSARScore v.1 Minimum");
			counter=0;
			while (str[counter]!='\0')
			{
				if (str[counter]==' ')
				{
					str[counter]='\0';
				}
				fwrite(&str[counter], sizeof(char), 1, new_file);
				counter++;
			}
			fwrite(&str[counter], sizeof(char), 1, new_file);
			x=300;
			fwrite(&x, sizeof(int), 1, new_file);
			x=1;
			fwrite(&x, sizeof(int), 1, new_file);
			x=0;
			fwrite(&x, sizeof(int), 1, new_file);
			counter=0;
			while (chr_name[m][counter]!='\0')
			{
				fwrite(&chr_name[m][counter], sizeof(char), 1, new_file);
				counter++;
			}
			fwrite(&chr_name[m][counter], sizeof(char), 1, new_file);
			fwrite(&chrL[m], sizeof(int), 1, new_file);
			start_c=atoi(start[m]);
			step_c=atoi(step[m]);
			span_c=atoi(span[m]);
			if(m+1<j)
			{
				i=1;
				y=0;
				z=0;
				while (n<tab[m+1])
				{
					ret = fscanf(file,"%f",&y);
					if (ret == -1) { break; }
					z=0;
					while(i<start_c)
					{
						i++;
						fwrite(&z, sizeof(int), 1, new_file);
					}
					z=(int)(y+0.5);
					for(i=i;i<start_c+span_c;i++){
						fwrite(&z, sizeof(int), 1, new_file);
					}
					start_c=start_c+step_c;
					n=n+1;
				}
				z=0;
				for(i=i;i<=chrL[m];i++)
				{
					fwrite(&z, sizeof(int), 1, new_file);
				}
			}
			else
			{
				i=1;
				y=0;
				z=0;
				while (!feof(file))
				{
					ret = fscanf(file,"%f",&y);
					if (ret == -1) { break; }
					z=0;
					while(i<start_c)
					{
						i++;
						fwrite(&z, sizeof(int), 1, new_file);
					}
					z=(int)(y+0.5);
					for(i=i;i<start_c+span_c;i++){
						fwrite(&z, sizeof(int), 1, new_file);
					}
					start_c=start_c+step_c;
					n=n+1;
				}
				z=0;
				for(i=i;i<=chrL[m];i++)
				{
					fwrite(&z, sizeof(int), 1, new_file);
				}
			}
		fclose(new_file);
		Rprintf(" OK\n");
		}
	}
	fclose(file);
return 0;
}

