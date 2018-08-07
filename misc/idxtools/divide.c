#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <assert.h>


#define MAX_CHR 1024                       //maximum number of chromosomes
#define MAX_CHR_NAME 256                   //maximum number of characters in a chromosome name
#define OUTPUT_FILE_FORMAT "part%d.fa"     //output file format


int NUM_PARTS=0;                      //maximum number of parts

/*Die on error. Print the error and exit if the return value of the previous function NULL*/
#define errorCheckNULL(ret) ({\
    if (ret==NULL){ \
        fprintf(stderr,"Error at File %s line number %d : %s\n",__FILE__, __LINE__,strerror(errno));\
        exit(EXIT_FAILURE);\
    }\
    })

/*Die on error. Print the error and exit if the return value of the previous function is -1*/
#define errorCheck(ret) ({\
    if (ret<0){ \
        fprintf(stderr,"Error at File %s line number %d : %s\n",__FILE__, __LINE__,strerror(errno));\
        exit(EXIT_FAILURE);\
    }\
    })
    
    
struct chromsomeData{
    char chr_name [MAX_CHR_NAME];   //name of chromosome
    long chr_len;                   //length of the chromosome in bases
    long chr_len_no_N ;             //length of the chromosome in bases (no N)
};    
    
    
//compare function for descending sort    
int cmpfunc (const void * a, const void * b) {
   struct chromsomeData *A = (struct chromsomeData *)a;
   struct chromsomeData *B = (struct chromsomeData *)b;
   return ( (B->chr_len_no_N) - (A->chr_len_no_N) );
}    

//find the index of the minimum value in an array
int minimum_index(long *list, int listsize){
    long minimum=LONG_MAX;
    int min_index=-1;
    int j;
    for(j=0;j<listsize;j++){
        if(list[j]<minimum){
            min_index=j;
            minimum=list[j];
        }
    }    
    assert(min_index>=0);
    return min_index;
    
}
  
//find to which part a certain chromosome(chr_name) is belonging to
int belong_to_which_part(char *chr_name, struct chromsomeData chr_per_part[][MAX_CHR],int *numchr_per_part){
    int i,j;
    for(i=0;i<NUM_PARTS;i++){
        for(j=0;j<numchr_per_part[i];j++){
            if(strcmp(chr_per_part[i][j].chr_name,chr_name)==0){
                return i; 
            }
        }
    }
    assert(0);
    return -1;

}    
    
int main(int argc, char **argv){
	
    if(argc!=3){
        fprintf(stderr,"Usage ./%s <reference.fa> <numparts>\n",argv[0]);
        exit(EXIT_FAILURE);
    }
    
    //open the fasta
    FILE *fasta = fopen(argv[1],"r");
    errorCheckNULL(fasta);
	
	//number of parts
	NUM_PARTS=atoi(argv[2]);
	if(NUM_PARTS<2){
		fprintf(stderr,"Number of parts should be greater than 2\n");
		exit(EXIT_FAILURE);
	}

	
    
    //for getline
    size_t bufferSize = 1000;
	char *buffer = malloc(sizeof(char)*bufferSize);
    errorCheckNULL(buffer);
    size_t readlinebytes = 0;
    
    int num_chr=0;  //keep track of number of chromosomes in the fasta
    struct chromsomeData theGenome[MAX_CHR];  //information for the whole genome
    
    //variables for current chromosome being processed
    long chr_len;
    long chr_len_no_N;
    char chr_name[MAX_CHR_NAME];
  
  
    //read through the file
    while(1){
        
        readlinebytes=getline(&buffer, &bufferSize, fasta); 
        if(readlinebytes==-1){  //EOF
            if (num_chr!=0){    //should save the data for the previous chromosome    
                theGenome[num_chr-1].chr_len=chr_len;
                theGenome[num_chr-1].chr_len_no_N=chr_len_no_N;
                strcpy(theGenome[num_chr-1].chr_name,chr_name);
                fprintf(stderr,"Processed %s\n",theGenome[num_chr-1].chr_name);
            }            
            break;
        }
        if(readlinebytes==0){
            fprintf(stderr,"We read nothing. Why is that?\n");
            exit(EXIT_FAILURE);
        }
        
        //a new chromosme
        if(buffer[0]=='>'){     
            if (num_chr!=0){    //if not the first chromosome in the file, save the data of the previous chromosome that was processed    
                theGenome[num_chr-1].chr_len=chr_len;
                theGenome[num_chr-1].chr_len_no_N=chr_len_no_N;
                strcpy(theGenome[num_chr-1].chr_name,chr_name);
                fprintf(stderr,"Processed %s\n",theGenome[num_chr-1].chr_name);
            }
            
            num_chr++;
            
            //reset lengths for the new chromosome
            chr_len=0;
            chr_len_no_N=0;
            
            if(readlinebytes-1>MAX_CHR_NAME){
                fprintf(stderr,"Chromosome name too large, Increase MAX_CHR_NAME\n");
                exit(EXIT_FAILURE);   
            }
            strcpy(chr_name, &buffer[1]); //copy the chromosome name except the >
            if(chr_name[strlen(chr_name)-1]=='\n'){
                chr_name[strlen(chr_name)-1]='\0';
            }
            
            else{
                fprintf(stderr,"Bad new line character\n");
            }
           
            if(num_chr>MAX_CHR){
                fprintf(stderr,"So many chromosomes, Increase MAX_CHR\n");
                exit(EXIT_FAILURE);
            }
        }
        
        //going through the chromosome
        else{
            int i=0;
            //go through all bases
            for(i=0;i<readlinebytes;i++){
                if(buffer[i]=='A' || buffer[i]=='C' || buffer[i]=='G' || buffer[i]=='T' || buffer[i]=='a' || buffer[i]=='c' || buffer[i]=='g' || buffer[i]=='t' ){
                    chr_len++;
                    chr_len_no_N++;
                }
				else if (buffer[i]=='R' || buffer[i]=='Y' || buffer[i]=='K' || buffer[i]=='M' || buffer[i]=='S' || buffer[i]=='W' || buffer[i]=='B' || buffer[i]=='D' || buffer[i]=='H' || buffer[i]=='V' || \
						 buffer[i]=='r' || buffer[i]=='y' || buffer[i]=='k' || buffer[i]=='m' || buffer[i]=='s' || buffer[i]=='w' || buffer[i]=='b' || buffer[i]=='d' || buffer[i]=='h' || buffer[i]=='v'){
                    chr_len++;
                    chr_len_no_N++;					
				}
                else if(buffer[i]=='N' || buffer[i]=='n'){
                    chr_len++;
                }                
                else if(buffer[i]!='\n'){
                    fprintf(stderr,"Funny character in %s : |%c|\n",chr_name,buffer[i]);
                    chr_len++;
                }
            }
            
        }
    }
    
    fclose(fasta); 
    
    long sum=0;
    fprintf(stderr,"\nNumber of chromosomes : %d\n",num_chr);
    
    //print stats
    FILE *stat=fopen("stat.csv","w");
    errorCheckNULL(stat);    
   
    fprintf(stat,"Chromosome name,Chromosome length,Chromosome length (without N)\n");
    int i;
    for(i=0;i<num_chr;i++){
        fprintf (stat,"\"%s\",%ld,%ld\n",theGenome[i].chr_name, theGenome[i].chr_len, theGenome[i].chr_len_no_N);
        sum+=theGenome[i].chr_len_no_N;
    }
    fclose(stat);
    
    //sort based on lengths
    //long per_part = sum/NUM_PARTS;
    qsort(theGenome,num_chr,sizeof(struct chromsomeData),cmpfunc);
    
    //for(i=0;i<num_chr;i++){
    //    fprintf (stderr,"\"%s\",%ld,%ld\n",theGenome[i].chr_name, theGenome[i].chr_len, theGenome[i].chr_len_no_N);

    //}    
    
    
    /*********************** Now do the partitioning******************************************/
    
    struct chromsomeData chr_per_part[NUM_PARTS][MAX_CHR]; //per each part we need to store the chromomes to be processed
    int numchr_per_part[NUM_PARTS];                        //number of chromosomes for each part        
    long length_per_part[NUM_PARTS];                       //the total length of chromomes (only ACGT) for each part        
    
    for(i=0;i<NUM_PARTS;i++){
        numchr_per_part[i]=0;
        length_per_part[i]=0;
        
    }
    
    int j;
    int min_index=0;

    //go through all the chromosomes
    //add the current chromosome to the list with lowest sum
    for(i=0;i<num_chr;i++){
        min_index = minimum_index(length_per_part,NUM_PARTS);
        chr_per_part[min_index][numchr_per_part[min_index]]=theGenome[i];
        numchr_per_part[min_index]++;
        length_per_part[min_index]+=theGenome[i].chr_len_no_N;
    }
    
    char filename[1024];
    
    //print the stats for each part
    for(i=0;i<NUM_PARTS;i++){
        fprintf(stderr,"\nFor part %d we have %d chromosomes with a sum length of %ld neucleotides (only ACGT)\n",i,numchr_per_part[i],length_per_part[i]);
        sprintf(filename,"stat_part%d.csv",i);
        stat=fopen(filename,"w");
        errorCheckNULL(stat); 
        fprintf(stat,"Chromosome name,Chromosome length,Chromosome length (without N)\n");
        for(j=0;j<numchr_per_part[i];j++){
            fprintf(stat,"%s,%ld,%ld\n",chr_per_part[i][j].chr_name, chr_per_part[i][j].chr_len, chr_per_part[i][j].chr_len_no_N);
        }
        fclose(stat);
        
    }
    
    /*********************** Now write the fastas for each part******************************************/

    fasta = fopen(argv[1],"r");
    errorCheckNULL(fasta);    
    
    FILE *outputs[NUM_PARTS];
    
    for(i=0;i<NUM_PARTS;i++){
        sprintf(filename,OUTPUT_FILE_FORMAT,i);
        outputs[i]=fopen(filename,"w");
        errorCheckNULL(outputs[i]);
    }
 
    int part_index=-1;
 
    //read each line of the input file
    while(1){
        readlinebytes=getline(&buffer, &bufferSize, fasta); 
        if(readlinebytes==-1){
            break;
        }
        if(readlinebytes==0){
            fprintf(stderr,"We read nothing. Why is that?\n");
            exit(EXIT_FAILURE);
        }
        
        if(buffer[0]=='>'){

            if(readlinebytes-1>MAX_CHR_NAME){
                fprintf(stderr,"Chromosome name too large, Increase MAX_CHR_NAME\n");
                exit(EXIT_FAILURE);   
            }
            strcpy(chr_name, &buffer[1]);
            if(chr_name[strlen(chr_name)-1]=='\n'){
                chr_name[strlen(chr_name)-1]='\0';
            }
            else{
                fprintf(stderr,"Bad new line character\n");
            }
            
            //find to which part this chromosome must be allocated to
            part_index=belong_to_which_part(chr_name, chr_per_part,numchr_per_part);
            fprintf(outputs[part_index],"%s",buffer); 
            fprintf(stderr,"Printing chromosome *%s* to part %d\n",chr_name,part_index);
                        

        }
        else{
            //write to the correct file
            fprintf(outputs[part_index],"%s",buffer);  
            
        }
    }
    
    for(i=0;i<NUM_PARTS;i++){
        fclose(outputs[i]);
    }
    fclose(fasta); 
    free(buffer);

    
	return 0;
}
