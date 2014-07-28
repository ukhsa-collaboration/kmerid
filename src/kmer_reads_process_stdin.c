/* ***************************************************************

Creates a sorted list of unique kmers in a fastq file. Only considers
the alphabetically 'smaller' one between the forward and the reverse
complement of each kmer.

Reads from stdin. Only pipe in actual reads. E.g.:

cat file.fq | sed -n '2~4p' | kmer reads_process_stdin 18 > outfile.txt

Memory requirements for large fastq file might be an issue.
(> 800MB on 4GB RAM is a [soft] limit) 

Author: ulf.schaefer@phe.gov.uk 24Jun2013

*************************************************************** */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define INILINELEN 1000
#define ININOFKMERS 1000000
#define ININOFREADS 500000

void q_sort(long long *a, long long llL, long long llR);
void convert_to_numerical(char *sLine, long lLen, int **ipSeq);
long rmdup(long long *a, long long lLen);
int compare (const void * a, const void * b);

// --------------------------------------------------------------------------------------------------------

int main(int argv, const char **args)
{
    time_t start;
    start = time(NULL);

    if (argv != 2)
    {
        printf("\nUsage: kmer_reads_process [kmerlen]\n\n");
        exit(1);
    }

    int KMERLEN=atoi(args[1]);
    int KMERLENMINUSONE = KMERLEN-1;
    long nNofKmers=(long)pow(4.0, KMERLEN);

    char sLine[INILINELEN];
    memset(sLine, '\0', sizeof(char)*INILINELEN);

    long long *llpKmers=0, *llpKmers2=0;

    long iAvailReadSpace=0, iNofReads=0, i=0;
    long long lTotalKmers=0;
    
    int iLineLen=0;
    char** aReads;
    char** aReads2=0;
    if ((aReads=(char**)malloc(sizeof(char*)*ININOFREADS)) == NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(2);
    }           
    memset(aReads, '\0', sizeof(char*)*ININOFREADS);
    iAvailReadSpace = ININOFREADS;
    while (fgets(sLine, INILINELEN, stdin))
    {
        iLineLen = strlen(sLine)-1;
        lTotalKmers += (iLineLen - KMERLENMINUSONE);
    
        if ((iNofReads + 1) > iAvailReadSpace)        
        {
            if ((aReads2=(char**)malloc(sizeof(char*)*(iAvailReadSpace*2))) == NULL)
            {
                fprintf(stderr, "Memory allocation failed\n");
                exit(2);
            }
            memset(aReads2, '\0', sizeof(char*)*(iAvailReadSpace*2));
            memcpy(aReads2, aReads, sizeof(char*)*(iNofReads));    
            free(aReads);
            aReads = aReads2;
            aReads2 = 0;
            iAvailReadSpace *= 2;
        }
        
        if ((aReads[iNofReads]=(char*)malloc(sizeof(char)*(iLineLen+1))) == NULL)
        {
            fprintf(stderr, "Memory allocation failed\n");
            exit(2);
        }        
        memset(aReads[iNofReads], '\0', sizeof(char)*(iLineLen+1));
        memcpy(aReads[iNofReads], sLine, sizeof(char)*(iLineLen+1));    
        
        iNofReads++;
    }
    
    // printf("%lld\n", lTotalKmers);
    // printf("%li\n", sizeof(long long));   
    
    if ((llpKmers=(long long*)malloc(sizeof(long long)*lTotalKmers)) == NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(2);
    }
    memset(llpKmers, 0, sizeof(long long)*lTotalKmers);    
    
    long long iaPowFour[KMERLEN];
    long k=0;    
    for (k = 0; k < KMERLEN; k++)
    {
        iaPowFour[k] = (long long)pow(4.0, k);
    }
    
    long p=0;
    int j=0, iReadLen=0, x=0, iValid = 1;
    long long q=0, nKmerFWD=0, nKmerRVC=0;
    int *ipSeq;
    // while (fgets(sLine, INILINELEN, stdin))
    for (i=0; i<iNofReads; i++)
    {
        iReadLen = strlen(aReads[i])-1;
        ipSeq=(int*)malloc(sizeof(int)*iReadLen);
        memset(ipSeq, 0, sizeof(int)*iReadLen);
        convert_to_numerical(aReads[i], iReadLen, &ipSeq);

        for (j=0; j<iReadLen-KMERLENMINUSONE; j++)
        { 
            nKmerFWD=0;
            nKmerRVC=0;
            iValid = 1;                
            for (p = KMERLENMINUSONE; p >= 0; p--)
            {
                // check for invalid chars go to next k-mer if found                    
                if (ipSeq[j+p] == -1) // N
                {
                    iValid = 0;
                    break;
                }
                // convert kmer from base 4 to base 10 index for counting
                nKmerFWD += (ipSeq[j+p] * iaPowFour[KMERLENMINUSONE-p]);
            }
       
            // reverse complement strand
            if (iValid!=0)
            {
                for (p = 0; p <= KMERLENMINUSONE; p++)
                {
                    // convert kmer from base 4 to base 10 index for counting
                    nKmerRVC += (3-ipSeq[j+p]) * iaPowFour[p];
                }
                
                llpKmers[q] = (nKmerFWD <= nKmerRVC) ? nKmerFWD : nKmerRVC;
                q++;
            }           
        }
        free(ipSeq);
                   
        // memset(sLine, '\0', sizeof(char)*INILINELEN);
    }

    // printf("\nProcessing time excluding sort: %ld secs\n", time(NULL)-start);

    // use stdlib qsort instead of own nonsense    
    // q_sort(llpKmers, 0, q-1);
    qsort (llpKmers, q, sizeof(long long), compare);
        
    long long a=0, b=0;
    for (a=0; a<q-1; a++)
    {
        if(llpKmers[a] == llpKmers[a+1])
        {
            // printf("%lld\t", llpKmers[a]);
            a++;
            b++;     
        }    
    }    
    
    long long *llpNonUniqKmers=0;
    if ((llpNonUniqKmers=(long long*)malloc(sizeof(long long)*b)) == NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(2);
    }
    memset(llpNonUniqKmers, 0, sizeof(long long)*b);  
    
    b=0;
    for (a=0; a<q-1; a++)
    {
        if(llpKmers[a] == llpKmers[a+1])
        {
            llpNonUniqKmers[b] = llpKmers[a];
            a++;
            b++;     
        }    
    }
    free(llpKmers);
    
    long lNewSize=0;
    lNewSize = rmdup(llpNonUniqKmers, b);

    // output kmer
    // printf("\t%lld", llpNonUniqKmers[0]);
    for (a=0; a<lNewSize; a++)
        printf("%lld\n", llpNonUniqKmers[a]);
    // printf("\n");
    
    // printf("NewSize: %ld\n", lNewSize );  
        
    for (i=0; i<iNofReads; i++)
        free(aReads[i]);    
    free(aReads);
    
    // printf("\nTotal processing time: %ld secs\n", time(NULL)-start);
        
    return 0;
}

// ----------------------------------------------------------------------------

// helper function for the stdlib qsort
int compare (const void * a, const void * b)
{
  if ( *(long long*)a <  *(long long*)b ) return -1;
  if ( *(long long*)a == *(long long*)b ) return 0;
  if ( *(long long*)a >  *(long long*)b ) return 1;
}

// ----------------------------------------------------------------------------

long rmdup(long long *a, long long lLen)
{
    long i=1, j=0;
    
    for (; i < lLen; i++) 
    { 
        if (a[i] != a[j]) 
        { 
            j++; 
            a[j] = a[i]; 
        } 
    }  

    return j+1;
}
// --------------------------------------------------------------------------------------------------------

void convert_to_numerical(char *sLine, long lLen, int **ipSeq)
{
    long j=0;
    // convert read string to int array            
    for (j = 0; j < lLen; j++)
    {
        switch (toupper(sLine[j]))
        {
            case 'A':
                **ipSeq = 0;
                break;
            case 'C':
                **ipSeq = 1;
                break;
            case 'G':
                **ipSeq = 2;
                break;
            case 'T':
                **ipSeq = 3;
                break;
            default:
                **ipSeq = -1;
                break;
        }
        *ipSeq += 1;
    }

    *ipSeq -= lLen;

    return;

}


// eof
