/* ***************************************************************

Creates a sorted list of unique 18mers in a fasta file. All sequences 
in the file are considered to be one big contiguous sequence. Only counts
the alphabetically 'smaller' one between the forward and the reverse
complement of each 18mer.

Author: ulf.schaefer@phe.gov.uk 26Jun2013

*************************************************************** */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define INILINELEN 1000
#define INISEQLEN 10000

void q_sort(long long *a, long long llL, long long llR);
long rmdup(long long *a, long lLen);
long long median3(long long a, long long b, long long c);
int compare (const void * a, const void * b); 
// --------------------------------------------------------------------------------------------------------

int main(int argv, const char **args)
{
    time_t start;
    start = time(NULL);

    if (argv != 3)
    {
        printf("\nUsage: kmer_refset_process [kmerlen] [file.fa]\n\n");
        exit(1);
    }

    int KMERLEN=atoi(args[1]);
    int KMERLENMINUSONE = KMERLEN-1;
    long long nNofKmers=(long long)pow(4.0, KMERLEN);

    FILE *fFastaFile;

    if ((fFastaFile = fopen(args[2], "r")) == NULL)
    {
        fprintf(stderr, "Can't open file: %s\n\n", args[2]);
        exit(1);
    }

    char sLine[INILINELEN];
    memset(sLine, '\0', sizeof(char)*INILINELEN);

    char *cpSeq=0, *cpSeq2=0;
    int iAvailLen = INISEQLEN;    
    if ((cpSeq=(char*)malloc(sizeof(char)*INISEQLEN)) == NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(2);
    }
    memset(cpSeq, '\0', sizeof(char)*INISEQLEN);

    long i=0, j=0, x=0, iCurLen=0;
    int *ipSeq=0;
    int iLineLen=0;
    while (fgets(sLine, INILINELEN, fFastaFile))
    {
        if (sLine[0] == '>')
        {
            memset(sLine, '\0', sizeof(char)*INILINELEN);
            continue;
        }

        iLineLen = strlen(sLine)-1;
        
        if (iCurLen + iLineLen > iAvailLen)
        {
            if ((cpSeq2=(char*)malloc(sizeof(char)*(iAvailLen*2))) == NULL)
            {
                fprintf(stderr, "Memory allocation failed\n");
                exit(2);
            }
            iAvailLen *= 2;      
            memset(cpSeq2, '\0', sizeof(char)*iAvailLen);
            cpSeq -= iCurLen;
            memcpy(cpSeq2, cpSeq, sizeof(char)*iCurLen);
            free(cpSeq);
            cpSeq = cpSeq2;
            cpSeq += iCurLen;                       
        }
        
        memcpy(cpSeq, sLine, sizeof(char)*iLineLen);
        cpSeq += iLineLen;
        iCurLen += iLineLen;
        
        memset(sLine, '\0', sizeof(char)*INILINELEN);
    }        
    
    fclose(fFastaFile);    

    cpSeq -= iCurLen;
    
    long lNofKmerPos = (iCurLen - KMERLEN) +  1;   
    
    if ((ipSeq=(int*)malloc(sizeof(int)*iCurLen)) == NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(2);
    }
    memset(ipSeq, 0, sizeof(int)*iCurLen);
    
    // printf("\n%ld\n", iCurLen);    
    
    // convert read string to int array            
    for (j = 0; j < iCurLen; j++)
    {
        switch (toupper(cpSeq[j]))
        {
            case 'A':
                *ipSeq = 0;
                break;
            case 'C':
                *ipSeq = 1;
                break;
            case 'G':
                *ipSeq = 2;
                break;
            case 'T':
                *ipSeq = 3;
                break;
            default:
                *ipSeq = -1;
                break;
        }
        ipSeq += 1;
    }
    ipSeq -= iCurLen;

    free(cpSeq);    

    long long *llpKmers = 0;
    if ((llpKmers=(long long*)malloc(sizeof(long long)*lNofKmerPos)) == NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(2);
    }

    memset(llpKmers, 0, sizeof(long long)*lNofKmerPos);

    int iValid = 1;
    long double nKmerFWD=0.0, nKmerRVC=0.0;
    
    long long llThis=0, llPrev=-1;    
    
    i=0;    
    for (j = 0; j < lNofKmerPos; j++)
    {
        nKmerFWD=0.0;
        nKmerRVC=0.0;
        iValid = 1;
        // forward strand        
        for (x = KMERLENMINUSONE; x >= 0; x--)
        {
            // check for invalid chars go to next k-mer if found                    
            if (ipSeq[j+x] == -1) // N
            {
                iValid = 0;
                break;
            }
            // convert kmer from base 4 to base 10 index for counting
            nKmerFWD += ipSeq[j+x] * pow(4.0, KMERLENMINUSONE-x);
        }
        // reverse complement strand
        if (iValid!=0)
        {
            for (x = 0; x <= KMERLENMINUSONE; x++)
            {
                // convert kmer from base 4 to base 10 index for counting
                nKmerRVC += (3-ipSeq[j+x]) * pow(4.0, x);
            }
        }

        if (iValid!=0)
        {
            llThis = (long long)((nKmerFWD <= nKmerRVC) ? nKmerFWD : nKmerRVC);            
            if (llThis != llPrev)
            {            
                llpKmers[i] = llThis;
                i++;
                llPrev = llThis;
            }
        }
    }
    free(ipSeq);
    
    lNofKmerPos  = i;
    
    // printf("\nProcessing time until sort: %ld secs\n", time(NULL)-start);

    // use stdlib qsort instead of own nonsense    
    // q_sort(llpKmers, 0, i-1);
    qsort (llpKmers, i, sizeof(long long), compare);
    
    // printf("\nProcessing time including sort: %ld secs\n", time(NULL)-start);
    // printf("OldSize: %ld\n", i);
    
    long lNewSize=0;
    lNewSize = rmdup(llpKmers, lNofKmerPos);

    // output kmer
    for (i=0; i<lNewSize; i++)
        printf("%lld\n",llpKmers[i]);
    
    
    // printf("NewSize: %ld\n", lNewSize );    
    // printf("Total processing time: %ld secs\n", time(NULL)-start);

    free(llpKmers);

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

long rmdup(long long *a, long lLen)
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

// ----------------------------------------------------------------------------

void q_sort(long long *a, long long llL, long long llR)
{
    int x=0; long y=0;
    for (y=llL; y<=llR; y++)
        if (a[y] != 0)
            x =1;  
    
    if (x==0) return;    
    
    long long llPivot=0,
              llTmp1=0,
              llTmp2=0;
 
    llTmp1 = llL;
    llTmp2 = llR;
    llPivot = a[llL];
    while (llL < llR)
    {
        while ((a[llR] >= llPivot) && (llL < llR))
        {
            llR--;
        }        
        if (llL != llR)
        {
            a[llL] = a[llR];
            llL++;
        }
        while ((a[llL] <= llPivot) && (llL < llR))
        {
            llL++;
        }
        if (llL != llR)
        {
            a[llR] = a[llL];
            llR--;
        }
    }
    a[llL] = llPivot;
    llPivot = llL;
    llL = llTmp1;
    llR = llTmp2;
    if (llL < llPivot)
    {
        q_sort(a, llL, llPivot-1);
    }
    if (llR > llPivot)
    {
        return q_sort(a, llPivot+1, llR);
    }
    
    return;
}

// ----------------------------------------------------------------------------

// eof
