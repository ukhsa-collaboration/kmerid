/* ***************************************************************

Takes 2 input files (numerically sorted lists of numbers)
and the Jaccard index.

Author: ulf.schaefer@phe.gov.uk 26Jul2013

*************************************************************** */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <glob.h>

#define INILISTLEN 20000000

// --------------------------------------------------------------------------------------------------------

int main(int argv, const char **args)
{
    time_t start;
    start = time(NULL);

    if (argv != 3)
    {
        printf("\nUsage: intersect_kmer_lists [readkmerlist1] [refkmerlist2]\n\n");
        exit(1);
    }

    FILE *fListFile1;
    FILE *fListFile2;

    if ((fListFile1 = fopen(args[1], "r")) == NULL)
    {
        fprintf(stderr, "Can't open file: %s\n\n", args[1]);
        exit(1);
    }
    
    if ((fListFile2 = fopen(args[2], "r")) == NULL)
    {
        fprintf(stderr, "Can't open file: %s\n\n", args[2]);
        exit(1);
    }

    long long *laList1;
    long long *laList2;
    if ((laList1=(long long*)malloc(sizeof(long long)*INILISTLEN)) == NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(2);
    }
    if ((laList2=(long long*)malloc(sizeof(long long)*INILISTLEN)) == NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(2);
    }
    memset(laList1, 0, sizeof(long long)*INILISTLEN);
    memset(laList2, 0, sizeof(long long)*INILISTLEN);

    long long llLen1=0, llLen2=0;
    while(!feof(fListFile1))
    {
        fscanf(fListFile1,"%lld\n",&laList1[llLen1]);
        llLen1++;
    }

    while(!feof(fListFile2))
    {
        fscanf(fListFile2,"%lld\n",&laList2[llLen2]);
        llLen2++;
    }    

    long long i=0, j=0, c=0, u=0;
    while (i<llLen1 && j<llLen2)
    {
        u++;
        if (laList1[i] == laList2[j])
        {
            i++; j++; c++; continue;
        }
        if (laList1[i] > laList2[j])
        {
            j++;
        }
        else
        {
            i++;
        }
    }
    
    u += (llLen1-i);    
    u += (llLen2-j);    

    long double flJacc=0.0;    
    flJacc = (long double)c / (long double)u; 
    
    printf("%Lf\t%s\t%s\n", flJacc, args[1], args[2]);

    free(laList1);
    free(laList2);
    fclose(fListFile1);
    fclose(fListFile2);

    // printf("Total processing time: %ld secs\n", time(NULL)-start);

    return 0;
}
// ------------------------------------------------------------------

// eof
