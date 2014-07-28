/* ***************************************************************

Takes a list of input kmer list files and prints the similarity between the
first list (assumed to be from reads) and all other lists (assumed to be from a 
set of reference genomes).

Author: ulf.schaefer@phe.gov.uk 31Jul2013

*************************************************************** */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <glob.h>
#include <dirent.h>

#define INILISTLEN 20000000

// --------------------------------------------------------------------------------------------------------

int main(int argv, const char **args)
{
    time_t start;
    start = time(NULL);

    if (argv < 3 )
    {
        printf("\nUsage: intersect_kmer_lists_filelist [readkmerlist] [refkmerlist_1] [refkmerlist_2] ... [refkmerlist_n]\n\n");
        exit(1);
    }
      
    FILE *fListFile1;

    if ((fListFile1 = fopen(args[1], "r")) == NULL)
    {
        fprintf(stderr, "Can't open file: %s\n\n", args[1]);
        exit(1);
    }
        
    long long *laList1, *laList2;
    if ((laList1=(long long*)malloc(sizeof(long long)*INILISTLEN)) == NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(2);
    }
    memset(laList1, 0, sizeof(long long)*INILISTLEN);

    if ((laList2=(long long*)malloc(sizeof(long long)*INILISTLEN)) == NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(2);
    }
    memset(laList2, 0, sizeof(long long)*INILISTLEN);
    
    long long llLen1=0, llLen2=0;
    while(!feof(fListFile1))
    {
        fscanf(fListFile1,"%lld\n",&laList1[llLen1]);
        llLen1++;
    }

    int k=0;    
    long long i=0, j=0, c=0;
    float flSim=0.0, flDist=0.0;
    FILE *fListFile2;
    for (k=2; k<argv; k++)
    {                  
        if ((fListFile2 = fopen(args[k], "r")) == NULL)
        {
            fprintf(stderr, "Can't open file: %s\n\n", args[k]);
            exit(1);
        }

        llLen2 = 0;
        while(!feof(fListFile2))
        {
            fscanf(fListFile2,"%lld\n",&laList2[llLen2]);
            llLen2++;
        }        

        i=0, j=0, c=0;
        while (i<llLen1 && j<llLen2)
        {
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

        // Richa: "Similarity is simply percentage of 18mers in reference seen in read set as well."
        flSim = (float)c / ( (float)llLen2 / 100.0);
        
        // Similarity = percentage of kmers in reads seen in reference as well     
        // flSim = (float)c / ( (float)llLen1 / 100.0); 
                
        flDist = 100.0 - flSim;

        printf("%f\t%f\t%s\n", flSim, flDist, args[k]);
        
        fclose(fListFile2);
        memset(laList2, 0, sizeof(long long)*INILISTLEN);
    }
    
    free(laList1);
    free(laList2);
    fclose(fListFile1);

    return 0;
}

// ------------------------------------------------------------------

// eof
