/****************************************************************

Takes a list of input kmer list files and prints the similarity between the
first list (assumed to be from reads) and all other lists (assumed to be from a 
set of reference genomes).

Author: ulf.schaefer@phe.gov.uk 31Jul2013
Modified: tobingallop@yahoo.co.uk 20Nov2018

****************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <glob.h>
#include <dirent.h>
#include <libgen.h>

#define VERSION 0.2

void displayUsage(char*);
long long countFileLines(FILE*);

//---------------------------------------------------------------
int main(int argc,  char *argv[])
{
 if (argc < 3 ) {
  displayUsage(argv[0]);
  exit(1);
 }

 int k = 0;    
 long long i = 0, j = 0, c = 0;
 long long llLen1 = 0, llLen2 = 0, lines = 0;
 float flSim = 0.0, flDist = 0.0;
 FILE *fListFile1; 
 FILE *fListFile2;
 long long *laList1, *laList2;

 if ((fListFile1 = fopen(argv[1], "r")) == NULL) {
  fprintf(stderr, "Can't open file: %s\n", argv[1]);
  exit(1);
 }
 
 lines = countFileLines(fListFile1);
 //fprintf(stdout, "%s contains %lld lines\n", argv[1], lines);
 
 if ((laList1 = (long long*)malloc(sizeof(long long) * lines)) == NULL) {
  fprintf(stderr, "Memory allocation failed\n");
  exit(2);
 }
 memset(laList1, 0, sizeof(long long) * lines);
 
 fseek(fListFile1, 0, SEEK_SET);
 while(!feof(fListFile1)) {
  fscanf(fListFile1, "%lld\n", &laList1[llLen1]);
  llLen1++;
 }

 for (k = 2; k < argc; k++) {                  
  if ((fListFile2 = fopen(argv[k], "r")) == NULL) {
   fprintf(stderr, "Can't open file: %s\n", argv[k]);
   exit(1);
  }
  
  lines = countFileLines(fListFile2);
  //fprintf(stdout, "%s contains %lld lines\n", argv[k], lines);
  
  if ((laList2 = (long long*)malloc(sizeof(long long) * lines)) == NULL) {
   fprintf(stderr, "Memory allocation failed\n");
   exit(2);
  }
  memset(laList2, 0, sizeof(long long) * lines);

  llLen2 = 0;
  fseek(fListFile2, 0, SEEK_SET);
  while(!feof(fListFile2)) {
   fscanf(fListFile2, "%lld\n", &laList2[llLen2]);
   llLen2++;
  }        

  i = 0, j = 0, c = 0;
  while (i < llLen1 && j < llLen2) {
   if (laList1[i] == laList2[j]) {
    i++; j++; c++; continue;
   }
   if (laList1[i] > laList2[j]) {
    j++;
   } else {
    i++;
   }
  }

  // Richa: "Similarity is simply percentage of 18mers in reference seen in read set as well."
  flSim = (float)c / ( (float)llLen2 / 100.0);

  // Similarity = percentage of kmers in reads seen in reference as well     
  // flSim = (float)c / ( (float)llLen1 / 100.0); 

  flDist = 100.0 - flSim;
  fprintf(stdout, "%f\t%f\t%s\n", flSim, flDist, argv[k]);
  fclose(fListFile2);
  free(laList2);
 }

 fclose(fListFile1);
 free(laList1);
 return 0;
}

//---------------------------------------------------------------
void displayUsage(char* parent)
{
 char *app;
 char *path = strdup(parent);
 app = basename(path);
 fprintf(stdout, "\n%s v%0.1f\n", app, VERSION);
 fprintf(stdout, "Usage: %s [readkmerlist] [refkmerlist_1] [refkmerlist_2] ... [refkmerlist_n]\n", app);
 fprintf(stdout, " [refkmerlist]       - File containing a list of sorted kmers. This is generate off of a fastq\n");
 fprintf(stdout, "                     - file (reads) and used to investigate similarities against the set of\n");
 fprintf(stdout, "                     - reference genomes\n");
 fprintf(stdout, " [refkmerlist_1,2,n] - List of files containing sorted kmers. These files are the reference\n");
 fprintf(stdout, "                     - genomes used to compare against the first kmer list (reads)\n");
}

//---------------------------------------------------------------
long long countFileLines(FILE *fp)
{
 char ch;
 long long linecount = 0;
 
 // check if file is open
 if (fp) {
  // read character from file until EOF
  while ((ch = getc(fp)) != EOF) {
   if (ch == '\n')
    ++linecount;
  }
 } else {
  fprintf(stderr, "File is not open\n");
 }
 
 return linecount;
}