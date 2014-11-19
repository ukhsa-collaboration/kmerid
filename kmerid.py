#!/usr/bin/env python
"""

"""
import sys, argparse, subprocess, os, operator
import ConfigParser
import tempfile

__version__= '0.1'
__date__= '12Feb2014'
__author__ = 'ulf.schaefer@phe.gov.uk'

# ---------------------------------------------------------------

def parse_args():
    
    sDescription = 'version %s, date %s, author %s' %(__version__, __date__, __author__) 
 
    oParser = argparse.ArgumentParser(description=sDescription)
 
    oParser.add_argument('-f', '--fastq',
                         metavar='FILE',
                         dest='fastq',
                         required=True,
                         help='REQUIRED: Investigate this fastq file.') 

    oParser.add_argument('-c', '--config',
                         metavar='FILE',
                         dest='config',
                         required=True,
                         help='REQUIRED: Configuration file. Usually config/config.cnf.')
    
    oParser.add_argument('-n', '--nomix',
                         action='store_true',
                         dest='nomix',
                         help='Do not investigate sample for mixing. [default: Investigate. (Takes about 2 minutes.)]')    

    oArgs = oParser.parse_args()
    return oArgs, oParser

# ---------------------------------------------------------------

def main():
    oArgs, oParser = parse_args()

    oConf = ConfigParser.RawConfigParser()
    oConf.read(oArgs.config)
        
    # create kmer list for sample reads
    fTmpFile = tempfile.NamedTemporaryFile()
    createReadKmerList(os.path.abspath(oArgs.fastq), fTmpFile)
 
    dTestGenera = determineTestGenera(fTmpFile, oConf)     
    
    aResults = determineExactMatch(dTestGenera, fTmpFile, oConf)
    
    sys.stdout.write("#Kmer based similarities\n#similarity\tgroups\tfile\n")
    for aRes in aResults:
        sys.stdout.write("%f\t%s\t%s\n" % (aRes[0], aRes[3], aRes[4]))

    if oArgs.nomix == False:
        checkMixing(aResults, fTmpFile, oConf)
    
    fTmpFile.close()
    
    return

# end of main ---------------------------------------------------------------

def createReadKmerList(sFastq, fFile):
    sCmd = "cat %s | sed -n '2~4p' | bin/kmer_reads_process_stdin 18 > %s" % (sFastq, fFile.name)
    if sFastq.endswith('.gz') == True:
        sCmd = "z" + sCmd 
    p = subprocess.Popen(sCmd, shell=True, stdin=None,stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    p.wait()
    return

# ---------------------------------------------------------------

def determineTestGenera(fFile, oConf):

    aGenusResults = []    
    sCmd2 = "bin/intersect_kmer_lists_filelist " + fFile.name
    
    aGenera = oConf.options('group_folders')
    dFileToGroup = {}
    for sGen in aGenera:
        sCentSec = '%s_centroids' % sGen
        aCents  = oConf.options(sCentSec)
        sFolder = oConf.get('group_folders', sGen)
        for sCenNum in aCents:
            sCent= oConf.get(sCentSec, sCenNum)
            sCmd2 += " %s%s%s_kmers.txt" % (sFolder, os.sep, sCent)
            dFileToGroup["%s%s%s_kmers.txt" % (sFolder, os.sep, sCent)] = sGen

    p = subprocess.Popen(sCmd2, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    aOutLines = p.stdout.readlines()        
    for sLine in aOutLines:
        sLine = sLine.strip()
        aCols = [x.strip() for x in sLine.split("\t")]
        aGenusResults.append([float(aCols[0]), aCols[1], aCols[2]])
    p.stdout.close()    

    # sort results descendingly by similarity value
    aGenusResults.sort(key=operator.itemgetter(0))
    aGenusResults.reverse()
    
    dTestGenera = {}    
    # for each genus in the x highest hits
    for aHighHitGenusResult in aGenusResults[0:5]:
        # get name of similarity file
        sHighHitFolder = aHighHitGenusResult[2]
        dTestGenera[dFileToGroup[sHighHitFolder]] = 1    

    return dTestGenera

# ------------------------------------------------------------------------------

def determineExactMatch(dTestGenera, fFile, oConf):
    
    aResults = []
    dFileToGroup = {}
    # iterate over the genera to determine closest genome
    for sGen in dTestGenera.keys():
        
        sCmd3 = "bin/intersect_kmer_lists_filelist " + fFile.name
        sRefSec = '%s_refset' % sGen
        aRefs  = oConf.options(sRefSec)
        sFolder = oConf.get('group_folders', sGen)
        for sRefNum in aRefs:
            sRef= oConf.get(sRefSec, sRefNum)
            sCmd3 += " %s%s%s_kmers.txt" % (sFolder, os.sep, sRef)
            dFileToGroup["%s%s%s_kmers.txt" % (sFolder, os.sep, sRef)] = sGen
        p = subprocess.Popen(sCmd3, shell=True, stdin=None, stdout=subprocess.PIPE, close_fds=True)
        aOutLines = p.stdout.readlines()
        for sLine in aOutLines:
            sLine = sLine.strip()
            aCols = [x.strip() for x in sLine.split("\t")]
            aResults.append([float(aCols[0]), aCols[1], aCols[2]])
        p.stdout.close()
        
    # sort results array and write out results
    aResults.sort(key=operator.itemgetter(0))
    aResults.reverse()
    
    for aRes in aResults:
        aRes.append(dFileToGroup[aRes[2]])
        aRes.append(os.path.basename(aRes[2]).replace("_kmers.txt", ""))
    
    return aResults

# ------------------------------------------------------------------------------

def checkMixing(aRes, fFile, oConf):

    oOut = sys.stdout

    aTopHit = aRes[0]
    sTopHitKmerList = aTopHit[2]
    sTopHitGenus = aTopHit[3]
    sTopHitGenome = aTopHit[4]
    
    sOutput = "\n#Mixing analysis:\n#Top hit - Group: %s\tFile: %s\tSimilarity: %f%%\n" % \
              (sTopHitGenus, sTopHitGenome, aTopHit[0])    
    
    dKmerListFiles = {}
    dOrigSim = {}
    
    dFile2Group = {}
    for [flSim, dummy1, sFileName, sGroup, sFileBase] in aRes[1:]:
        dKmerListFiles[sFileName] = 1
        dOrigSim[sFileName] = flSim
        dFile2Group[sFileBase] = sGroup
    sCmd = "bin/intersect_kmer_lists_filelist %s" % sTopHitKmerList

    for k in dKmerListFiles.keys():
        sCmd += " %s" % k

    dCompResults = {}
    p = subprocess.Popen(sCmd, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    aOutLines = p.stdout.readlines()        
    for sLine in aOutLines:
        sLine = sLine.strip()
        aCols = [x.strip() for x in sLine.split("\t")]
        dCompResults[aCols[2]] = float(aCols[0])
    p.stdout.close()    
    
    aMixResults = []    
    for sF in dOrigSim.keys():
        aMixResults.append([abs(dOrigSim[sF] - dCompResults[sF]), 
                            dOrigSim[sF] - dCompResults[sF], os.path.basename(sF)]) 
    
    # sort results descendingly by similarity value
    aMixResults.sort(key=operator.itemgetter(0))
    aMixResults.reverse()    
    
    # write results 
    sOutput += "\n#Comparison of results:\n#sim diff absolute\tsim(reads,thisfile)-sim(tophit,thisfile)\tgroup\tfile\n"    
    for s in aMixResults:
        x = s[2].replace("_kmers.txt", "")
        sOutput += "%s\t%s\t%s\t%s\n" % (s[0], s[1], dFile2Group[x], x)
    oOut.write(sOutput)
    
    return

# ------------------------------------------------------------------------------

if __name__=='__main__':
    main()
