#!/usr/bin/env python
"""


"""
import sys, argparse, os, glob, subprocess
import ConfigParser
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import squareform

__version__= '0.1'
__date__= '12Feb2014'
__author__ = 'ulf.schaefer@phe.gov.uk'

# ---------------------------------------------------------------

def parse_args():
    
    sDescription = 'version %s, date %s, author %s' %(__version__, __date__, __author__) 
 
    oParser = argparse.ArgumentParser(description=sDescription)
 
    oParser.add_argument('-f', '--folder',
                         metavar='FILE',
                         dest='folder',
                         required=True,
                         help='REQUIRED: Folder that contains a set of fasta files with one reference genome each.') 

    oParser.add_argument('-n', '--name',
                         metavar='FILE',
                         dest='name',
                         required=True,
                         help='REQUIRED: Unique name for this set of references. Case insensitive. [e.g. salmonella]')

    oArgs = oParser.parse_args()
    return oArgs, oParser

# ---------------------------------------------------------------

def main():
    oArgs, oParser = parse_args()
    oConf = ConfigParser.RawConfigParser()

    sFolder = os.path.abspath(oArgs.folder)
    oArgs.name = oArgs.name.lower()
    sConfFile = "config/config.cnf"

    oConf.read(sConfFile)
    try:
        oConf.add_section('group_folders')
    except ConfigParser.DuplicateSectionError:
        pass
    oConf.set('group_folders', oArgs.name, sFolder)

    aFileEndings = ["fa", "fna", "fas"]
    aFileList = []
    for sFileEnd in aFileEndings:
        aFileList += glob.glob(os.path.join(sFolder, "*." + sFileEnd))
        aFileList += glob.glob(os.path.join(sFolder, "*." + sFileEnd + ".gz"))

    gZipped = False
    aKmerLists = []
    for sFile in aFileList:
        if sFile.endswith(".gz"):
            gZipped = True
            p = subprocess.Popen("gunzip " + sFile, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
            stdout_write("gunzipping %s ..." % sFile)
            p.wait()
            sFile = sFile[:-3]
        
        k = sFile.rfind(".")
        sKmerList = sFile[:k] + "_kmers.txt"
        if os.path.exists(sKmerList) != True:
            sCmd = "bin/kmer_refset_process %i %s > %s" % (18, sFile, sKmerList)
            p = subprocess.Popen(sCmd, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
            stdout_write("Calculating kmer list for %s ..." % sFile)
            p.wait()
        
        aKmerLists.append(sKmerList)
        if gZipped == True:
            p = subprocess.Popen("gzip " + sFile, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
            stdout_write("zipping %s ..." % sFile)
            p.wait()
        
    sSimMatFile = "config%s%s_simmat.tsv" % (os.sep, oArgs.name)
    if os.path.exists(sSimMatFile) != True:
        stdout_write("creating similarity matrix for reference group %s ..." % oArgs.name)
        create_sim_matrix(aKmerLists, sSimMatFile)
    
    try:
        oConf.add_section('matrices')
    except ConfigParser.DuplicateSectionError:
        pass
    oConf.set('matrices', oArgs.name, sSimMatFile)
    
    d3Cl = cluster_group(sSimMatFile, 3)
    d3Cen = get_centroids(d3Cl, sSimMatFile)
    
    try:
        oConf.add_section('%s_centroids' % oArgs.name)
    except ConfigParser.DuplicateSectionError:
        pass
    for k in d3Cen.keys():
        oConf.set('%s_centroids' % oArgs.name, str(k), d3Cen[k])
    
    if len(aFileList) > 40:
        d40Cl = cluster_group(sSimMatFile, 40)
        d40Cen = get_centroids(d40Cl, sSimMatFile)
        try:
            oConf.add_section('%s_refset' % oArgs.name)
        except ConfigParser.DuplicateSectionError:
            pass
        for k in d40Cen.keys():
            oConf.set('%s_refsets' % oArgs.name, str(k), d40Cen[k])
    else:
        aGenomes = [os.path.basename(s).replace("_kmers.txt", "") for s in aKmerLists]
        try:
            oConf.add_section('%s_refset' % oArgs.name)
        except ConfigParser.DuplicateSectionError:
            pass
        for i in range(1, len(aGenomes)+1):
            oConf.set('%s_refset' % oArgs.name, str(i),aGenomes[i-1])
    
    fCnf = open(sConfFile, 'w')
    oConf.write(fCnf)
    fCnf.close()
    
    return

# end of main ---------------------------------------------------

def stdout_write(s):
    sys.stdout.write("%s\n" % s)
    sys.stdout.flush()
    return

# ---------------------------------------------------------------

def create_sim_matrix(aFiles, sSimMat):

    fOut = open(sSimMat, 'w')
    
    iNofFiles = len(aFiles)
    aNames = [os.path.basename(x).replace("_kmers.txt", "") for x in aFiles]    

    for s in aNames:
        fOut.write("\t%s" % (s))
    fOut.write("\n")
    
    d = {}    
    for i in range(0, iNofFiles):
        fOut.write(aNames[i])
        for j in range(0, iNofFiles):
            if i==j:
                fOut.write("\t1.0")
            else:
                flSim = 0.0
                k1 = str(i) + "-" + str(j)
                k2 = str(j) + "-" + str(i)
                try:
                    fOut.write("\t%f" % (d[k2]))                
                except KeyError:
                    sCmd = "bin/kmer_jaccard_index %s %s" % (aFiles[i], aFiles[j])
                    p = subprocess.Popen(sCmd, shell=True, stdin=None, stdout=subprocess.PIPE, close_fds=True)
                    # get output from subprocess
                    sOutLine = p.stdout.readline()
                    flSim = float(sOutLine.split("\t")[0])
                    fOut.write("\t%f" % (flSim))
                    d[k1] = flSim                
        
        fOut.write("\n")

    fOut.close()
    return

# ---------------------------------------------------------------

def cluster_group(sFile, iClusters):
    
    f = open(sFile, 'r')
    a = []
    aLines = []
    for s in f:
        aLines.append(s.strip())  
    f.close()
    
    aNames = aLines[0].strip().split("\t")
    for sL in aLines[1:]:
        a.append([1.0-float(x) for x in sL.split("\t")[1:]])
    
    linkage_matrix = linkage(squareform(a), method='average')    
    cluster = fcluster(linkage_matrix, float(iClusters), criterion='maxclust')
    
    dClus = {}    
    for x in range(0, len(cluster)):
        clusnum = cluster[x]        
        try:        
            dClus[clusnum].append(aNames[x])
        except KeyError:
            dClus[clusnum] = [aNames[x]]

    return dClus

# ---------------------------------------------------------------

def get_centroids(dClusters, sSimMatFile):

    fMat = open(sSimMatFile, 'r')
    aaSimMatrix = []
    for sLine in fMat:
        aaSimMatrix.append([x.strip() for x in sLine.split("\t")])     
    fMat.close()

    aNames = aaSimMatrix[0][1:]
    dSimMatrix = {}    
    for aLine in aaSimMatrix[1:]:
        sName = aLine[0]        
        dSimMatrix[sName] = {}        
        for i in range(1, len(aLine)):
            dSimMatrix[sName][aNames[i-1]] = float(aLine[i])

    dCentroids = {}    
    
    for k in dClusters.keys():
        if len(dClusters[k]) <= 2:
            dCentroids[k] = dClusters[k][0]
        else:
            dAvgDists = {}
            for x in dClusters[k]:
                flAvgDist = 0.0
                for y in dClusters[k]:
                    flAvgDist += dSimMatrix[x][y]
                flAvgDist /= float(len(dClusters[k]))
                dAvgDists[flAvgDist] = x
            flMin = max(dAvgDists.keys())
            dCentroids[k] = dAvgDists[flMin]

    return dCentroids

# ---------------------------------------------------------------

if __name__=='__main__':
    main()
