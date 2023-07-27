import os, time, sys
from os import system, remove, path

def main():
#######################################
# Option to create a case             #
#######################################
#######################################
    if(len(sys.argv)<4):
        print 'Please use option -h to see all available options'
        return

    if(len(sys.argv)==4):
        nameString = '%s' %sys.argv[2]
        modeID = '%s' %sys.argv[3]
        print sys.argv[1]
        print 'SUGAR post-processing tool operating...' 
        print 'For %s th modeID' %(modeID)
        bfn = 'SpatialMode_%s_'%(nameString) #baselineFileName
        rf = 'TecFile_SpatialMode_%s_'%(nameString) #readlineFile

#        myRange = range(0, int(MPIRanks))
        removeMainFile = 'rm %s%s.tec' %(rf,modeID)
        system(removeMainFile)
        appendCommand = 'cat my1DMesh.dat >> %s%s.tec' %(rf,modeID)
        system(appendCommand)

        appendCommand = 'cat %s%s.dat >> %s%s.tec' %(bfn, modeID, rf, modeID)
        system(appendCommand)
#        print 'File #i is appended to %s%s.tec' %(bfn,planeID)

#        appendCommand = 'cat ../2DMeshGenerator/connectivity.dat >> %s%s.tec' %(rf,modeID)
#        system(appendCommand)
        if(sys.argv[1]=='keep'):
            print 'Partitioned files are kept since KEEP option was used'
        elif(sys.argv[1]=='remove'):
            print 'Partitioned files are removed since REMOVE option was used'
        else:
            print 'Unknown argument used, partitioned files are kept for safety'
        
    elif(sys.argv[1]=='-h'):
        print 'help option is invoked.'
        print 'LiberLocus multiscale multiphysics solver version 1.0.'
        print 'To post-process a case use 8000 32, the planeID and total number of MPI ranks'
if __name__ == "__main__":
    main()
