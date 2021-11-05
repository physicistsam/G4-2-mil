# this converts gadget 4 output trees to millennium format
#to run this $$ python3.8 conv_G4_mil.py <i/p treefile name> <# snapshots> <o/p number of treefile>
# i/p tree file := Gadged 4 final tree file
import os, sys, time
import math as m
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import h5py 

# check if all the arguments are there :
I = ["<i/p treefile name>","<# snapshots>"," <# treefilenum>"]
for ii in range (1, 4):
    if (sys.argv[ii] == None):
        print(I[ii-1],"is missing.", "\n", "Ending Now !")
        exit()


# Read final tree file :
in_path ="/data/chatterjee/gadget_run_sept2021/DM-L102-Nbody-512_seed1/output/"
out_path="./data_512/"
# code assumes there is a ./data/ directory which has all the input files
#
#? call first snapshot file for mass of the particles
snapfile=in_path+str('snapshot_000.hdf5')
snapf=h5py.File(snapfile) # read tree file
partmass = snapf["Header"].attrs["MassTable"][1]
print("particle mass =", partmass)
snapf.close()
tfname = in_path+str(sys.argv[1]) # Gadget 4 final tree file.
treef=h5py.File(tfname) # read tree file
print("Tree file keys :")
print(list(treef.keys()), "\n\n")
#
nsnap = int(sys.argv[2]) # number of napshot files, (0, nsnap-1)
print("fof_subhalo keys:")
#
#now read fof_subhalo_tab_###.hdf5
fshf = [] # list of fof_subhalo_tab files
for ii in range (nsnap):
    fsfname=in_path+'fof_subhalo_tab_%03d.hdf5'%ii
    #fshf.append(h5py.File('data/fof_subhalo_tab_%03d.hdf5'%ii))
    fshf.append(h5py.File(fsfname))
    print("File%03d"%ii,"\n", list(fshf[ii].keys()),"\n")
#
print("reading HDF5 files done.")

# Make list of the quantities required from fof_subhalo files :
Nsub = []
GSGN = [] # for /Group/Group_M_Crit200
GSGM = [] # for /Group/Group_M_Mean200
GSGT = [] # for /Group/Group_M_TopHat200

for ii in range (nsnap):
    GSGN.append(ii)
    GSGM.append(ii)
    GSGT.append(ii)
    #
    Nsubtf= fshf[ii]["Header"].attrs["Nsubhalos_ThisFile"]
    print("Nsubtf", Nsubtf)
    if (Nsubtf > 0):
    	GSGN[ii] = list(fshf[ii]['/Group/Group_M_Crit200'])
    	GSGM[ii] = list(fshf[ii]['/Group/Group_M_Mean200'])
    	GSGT[ii] = list(fshf[ii]['/Group/Group_M_TopHat200'])
    
    	Nsub.append(len(list(fshf[ii]['/Subhalo/SubhaloLen'])))
    else :
    	GSGN[ii] = [0]
    	GSGM[ii] = [0]
    	GSGT[ii] = [0]
    
    	Nsub.append(0)
    #print((GSGN[ii][100]))
print("mass collection from fof_subhalo files ... done.\n\n")
print("Start making list for trees ..... \n")
#
HID=list(treef['TreeHalos/TreeID'])
TINX=list(treef['TreeHalos/TreeIndex'])
Len = list(treef['TreeTable/Length'])
DEC = list(treef['TreeHalos/TreeDescendant'])
FP = list(treef['TreeHalos/TreeFirstProgenitor'])
NP = list(treef['TreeHalos/TreeNextProgenitor'])
FH = list(treef['TreeHalos/TreeFirstHaloInFOFgroup'])
NH = list(treef['TreeHalos/TreeNextHaloInFOFgroup'])
NPOS = list(treef['TreeHalos/SubhaloPos'])
POS = [element * 1000 for element in NPOS]
VEL = list(treef['TreeHalos/SubhaloVel'])
VELD = list(treef['TreeHalos/SubhaloVelDisp'])
VMAX = list(treef['TreeHalos/SubhaloVmax'])
NSPIN = list(treef['TreeHalos/SubhaloSpin'])
SPIN = [element * 1000 for element in NSPIN]
MBID = list(treef['TreeHalos/SubhaloIDMostbound'])
MVIR = list(treef['TreeHalos/Group_M_Crit200'])
SNUM = list(treef['TreeHalos/SnapNum'])
SLEN = list(treef['TreeHalos/SubhaloLen'])
GGN = list(treef['TreeHalos/GroupNr'])
GSN= list(treef['TreeHalos/SubhaloNr'])
#FNR = list(np.ones(np.sum(Len), dtype=int))
#FNR = list(np.ones((len(list(treef['TreeTable/TreeID']))), dtype=int))
RSD = list(treef['TreeTimes/Redshift'])
#
print("Done making list ! \n")

#print(len(FNR), FNR[10:20])
print("Start making mass list from fof_subhalo data ..... \n")
GMM200 = []# for /Group/Group_M_Mean200
GMT200 = []# for /Group/Group_M_TopHat200
NT = len(list(treef['TreeTable/TreeID']))
NHC = 0
for ii in range (0, NT):
    TID = HID[ii]
    LL = Len[ii]

    for jj in range(0, LL):
        ssn=int(SNUM[NHC])
        if (MVIR[NHC] > 0 ):
            Nc = GSGN[ssn][:].index(MVIR[NHC])
            gmm200 = GSGM[ssn][Nc]
            gmt200 = GSGT[ssn][Nc]
            
            GMM200.append(gmm200)
            GMT200.append(gmt200)
        else :
            GMM200.append(0)
            GMT200.append(0)

        NHC += 1
print("Done making mass list ! \n")     

print("write data into a file .... \n")
#

Nfile = int(sys.argv[3])
Cnhpf = np.zeros(Nfile)
Cntpf = np.zeros(Nfile)
i_count = 0
op =[]
tnh = np.empty(Nfile, dtype=np.object)

for ff in range (Nfile):
    opfname = out_path+"trees_%d.%d.hdf5"%(nsnap,ff)            #str(sys.argv[3])
    op.append(h5py.File(opfname, 'w'))
    tnh[ff] = []

    #op[ff].close()

    # TotNsubhalos = grp.create_dataset("TotNsubhalos", (1,))
    # TotNsubhalos[:] = [0]

    # TreeNHalos = grp.create_dataset("TreeNHalos", (NT,))
    # TreeNHalos[:] = Len[:]
NT = len(list(treef['TreeTable/TreeID']))
NHC = 0 ## subhalo number count
#'Descendant', 'FileNr', 'FirstHaloInFOFGroup', 'FirstProgenitor', '
#Group_M_Crit200', 'Group_M_Mean200', 'Group_M_TopHat200', 'NextHaloInFOFGroup', 
# 'NextProgenitor', 'SnapNum', 'SubhaloGrNr', 'SubhaloHalfmassRad', 'SubhaloHalfmassRadType', 
# 'SubhaloIDMostBound', 'SubhaloLen', 'SubhaloLenType', 'SubhaloMassInRadType', 
# 'SubhaloMassType', 'SubhaloNumber', 'SubhaloOffsetType', 'SubhaloPos', 'SubhaloSpin', 
# 'SubhaloVMax', 'SubhaloVel', 'SubhaloVelDisp' *
for ii in range (NT):
    TID = HID[ii]
    LL = Len[ii]
    
    #? choose file number for the 
    fnum = int(ii%Nfile)
    TGname = "Tree%d"%Cntpf[fnum]#+str(ii)
    grp =  op[fnum].create_group(TGname)
    ## adding data sets:
    Descendant = grp.create_dataset("Descendant", data=np.int32(DEC[NHC:(NHC+LL)]), dtype='i4')
    FirstProgenitor = grp.create_dataset("FirstProgenitor", data=np.int32(FP[NHC:(NHC+LL)]), dtype='i4')
    NextProgenitor = grp.create_dataset("NextProgenitor", data=np.int32(NP[NHC:(NHC+LL)]), dtype='i4')
    FirstHaloInFOFGroup = grp.create_dataset("FirstHaloInFOFGroup", data= np.int32(FH[NHC:(NHC+LL)]), dtype='i4')
    NextHaloInFOFgroup = grp.create_dataset("NextHaloInFOFGroup", data= np.int32(NH[NHC:(NHC+LL)]), dtype='i4')
    SubhaloLen = grp.create_dataset("SubhaloLen", data=np.int32(SLEN[NHC:(NHC+LL)]), dtype='i4')
    Group_M_Crit200 = grp.create_dataset("Group_M_Crit200", data=np.float32(MVIR[NHC:(NHC+LL)]), dtype='f4')
    SubhaloPos = grp.create_dataset("SubhaloPos", data=np.float32(POS[NHC:(NHC+LL)]), dtype='f4') # needs to be converted to kpc units, required by Sage-model)
    SubhaloVel= grp.create_dataset("SubhaloVel", data=np.float32(VEL[NHC:(NHC+LL)]), dtype='f4')
    SubhaloVelDisp = grp.create_dataset("SubhaloVelDisp", data=np.float32(VELD[NHC:(NHC+LL)]), dtype='f4')
    SubhaloVMax =  grp.create_dataset("SubhaloVMax", data=np.float32(VMAX[NHC:(NHC+LL)]), dtype='f4')
    SubhaloSpin  =  grp.create_dataset("SubhaloSpin", data=np.float32(SPIN[NHC:(NHC+LL)]), dtype='f4')
    SubhaloIDMostBound =  grp.create_dataset("SubhaloIDMostBound", data=np.ulonglong(MBID[NHC:(NHC+LL)]), dtype='u8')
    SnapNum =  grp.create_dataset("SnapNum", data=np.int32(SNUM[NHC:(NHC+LL)]), dtype='i4')
    #FileNr = grp.create_dataset("FileNr", data=np.int32(FNR[NHC:(NHC+LL)]))
    FNR = [fnum]*LL
    FileNr = grp.create_dataset("FileNr", data=np.int32(FNR[:]), dtype='i4')
    ## these following quantities need to be changed:
    Group_M_Mean200 = grp.create_dataset("Group_M_Mean200", data=np.float32(GMM200[NHC:(NHC+LL)]), dtype='f4')
    Group_M_TopHat200 = grp.create_dataset("Group_M_TopHat200", data=np.float32(GMT200[NHC:(NHC+LL)]), dtype='f4')
    ## increase Halo count
    NHC += LL
    Cntpf[fnum] += 1
    Cnhpf[fnum] += LL
    tnh[fnum].append(LL)

for ff in range (Nfile):
    grp = op[ff].create_group("Header")
    grp.attrs['AlphaExponent'] =0.8 #may not be required for Sage
    grp.attrs['FirstSnapshotNr'] =np.int32(0)
    grp.attrs['LastSnapshotNr'] =np.int32(nsnap-1)
    #grp.attrs['totNHalos'] =np.int32(np.sum(Len))
    # grp.attrs['NhalosPerFile'] =np.sum(Len)
    # grp.attrs['Ntrees'] =NT
    grp.attrs['NumberOfOutputFiles'] = np.int32(Nfile) # may be pass as an input for general case
    grp.attrs['ParticleMass'] =np.float64(partmass) #check ?
    grp.attrs['RunOutputDir'] =np.str("./output/") #not required for now
    grp.attrs['SnapSkipFac'] =np.int32(1) #neeed to be checked
    grp.attrs['SnapshotFileBase'] =np.str("snap")#not required for now
    grp.attrs['TreeOutputDir'] =np.str("./output/")#not required for now

    Redshifts = grp.create_dataset("Redshifts", (nsnap,), dtype='f8')
    Redshifts[:] = np.float64(RSD[:])#[3.04356, 0.978456, 0.575981, 4.44089e-16]
    #grp.attrs['totNHalos'] =np.sum(Len)
    grp.attrs['NhalosPerFile'] = np.int32(Cnhpf[ff])
    grp.attrs['NtreesPerFile'] =np.int32(Cntpf[ff])

    TotNsubhalos = grp.create_dataset("TotNsubhalos", (nsnap,), dtype='i4')
    TotNsubhalos[:] = np.int32(Nsub[:])
    print("Match ?",Cntpf[ff], len(tnh[ff]), "and", Cnhpf[ff], np.sum(tnh[ff][:]))
    TreeNHalos = grp.create_dataset("TreeNHalos", (Cntpf[ff],), dtype='i4')
    TreeNHalos[:] = np.int32(tnh[ff][:])
    op[ff].close()

print("Finished writing data in %s\n\n"%opfname)
print(".. Bye ..")
