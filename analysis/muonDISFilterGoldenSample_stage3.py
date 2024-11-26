import os
# Code snippet from Simona. This should go somewhere where it can be used by several different pieces of code (monitoring, analysis, etc)
import pickle

import numpy as np
import ROOT
import rootUtils as ut
import SndlhcGeo

p = open("/eos/experiment/sndlhc/convertedData/commissioning/TI18/FSdict.pkl",'rb')
FSdict = pickle.load(p)

def GetAvgScifiPos(DigiScifiHits):
    n_sf_hits_x = [0]*5
    n_sf_hits_y = [0]*5
    avg_sf_x = [0.]*5
    avg_sf_y = [0.]*5
    a, b = ROOT.TVector3(), ROOT.TVector3()
    for hit in DigiScifiHits :
        if not hit.isValid() :
            continue
        plane = hit.GetStation() - 1
        if hit.isVertical() :
            n_sf_hits_x[plane] += 1
            avg_sf_x[plane] += (a.X() + b.X())/2.
        else :
            n_sf_hits_y[plane] += 1
            avg_sf_y[plane] += (a.Y() + b.Y())/2.
    for i_plane in range(5) :
        if n_sf_hits_x[i_plane] :
            avg_sf_x[i_plane] /= n_sf_hits_x[i_plane]
        if n_sf_hits_y[i_plane] :
            avg_sf_y[i_plane] /= n_sf_hits_y[i_plane]
    return avg_sf_x, avg_sf_y

def getAvgScifiPos(event, scifiDet):
    n_sf_hits_x ={1:0, 2:0, 3:0, 4:0, 5:0}
    n_sf_hits_y ={1:0, 2:0, 3:0, 4:0, 5:0}
    avg_sf_x = {1:0, 2:0, 3:0, 4:0, 5:0}
    avg_sf_y = {1:0, 2:0, 3:0, 4:0, 5:0}
    a, b = ROOT.TVector3(), ROOT.TVector3()
    for aHit in event.Digi_ScifiHits:
        if not aHit.isValid(): continue
        plane = aHit.GetStation()
        detID = aHit.GetDetectorID()
        scifiDet.GetSiPMPosition(detID, a, b)
        if aHit.isVertical():
            n_sf_hits_x[plane]+=1
            avg_sf_x[plane]+= (a.X() + b.X())/2.
        else:
            n_sf_hits_y[plane]+=1
            avg_sf_y[plane]+= (a.Y() + b.Y())/2.
    for iplane in range(1, 5):
        if n_sf_hits_x[iplane]:
            avg_sf_x[plane]/=n_sf_hits_x[iplane]
        if n_sf_hits_y[iplane]:
            avg_sf_y[iplane]/=n_sf_hits_y[iplane]
    return avg_sf_x, avg_sf_y

def getTimeCorrectedRange(event, scifiDet):
    import rootUtils as ut
    rangePerStation = {1:[], 2:[], 3:[], 4:0, 5:[]}
    avg_sf_x, avg_sf_y = getAvgScifiPos(event, scifiDet)
    hist = {}
    for iplane in range(1, 6):
        a, b = ROOT.TVector3(), ROOT.TVector3()
        ut.bookHist(hist, 'ScifiHitTime_'+str(iplane), "Scifihittime corrected station "+str(iplane), 20, 0, 50)
        for aHit in event.Digi_ScifiHits:
            if not aHit.isValid(): continue
            if not aHit.GetStation()==iplane: continue
            scifiDet.GetSiPMPosition(aHit.GetDetectorID(), a, b)
            L = None
            if aHit.isVertical(): L = b.Y()-avg_sf_y[iplane]
            else: L = avg_sf_x[iplane]-a.X()
            hit_time = scifiDet.GetCorrectedTime(aHit.GetDetectorID(), aHit.GetTime()*TDC2ns, L)
            hist['ScifiHitTime_'+str(iplane)].Fill(hit_time)
        ibin = -1
        ibin = hist['ScifiHitTime_'+str(iplane)].GetMaximumBin()
        rangePerStation[iplane] = [hist['ScifiHitTime_'+str(iplane)].GetBinLowEdge(ibin),  hist['ScifiHitTime_'+str(iplane)].GetBinLowEdge(ibin+3)]
    return rangePerStation

def isInTimeRange(hit_time, time_low, time_up):
    if hit_time > time_low and hit_time <= time_up: return True
    else: return False

def scifiCluster(DigiScifiBranch, scifiDet):
    clusters = []
    hitDict = {}
    clusScifi   = ROOT.TObjArray(100)
    for k in range(DigiScifiBranch.GetEntries()):
        d = DigiScifiBranch[k]
        if not d.isValid(): continue 
        hitDict[d.GetDetectorID()] = k
    hitList = list(hitDict.keys())
    if len(hitList)>0:
            hitList.sort()
            tmp = [ hitList[0] ]
            cprev = hitList[0]
            ncl = 0
            last = len(hitList)-1
            hitvector = ROOT.std.vector("sndScifiHit*")()
            for i in range(len(hitList)):
                if i==0 and len(hitList)>1: continue
                c=hitList[i]
                neighbour = False
                if (c-cprev)==1:    # does not account for neighbours across sipms
                    neighbour = True
                    tmp.append(c)
                if not neighbour  or c==hitList[last]:
                    first = tmp[0]
                    N = len(tmp)
                    hitvector.clear()
                    for aHit in tmp: hitvector.push_back( DigiScifiBranch[hitDict[aHit]])
                    aCluster = ROOT.sndCluster(first,N,hitvector,scifiDet,False)
                    clusters.append(aCluster)
                    if c!=hitList[last]:
                            ncl+=1
                            tmp = [c]
                    elif not neighbour :   # save last channel
                        hitvector.clear()
                        hitvector.push_back( DigiScifiBranch[hitDict[c]])
                        aCluster = ROOT.sndCluster(c,1,hitvector,scifiDet,False)
                        clusters.append(aCluster)
                cprev = c
    clusScifi.Delete()            
    for c in clusters:  
        clusScifi.Add(c)
    return clusScifi

def CorrectScifi(event):
  nsf_statID = {1:0, 2:0, 3:0, 4:0, 5:0}
  nsf_statID_corr = {1:0, 2:0, 3:0, 4:0, 5:0}
  time_plane = {1:0, 2:0, 3:0, 4:0, 5:0}
  Nsf = 0
  Nsf_corr = 0
  hist = {}
  a, b = ROOT.TVector3(), ROOT.TVector3()
  avg_sf_x, avg_sf_y = GetAvgScifiPos(event.Digi_ScifiHits)
  for st in range(1, 6):
    ut.bookHist(hist, 'ScifiHitTime_'+str(st), "Scifihittime corrected station "+str(st), 20, 0, 50)
  for aHit in event.Digi_ScifiHits:
    if not aHit.isValid(): continue
    scifiDet.GetSiPMPosition(aHit.GetDetectorID(), a, b)
    station = aHit.GetStation()
    if aHit.isVertical() :
        L = b.Y() - avg_sf_y[station-1]      
    else :
        L = avg_sf_x[station-1] - a.X()
    Nsf+=1
    nsf_statID[station]+=1
    time_plane[station] = aHit.GetTime()*TDC2ns
    time_plane[station] = scifiDet.GetCorrectedTime(aHit.GetDetectorID(), time_plane[station], L)
    hist['ScifiHitTime_'+str(station)].Fill(time_plane[station])
  for station, sfhit in nsf_statID.items():
    ibin = hist['ScifiHitTime_'+str(station)].GetMaximumBin()
    if sfhit > 40:
      Nsf_corr+= hist['ScifiHitTime_'+str(station)].Integral(ibin, ibin+2)
      nsf_statID_corr[station] = hist['ScifiHitTime_'+str(station)].Integral(ibin, ibin+2)
    else:
      Nsf_corr+=sfhit
      nsf_statID_corr[station] = sfhit
  del hist
  return Nsf_corr, Nsf, nsf_statID_corr, nsf_statID

def getScifiHitDensity(hitcoll, width=1.):
    if len(hitcoll) == 0: return 0
    weights = []
    for i in hitcoll:
        neighbour_no_of_hits = 0
        for j in hitcoll:
            if i == j: continue
            if ROOT.TMath.Abs(i-j) <= width: neighbour_no_of_hits += 1
        weights.append(neighbour_no_of_hits)
    sum_weights = sum(weights)
    if sum_weights: return sum_weights
    else: return 0

def GetVetoBar(detID):
    plane = int((detID/1000)%10)
    bar = int((detID%10000)%1000)
    return plane, bar

def runTracking():
    print('Run tracking')
    tmp = args.inputFile.split('/')
    infile = tmp[len(tmp)-1]
    tmp2 = infile.split('.')
    outfile = '/'.join(tmp[:len(tmp)-1])+'/'+'.'.join(tmp2[:len(tmp2)-1])+'_track.root'
    print('Writing to', outfile)
    os.system('export EOSSHIP=root://eosuser.cern.ch/')
    command= 'python $SNDSW_ROOT/shipLHC/run_muonRecoSND.py -f'+str(args.inputFile)+' -g'+args.geoFile+' -c passing_mu_DS -sc 1 -s '+outfile+' -hf linearSlopeIntercept -o'
    os.system(command)
    return

def bunchXtype(eventTime, runN):
    if runN in FSdict: 
       fsdict = FSdict[runN]      
    else: fsdict = False
    if fsdict:
             bunchNumber = int(eventTime%(4*3564)/4+0.5)
             nb1 = (3564 + bunchNumber - fsdict['phaseShift1'])%3564
             nb2 = (3564 + bunchNumber - fsdict['phaseShift1']- fsdict['phaseShift2'])%3564
             if not "B1" in fsdict: b1 = False
             else: b1 = nb1 in fsdict['B1']
             if not "B2" in fsdict: b2 = False
             else: b2 = nb2 in fsdict['B2']
             IP1 = False
             IP2 = False
             B2noB1 = False
             b1Only = False
             noBeam = False
             if b1:
                IP1 =  fsdict['B1'][nb1]['IP1']
             if b2:
                IP2 =  fsdict['B2'][nb2]['IP2']
             if b2 and not b1:
                B2noB1 = True             
             if b1 and not b2 and not IP1:
                b1Only = True
             if not b1 and not b2: noBeam = True
             return IP1, IP2, b1Only, B2noB1, noBeam
# End snippet

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-f", dest = "inputFile", required = False)
parser.add_argument("-o", dest = "outputFile", required = False)
parser.add_argument("-t", dest = "trackFile", required = False)
parser.add_argument("-g", dest = "geoFile", required = False)
parser.add_argument("--pmu", dest="pmu", required=False, default=False, action='store_true')
parser.add_argument("-npmu", dest="npmu", required=False, default=None, type=int)

args = parser.parse_args()

geofiles    = {'DIS': '/eos/experiment/sndlhc/users/dancc/MuonDIS/ecut1.0_z_2.85_3.55m_Ioni_latelateFLUKA/muonDis_201/1/geofile_full.muonDIS-TGeant4-muonDis_201.root',
               'PMU': '/eos/experiment/sndlhc/MonteCarlo/MuonBackground/muons_down/scoring_1.8_Bfield_4xstat/geofile_full.Ntuple-TGeant4.root',
               'DATA2022': '/eos/experiment/sndlhc/convertedData/physics/2022/geofile_sndlhc_TI18_V5_14August2022.root',
               'DATA2023': '/eos/experiment/sndlhc/convertedData/physics/2023/geofile_sndlhc_TI18_V3_2023.root'}

# Set up TTrees
isMC = False
treeName = "rawConv"
TDC2ns = 1E9/160.316E6
ch = ROOT.TChain(treeName)
ch.Add(args.inputFile)
tag = 'DATA'

if ch.GetEntries() == 0 :
    treeName = "cbmsim"
    isMC = True
    TDC2ns = 1. # In the simulation hit times are always in ns
    del ch
    ch = ROOT.TChain(treeName)
    ch.Add(args.inputFile)
    tag='DIS'
    if args.pmu: 
        tag='PMU'
        if args.npmu: print('N. of samples', int(args.npmu))

if ch.GetEntries() == 0 :
    print("Chain is empty. Exitting")
    exit(-1)

if not args.geoFile:
    if not isMC:
        rindex = args.inputFile.find('run')
        run = args.inputFile[rindex+3:rindex+7]
        year=''
        if int(run) < 5562: year='2022'
        else: year='2023'
        tag=tag+year
    geoFile = geofiles[tag]
    print('Input geoFile not provided, using', geoFile)
else:
    geoFile = args.geoFile

snd_geo = SndlhcGeo.GeoInterface(geoFile)
scifiDet = ROOT.gROOT.GetListOfGlobals().FindObject('Scifi')
muFilterDet = ROOT.gROOT.GetListOfGlobals().FindObject('MuFilter')

ch_tracks = ROOT.TChain(treeName)
ch_tracks.Add(args.trackFile)
ch.AddFriend(ch_tracks)

#ADD INITEVENT

# Set up cuts
cuts = []
################################################################################
# At least 35 SciFi hits
################################################################################
#min_scifi_hits_cut = 35
min_scifi_hits_cut = 60
def min_scifi_hits(event) :
    n_hits = 0
    n_hits_corr = 0
    ret = False
    n_hits_corr, n_hits, nsf_statID_corr, nsf_statID = CorrectScifi(event)
    if isMC:
        if n_hits > min_scifi_hits_cut: ret = True
        return ret, n_hits
    else:
        if n_hits_corr > min_scifi_hits_cut: ret = True
        return ret, n_hits_corr
cuts.append(["More than {0} SciFi hits (CORRECTED)".format(min_scifi_hits_cut), min_scifi_hits, "scifi_nhits", 100, 0, 3000])
################################################################################
# Top and bottom Veto bars are not fired
################################################################################
def VetoBarsCut(event):
    ## add veto bars cut
    VetoBarFired    = {0: list(), 1: list()}
    for aHit in event.Digi_MuFilterHits:
        if aHit.GetSystem()!=1: continue
        plane, bar = GetVetoBar(aHit.GetDetectorID())
        VetoBarFired[plane].append(bar)
    ret = True
    value = 1
    for plane, barlist in VetoBarFired.items():
        for bar in barlist:
            if bar == 0 or bar == 6: 
                ret = False
                value = 0
    
    return ret, value
cuts.append(["Top and Bottom Veto bars not fired", VetoBarsCut, "VetoBars", 2, -0.5, 1.5])
################################################################################
# Single Veto bar per plane fired
################################################################################
def SingleVetoBarsCut(event):
    ## add veto bars cut
    VetoBarFired    = {0: list(), 1: list()}
    for aHit in event.Digi_MuFilterHits:
        if aHit.GetSystem()!=1: continue
        plane, bar = GetVetoBar(aHit.GetDetectorID())
        VetoBarFired[plane].append(bar)
    ret = False
    value = 1
    if len(VetoBarFired[0]) == len(VetoBarFired[1]) == 1:
        if VetoBarFired[1][0] >= VetoBarFired[0][0] -1 and VetoBarFired[1][0] <= VetoBarFired[0][0]+1:
            value = VetoBarFired[0][0]
            ret = True
    else:
        ret = False
    return ret, value
cuts.append(["Ask for single Veto Bars", SingleVetoBarsCut, "SingleVetoBars", 9, -0.5, 8.5])
################################################################################
# Single Cluster in Scifi 1 and 2
################################################################################
uptoiplane = 2
def SingleScifiCluster(event):
    NsfClPl_H       =  {1:0, 2:0, 3:0, 4:0, 5:0}
    NsfClPl_V       =  {1:0, 2:0, 3:0, 4:0, 5:0}
    Nsfcl= 0
    if not isMC:
        DATA_scifiCluster = scifiCluster(event.Digi_ScifiHits, scifiDet)
        for aCl in DATA_scifiCluster:
            detID = aCl.GetFirst()
            Nsfcl+=1
            station = int(detID/1e+6)
            if int(detID/100000)%10 == 1: 
                NsfClPl_V[int(detID/1e+6)]+=1
            else:
                NsfClPl_H[int(detID/1e+6)]+=1
    else:
        for aCl in event.Cluster_Scifi: # MonteCarlo clustering
            detID = aCl.GetFirst()
            Nsfcl+=1
            station = int(detID/1e+6)
            if int(detID/100000)%10 == 1: 
                NsfClPl_V[int(detID/1e+6)]+=1
            else:
                NsfClPl_H[int(detID/1e+6)]+=1
    ret = True
    value = Nsfcl
    for key in NsfClPl_H.keys():
        if key > uptoiplane: continue
        if NsfClPl_H[key] > 1 or NsfClPl_V[key] > 1:
            ret = False
    return ret, value
cuts.append(["Single H/V Cluster in first {0} Scifi planes".format(uptoiplane), SingleScifiCluster, "SingleSFClusters", 20, 0, 20])
################################################################################
# More than N cluster in Scifi 5
################################################################################
#nclus5 = 2
nclus5 = 6
clus_iplane = 5
def ScifiCluster(event):
    NsfClPl_H       =  {1:0, 2:0, 3:0, 4:0, 5:0}
    NsfClPl_V       =  {1:0, 2:0, 3:0, 4:0, 5:0}
    Nsfcl= 0
    if not isMC:
        DATA_scifiCluster = scifiCluster(event.Digi_ScifiHits, scifiDet)
        for aCl in DATA_scifiCluster:
            detID = aCl.GetFirst()
            Nsfcl+=1
            station = int(detID/1e+6)
            if int(detID/100000)%10 == 1: 
                NsfClPl_V[int(detID/1e+6)]+=1
            else:
                NsfClPl_H[int(detID/1e+6)]+=1
    else:
        for aCl in event.Cluster_Scifi: # MonteCarlo clustering
            detID = aCl.GetFirst()
            Nsfcl+=1
            station = int(detID/1e+6)
            if int(detID/100000)%10 == 1: 
                NsfClPl_V[int(detID/1e+6)]+=1
            else:
                NsfClPl_H[int(detID/1e+6)]+=1
    ret = False
    value = NsfClPl_H[clus_iplane]+NsfClPl_V[clus_iplane]
    if NsfClPl_V[clus_iplane] > nclus5 or NsfClPl_H[clus_iplane] > nclus5:
        ret = True
    return ret, value
cuts.append(["More than {0} clusters in SciFi {1}".format(nclus5, clus_iplane), ScifiCluster, "ScifiCluster", 201, -0.5, 200.5])
################################################################################
# Ask for shower developing
################################################################################
clusIncrease_plane = 5
def ScifiClusIncrease(event):
    NsfClPl_H       =  {1:0, 2:0, 3:0, 4:0, 5:0}
    NsfClPl_V       =  {1:0, 2:0, 3:0, 4:0, 5:0}
    NsfClPl         = {1:0, 2:0, 3:0, 4:0, 5:0}
    Nsfcl= 0
    if not isMC:
        DATA_scifiCluster = scifiCluster(event.Digi_ScifiHits, scifiDet)
        for aCl in DATA_scifiCluster:
            detID = aCl.GetFirst()
            Nsfcl+=1
            station = int(detID/1e+6)
            NsfClPl[int(detID/1e+6)]+=1
            if int(detID/100000)%10 == 1: 
                NsfClPl_V[int(detID/1e+6)]+=1
            else:
                NsfClPl_H[int(detID/1e+6)]+=1
    else:
        for aCl in event.Cluster_Scifi: # MonteCarlo clustering
            detID = aCl.GetFirst()
            Nsfcl+=1
            station = int(detID/1e+6)
            NsfClPl[int(detID/1e+6)]+=1
            if int(detID/100000)%10 == 1: 
                NsfClPl_V[int(detID/1e+6)]+=1
            else:
                NsfClPl_H[int(detID/1e+6)]+=1
    ret = False
    value = NsfClPl[clusIncrease_plane]-NsfClPl[clusIncrease_plane-1]
    if NsfClPl[clusIncrease_plane] >= NsfClPl[clusIncrease_plane-1]:
        ret = True
    return ret, value
cuts.append(["Scifi {0} clusters > SciFi {1} clusters".format(clusIncrease_plane, clusIncrease_plane-1), ScifiClusIncrease, "ScifiClusIncrease", 201, -100.5,100.5])
################################################################################
# Shower in US
################################################################################
def US_shower(event):
    nus_statID = {1:0, 2:0, 3:0, 4:0, 5:0}
    ret = False
    for hit in event.Digi_MuFilterHits :
        if hit.GetSystem() != 2 : continue
        if not hit.isValid() : continue
        MuFistation = hit.GetPlane()
        nus_statID[MuFistation+1]+=1
    avg_1 = float((nus_statID[1]+nus_statID[2]+nus_statID[3])/3)
    avg_2 = float((nus_statID[3]+nus_statID[4]+nus_statID[5])/3)
    if avg_1 >= avg_2: ret = True
    return ret, (avg_1-avg_2)
cuts.append(["Avg US1-2-3 bars > Avg US3-4-5 bars", US_shower, "US_avg", 20, -10, 10])
################################################################################
# Ask for at least J bars in US1 & US2
################################################################################
nbars = 2
def min_US1_US2_bars(event):
    nus_statID = {1:0, 2:0, 3:0, 4:0, 5:0}
    ret = False
    US_hits = 0.
    for hit in event.Digi_MuFilterHits :
        if hit.GetSystem() != 2 : continue
        if not hit.isValid() : continue
        US_hits+=1
        MuFistation = hit.GetPlane()
        nus_statID[MuFistation+1]+=1
    if nus_statID[1] >= nbars and nus_statID[2] >= nbars:
        ret = True
    return ret, US_hits
cuts.append(["More than {0} bars fired in US1 and US2 stations".format(nbars), min_US1_US2_bars, "US1_US2_bars", 51, -0.5, 50.5])
################################################################################
# END CUT DEFINITONS
################################################################################
obj = {}
h = {}

if isMC and not args.pmu:
    obj['crsecfile'] = ROOT.TFile('/eos/experiment/sndlhc/MonteCarlo/Pythia6/MuonDIS/muDIScrossSec.root')
    h['g_13'] = obj['crsecfile'].Get('g_13').Clone('g_13')
    h['g_-13'] = obj['crsecfile'].Get('g_-13').Clone('g_-13')

ch.GetEntry(0)
f = ch.GetFile()

# Set up output file
output_file = ROOT.TFile(args.outputFile, "RECREATE")
output_tree = ch.CloneTree(0)

# Copy branch list
branch_list = f.Get("BranchList")
branch_list_copy = branch_list.Clone()
branch_list_copy.Write("BranchList", 1)

# Set up cut flow histogram
cut_flow = f.Get("cutFlow_extended")
cut_flow_extended = ROOT.TH1D(cut_flow.GetName()+"_extended_more", cut_flow.GetTitle(), cut_flow.GetNbinsX()+len(cuts), 0, cut_flow.GetNbinsX()+len(cuts))
cut_flow2 = f.Get("cutFlow2_extended2")
cut_flow_extended2 = ROOT.TH1D(cut_flow2.GetName()+"_extended_more2", cut_flow2.GetTitle(), cut_flow2.GetNbinsX()+len(cuts), 0, cut_flow2.GetNbinsX()+len(cuts))

for i in range(1, cut_flow.GetNbinsX()+1) :
    cut_flow_extended.SetBinContent(i, cut_flow.GetBinContent(i))
    cut_flow_extended.GetXaxis().SetBinLabel(i, cut_flow.GetXaxis().GetBinLabel(i))
    cut_flow_extended2.SetBinContent(i, cut_flow2.GetBinContent(i))
    cut_flow_extended2.GetXaxis().SetBinLabel(i, cut_flow2.GetXaxis().GetBinLabel(i))

for i in range(len(cuts)) :
    cut_flow_extended.GetXaxis().SetBinLabel(i+1+cut_flow.GetNbinsX(), cuts[i][0])
    cut_flow_extended2.GetXaxis().SetBinLabel(i+1+cut_flow2.GetNbinsX(), cuts[i][0])

# Cut-by-cut histograms
cut_by_cut_var_histos = []
for i_cut in range(-1, len(cuts)) :
    this_cut_by_cut_var_histos = []
    for this_cut_name, cut_function, short_name, nbins, range_start, range_end in cuts :
        print("Initializing", short_name, nbins, range_start, range_end)
        this_cut_by_cut_var_histos.append(ROOT.TH1D("_"+str(i_cut+1)+"_"+short_name+"_0",
                                                    short_name,
                                                    nbins, range_start, range_end))
    cut_by_cut_var_histos.append(this_cut_by_cut_var_histos)

if isMC :
    suffix_dict = {}
    suffix_dict[0] = 'muDIS'
    if args.pmu: suffix_dict[0] = 'Pmu'

    cut_by_cut_truth_histos = []
    for i_cut in range(-1, len(cuts)) :
        this_cut_by_cut_truth_histos = []
        this_cut_by_cut_truth_histos.append(ROOT.TH1D(suffix_dict[0]+"_"+str(i_cut+1)+"_Emu", "Emu", 300, 0, 3000))
        this_cut_by_cut_truth_histos.append(ROOT.TH1D(suffix_dict[0]+"_"+str(i_cut+1)+"_vtxX", "vtxX", 200, -100, 0))
        this_cut_by_cut_truth_histos.append(ROOT.TH1D(suffix_dict[0]+"_"+str(i_cut+1)+"_vtxY", "vtxY", 200, 0, 100))
        this_cut_by_cut_truth_histos.append(ROOT.TH1D(suffix_dict[0]+"_"+str(i_cut+1)+"_vtxZ", "vtxZ", 200, 280, 380))
        cut_by_cut_truth_histos.append(this_cut_by_cut_truth_histos)

# N-1
n_minus_1_var_histos = []
for this_cut_name, cut_function, short_name, nbins, range_start, range_end in cuts :
    n_minus_1_var_histos.append(ROOT.TH1D("n_minus_1_"+short_name+"_0", short_name,
                                          nbins, range_start, range_end))

i_pass = 0
passes_cut = [False]*len(cuts)
cut_var = [0.]*len(cuts)
for event in ch:
    n_cuts_passed = 0
    accept_event = True
    WEIGHT=1.
    if isMC:
        if args.pmu: 
            WEIGHT = 8E8/2E8*event.MCTrack[0].GetWeight()*1E5
            if args.npmu: WEIGHT = WEIGHT/int(args.npmu)
        else:
            W = 8E8/2E8*event.MCTrack[0].GetWeight()
            wLHC = W/10/2  #I am using the FLUKA sample twice, mu->p & mu->n
            wInter = event.MCTrack[2].GetWeight()
            PID = event.MCTrack[0].GetPdgCode()
            wDIS = 0.6E-3*h["g_"+str(PID)].Eval(event.MCTrack[0].GetEnergy())
            WEIGHT=wLHC*wInter*wDIS*1E5
    else:
        scifiDet.InitEvent(event.EventHeader)
        muFilterDet.InitEvent(event.EventHeader)
    if not isMC:
        avg_sf_x, avg_sf_y = getAvgScifiPos(event, scifiDet)
        rangePerStation = getTimeCorrectedRange(event, scifiDet)
    for i_cut, cut in enumerate(cuts) :
        #if i_cut != len(cuts)-2: 
        this_cut_passed, this_cut_var = cut[1](event)
        """else: 
            this_cut_passed, this_cut_var, lastdeltas, maxdelta, maxkey = cut[1](event)
            if this_cut_passed:
                print(ch.EventHeader.GetEventNumber(), "Event passed, deltahits", lastdeltas, "maxdelta", maxdelta,"maxkey", maxkey)"""
        passes_cut[i_cut] = this_cut_passed
        if this_cut_passed :
            n_cuts_passed += 1
        cut_var[i_cut] = this_cut_var
        if accept_event and this_cut_passed :
            cut_flow_extended.Fill(cut_flow.GetNbinsX()+i_cut)
            cut_flow_extended2.Fill(cut_flow.GetNbinsX()+i_cut, WEIGHT)
        else :
            accept_event = False
    if accept_event :
        output_tree.Fill()
    
    # Fill histograms
    # Sequential
    for seq_cut in range(-1, len(passes_cut)) :
        if seq_cut >= 0 :
            if not passes_cut[seq_cut] :
                break
#        for i_hist in range(len(cut_by_cut_var_histos[seq_cut+1])) :
        for i_cut_var, this_cut_var in enumerate(cut_var) :
            cut_by_cut_var_histos[seq_cut+1][i_cut_var].Fill(this_cut_var)
        if isMC :

            cut_by_cut_truth_histos[seq_cut+1][0].Fill(event.MCTrack[0].GetEnergy())
            cut_by_cut_truth_histos[seq_cut+1][1].Fill(event.MCTrack[0].GetStartX())
            cut_by_cut_truth_histos[seq_cut+1][2].Fill(event.MCTrack[0].GetStartY())
            cut_by_cut_truth_histos[seq_cut+1][3].Fill(event.MCTrack[0].GetStartZ())
    
    # N-1
    for i_cut in range(len(passes_cut)) :
        if ((not passes_cut[i_cut] and n_cuts_passed == (len(passes_cut)-1))) or (n_cuts_passed == len(passes_cut)) :
            n_minus_1_var_histos[i_cut].Fill(cut_var[i_cut])
    
    if (n_cuts_passed == len(passes_cut)) :
        #print("EVENT {0}".format(i_pass))
        i_pass +=1

cut_flow_extended.Write()
cut_flow_extended2.Write()
output_file.Write()
output_file.Close()

############################################################################
######################## UNUSED CUT DEFINITIONS ############################
############################################################################


################################################################################
# Kill EM showers in upstream SciFi
################################################################################
def NoEMShower(event):
    nsf_stat = {1:0, 2:0, 3:0, 4:0, 5:0}
    nsf_stat_corr = {1:0, 2:0, 3:0, 4:0, 5:0}
    nsf_statID = {1:0, 2:0, 3:0, 4:0, 5:0}
    Nsf_corr, Nsf, nsf_stat_corr, nsf_stat = CorrectScifi(event)
    if isMC: nsf_statID=nsf_stat
    else: nsf_statID=nsf_stat_corr
    FirstSFHit = -1
    deltahits = {2:0, 3:0, 4:0, 5:0}
    FirstSFHit = next((i for i, x in enumerate(nsf_statID.values()) if x), None)+1
    for detID in range(1, 6):
        if detID > FirstSFHit and nsf_statID[detID-1] and detID > 1:
            deltahits[detID] = float((nsf_statID[detID]-nsf_statID[detID-1])/nsf_statID[detID-1])
    ret = False
    delta2132 = [value for key, value in deltahits.items() if key < 4]
    if all(value > -0.05 for value in delta2132): ret = True
    #if all(value for key, value in deltahits.items() if value > -0.05 and key<4): ret = True
    return ret, deltahits[2]
#cuts.append(["Delta2_3_positive", NoEMShower, "Delta2_3_positive", 180, -3, 15])
################################################################################
# Ask for at least 6 US hits
################################################################################
min_US_hits_cut = 8
def min_US_hits(event):
    US_hits = 0.
    ret = True
    for hit in event.Digi_MuFilterHits :
        if hit.GetSystem() != 2 : continue
        if not hit.isValid() : continue
        US_hits += 1.
    if US_hits < min_US_hits_cut :
        ret =  False
    return ret, US_hits
#cuts.append(["Number of US hits is > {0}".format(min_US_hits_cut), min_US_hits, "US_hits", 51, -0.5, 50.5])
################################################################################
# Shower forming
################################################################################
def FollowShower(event):
    nsf_stat = {1:0, 2:0, 3:0, 4:0, 5:0}
    nsf_stat_corr = {1:0, 2:0, 3:0, 4:0, 5:0}
    nsf_statID = {1:0, 2:0, 3:0, 4:0, 5:0}
    Nsf_corr, Nsf, nsf_stat_corr, nsf_stat = CorrectScifi(event)
    if isMC: nsf_statID=nsf_stat
    else: nsf_statID=nsf_stat_corr
    FirstSFHit = -1
    deltahits = {2:0, 3:0, 4:0, 5:0}
    FirstSFHit = next((i for i, x in enumerate(nsf_statID.values()) if x), None)+1
    for detID in range(1, 6):
        if detID > FirstSFHit and nsf_statID[detID-1] and detID > 1:
            deltahits[detID] = float((nsf_statID[detID]-nsf_statID[detID-1])/nsf_statID[detID-1])
    ret = False
    if deltahits[4] > 1. and deltahits[5] > 0: ret = True
    return ret, deltahits[4]
#cuts.append(["Delta4>1 and Delta5>0", FollowShower, "FollowShower", 180, -3, 15])
################################################################################
# Hit density in Scifi
################################################################################
ScifiHitDens_min = 100.
def HitDensityScifi(event):
    vLeft, vRight = ROOT.TVector3(), ROOT.TVector3()
    Scifi_HitCollectionX = {1:[], 2:[], 3:[], 4:[], 5:[]}
    Scifi_HitCollectionY = {1:[], 2:[], 3:[], 4:[], 5:[]}
    ScifiDensityMean     = {1:0, 2:0, 3:0, 4:0, 5:0}
    ret = False
    for aHit in event.Digi_ScifiHits:
        if not aHit.isValid(): continue
        station = aHit.GetStation()
        detID = aHit.GetDetectorID()
        scifiDet.GetSiPMPosition(detID, vLeft, vRight)
        if not isMC:
            L = None
            if aHit.isVertical(): L = vRight.Y()-avg_sf_y[station]
            else: L = avg_sf_x[station]-vLeft.X()
            hit_time = scifiDet.GetCorrectedTime(detID, aHit.GetTime()*TDC2ns, L)
            if not isInTimeRange(hit_time, rangePerStation[station][0], rangePerStation[station][1]): continue
        if aHit.isVertical():
            Scifi_HitCollectionX[station].append(vLeft.X())
        else:
            Scifi_HitCollectionY[station].append(vRight.Y())
    for plane in range(1, 6):
        ScifiDensityMean[plane] = float((getScifiHitDensity(Scifi_HitCollectionX[plane])+getScifiHitDensity(Scifi_HitCollectionY[plane]))/2.)
    if ScifiDensityMean[5] >= ScifiHitDens_min or ScifiDensityMean[4] >= ScifiHitDens_min or ScifiDensityMean[3] >= ScifiHitDens_min: ret = True
    return ret, ScifiDensityMean[5]
#cuts.append(["ScifiHitDensity > {} in at least one of the 3 last planes".format(ScifiHitDens_min), HitDensityScifi, "SFHitDens", 250, 0, 25000])
################################################################################
# More hit variation in Scifi planes
################################################################################
min_maxdelta = 7.35
def MaxDelta(event):
    n_hits_corr, n_hits, nsf_statID_corr, nsf_statID = CorrectScifi(event)
    FirstSFHit = -1
    deltahits = {2:0, 3:0, 4:0, 5:0}
    if not isMC:
        nsf_statID = nsf_statID_corr
    FirstSFHit = next((i for i, x in enumerate(nsf_statID.values()) if x), None)+1
    for detID in range(1, 6):
        if detID > FirstSFHit and nsf_statID[detID-1] and detID > 1:
            deltahits[detID] = float((nsf_statID[detID]-nsf_statID[detID-1])/nsf_statID[detID-1])
    MaxDelta = max([value for key, value in deltahits.items() if key > 2], key=abs)
    MaxKey = [key for key, value in deltahits.items() if value == MaxDelta]
    ret = False
    value = MaxKey[0]
    """if deltahits[MaxKey[0]] < 0: return ret, value  
    lastdeltas = []
    lastdeltas = [value for key, value in deltahits.items() if key > MaxKey[0]]
    if not all(_value > -0.05 for _value in lastdeltas): ret = False
    else: ret = True"""
    if deltahits[value] > min_maxdelta: ret = True
    return ret, value
#cuts.append(["Ask for maximum variation > {}".format(min_maxdelta), MaxDelta, "MinMaxDelta", 5, 0.5, 5.5,])
################################################################################
# Ask for at least J bars in K different planes
################################################################################
nplanes = 3
nus_thresh = 2
def min_US_bars_nplanes(event):
    nus_statID = {1:0, 2:0, 3:0, 4:0, 5:0}
    ret = False
    US_hits = 0.
    for hit in event.Digi_MuFilterHits :
        if hit.GetSystem() != 2 : continue
        if not hit.isValid() : continue
        US_hits+=1
        MuFistation = hit.GetPlane()
        nus_statID[MuFistation+1]+=1
    if sum(1 for x in nus_statID.values() if x >= nus_thresh) >= nplanes:
        ret = True
    return ret, US_hits
#cuts.append(["More than {0} bars fired in at least {1} US planes".format(nus_thresh,nplanes), min_US_bars_nplanes, "USbar_planes", 51, -0.5, 50.5])