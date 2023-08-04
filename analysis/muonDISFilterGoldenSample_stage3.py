import ROOT
import numpy as np
import os
import SndlhcGeo

# Code snippet from Simona. This should go somewhere where it can be used by several different pieces of code (monitoring, analysis, etc)
import pickle
p = open("/eos/experiment/sndlhc/convertedData/commissioning/TI18/FSdict.pkl",'rb')
FSdict = pickle.load(p)

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

args = parser.parse_args()

    
snd_geo = SndlhcGeo.GeoInterface(args.geoFile)
scifiDet = ROOT.gROOT.GetListOfGlobals().FindObject('Scifi')
muFilterDet = ROOT.gROOT.GetListOfGlobals().FindObject('MuFilter')

# Set up TTrees
isMC = False
treeName = "rawConv"

ch = ROOT.TChain(treeName)
ch.Add(args.inputFile)

if ch.GetEntries() == 0 :
    treeName = "cbmsim"
    isMC = True
    TDC2ns = 1. # In the simulation hit times are always in ns
    del ch
    ch = ROOT.TChain(treeName)
    ch.Add(args.inputFile)

if ch.GetEntries() == 0 :
    print("Chain is empty. Exitting")
    exit(-1)

ch_tracks = ROOT.TChain(treeName)
ch_tracks.Add(args.trackFile)
ch.AddFriend(ch_tracks)

# Set up cuts
cuts = []

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
cuts.append(["Top and Bottom Veto bars fired", VetoBarsCut, "VetoBars", 2, -0.5, 1.5])
################################################################################
# Ask for shower developing
################################################################################
"""def MaxDelta(event):
    nsf_statID = {1:0, 2:0, 3:0, 4:0, 5:0}
    for aHit in event.Digi_ScifiHits:
        if not aHit.isValid(): continue
        station = aHit.GetStation()
        nsf_statID[station]+=1
    FirstSFHit = -1
    deltahits = {2:0, 3:0, 4:0, 5:0}
    FirstSFHit = next((i for i, x in enumerate(nsf_statID.values()) if x), None)+1
    for detID in range(1, 6):
        if detID > FirstSFHit and nsf_statID[detID-1] and detID > 1:
            deltahits[detID] = float((nsf_statID[detID]-nsf_statID[detID-1])/nsf_statID[detID-1])
    MaxDelta = max([value for key, value in deltahits.items() if key > 2], key=abs)
    MaxKey = [key for key, value in deltahits.items() if value == MaxDelta]
    ret = False
    value = MaxKey[0]
    if deltahits[MaxKey[0]] < 0: return ret, value, deltahits, MaxDelta, MaxKey  
    lastdeltas = []
    lastdeltas = [value for key, value in deltahits.items() if key > MaxKey[0]]
    if not all(_value > -0.05 for _value in lastdeltas): ret = False
    else: ret = True
    return ret, value, deltahits, MaxDelta, MaxKey
cuts.append(["Ask for HAD shower develop", MaxDelta, "HShowerDev", 5, 0.5, 5.5,])"""
################################################################################
# Ask for at least 6 US hits
################################################################################
min_US_hits_cut = 7
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
cuts.append(["Number of US hits is > {0}".format(min_US_hits_cut), min_US_hits, "US_hits", 51, -0.5, 50.5])
################################################################################
# Shower forming
################################################################################
def FollowShower(event):
    nsf_statID = {1:0, 2:0, 3:0, 4:0, 5:0}
    for aHit in event.Digi_ScifiHits:
        if not aHit.isValid(): continue
        station = aHit.GetStation()
        nsf_statID[station]+=1
    FirstSFHit = -1
    deltahits = {2:0, 3:0, 4:0, 5:0}
    FirstSFHit = next((i for i, x in enumerate(nsf_statID.values()) if x), None)+1
    for detID in range(1, 6):
        if detID > FirstSFHit and nsf_statID[detID-1] and detID > 1:
            deltahits[detID] = float((nsf_statID[detID]-nsf_statID[detID-1])/nsf_statID[detID-1])
    ret = False
    if deltahits[4] > 1. and deltahits[5] > 0: ret = True
    return ret, deltahits[4]
cuts.append(["Delta4>1 and Delta5>0", FollowShower, "FollowShower", 180, -3, 15])
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
    if nus_statID[1] > 2. and nus_statID[2] > 2 and nus_statID[3] > 2: ret = True
    return ret, nus_statID[1]
cuts.append(["US1-2-3 hits > 2", US_shower, "ThreeUS_2", 11, -0.5, 10.5])
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
        if args.pmu: WEIGHT = 8E8/2E8*event.MCTrack[0].GetWeight()*1E5
        else:
            W = 8E8/2E8*event.MCTrack[0].GetWeight()
            wLHC = W/10/2  #I am using the FLUKA sample twice, mu->p & mu->n
            wInter = event.MCTrack[2].GetWeight()
            PID = event.MCTrack[0].GetPdgCode()
            wDIS = 0.6E-3*h["g_"+str(PID)].Eval(event.MCTrack[0].GetEnergy())
            WEIGHT=wLHC*wInter*wDIS*1E5
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
