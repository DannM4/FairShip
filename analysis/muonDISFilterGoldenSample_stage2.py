import os
# Code snippet from Simona. This should go somewhere where it can be used by several different pieces of code (monitoring, analysis, etc)
import pickle

import numpy as np
import ROOT
import SndlhcGeo

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
    outfile = './'.join(tmp[:len(tmp)-1])+'/'+'.'.join(tmp2[:len(tmp2)-1])+'_track.root'
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
parser.add_argument("--tracking", dest="tracking", required=False, default=False, action='store_true')
parser.add_argument("--pmu", dest="pmu", required=False, default=False, action='store_true')
parser.add_argument("-npmu", dest="npmu", required=False, default=None, type=int)

args = parser.parse_args()

geofiles    = {'DIS': '/eos/experiment/sndlhc/users/dancc/MuonDIS/ecut1.0_z_2.85_3.55m_Ioni_latelateFLUKA/muonDis_201/1/geofile_full.muonDIS-TGeant4-muonDis_201.root',
               'PMU': '/eos/experiment/sndlhc/MonteCarlo/MuonBackground/muons_down/scoring_1.8_Bfield_4xstat/geofile_full.Ntuple-TGeant4.root',
               'DATA2022': '/eos/experiment/sndlhc/convertedData/physics/2022/geofile_sndlhc_TI18_V5_14August2022.root',
               'DATA2023': '/eos/experiment/sndlhc/convertedData/physics/2023/geofile_sndlhc_TI18_V3_2023.root'}

if args.tracking and args.inputFile and args.geoFile:
    runTracking()
    exit()

# Set up TTrees
isMC = False
treeName = "rawConv"

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


snd_geo = SndlhcGeo.GeoInterface(geoFile)
scifiDet = ROOT.gROOT.GetListOfGlobals().FindObject('Scifi')
muFilterDet = ROOT.gROOT.GetListOfGlobals().FindObject('MuFilter')

ch_tracks = ROOT.TChain(treeName)
ch_tracks.Add(args.trackFile)
ch.AddFriend(ch_tracks)

# Set up cuts
cuts = []

# Event in time with IP1 bunch crossing
def eventInBunchCrossing(event) :
    IP1, IP2, b1Only, B2noB1, noBeam = bunchXtype(event.EventHeader.GetEventTime(), event.EventHeader.GetRunId())
    return IP1
#cuts.append(["Event in time with IP1 collision", eventInBunchCrossing])

################################################################################
# Event direction cut
################################################################################
def direction(event) :
    if not isMC :
        scifiDet.InitEvent(event.EventHeader)
        muFilterDet.InitEvent(event.EventHeader)

    t_scifi = []
    t_muon = []

    for hit in event.Digi_ScifiHits :
        if not hit.isValid() : 
            continue

        t = hit.GetTime()*TDC2ns
#        if not isMC :
#            t = scifiDet.GetCorrectedTime(hit.GetDetectorID(), t, 0)
        
        t_scifi.append(t)

    for hit in event.Digi_MuFilterHits :
        if hit.GetSystem() != 3 :
            continue
        if not hit.isValid() :
            continue

        if hit.isVertical() :
            t = hit.GetTime(0)*TDC2ns
#            if not isMC :
#                t = muFilterDet.GetCorrectedTime(hit.GetDetectorID(), 0, t, 0)
        else :
            t = [hit.GetTime(i)*TDC2ns for i in range(2)]

 #           if not isMC :
 #               for i in range(2) :
 #                   t[i] = muFilterDet.GetCorrectedTime(hit.GetDetectorID(), i, t[i], 0)
            t = np.mean(t)
            t_muon.append(t)

    delta_t = np.sort(t_scifi)[0] - np.sort(t_muon)[-1]
    ret = delta_t < 0

    return ret, delta_t
#cuts.append(["Event direction", direction])

################################################################################
# Event has one reconstructed DS track
################################################################################
def eventHasOneTrack(event) :
    ret = False
    n_tracks = event.Reco_MuonTracks.GetEntries()
    if n_tracks >= 1 :
        ret =  True
    return ret, n_tracks 
cuts.append(["Event has one reconstructed DS track", eventHasOneTrack, "n_DS_tracks", 3, 0, 3])

################################################################################
# Track intercepts first layer < 5 cm from shower center
################################################################################
def trackInterceptShower(event, thresh) : 
    n_ver = [0]*5
    n_hor = [0]*5
    x_sta = [0.]*5
    y_sta = [0.]*5

    a = ROOT.TVector3()
    b = ROOT.TVector3()

    for hit in event.Digi_ScifiHits :
        if not hit.isValid() :
            continue
            
        scifiDet.GetSiPMPosition(hit.GetDetectorID(), a, b)

        if hit.isVertical() :
            n_ver[hit.GetStation()-1] += 1
            x_sta[hit.GetStation()-1] += (a.X() + b.X())/2.
        else :
            n_hor[hit.GetStation()-1] += 1
            y_sta[hit.GetStation()-1] += (a.Y() + b.Y())/2.
    
    fracsum = np.cumsum(np.add(n_ver, n_hor)/(np.sum(n_ver)+np.sum(n_hor)))
    station = next(x[0] for x in enumerate(fracsum) if x[1] > thresh)
    
    x_first_sta = None
    y_first_sta = None

    for i in range(station, 5) :
        if n_ver[i] >= 2 :
            x_first_sta = x_sta[i]/n_ver[i]
            break

    for i in range(station, 5) :
        if n_hor[i] >= 2 :
            y_first_sta = y_sta[i]/n_hor[i]
            break
    
    track = event.Reco_MuonTracks.At(0)

    track_start = track.getStart()
    track_stop = track.getStop()

    slope_X = (track_stop.X() - track_start.X())/(track_stop.Z() - track_start.Z())
    slope_Y = (track_stop.Y() - track_start.Y())/(track_stop.Z() - track_start.Z())
    
    station_z = 298.97+13*station
        
    station_vertex = ROOT.TVector3(track_start.X() + slope_X*(station_z - track_start.Z()), track_start.Y() + slope_Y*(station_z - track_start.Y()), station_z)
    try :
        if ((station_vertex.X() - x_first_sta)**2 + (station_vertex.Y() - y_first_sta)**2)**0.5 < 5 :
            return True
        else :
            return False
    except :
        print(n_ver)
        print(n_hor)
        exit()
#cuts.append(["Track intercepts first layer < 5 cm from shower center", trackInterceptShower])

################################################################################
# Track intercepts first SciFi plane within 5 cm of the edge
################################################################################
d_scifi_fiducial = 5
def trackInScifiFiducial(event) :
    if len(event.Reco_MuonTracks) < 1 :
        return False, -999
    track = event.Reco_MuonTracks.At(0)

    track_start = track.getStart()
    track_stop = track.getStop()

    slope_X = (track_stop.X() - track_start.X())/(track_stop.Z() - track_start.Z())
    slope_Y = (track_stop.Y() - track_start.Y())/(track_stop.Z() - track_start.Z())

    z_first_scifi = 300.

    track_first_scifi_extrap = ROOT.TVector3(track_start.X() + slope_X*(z_first_scifi - track_start.Z()), track_start.Y() + slope_Y*(z_first_scifi - track_start.Z()), z_first_scifi)

    min_d = 1e6
    if -15.5 + track_first_scifi_extrap.Y() < min_d :
        min_d = -15.5 + track_first_scifi_extrap.Y() < min_d
    if + 15.5+39 - track_first_scifi_extrap.Y()  < min_d :
        min_d = + 15.5+39 - track_first_scifi_extrap.Y()
    if -8 - track_first_scifi_extrap.X() < min_d:
        min_d = -8  - track_first_scifi_extrap.X()
    if +8 +39 + track_first_scifi_extrap.X() < min_d :
        min_d = +8 +39 + track_first_scifi_extrap.X()
    ret = True
    if track_first_scifi_extrap.Y() < 15.5 + d_scifi_fiducial :
        ret =  False
    if track_first_scifi_extrap.Y() > 15.5+39 - d_scifi_fiducial :
        ret = False
    if track_first_scifi_extrap.X() > -8 - d_scifi_fiducial :
        ret = False
    if track_first_scifi_extrap.X() < -8 -39 + d_scifi_fiducial:
        ret = False
    return ret, min_d
cuts.append(["Track intercepts first SciFi plane < {0} cm from edge".format(d_scifi_fiducial), trackInScifiFiducial, "ds_track_extrap_fiducial", 100, -50, 50])

################################################################################
# SciFi hit to DS track DOCA per plane < 3 cm for both projections
################################################################################
sum_min_dca_cut = 3
def sum_min_dca(event) :
    if len(event.Reco_MuonTracks) < 1 :
        return False, 999
    
    track = event.Reco_MuonTracks.At(0)

    track_start = track.getStart()
    track_stop = track.getStop()

    slope_X = (track_stop.X() - track_start.X())/(track_stop.Z() - track_start.Z())
    slope_Y = (track_stop.Y() - track_start.Y())/(track_stop.Z() - track_start.Z())

    a = ROOT.TVector3()
    b = ROOT.TVector3()

    min_dca_ver = [1e10]*5
    min_dca_hor = [1e10]*5
    
    for hit in event.Digi_ScifiHits :
        if not hit.isValid() :
            continue
            
        scifiDet.GetSiPMPosition(hit.GetDetectorID(), a, b)

        if hit.isVertical() :
            slope_XY = slope_X
            track_start_XY = track_start.X()
            track_stop_XY = track_stop.X()
            hit_pos = [(a.Z() + b.Z())/2., (a.X() + b.X())/2.]
        else :
            slope_XY = slope_Y
            track_start_XY = track_start.Y()
            track_stop_XY = track_stop.Y()
            hit_pos = [(a.Z() + b.Z())/2., (a.Y() + b.Y())/2.]

        if hit.isVertical() :
            this_d = np.abs(track_start_XY + slope_XY*(hit_pos[0]-track_start.Z()) - hit_pos[1])
            if this_d < min_dca_ver[hit.GetStation()-1] :
                min_dca_ver[hit.GetStation()-1] = this_d
        else :
            this_d = np.abs(track_start_XY + slope_XY*(hit_pos[0]-track_start.Z()) - hit_pos[1])
            if this_d < min_dca_hor[hit.GetStation()-1] :
                min_dca_hor[hit.GetStation()-1] = this_d

    min_dca_ver = np.array(min_dca_ver)
    min_dca_hor = np.array(min_dca_hor)
    
    n_ver = 5 - np.sum(min_dca_ver > 9.9e9)
    n_hor = 5 - np.sum(min_dca_hor > 9.9e9)

    min_dca_ver[min_dca_ver > 9.9e9] = 0
    min_dca_hor[min_dca_hor > 9.9e9] = 0

    avg_ver = np.sum(min_dca_ver)/n_ver
    avg_hor = np.sum(min_dca_hor)/n_hor

    ret_value = np.max([avg_ver, avg_hor])
    ret = ret_value <= sum_min_dca_cut

    return ret, ret_value
cuts.append(["Sum of min DOCA per station < {0} cm".format(sum_min_dca_cut), sum_min_dca, "doca", 30, 0, 30])

################################################################################
# At least 35 SciFi hits
################################################################################
min_scifi_hits_cut = 35
def min_scifi_hits(event) :
    n_hits = 0
    ret = False
    for hit in event.Digi_ScifiHits :
        if hit.isValid() and hit.GetStation!=1:
            n_hits += 1
            if n_hits > 35 :
                ret = True
    return ret, n_hits
#cuts.append(["More than {0} SciFi hits".format(min_scifi_hits_cut), min_scifi_hits, "scifi_nhits", 100, 0, 3000])

################################################################################
# Min QDC
################################################################################
min_QDC_data = 600
min_QDC_MC = 700
def min_US_QDC(event) :
    US_QDC = 0
    ret = False
    for hit in event.Digi_MuFilterHits :
        if hit.GetSystem() != 2 :
            continue
        if not hit.isValid() :
            continue
        for key, value in hit.GetAllSignals() :
            US_QDC += value
            if isMC and (US_QDC > min_QDC_MC) :
                ret = True
            if (not isMC) and (US_QDC > min_QDC_data) :
                ret = True
    if isMC :
        return ret, US_QDC - 100
    else :
        return ret, US_QDC
cuts.append(["US QDC larger than {0} ({1}) for data (MC)".format(min_QDC_data, min_QDC_MC), min_US_QDC, "US_QDC", 100, 0, 4000])

################################################################################
# Max DS activity < 10
################################################################################
max_DS_hits_cut = 10
def max_DS_hits(event) :
    DS_hits = 0.
    ret = True
    for hit in event.Digi_MuFilterHits :
        if hit.GetSystem() != 3 :
            continue
        if not hit.isValid() :
            continue
        
        if hit.GetPlane() == 3 :
            DS_hits += 1.
        else :
            DS_hits += 0.5
    if DS_hits > max_DS_hits_cut :
        ret =  False
    return ret, DS_hits
cuts.append(["Number of DS hits per projection is < {0}".format(max_DS_hits_cut), max_DS_hits, "DS_hits", 100, 0, 100])

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
cut_flow = f.Get("cutFlow")
cut_flow_extended = ROOT.TH1D(cut_flow.GetName()+"_extended", cut_flow.GetTitle(), cut_flow.GetNbinsX()+len(cuts), 0, cut_flow.GetNbinsX()+len(cuts))
cut_flow2 = f.Get("cutFlow2")
cut_flow_extended2 = ROOT.TH1D(cut_flow2.GetName()+"_extended2", cut_flow2.GetTitle(), cut_flow2.GetNbinsX()+len(cuts), 0, cut_flow2.GetNbinsX()+len(cuts))

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
for i_event, event in enumerate(ch):
    n_cuts_passed = 0
    accept_event = True
    WEIGHT=1.
    if i_event%10000 == 0: print("Sanity check, current event ", i_event)
    if isMC:
        if args.pmu: 
            WEIGHT = 8E8/2E8*event.MCTrack[0].GetWeight()*1E5
            if args.npmu: WEIGHT = WEIGHT/int(args.npmu)
        else:
            W = 8E8/2E8*event.MCTrack[0].GetWeight()
            wLHC = W/10/2. # I am using the same FLUKA sample twice, mu->p & mu->n
            wInter = event.MCTrack[2].GetWeight()
            PID = event.MCTrack[0].GetPdgCode()
            wDIS = 0.6E-3*h["g_"+str(PID)].Eval(event.MCTrack[0].GetEnergy())
            WEIGHT=wLHC*wInter*wDIS*1E5
    else:
        scifiDet.InitEvent(event.EventHeader)
        muFilterDet.InitEvent(event.EventHeader)
    for i_cut, cut in enumerate(cuts) : 
        this_cut_passed, this_cut_var = cut[1](event)
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
        if not i_pass %10000: print("EVENT {0}".format(i_pass))
        i_pass +=1

cut_flow_extended.Write()
cut_flow_extended2.Write()
output_file.Write()
output_file.Close()
