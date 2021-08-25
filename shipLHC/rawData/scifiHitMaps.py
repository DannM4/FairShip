import ROOT,os
import rootUtils as ut
h={}
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-r", "--runNumber", dest="runNumber", help="run number", type=int,required=True)
options = parser.parse_args()

f=ROOT.TFile('sndsw_raw_'+str(options.runNumber).zfill(6)+'.root')

def xPos(detID):
        nStation = detID//100000
        chan = detID - nStation*100000
        mat   = chan//512
        X = chan-mat*512
        return [nStation,mat,X]   # even numbers are Y, odd numbers X

def makeDisplayHit(x,y,col):
      dHit = ROOT.TGraph()
      dHit.SetPoint(0,x-5,y)
      dHit.SetPoint(1,x+5,y)
      dHit.SetPoint(2,x,y)
      dHit.SetPoint(3,x,y+0.1)
      dHit.SetPoint(4,x,y-0.1)
      dHit.SetLineColor(col)
      return dHit

def loopEvents(start=0,save=False):
 ut.bookHist(h,'xy',' ',200,-0.5,3*512.,10,-0.5,5.5)
 ut.bookCanvas(h,'event','event',1024,768,1,2)
 h['xy'].SetStats(0)

 N = -1
 Tprev = -1
 for sTree in f.rawConv:
    color={0:ROOT.kBlue,1:ROOT.kGreen}
    ptext={0:'   Y projection',1:'   X projection'}
    text = ROOT.TLatex()
    N+=1
    if N<start: continue
    T = sTree.EventHeader.GetEventTime()
    dT = 0
    if Tprev >0: dT = T-Tprev
    print('event ',N,':',dT/160.E6)
    Tprev = T
    h['hits0']=[]
    h['hits1']=[]
    for aHit in sTree.Digi_SciFiHits:
        X =  xPos(aHit.GetDetectorID())
        proj = X[0]%2
        h['hits'+str(proj)].append(makeDisplayHit(X[2]+512*X[1],int(X[0]/2+0.5),color[proj]))
        print('hit:',X)
    for p in range(2):
       tc = h['event'].cd(p+1)
       h['xy'].Draw()
       for x in h['hits'+str(p)]: x.Draw('same')
       text.DrawLatex(20,5.8,'event '+str(N)+ptext[p])
    h['event'].Update()
    if save: h['event'].Print('event_'+"{:04d}".format(N)+'.png')
    rc = input("hit return for next event: ")
    print(rc)
    if rc=='q': break
 if save: os.system("convert -delay 60 -loop 0 *.png animated.gif")

def hitMaps():
 for mat in range(30):
    ut.bookHist(h,'mat_'+str(mat),'hit map / mat',512,-0.5,511.5)
    ut.bookHist(h,'sig_'+str(mat),'signal / mat',100,0.0,10.)
 N=-1
 for sTree in f.rawConv:
    N+=1
    for aHit in sTree.Digi_SciFiHits:
        X =  xPos(aHit.GetDetectorID())
        rc = h['mat_'+str(X[0]*3+X[1])].Fill(X[2])
        rc  = h['sig_'+str(X[0]*3+X[1])].Fill(aHit.GetSignal(0))
 ut.bookCanvas(h,'hitmaps',' ',1024,768,6,5)
 ut.bookCanvas(h,'signal',' ',1024,768,6,5)
 for mat in range(30):
    tc = h['hitmaps'].cd(mat+1)
    A = h['mat_'+str(mat)].GetSumOfWeights()/512.
    if h['mat_'+str(mat)].GetMaximum()>10*A: h['mat_'+str(mat)].SetMaximum(10*A)
    h['mat_'+str(mat)].Draw()
    tc = h['signal'].cd(mat+1)
    h['sig_'+str(mat)].Draw()
def eventTime():
 Tprev = -1
 ut.bookHist(h,'Etime','delta event time',100,0.0,1.)
 ut.bookCanvas(h,'E',' ',1024,768,1,1)
 for sTree in f.rawConv:
    T = sTree.EventHeader.GetEventTime()
    dT = 0
    if Tprev >0: dT = T-Tprev
    Tprev = T
    rc = h['Etime'].Fill(dT/160.E6)
 tc = h['E'].cd()
 rc = h['Etime'].Fit('expo')
 stats = h['Etime'].FindObject("stats")
 stats.SetOptFit(111)



