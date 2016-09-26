#!/usr/bin/env python

############ INIT ############

from ROOT import *
from math import *
from array import *

### Analysis settings
debug       = 1
escale      = true
save_plots  = false

Edec=[4,5,6,7]             # Energy decade (MeV)
NEbins=[5,5,5]             # Number of bins between 10**Edec

AltRef      = 550.           # Reference altitude in Kamae Limb model
LimbPeak    = 67.798         # Corresponding angle of Limb peak
thShiftFact = 0.0211

ThetaCut  = 70.            # Selection on Theta_LAT for Limb
ZenithCut = 100.           # Selection on Theta_nadir for Limb
Cut="FT1Theta<%f&&FT1ZenithTheta>%f"%(ThetaCut,ZenithCut)

phbw  = [72,  0., 360.]    # Phi   [N bins, phi1, phi2]
thbw  = [400, 0., 80.]     # Theta [N bins, theta1, theta2]
nDraw = [0,4,8,12]         # Plot count maps for these i^th bin of Ebins

### Root settings
gROOT.SetStyle('Plain')
gStyle.SetPalette(1)
gStyle.SetOptStat(10)
gStyle.SetErrorX(0)

### Data sets
datV='P301v1'              # Data set version
irfV='P8V6ULTRACLEANVETO'  # IRF version
inFile='SkimMerge.root'

### Output
ofileroot="CntMap%s_%stest963_thetashiftpeak67.798-2.root"%(irfV,datV)
fout=TFile(ofileroot,"new")


Ebins=[]
for i in range(len(Edec)-1):
    for j in range(NEbins[i]):
        Ebins.append(10**(Edec[i]+float(j)/NEbins[i]))
Ebins.append(10**Edec[-1])
EbinsArray=array('d',Ebins)

inM=TChain("MeritTuple")
inM.Add(inFile)

cntSpec=TH1F("cntSpec","Raw Count",len(Ebins)-1,EbinsArray)
altDist=TH1F("altDist","LAT altitude",100, 525., 575.)

c1=TCanvas("c1")
c1.SetLogx()
c1.SetLogy()
inM.Draw("0.963*FT1Energy>>cntSpec",Cut) #test E scale 0.963
cntSpec.Draw("E1same")
c1.SaveAs("CntMap_rawspec963_thetashift.gif")
cntSpec.Write()

c2=TCanvas("c2")
inM.Draw("PtAlt>>altDist",Cut)
altDist.Draw()
c2.SaveAs("altDist_rawspec963_thetashift.gif")
altDist.Write()

cntMap=[]
corrTh=[]
cntMapNocorr=[]
cntDist=[]
cntDistNocorr=[]
shiftMap=[]
shiftMapNocorr=[]
for n in range(len(Ebins)-1):
        Etitle="%.1f-%.1f GeV"%(Ebins[n]/1000.,Ebins[n+1]/1000.)
        cntMapName="cntmap%d"%(n)
        cntMaptmp=TH2F(cntMapName,Etitle,phbw[0],phbw[1],phbw[2],thbw[0],thbw[1],thbw[2])
        cntMap.append(cntMaptmp)
        cntMapNocorrName="cntmapnocorr%d"%(n)
        cntMapNocorrtmp=TH2F(cntMapNocorrName,Etitle,phbw[0],phbw[1],phbw[2],thbw[0],thbw[1],thbw[2])
        cntMapNocorr.append(cntMapNocorrtmp)
        cntDistName="cntdist%d"%(n)
        cntDisttmp=TH2F(cntDistName,Etitle,phbw[0],phbw[1],phbw[2],thbw[0],-5.,5.)
        cntDist.append(cntDisttmp)
        cntDistNocorrName="cntdistnocorr%d"%(n)
        cntDistNocorrtmp=TH2F(cntDistNocorrName,Etitle,phbw[0],phbw[1],phbw[2],thbw[0],-5.,5.)
        cntDistNocorr.append(cntDistNocorrtmp)
        corrThName="corrth%d"%(n)
        corrThtmp=TH2F(corrThName,Etitle,100,63.,80.,100,63.,80.)
        corrTh.append(corrThtmp)
        shiftMapName="shiftmap%d"%(n)
        shiftMaptmp=TH2F(shiftMapName,Etitle,100, 525., 575.,100,63.,80.)
        shiftMap.append(shiftMaptmp)
        shiftMapNocorrName="shiftmapnocorr%d"%(n)
        shiftMapNocorrtmp=TH2F(shiftMapNocorrName,Etitle,100, 525., 575.,100,63.,80.)
        shiftMapNocorr.append(shiftMapNocorrtmp)

for i in range(inM.GetEntries()):
#for i in range(100):
    if (not(i%1000)): print "Event ", i
#    print "Event ", i
    inM.GetEntry(i)
    DeltaTheta = thShiftFact*(AltRef - getattr(inM,"PtAlt"))
    if (not(i%1000)):
        print "  PtAlt =", getattr(inM,"PtAlt")," DeltaTheta =", DeltaTheta
        print "  Ft1Theta = ", getattr(inM,"FT1Theta"), " FT1ZenithTheta = ", getattr(inM,"FT1ZenithTheta"), " corrected FT1ZenithTheta = ", getattr(inM,"FT1ZenithTheta")+DeltaTheta
#   if getattr(inM,"FT1Theta")<ThetaCut and getattr(inM,"FT1ZenithTheta")>ZenithCut:   # va corretto?... eventualmente ==> getattr(inM,"FT1ZenithTheta")+DeltaTheta>ZenithCut
    if getattr(inM,"FT1Theta")<ThetaCut and (getattr(inM,"FT1ZenithTheta")+DeltaTheta)>ZenithCut:
        for j in range(len(Ebins)-1):
#           if Ebins[j]<=getattr(inM,"FT1Energy") and getattr(inM,"FT1Energy")<Ebins[j+1]:                        # no E scale
            if Ebins[j]<=0.963*getattr(inM,"FT1Energy") and 0.963*getattr(inM,"FT1Energy")<Ebins[j+1]:            # test E scale 0.963
                cntMapNocorr[j].Fill(getattr(inM,"FT1EarthAzimuth"),180.-getattr(inM,"FT1ZenithTheta"))           # no Theta shift
                cntMap[j].Fill(getattr(inM,"FT1EarthAzimuth"),180.-getattr(inM,"FT1ZenithTheta")-DeltaTheta)      # Theta -> Theta-DeltaTheta
                #print "  Filling cntmap Ebin ", Ebins[j]," with ",getattr(inM,"FT1EarthAzimuth")," - ",180.-getattr(inM,"FT1ZenithTheta")-DeltaTheta
                cntDistNocorr[j].Fill(getattr(inM,"FT1EarthAzimuth"),LimbPeak-(180.-getattr(inM,"FT1ZenithTheta")))           # no Theta shift
                cntDist[j].Fill(getattr(inM,"FT1EarthAzimuth"),LimbPeak-(180.-getattr(inM,"FT1ZenithTheta")-DeltaTheta))      # Theta -> Theta-DeltaTheta
                #print "  Filling cntdist Ebin ", Ebins[j]," with ",getattr(inM,"FT1EarthAzimuth")," - ",68.-(180.-getattr(inM,"FT1ZenithTheta")-DeltaTheta)
                shiftMapNocorr[j].Fill(getattr(inM,"PtAlt"),180.-getattr(inM,"FT1ZenithTheta"))           # no Theta shift
                shiftMap[j].Fill(getattr(inM,"PtAlt"),180.-getattr(inM,"FT1ZenithTheta")-DeltaTheta)      # Theta -> Theta-DeltaTheta
                corrTh[j].Fill(180.-getattr(inM,"FT1ZenithTheta"),180.-getattr(inM,"FT1ZenithTheta")-DeltaTheta)
                
                break

for i in range(len(Ebins)-1):
    cntMap[i].Write()
    cntMapNocorr[i].Write()
    cntDist[i].Write()
    cntDistNocorr[i].Write()
    shiftMap[i].Write()
    shiftMapNocorr[i].Write()
    corrTh[i].Write()

"""
Earth=TEllipse(0,0,.485,.485)
Earth.SetFillStyle(0)
Earth.SetLineWidth(1)
Earth.SetLineStyle(2)
c2=TCanvas("c2","c2",300*len(nDraw),300)
c2.Divide(len(nDraw),1)
nd=1
for i in nDraw:
    c2.cd(nd)
    c2.cd(nd).SetLogz()
    gPad.SetTheta(-90)
    gPad.SetPhi(-90)
    cntMap[i].Draw("SURF2POLZ")
    Earth.Draw()
    c2.cd(nd).Update()
    nd+=1
"""

fout.Close()

print "Done!"

raw_input()
