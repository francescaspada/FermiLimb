#!/usr/bin/env python

############ INIT ############

import numpy
import sys
from ROOT import *
from math import *
from array import *
from copy import *

### Analysis settings
debug = 1

escale    = true
subtract_bg = true
save_plots = true

Edec=[1,2,3,4]         # Energy decade (GeV) [DATA]
NEbins=[5,5,5]         # Number of bins between 10**Edec [DATA]
thbw=[400,0.,80.]      # Theta [N bins, theta1, theta2] [DATA]
#thrn=[66.,70.]         # Spectrum for this theta range [DATA]
thrn=[68.4,70.]       # Spectrum for this theta range [DATA] thin target
thBGrn=[77.,80.]       # Diff. BG from this theta range [DATA]
nDraw=[0,3,6,9,10]     # Plot maps for these i^th bin of Ebins [DATA]

### Constants
if escale: 
    NormConst=8.954e-5
else:
    NormConst=5.500e-5      # range 66-70 coincidence ad minchiam

SpecIndex = 2.75            # spectral index


### Root settings
gROOT.SetStyle('Plain')
gStyle.SetPalette(1)
gStyle.SetErrorX(0)
gStyle.SetOptFit(0)
gStyle.SetOptStat(0)

kamaeMarker = 24
limbMarker  = 25
otherMarker = 24

myMarkerColor  =TColor.GetColor("#000000")
myLightGreen   =TColor.GetColor("#55CCCC")
myDarkGreen    =TColor.GetColor("#338888")
myLightBlue    =TColor.GetColor("#55AADD")
myDarkBlue     =TColor.GetColor("#3366AA")
myLightMagenta =TColor.GetColor("#FF99DD")
myDarkMagenta  =TColor.GetColor("#AA4488")
myLightRosso   =TColor.GetColor("#FF9999")
myDarkRosso    =TColor.GetColor("#AA5555")

### Data sets
datV='P301v1'
irfV='P8V6ULTRACLEANVETO' #IRF version
if escale:
    #cmFile=TFile('CntMap'+irfV+'_'+datV+'test963SLAC.root')  # WITH E*0.963-SCALE
    cmFile=TFile('CntMap'+irfV+'_'+datV+'test963_thetashiftpeak67.798-2.root')  # WITH E*0.963-SCALE
    escFact='0.963ESCALE'
else:
    cmFile=TFile('CntMap'+irfV+'_'+datV+'.root')             # NO E-SCALE
    escFact='NOESCALE'
emFile=TFile('ExpMap'+irfV+'_'+datV+'.root')
kaFile=open("gammaspectrum.txt")
kaFileH=open("gammaspectrumHe.txt")
kaFileM=open("gammaspectrumM.txt")
kaFile10p=open("gammaspectrumEp.txt")
kaFile10m=open("gammaspectrumEm.txt")

### Theta bin width, bins in 0.2 deg, used bins
thw=(thbw[2]-thbw[1])/thbw[0]  
thnb=int(round(thbw[0]/(thbw[2]-thbw[1])))  ## how many bins for 1 degree

### Phi bin width 5.0 deg
phw = 5.0

### n. of BG bins
thBGbr=[int(round(thBGrn[0]*thnb+1)),int(round(thBGrn[1]*thnb))]
nthBGb=thBGbr[1]-thBGbr[0]


############ KAMAE ############

### Create energy bins for Kamae flux and get flux from AMS parameters
if debug:
    print "===== KAMAE ====="
Kbins=[]
KbinsErr=[]
KFlux=[]
KFluxH=[]
KFluxM=[]
KFlux10p=[]
KFlux10m=[]
KFluxErr=[]
KFluxN=[]
KFluxNH=[]
KFluxNM=[]
KFluxNErr=[]
KFluxEind=[]
KFluxEindH=[]
KFluxEindM=[]
KFlux10pEind=[]
KFlux10mEind=[]

p=0
if debug:
    print "Reading input file -+-+-", kaFile.name
if debug>1:
    print " line\tenergy\t  flx\t  flx/energy**SpecIndex"
while 1:
    line=kaFile.readline()
    if not line:
        break;
    energy=float(line.split(' ')[0])     #GeV
    flx=float(line.split(' ')[1])        #get flux as E^2.75*dN/dE
    Kbins.append(energy)
    KFlux.append(flx/energy**SpecIndex)  #flux as dN/dE
    KFluxEind.append(flx)                #flux as E^2.75*dN/dE
    p+=1
    if debug>1:
        print " %2d %10.4f %10.1f %12.6f"%(p, energy, flx, flx/energy**SpecIndex)
kaFile.close()
if debug:
    print "  Filling vector KFlux[",len(KFlux),"]"
if debug>1:
    print "len(KFlux) =", len(KFlux)
    print "KFlux = ", KFlux

r=0
if debug:
    print "Reading input file -+-+-", kaFileM.name
if debug>1:
    print " line\tenergy\t  flx\t  flx/energy**SpecIndex"
while 1:
    line=kaFileM.readline()
    if not line:
        break;
    energy=float(line.split(' ')[0])     #GeV
    flx=float(line.split(' ')[1])        #get flux as E^2.75*dN/dE
    KFluxM.append(flx/energy**SpecIndex)  #flux as dN/dE
    KFluxEindM.append(flx)                #flux as E^2.75*dN/dE
    r+=1
    if debug>1:
        print " %2d %10.4f %10.1f %12.6f"%(r, energy, flx, flx/energy**SpecIndex)
kaFileM.close()
if debug:
    print "  Filling vector KFluxM[",len(KFluxM),"]"
if debug>1:
    print "len(KFluxM) =", len(KFluxM)
    print "KFluxM = ", KFluxM

q=0
if debug:
    print "Reading input file -+-+-", kaFileH.name
if debug>1:
    print " line\tenergy\t  flx\t  flx/energy**SpecIndex"
while 1:
    line=kaFileH.readline()
    if not line:
        break;
    energy=float(line.split(' ')[0])     #GeV
    flx=float(line.split(' ')[1])        #get flux as E^2.75*dN/dE
    KFluxH.append(flx/energy**SpecIndex)  #flux as dN/dE
    KFluxEindH.append(flx)                #flux as E^2.75*dN/dE
    q+=1
    if debug>1:
        print " %2d %10.4f %10.1f %12.6f"%(q, energy, flx, flx/energy**SpecIndex)
kaFileH.close()
if debug:
    print "  Filling vector KFluxH[",len(KFluxH),"]"
if debug>1:
    print "len(KFluxH) =", len(KFluxH)
    print "KFluxH = ", KFluxH

### Get Kamae flux with dgamma+1%
q=0
if debug:
    print "Reading input file -+-+-", kaFile10p.name
if debug>1:
    print " line\tenergy\t  flx\t  flx/energy**SpecIndex"
while 1:
    line=kaFile10p.readline()
    if not line:
        break;
    energy10p=float(line.split(' ')[0])    #GeV
    flx10p=float(line.split(' ')[1])       #get flux as E^2.75*E*dN/dE
    KFlux10p.append(flx10p/energy10p**SpecIndex) #flux as dN/dE
    KFlux10pEind.append(flx10p)            #flux as E^2.75*dN/dE
    q+=1
    if debug>1:
        print " %2d %10.4f %10.1f %12.6f"%(q, energy10p, flx10p, flx10p/energy10p**SpecIndex)
kaFile10p.close()
if debug:
    print "  Filling vector KFlux10p[",len(KFlux10p),"]"
if debug>1:    
    print "len(KFlux10p) =", len(KFlux10p)
    print "KFlux10p = ", KFlux10p

### Get Kamae flux with dgamma-1%
r=0
if debug:
    print "Reading input file -+-+-", kaFile10m.name
if debug>1:
    print " line\tenergy\t  flx\t  flx/energy**SpecIndex"
while 1:
    line=kaFile10m.readline()
    if not line:
        break;
    energy10m=float(line.split(' ')[0])    #GeV
    flx10m=float(line.split(' ')[1])       #get flux as E^2.75*E*dN/dE
    KFlux10m.append(flx10m/energy10m**SpecIndex) #flux as dN/dE
    KFlux10mEind.append(flx10m)            #flux as E^2.75*dN/dE
    r+=1
    if debug>1:
        print " %2d %10.4f %10.1f %12.6f"%(r, energy10m, flx10m, flx10m/energy10m**SpecIndex)
kaFile10m.close()
if debug:
    print "  Filling vector KFlux10m[",len(KFlux10m),"]"
if debug>1:    
    print "len(KFlux10m) =", len(KFlux10m)
    print "KFlux10m = ", KFlux10m

if debug>1:
    print "Kbins created from file -+-+- len(Kbins) = ", len(Kbins)
    for mn in range(len(Kbins)): print "Kbins[%2d] = %10.4f"%(mn, Kbins[mn])

### Create Kamae energy bins arrays 
for j in range(len(Kbins)-1):
    dEK=Kbins[j+1]-Kbins[j]
    #KbinsErr.append(dEK/2.)   # removed for the fit's sake
    KbinsErr.append(0.)
KbinsErr.append(0.)

KbinsArray=array('d',Kbins)
KbinsErrArray=array('d',KbinsErr)

if debug:
    print "  len(Kbins)       =",len(Kbins)," -+-+-  len(KbinsErr) =",len(KbinsErr)
    if debug>1:
        print "KbinsErr =", KbinsErr
    print "  len(KFluxEind)   =",len(KFluxEind)
    print "-+-+- Kamae flux arrays:"
    print "ind |   Kbins    |   KFluxEind  |  KFluxEindH  |  KFluxEindM  "
    for j in range(len(Kbins)):
        print '%2d  | %9.4f  | %10.1f   | %10.1f   | %10.1f   '%(j, Kbins[j], KFluxEind[j], KFluxEindH[j], KFluxEindM[j])



############ DATA ############
if debug:
    print "===== DATA ====="
    print "Maps   from file -+-+-", emFile.GetName()
    print "Counts from file -+-+-", cmFile.GetName()

### Set energy bins and load Count and Exposure map
Ebins=[]
cntMap=[]
expMap=[]
mn=0
for i in range(len(Edec)-1):
    for j in range(NEbins[i]):
        Ebins.append(10**(Edec[i]+float(j)/NEbins[i]))
        cntname='cntmap%d'%(mn)
        expname='totexpmap%d'%(mn)
        cntMap.append(cmFile.Get(cntname))
        expMap.append(emFile.Get(expname))
        mn+=1
Ebins.append(10**Edec[-1])

if debug:
    print "Ebins created from settings -+-+- len(Ebins) = ", len(Ebins)
if debug>1:
    for mn in range(len(Ebins)): print "Ebins[%2d] = %10.4f"%(mn, Ebins[mn])

if debug:
    print "  len(cntMap) =", len(cntMap)
    print "  len(expMap) =", len(expMap)
if (len(cntMap) != len(Ebins)-1):
    print "ERROR! Maps and energy bins do not correspond"
    sys.exit()
if (len(cntMap) != len(expMap)):
    print "ERROR! Count maps and exposure maps do not correspond"
    sys.exit()

## Solid angle
pi=3.14159 
strMap=expMap[0].Clone()
strMap.SetName("strMap")
for j in range(strMap.GetNbinsY()):
    str=(phw*pi/180.)*(thw*pi/180.)*sin(pi*(j+0.5)*thw/180.)  #Solid angle d(Phi)d(Theta)sin(Theta) #Corrected
    for i in range(strMap.GetNbinsX()):
        strMap.SetBinContent(i+1,j+1,str)

flxMap=[]
invExpMap=[]
##  Differential flux map from raw counts and exposure
for i in range(len(cntMap)):
    ftmp=cntMap[i].Clone()
    flxMap.append(ftmp)
    flxMap[i].Divide(flxMap[i], expMap[i])
    flxMap[i].Scale(1e4)        #cm^2-->m^2

## 1/exposure map for bg subtraction
for i in range(len(expMap)):
    etpm=expMap[i].Clone()
    invExpMap.append(etpm)
    for j in range(invExpMap[i].GetNbinsY()):
        for l in range(invExpMap[i].GetNbinsX()):
            invExpMap[i].SetBinContent(l+1,j+1,1.)
    invExpMap[i].Divide(invExpMap[i], expMap[i])
    #invExpMap[i].Scale(1e4)        #cm^2-->m^2

## Exposure map in the right units
    expMap[i].Scale(1e-4)        #cm^2-->m^2    


## Raw counts and Background from diffuse emission 
c5=TCanvas("c5","c5",700,800)
c5.Divide(1,3)
c4=TCanvas("c4","c4",500,300)
c6=TCanvas("c6","c6",700,800)
c6.Divide(1,3)

thbr=[int(round(thrn[0]*thnb+1)),int(round(thrn[1]*thnb))]
nthb=thbr[1]-thbr[0]

if debug:
    print "  SIGNAL theta range = %5.1f - %5.1f deg -+-+- bin width = %3.1f deg, used [%3d, %3d]"%(thrn[0],thrn[1],thw,thbr[0],thbr[1])
    print "  BG theta range     = %5.1f - %5.1f deg -+-+- used [%3d, %3d]"%(thBGrn[0],thBGrn[1],thBGbr[0],thBGbr[1])

rawcnts=[]
sigevs=[]
bgevs=[]
nbgevs=[]

for i in range(len(Ebins)-1):
    if debug>1:
        print 'i=%d  All Counts=%d  Bg Events=%d'%(i,cntMap[i].GetEntries(),cntMap[i].Integral(1,cntMap[i].GetNbinsX(),thBGbr[0],thBGbr[1]))
    rawcnts.append(cntMap[i].GetEntries())
    sigevs.append(cntMap[i].Integral(1,cntMap[i].GetNbinsX(),thbr[0],thbr[1]))
    bgevs.append(cntMap[i].Integral(1,cntMap[i].GetNbinsX(),thBGbr[0],thBGbr[1]))
    nbgevs.append(round(bgevs[i]/nthBGb*nthb))
if debug:
    print "   Raw counts before theta selection:   ", rawcnts
    print "   Events in signal theta range:        ", sigevs
    print "   Events in bg theta range:            ", bgevs, " in ",nthBGb," th bins and ", 360./phw, " phi bins."
    print "   ==> Bg events in signal theta range: ", nbgevs, " in ",nthb," th bins and ", 360./phw, " phi bins."

LFluxEind=[]
LFluxEindErr=[]
KLRatio=[]
KLRatioH=[]
KLRatioM=[]
KLRatioErr=[]
## Array of energy bin CENTERS!
CEbins = []
CEbinsErr=[]

bb=0
if debug:
    print "  -+-+- Limb flux arrays:"
    print "  ind |   Ebins   | CntPerBin |     CEbins  [  Kbins  ]  |  KLRatio   |  KLRatioH  |  KLRatioM  |  LFluxEind "

dOmega=strMap.Integral(1,strMap.GetNbinsX(),thbr[0],thbr[1])

for i in range(len(Ebins)-1):

## Energy bin width
    dE=Ebins[i+1]-Ebins[i]

## Raw counts    
    CntPerBin=cntMap[i].Integral(1,cntMap[i].GetNbinsX(),thbr[0],thbr[1])
    
## Exposure
    AvgExp=expMap[i].Integral(1,expMap[i].GetNbinsX(),thbr[0],thbr[1])

## Limb DATA flux as dN/dE
    if not subtract_bg:
        thisLFlux=flxMap[i].Integral(1,flxMap[i].GetNbinsX(),thbr[0],thbr[1])/dOmega/dE
    else:
        thisLFlux=(flxMap[i].Integral(1,flxMap[i].GetNbinsX(),thbr[0],thbr[1]) - nbgevs[i]*invExpMap[i].Integral(1,invExpMap[i].GetNbinsX(),thbr[0],thbr[1]))/dOmega/dE
        thisLFluxErr=thisLFlux/(CntPerBin+0.5)**0.5
    #if (CntPerBin == 1): thisLFluxErr=0.999**0.5/dE/AvgExp                 # avoid error going to y=0
    if debug:
        print '  %2d  | %9.4f | %8d  |'%(i,Ebins[i],CntPerBin),
    if (CntPerBin == 0): print " "
    if (CntPerBin != 0):
        thisEC = Ebins[i]+dE/2.
        CEbins.append(thisEC)
        CEbinsErr.append(0.)

## Limb DATA flux as E^2.75*dN/dE
        thisLFluxEind=thisEC**SpecIndex*thisLFlux
        thisLFluxEindErr=thisEC**SpecIndex*thisLFluxErr
        LFluxEind.append(thisLFluxEind)
        LFluxEindErr.append(thisLFluxEindErr)
    
        for j in range(len(Kbins)-1):
            if (abs(CEbins[i]-Kbins[j])/CEbins[i] < 0.05):
                bb+=1
                
# Ratio DATA/Kamae fluxes as (dN/dE)/(dN(kamae)/dE) to get NormConst
                KLRatio.append(thisLFlux/KFlux[j])
                KLRatioH.append(thisLFlux/KFluxH[j])
                KLRatioM.append(thisLFlux/KFluxM[j])
                KLRatioErr.append(thisLFluxErr/KFlux[j])
                if debug:
                    print '  %9.4f [%9.4f]  | %5.3e  | %5.3e  | %5.3e  |  %7.5f  '%(CEbins[i],Kbins[j],KLRatio[i],KLRatioH[i],KLRatioM[i],LFluxEind[i])

if debug:
    print "  len(CEbins) = ", len(CEbins)," -+-+- len(CEbinsErr) = ", len(CEbinsErr) 
if debug>1:
    print "CEbins = ",CEbins
    print "CEbinsErr = ",CEbinsErr

CEbinsArray=array('d',CEbins)
CEbinsErrArray=array('d',CEbinsErr)
LFluxEindArray = array('d',LFluxEind)
LFluxEindErrArray = array('d',LFluxEindErr)
KLRatioArray = array('d',KLRatio)
KLRatioArrayH = array('d',KLRatioH)
KLRatioArrayM = array('d',KLRatioM)
KLRatioErrArray = array('d',KLRatioErr)
flxSpecEindTG=TGraphErrors(len(CEbins),CEbinsArray,LFluxEindArray,CEbinsErrArray,LFluxEindErrArray)
pltit='Flux*E^%4.2f #theta_{Nadir} = %4.1f-%4.1f [bins %d-%d]'%(SpecIndex,thrn[0],thrn[1],thbr[0],thbr[1])
flxSpecEindTG.SetTitle(pltit)
plname='flxSpecEindTG'
flxSpecEindTG.SetName(plname)

flxRatioTG=TGraphErrors(len(CEbins),CEbinsArray,KLRatioArray,CEbinsErrArray,KLRatioErrArray)
pltit='Ratio Flux/Kamae #theta_{Nadir} = %4.1f-%4.1f [protons]'%(thrn[0],thrn[1])
flxRatioTG.SetTitle(pltit)
plname='flxRatioTG'
flxRatioTG.SetName(plname)
flxRatioTGH=TGraphErrors(len(CEbins),CEbinsArray,KLRatioArrayH,CEbinsErrArray,KLRatioErrArray)
pltit='Ratio Flux/Kamae #theta_{Nadir} = %4.1f-%4.1f [helium]'%(thrn[0],thrn[1])
flxRatioTGH.SetTitle(pltit)
plname='flxRatioTGH'
flxRatioTGH.SetName(plname)
flxRatioTGM=TGraphErrors(len(CEbins),CEbinsArray,KLRatioArrayM,CEbinsErrArray,KLRatioErrArray)
pltit='Ratio Flux/Kamae #theta_{Nadir} = %4.1f-%4.1f [p+0.085He]'%(thrn[0],thrn[1])
flxRatioTGM.SetTitle(pltit)
plname='flxRatioTGM'
flxRatioTGM.SetName(plname)

### Init fit function

fitRatio = TF1( "fitRatio", "[0]+pow(10,-30)*x", Ebins[0], Ebins[4])  #this one for coincidence in low-en bins
fitRatioH = TF1( "fitRatioH", "[0]+pow(10,-30)*x", Ebins[0], Ebins[4])  #this one for coincidence in low-en bins
fitRatioM = TF1( "fitRatioM", "[0]+pow(10,-30)*x", Ebins[0], Ebins[4])  #this one for coincidence in low-en bins

############ DRAW ############

# Ratio of flux/kamae to find normalization constant
c6.cd(1)
c6.cd(1).SetLogx()
flxRatioTG.SetMarkerStyle(otherMarker)
flxRatioTG.SetMarkerSize(.6)
flxRatioTG.GetXaxis().SetTitle("E [GeV]")
flxRatioTG.Draw("APZ")
fitRatio.SetLineColor(kRed-3);
fitRatio.SetLineWidth(1);
if debug:
    print "\n ",c6.GetName()," -+-+- Fit ",flxRatioTG.GetName(), "with f(x) =",fitRatio.GetExpFormula()
    print ' ',
flxRatioTG.Fit("fitRatio","R")
if debug:
    print "  CHISQ/NDOF= %10.4f / %4d"%(fitRatio.GetChisquare(), fitRatio.GetNDF())
pars = fitRatio.GetParameters()
dpars = fitRatio.GetParErrors()
chisq = fitRatio.GetChisquare()
prob = fitRatio.GetProb()
ndof = fitRatio.GetNDF()
txt=("#splitline{Normalization constant = %8.3e}{#chi^{2}/ndf = %8.2f / %d}")%(pars[0],chisq,ndof)
testo = []
testo.append(TLatex(200,0.000025,txt))
testo[-1].SetTextColor(kRed+2)
testo[-1].SetTextSize(0.04)
testo[-1].Draw()
gPad.Update()
NormConst = pars[0]
NormConst += 0.1*NormConst

c6.cd(2)
c6.cd(2).SetLogx()
flxRatioTGH.SetMarkerStyle(otherMarker)
flxRatioTGH.SetMarkerSize(.6)
flxRatioTGH.GetXaxis().SetTitle("E [GeV]")
flxRatioTGH.Draw("APZ")
fitRatioH.SetLineColor(kRed-3);
fitRatioH.SetLineWidth(1);
if debug:
    print "\n ",c6.GetName()," -+-+- Fit ",flxRatioTGH.GetName(), "with f(x) =",fitRatioH.GetExpFormula()
    print ' ',
flxRatioTGH.Fit("fitRatioH","R")
if debug:
    print "  CHISQ/NDOF= %10.4f / %4d"%(fitRatioH.GetChisquare(), fitRatioH.GetNDF())
parsH = fitRatioH.GetParameters()
dparsH = fitRatioH.GetParErrors()
chisqH = fitRatioH.GetChisquare()
probH = fitRatioH.GetProb()
ndofH = fitRatioH.GetNDF()
txt=("#splitline{Normalization constant = %8.3e}{#chi^{2}/ndf = %8.2f / %d}")%(parsH[0],chisqH,ndofH)
testoH = []
testoH.append(TLatex(200,0.000025,txt))
testoH[-1].SetTextColor(kRed+2)
testoH[-1].SetTextSize(0.04)
testoH[-1].Draw()
gPad.Update()
NormConstH = parsH[0]
NormConstH += 0.1*NormConstH

c6.cd(3)
c6.cd(3).SetLogx()
flxRatioTGM.SetMarkerStyle(otherMarker)
flxRatioTGM.SetMarkerSize(.6)
flxRatioTGM.GetXaxis().SetTitle("E [GeV]")
flxRatioTGM.Draw("APZ")
fitRatioM.SetLineColor(kRed-3);
fitRatioM.SetLineWidth(1);
if debug:
    print "\n ",c6.GetName()," -+-+- Fit ",flxRatioTGM.GetName(), "with f(x) =",fitRatioM.GetExpFormula()
    print ' ',
flxRatioTGM.Fit("fitRatioM","R")
if debug:
    print "  CHISQ/NDOF= %10.4f / %4d"%(fitRatioM.GetChisquare(), fitRatioM.GetNDF())
parsM = fitRatioM.GetParameters()
dparsM = fitRatioM.GetParErrors()
chisqM = fitRatioM.GetChisquare()
probM = fitRatioM.GetProb()
ndofM = fitRatioM.GetNDF()
txt=("#splitline{Normalization constant = %8.3e}{#chi^{2}/ndf = %8.2f / %d}")%(parsM[0],chisqM,ndofM)
testoM = []
testoM.append(TLatex(200,0.000025,txt))
testoM[-1].SetTextColor(kRed+2)
testoM[-1].SetTextSize(0.04)
testoM[-1].Draw()
gPad.Update()
NormConstM = parsM[0]
NormConstM += 0.1*NormConstM
gfnam="FluxRatio_p-He_THETA%4.1f-%4.1f_%s_%s.gif"%(thrn[0],thrn[1],escFact,irfV)
#if save_plots:
#    c6.SaveAs(gfnam)

chisqList=[chisq/float(ndof),chisqH/float(ndofH),chisqM/float(ndofM)]
scanList=[1.,2.,3.]
scanLabels=['p','He','p+8.5He']
chisqArray = array('d',chisqList)
scanArray = array('d',scanList)
chisqScanTG=TGraph(3,scanArray,chisqArray)
c4.cd()
chisqScanTG.SetMarkerStyle(20)
chisqScanTG.SetMarkerSize(.6)
chisqScanTG.GetYaxis().SetTitle("#chi^{2}/NDoF")
xax = chisqScanTG.GetXaxis()
xax.LabelsOption("h")        # NON FUNGE...
i = 1
while i <= xax.GetXmax():
    bin_index = xax.FindBin(i)
#    xax.SetBinLabel(bin_index,scanLabels[bin_index-1])
    if (bin_index == 9): 
        xax.SetBinLabel(bin_index,scanLabels[0])
    elif (bin_index == 50): 
        xax.SetBinLabel(bin_index,scanLabels[1])
    elif (bin_index == 92): 
        xax.SetBinLabel(bin_index,scanLabels[2])
    print "i=",i," bin_index=",bin_index
    i+=1
chisqScanTG.Draw("APZ")
gPad.Update()


### Create Kamae normalized flux arrays
KFluxEindN=[]
KFluxEindNH=[]
KFluxEindNM=[]
KFluxEindNErr=[]
for j in range(len(Kbins)-1):
### Error on Kamae flux is |KFlux10p-KFlux10m|/2.
    KFluxEindN.append(NormConst*KFluxEind[j])
    KFluxEindNH.append(NormConstH*KFluxEindH[j])
    KFluxEindNM.append(NormConstM*KFluxEindM[j])
    KFluxEindNErr.append(NormConst*abs(KFlux10pEind[j]-KFlux10mEind[j])/2.)    
KFluxEindN.append(NormConst*KFluxEind[-1])
KFluxEindNH.append(NormConstH*KFluxEindH[-1])
KFluxEindNM.append(NormConstM*KFluxEindM[-1])
KFluxEindNErr.append(NormConst*abs(KFlux10pEind[-1]-KFlux10mEind[-1])/2.)

if debug:
    print "  len(KFluxEindN)  =",len(KFluxEindN)," -+-+-  len(KFluxEindNErr) =",len(KFluxEindNErr)
    print "  len(KFluxEindNH)  =",len(KFluxEindNH)
    print "  len(KFluxEindNM)  =",len(KFluxEindNM)
    
    print "  -+-+- Kamae normalized flux arrays:"
    print "  ind |   Kbins    |  KFluxEindN  |  KFluxEindNH |  KFluxEindNM [KFluxEindErr]"
    for j in range(len(Kbins)):
        print '  %2d  | %9.4f  |  %10.7f  |  %10.7f  |  %10.7f  [ %10.7f ]'%(j, Kbins[j], KFluxEindN[j], KFluxEindNH[j], KFluxEindNM[j], KFluxEindNErr[j])


### Fill normalized graphs for Kamae flux
KFluxEindNArray = array('d',KFluxEindN)
KFluxEindNArrayH = array('d',KFluxEindNH)
KFluxEindNArrayM = array('d',KFluxEindNM)
KFluxEindNErrArray = array('d',KFluxEindNErr)
flxKamaeEindNTG = TGraphErrors(len(Kbins),KbinsArray,KFluxEindNArray,KbinsErrArray,KFluxEindNErrArray)
pltit = 'Kamae(2006) flux*E^%4.2f from protons'%(SpecIndex)
flxKamaeEindNTG.SetTitle(pltit)
plname='flxKamaeEindNTG'
flxKamaeEindNTG.SetName(plname)
flxKamaeEindNTGH = TGraphErrors(len(Kbins),KbinsArray,KFluxEindNArrayH,KbinsErrArray,KFluxEindNErrArray)
pltit = 'Kamae(2006) flux*E^%4.2f from He'%(SpecIndex)
flxKamaeEindNTGH.SetTitle(pltit)
plname='flxKamaeEindNTGH'
flxKamaeEindNTGH.SetName(plname)
flxKamaeEindNTGM = TGraphErrors(len(Kbins),KbinsArray,KFluxEindNArrayM,KbinsErrArray,KFluxEindNErrArray)
pltit = 'Kamae(2006) flux*E^%4.2f from p+0.085He'%(SpecIndex)
flxKamaeEindNTGM.SetTitle(pltit)
plname='flxKamaeEindNTGM'
flxKamaeEindNTGM.SetName(plname)


# Limb flux over Kamae flux
c5.cd(1)
c5.cd(1).SetLogx()
c5.cd(1).SetGridx()
c5.cd(1).SetGridy()
flxSpecEindTG.GetXaxis().SetTitle("E [GeV]")
flxSpecEindTG.GetYaxis().SetTitle("Flux #times E^{2.75} (m^{-2} s^{-1} sr^{-1} GeV^{1.75})")
flxSpecEindTG.SetMarkerStyle(limbMarker)
flxSpecEindTG.SetMarkerColor(myDarkBlue)
flxSpecEindTG.SetLineColor(myDarkBlue)
flxSpecEindTG.SetFillColor(0)
flxSpecEindTG.SetMarkerSize(0.6)
flxSpecEindTG.Draw("APZ");
flxKamaeEindNTG.SetFillColor(myLightRosso)
flxKamaeEindNTG.SetFillStyle(3002)
flxKamaeEindNTG.SetLineColor(0)
flxKamaeEindNTG.SetLineStyle(3)
flxKamaeEindNTG.SetMarkerColor(myLightRosso)
flxKamaeEindNTG.Draw("3,SAME")
leg = TLegend(0.15,0.70,0.6,0.85)
if escale:
    leg.AddEntry(flxSpecEindTG,"Gamma-ray flux from Earth Limb ESCALE [ULTRACLEANVETO]")
    leg.AddEntry(flxKamaeEindNTG,"Gamma-ray flux from Kamae(2006) protons")
else:
    leg.AddEntry(flxSpecEindTG,"Gamma-ray flux from Earth Limb [ULTRACLEANVETO]")
    leg.AddEntry(flxKamaeEindNTG,"Gamma-ray flux from Kamae(2006) protons")
leg.SetTextSize(.03)
leg.SetFillColor(0)
leg.SetLineColor(0)
leg.Draw()
gPad.Update()
        
c5.cd(2)
c5.cd(2).SetLogx()
c5.cd(2).SetGridx()
c5.cd(2).SetGridy()
flxSpecEindTG.Draw("APZ");
flxKamaeEindNTGH.SetFillColor(myLightBlue)
flxKamaeEindNTGH.SetFillStyle(3002)
flxKamaeEindNTGH.SetLineColor(0)
flxKamaeEindNTGH.SetLineStyle(3)
flxKamaeEindNTGH.SetMarkerColor(myLightBlue)
flxKamaeEindNTGH.Draw("3,SAME")
legH = TLegend(0.15,0.70,0.6,0.85)
if escale:
    legH.AddEntry(flxSpecEindTG,"Gamma-ray flux from Earth Limb ESCALE [ULTRACLEANVETO]")
    legH.AddEntry(flxKamaeEindNTGH,"Gamma-ray flux from Kamae(2006) Helium")
else:
    legH.AddEntry(flxSpecEindTG,"Gamma-ray flux from Earth Limb [ULTRACLEANVETO]")
    legH.AddEntry(flxKamaeEindNTGH,"Gamma-ray flux from Kamae(2006) Helium")
legH.SetTextSize(.03)
legH.SetFillColor(0)
legH.SetLineColor(0)
legH.Draw()
gPad.Update()
        
c5.cd(3)
c5.cd(3).SetLogx()
c5.cd(3).SetGridx()
c5.cd(3).SetGridy()
flxSpecEindTG.Draw("APZ");
flxKamaeEindNTGM.SetFillColor(kMagenta-5)
flxKamaeEindNTGM.SetFillStyle(3002)
flxKamaeEindNTGM.SetLineColor(0)
flxKamaeEindNTGM.SetLineStyle(3)
flxKamaeEindNTGM.SetMarkerColor(kMagenta-5)
flxKamaeEindNTGM.Draw("3,SAME")
legM = TLegend(0.15,0.70,0.6,0.85)
if escale:
    legM.AddEntry(flxSpecEindTG,"Gamma-ray flux from Earth Limb ESCALE [ULTRACLEANVETO]")
    legM.AddEntry(flxKamaeEindNTG,"Gamma-ray flux from Kamae(2006) p+0.085He")
else:
    legM.AddEntry(flxSpecEindTG,"Gamma-ray flux from Earth Limb [ULTRACLEANVETO]")
    legM.AddEntry(flxKamaeEindNTG,"Gamma-ray flux from Kamae(2006) p+0.085He")
legM.SetTextSize(.03)
legM.SetFillColor(0)
legM.SetLineColor(0)
legM.Draw()
gPad.Update()
gfnam="LimbFlux_p-He_OnKamaeBand_THETA%4.1f-%4.1f_%s_%s.gif"%(thrn[0],thrn[1],escFact,irfV)
if save_plots:
    print "\n",
    c5.SaveAs(gfnam)
    print "\n",
        

    
print "  ==== BIAS SUMMARY ===="
print "  Background subtraction (total %d/%d events) makes a %5.3f"%(sum(bgevs), sum(sigevs), sum(bgevs)/float(sum(sigevs))*100.),"% effect"
print "  Normalization Constant for p:         [from Ratio fit =", pars[0], "], used =", NormConst
print "  Normalization Constant for He:        [from Ratio fit =", parsH[0], "], used =", NormConstH
print "  Normalization Constant for p+0.085He: [from Ratio fit =", parsM[0], "], used =", NormConstM

if escale:
    print "  Using constant energy pre-scaling 0.963"
else:
    print "Using NO constant energy pre-scaling"
        
print "\nDone!"

raw_input()
