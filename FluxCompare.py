#!/usr/bin/env python

############ INIT ############

import numpy
import sys
from ROOT import *
from math import *
from array import *
from copy import *

### User settings
debug = 2

escale      = true    # use 0.963 constant E scaling
escale500   = false   # (!) usa in cascata a escale (!)
escale750   = false   # (!) usa in cascata a escale (!)

subtract_bg = true    # subtract BG from diffuse emission
thin_target = true    # thin target approximation

use_kamae   = false   # use Kamae model - otherwise: Ostapchenko
use_onlyP   = true    # use gammas from proton interaction only
use_onlyHe  = false   # use gammas from Helium interaction only

use_glfit   = false   # use parameters from the total data fit
                      # (otherwise pars from the low- and high-en parts separately)
save_plots  = false
nDraw=[0,3,6,9,10]    # Plot maps for these i^th bin of Ebins [DATA]


### Fixed analysis settings
Edec=[1,2,3,4]         # Energy decade (GeV) [DATA]
NEbins=[5,5,5]         # Number of bins between 10**Edec [DATA]
thbw=[400,0.,80.]      # Theta [N bins, theta1, theta2] [DATA]
if thin_target:
    thrn=[68.4,70.]        # Spectrum for this theta range [DATA] thin target
else:
    thrn=[66.,70.]         # Spectrum for this theta range [DATA]
thBGrn=[77.,80.]       # Diff. BG from this theta range [DATA]

### Constants
SpecIndex = 2.75       # spectral index
NormConst=1.688e-05     

### Root settings
gROOT.SetStyle('Plain')
gStyle.SetPalette(1)
gStyle.SetErrorX(0)
gStyle.SetOptFit(0)
gStyle.SetOptStat(0)

modelMarker = 24
limbMarker  = 25
otherMarker = 26

myMarkerColor  =TColor.GetColor("#000000")
myLightGreen   =TColor.GetColor("#55CCCC")
myDarkGreen    =TColor.GetColor("#338888")
myLightBlue    =TColor.GetColor("#55AADD")
myDarkBlue     =TColor.GetColor("#3366AA")
myLightMagenta =TColor.GetColor("#FF99DD")
myDarkMagenta  =TColor.GetColor("#AA4488")
myLightRosso   =TColor.GetColor("#FF9999")
myDarkRosso    =TColor.GetColor("#AA5555")
kamaeColor     =TColor.GetColor("#6677FF")
ostapColor     =TColor.GetColor("#77DDBB")

### Data sets
datV='P301v1'
irfV='P8V6ULTRACLEANVETO' #IRF version
if escale:
    cmFile=TFile('CntMap'+irfV+'_'+datV+'test963_thetashiftpeak67.798-2.root')  # WITH E*0.963-SCALE
    escFact='0.963ESCALE'
    if escale500:
        cmFile=TFile('CntMap'+irfV+'_'+datV+'test500.root')  # WITH E*0.500-SCALE
        escFact='0.500ESCALE'
    elif escale750:
        cmFile=TFile('CntMap'+irfV+'_'+datV+'test750.root')  # WITH E*0.750-SCALE
        escFact='0.750ESCALE'
else:
    cmFile=TFile('CntMap'+irfV+'_'+datV+'.root')             # NO E-SCALE
    escFact='NOESCALE'
emFile=TFile('ExpMap'+irfV+'_'+datV+'.root')

if use_kamae:
### Kamae Model
    if use_onlyP:
        kaFile=open("gammaspectrum.txt")
    elif use_onlyHe:
        kaFile=open("gammaspectrumHe.txt")
    else:
        kaFile=open("gammaspectrumM.txt")
    kaFile10p=open("gammaspectrumEp.txt")
    kaFile10m=open("gammaspectrumEm.txt")
else:
### Ostapchenko model
    if use_onlyP:
        kaFile=open("p+p-gamma.txt")
    elif use_onlyHe:
        kaFile=open("gammaspectrumHe.txt")
    else:
        kaFile=open("gammaspectrumM.txt")
    kaFile10p=open("p+p-gamma-Ep.txt")
    kaFile10m=open("p+p-gamma-Em.txt")

### Theta bin width, bins in 0.2 deg, used bins
thw=(thbw[2]-thbw[1])/thbw[0]  
thnb=int(round(thbw[0]/(thbw[2]-thbw[1])))  ## how many bins for 1 degree
thbr=[int(round(thrn[0]*thnb+1)),int(round(thrn[1]*thnb))]
nthb=thbr[1]-thbr[0]
thBGbr=[int(round(thBGrn[0]*thnb+1)),int(round(thBGrn[1]*thnb))]
nthBGb=thBGbr[1]-thBGbr[0]

### Phi bin width 5.0 deg
phw = 5.0

print "\nExecuting %s with settings:\n"%(__file__)
if debug:
    print "SIGNAL theta range = %5.1f - %5.1f deg -+-+- bin width = %3.1f deg, used [%3d, %3d]"%(thrn[0],thrn[1],thw,thbr[0],thbr[1])
    print "BG theta range     = %5.1f - %5.1f deg -+-+- used [%3d, %3d]"%(thBGrn[0],thBGrn[1],thBGbr[0],thBGbr[1])
    print "SpecIndex = %5.3f\n"%(SpecIndex)


############ KAMAE/OSTAP ############
### Lists for Model flux and energy bins
Kbins=[]
KbinsErr=[]
KFlux=[]
KFlux10p=[]
KFlux10m=[]
KFluxEind=[]
KFlux10pEind=[]
KFlux10mEind=[]

### Get input from files
p=0
if debug:
    print "Reading input file -+-+-", kaFile.name
if debug>2:
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
    if debug>2:
        print " %2d %10.4f %10.1f %12.6f"%(p, energy, flx, flx/energy**SpecIndex)
kaFile.close()
if debug>1:
    print "  Filling vector KFlux[",len(KFlux),"]"
if debug>2:
    print "KFlux = ", KFlux

### Get Model flux with dgamma+1%
q=0
if debug:
    print "Reading input file -+-+-", kaFile10p.name
if debug>2:
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
    if debug>2:
        print " %2d %10.4f %10.1f %12.6f"%(q, energy10p, flx10p, flx10p/energy10p**SpecIndex)
kaFile10p.close()
if debug>1:
    print "  Filling vector KFlux10p[",len(KFlux10p),"]"
if debug>2:    
    print "KFlux10p = ", KFlux10p

### Get Model flux with dgamma-1%
r=0
if debug:
    print "Reading input file -+-+-", kaFile10m.name
if debug>2:
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
    if debug>2:
        print " %2d %10.4f %10.1f %12.6f"%(r, energy10m, flx10m, flx10m/energy10m**SpecIndex)
kaFile10m.close()
if debug>1:
    print "  Filling vector KFlux10m[",len(KFlux10m),"]"
if debug>2:    
    print "KFlux10m = ", KFlux10m

if debug>1:
    print "Kbins created from file     -+-+- len(Kbins)   	     =",len(Kbins)
if debug>2:
    for mn in range(len(Kbins)): print "Kbins[%2d] = %10.4f"%(mn, Kbins[mn])

### Create Model energy bins arrays 
for j in range(len(Kbins)-1):
    dEK=Kbins[j+1]-Kbins[j]
    #KbinsErr.append(dEK/2.)   # removed for the fit's sake
    KbinsErr.append(0.)
KbinsErr.append(0.)

if debug>1:
    print "                            -+-+- len(KbinsErr)      =",len(KbinsErr)
if debug>2:
    print "KbinsErr =", KbinsErr


############ DATA ############
if debug:
    print "Maps   from file   -+-+-", emFile.GetName()
    print "Counts from file   -+-+-", cmFile.GetName()

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

if debug>1:
    print "Ebins created from settings -+-+- len(Ebins)   	     =",len(Ebins)
if debug>2:
    for mn in range(len(Ebins)): print "Ebins[%2d] = %10.4f"%(mn, Ebins[mn])

if debug>1:
    print "                            -+-+- len(cntMap)  	     =",len(cntMap)
    print "                            -+-+- len(expMap)  	     =",len(expMap)
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

##  Differential flux map from raw counts and exposure
flxMap=[]
for i in range(len(cntMap)):
    ftmp=cntMap[i].Clone()
    flxMap.append(ftmp)
    flxMap[i].Divide(flxMap[i], expMap[i])
    flxMap[i].Scale(1e4)        #cm^2-->m^2

invExpMap=[]
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

rawcnts=[]
sigevs=[]
bgevs=[]
nbgevs=[]

for i in range(len(Ebins)-1):
    if debug>2:
        print 'i=%3d  All Counts=%6d  Bg Events=%4d'%(i,cntMap[i].GetEntries(),cntMap[i].Integral(1,cntMap[i].GetNbinsX(),thBGbr[0],thBGbr[1]))
    rawcnts.append(cntMap[i].GetEntries())
    sigevs.append(cntMap[i].Integral(1,cntMap[i].GetNbinsX(),thbr[0],thbr[1]))
    bgevs.append(cntMap[i].Integral(1,cntMap[i].GetNbinsX(),thBGbr[0],thBGbr[1]))
    nbgevs.append(round(bgevs[i]/nthBGb*nthb))
if debug:
    print "Raw counts before theta selection:   ", rawcnts
    print "Events in signal theta range:        ", sigevs
    print "Events in bg theta range:            ", bgevs, " in ",nthBGb," th bins and ", 360./phw, " phi bins."
    print "==> Bg events in signal theta range: ", nbgevs, " in ",nthb," th bins and ", 360./phw, " phi bins."


############ RATIO ############
LFlux=[]
LFluxErr=[]
LFluxEind=[]
LFluxEindErr=[]
KLRatio=[]
KLRatioErr=[]

## Array of energy bin CENTERS!
CEbins = []
CEbinsErr=[]

MatchInd=[]

if debug:
    print "\nGet NormConst -+-+- Ratio arrays:"
    print "  ind |   Ebins   | CntPerBin |     CEbins  [  Kbins  ]  |  KLRatio   "

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

    if debug:
        print '  %2d  | %9.4f | %8d  |'%(i,Ebins[i],CntPerBin),
    if (CntPerBin == 0): print " "
    if (CntPerBin != 0):
        thisEC = Ebins[i]+dE/2.
        CEbins.append(thisEC)
        CEbinsErr.append(0.)

        LFlux.append(thisLFlux)
        LFluxErr.append(thisLFluxErr)
        
## Limb DATA flux as E^2.75*dN/dE
        thisLFluxEind=thisEC**SpecIndex*thisLFlux
        thisLFluxEindErr=thisEC**SpecIndex*thisLFluxErr
        LFluxEind.append(thisLFluxEind)
        LFluxEindErr.append(thisLFluxEindErr)
    
        for j in range(len(Kbins)-1):
            if (abs(CEbins[i]-Kbins[j])/CEbins[i] < 0.05):
                MatchInd.append(j)
# Ratio DATA/Model fluxes as (dN/dE)/(dN(model)/dE) to get NormConst
                KLRatio.append(thisLFlux/KFlux[j])
                KLRatioErr.append(thisLFluxErr/KFlux[j])
                if debug:
                    print '  %9.4f [%9.4f]  | %5.3e'%(CEbins[i],Kbins[j],KLRatio[i]),
                    print " j =",MatchInd[i]

if debug>1:
    print "\nCEbins created from Ebins   -+-+- len(CEbins)        =",len(CEbins)
    print "                            -+-+- len(CEbinsErr)     =",len(CEbinsErr) 
if debug>2:
    print "CEbins    = ",CEbins
    print "CEbinsErr = ",CEbinsErr

if debug>1:
    print "KLRatio created from input  -+-+- len(KLRatio)       =",len(KLRatio)

CEbinsArray=array('d',CEbins)
CEbinsErrArray=array('d',CEbinsErr)
KLRatioArray = array('d',KLRatio)
KLRatioErrArray = array('d',KLRatioErr)
flxRatioTG=TGraphErrors(len(CEbins),CEbinsArray,KLRatioArray,CEbinsErrArray,KLRatioErrArray)
if use_onlyP:
    pltit='Ratio Flux/Model #theta_{Nadir} = %4.1f-%4.1f [protons]'%(thrn[0],thrn[1])
elif use_onlyHe:
    pltit='Ratio Flux/Model #theta_{Nadir} = %4.1f-%4.1f [Helium]'%(thrn[0],thrn[1])
else:
    pltit='Ratio Flux/Model #theta_{Nadir} = %4.1f-%4.1f [p+0.085He]'%(thrn[0],thrn[1])
flxRatioTG.SetTitle(pltit)
plname='flxRatioTG'
flxRatioTG.SetName(plname)

### Init fit function
fitRatio = TF1( "fitRatio", "[0]+pow(10,-30)*x", Ebins[0], Ebins[4])  #this one for coincidence in low-en bins

### DRAW Ratio of flux/model to find normalization constant
c6=TCanvas("c6","c6 - Ratio")
c6.cd()
c6.SetLogx()
flxRatioTG.SetMarkerStyle(otherMarker)
flxRatioTG.SetMarkerSize(.7)
flxRatioTG.GetXaxis().SetTitle("E [GeV]")
flxRatioTG.Draw("APZ")
fitRatio.SetLineColor(kRed-3);
fitRatio.SetLineWidth(1);
if debug:
    print "\n",c6.GetName()," -+-+- Fit ",flxRatioTG.GetName(), "with f(x) =",fitRatio.GetExpFormula()
    print ' ',
flxRatioTG.Fit("fitRatio","R")
if debug:
    print "  CHISQ/NDOF= %10.4f / %4d"%(fitRatio.GetChisquare(), fitRatio.GetNDF())
pars3 = fitRatio.GetParameters()
dpars3 = fitRatio.GetParErrors()
chisq3 = fitRatio.GetChisquare()
prob3 = fitRatio.GetProb()
ndof3 = fitRatio.GetNDF()
txt=("#splitline{Normalization constant = %8.3e}{#chi^{2}/ndf = %8.2f / %d}")%(pars3[0],chisq3,ndof3)
NormConst = pars3[0]
NormConst += 0.08*NormConst
testo3 = []
testo3.append(TLatex(20,NormConst/2.,txt))
testo3[-1].SetTextColor(kRed+2)
testo3[-1].SetTextSize(0.04)
testo3[-1].Draw()
gPad.Update()
gfnam="FluxRatio_THETA%4.1f-%4.1f_%s_%s.gif"%(thrn[0],thrn[1],escFact,irfV)
if save_plots:
    c6.SaveAs(gfnam)

### Create normalized Model fluxes arrays
KFluxN=[]
KFluxNErr=[]
KFluxEindN=[]
KFluxEindNErr=[]
for j in range(len(Kbins)-1):
### Error on Model flux is |KFlux10p-KFlux10m|/2.
    KFluxN.append(NormConst*KFlux[j])
    KFluxNErr.append(NormConst*abs(KFlux10p[j]-KFlux10m[j])/2.)
    KFluxEindN.append(NormConst*KFluxEind[j])
    KFluxEindNErr.append(NormConst*abs(KFlux10pEind[j]-KFlux10mEind[j])/2.)
KFluxN.append(NormConst*KFlux[-1])
KFluxNErr.append(NormConst*abs(KFlux10p[-1]-KFlux10m[-1])/2.)
KFluxEindN.append(NormConst*KFluxEind[-1])
KFluxEindNErr.append(NormConst*abs(KFlux10pEind[-1]-KFlux10mEind[-1])/2.)

if debug>1:
    print "KFlux created from input    -+-+- len(KFlux)         =",len(KFlux)
    print "                            -+-+- len(KFluxN)        =",len(KFluxN)
    print "                            -+-+- len(KFluxNErr)     =",len(KFluxNErr)
    print "                            -+-+- len(KFluxEind)     =",len(KFluxEind)
    print "                            -+-+- len(KFluxEindN)    =",len(KFluxEindN)
    print "                            -+-+- len(KFluxEindNErr) =",len(KFluxEindNErr)

if debug:
    print "-+-+- Model flux arrays:"
    print "ind |   Kbins   |     KFlux    |    KFluxN    [ KFluxNErr]  |  KFluxEind   |  KFluxEindN  [KFluxEindNErr]"
    for j in range(len(Kbins)):
        print '%2d  | %9.4f |  %5.4e  |  %5.4e  [%5.4e]  | %10.1f   |  %10.7f  [ %10.7f  ]'%(j, Kbins[j], KFlux[j], KFluxN[j], KFluxNErr[j], KFluxEind[j], KFluxEindN[j], KFluxEindNErr[j])

if debug>1:
    print "\nLFlux created from input    -+-+- len(LFlux)         =",len(LFlux)
    print "                            -+-+- len(LFluxErr)      =",len(LFluxErr)    
    print "                            -+-+- len(LFluxEind)     =",len(LFluxEind)    
    print "                            -+-+- len(LFluxEindErr)  =",len(LFluxEindErr)

if debug:
    print "\nNormConst for Model fluxes = %8.3e [%8.3e]\n"%(NormConst, pars3[0])

if debug:
    print "-+-+- Limb flux arrays:"
    print "ind |   CEbins  [  Kbins  ]  |    KFluxN    |     LFlux    [  LFluxErr  ]  |  KFluxEindN  |  LFluxEind  [ LFluxErr ]"
    for j in range(len(CEbins)):
         print '%2d  | %9.4f [%9.4f]  |  %5.4e  |  %5.4e  [ %5.4e ]  |    %7.4f   |   %7.4f   [ %7.4f  ]'%(j,CEbins[j],Kbins[MatchInd[j]],KFluxN[MatchInd[j]],LFlux[j],LFluxErr[j],KFluxEindN[MatchInd[j]],LFluxEind[j],LFluxEindErr[j])

############ DRAW ############
### Fill graphs
KbinsArray=array('d',Kbins)
KbinsErrArray=array('d',KbinsErr)

KFluxNArray = array('d',KFluxN)
KFluxNErrArray = array('d',KFluxNErr)
flxModelTG = TGraphErrors(len(Kbins),KbinsArray,KFluxNArray,KbinsErrArray,KFluxNErrArray)
flxModelTG.SetName("flxModelTG")
if use_kamae:
    pltit = 'Kamae (2006) flux'
else:
    pltit = 'Kachelriess/Ostapchenko (2011) flux'
flxModelTG.SetTitle(pltit)

KFluxEindNArray = array('d',KFluxEindN)
KFluxEindNErrArray = array('d',KFluxEindNErr)
flxModelEindTG = TGraphErrors(len(Kbins),KbinsArray,KFluxEindNArray,KbinsErrArray,KFluxEindNErrArray)
flxModelEindTG.SetName("flxModelEindTG")
pltit = 'Model(2006) flux*E^%4.2f'%(SpecIndex)
if use_kamae:
    pltit = 'Kamae (2006) flux*E^%4.2f'%(SpecIndex)
else:
    pltit = 'Kachelriess/Ostapchenko (2011) flux*E^%4.2f'%(SpecIndex)
flxModelEindTG.SetTitle(pltit)

LFluxArray = array('d',LFlux)
LFluxErrArray = array('d',LFluxErr)
flxSpecTG=TGraphErrors(len(CEbins),CEbinsArray,LFluxArray,CEbinsErrArray,LFluxErrArray)
pltit='Flux #theta_{Nadir} = %4.1f-%4.1f'%(thrn[0],thrn[1])
flxSpecTG.SetTitle(pltit)
flxSpecTG.SetName("flxSpecTG")

LFluxEindArray = array('d',LFluxEind)
LFluxEindErrArray = array('d',LFluxEindErr)
flxSpecEindTG=TGraphErrors(len(CEbins),CEbinsArray,LFluxEindArray,CEbinsErrArray,LFluxEindErrArray)
pltit='Flux*E^%4.2f #theta_{Nadir} = %4.1f-%4.1f [bins %d-%d]'%(SpecIndex,thrn[0],thrn[1],thbr[0],thbr[1])
flxSpecEindTG.SetTitle(pltit)
plname='flxSpecEindTG'
flxSpecEindTG.SetName(plname)


### Init fit function
fitFlux1 = TF1( "fitFlux1", "[0]*pow(x,[1])", Ebins[0], Ebins[5])
fitFlux2 = TF1( "fitFlux2", "[0]*pow(x,[1])", Ebins[4], Ebins[10])
fitFlux1.SetParameters(3e-1,-2.8)
fitFlux2.SetParameters(3e-1,-2.8)
fitPowK1 = TF1( "fitPowK1", "[0]*pow(x,[1])", Kbins[25], Kbins[40])
fitPowK2 = TF1( "fitPowK2", "[0]*pow(x,[1])", Kbins[40], Kbins[-10])
fitPowL1 = TF1( "fitPowL1", "[0]*pow(x,[1])", Ebins[0], Ebins[5])
fitPowL2 = TF1( "fitPowL2", "[0]*pow(x,[1])*([2]+[3]*x)", Ebins[5], Ebins[11])
fitPowL2.SetParameters(100.,-2.83,0.4197,0.0004271)
fitTotalL = TF1("fitTotalL","fitPowL1+fitPowL2",Ebins[0],Ebins[11])
fitUnBiased = TF1( "fitUnBiased", "[0]*pow(x,[1])", Ebins[4], Ebins[12])


# All-Theta counts and exposure
c1=TCanvas("c1","c1 - Counts",600,600)
c1.Divide(2,2)
ndd=1
preli = []
for i in nDraw:
    c1.cd(ndd)
    c1.cd(ndd).SetRightMargin(0.13);
    cntMap[i].GetXaxis().SetTitle("#phi [deg]")
    cntMap[i].GetYaxis().SetTitle("#theta_{Nadir} [deg]")
    cntMap[i].GetYaxis().SetRangeUser(63.,80.)
    cntMap[i].Draw("COLZ")
    preli.append(TLatex(40,64,"PRELIMINARY"))
    preli[ndd-1].SetTextColor(kViolet+6)
    preli[ndd-1].Draw()
    ndd+=1
gPad.Update()
gfnam="RawCountMap_%s_%s.gif"%(escFact,irfV)
if save_plots:
    c1.SaveAs(gfnam)

c2=TCanvas("c2","c2 - Exposure",600,600)
c2.Divide(2,2)
nddd=1
for i in nDraw:
    c2.cd(nddd)
    c2.cd(nddd).SetRightMargin(0.15);
    expMap[i].GetXaxis().SetTitle("#phi [deg]")
    expMap[i].GetYaxis().SetTitle("#theta_{Nadir} [deg]")
    expMap[i].GetZaxis().SetTitle("[m^{2} s]")
    expMap[i].GetYaxis().SetRangeUser(63.,80.);
    expMap[i].Draw("COLZ")
    nddd+=1
gPad.Update()
gfnam="ExposureMap_%s.gif"%(irfV)
if save_plots:
    c2.SaveAs(gfnam)


# Limb flux over Model flux
c5=TCanvas("c5","c5 - Comparison")
c5.SetLogx()
c5.SetGridx()
c5.SetGridy()
#c5.SetLogy()
flxSpecEindTG.GetXaxis().SetTitle("E [GeV]")
flxSpecEindTG.GetYaxis().SetTitle("Flux #times E^{2.75} (m^{-2} s^{-1} sr^{-1} GeV^{1.75})")
flxSpecEindTG.SetMarkerStyle(limbMarker)
flxSpecEindTG.SetMarkerColor(myMarkerColor)
flxSpecEindTG.SetLineColor(myMarkerColor)
flxSpecEindTG.SetFillColor(0)
flxSpecEindTG.SetMarkerSize(0.6)
flxSpecEindTG.Draw("APZ");
if use_kamae:
    flxModelEindTG.SetFillColor(kamaeColor)
    flxModelEindTG.SetMarkerColor(kamaeColor)
else:
    flxModelEindTG.SetFillColor(ostapColor)
    flxModelEindTG.SetMarkerColor(ostapColor)    
flxModelEindTG.SetFillStyle(3002)
flxModelEindTG.SetLineColor(0)
flxModelEindTG.SetLineStyle(3)
flxModelEindTG.Draw("3,SAME")
leg = TLegend(0.15,0.70,0.80,0.85)
if escale:
    leg.AddEntry(flxSpecEindTG,"Gamma-ray flux from Earth Limb ESCALE [ULTRACLEANVETO]")
else:
    leg.AddEntry(flxSpecEindTG,"Gamma-ray flux from Earth Limb [ULTRACLEANVETO]")
if use_kamae:
    leg.AddEntry(flxModelEindTG,"Gamma-ray flux from Kamae (2006)")
else:
    leg.AddEntry(flxModelEindTG,"Gamma-ray flux from Kachelriess/Ostapchenko (2011)")
leg.SetTextSize(.03)
leg.SetFillColor(0)
leg.SetLineColor(0)
leg.Draw()
gPad.Update()
if use_kamae:
    gfnam="LimbFluxOnKamaeModelBand_THETA%4.1f-%4.1f_%s_%s.gif"%(thrn[0],thrn[1],escFact,irfV)
else:
    gfnam="LimbFluxOnOstapModelBand_THETA%4.1f-%4.1f_%s_%s.gif"%(thrn[0],thrn[1],escFact,irfV)
if save_plots:
    c5.SaveAs(gfnam)

# All fluxes
# Model (again)
c8=TCanvas("c8","c8 - All fluxes",1400,370)
c8.Divide(3,1)
c8.cd(1)
c8.cd(1).SetGridx()
c8.cd(1).SetGridy()
c8.cd(1).SetLogx()
c8.cd(1).SetLogy()
flxModelTG.SetMarkerStyle(modelMarker)
flxModelTG.SetMarkerSize(.4)
flxModelTG.SetMarkerColor(myMarkerColor)
flxModelTG.SetLineColor(myMarkerColor)
xaxis = flxModelTG.GetXaxis();
xaxis.SetLimits(9.,11000.);
flxModelTG.GetHistogram().SetMaximum(0.05);
flxModelTG.GetHistogram().SetMinimum(5e-12);
flxModelTG.Draw("APZ")
fitPowK1.SetLineWidth(1);
fitPowK1.SetLineColor(myLightRosso)
if debug:
    print "\n",c8.GetName(),"-+-+- Fit ",flxModelTG.GetName(), "with f(x) =",fitPowK1.GetExpFormula()
flxModelTG.Fit("fitPowK1","R")   #R
if debug:
    print "  CHISQ/NDOF= %10.4f / %4d"%(fitPowK1.GetChisquare(), fitPowK1.GetNDF())
fitPowK2.SetLineWidth(1);
fitPowK2.SetLineColor(myLightRosso)
if debug:
    print "\n",c8.GetName(),"-+-+- Fit ",flxModelTG.GetName(), "with f(x) =",fitPowK2.GetExpFormula()
flxModelTG.Fit("fitPowK2","R+")   #R+
if debug:
    print "  CHISQ/NDOF= %10.4f / %4d"%(fitPowK2.GetChisquare(), fitPowK2.GetNDF())
parsK1 = fitPowK1.GetParameters()
parsK2 = fitPowK2.GetParameters()
txt=("#splitline{C_{1} = %8.3e, #gamma_{1} = %8.3f}{C_{2} = %8.3e, #gamma_{2} = %8.3f}")%(parsK1[0],parsK1[1],parsK2[0],parsK2[1])
testoK = []
testoK.append(TLatex(120,0.001,txt))
testoK[-1].SetTextColor(myDarkRosso)
testoK[-1].SetTextSize(0.04)
testoK[-1].Draw()
gPad.Update()
# Limb (again)
c8.cd(2)
c8.cd(2).SetGridx()
c8.cd(2).SetGridy()
c8.cd(2).SetLogx()
c8.cd(2).SetLogy()
xaxis = flxSpecTG.GetXaxis();
xaxis.SetLimits(9.,11000.);
flxSpecTG.GetHistogram().SetMaximum(0.05);
flxSpecTG.GetHistogram().SetMinimum(5e-12);
flxSpecTG.SetMarkerStyle(limbMarker)
flxSpecTG.SetMarkerSize(.4)
flxSpecTG.SetMarkerColor(myMarkerColor)
flxSpecTG.SetLineColor(myMarkerColor)
flxSpecTG.Draw("APZ")
fitPowL1.SetLineWidth(1);
fitPowL1.SetLineColor(myLightBlue)
if debug:
    print "\n",c8.GetName(),"-+-+- Fit ",flxSpecTG.GetName(), "with f(x) =",fitPowL1.GetExpFormula()
flxSpecTG.Fit("fitPowL1","NR")   #R
if debug:
    print "  CHISQ/NDOF= %10.4f / %4d"%(fitPowL1.GetChisquare(), fitPowL1.GetNDF())
fitPowL2.SetLineWidth(1);
fitPowL2.SetLineColor(myLightBlue)
if debug:
    print "\n",c8.GetName(),"-+-+- Fit ",flxSpecTG.GetName(), "with f(x) =",fitPowL2.GetExpFormula()
flxSpecTG.Fit("fitPowL2","NR+")   #R+
if debug:
    print "  CHISQ/NDOF= %10.4f / %4d"%(fitPowL2.GetChisquare(), fitPowL2.GetNDF())
parsL1 = fitPowL1.GetParameters()
parsL2 = fitPowL2.GetParameters()
NPL1 = fitPowL1.GetNpar()
NPL2 = fitPowL2.GetNpar()
parsLT = []
for i in range(NPL1): 
    parsLT.append(parsL1[i])
for i in range(NPL2): 
    parsLT.append(parsL2[i])
parsLTArray=array('d',parsLT)
if debug:
    print "\n parsLT SET = ", parsLT
fitTotalL.SetParameters(parsLTArray);
fitTotalL.SetLineWidth(1);
fitTotalL.SetLineColor(myLightBlue)
if debug:
    print "\n",c8.GetName(),"-+-+- Fit ",flxSpecTG.GetName(), "with f(x) =",fitTotalL.GetExpFormula()
flxSpecTG.Fit(fitTotalL,"R+"); 
if debug:
    print "  CHISQ/NDOF= %10.4f / %4d"%(fitTotalL.GetChisquare(), fitTotalL.GetNDF())
parsLT = fitTotalL.GetParameters()
if debug:
    print "\n parsLT NEW = [", parsLT[0], parsLT[1], parsLT[2], parsLT[3], parsLT[4], parsLT[5], "]\n"
txt=("#splitline{C_{1} = %8.3e, #gamma_{1} = %8.3f}{#splitline{C_{2} = %8.3e, #gamma_{2} = %8.3f}{A = %8.3e, B = %8.3e}}")%(parsL1[0],parsL1[1],parsL2[0],parsL2[1],parsL2[2],parsL2[3])
testoL = []
testoL.append(TLatex(120,0.001,txt))
testoL[-1].SetTextColor(myDarkBlue)
testoL[-1].SetTextSize(0.04)
testoL[-1].Draw()
txt=("#splitline{Total:}{#splitline{C_{1} = %8.3e, #gamma_{1} = %8.3f}{C_{2} = %8.3e, #gamma_{2} = %8.3f}}")%(parsLT[0],parsLT[1],parsLT[2],parsLT[3])
testoL2 = []
testoL2.append(TLatex(20,0.5e-9,txt))
testoL2[-1].SetTextColor(myDarkBlue)
testoL2[-1].SetTextSize(0.04)
#testoL2[-1].Draw()
gPad.Update()

if use_glfit:
    if debug:
        print "==== UNBIAS = USING TOTAL FIT PARAMETERS ====\n"
    parsL1[0] = parsLT[0]
    parsL1[1] = parsLT[1]
    parsL2[0] = parsLT[2]
    parsL2[1] = parsLT[3]
    parsL2[2] = parsLT[4]
    parsL2[3] = parsLT[5]
else:
    if debug:
        print "==== UNBIAS = USING SEPARATE FIT PARAMETERS ====\n"


# Limb flux from "unbiased" E
c8.cd(3)
c8.cd(3).SetGridx()
c8.cd(3).SetGridy()
c8.cd(3).SetLogx()
c8.cd(3).SetLogy()
CM1 = parsK1[0]
KM1 = parsK1[1]
CM2 = parsK2[0]
KM2 = parsK2[1]
CD1 = parsL1[0]
KD1 = parsL1[1]
CD2 = parsL2[0]
KD2 = parsL2[1]
AD  = parsL2[2]
BD  = parsL2[3]

CDT1 = parsLT[0]
KDT1 = parsLT[1]
CDT2 = parsLT[2]
KDT2 = parsLT[3]
ADT  = parsLT[4]
BDT  = parsLT[5]

EMbins=[]
EMbinsErr=[]
FD=[]
FDFit=[]
FDErr=[]
if debug:
    print "-+-+- Unbiased flux arrays:"
    print "ind |    FDfit    [    FD    ]  |      ED     ===>       EMM     [    EM    ]"
for i in range(len(CEbins)):
    ED = CEbins[i]                                 #E from [D]ata
    thisFD = LFlux[i]                              #Flux for this ED (measured)
    thisFDErr = LFluxErr[i]
    if (ED<100):
        thisFDFit = CD1*pow(ED,KD1)             #Flux at this ED from fit  
        EM = pow(thisFD/CM1,1./KM1)             #E corresponding to this flux value from [M]odel
        EMM = pow(thisFDFit/CM1,1./KM1)         #E corresponding to this flux fit-value from [M]odel
        FDFit.append(thisFDFit)
        EMbins.append(EMM)
    else:                                                  
        thisFDFit = CD2*pow(ED,KD2)*(AD+BD*ED)  #PLUS Linear over 100 GeV
        EM = pow(thisFD/CM2,1./KM2)            
        EMM = pow(thisFDFit/CM2,1./KM2)
        FDFit.append(thisFDFit)
        EMbins.append(EMM)
    EMbinsErr.append(0)
    FD.append(thisFD)
    FDErr.append(thisFDErr)
    if debug:
        print '%2d  | %10.4e  [%10.4e]  | %10.4f  ===>  %10.4f   [%10.4f]'%(i,thisFDFit,thisFD,ED,EMM,EM)
if debug>1:
    print "\nEMbins created from unbias  -+-+- len(EMbins)        =", len(EMbins)
if debug>2:
    print "EMbins = ", EMbins
if debug>1:
    print "                            -+-+- len(FDFit)         =", len(FDFit)
    print "                            -+-+- len(FDErr)         =",  len(FDErr)
if debug>2:
    print "FDFit  = ",FDFit
    print "FDErr  = ",FDErr
EMbinsArray=array('d',EMbins)
EMbinsErrArray=array('d',EMbinsErr)
FDArray=array('d',FDFit)                                   #Use FIT!
FDErrArray=array('d',FDErr)
pltit='Flux Unbiased #theta_{Nadir} = %4.1f-%4.1f'%(thrn[0],thrn[1])
flxUnBiasedTG=TGraphErrors(len(EMbins),EMbinsArray,FDArray,EMbinsErrArray,FDErrArray)
flxUnBiasedTG.SetTitle(pltit)
xaxis = flxUnBiasedTG.GetXaxis();
xaxis.SetLimits(9.,11000.);
flxUnBiasedTG.GetHistogram().SetMaximum(0.05);
flxUnBiasedTG.GetHistogram().SetMinimum(5e-12);
flxUnBiasedTG.SetName("flxUnBiasedTG")
flxUnBiasedTG.SetMarkerStyle(otherMarker)
flxUnBiasedTG.SetMarkerSize(.4)
flxUnBiasedTG.SetMarkerColor(myMarkerColor)
flxUnBiasedTG.Draw("APZ")
fitUnBiased.SetLineColor(myLightGreen)
fitUnBiased.SetLineWidth(1)
if debug:
    print "\n",c8.GetName(),"-+-+- Fit ",flxUnBiasedTG.GetName(), "with f(x) =",fitUnBiased.GetExpFormula()
flxUnBiasedTG.Fit("fitUnBiased","R")   #R
parsUB = fitUnBiased.GetParameters()
if debug:
    print "  CHISQ/NDOF= %10.4f / %4d"%(fitUnBiased.GetChisquare(), fitUnBiased.GetNDF())
txt=("C = %8.3e, #gamma = %8.3f")%(parsUB[0],parsUB[1])
testoUB = []
testoUB.append(TLatex(120,0.001,txt))
testoUB[-1].SetTextColor(myDarkGreen)
testoUB[-1].SetTextSize(0.04)
testoUB[-1].Draw()
gPad.Update()
gfnam="AllFluxesWithUnbias_THETA%4.1f-%4.1f_%s_%s.gif"%(thrn[0],thrn[1],escFact,irfV)
if save_plots:
    c8.SaveAs(gfnam)


# Unbiased over Model
c9=TCanvas("c9","c9 - K + Unb")
c9.SetLogx()
c9.SetLogy()
c9.SetGridx()
c9.SetGridy()
pltit='Flux #theta_{Nadir} = %4.1f-%4.1f'%(thrn[0],thrn[1])
flxModelTG.SetTitle(pltit)
flxModelTG.SetMarkerSize(.6)
flxModelTG.SetFillColor(0)
flxModelTG.SetMarkerColor(myDarkRosso)
flxModelTG.SetLineColor(myDarkRosso)
fitPowK1.SetLineColor(0)
fitPowK2.SetLineColor(0)
flxModelTG.GetXaxis().SetTitle("E [GeV]")
flxModelTG.GetYaxis().SetTitle("Flux (m^{-2} s^{-1} sr^{-1} GeV^{-1})")
flxModelTG.Draw("APZ")
flxUnBiasedTG.SetMarkerSize(.7)
flxUnBiasedTG.SetFillColor(0)
flxUnBiasedTG.SetMarkerColor(myDarkGreen)
flxUnBiasedTG.SetLineColor(myDarkGreen)
flxUnBiasedTG.Draw("PZ,SAME")
legk = TLegend(0.4,0.75,0.85,0.85);
if use_kamae:
    legk.AddEntry(flxModelTG,"Kamae (2006)");
else:
    legk.AddEntry(flxModelTG,"Kachelriess/Ostapchenko (2011)");
legk.AddEntry(flxUnBiasedTG,"ULTRACLEANVETO Unbiased");
legk.SetTextSize(.03)
legk.SetFillColor(0)
legk.SetLineColor(0)
legk.Draw()
gPad.Update()
if use_kamae:
    gfnam="FluxUnbiasedOnKamaeModel_THETA%4.1f-%4.1f_%s_%s.gif"%(thrn[0],thrn[1],escFact,irfV)
else:
    gfnam="FluxUnbiasedOnOstapModel_THETA%4.1f-%4.1f_%s_%s.gif"%(thrn[0],thrn[1],escFact,irfV)
if save_plots:
    c9.SaveAs(gfnam)

# Unbiased over Limb
c11=TCanvas("c11","c11 - L + Unb")
c11.SetLogx()
c11.SetLogy()
c11.SetGridx()
c11.SetGridy()
pltit='Flux #theta_{Nadir} = %4.1f-%4.1f'%(thrn[0],thrn[1])
flxSpecTG.SetTitle(pltit)
flxSpecTG.SetMarkerSize(.6)
flxSpecTG.SetMarkerColor(myDarkBlue)
flxSpecTG.SetLineColor(myDarkBlue)
flxSpecTG.SetFillColor(0)
flxSpecTG.GetXaxis().SetTitle("E [GeV]")
flxSpecTG.GetYaxis().SetTitle("Flux (m^{-2} s^{-1} sr^{-1} GeV^{-1})")
flxSpecTG.Draw("APZ0")
flxUnBiasedTG.SetMarkerSize(.7)
flxUnBiasedTG.SetFillColor(0)
flxUnBiasedTG.SetMarkerColor(myDarkGreen)
flxUnBiasedTG.SetLineColor(myDarkGreen)
#flxUnBiasedTG.GetXaxis().SetRangeUser(200.,10000.)
flxUnBiasedTG.Draw("PZ0,SAME")
leg2 = TLegend(0.4,0.75,0.85,0.85); 
leg2.AddEntry(flxSpecTG,"ULTRACLEANVETO");
leg2.AddEntry(flxUnBiasedTG,"ULTRACLEANVETO Unbiased")
leg2.SetTextSize(.03)
leg2.SetFillColor(0)
leg2.SetLineColor(0)
leg2.Draw()
gPad.Update()
gfnam="FluxUnbiasedOnBiased_THETA%4.1f-%4.1f_%s_%s.gif"%(thrn[0],thrn[1],escFact,irfV)
if save_plots:
    c11.SaveAs(gfnam)

# Plot bias function b(E) (--> gives Ebias from Ereal)
c10=TCanvas("c10","c10 - Ebias",800,500)
c10.SetGridx()
c10.SetGridy()
#Bias function at E<100 GeV, use only low-energy parts of fits
bias1= TF1( "bias1", "pow([0]*[1]/[2]*pow(x,[3]),1./[4])", 10, 2000);
bias1.SetTitle("Energy BIAS function")
bias1.SetParameters(1,CM1,CD1,KM1,KD1)              #NormConst=1
bias1.SetLineWidth(1)
bias1.SetLineColor(kBlue)
bias1.Draw()
bias1.GetXaxis().SetTitle("E_{real}")
bias1.GetYaxis().SetTitle("E_{meas}")
#bias2= TF1( "bias2", "pow([0]*[1]/[2]*pow(x,[3])*(parsL2[2]+parsL2[3]*x),1./[4])", 100, 2000);
#bias2.SetParameters(1,CM2,CD2,KM2,KD2)      #NormConst=1
#bias2.SetLineWidth(1)
#bias2.SetLineColor(kRed)
#bias2.Draw("SAME")
txt1=("low-en b(E) = [C_{M1}/C_{D1}]^{1/#gamma_{D1}} E^{#gamma_{M1}/#gamma_{D1}} = %8.3f E^{%8.3f}")%(pow(1*CM1/CD1,1./KD1), KM1/KD1)
testoB1 = []
testoB1.append(TLatex(200,1800,txt1))
testoB1[-1].SetTextColor(kBlue)
testoB1[-1].SetTextSize(0.035)
testoB1[-1].Draw()
#txt2=("b(E) = %8.5f E^{%8.5f} (%8.5f + %8.5f E)^{%8.5f}")%(pow(1*CM2/CD2,1./KD2), KM2/KD2, AD, parsL2[3], 1./KD2 )
gPad.Update()
gfnam="BiasFunction_THETA%4.1f-%4.1f_%s_%s.gif"%(thrn[0],thrn[1],escFact,irfV)
if save_plots:
    c10.SaveAs(gfnam)


# Plot inverse of bias function b-1(E) (--> gives Ereal from Ebias==Emeasured)
c12=TCanvas("c12","c12 - Unbias",800,500)
c12.SetGridx()
c12.SetGridy()
### InvBias function at E<100 GeV, use only low-energy parts of fits
unbiasS1=TF1("unbiasS1","[0]*pow(x,[1])",10,2000)
unbiasS1.SetTitle("Inverse of BIAS function")
unbiasS1.SetParameters(pow(CD1/CM1,1./KM1), KD1/KM1)
unbiasS1.SetLineWidth(1)
unbiasS1.SetLineColor(myLightMagenta)
unbiasS1.GetYaxis().SetTitle("E_{real}")
unbiasS1.GetXaxis().SetTitle("E_{meas}")
unbiasS1.Draw()
txts1=("b^{-1}(E) = (C_{D1}/C_{M1})^{1/#gamma_{M1}} E^{#gamma_{D1}/#gamma_{M1}} = %8.3f E^{%8.3f}")%(pow(CD1/CM1,1./KM1),KD1/KM1)
testoUs1 = []
testoUs1.append(TLatex(60,1700,txts1))
testoUs1[-1].SetTextColor(myLightMagenta)
testoUs1[-1].SetTextSize(0.035)
testoUs1[-1].Draw()
### InvBias function at E>100 GeV, use only high-energy parts of fits
unbiasS2=TF1("unbiasS2","[0]*pow(x,[1])*pow([2]+[3]*x,[4])",10,2000)
unbiasS2.SetParameters(pow(CD2/CM2,1./KM2), KD2/KM2, AD, BD, 1./KM2)
#print "Magenta a 500 = ", unbiasS2.Eval(500)
unbiasS2.SetLineWidth(1)
unbiasS2.SetLineColor(myDarkMagenta)
unbiasS2.Draw("SAME")
txts2=("b^{-1}(E) = (C_{D2}/C_{M2})^{1/#gamma_{M2}} E^{#gamma_{D2}/#gamma_{M2}} (A_{D}+B_{D}E)^{1/#gamma_{M2}} = %8.3f E^{%8.3f} (%7.3f + %5.3f E)^{%8.3f}")%(pow(CD2/CM2,1./KM2), KD2/KM2, AD, BD, 1./KM2)
testoUs2 = []
testoUs2.append(TLatex(60,1500,txts2))
testoUs2[-1].SetTextColor(myDarkMagenta)
testoUs2[-1].SetTextSize(0.035)
testoUs2[-1].Draw()
### InvBias function at E<100 GeV, use only low-energy parts of fits FOR MODEL only, TOTAL for DATA
unbias1=TF1("unbias1","pow([0]*([1]*pow(x,[2])+[3]*pow(x,[4])*([5]+[6]*x)),[7])",10,2000)
unbias1.SetParameters(1./CM1,CDT1,KDT1,CDT2,KDT2,ADT,BDT,1./KM1)
unbias1.SetLineWidth(1)
unbias1.SetLineColor(myLightBlue)
unbias1.Draw("SAME")
txt1=("#splitline{b^{-1}(E) = (1/C_{M1})^{1/#gamma_{M1}} [C'_{D1} E^{#gamma'_{D1}} +  C'_{D2} E^{#gamma'_{D2}} (A'_{D}+B'_{D}E)^{1/#gamma_{M1}}}{         = %8.3f [ %8.3f E^{%8.3f} + %8.3f E^{%8.3f} (%7.3f + %5.3e E)]^{%8.3f}}")%(pow(1./CM1,1./KM1), CDT1, KDT1, CDT2, KDT2, ADT, BDT, 1./KM1)
testoU1 = []
testoU1.append(TLatex(60,1200,txt1))
testoU1[-1].SetTextColor(myLightBlue)
testoU1[-1].SetTextSize(0.035)
testoU1[-1].Draw()
### InvBias function at E>100 GeV, use only high-energy parts of fits FOR MODEL only, TOTAL for DATA
unbias2=TF1("unbias2","pow([0]*([1]*pow(x,[2])+[3]*pow(x,[4])*([5]+[6]*x)),[7])",10,2000)
unbias2.SetParameters(1./CM2,CDT1,KDT1,CDT2,KDT2,ADT,BDT,1./KM2)
unbias2.SetLineWidth(1)
unbias2.SetLineColor(myDarkBlue)
unbias2.Draw("SAME")
txt2=("#splitline{b^{-1}(E) = (1/C_{M2})^{1/#gamma_{M2}} [C'_{D1} E^{#gamma'_{D1}} +  C'_{D2} E^{#gamma'_{D2}} (A'_{D}+B'_{D}E)]^{1/#gamma_{M2}}}{         = %8.3f [ %8.3f E^{%8.3f} + %8.3f E^{%8.3f} (%7.3f + %5.3e E)]^{%8.3f}}")%(pow(1./CM2,1./KM2), CDT1, KDT1, CDT2, KDT2, ADT, BDT, 1./KM2)
testoU2 = []
testoU2.append(TLatex(60,900,txt2))
testoU2[-1].SetTextColor(myDarkBlue)
testoU2[-1].SetTextSize(0.035)
testoU2[-1].Draw()
gPad.Update()
gfnam="UnbiasFunctions_THETA%4.1f-%4.1f_%s_%s.gif"%(thrn[0],thrn[1],escFact,irfV)
if save_plots:
    c12.SaveAs(gfnam)

# Plot inverse of bias function b-1(E) (--> gives Ereal from Ebias==Emeasured)
c3=TCanvas("c3","c3 - Unbias func Total",800,500)
c3.SetGridx()
c3.SetGridy()
### InvBias function at E>100 GeV, use only high-energy parts of fits FOR MODEL only, TOTAL for DATA
unbias2.SetLineWidth(1)
unbias2.SetLineColor(myDarkBlue)
unbias2.SetTitle("Inverse of BIAS function")
unbias2.Draw()
testo3 = []
testo3.append(TLatex(60,1300,txt2))
testo3[-1].SetTextColor(myDarkBlue)
testo3[-1].SetTextSize(0.035)
testo3[-1].Draw()
gPad.Update()
gfnam="UnbiasFunctionHE_THETA%4.1f-%4.1f_%s_%s.gif"%(thrn[0],thrn[1],escFact,irfV)
if save_plots:
    c3.SaveAs(gfnam)


print "\n==== BIAS SUMMARY ===="
print "Background subtraction (total %d/%d events) makes a %5.3f"%(sum(bgevs), sum(sigevs),  sum(bgevs)/float(sum(sigevs))*100.),"% effect"
if escale:
    print "Using constant energy pre-scaling ",
    if escale500:
        print " 0.500"
    elif escale750:
        print " 0.750"
    else:
        print " 0.963"
else:
    print "Using NO constant energy pre-scaling"
if use_kamae:
    print "Using Kamae model from",
else:
    print "Using Ostap model from",
if use_onlyP:
    print " protons"
if use_onlyHe:
    print " Helium"
if not(use_onlyP) and not(use_onlyHe):
    print " protons + 0.085 Helium"
print "Normalization Constant: [from Ratio fit = %8.3e], used = %8.3e"%(pars3[0], NormConst)
if use_glfit:
    print "Using parameters from total fit for unbias"
else:
    print "Using parameters from separate fits for unbias"
print "LE spectral index - Model: %6.3f - DATA:  %6.3f\n)"%(parsK1[1],parsL1[1])

BT1 = pow(1./CM2,1./KM2)
BT2 = CDT1
BT3 = KDT1
BT4 = CDT2
B5 = 1./KM2
BT5 = KDT2
BT6 = ADT
BT7 = BDT
for n in [10.,50.]:
    print "Meas E = %4d ==> b-1(E) = %6.1f (bias[LE] = %6.2f"%(n,unbias1.Eval(n),(n-unbias1.Eval(n))/n*100.), "%)"
#print "\nHIGH-EN (TOTAL FIT) b-1(E') = %6.3f [%5.3f E'^{%6.3f} + %6.3f E'^{%6.3f}(%5.3f + %6.3f E')]^{%6.3f}"%(BT1,BT2,BT3,BT4,BT5,BT6,BT7,B5)
for n in [50.,100.,200.,500.,1000.]:
    print "Meas E = %4d ==> b-1(E) = %6.1f (bias[HE] = %6.2f"%(n,unbias2.Eval(n),(n-unbias2.Eval(n))/n*100.), "%)"
  
print "\nESCALEFUNCTION = '%6.3f*pow(%5.3f*pow((EvtJointEnergy/1000.),%6.3f)+%6.3f*pow((EvtJointEnergy/1000.),%6.3f)*(%5.3f+%6.3f*(EvtJointEnergy/1000.)),%6.3f)'"%(BT1,BT2,BT3,BT4,BT5,BT6,BT7,B5) 
print "ESCALEFACTOR = '%6.3f/(EvtJointEnergy/1000.)*pow(%5.3f*pow((EvtJointEnergy/1000.),%6.3f)+%6.3f*pow((EvtJointEnergy/1000.),%6.3f)*(%5.3f+%6.3f*(EvtJointEnergy/1000.)),%6.3f)'"%(BT1,BT2,BT3,BT4,BT5,BT6,BT7,B5) 
print "\nDone!"

raw_input()
