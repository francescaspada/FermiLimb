#!/usr/bin/env python

############ INIT ############

import numpy
import sys
from ROOT import *
from math import *
from array import *
from copy import *

### Analysis settings
debug       = 1
escale      = true
subtract_bg = true
save_plots  = false

Edec=[1,2,3,4]         # Energy decade (GeV) [DATA]
NEbins=[5,5,5]         # Number of bins between 10**Edec [DATA]
thbw=[400,0.,80.]      # Theta [N bins, theta1, theta2] [DATA]
#thrn=[66.,70.]         # Spectrum for this theta range [DATA]
thrn=[68.4,70.]        # Spectrum for this theta range [DATA] thin target
thBGrn=[77.,80.]       # Diff. BG from this theta range [DATA]
nDraw=[0,3,6,9,10]     # Plot maps for these i^th bin of Ebins [DATA]

### Constants
if escale: 
    NormConst=4.9e-5      # range 66-70
else:
    NormConst=5.500e-5      # range 66-70

SpecIndex  = 2.75      # spectral index
AeffScFact = 1e-3      # slope of the Aeff scaling function with E
maxESc     = 110       # max energy for Aeff scaling

### Root settings
gROOT.SetStyle('Plain')
gStyle.SetPalette(1)
gStyle.SetErrorX(0)
gStyle.SetOptFit(0)
gStyle.SetOptStat(0)

kamaeMarker = 24
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

### Data sets
datV='P301v1'
irfV='P8V6ULTRACLEANVETO' #IRF version
if escale:
    cmFile=TFile('CntMap'+irfV+'_'+datV+'test963_thetashiftpeak67.798-2.root')  # WITH E*0.963-SCALE
    escFact='0.963ESCALE'
else:
    cmFile=TFile('CntMap'+irfV+'_'+datV+'.root')             # NO E-SCALE
    escFact='NOESCALE'
emFile=TFile('ExpMap'+irfV+'_'+datV+'.root')
kaFile=open("gammaspectrum.txt")
kaFile10p=open("gammaspectrumEp.txt")
kaFile10m=open("gammaspectrumEm.txt")

### Theta bin width, bins in 0.2 deg, used bins
thw=(thbw[2]-thbw[1])/thbw[0]  
thnb=int(round(thbw[0]/(thbw[2]-thbw[1])))  ## how many bins for 1 degree
thbr=[int(round(thrn[0]*thnb+1)),int(round(thrn[1]*thnb))]
nthb=thbr[1]-thbr[0]
thBGbr=[int(round(thBGrn[0]*thnb+1)),int(round(thBGrn[1]*thnb))]
nthBGb=thBGbr[1]-thBGbr[0]

### Phi bin width 5.0 deg
phw = 5.0

if debug:
    print "SIGNAL theta range = %5.1f - %5.1f deg -+-+- bin width = %3.1f deg, used [%3d, %3d]"%(thrn[0],thrn[1],thw,thbr[0],thbr[1])
    print "BG theta range     = %5.1f - %5.1f deg -+-+- used [%3d, %3d]"%(thBGrn[0],thBGrn[1],thBGbr[0],thBGbr[1])
    print "Using NormConst = %6.5e, SpecIndex = %5.3f"%(NormConst, SpecIndex)


############ KAMAE ############

### Create energy bins for Kamae flux and get flux from AMS parameters
if debug:
    print "===== KAMAE ====="
Kbins=[]
KbinsErr=[]
KFlux=[]
KFlux10p=[]
KFlux10m=[]
KFluxErr=[]
KFluxN=[]
KFluxNErr=[]
KFluxEind=[]
KFluxEindN=[]
KFlux10pEind=[]
KFlux10mEind=[]
KFluxEindNErr=[]

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

### Create Kamae fluxes arrays 
for j in range(len(Kbins)-1):
    dEK=Kbins[j+1]-Kbins[j]
    #KbinsErr.append(dEK/2.)   # removed for the fit's sake
    KbinsErr.append(0.)
### Error on Kamae flux is |KFlux10p-KFlux10m|/2.
    KFluxN.append(NormConst*KFlux[j])
    KFluxErr.append(abs(KFlux10p[j]-KFlux10m[j])/2.)
    KFluxNErr.append(NormConst*abs(KFlux10p[j]-KFlux10m[j])/2.)
    KFluxEindN.append(NormConst*KFluxEind[j])
    KFluxEindNErr.append(NormConst*abs(KFlux10pEind[j]-KFlux10mEind[j])/2.)
    
#KbinsErr.append(Kbins[-1]-(Kbins[-2]+KbinsErr[-2]))  # removed for the fit's sake
KbinsErr.append(0.)
KFluxN.append(NormConst*KFlux[-1])
KFluxErr.append(abs(KFlux10p[-1]-KFlux10m[-1])/2.)
KFluxNErr.append(NormConst*abs(KFlux10p[-1]-KFlux10m[-1])/2.)
KFluxEindN.append(NormConst*KFluxEind[-1])
KFluxEindNErr.append(NormConst*abs(KFlux10pEind[-1]-KFlux10mEind[-1])/2.)

if debug:
    print "  len(Kbins)       =",len(Kbins)," -+-+-  len(KbinsErr) =",len(KbinsErr)
    if debug>1:
        print "KbinsErr =", KbinsErr
    print "  len(KFlux)       =",len(KFlux)," -+-+-  len(KFluxErr)  =",len(KFluxErr)
    print "  len(KFluxN)      =",len(KFluxN)," -+-+-  len(KFluxNErr) =",len(KFluxNErr)
    print "  len(KFluxEind)   =",len(KFluxEind)
    print "  len(KFluxEindN)  =",len(KFluxEindN)," -+-+-  len(KFluxEindNErr) =",len(KFluxEindNErr)

    print "-+-+- Kamae flux arrays:"
    print "ind |   Kbins   |    KFlux    [ KFluxErr ]  |   KFluxN    [ KFluxNErr]  |   KFluxEind  |  KFluxEindN  [KFluxEindErr]"
    for j in range(len(Kbins)):
        print '%2d  | %9.4f | %5.4e  [%5.4e]  | %5.4e  [%5.4e]  | %10.1f   |  %10.7f  [ %10.7f ]'%(j, Kbins[j], KFlux[j], KFluxErr[j], KFluxN[j], KFluxNErr[j], KFluxEind[j], KFluxEindN[j], KFluxEindNErr[j])


### Fill graphs for Kamae flux
KbinsArray=array('d',Kbins)
KbinsErrArray=array('d',KbinsErr)
KFluxArray = array('d',KFlux)
KFluxNArray = array('d',KFluxN)
KFluxNErrArray = array('d',KFluxNErr)
KFluxEindNArray = array('d',KFluxEindN)
KFluxErrArray = array('d',KFluxErr)
KFluxEindNErrArray = array('d',KFluxEindNErr)
pltit = 'Kamae(2006) flux'
flxKamaeTG = TGraphErrors(len(Kbins),KbinsArray,KFluxNArray,KbinsErrArray,KFluxNErrArray)
flxKamaeTG.SetTitle(pltit)
flxKamaeTG.SetName("flxKamaeTG")
pltit = 'Kamae(2006) flux*E^%4.2f'%(SpecIndex)
flxKamaeEindTG = TGraphErrors(len(Kbins),KbinsArray,KFluxEindNArray,KbinsErrArray,KFluxEindNErrArray)
flxKamaeEindTG.SetTitle(pltit)
flxKamaeEindTG.SetName("flxKamaeEindTG")


############ DATA ############

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

    
### Fill graphs for Limb flux

if debug:
    print "===== DATA ====="
    print "Maps   from file -+-+-", emFile.GetName()
    print "Counts from file -+-+-", cmFile.GetName()

EbinsArray=array('d',Ebins)

pltit='Raw Count %4.1f-%4.1f'%(thrn[0],thrn[1])
cntSpec=TH1F("cntSpec",pltit,len(Ebins)-1,EbinsArray)
pltit='Exposure %4.1f-%4.1f'%(thrn[0],thrn[1])
expSpec=TH1F("expSpec",pltit,len(Ebins)-1,EbinsArray)
pltit='Flux #theta_{Nadir} = %4.1f-%4.1f'%(thrn[0],thrn[1])
flxSpec=TH1F("flxSpec",pltit,len(Ebins)-1,EbinsArray)

if debug:
    print "  len(cntMap) =", len(cntMap)
    print "  len(expMap) =", len(expMap)
if (len(cntMap) != len(Ebins)-1):
    print "ERROR! Maps and energy bins do not correspond"
    sys.exit()
if (len(cntMap) != len(expMap)):
    print "ERROR! Count maps and exposure maps do not correspond"
    sys.exit()
    
## Raw counts and Background from diffuse emission
rawcnts=[]
sigevs=[]
bgevs=[]
nbgevs=[]
for i in range(len(cntMap)):
    if debug>1:
        print 'i=d  All Counts=%d  Bg Events=%d'%(i,cntMap[i].GetEntries(),cntMap[i].Integral(1,cntMap[i].GetNbinsX(),thBGbr[0],thBGbr[1]))
    rawcnts.append(cntMap[i].GetEntries())
    sigevs.append(cntMap[i].Integral(1,cntMap[i].GetNbinsX(),thbr[0],thbr[1]))
    bgevs.append(cntMap[i].Integral(1,cntMap[i].GetNbinsX(),thBGbr[0],thBGbr[1]))
    nbgevs.append(round(bgevs[i]/nthBGb*nthb))
if debug:
    print "  Raw counts before theta selection:   ", rawcnts
    print "  Events in signal theta range:        ", sigevs
    print "  Events in bg theta range:            ", bgevs, " in ",nthBGb," th bins and ", 360./phw, " phi bins."
    print "  ==> Bg events in signal theta range: ", nbgevs, " in ",nthb," th bins and ", 360./phw, " phi bins."

flxMap=[]
flxMapNS=[]
scExpMap=[]
invExpMap=[]
##  Differential flux map from raw counts and exposure
if debug:  print "-+-+- Scaling the Exposure:"
for i in range(len(cntMap)):
    ftmpNS=cntMap[i].Clone()
    flxMapNS.append(ftmpNS)
    flxMapNS[i].Divide(flxMapNS[i], expMap[i])
    flxMapNS[i].Scale(1e4)        #cm^2-->m^2

    ftmp=cntMap[i].Clone()
    flxMap.append(ftmp)
    scmp=expMap[i].Clone()
    scExpMap.append(scmp)
    dE=Ebins[i+1]-Ebins[i]
    EC = Ebins[i]+dE/2.
    if (EC < maxESc):
        scFact=1+AeffScFact*EC
    else:
        scFact=1+AeffScFact*maxESc
    if debug:
        print "    Bin %2d | EC = %7.2f | scFact = %f"%(i, EC, scFact)
    scExpMap[i].Scale(scFact)
    flxMap[i].Divide(flxMap[i], scExpMap[i])
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

## Solid angle
pi=3.14159 
strMap=expMap[0].Clone()
strMap.SetName("strMap")
for j in range(strMap.GetNbinsY()):
    #str=(2*pi)*(thw*pi/180.)*sin(pi*(j+0.5)*thw/180.)  #Solid angle (2pi)d(Theta)sin(Theta)  WRONG!
    str=(phw*pi/180.)*(thw*pi/180.)*sin(pi*(j+0.5)*thw/180.)  #Solid angle d(Phi)d(Theta)sin(Theta) #Corrected
    for i in range(strMap.GetNbinsX()):
        strMap.SetBinContent(i+1,j+1,str)
dOmega=strMap.Integral(1,strMap.GetNbinsX(),thbr[0],thbr[1])


CEbins = []
CEbinsErr=[]
LFlux=[]
LFluxErr=[]
LFluxEind=[]
LFluxEindNS=[]
LFluxEindErr=[]
KLRatio=[]
KLRatioErr=[]
KLDiff=[]
KLDiffErr=[]

bb=0
if debug:
    print "-+-+- Limb flux arrays:"
    print "ind |   Ebins   | CntPerBin |   CEbins  [  Kbins  ]  |     LFlux    |   KFluxN   |  KLRatio  |   KLDiff   | LFluxEind | LFluxEindNS | KFluxEindN"

for i in range(len(Ebins)-1):

## Energy bin width
    dE=Ebins[i+1]-Ebins[i]

## Raw counts    
    CntPerBin=cntMap[i].Integral(1,cntMap[i].GetNbinsX(),thbr[0],thbr[1])
    cntSpec.SetBinContent(i+1,CntPerBin/dE)
    cntSpec.SetBinError(i+1,CntPerBin**0.5/dE)
    
## Exposure
    expMap[i].Scale(1e-4)        #cm^2-->m^2    
    AvgExp=expMap[i].Integral(1,expMap[i].GetNbinsX(),thbr[0],thbr[1])
    expSpec.SetBinContent(i+1,AvgExp)

    scExpMap[i].Scale(1e-4)        #cm^2-->m^2    

## Limb DATA flux as dN/dE
    if not subtract_bg:
        thisLFlux=flxMap[i].Integral(1,flxMap[i].GetNbinsX(),thbr[0],thbr[1])/dOmega/dE
        thisLFluxNS=flxMapNS[i].Integral(1,flxMapNS[i].GetNbinsX(),thbr[0],thbr[1])/dOmega/dE
    else:
        thisLFlux=(flxMap[i].Integral(1,flxMap[i].GetNbinsX(),thbr[0],thbr[1]) - nbgevs[i]*invExpMap[i].Integral(1,invExpMap[i].GetNbinsX(),thbr[0],thbr[1]))/dOmega/dE
        thisLFluxNS=(flxMapNS[i].Integral(1,flxMapNS[i].GetNbinsX(),thbr[0],thbr[1]) - nbgevs[i]*invExpMap[i].Integral(1,invExpMap[i].GetNbinsX(),thbr[0],thbr[1]))/dOmega/dE
    thisLFluxErr=thisLFlux/(CntPerBin+0.5)**0.5
    #if (CntPerBin == 1): thisLFluxErr=0.999**0.5/dE/AvgExp                 # avoid error going to y=0
    flxSpec.SetBinContent(i+1,thisLFlux)
    flxSpec.SetBinError(i+1,thisLFluxErr)
    if debug:
        print '%2d  | %9.4f | %8d  |'%(i,Ebins[i],CntPerBin),
    if (CntPerBin == 0): print " "
    if (CntPerBin != 0):
## Array of energy bin CENTERS!
        thisEC = Ebins[i]+dE/2.
        CEbins.append(thisEC)
        #CEbinsErr.append(dE/2.)                                           # removed for the fit's sake
        CEbinsErr.append(0.)
        LFlux.append(thisLFlux)
        LFluxErr.append(thisLFluxErr)
## Limb DATA flux as E^2.75*dN/dE
        thisLFluxEind=thisEC**SpecIndex*thisLFlux
        thisLFluxEindErr=thisEC**SpecIndex*thisLFluxErr
        LFluxEind.append(thisLFluxEind)
        LFluxEindErr.append(thisLFluxEindErr)
## Limb DATA flux as E^2.75*dN/dE without AEff scale
        thisLFluxEindNS=thisEC**SpecIndex*thisLFluxNS
        LFluxEindNS.append(thisLFluxEindNS)
    
        for j in range(len(Kbins)-1):
            if (abs(CEbins[i]-Kbins[j])/CEbins[i] < 0.05):
                bb+=1
# Ratio DATA/Kamae fluxes as (dN/dE)/(dN(kamae)/dE) to get NormConst
                KLRatio.append(thisLFlux/KFlux[j])
                KLRatioErr.append(thisLFluxErr/KFlux[j])
# Difference DATA-Kamae fluxes as (dN/dE)-(NormConst*dN(kamae)/dE)
                KLDiff.append((thisLFlux-KFluxN[j])/KFluxN[j])
                KLDiffErr.append((thisLFluxErr+KFluxNErr[j])/KFluxN[j])
                if debug:
                    print '%9.4f [%9.4f]  |  %5.4e  | %5.4e | %5.3e |%11.3e | %7.5f  |   %7.5f   |  %7.5f '%(CEbins[i],Kbins[j],LFlux[i],KFluxN[j],KLRatio[i],KLDiff[i],LFluxEind[i],LFluxEindNS[i],KFluxEindN[j])

#sys.exit()


if debug:
    print "  len(CEbins) = ", len(CEbins)," -+-+- len(CEbinsErr) = ", len(CEbinsErr) 
if debug>1:
    print "CEbins = ",CEbins
    print "CEbinsErr = ",CEbinsErr
if debug:
    print "  len(LFlux)  = ", len(LFlux), " -+-+- len(LFluxErr)  = ", len(LFluxErr)
if debug>1:
    print "LFlux = ",LFlux
    print "LFluxErr = ",LFluxErr

CEbinsArray=array('d',CEbins)
CEbinsErrArray=array('d',CEbinsErr)
LFluxArray = array('d',LFlux)
LFluxErrArray = array('d',LFluxErr)
LFluxEindArray = array('d',LFluxEind)
LFluxEindNSArray = array('d',LFluxEindNS)
LFluxEindErrArray = array('d',LFluxEindErr)
KLRatioArray = array('d',KLRatio)
KLRatioErrArray = array('d',KLRatioErr)
KLDiffArray = array('d',KLDiff)
KLDiffErrArray = array('d',KLDiffErr)
pltit='Flux #theta_{Nadir} = %4.1f-%4.1f'%(thrn[0],thrn[1])
flxSpecTG=TGraphErrors(len(CEbins),CEbinsArray,LFluxArray,CEbinsErrArray,LFluxErrArray)
flxSpecTG.SetTitle(pltit)
flxSpecTG.SetName("flxSpecTG")
pltit='Flux*E^%4.2f #theta_{Nadir} = %4.1f-%4.1f'%(SpecIndex,thrn[0],thrn[1])
flxSpecEindTG=TGraphErrors(len(CEbins),CEbinsArray,LFluxEindArray,CEbinsErrArray,LFluxEindErrArray)
flxSpecEindTG.SetTitle(pltit)
flxSpecEindTG.SetName("flxSpecEindTG")
pltit='Flux*E^%4.2f #theta_{Nadir} = %4.1f-%4.1f NO AeffScale'%(SpecIndex,thrn[0],thrn[1])
flxSpecEindNSTG=TGraphErrors(len(CEbins),CEbinsArray,LFluxEindNSArray,CEbinsErrArray,LFluxEindErrArray)
flxSpecEindNSTG.SetTitle(pltit)
flxSpecEindNSTG.SetName("flxSpecEindNSTG")
pltit='Ratio Flux/Kamae #theta_{Nadir} = %4.1f-%4.1f'%(thrn[0],thrn[1])
flxRatioTG=TGraphErrors(len(CEbins),CEbinsArray,KLRatioArray,CEbinsErrArray,KLRatioErrArray)
flxRatioTG.SetTitle(pltit)
flxRatioTG.SetName("flxRatioTG")
pltit='Diff Flux-Kamae*Const #theta_{Nadir} = %4.1f-%4.1f'%(thrn[0],thrn[1])
flxDiffTG=TGraphErrors(len(CEbins),CEbinsArray,KLDiffArray,CEbinsErrArray,KLDiffErrArray)
flxDiffTG.SetTitle(pltit)
flxDiffTG.SetName("flxDiffTG")

### Init fit function

fitFlux1 = TF1( "fitFlux1", "[0]*pow(x,[1])", Ebins[0], Ebins[5])
fitFlux2 = TF1( "fitFlux2", "[0]*pow(x,[1])", Ebins[4], Ebins[10])
fitFlux1.SetParameters(3e-1,-2.8)
fitFlux2.SetParameters(3e-1,-2.8)
if escale:
    fitRatio = TF1( "fitRatio", "[0]+pow(10,-30)*x", Ebins[0], Ebins[5])  #this one for coincidence in low-en bins
else:
    #fitRatio = TF1( "fitRatio", "[0]+pow(10,-30)*x", Ebins[5], Ebins[9])  #this one for coincidence in mid-en bins
    fitRatio = TF1( "fitRatio", "[0]+pow(10,-30)*x", Ebins[2], Ebins[8])  #this one for coincidence in low-en bins
fitBias = TF1( "fitBias", "[0]+[2]*log(x)", Ebins[0], Ebins[12])
fitPowK1 = TF1( "fitPowK1", "[0]*pow(x,[1])", Kbins[25], Kbins[40])
fitPowK2 = TF1( "fitPowK2", "[0]*pow(x,[1])", Kbins[40], Kbins[-10])
fitPowL1 = TF1( "fitPowL1", "[0]*pow(x,[1])", Ebins[0], Ebins[5])
fitPowL2 = TF1( "fitPowL2", "[0]*pow(x,[1])*([2]+[3]*x)", Ebins[5], Ebins[11])
fitPowL2.SetParameters(100.,-2.83,0.4197,0.0004271)
fitTotalL = TF1("fitTotalL","fitPowL1+fitPowL2",Ebins[0],Ebins[11])
fitUnBiased = TF1( "fitUnBiased", "[0]*pow(x,[1])", Ebins[4], Ebins[12])

############ DRAW ############

# All-Theta counts and exposure
c1=TCanvas("c1","c1",600,600)
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

c2=TCanvas("c2","c2",600,600)
c2.Divide(2,2)
nddd=1
for i in nDraw:
    c2.cd(nddd)
    c2.cd(nddd).SetRightMargin(0.15);
#    expMap[i].GetXaxis().SetTitle("#phi [deg]")
#    expMap[i].GetYaxis().SetTitle("#theta_{Nadir} [deg]")
#    expMap[i].GetYaxis().SetRangeUser(63.,80.);
#    expMap[i].Draw("COLZ")
    scExpMap[i].GetXaxis().SetTitle("#phi [deg]")
    scExpMap[i].GetYaxis().SetTitle("#theta_{Nadir} [deg]")
    scExpMap[i].GetZaxis().SetTitle("[m^{2} s]")
    scExpMap[i].GetYaxis().SetRangeUser(63.,80.);
    scExpMap[i].Draw("COLZ")
    nddd+=1
gPad.Update()
gfnam="ExposureMap_%s.gif"%(irfV)
if save_plots:
    c2.SaveAs(gfnam)


# Flux
c4=TCanvas("c4")
c4.SetLogx()
c4.SetLogy()
c4.SetGridx()
c4.SetGridy()
flxSpec.GetXaxis().SetTitle("Energy (GeV)")
flxSpec.GetYaxis().SetTitle("Flux (m^{-2} s^{-1} sr^{-1} GeV^{-1})")
flxSpec.SetMarkerStyle(7)
flxSpec.Draw("E1")
fitFlux1.SetLineColor(kCyan-3);
fitFlux1.SetLineWidth(1);
if debug:
    print "\n",c4.GetName(),"-+-+- Fit ",flxSpec.GetName(), "with",fitFlux1.GetExpFormula()
#flxSpec.Fit("fitFlux1","R")
if debug:
    print "  CHISQ/NDOF= %10.4f / %4d"%(fitFlux1.GetChisquare(), fitFlux1.GetNDF())
fitFlux2.SetLineColor(kMagenta-8);
fitFlux2.SetLineWidth(1);
if debug:
    print "\n",c4.GetName(),"-+-+- Fit ",flxSpec.GetName(), "with",fitFlux2.GetExpFormula()
#flxSpec.Fit("fitFlux2","R+")
if debug:
    print "  CHISQ/NDOF= %10.4f / %4d"%(fitFlux2.GetChisquare(), fitFlux2.GetNDF())
pars = fitFlux1.GetParameters()
dpars = fitFlux1.GetParErrors()
chisq = fitFlux1.GetChisquare()
prob = fitFlux1.GetProb()
ndof = fitFlux1.GetNDF()
pars2 = fitFlux2.GetParameters()
dpars2 = fitFlux2.GetParErrors()
chisq2 = fitFlux2.GetChisquare()
prob2 = fitFlux2.GetProb()
ndof2 = fitFlux2.GetNDF()
txt=("#splitline{slope = %8.3f #pm %6.3f}{#chi^{2}/ndf = %8.2f / %d}")%(pars[1],dpars[1],chisq,ndof)
testo1 = []
testo2 = []
testo1.append(TLatex(400,0.00001,txt))
testo1[-1].SetTextColor(kCyan+2)
testo1[-1].SetTextSize(0.04)
#testo1[-1].Draw()
txt=("#splitline{slope = %8.3f #pm %6.3f}{#chi^{2}/ndf = %8.2f / %d}")%(pars2[1],dpars2[1],chisq2,ndof2)
testo2.append(TLatex(400,0.0000005,txt))
testo2[-1].SetTextColor(kMagenta-5)
testo2[-1].SetTextSize(0.04)
#testo2[-1].Draw()
gPad.Update()
gfnam="LimbFlux_THETA%4.1f-%4.1f_%s_%s.gif"%(thrn[0],thrn[1],escFact,irfV)
if save_plots:
    c4.SaveAs(gfnam)

# Ratio of flux/kamae to find normalization constant
c6=TCanvas("c6")
c6.SetLogx()
flxRatioTG.SetMarkerStyle(otherMarker)
flxRatioTG.SetMarkerSize(.7)
flxRatioTG.GetXaxis().SetTitle("E [GeV]")
#flxRatio.Draw("E1")
flxRatioTG.Draw("APZ")
fitRatio.SetLineColor(kRed-3);
fitRatio.SetLineWidth(1);
if debug:
    print "\n",c6.GetName(),"-+-+- Fit ",flxRatioTG.GetName(), "with f(x) =",fitRatio.GetExpFormula()
flxRatioTG.Fit("fitRatio","R")
if debug:
    print "  CHISQ/NDOF= %10.4f / %4d"%(fitRatio.GetChisquare(), fitRatio.GetNDF())
pars3 = fitRatio.GetParameters()
dpars3 = fitRatio.GetParErrors()
chisq3 = fitRatio.GetChisquare()
prob3 = fitRatio.GetProb()
ndof3 = fitRatio.GetNDF()
txt=("#splitline{Normalization constant = %8.3e}{#chi^{2}/ndf = %8.2f / %d}")%(pars3[0],chisq3,ndof3)
testo3 = []
testo3.append(TLatex(20,0.000001,txt))
testo3[-1].SetTextColor(kRed+2)
testo3[-1].SetTextSize(0.03)
testo3[-1].Draw()
gPad.Update()
gfnam="FluxRatio_THETA%4.1f-%4.1f_%s_%s.gif"%(thrn[0],thrn[1],escFact,irfV)
if save_plots:
    c6.SaveAs(gfnam)


# Limb flux over Kamae flux 
c5=TCanvas("c5")
c5.SetLogx()
c5.SetGridx()
c5.SetGridy()
#c5.SetLogy()
flxSpecEindTG.GetXaxis().SetTitle("E [GeV]")
flxSpecEindTG.GetYaxis().SetTitle("Flux #times E^{2.75} (m^{-2} s^{-1} sr^{-1} GeV^{1.75})")
flxSpecEindTG.SetMarkerStyle(limbMarker)
flxSpecEindTG.SetMarkerColor(myDarkMagenta)
flxSpecEindTG.SetLineColor(myDarkMagenta)
flxSpecEindTG.SetFillColor(0)
flxSpecEindTG.SetMarkerSize(0.7)
flxSpecEindTG.Draw("APZ");
flxSpecEindNSTG.SetMarkerStyle(24)
flxSpecEindNSTG.SetMarkerColor(myDarkBlue)
flxSpecEindNSTG.SetLineColor(myDarkBlue)
flxSpecEindNSTG.SetFillColor(0)
flxSpecEindNSTG.SetMarkerSize(0.7)
flxSpecEindNSTG.Draw("PZ");
flxKamaeEindTG.SetFillColor(myLightRosso)
flxKamaeEindTG.SetFillStyle(3002)
flxKamaeEindTG.SetLineColor(0)
flxKamaeEindTG.SetLineStyle(3)
flxKamaeEindTG.SetMarkerColor(myLightRosso)
flxKamaeEindTG.Draw("3,SAME")
leg = TLegend(0.15,0.70,0.80,0.85)
if escale:
    leg.AddEntry(flxSpecEindTG,"Gamma-ray flux from Earth Limb ESCALE AEFFSCALE")
    leg.AddEntry(flxSpecEindNSTG,"Gamma-ray flux from Earth Limb ESCALE")
    leg.AddEntry(flxKamaeEindTG,"Gamma-ray flux from Kamae(2006) normal. lowEn-bins")
else:
    leg.AddEntry(flxSpecEindTG,"Gamma-ray flux from Earth Limb")
    leg.AddEntry(flxKamaeEindTG,"Gamma-ray flux from Kamae(2006) normal. lowEn-bins")
#    leg.AddEntry(flxKamaeEindTG,"Gamma-ray flux from Kamae(2006) normal. midEn-bins")
leg.SetTextSize(.03)
leg.SetFillColor(0)
leg.SetLineColor(0)
leg.Draw()
gPad.Update()
gfnam="LimbFluxOnKamaeBand_THETA%4.1f-%4.1f_%s_%s.gif"%(thrn[0],thrn[1],escFact,irfV)
if save_plots:
    c5.SaveAs(gfnam)


# Diff of flux-kamae to find bias function
c7=TCanvas("c7")
c7.SetLogx()
c7.SetGridx()
c7.SetGridy()
fitBias.SetLineColor(myLightGreen);
fitBias.SetLineWidth(1);
flxDiffTG.SetMarkerStyle(otherMarker)
flxDiffTG.SetMarkerSize(.7)
flxDiffTG.SetMarkerColor(myMarkerColor)
flxDiffTG.SetLineColor(myMarkerColor)
flxDiffTG.GetXaxis().SetTitle("E [GeV]")
flxDiffTG.Draw("APZ")
if debug:
    print "\n",c7.GetName(),"-+-+- Fit ",flxDiffTG.GetName(), "with f(x) =",fitBias.GetExpFormula()
flxDiffTG.Fit("fitBias","R")
if debug:
    print "  CHISQ/NDOF= %10.4f / %4d"%(fitBias.GetChisquare(), fitBias.GetNDF())
pars4 = fitBias.GetParameters()
dpars4 = fitBias.GetParErrors()
chisq4 = fitBias.GetChisquare()
prob4 = fitBias.GetProb()
ndof4 = fitBias.GetNDF()
txt=("#splitline{#splitline{P0 = %8.3e #pm %6.3e}{P1 = %8.3e #pm %6.3e}}{#splitline{P2 = %8.3e #pm %6.3e}{#chi^{2}/ndf = %10.4f / %d}}")%(pars4[0],dpars4[0],pars4[1],dpars4[1],pars4[2],dpars4[2],chisq4,ndof4)
testo4 = []
testo4.append(TLatex(20,-.6,txt))
testo4[-1].SetTextColor(myDarkGreen)
testo4[-1].SetTextSize(0.035)
testo4[-1].Draw()
gPad.Update()
gfnam="FluxDifference_THETA%4.1f-%4.1f_%s_%s.gif"%(thrn[0],thrn[1],escFact,irfV)
if save_plots:
    c7.SaveAs(gfnam)


# All fluxes
# Kamae (again)
c8=TCanvas("c8","c8",1400,370)
c8.Divide(3,1)
c8.cd(1)
c8.cd(1).SetGridx()
c8.cd(1).SetGridy()
c8.cd(1).SetLogx()
c8.cd(1).SetLogy()
flxKamaeTG.SetMarkerStyle(kamaeMarker)
flxKamaeTG.SetMarkerSize(.4)
flxKamaeTG.SetMarkerColor(myMarkerColor)
flxKamaeTG.SetLineColor(myMarkerColor)
xaxis = flxKamaeTG.GetXaxis();
xaxis.SetLimits(9.,11000.);
flxKamaeTG.GetHistogram().SetMaximum(0.05);
flxKamaeTG.GetHistogram().SetMinimum(5e-12);
flxKamaeTG.Draw("APZ")
fitPowK1.SetLineWidth(1);
fitPowK1.SetLineColor(myLightRosso)
if debug:
    print "\n",c8.GetName(),"-+-+- Fit ",flxKamaeTG.GetName(), "with f(x) =",fitPowK1.GetExpFormula()
flxKamaeTG.Fit("fitPowK1","R")   #R
if debug:
    print "  CHISQ/NDOF= %10.4f / %4d"%(fitPowK1.GetChisquare(), fitPowK1.GetNDF())
fitPowK2.SetLineWidth(1);
fitPowK2.SetLineColor(myLightRosso)
if debug:
    print "\n",c8.GetName(),"-+-+- Fit ",flxKamaeTG.GetName(), "with f(x) =",fitPowK2.GetExpFormula()
flxKamaeTG.Fit("fitPowK2","R+")   #R+
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

if debug:
    print "==== UNBIAS = USING SEPARATE FIT PARAMETERS ===="

# UNCOMMENT THIS TO USE TOTAL FIT PARS IN THE FOLLOWING
#if debug:
#    print "==== UNBIAS = USING TOTAL FIT PARAMETERS ===="
#parsL1[0] = parsLT[0]
#parsL1[1] = parsLT[1]
#parsL2[0] = parsLT[2]
#parsL2[1] = parsLT[3]
#parsL2[2] = parsLT[4]
#parsL2[3] = parsLT[5]

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
if debug:
    print "  len(EMbins) = ", len(EMbins)
if debug>1:
    print "EMbins = ", EMbins
if debug:
    print "  len(FDFit)  = ", len(FDFit), "  -+-+-  len(FDErr) =",  len(FDErr)
if debug>1:
    print "FDFit = ",FDFit
    print "FDErr = ",FDErr
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


# Unbiased over Kamae
c9=TCanvas("c9")
c9.SetLogx()
c9.SetLogy()
c9.SetGridx()
c9.SetGridy()
pltit='Flux #theta_{Nadir} = %4.1f-%4.1f'%(thrn[0],thrn[1])
flxKamaeTG.SetTitle(pltit)
flxKamaeTG.SetMarkerSize(.6)
flxKamaeTG.SetFillColor(0)
flxKamaeTG.SetMarkerColor(myDarkRosso)
flxKamaeTG.SetLineColor(myDarkRosso)
fitPowK1.SetLineColor(0)
fitPowK2.SetLineColor(0)
flxKamaeTG.GetXaxis().SetTitle("E [GeV]")
flxKamaeTG.GetYaxis().SetTitle("Flux (m^{-2} s^{-1} sr^{-1} GeV^{-1})")
flxKamaeTG.Draw("APZ")
flxUnBiasedTG.SetMarkerSize(.7)
flxUnBiasedTG.SetFillColor(0)
flxUnBiasedTG.SetMarkerColor(myDarkGreen)
flxUnBiasedTG.SetLineColor(myDarkGreen)
flxUnBiasedTG.Draw("PZ,SAME")
legk = TLegend(0.4,0.75,0.85,0.85); 
legk.AddEntry(flxKamaeTG,"Kamae(2006)");
legk.AddEntry(flxUnBiasedTG,"ULTRACLEANVETO Unbiased");
legk.SetTextSize(.03)
legk.SetFillColor(0)
legk.SetLineColor(0)
legk.Draw()
gPad.Update()
gfnam="FluxUnbiasedOnKamae_THETA%4.1f-%4.1f_%s_%s.gif"%(thrn[0],thrn[1],escFact,irfV)
if save_plots:
    c9.SaveAs(gfnam)

# Unbiased over Limb
c11=TCanvas("c11")
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
c10=TCanvas("c10","c10",800,500)
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
c12=TCanvas("c12","c12",800,500)
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
### InvBias function at E<100 GeV, use only low-energy parts of fits FOR KAMAE only, TOTAL for DATA
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
### InvBias function at E>100 GeV, use only high-energy parts of fits FOR KAMAE only, TOTAL for DATA
unbias2=TF1("unbias2","pow([0]*([1]*pow(x,[2])+[3]*pow(x,[4])*([5]+[6]*x)),[7])",10,2000)
unbias2.SetParameters(1./CM2,CDT1,KDT1,CDT2,KDT2,ADT,BDT,1./KM2)
#print "Blu a 50 = ", unbias2.Eval(50.)
#print "Blu a 200 = ", unbias2.Eval(200.)
#print "Blu a 500 = ", unbias2.Eval(500.)
#print "Blu a 1000 = ", unbias2.Eval(1000.)
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
c3=TCanvas("c3","c3",800,500)
c3.SetGridx()
c3.SetGridy()
### InvBias function at E>100 GeV, use only high-energy parts of fits FOR KAMAE only, TOTAL for DATA
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
    print "Using constant energy pre-scaling  0.963"
else:
    print "Using NO constant energy pre-scaling"
print "Normalization Constant: [from Ratio fit =", pars3[0], "], used =", NormConst

# high energy part only
#unbiasS2=TF1("unbiasS2","[0]*pow(x,[1])*pow([2]+[3]*x,[4])",10,2000)
#unbiasS2=TF1("unbiasS2","[0]*pow(x,[1])*pow([2]+[3]*x,[5])",10,2000)
#unbiasS2.SetParameters(pow(CD2/CM2,1./KM2), KD2/KM2, AD, BD, 1./KM2)
#unbiasS2.SetParameters(pow(CD2/CM2,1./KM2), KD2/KM2, AD, BD, 1./KM2)

# TOTAL, high energy part
#unbias2=TF1("unbias2","pow([0]*([1]*pow(x,[2])+[3]*pow(x,[4])*([5]+[6]*x)),[7])",10,2000)
#unbias2.SetParameters(1./CM2,CDT1,KDT1,CDT2,KDT2,ADT,BDT,1./KM2)

B1 = pow(CD2/CM2,1./KM2)
B2 = KD2/KM2
B3 = AD
B4 = BD
B5 = 1./KM2
print "HIGH-EN b-1(E') = %6.3f E'^{%6.3f} (%6.3f+%6.3fE')^{%6.3f}"%(B1, B2, B3, B4, B5)
print "meas E =   10 ==> b-1(E) = %6.1f (bias = %6.2f"%(B1*pow(10,B2)*pow(B3+B4*10,B5), (10-(B1*pow(10,B2)*pow(B3+B4*10,B5)))/10*100), "%)"
print "meas E =  100 ==> b-1(E) = %6.1f (bias = %6.2f"%(B1*pow(100,B2)*pow(B3+B4*100,B5), (100-(B1*pow(100,B2)*pow(B3+B4*100,B5)))/100*100), "%)"
print "meas E =  500 ==> b-1(E) = %6.1f (bias = %6.2f"%(B1*pow(500,B2)*pow(B3+B4*500,B5), (500-(B1*pow(500,B2)*pow(B3+B4*500,B5)))/500*100), "%)"
print "meas E = 1000 ==> b-1(E) = %6.1f (bias = %6.2f"%(B1*pow(1000,B2)*pow(B3+B4*1000,B5), (1000-(B1*pow(1000,B2)*pow(B3+B4*1000,B5)))/1000*100), "%)"

BT1 = pow(1./CM2,1./KM2)
BT2 = CDT1
BT3 = KDT1
BT4 = CDT2
BT5 = KDT2
BT6 = ADT
BT7 = BDT
print "\nHIGH-EN (TOTAL FIT) b-1(E') = %6.3f [%5.3f E'^{%6.3f} + %6.3f E'^{%6.3f}(%5.3f + %6.3f E')]^{%6.3f}"%(BT1,BT2,BT3,BT4,BT5,BT6,BT7,B5) 
print "meas E =    10 ==> b-1(E) = %6.1f (bias = %6.2f"%(BT1*pow(BT2*pow(10,BT3)+BT4*pow(10,BT5)*(BT6+BT7*10),B5), (10-BT1*pow(BT2*pow(10,BT3)+BT4*pow(10,BT5)*(BT6+BT7*10),B5))/10.*100), "%)"
print "meas E =   100 ==> b-1(E) = %6.1f (bias = %6.2f"%(BT1*pow(BT2*pow(100,BT3)+BT4*pow(100,BT5)*(BT6+BT7*100),B5), (100-BT1*pow(BT2*pow(100,BT3)+BT4*pow(100,BT5)*(BT6+BT7*100),B5))/100.*100), "%)"
print "meas E =   500 ==> b-1(E) = %6.1f (bias = %6.2f"%(BT1*pow(BT2*pow(500,BT3)+BT4*pow(500,BT5)*(BT6+BT7*500),B5), (500-BT1*pow(BT2*pow(500,BT3)+BT4*pow(500,BT5)*(BT6+BT7*500),B5))/500.*100), "%)"
print "meas E =  1000 ==> b-1(E) = %6.1f (bias = %6.2f"%(BT1*pow(BT2*pow(1000,BT3)+BT4*pow(1000,BT5)*(BT6+BT7*1000),B5), (1000-BT1*pow(BT2*pow(1000,BT3)+BT4*pow(1000,BT5)*(BT6+BT7*1000),B5))/1000.*100), "%)"

print "\nESCALEFUNCTION = '%6.3f*pow(%5.3f*pow((EvtJointEnergy/1000.),%6.3f)+%6.3f*pow((EvtJointEnergy/1000.),%6.3f)*(%5.3f+%6.3f*(EvtJointEnergy/1000.)),%6.3f)'"%(BT1,BT2,BT3,BT4,BT5,BT6,BT7,B5) 
print "ESCALEFACTOR = '%6.3f/(EvtJointEnergy/1000.)*pow(%5.3f*pow((EvtJointEnergy/1000.),%6.3f)+%6.3f*pow((EvtJointEnergy/1000.),%6.3f)*(%5.3f+%6.3f*(EvtJointEnergy/1000.)),%6.3f)'"%(BT1,BT2,BT3,BT4,BT5,BT6,BT7,B5) 
print "\nDone!"

raw_input()
