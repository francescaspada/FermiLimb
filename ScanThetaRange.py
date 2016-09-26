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
save_plots = false

Edec=[1,2,3,4]         # Energy decade (GeV) [DATA]
NEbins=[5,5,5]         # Number of bins between 10**Edec [DATA]
thbw=[400,0.,80.]      # Theta [N bins, theta1, theta2] [DATA]
#thrn=[66.,70.]         # Spectrum for this theta range [DATA]
thrn=[68.4,70.]        # Spectrum for this theta range [DATA] thin target
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
KFlux10p=[]
KFlux10m=[]
KFluxErr=[]
KFluxN=[]
KFluxNErr=[]
KFluxEind=[]
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
    print "ind |   Kbins    |   KFluxEind  "
    for j in range(len(Kbins)):
        print '%2d  | %9.4f  | %10.1f   '%(j, Kbins[j], KFluxEind[j])



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

flxKamaeEindNTG=[]
flxSpecEindTG=[]
flxRatioTG=[]

## Raw counts and Background from diffuse emission  ===  SCAN in theta range
Nsteps = 5
c4=TCanvas("c4","c4",500,600)
c4.SetLogx()
c5=TCanvas("c5","c5",1200,700)
c5.Divide(2,int(Nsteps+1)/2)
c6=TCanvas("c6","c6",1200,700)
c6.Divide(2,int(Nsteps+1)/2)

leg5=[]
leg4 = TLegend(0.15,0.70,0.80,0.85)


for step in range(Nsteps):
    thrnst1=(thrn[0]+0.5*step)            # from thrn[0] in steps of 0. degs, width 1 deg
    thrnst2=(thrn[1]-(0.5*((Nsteps-1)-step)))
    if debug:
        print "\nSCAN step %d === range %5.1f - %5.1f deg"%(step,thrnst1,thrnst2)

    thbr=[int(round(thrnst1*thnb+1)),int(round(thrnst2*thnb))]
    nthb=thrnst2-thrnst1
    if debug:
        print "  SIGNAL theta range = %5.1f - %5.1f deg -+-+- bin width = %3.1f deg, used [%3d, %3d]"%(thrnst1,thrnst2,thw,thbr[0],thbr[1])
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
    KLRatioErr=[]
## Array of energy bin CENTERS!
    CEbins = []
    CEbinsErr=[]

    bb=0
    if debug:
        print "  -+-+- Limb flux arrays:"
        print "  ind |   Ebins   | CntPerBin |     CEbins  [  Kbins  ]  |  KLRatio   |  LFluxEind | KFluxEindN"

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
                    KLRatioErr.append(thisLFluxErr/KFlux[j])
                    if debug:
                        print '  %9.4f [%9.4f]  | %5.3e  |  %7.5f  '%(CEbins[i],Kbins[j],KLRatio[i],LFluxEind[i])

#sys.exit()

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
    KLRatioErrArray = array('d',KLRatioErr)
    _flxSpecEindTG=TGraphErrors(len(CEbins),CEbinsArray,LFluxEindArray,CEbinsErrArray,LFluxEindErrArray)
    pltit='Flux*E^%4.2f #theta_{Nadir} = %4.1f-%4.1f [bins %d-%d]'%(SpecIndex,thrnst1,thrnst2,thbr[0],thbr[1])
    _flxSpecEindTG.SetTitle(pltit)
    plname='flxSpecEindTG%d'%step
    _flxSpecEindTG.SetName(plname)
    _flxRatioTG=TGraphErrors(len(CEbins),CEbinsArray,KLRatioArray,CEbinsErrArray,KLRatioErrArray)
    pltit='Ratio Flux/Kamae #theta_{Nadir} = %4.1f-%4.1f [bins %d-%d]'%(thrnst1,thrnst2,thbr[0],thbr[1])
    _flxRatioTG.SetTitle(pltit)
    plname='flxRatioTG%d'%step
    _flxRatioTG.SetName(plname)

### Init fit function

    fitRatio = TF1( "fitRatio", "[0]+pow(10,-30)*x", Ebins[0], Ebins[4])  #this one for coincidence in low-en bins

############ DRAW ############

# Ratio of flux/kamae to find normalization constant
    c6.cd(step+1)
    c6.cd(step+1).SetLogx()
    flxRatioTG.append(_flxRatioTG)
    flxRatioTG[step].SetMarkerStyle(otherMarker)
    flxRatioTG[step].SetMarkerSize(.6)
    flxRatioTG[step].GetXaxis().SetTitle("E [GeV]")
#flxRatio.Draw("E1")
    flxRatioTG[step].Draw("APZ")
    fitRatio.SetLineColor(kRed-3);
    fitRatio.SetLineWidth(1);
    if debug:
        print "\n ",c6.GetName()," -+-+- Fit ",_flxRatioTG.GetName(), "with f(x) =",fitRatio.GetExpFormula()
        print ' ',
    flxRatioTG[step].Fit("fitRatio","R")
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

    NormConst = pars3[0]
    if (step!=0):
        NormConst += 0.1*NormConst
    
### Create Kamae normalized flux arrays
    KFluxEindN=[]
    KFluxEindNErr=[]
    for j in range(len(Kbins)-1):
### Error on Kamae flux is |KFlux10p-KFlux10m|/2.
        KFluxEindN.append(NormConst*KFluxEind[j])
        KFluxEindNErr.append(NormConst*abs(KFlux10pEind[j]-KFlux10mEind[j])/2.)    
    KFluxEindN.append(NormConst*KFluxEind[-1])
    KFluxEindNErr.append(NormConst*abs(KFlux10pEind[-1]-KFlux10mEind[-1])/2.)

    if debug:
        print "  len(KFluxEindN)  =",len(KFluxEindN)," -+-+-  len(KFluxEindNErr) =",len(KFluxEindNErr)

        print "  -+-+- Kamae normalized flux arrays:"
        print "  ind |   Kbins    |  KFluxEindN  [KFluxEindErr]"
        for j in range(len(Kbins)):
            print '  %2d  | %9.4f  |  %10.7f  [ %10.7f ]'%(j, Kbins[j], KFluxEindN[j], KFluxEindNErr[j])

### Fill normalized graphs for Kamae flux
    KFluxEindNArray = array('d',KFluxEindN)
    KFluxEindNErrArray = array('d',KFluxEindNErr)
    _flxKamaeEindNTG = TGraphErrors(len(Kbins),KbinsArray,KFluxEindNArray,KbinsErrArray,KFluxEindNErrArray)
    pltit = 'Kamae(2006) flux*E^%4.2f'%(SpecIndex)
    _flxKamaeEindNTG.SetTitle(pltit)
    plname='flxKamaeEindNTG%d'%step
    _flxKamaeEindNTG.SetName(plname)


# Limb flux over Kamae flux
    c5.cd(step+1)
    c5.cd(step+1).SetLogx()
    c5.cd(step+1).SetGridx()
    c5.cd(step+1).SetGridy()
#c5.SetLogy()
    flxSpecEindTG.append(_flxSpecEindTG)
    flxSpecEindTG[step].GetXaxis().SetTitle("E [GeV]")
    flxSpecEindTG[step].GetYaxis().SetTitle("Flux #times E^{2.75} (m^{-2} s^{-1} sr^{-1} GeV^{1.75})")
    flxSpecEindTG[step].SetMarkerStyle(limbMarker)
    flxSpecEindTG[step].SetMarkerColor(myLightBlue+step)
    flxSpecEindTG[step].SetLineColor(myLightBlue+step)
    flxSpecEindTG[step].SetFillColor(0)
    flxSpecEindTG[step].SetMarkerSize(0.6)
    flxSpecEindTG[step].Draw("APZ");
    flxKamaeEindNTG.append(_flxKamaeEindNTG)
    flxKamaeEindNTG[step].SetFillColor(myLightRosso)
    flxKamaeEindNTG[step].SetFillStyle(3002)
    flxKamaeEindNTG[step].SetLineColor(0)
    flxKamaeEindNTG[step].SetLineStyle(3)
    flxKamaeEindNTG[step].SetMarkerColor(myLightRosso)
    flxKamaeEindNTG[step].Draw("3,SAME")
    _leg = TLegend(0.15,0.70,0.80,0.85)
    if escale:
        _leg.AddEntry(flxSpecEindTG[step],"Gamma-ray flux from Earth Limb ESCALE [ULTRACLEANVETO]")
        _leg.AddEntry(flxKamaeEindNTG[step],"Gamma-ray flux from Kamae(2006) normal. lowEn-bins")
    else:
        _leg.AddEntry(flxSpecEindTG[step],"Gamma-ray flux from Earth Limb [ULTRACLEANVETO]")
        _leg.AddEntry(flxKamaeEindNTG[step],"Gamma-ray flux from Kamae(2006) normal. lowEn-bins")
    leg5.append(_leg)
    leg5[step].SetTextSize(.03)
    leg5[step].SetFillColor(0)
    leg5[step].SetLineColor(0)
    leg5[step].Draw()
    gPad.Update()
    gfnam="LimbFluxOnKamaeBand_THETA%4.1f-%4.1f_%s_%s.gif"%(thrn[0],thrn[1],escFact,irfV)
    if save_plots:
        c5.SaveAs(gfnam)
        

    c4.cd()
    if (step==0):
        flxSpecEindTG[step].GetYaxis().SetRangeUser(0.,180.)
        flxSpecEindTG[step].Draw("APZ");
    else:
        flxSpecEindTG[step].Draw("PZ");
    _text='#theta_{Nadir} = %4.1f-%4.1f [bins %d-%d]'%(thrnst1,thrnst2,thbr[0],thbr[1])
    leg4.AddEntry(flxSpecEindTG[step],_text)
    if (step==Nsteps-1):
        leg4.SetTextSize(.03)
        leg4.SetFillColor(0)
        leg4.SetLineColor(0)
        leg4.Draw()

    gPad.Update()

    print "  ==== BIAS SUMMARY ===="
    print "  Background subtraction (total %d/%d events) makes a %5.3f"%(sum(bgevs), sum(sigevs),  sum(bgevs)/float(sum(sigevs))*100.),"% effect"
    print "  Normalization Constant: [from Ratio fit =", pars3[0], "], used =", NormConst


    if escale:
        print "  Using constant energy pre-scaling 0.963"
    else:
        print "Using NO constant energy pre-scaling"
        
print "\nDone!"

raw_input()
