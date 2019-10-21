import ROOT, sys

savefits = False

samples = []
methods = []
labels = []
colors = []

#samples.append("prodgenie_bnb_nu_cosmic_uboone_mcc9.0_pre_oct_reco_1d")
#methods.append("12")
#labels.append("Tag0 MC")
#colors.append(ROOT.kMagenta)

#samples.append("data_bnb_rolling_unblind_mcc9_pre_oct_reco_1d")
#methods.append("12")
#labels.append("Tag0 rolling data")
#colors.append(ROOT.kRed)

#samples.append("data_bnb_run1_unblind_mcc9_pre_oct_reco_1d")
#methods.append("12")
#labels.append("Tag0 run1 data")
#colors.append(ROOT.kCyan+1)

#samples.append("prodgenie_bnb_nu_cosmic_uboone_mcc9.0_beta1_oct_reco_2d_dic")
#methods.append("12")
#labels.append("Tag1 MC")
#colors.append(ROOT.kBlue)

samples.append("data_bnb_rolling_unblind_mcc9_beta1_oct_reco_2d")
methods.append("12")
labels.append("Tag1 rolling data")
colors.append(ROOT.kBlack)

samples.append("data_bnb_run1_unblind_mcc9_beta1_oct_reco_2d")
methods.append("12")
labels.append("Tag1 run1 data")
colors.append(ROOT.kGray+1)

#samples.append("prodgenie_bnb_nu_cosmic_uboone_mcc9.0_beta2_oct_reco_2d_wc")
#methods.append("12")
#labels.append("Tag2 MC")
#colors.append(ROOT.kGreen+2)

sample1 = samples[0]
sample2 = samples[1]
method1 = methods[0]
method2 = methods[1]
label1 = labels[0]
label2 = labels[1]
color1 = colors[0]
color2 = colors[1]

plot = sys.argv[1]

labs = plot.split("_vs_")
xlab = labs[1]
ylab = labs[0]

minentries = 30
nrebinx = 1
nrebiny = 2
#if  "mommumcs" in xlab: nrebinx = 2

maxy = 0.3
if "dp0" in ylab: maxy = 0.5
if "dp1" in ylab: maxy = 0.5
if "dp2" in ylab: maxy = 0.15
if "dp3" in ylab: maxy = 0.15
if "dx" in ylab: maxy = 0.5
if "dy" in ylab: maxy = 0.5
if "dz" in ylab: maxy = 0.5
if "du" in ylab: maxy = 0.15
if "dphi" in ylab: maxy = 0.2
if "dtheta" in ylab: maxy = 0.15
if ylab=="dx" or ylab=="dy" or ylab=="dz" or ylab=="dp0" or ylab=="dp1":
    ylab = ylab+" [cm]"
    
ylab = ylab.replace("d","#sigma")
ylab = ylab.replace("phi","#phi")
ylab = ylab.replace("theta","#theta")
ylab = ylab.replace("p0","p_{0}")
ylab = ylab.replace("p1","p_{1}")
ylab = ylab.replace("p2","p_{2}")
ylab = ylab.replace("p3","p_{3}")
if "theta" in ylab or "phi" in ylab:
    ylab = ylab+" [rad]"
if xlab=="len":
    xlab = "Track length [cm]"
if xlab=="x" or xlab=="y" or xlab=="z":
    xlab = xlab+" [cm]"
if xlab=="mommumcs":
    xlab = "p_{MCS}^{#mu} [GeV]"

algo = "KalmanTrack"
split = "Half10"
tag = "pandoraNu"+algo

file1 = ROOT.TFile.Open("plots_"+algo+split+"_"+sample1+".root")
file2 = ROOT.TFile.Open("plots_"+algo+split+"_"+sample2+".root")
ROOT.gStyle.SetOptStat(111110)
ROOT.gStyle.SetTitleYOffset(1.3)
ROOT.gStyle.SetLabelSize(ROOT.gStyle.GetLabelSize()*0.9,"xyz")

h1 = file1.Get("h"+method1+plot)
h2 = file2.Get("h"+method2+plot)
h1.RebinX(nrebinx)
h2.RebinX(nrebinx)

if len(h1.GetXaxis().GetXbins())>0:
    h1p = ROOT.TH1F("h1p","h1p",h1.GetNbinsX(),h1.GetXaxis().GetXbins().GetArray())
else:
    h1p = ROOT.TH1F("h1p","h1p",h1.GetNbinsX(),h1.GetXaxis().GetXmin(),h1.GetXaxis().GetXmax())
h1p.SetLineColor(color1)
#h1p.SetFillColor(color1)
#h1p.SetLineWidth(2)
h1p.SetMarkerColor(color1)
h1p.SetMarkerStyle(21)
h1p.SetTitle(label1)
if len(h2.GetXaxis().GetXbins())>0:
    h2p = ROOT.TH1F("h2p","h2p",h2.GetNbinsX(),h2.GetXaxis().GetXbins().GetArray())
else:
    h2p = ROOT.TH1F("h2p","h2p",h2.GetNbinsX(),h2.GetXaxis().GetXmin(),h2.GetXaxis().GetXmax())
h2p.SetLineColor(color2)
#h2p.SetLineWidth(2)
h2p.SetTitle(label2)
h2p.SetMarkerColor(color2)
h2p.SetMarkerStyle(20)

label = ROOT.TText()
label.SetNDC()
label.SetTextFont(1)
label.SetTextColor(1)
label.SetTextSize(0.06)
label.SetTextAlign(11)
label.SetTextAngle(0)
    
ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gStyle.SetLegendFont(43)  # Absolute font size in pixel (precision 3)
ROOT.gStyle.SetLegendTextSize(25)
ROOT.gStyle.SetOptTitle(0)
c = ROOT.TCanvas("c","c",600,600)

for i in range (0,h1.GetNbinsX()):
    h1pi = h1.ProjectionY( ("bin%d" % (i+1)),i+1,i+1)
    if method1 is "12":
        h1pi.GetXaxis().SetTitle(ylab.replace("sigma","sqrt{1/2}#upoint#Delta"))
    else:
        h1pi.GetXaxis().SetTitle(ylab.replace("sigma","Delta"))
    if i == 0:
        h1pi.Rebin(nrebiny)
    elif i<3:
        h1pi.Rebin(nrebiny)
    else:
        h1pi.Rebin(nrebiny)
    if h1pi.GetEntries() < minentries: continue
    h1pi.Draw()
    h1pi.Fit("gaus","MLE","",h1pi.GetMean()-h1pi.GetRMS()*1.5,h1pi.GetMean()+h1pi.GetRMS()*1.5)
    #h1pi.Fit("gaus","MLE","",h1pi.GetMean()-h1pi.GetRMS()*2,h1pi.GetMean()+h1pi.GetRMS()*2)
    if savefits:
        c.SaveAs("mcbin"+str(i+1)+"_"+tag+"_"+split+"_"+plot+".png")
    fr = h1pi.GetFunction("gaus")
    if fr and fr.IsValid():
        #h1p.SetBinContent(i+1,h1pi.GetRMS())
        #h1p.SetBinError(i+1,h1pi.GetRMSError())
        h1p.SetBinContent(i+1,abs(fr.GetParameter(2)))
        h1p.SetBinError(i+1,abs(fr.GetParError(2)))
for i in range (0,h2.GetNbinsX()):
    h2pi = h2.ProjectionY( ("bin%d" % (i+1)),i+1,i+1)
    if method1 is "12":
        h2pi.GetXaxis().SetTitle(ylab.replace("sigma","sqrt{1/2}#upoint#Delta"))
    else:
        h2pi.GetXaxis().SetTitle(ylab.replace("sigma","Delta"))
    if i == 0:
        h2pi.Rebin(nrebiny)
    elif i<3:
        h2pi.Rebin(nrebiny)
    else:
        h2pi.Rebin(nrebiny)
    if h2pi.GetEntries() < minentries: continue
    h2pi.Draw()
    h2pi.Fit("gaus","MLE","",h2pi.GetMean()-h2pi.GetRMS()*1.5,h2pi.GetMean()+h2pi.GetRMS()*1.5)
    #h2pi.Fit("gaus","MLE","",h2pi.GetMean()-h2pi.GetRMS()*2,h2pi.GetMean()+h2pi.GetRMS()*2)
    if savefits: c.SaveAs("ddbin"+str(i+1)+"_"+tag+"_"+split+"_"+plot+".png")
    fr = h2pi.GetFunction("gaus")
    if fr and fr.IsValid():
        #h2p.SetBinContent(i+1,h2pi.GetRMS())
        #h2p.SetBinError(i+1,h2pi.GetRMSError())
        h2p.SetBinContent(i+1,abs(fr.GetParameter(2)))
        h2p.SetBinError(i+1,abs(fr.GetParError(2)))

c.Clear()
ROOT.gStyle.SetOptStat(0)

pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
pad1.SetBottomMargin(0.025) # Leave some room so that y axis labels are not truncated
pad1.SetGridx()
pad1.SetGridy()
pad1.Draw()
pad1.cd()

h1p.GetXaxis().SetTitle(xlab)
h1p.GetYaxis().SetTitle(ylab)
h1p.GetYaxis().SetRangeUser(0,maxy)

h1p.GetYaxis().SetTitleSize(25)
h1p.GetYaxis().SetTitleFont(43)
h1p.GetYaxis().SetTitleOffset(1.12)
h1p.GetYaxis().SetLabelFont(43)
h1p.GetYaxis().SetLabelSize(18)
# Do not draw the X axis label and title
h1p.GetXaxis().SetLabelSize(0.);
h1p.GetXaxis().SetTitleSize(0.);

h1p.Draw("E1")
h2p.Draw("E1,same")

pad1.BuildLegend(0.581,0.552,0.895,0.892,"Exiting Tracks","")

#label.DrawText(0.1, 0.93, "MicroBooNE Preliminary")
#label.DrawText(0.68, 0.93,  "1.62e20 POT")
#label.DrawText(0.1, 0.93, "MicroBooNE Simulation Preliminary")

c.cd() # Go back to the main canvas before defining pad2
pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
pad2.SetTopMargin(0.025)
pad2.SetBottomMargin(0.38)
pad2.SetGridx()
pad2.SetGridy()
pad2.Draw()
pad2.cd()

# Define the ratio plot
h3 = h1p.Clone("h3")
h3.Divide(h2p)
for bin in range(1,h3.GetNbinsX()+1):
    if h3.GetBinError(bin)>0.5*h3.GetBinContent(bin): 
        h3.SetBinContent(bin,0)
        h3.SetBinError(bin,0)
h3.SetLineColor(ROOT.kBlack)
h3.GetYaxis().SetRangeUser(0.5,1.5)
h3.SetMarkerStyle(0)
h3.GetXaxis().SetTitle(xlab)
h3.GetYaxis().SetTitle("ratio")
h3.SetTitle("")
h3.Draw("EP1")

# Y axis ratio plot settings
h3.GetYaxis().SetNdivisions(505)
h3.GetYaxis().SetTitleSize(25)
h3.GetYaxis().SetTitleFont(43)
h3.GetYaxis().SetTitleOffset(1.12)
h3.GetYaxis().SetLabelFont(43)
h3.GetYaxis().SetLabelSize(16)

# X axis ratio plot settings
h3.GetXaxis().SetTitleSize(25)
h3.GetXaxis().SetTitleFont(43)
h3.GetXaxis().SetTitleOffset(3.5)
h3.GetXaxis().SetLabelFont(43)
h3.GetXaxis().SetLabelSize(20)

c.SaveAs("data_"+sample1+"_"+method1+"_"+sample2+"_"+method2+"_"+plot+".png")

ROOT.gStyle.SetOptStat(110010)
