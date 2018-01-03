import ROOT, sys

plot = sys.argv[1]

labs = plot.split("_vs_")
xlab = labs[1]
ylab = labs[0]

maxy = 0.2
if "du" in ylab: maxy = 0.1

algo = "KalmanTrack"
#algo = "PMA"
split = "Half10"
#split = "Half15"
sample = "bnb_cosmics_mcc8.4"
tag = "pandoraNu"+algo+"_cutsv1"
#tag = tag+"_scale001"

file = ROOT.TFile.Open("plots_"+algo+split+"_"+sample+"_"+tag+".root")
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTitleYOffset(1.3)
ROOT.gStyle.SetLabelSize(ROOT.gStyle.GetLabelSize()*0.9,"xyz")

h12 = file.Get("h12"+plot)
h12.FitSlicesY()
h12p = ROOT.gDirectory.Get("h12"+plot+"_2")
h12p.SetLineColor(ROOT.kRed)

h1m = file.Get("h1m"+plot)
h1m.FitSlicesY()
h1mp = ROOT.gDirectory.Get("h1m"+plot+"_2")
h1mp.SetLineColor(ROOT.kBlue)

htm = file.Get("htm"+plot)
htm.FitSlicesY()
htmp = ROOT.gDirectory.Get("htm"+plot+"_2")
htmp.SetLineColor(ROOT.kBlack)

c = ROOT.TCanvas("c","c",600,600)

if ylab=="dx" or ylab=="dy" or ylab=="dz":
    ylab = ylab+" [cm]"

ylab = ylab.replace("d","#sigma")

if xlab=="len" or xlab=="x" or xlab=="y" or xlab=="z":
    xlab = xlab+" [cm]"

h12p.SetTitle("")
h12p.GetXaxis().SetTitle(xlab)
h12p.GetYaxis().SetTitle(ylab)
h12p.GetYaxis().SetRangeUser(0,maxy)

h12p.Draw()
h1mp.Draw("same")
htmp.Draw("same")

leg = ROOT.TLegend(0.1,0.9,0.9,0.99)
leg.SetNColumns(3)
leg.AddEntry(h12p,"split track 1,2","L")
leg.AddEntry(h1mp,"split track MC","L")
leg.AddEntry(htmp,"full track MC","L")
leg.Draw()

c.SaveAs(tag+"_"+split+"_"+plot+".png")
