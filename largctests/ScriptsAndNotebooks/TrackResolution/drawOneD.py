import ROOT, sys

fitdouble = False
blind = False
savefits = False

plot = sys.argv[1]

labs = plot.split("_vs_")
xlab = labs[1]
ylab = labs[0]

maxy = 0.20
if "du" in ylab: maxy = 0.1

algo = "KalmanTrack"
split = "Half10"
#
sample1 = "mc_mcc8.7_halfcuts3_full_isdata"
label1 = "BNB+Cosmics"
sample2 = "data_mcc8.8_halfcuts3_newbins"
label2 = "Data"
#
tag = "pandoraNu"+algo

file1 = ROOT.TFile.Open("plots_"+algo+split+"_"+sample1+".root")
file2 = ROOT.TFile.Open("plots_"+algo+split+"_"+sample2+".root")
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTitleYOffset(1.3)
ROOT.gStyle.SetLabelSize(ROOT.gStyle.GetLabelSize()*0.9,"xyz")

hmc = file1.Get("h12"+plot)
hdd = file2.Get("h12"+plot)

hmcp = hmc.ProjectionX().Clone()
hddp = hdd.ProjectionX().Clone()

hmcp.Sumw2()
hddp.Sumw2()

hmcp.SetLineColor(ROOT.kRed)
hddp.SetLineColor(ROOT.kMagenta)
if label2 is "Data":
    hddp.SetLineColor(ROOT.kBlack)
    hddp.SetMarkerColor(ROOT.kBlack)
    hddp.SetMarkerStyle(20)

c = ROOT.TCanvas("c","c",600,600)


xtit = xlab
if xlab=="len" or xlab=="x" or xlab=="y" or xlab=="z":
    xtit = xlab+" [cm]"

hmcp.SetTitle("")
hmcp.GetXaxis().SetTitle(xtit)

#hmcp.Sumw2()
#hddp.Sumw2()

hmcp.Scale(1./float(hmcp.GetEntries()))
hddp.Scale(1./float(hddp.GetEntries()))
hmcp.GetYaxis().SetTitle("Fraction of events")

hmcp.GetYaxis().SetRangeUser(0,1.1*max(hmcp.GetMaximum(),hddp.GetMaximum()))

c.SetGridy()
hmcp.Draw("E1")
hddp.Draw("E1,same")

leg = ROOT.TLegend(0.1,0.9,0.9,0.99)
leg.SetNColumns(4)
leg.AddEntry(hmcp,label1,"L")
leg.AddEntry(hddp,label2,"L")
leg.Draw()

c.SaveAs("one_"+xlab+"_"+tag+split+".png")
