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
sample1 = "bnb_cosmics_mcc8.4"
#sample1 = "bnb_mcc8.4"
sample2 = "data_mcc8.4" #sample1 to check biases if mc cuts
#sample2 = sample1 #sample1 to check biases if mc cuts
#sample2 = "bnb_cosmics_mcc8.4"
tag = "pandoraNu"+algo #+"_singlebin"

file1 = ROOT.TFile.Open("plots_"+algo+split+"_"+sample1+"_"+tag+".root")
file2 = ROOT.TFile.Open("plots_"+algo+split+"_"+sample2+"_"+tag+".root")
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTitleYOffset(1.3)
ROOT.gStyle.SetLabelSize(ROOT.gStyle.GetLabelSize()*0.9,"xyz")

hmc = file1.Get("h12"+plot)
h1m = file1.Get("h1m"+plot)
htm = file1.Get("htm"+plot)
hdd = file2.Get("h12"+plot)

hmcp = ROOT.TH1F("hmcp","hmcp",hmc.GetNbinsX(),hmc.GetXaxis().GetXmin(),hmc.GetXaxis().GetXmax())
hmcp.SetLineColor(ROOT.kRed)

h1mp = ROOT.TH1F("h1mp","h1mp",h1m.GetNbinsX(),h1m.GetXaxis().GetXmin(),h1m.GetXaxis().GetXmax())
h1mp.SetLineColor(ROOT.kBlue)

htmp = ROOT.TH1F("htmp","htmp",htm.GetNbinsX(),htm.GetXaxis().GetXmin(),htm.GetXaxis().GetXmax())
htmp.SetLineColor(ROOT.kBlack)

hddp = ROOT.TH1F("hddp","hddp",hdd.GetNbinsX(),hdd.GetXaxis().GetXmin(),hdd.GetXaxis().GetXmax())
hddp.SetLineColor(ROOT.kMagenta)

c = ROOT.TCanvas("c","c",600,600)

for i in range (0,hmc.GetNbinsX()):
    hmcpi = hmc.ProjectionY( ("bin%d" % (i+1)),i+1,i+1)
    hmcpi.Rebin(4)
    if hmcpi.GetEntries() < 10: continue
    hmcpi.Draw()
    if fitdouble:
        r1 = hmcpi.Fit("gaus","S","",-0.2,0.2)
        r2 = hmcpi.Fit("gaus","S","",-0.4,0.4)
        f1 = ROOT.TF1("fitFunc","gaus(0)+gaus(3)")
        f1.SetParameter(0,r1.Parameter(0))
        f1.SetParameter(1,r1.Parameter(1))
        f1.SetParameter(2,r1.Parameter(2))
        f1.SetParameter(3,r2.Parameter(0))
        f1.SetParameter(4,r2.Parameter(1))
        f1.SetParameter(5,r2.Parameter(2))
        hmcpi.Fit(f1,"","",-0.4,0.4)
        if not f1 or f1.IsValid() == False: continue
        fa = ROOT.TF1("fa","gaus",-0.2,0.2)
        fa.SetParameter(0,f1.GetParameter(0))
        fa.SetParameter(1,f1.GetParameter(1))
        fa.SetParameter(2,f1.GetParameter(2))
        fa.SetParError(0,f1.GetParError(0))
        fa.SetParError(1,f1.GetParError(1))
        fa.SetParError(2,f1.GetParError(2))
        fa.Draw("same")
        fb = ROOT.TF1("fb","gaus",-0.4,0.4)
        fb.SetParameter(0,f1.GetParameter(3))
        fb.SetParameter(1,f1.GetParameter(4))
        fb.SetParameter(2,f1.GetParameter(5))
        fb.SetParError(0,f1.GetParError(3))
        fb.SetParError(1,f1.GetParError(4))
        fb.SetParError(2,f1.GetParError(5))
        fb.Draw("same")
        if savefits: c.SaveAs("mcbin"+str(i+1)+"_"+tag+"_"+split+"_"+plot+".png")
        if fb and fb.IsValid():
            hmcp.SetBinContent(i+1,abs(fb.GetParameter(2)))
            hmcp.SetBinError(i+1,abs(fb.GetParError(2)))
    else:
        hmcpi.Fit("gaus","","",-0.2,0.2)
        if savefits: c.SaveAs("mcbin"+str(i+1)+"_"+tag+"_"+split+"_"+plot+".png")
        fr = hmcpi.GetFunction("gaus")
        if fr and fr.IsValid():
            hmcp.SetBinContent(i+1,abs(fr.GetParameter(2)))
            hmcp.SetBinError(i+1,abs(fr.GetParError(2)))
for i in range (0,h1m.GetNbinsX()):
    if not blind: continue
    h1mpi = h1m.ProjectionY( ("bin%d" % (i+1)),i+1,i+1)
    h1mpi.Rebin(4)
    if h1mpi.GetEntries() < 10: continue
    h1mpi.Draw()
    if fitdouble:
        r1 = h1mpi.Fit("gaus","S","",-0.2,0.2)
        r2 = h1mpi.Fit("gaus","S","",-0.4,0.4)
        f1 = ROOT.TF1("fitFunc","gaus(0)+gaus(3)")
        f1.SetParameter(0,r1.Parameter(0))
        f1.SetParameter(1,r1.Parameter(1))
        f1.SetParameter(2,r1.Parameter(2))
        f1.SetParameter(3,r2.Parameter(0))
        f1.SetParameter(4,r2.Parameter(1))
        f1.SetParameter(5,r2.Parameter(2))
        h1mpi.Fit(f1,"","",-0.4,0.4)
        if not f1 or f1.IsValid() == False: continue
        fa = ROOT.TF1("fa","gaus",-0.2,0.2)
        fa.SetParameter(0,f1.GetParameter(0))
        fa.SetParameter(1,f1.GetParameter(1))
        fa.SetParameter(2,f1.GetParameter(2))
        fa.SetParError(0,f1.GetParError(0))
        fa.SetParError(1,f1.GetParError(1))
        fa.SetParError(2,f1.GetParError(2))
        fa.Draw("same")
        fb = ROOT.TF1("fb","gaus",-0.4,0.4)
        fb.SetParameter(0,f1.GetParameter(3))
        fb.SetParameter(1,f1.GetParameter(4))
        fb.SetParameter(2,f1.GetParameter(5))
        fb.SetParError(0,f1.GetParError(3))
        fb.SetParError(1,f1.GetParError(4))
        fb.SetParError(2,f1.GetParError(5))
        fb.Draw("same")
        #if savefits: c.SaveAs("1mbin"+str(i+1)+"_"+tag+"_"+split+"_"+plot+".png")
        if fb and fb.IsValid():
            h1mp.SetBinContent(i+1,abs(fb.GetParameter(2)))
            h1mp.SetBinError(i+1,abs(fb.GetParError(2)))
    else:
        h1mpi.Fit("gaus","","",-0.2,0.2)
        #if savefits: c.SaveAs("1mbin"+str(i+1)+"_"+tag+"_"+split+"_"+plot+".png")
        fr = h1mpi.GetFunction("gaus")
        if fr and fr.IsValid():
            h1mp.SetBinContent(i+1,abs(fr.GetParameter(2)))
            h1mp.SetBinError(i+1,abs(fr.GetParError(2)))
for i in range (0,htm.GetNbinsX()):
    if not blind: continue
    htmpi = htm.ProjectionY( ("bin%d" % (i+1)),i+1,i+1)
    htmpi.Rebin(4)
    if htmpi.GetEntries() < 10: continue
    htmpi.Draw()
    if fitdouble:
        r1 = htmpi.Fit("gaus","S","",-0.2,0.2)
        r2 = htmpi.Fit("gaus","S","",-0.4,0.4)
        f1 = ROOT.TF1("fitFunc","gaus(0)+gaus(3)")
        f1.SetParameter(0,r1.Parameter(0))
        f1.SetParameter(1,r1.Parameter(1))
        f1.SetParameter(2,r1.Parameter(2))
        f1.SetParameter(3,r2.Parameter(0))
        f1.SetParameter(4,r2.Parameter(1))
        f1.SetParameter(5,r2.Parameter(2))
        htmpi.Fit(f1,"","",-0.4,0.4)
        if not f1 or f1.IsValid() == False: continue
        fa = ROOT.TF1("fa","gaus",-0.2,0.2)
        fa.SetParameter(0,f1.GetParameter(0))
        fa.SetParameter(1,f1.GetParameter(1))
        fa.SetParameter(2,f1.GetParameter(2))
        fa.SetParError(0,f1.GetParError(0))
        fa.SetParError(1,f1.GetParError(1))
        fa.SetParError(2,f1.GetParError(2))
        fa.Draw("same")
        fb = ROOT.TF1("fb","gaus",-0.4,0.4)
        fb.SetParameter(0,f1.GetParameter(3))
        fb.SetParameter(1,f1.GetParameter(4))
        fb.SetParameter(2,f1.GetParameter(5))
        fb.SetParError(0,f1.GetParError(3))
        fb.SetParError(1,f1.GetParError(4))
        fb.SetParError(2,f1.GetParError(5))
        fb.Draw("same")
        #if savefits: c.SaveAs("tmbin"+str(i+1)+"_"+tag+"_"+split+"_"+plot+".png")
        if fb and fb.IsValid():
            htmp.SetBinContent(i+1,abs(fb.GetParameter(2)))
            htmp.SetBinError(i+1,abs(fb.GetParError(2)))
    else:
        htmpi.Fit("gaus","","",-0.2,0.2)
        #if savefits: c.SaveAs("tmbin"+str(i+1)+"_"+tag+"_"+split+"_"+plot+".png")
        fr = htmpi.GetFunction("gaus")
        if fr and fr.IsValid():
            htmp.SetBinContent(i+1,abs(fr.GetParameter(2)))
            htmp.SetBinError(i+1,abs(fr.GetParError(2)))
for i in range (0,hdd.GetNbinsX()):
    if blind: continue
    hddpi = hdd.ProjectionY( ("bin%d" % (i+1)),i+1,i+1)
    hddpi.Rebin(4)
    if hddpi.GetEntries() < 10: continue
    hddpi.Draw()
    if fitdouble:
        r1 = hddpi.Fit("gaus","S","",-0.2,0.2)
        r2 = hddpi.Fit("gaus","S","",-0.4,0.4)
        f1 = ROOT.TF1("fitFunc","gaus(0)+gaus(3)")
        f1.SetParameter(0,r1.Parameter(0))
        f1.SetParameter(1,r1.Parameter(1))
        f1.SetParameter(2,r1.Parameter(2))
        f1.SetParameter(3,r2.Parameter(0))
        f1.SetParameter(4,r2.Parameter(1))
        f1.SetParameter(5,r2.Parameter(2))
        hddpi.Fit(f1,"","",-0.4,0.4)
        if not f1 or f1.IsValid() == False: continue
        fa = ROOT.TF1("fa","gaus",-0.2,0.2)
        fa.SetParameter(0,f1.GetParameter(0))
        fa.SetParameter(1,f1.GetParameter(1))
        fa.SetParameter(2,f1.GetParameter(2))
        fa.SetParError(0,f1.GetParError(0))
        fa.SetParError(1,f1.GetParError(1))
        fa.SetParError(2,f1.GetParError(2))
        fa.Draw("same")
        fb = ROOT.TF1("fb","gaus",-0.4,0.4)
        fb.SetParameter(0,f1.GetParameter(3))
        fb.SetParameter(1,f1.GetParameter(4))
        fb.SetParameter(2,f1.GetParameter(5))
        fb.SetParError(0,f1.GetParError(3))
        fb.SetParError(1,f1.GetParError(4))
        fb.SetParError(2,f1.GetParError(5))
        fb.Draw("same")
        if savefits: c.SaveAs("ddbin"+str(i+1)+"_"+tag+"_"+split+"_"+plot+".png")
        if fb and fb.IsValid():
            hddp.SetBinContent(i+1,abs(fb.GetParameter(2)))
            hddp.SetBinError(i+1,abs(fb.GetParError(2)))
    else:
        hddpi.Fit("gaus","","",-0.2,0.2)
        if savefits: c.SaveAs("ddbin"+str(i+1)+"_"+tag+"_"+split+"_"+plot+".png")
        fr = hddpi.GetFunction("gaus")
        if fr and fr.IsValid():
            hddp.SetBinContent(i+1,abs(fr.GetParameter(2)))
            hddp.SetBinError(i+1,abs(fr.GetParError(2)))

if ylab=="dx" or ylab=="dy" or ylab=="dz":
    ylab = ylab+" [cm]"

ylab = ylab.replace("d","#sigma")

if xlab=="len" or xlab=="x" or xlab=="y" or xlab=="z":
    xlab = xlab+" [cm]"

hmcp.SetTitle("")
hmcp.GetXaxis().SetTitle(xlab)
hmcp.GetYaxis().SetTitle(ylab)
hmcp.GetYaxis().SetRangeUser(0,maxy)

c.SetGridy()
hmcp.Draw("E1")
h1mp.Draw("E1,same")
htmp.Draw("E1,same")
if not blind: hddp.Draw("E1,same")

leg = ROOT.TLegend(0.1,0.9,0.9,0.99)
leg.SetNColumns(4)
leg.AddEntry(hmcp,"split 1,2 MC","L")#bnb+cosmics
if blind:
    leg.AddEntry(h1mp,"split MC","L")
    leg.AddEntry(htmp,"full MC","L")
if not blind: leg.AddEntry(hddp,"split 1,2 data","L")
leg.Draw()

c.SaveAs("data_"+tag+"_"+split+"_"+plot+".png")

if not blind:
    hmcp2 = hmcp.Clone()
    hmcp2.Divide(hddp)
    for bin in range(1,hmcp2.GetNbinsX()+1):
        if hmcp2.GetBinError(bin)>0.5*hmcp2.GetBinContent(bin): 
            hmcp2.SetBinContent(bin,0)
            hmcp2.SetBinError(bin,0)
    hmcp2.SetLineColor(ROOT.kBlack)
    hmcp2.GetYaxis().SetRangeUser(0.,2.)
    hmcp2.Draw("E1")
    hmcp2.GetXaxis().SetTitle(xlab)
    hmcp2.GetYaxis().SetTitle("MC/data")
    c.SetGridy()
    c.SaveAs("ratio_"+tag+"_"+split+"_"+plot+".png")

ROOT.gStyle.SetOptStat(110010)
