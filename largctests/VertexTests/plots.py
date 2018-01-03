import ROOT as r

f = r.TFile.Open("/uboone/data/users/cerati/fittedVertexNtuple_bnb_pandoraNuKalmanTrack_v06_48_00_v2.root","READ")

t1 = f.Get("fittedvertexntuplizer/tree")
t2 = f.Get("fittedvertexntuplizerscaled/tree")

var = "y"
xtitle = "#Delta"+var+" [cm]"
den = ""
xarange = 2

doPull = False
if doPull:
    den = "/sqrt(vtx_fit_c"+var+var+")"
    xtitle = "#Delta"+var+"/#sigma"+var
    xarange = 10

#cut = "vtx_fit_tkdists_max<0"
cut = "(vtx_fit_chi2/vtx_fit_ndof)>10."

h11 = r.TH1F("h11","h11",100,-xarange,xarange)
h11.SetLineColor(r.kBlue)
h12 = r.TH1F("h12","h12",100,-xarange,xarange)
h12.SetLineColor(r.kRed)
h21 = r.TH1F("h21","h21",100,-xarange,xarange)
h21.SetLineColor(r.kBlack)

t1.Draw("(vtx_fit_"+var+" - vtx_gen_"+var+")"+den+">>h11",cut,"goff")
t2.Draw("(vtx_fit_"+var+" - vtx_gen_"+var+")"+den+">>h12",cut,"goff")
t1.Draw("(vtx_pan_"+var+" - vtx_gen_"+var+")>>h21",cut,"goff")

c = r.TCanvas("c","c",1200,1200)
c.Divide(2,2)

c.cd(4)
h11c = h11.Clone()
h12c = h12.Clone()
h21c = h21.Clone()
h11c.SetTitle("comparison")
h11c.Scale(1./h11c.Integral())
h12c.Scale(1./h12c.Integral())
h21c.Scale(1./h21c.Integral())
h11c.GetXaxis().SetTitle(xtitle)
h11c.GetYaxis().SetRangeUser(0,0.25)
h11c.Draw("hist")
h12c.Draw("hist,same")
h21c.Draw("hist,same")

c.cd(1)
h11.GetXaxis().SetTitle(xtitle)
h11.SetTitle("pandoraNuKalmanTrack")
h11.Draw()
h11.Fit("gaus")#,"","",-0.5,0.5

c.cd(2)
h12.GetXaxis().SetTitle(xtitle)
h12.SetTitle("pandoraNuKalmanTrack - scaled #sigma_{hit}")
h12.Draw()
h12.Fit("gaus")#,"","",-0.5,0.5

c.cd(3)
h21.GetXaxis().SetTitle(xtitle)
h21.SetTitle("pandoraNu")
h21.Draw()
h21.Fit("gaus")#,"","",-0.5,0.5
