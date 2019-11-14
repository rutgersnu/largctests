import ROOT
import math

alg = "KalmanTrack"
algnoscale = alg.replace("_scale001","")

tag = "_data_bnb_optfilter_C1_5e19_goodruns_v08_00_00_12_reco2"
tag = "_prodgenie_bnb_nu_uboone_overlay_mcc9_v12"
tag = "_prodgenie_bnb_nu_only_mcc9_v12_reco2"

isdata = True
if "data" in tag: isdata = True
#isdata = True

fnme = algnoscale+"Half10"

file = ROOT.TFile.Open("/Users/cerati/Notebooks/mcc9.0-files/allstudies"+tag+".root")
#ROOT.gStyle.SetOptStat(0)

print fnme
print "Comparison at midpoint!"
print 'isdata=', isdata    

tree = file.Get("splitTrackNtuplizer"+fnme+"/tree")
#tree = file.Get("tree")

tk_gpar = ROOT.vector('float')()
tree.SetBranchAddress("tk_vtx_gpar", tk_gpar)

tk1_gpar_vtx = ROOT.vector('float')()
tree.SetBranchAddress("tk1_vtx_gpar", tk1_gpar_vtx)
tk1_gpar_end = ROOT.vector('float')()
tree.SetBranchAddress("tk1_end_gpar", tk1_gpar_end)
tk1_tk_gpar_vtx = ROOT.vector('float')()
tree.SetBranchAddress("tk1_vtxtk_gpar", tk1_tk_gpar_vtx)
tk1_tk_gpar_mid = ROOT.vector('float')()
tree.SetBranchAddress("tk1_midtk_gpar", tk1_tk_gpar_mid)
tk1_tk_gpar_end = ROOT.vector('float')()
tree.SetBranchAddress("tk1_endtk_gpar", tk1_tk_gpar_end)

tk2_gpar_vtx = ROOT.vector('float')()
tree.SetBranchAddress("tk2_vtx_gpar", tk2_gpar_vtx)
tk2_gpar_end = ROOT.vector('float')()
tree.SetBranchAddress("tk2_end_gpar", tk2_gpar_end)
tk2_tk1_gpar_vtx = ROOT.vector('float')()
tree.SetBranchAddress("tk2_vtxtk1_gpar", tk2_tk1_gpar_vtx)
tk2_tk1_gpar_end = ROOT.vector('float')()
tree.SetBranchAddress("tk2_endtk1_gpar", tk2_tk1_gpar_end)
tk2_tk_gpar_vtx = ROOT.vector('float')()
tree.SetBranchAddress("tk2_vtxtk_gpar", tk2_tk_gpar_vtx)
tk2_tk_gpar_mid = ROOT.vector('float')()
tree.SetBranchAddress("tk2_midtk_gpar", tk2_tk_gpar_mid)
tk2_tk_gpar_end = ROOT.vector('float')()
tree.SetBranchAddress("tk2_endtk_gpar", tk2_tk_gpar_end)

mc_gpar_vtx = ROOT.vector('float')()
tree.SetBranchAddress("mc_vtx_gpar", mc_gpar_vtx)
mc_tk_gpar_vtx = ROOT.vector('float')()
tree.SetBranchAddress("mc_vtxtk_gpar", mc_tk_gpar_vtx)
mc_tk1_gpar_vtx = ROOT.vector('float')()
tree.SetBranchAddress("mc_vtxtk1_gpar", mc_tk1_gpar_vtx)
mc_tk2_gpar_vtx = ROOT.vector('float')()
tree.SetBranchAddress("mc_vtxtk2_gpar", mc_tk2_gpar_vtx)
mc_tk_gpar_mid = ROOT.vector('float')()
tree.SetBranchAddress("mc_midtk_gpar", mc_tk_gpar_mid)
mc_tk_gpar_end = ROOT.vector('float')()
tree.SetBranchAddress("mc_endtk_gpar", mc_tk_gpar_end)
mc_tk1_gpar_end = ROOT.vector('float')()
tree.SetBranchAddress("mc_endtk1_gpar", mc_tk1_gpar_end)
mc_tk2_gpar_end = ROOT.vector('float')()
tree.SetBranchAddress("mc_endtk2_gpar", mc_tk2_gpar_end)

tk1_mc_gpar_mid = ROOT.vector('float')()
tree.SetBranchAddress("tk1_midmc_gpar", tk1_mc_gpar_mid)
tk2_mc_gpar_mid = ROOT.vector('float')()
tree.SetBranchAddress("tk2_midmc_gpar", tk2_mc_gpar_mid)
tk_mc_gpar_vtx = ROOT.vector('float')()
tree.SetBranchAddress("tk_vtxmc_gpar", tk_mc_gpar_vtx)
tk_mc_gpar_end = ROOT.vector('float')()
tree.SetBranchAddress("tk_endmc_gpar", tk_mc_gpar_end)

vtxpos = ROOT.vector('float')()
tree.SetBranchAddress("vx_pos", vtxpos)

fout = ROOT.TFile("plots_"+fnme+tag+".root","RECREATE")

# split track resolution (data-driven)
h12dx = ROOT.TH1F("h12dx","h12dx",100,-0.5,0.5)
h12dy = ROOT.TH1F("h12dy","h12dy",100,-0.5,0.5)
h12dz = ROOT.TH1F("h12dz","h12dz",100,-0.5,0.5)
h12dux = ROOT.TH1F("h12dux","h12dux",100,-0.1,0.1)
h12duy = ROOT.TH1F("h12duy","h12duy",100,-0.1,0.1)
h12duz = ROOT.TH1F("h12duz","h12duz",100,-0.1,0.1)
h12duxr = ROOT.TH1F("h12duxr","h12duxr",100,-0.25,0.25)
h12duyr = ROOT.TH1F("h12duyr","h12duyr",100,-0.25,0.25)
h12duzr = ROOT.TH1F("h12duzr","h12duzr",100,-0.25,0.25)
h12dp = [h12dx, h12dy, h12dz, h12dux, h12duy, h12duz, h12duxr, h12duyr, h12duzr]
for h in h12dp: h.SetLineColor(ROOT.kRed)

# split track resolution (MC)
h1mdx = ROOT.TH1F("h1mdx","h1mdx",100,-0.5,0.5)
h1mdy = ROOT.TH1F("h1mdy","h1mdy",100,-0.5,0.5)
h1mdz = ROOT.TH1F("h1mdz","h1mdz",100,-0.5,0.5)
h1mdux = ROOT.TH1F("h1mdux","h1mdux",100,-0.1,0.1)
h1mduy = ROOT.TH1F("h1mduy","h1mduy",100,-0.1,0.1)
h1mduz = ROOT.TH1F("h1mduz","h1mduz",100,-0.1,0.1)
h1mduxr = ROOT.TH1F("h1mduxr","h1mduxr",100,-0.25,0.25)
h1mduyr = ROOT.TH1F("h1mduyr","h1mduyr",100,-0.25,0.25)
h1mduzr = ROOT.TH1F("h1mduzr","h1mduzr",100,-0.25,0.25)
h1mdp = [h1mdx, h1mdy, h1mdz, h1mdux, h1mduy, h1mduz, h1mduxr, h1mduyr, h1mduzr]
for h in h1mdp: h.SetLineColor(ROOT.kBlue)

# full track resolution
htmdx = ROOT.TH1F("htmdx","htmdx",100,-0.5,0.5)
htmdy = ROOT.TH1F("htmdy","htmdy",100,-0.5,0.5)
htmdz = ROOT.TH1F("htmdz","htmdz",100,-0.5,0.5)
htmdux = ROOT.TH1F("htmdux","htmdux",100,-0.1,0.1)
htmduy = ROOT.TH1F("htmduy","htmduy",100,-0.1,0.1)
htmduz = ROOT.TH1F("htmduz","htmduz",100,-0.1,0.1)
htmduxr = ROOT.TH1F("htmduxr","htmduxr",100,-0.25,0.25)
htmduyr = ROOT.TH1F("htmduyr","htmduyr",100,-0.25,0.25)
htmduzr = ROOT.TH1F("htmduzr","htmduzr",100,-0.25,0.25)
htmdp = [htmdx, htmdy, htmdz, htmdux, htmduy, htmduz, htmduxr, htmduyr, htmduzr]
for h in htmdp: h.SetLineColor(ROOT.kBlack)

testcase = ["12","1m","tm"]
if isdata: testcase = ["12"]
print 'testcase=', testcase

h_dp_vs_len = []
h_dp_vs_nh = []
h_dp_vs_x = []
h_dp_vs_y = []
h_dp_vs_z = []
h_dp_vs_ux = []
h_dp_vs_uy = []
h_dp_vs_uz = []
h_dp_vs_mcmom = []
h_dp_vs_coszx = []
h_dp_vs_coszy = []
h_dp_vs_ismu = []
h_dp_vs_ntkvtx = []
h_dp_vs_contain = []
h_dp_vs_mommumcs = []
h_dp_vs_mommurng = []
h_dp_vs_mompmcs = []
h_dp_vs_momprng = []
for tc in testcase:
    #print tc
    # track resolution vs track lenght
    dx_vs_len = ROOT.TH2F("h"+tc+"dx_vs_len","h"+tc+"dx_vs_len",10,0,500,100,-0.5,0.5)
    dy_vs_len = ROOT.TH2F("h"+tc+"dy_vs_len","h"+tc+"dy_vs_len",10,0,500,100,-0.5,0.5)
    dz_vs_len = ROOT.TH2F("h"+tc+"dz_vs_len","h"+tc+"dz_vs_len",10,0,500,100,-0.5,0.5)
    dux_vs_len = ROOT.TH2F("h"+tc+"dux_vs_len","h"+tc+"dux_vs_len",10,0,500,100,-0.1,0.1)
    duy_vs_len = ROOT.TH2F("h"+tc+"duy_vs_len","h"+tc+"duy_vs_len",10,0,500,100,-0.1,0.1)
    duz_vs_len = ROOT.TH2F("h"+tc+"duz_vs_len","h"+tc+"duz_vs_len",10,0,500,100,-0.1,0.1)
    dr_vs_len = ROOT.TH2F("h"+tc+"dr_vs_len","h"+tc+"dr_vs_len",10,0,500,100,-0.5,0.5)
    dur_vs_len = ROOT.TH2F("h"+tc+"dur_vs_len","h"+tc+"dur_vs_len",10,0,500,100,-0.1,0.1)
    dphi_vs_len = ROOT.TH2F("h"+tc+"dphi_vs_len","h"+tc+"dphi_vs_len",10,0,500,100,-0.3,0.3)
    dtheta_vs_len = ROOT.TH2F("h"+tc+"dtheta_vs_len","h"+tc+"dtheta_vs_len",10,0,500,100,-0.2,0.2)
    dp_vs_len = [dx_vs_len, dy_vs_len, dz_vs_len, dux_vs_len, duy_vs_len, duz_vs_len, dr_vs_len, dur_vs_len, dphi_vs_len, dtheta_vs_len]
    h_dp_vs_len.append(dp_vs_len)
    # track resolution vs track nhits
    dx_vs_nh = ROOT.TH2F("h"+tc+"dx_vs_nh","h"+tc+"dx_vs_nh",10,0,2000,100,-0.5,0.5)
    dy_vs_nh = ROOT.TH2F("h"+tc+"dy_vs_nh","h"+tc+"dy_vs_nh",10,0,2000,100,-0.5,0.5)
    dz_vs_nh = ROOT.TH2F("h"+tc+"dz_vs_nh","h"+tc+"dz_vs_nh",10,0,2000,100,-0.5,0.5)
    dux_vs_nh = ROOT.TH2F("h"+tc+"dux_vs_nh","h"+tc+"dux_vs_nh",10,0,2000,100,-0.1,0.1)
    duy_vs_nh = ROOT.TH2F("h"+tc+"duy_vs_nh","h"+tc+"duy_vs_nh",10,0,2000,100,-0.1,0.1)
    duz_vs_nh = ROOT.TH2F("h"+tc+"duz_vs_nh","h"+tc+"duz_vs_nh",10,0,2000,100,-0.1,0.1)
    dr_vs_nh = ROOT.TH2F("h"+tc+"dr_vs_nh","h"+tc+"dr_vs_nh",10,0,2000,100,-0.5,0.5)
    dur_vs_nh = ROOT.TH2F("h"+tc+"dur_vs_nh","h"+tc+"dur_vs_nh",10,0,2000,100,-0.1,0.1)
    dphi_vs_nh = ROOT.TH2F("h"+tc+"dphi_vs_nh","h"+tc+"dphi_vs_nh",10,0,2000,100,-0.3,0.3)
    dtheta_vs_nh = ROOT.TH2F("h"+tc+"dtheta_vs_nh","h"+tc+"dtheta_vs_nh",10,0,2000,100,-0.2,0.2)
    dp_vs_nh = [dx_vs_nh, dy_vs_nh, dz_vs_nh, dux_vs_nh, duy_vs_nh, duz_vs_nh, dr_vs_nh, dur_vs_nh, dphi_vs_nh, dtheta_vs_nh]
    h_dp_vs_nh.append(dp_vs_nh)
    # track resolution vs track x (MC)
    dx_vs_x = ROOT.TH2F("h"+tc+"dx_vs_x","h"+tc+"dx_vs_x",10,0,250,100,-0.5,0.5)
    dy_vs_x = ROOT.TH2F("h"+tc+"dy_vs_x","h"+tc+"dy_vs_x",10,0,250,100,-0.5,0.5)
    dz_vs_x = ROOT.TH2F("h"+tc+"dz_vs_x","h"+tc+"dz_vs_x",10,0,250,100,-0.5,0.5)
    dux_vs_x = ROOT.TH2F("h"+tc+"dux_vs_x","h"+tc+"dux_vs_x",10,0,250,100,-0.1,0.1)
    duy_vs_x = ROOT.TH2F("h"+tc+"duy_vs_x","h"+tc+"duy_vs_x",10,0,250,100,-0.1,0.1)
    duz_vs_x = ROOT.TH2F("h"+tc+"duz_vs_x","h"+tc+"duz_vs_x",10,0,250,100,-0.1,0.1)
    dr_vs_x = ROOT.TH2F("h"+tc+"dr_vs_x","h"+tc+"dr_vs_x",10,0,250,100,-0.5,0.5)
    dur_vs_x = ROOT.TH2F("h"+tc+"dur_vs_x","h"+tc+"dur_vs_x",10,0,250,100,-0.1,0.1)
    dphi_vs_x = ROOT.TH2F("h"+tc+"dphi_vs_x","h"+tc+"dphi_vs_x",10,0,250,100,-0.3,0.3)
    dtheta_vs_x = ROOT.TH2F("h"+tc+"dtheta_vs_x","h"+tc+"dtheta_vs_x",10,0,250,100,-0.2,0.2)
    dp_vs_x = [dx_vs_x, dy_vs_x, dz_vs_x, dux_vs_x, duy_vs_x, duz_vs_x, dr_vs_x, dur_vs_x, dphi_vs_x, dtheta_vs_x]
    h_dp_vs_x.append(dp_vs_x)
    # track resolution vs track y (MC)
    dx_vs_y = ROOT.TH2F("h"+tc+"dx_vs_y","h"+tc+"dx_vs_y",10,-125,125,100,-0.5,0.5)
    dy_vs_y = ROOT.TH2F("h"+tc+"dy_vs_y","h"+tc+"dy_vs_y",10,-125,125,100,-0.5,0.5)
    dz_vs_y = ROOT.TH2F("h"+tc+"dz_vs_y","h"+tc+"dz_vs_y",10,-125,125,100,-0.5,0.5)
    dux_vs_y = ROOT.TH2F("h"+tc+"dux_vs_y","h"+tc+"dux_vs_y",10,-125,125,100,-0.1,0.1)
    duy_vs_y = ROOT.TH2F("h"+tc+"duy_vs_y","h"+tc+"duy_vs_y",10,-125,125,100,-0.1,0.1)
    duz_vs_y = ROOT.TH2F("h"+tc+"duz_vs_y","h"+tc+"duz_vs_y",10,-125,125,100,-0.1,0.1)
    dr_vs_y = ROOT.TH2F("h"+tc+"dr_vs_y","h"+tc+"dr_vs_y",10,-125,125,100,-0.5,0.5)
    dur_vs_y = ROOT.TH2F("h"+tc+"dur_vs_y","h"+tc+"dur_vs_y",10,-125,125,100,-0.1,0.1)
    dphi_vs_y = ROOT.TH2F("h"+tc+"dphi_vs_y","h"+tc+"dphi_vs_y",10,-125,125,100,-0.3,0.3)
    dtheta_vs_y = ROOT.TH2F("h"+tc+"dtheta_vs_y","h"+tc+"dtheta_vs_y",10,-125,125,100,-0.2,0.2)
    dp_vs_y = [dx_vs_y, dy_vs_y, dz_vs_y, dux_vs_y, duy_vs_y, duz_vs_y, dr_vs_y, dur_vs_y, dphi_vs_y, dtheta_vs_y]
    h_dp_vs_y.append(dp_vs_y)
    # track resolution vs track z (MC)
    dx_vs_z = ROOT.TH2F("h"+tc+"dx_vs_z","h"+tc+"dx_vs_z",10,0,1100,100,-0.5,0.5)
    dy_vs_z = ROOT.TH2F("h"+tc+"dy_vs_z","h"+tc+"dy_vs_z",10,0,1100,100,-0.5,0.5)
    dz_vs_z = ROOT.TH2F("h"+tc+"dz_vs_z","h"+tc+"dz_vs_z",10,0,1100,100,-0.5,0.5)
    dux_vs_z = ROOT.TH2F("h"+tc+"dux_vs_z","h"+tc+"dux_vs_z",10,0,1100,100,-0.1,0.1)
    duy_vs_z = ROOT.TH2F("h"+tc+"duy_vs_z","h"+tc+"duy_vs_z",10,0,1100,100,-0.1,0.1)
    duz_vs_z = ROOT.TH2F("h"+tc+"duz_vs_z","h"+tc+"duz_vs_z",10,0,1100,100,-0.1,0.1)
    dr_vs_z = ROOT.TH2F("h"+tc+"dr_vs_z","h"+tc+"dr_vs_z",10,0,1100,100,-0.5,0.5)
    dur_vs_z = ROOT.TH2F("h"+tc+"dur_vs_z","h"+tc+"dur_vs_z",10,0,1100,100,-0.1,0.1)
    dphi_vs_z = ROOT.TH2F("h"+tc+"dphi_vs_z","h"+tc+"dphi_vs_z",10,0,1100,100,-0.3,0.3)
    dtheta_vs_z = ROOT.TH2F("h"+tc+"dtheta_vs_z","h"+tc+"dtheta_vs_z",10,0,1100,100,-0.2,0.2)
    dp_vs_z = [dx_vs_z, dy_vs_z, dz_vs_z, dux_vs_z, duy_vs_z, duz_vs_z, dr_vs_z, dur_vs_z, dphi_vs_z, dtheta_vs_z]
    h_dp_vs_z.append(dp_vs_z)
    # track resolution vs track ux (MC)
    dx_vs_ux = ROOT.TH2F("h"+tc+"dx_vs_ux","h"+tc+"dx_vs_ux",10,-1,1,100,-0.5,0.5)
    dy_vs_ux = ROOT.TH2F("h"+tc+"dy_vs_ux","h"+tc+"dy_vs_ux",10,-1,1,100,-0.5,0.5)
    dz_vs_ux = ROOT.TH2F("h"+tc+"dz_vs_ux","h"+tc+"dz_vs_ux",10,-1,1,100,-0.5,0.5)
    dux_vs_ux = ROOT.TH2F("h"+tc+"dux_vs_ux","h"+tc+"dux_vs_ux",10,-1,1,100,-0.1,0.1)
    duy_vs_ux = ROOT.TH2F("h"+tc+"duy_vs_ux","h"+tc+"duy_vs_ux",10,-1,1,100,-0.1,0.1)
    duz_vs_ux = ROOT.TH2F("h"+tc+"duz_vs_ux","h"+tc+"duz_vs_ux",10,-1,1,100,-0.1,0.1)
    dr_vs_ux = ROOT.TH2F("h"+tc+"dr_vs_ux","h"+tc+"dr_vs_ux",10,-1,1,100,-0.5,0.5)
    dur_vs_ux = ROOT.TH2F("h"+tc+"dur_vs_ux","h"+tc+"dur_vs_ux",10,-1,1,100,-0.1,0.1)
    dphi_vs_ux = ROOT.TH2F("h"+tc+"dphi_vs_ux","h"+tc+"dphi_vs_ux",10,-1,1,100,-0.3,0.3)
    dtheta_vs_ux = ROOT.TH2F("h"+tc+"dtheta_vs_ux","h"+tc+"dtheta_vs_ux",10,-1,1,100,-0.2,0.2)
    dp_vs_ux = [dx_vs_ux, dy_vs_ux, dz_vs_ux, dux_vs_ux, duy_vs_ux, duz_vs_ux, dr_vs_ux, dur_vs_ux, dphi_vs_ux, dtheta_vs_ux]
    h_dp_vs_ux.append(dp_vs_ux)
    # track resolution vs track uy (MC)
    dx_vs_uy = ROOT.TH2F("h"+tc+"dx_vs_uy","h"+tc+"dx_vs_uy",10,-1,1,100,-0.5,0.5)
    dy_vs_uy = ROOT.TH2F("h"+tc+"dy_vs_uy","h"+tc+"dy_vs_uy",10,-1,1,100,-0.5,0.5)
    dz_vs_uy = ROOT.TH2F("h"+tc+"dz_vs_uy","h"+tc+"dz_vs_uy",10,-1,1,100,-0.5,0.5)
    dux_vs_uy = ROOT.TH2F("h"+tc+"dux_vs_uy","h"+tc+"dux_vs_uy",10,-1,1,100,-0.1,0.1)
    duy_vs_uy = ROOT.TH2F("h"+tc+"duy_vs_uy","h"+tc+"duy_vs_uy",10,-1,1,100,-0.1,0.1)
    duz_vs_uy = ROOT.TH2F("h"+tc+"duz_vs_uy","h"+tc+"duz_vs_uy",10,-1,1,100,-0.1,0.1)
    dr_vs_uy = ROOT.TH2F("h"+tc+"dr_vs_uy","h"+tc+"dr_vs_uy",10,-1,1,100,-0.5,0.5)
    dur_vs_uy = ROOT.TH2F("h"+tc+"dur_vs_uy","h"+tc+"dur_vs_uy",10,-1,1,100,-0.1,0.1)
    dphi_vs_uy = ROOT.TH2F("h"+tc+"dphi_vs_uy","h"+tc+"dphi_vs_uy",10,-1,1,100,-0.3,0.3)
    dtheta_vs_uy = ROOT.TH2F("h"+tc+"dtheta_vs_uy","h"+tc+"dtheta_vs_uy",10,-1,1,100,-0.2,0.2)
    dp_vs_uy = [dx_vs_uy, dy_vs_uy, dz_vs_uy, dux_vs_uy, duy_vs_uy, duz_vs_uy, dr_vs_uy, dur_vs_uy, dphi_vs_uy, dtheta_vs_uy]
    h_dp_vs_uy.append(dp_vs_uy)
    # track resolution vs track uz (MC)
    dx_vs_uz = ROOT.TH2F("h"+tc+"dx_vs_uz","h"+tc+"dx_vs_uz",10,-1,1,100,-0.5,0.5)
    dy_vs_uz = ROOT.TH2F("h"+tc+"dy_vs_uz","h"+tc+"dy_vs_uz",10,-1,1,100,-0.5,0.5)
    dz_vs_uz = ROOT.TH2F("h"+tc+"dz_vs_uz","h"+tc+"dz_vs_uz",10,-1,1,100,-0.5,0.5)
    dux_vs_uz = ROOT.TH2F("h"+tc+"dux_vs_uz","h"+tc+"dux_vs_uz",10,-1,1,100,-0.1,0.1)
    duy_vs_uz = ROOT.TH2F("h"+tc+"duy_vs_uz","h"+tc+"duy_vs_uz",10,-1,1,100,-0.1,0.1)
    duz_vs_uz = ROOT.TH2F("h"+tc+"duz_vs_uz","h"+tc+"duz_vs_uz",10,-1,1,100,-0.1,0.1)
    dr_vs_uz = ROOT.TH2F("h"+tc+"dr_vs_uz","h"+tc+"dr_vs_uz",10,-1,1,100,-0.5,0.5)
    dur_vs_uz = ROOT.TH2F("h"+tc+"dur_vs_uz","h"+tc+"dur_vs_uz",10,-1,1,100,-0.1,0.1)
    dphi_vs_uz = ROOT.TH2F("h"+tc+"dphi_vs_uz","h"+tc+"dphi_vs_uz",10,-1,1,100,-0.3,0.3)
    dtheta_vs_uz = ROOT.TH2F("h"+tc+"dtheta_vs_uz","h"+tc+"dtheta_vs_uz",10,-1,1,100,-0.2,0.2)
    dp_vs_uz = [dx_vs_uz, dy_vs_uz, dz_vs_uz, dux_vs_uz, duy_vs_uz, duz_vs_uz, dr_vs_uz, dur_vs_uz, dphi_vs_uz, dtheta_vs_uz]
    h_dp_vs_uz.append(dp_vs_uz)
    # track resolution vs track momentum (MC)
    dx_vs_mcmom = ROOT.TH2F("h"+tc+"dx_vs_mcmom","h"+tc+"dx_vs_mcmom",10,0,5,100,-0.5,0.5)
    dy_vs_mcmom = ROOT.TH2F("h"+tc+"dy_vs_mcmom","h"+tc+"dy_vs_mcmom",10,0,5,100,-0.5,0.5)
    dz_vs_mcmom = ROOT.TH2F("h"+tc+"dz_vs_mcmom","h"+tc+"dz_vs_mcmom",10,0,5,100,-0.5,0.5)
    dux_vs_mcmom = ROOT.TH2F("h"+tc+"dux_vs_mcmom","h"+tc+"dux_vs_mcmom",10,0,5,100,-0.1,0.1)
    duy_vs_mcmom = ROOT.TH2F("h"+tc+"duy_vs_mcmom","h"+tc+"duy_vs_mcmom",10,0,5,100,-0.1,0.1)
    duz_vs_mcmom = ROOT.TH2F("h"+tc+"duz_vs_mcmom","h"+tc+"duz_vs_mcmom",10,0,5,100,-0.1,0.1)
    dr_vs_mcmom = ROOT.TH2F("h"+tc+"dr_vs_mcmom","h"+tc+"dr_vs_mcmom",10,0,5,100,-0.5,0.5)
    dur_vs_mcmom = ROOT.TH2F("h"+tc+"dur_vs_mcmom","h"+tc+"dur_vs_mcmom",10,0,5,100,-0.1,0.1)
    dphi_vs_mcmom = ROOT.TH2F("h"+tc+"dphi_vs_mcmom","h"+tc+"dphi_vs_mcmom",10,0,5,100,-0.3,0.3)
    dtheta_vs_mcmom = ROOT.TH2F("h"+tc+"dtheta_vs_mcmom","h"+tc+"dtheta_vs_mcmom",10,0,5,100,-0.2,0.2)
    dp_vs_mcmom = [dx_vs_mcmom, dy_vs_mcmom, dz_vs_mcmom, dux_vs_mcmom, duy_vs_mcmom, duz_vs_mcmom, dr_vs_mcmom, dur_vs_mcmom, dphi_vs_mcmom, dtheta_vs_mcmom]
    h_dp_vs_mcmom.append(dp_vs_mcmom)
    # track resolution vs track angle in zx plane
    dx_vs_coszx = ROOT.TH2F("h"+tc+"dx_vs_coszx","h"+tc+"dx_vs_coszx",8,-1,1,100,-0.5,0.5)
    dy_vs_coszx = ROOT.TH2F("h"+tc+"dy_vs_coszx","h"+tc+"dy_vs_coszx",8,-1,1,100,-0.5,0.5)
    dz_vs_coszx = ROOT.TH2F("h"+tc+"dz_vs_coszx","h"+tc+"dz_vs_coszx",8,-1,1,100,-0.5,0.5)
    dux_vs_coszx = ROOT.TH2F("h"+tc+"dux_vs_coszx","h"+tc+"dux_vs_coszx",8,-1,1,100,-0.1,0.1)
    duy_vs_coszx = ROOT.TH2F("h"+tc+"duy_vs_coszx","h"+tc+"duy_vs_coszx",8,-1,1,100,-0.1,0.1)
    duz_vs_coszx = ROOT.TH2F("h"+tc+"duz_vs_coszx","h"+tc+"duz_vs_coszx",8,-1,1,100,-0.1,0.1)
    dr_vs_coszx = ROOT.TH2F("h"+tc+"dr_vs_coszx","h"+tc+"dr_vs_coszx",8,-1,1,100,-0.5,0.5)
    dur_vs_coszx = ROOT.TH2F("h"+tc+"dur_vs_coszx","h"+tc+"dur_vs_coszx",8,-1,1,100,-0.1,0.1)
    dphi_vs_coszx = ROOT.TH2F("h"+tc+"dphi_vs_coszx","h"+tc+"dphi_vs_coszx",8,-1,1,100,-0.3,0.3)
    dtheta_vs_coszx = ROOT.TH2F("h"+tc+"dtheta_vs_coszx","h"+tc+"dtheta_vs_coszx",8,-1,1,100,-0.2,0.2)
    dp_vs_coszx = [dx_vs_coszx, dy_vs_coszx, dz_vs_coszx, dux_vs_coszx, duy_vs_coszx, duz_vs_coszx, dr_vs_coszx, dur_vs_coszx, dphi_vs_coszx, dtheta_vs_coszx]
    h_dp_vs_coszx.append(dp_vs_coszx)
    # track resolution vs track angle in zy plane
    dx_vs_coszy = ROOT.TH2F("h"+tc+"dx_vs_coszy","h"+tc+"dx_vs_coszy",8,-1,1,100,-0.5,0.5)
    dy_vs_coszy = ROOT.TH2F("h"+tc+"dy_vs_coszy","h"+tc+"dy_vs_coszy",8,-1,1,100,-0.5,0.5)
    dz_vs_coszy = ROOT.TH2F("h"+tc+"dz_vs_coszy","h"+tc+"dz_vs_coszy",8,-1,1,100,-0.5,0.5)
    dux_vs_coszy = ROOT.TH2F("h"+tc+"dux_vs_coszy","h"+tc+"dux_vs_coszy",8,-1,1,100,-0.1,0.1)
    duy_vs_coszy = ROOT.TH2F("h"+tc+"duy_vs_coszy","h"+tc+"duy_vs_coszy",8,-1,1,100,-0.1,0.1)
    duz_vs_coszy = ROOT.TH2F("h"+tc+"duz_vs_coszy","h"+tc+"duz_vs_coszy",8,-1,1,100,-0.1,0.1)
    dr_vs_coszy = ROOT.TH2F("h"+tc+"dr_vs_coszy","h"+tc+"dr_vs_coszy",8,-1,1,100,-0.5,0.5)
    dur_vs_coszy = ROOT.TH2F("h"+tc+"dur_vs_coszy","h"+tc+"dur_vs_coszy",8,-1,1,100,-0.1,0.1)
    dphi_vs_coszy = ROOT.TH2F("h"+tc+"dphi_vs_coszy","h"+tc+"dphi_vs_coszy",8,-1,1,100,-0.3,0.3)
    dtheta_vs_coszy = ROOT.TH2F("h"+tc+"dtheta_vs_coszy","h"+tc+"dtheta_vs_coszy",8,-1,1,100,-0.2,0.2)
    dp_vs_coszy = [dx_vs_coszy, dy_vs_coszy, dz_vs_coszy, dux_vs_coszy, duy_vs_coszy, duz_vs_coszy, dr_vs_coszy, dur_vs_coszy, dphi_vs_coszy, dtheta_vs_coszy]
    h_dp_vs_coszy.append(dp_vs_coszy)
    # track resolution vs ismu
    dx_vs_ismu = ROOT.TH2F("h"+tc+"dx_vs_ismu","h"+tc+"dx_vs_ismu",2,0,2,100,-0.5,0.5)
    dy_vs_ismu = ROOT.TH2F("h"+tc+"dy_vs_ismu","h"+tc+"dy_vs_ismu",2,0,2,100,-0.5,0.5)
    dz_vs_ismu = ROOT.TH2F("h"+tc+"dz_vs_ismu","h"+tc+"dz_vs_ismu",2,0,2,100,-0.5,0.5)
    dux_vs_ismu = ROOT.TH2F("h"+tc+"dux_vs_ismu","h"+tc+"dux_vs_ismu",2,0,2,100,-0.1,0.1)
    duy_vs_ismu = ROOT.TH2F("h"+tc+"duy_vs_ismu","h"+tc+"duy_vs_ismu",2,0,2,100,-0.1,0.1)
    duz_vs_ismu = ROOT.TH2F("h"+tc+"duz_vs_ismu","h"+tc+"duz_vs_ismu",2,0,2,100,-0.1,0.1)
    dr_vs_ismu = ROOT.TH2F("h"+tc+"dr_vs_ismu","h"+tc+"dr_vs_ismu",2,0,2,100,-0.5,0.5)
    dur_vs_ismu = ROOT.TH2F("h"+tc+"dur_vs_ismu","h"+tc+"dur_vs_ismu",2,0,2,100,-0.1,0.1)
    dphi_vs_ismu = ROOT.TH2F("h"+tc+"dphi_vs_ismu","h"+tc+"dphi_vs_ismu",10,0,2,100,-0.3,0.3)
    dtheta_vs_ismu = ROOT.TH2F("h"+tc+"dtheta_vs_ismu","h"+tc+"dtheta_vs_ismu",10,0,2,100,-0.2,0.2)
    dp_vs_ismu = [dx_vs_ismu, dy_vs_ismu, dz_vs_ismu, dux_vs_ismu, duy_vs_ismu, duz_vs_ismu, dr_vs_ismu, dur_vs_ismu, dphi_vs_ismu, dtheta_vs_ismu]
    h_dp_vs_ismu.append(dp_vs_ismu)
    # track resolution vs ntkvtx
    dx_vs_ntkvtx = ROOT.TH2F("h"+tc+"dx_vs_ntkvtx","h"+tc+"dx_vs_ntkvtx",5,1,6,100,-0.5,0.5)
    dy_vs_ntkvtx = ROOT.TH2F("h"+tc+"dy_vs_ntkvtx","h"+tc+"dy_vs_ntkvtx",5,1,6,100,-0.5,0.5)
    dz_vs_ntkvtx = ROOT.TH2F("h"+tc+"dz_vs_ntkvtx","h"+tc+"dz_vs_ntkvtx",5,1,6,100,-0.5,0.5)
    dux_vs_ntkvtx = ROOT.TH2F("h"+tc+"dux_vs_ntkvtx","h"+tc+"dux_vs_ntkvtx",5,1,6,100,-0.1,0.1)
    duy_vs_ntkvtx = ROOT.TH2F("h"+tc+"duy_vs_ntkvtx","h"+tc+"duy_vs_ntkvtx",5,1,6,100,-0.1,0.1)
    duz_vs_ntkvtx = ROOT.TH2F("h"+tc+"duz_vs_ntkvtx","h"+tc+"duz_vs_ntkvtx",5,1,6,100,-0.1,0.1)
    dr_vs_ntkvtx = ROOT.TH2F("h"+tc+"dr_vs_ntkvtx","h"+tc+"dr_vs_ntkvtx",5,1,6,100,-0.5,0.5)
    dur_vs_ntkvtx = ROOT.TH2F("h"+tc+"dur_vs_ntkvtx","h"+tc+"dur_vs_ntkvtx",5,1,6,100,-0.1,0.1)
    dphi_vs_ntkvtx = ROOT.TH2F("h"+tc+"dphi_vs_ntkvtx","h"+tc+"dphi_vs_ntkvtx",10,1,6,100,-0.3,0.3)
    dtheta_vs_ntkvtx = ROOT.TH2F("h"+tc+"dtheta_vs_ntkvtx","h"+tc+"dtheta_vs_ntkvtx",10,1,6,100,-0.2,0.2)
    dp_vs_ntkvtx = [dx_vs_ntkvtx, dy_vs_ntkvtx, dz_vs_ntkvtx, dux_vs_ntkvtx, duy_vs_ntkvtx, duz_vs_ntkvtx, dr_vs_ntkvtx, dur_vs_ntkvtx, dphi_vs_ntkvtx, dtheta_vs_ntkvtx]
    h_dp_vs_ntkvtx.append(dp_vs_ntkvtx)
    # track resolution vs contain
    dx_vs_contain = ROOT.TH2F("h"+tc+"dx_vs_contain","h"+tc+"dx_vs_contain",2,0,2,100,-0.5,0.5)
    dy_vs_contain = ROOT.TH2F("h"+tc+"dy_vs_contain","h"+tc+"dy_vs_contain",2,0,2,100,-0.5,0.5)
    dz_vs_contain = ROOT.TH2F("h"+tc+"dz_vs_contain","h"+tc+"dz_vs_contain",2,0,2,100,-0.5,0.5)
    dux_vs_contain = ROOT.TH2F("h"+tc+"dux_vs_contain","h"+tc+"dux_vs_contain",2,0,2,100,-0.1,0.1)
    duy_vs_contain = ROOT.TH2F("h"+tc+"duy_vs_contain","h"+tc+"duy_vs_contain",2,0,2,100,-0.1,0.1)
    duz_vs_contain = ROOT.TH2F("h"+tc+"duz_vs_contain","h"+tc+"duz_vs_contain",2,0,2,100,-0.1,0.1)
    dr_vs_contain = ROOT.TH2F("h"+tc+"dr_vs_contain","h"+tc+"dr_vs_contain",2,0,2,100,-0.5,0.5)
    dur_vs_contain = ROOT.TH2F("h"+tc+"dur_vs_contain","h"+tc+"dur_vs_contain",2,0,2,100,-0.1,0.1)
    dphi_vs_contain = ROOT.TH2F("h"+tc+"dphi_vs_contain","h"+tc+"dphi_vs_contain",10,0,2,100,-0.3,0.3)
    dtheta_vs_contain = ROOT.TH2F("h"+tc+"dtheta_vs_contain","h"+tc+"dtheta_vs_contain",10,0,2,100,-0.2,0.2)
    dp_vs_contain = [dx_vs_contain, dy_vs_contain, dz_vs_contain, dux_vs_contain, duy_vs_contain, duz_vs_contain, dr_vs_contain, dur_vs_contain, dphi_vs_contain, dtheta_vs_contain]
    h_dp_vs_contain.append(dp_vs_contain)
    # track resolution vs mommumcs
    dx_vs_mommumcs = ROOT.TH2F("h"+tc+"dx_vs_mommumcs","h"+tc+"dx_vs_mommumcs",10,0,2.5,100,-0.5,0.5)
    dy_vs_mommumcs = ROOT.TH2F("h"+tc+"dy_vs_mommumcs","h"+tc+"dy_vs_mommumcs",10,0,2.5,100,-0.5,0.5)
    dz_vs_mommumcs = ROOT.TH2F("h"+tc+"dz_vs_mommumcs","h"+tc+"dz_vs_mommumcs",10,0,2.5,100,-0.5,0.5)
    dux_vs_mommumcs = ROOT.TH2F("h"+tc+"dux_vs_mommumcs","h"+tc+"dux_vs_mommumcs",10,0,2.5,100,-0.1,0.1)
    duy_vs_mommumcs = ROOT.TH2F("h"+tc+"duy_vs_mommumcs","h"+tc+"duy_vs_mommumcs",10,0,2.5,100,-0.1,0.1)
    duz_vs_mommumcs = ROOT.TH2F("h"+tc+"duz_vs_mommumcs","h"+tc+"duz_vs_mommumcs",10,0,2.5,100,-0.1,0.1)
    dr_vs_mommumcs = ROOT.TH2F("h"+tc+"dr_vs_mommumcs","h"+tc+"dr_vs_mommumcs",10,0,2.5,100,-0.5,0.5)
    dur_vs_mommumcs = ROOT.TH2F("h"+tc+"dur_vs_mommumcs","h"+tc+"dur_vs_mommumcs",10,0,2.5,100,-0.1,0.1)
    dphi_vs_mommumcs = ROOT.TH2F("h"+tc+"dphi_vs_mommumcs","h"+tc+"dphi_vs_mommumcs",10,0,2.5,100,-0.3,0.3)
    dtheta_vs_mommumcs = ROOT.TH2F("h"+tc+"dtheta_vs_mommumcs","h"+tc+"dtheta_vs_mommumcs",10,0,2.5,100,-0.2,0.2)
    dp_vs_mommumcs = [dx_vs_mommumcs, dy_vs_mommumcs, dz_vs_mommumcs, dux_vs_mommumcs, duy_vs_mommumcs, duz_vs_mommumcs, dr_vs_mommumcs, dur_vs_mommumcs, dphi_vs_mommumcs, dtheta_vs_mommumcs]
    h_dp_vs_mommumcs.append(dp_vs_mommumcs)
    # track resolution vs mommurng
    dx_vs_mommurng = ROOT.TH2F("h"+tc+"dx_vs_mommurng","h"+tc+"dx_vs_mommurng",10,0,2.5,100,-0.5,0.5)
    dy_vs_mommurng = ROOT.TH2F("h"+tc+"dy_vs_mommurng","h"+tc+"dy_vs_mommurng",10,0,2.5,100,-0.5,0.5)
    dz_vs_mommurng = ROOT.TH2F("h"+tc+"dz_vs_mommurng","h"+tc+"dz_vs_mommurng",10,0,2.5,100,-0.5,0.5)
    dux_vs_mommurng = ROOT.TH2F("h"+tc+"dux_vs_mommurng","h"+tc+"dux_vs_mommurng",10,0,2.5,100,-0.1,0.1)
    duy_vs_mommurng = ROOT.TH2F("h"+tc+"duy_vs_mommurng","h"+tc+"duy_vs_mommurng",10,0,2.5,100,-0.1,0.1)
    duz_vs_mommurng = ROOT.TH2F("h"+tc+"duz_vs_mommurng","h"+tc+"duz_vs_mommurng",10,0,2.5,100,-0.1,0.1)
    dr_vs_mommurng = ROOT.TH2F("h"+tc+"dr_vs_mommurng","h"+tc+"dr_vs_mommurng",10,0,2.5,100,-0.5,0.5)
    dur_vs_mommurng = ROOT.TH2F("h"+tc+"dur_vs_mommurng","h"+tc+"dur_vs_mommurng",10,0,2.5,100,-0.1,0.1)
    dphi_vs_mommurng = ROOT.TH2F("h"+tc+"dphi_vs_mommurng","h"+tc+"dphi_vs_mommurng",10,0,2.5,100,-0.3,0.3)
    dtheta_vs_mommurng = ROOT.TH2F("h"+tc+"dtheta_vs_mommurng","h"+tc+"dtheta_vs_mommurng",10,0,2.5,100,-0.2,0.2)
    dp_vs_mommurng = [dx_vs_mommurng, dy_vs_mommurng, dz_vs_mommurng, dux_vs_mommurng, duy_vs_mommurng, duz_vs_mommurng, dr_vs_mommurng, dur_vs_mommurng, dphi_vs_mommurng, dtheta_vs_mommurng]
    h_dp_vs_mommurng.append(dp_vs_mommurng)
    # track resolution vs mompmcs
    dx_vs_mompmcs = ROOT.TH2F("h"+tc+"dx_vs_mompmcs","h"+tc+"dx_vs_mompmcs",10,0,2.5,100,-0.5,0.5)
    dy_vs_mompmcs = ROOT.TH2F("h"+tc+"dy_vs_mompmcs","h"+tc+"dy_vs_mompmcs",10,0,2.5,100,-0.5,0.5)
    dz_vs_mompmcs = ROOT.TH2F("h"+tc+"dz_vs_mompmcs","h"+tc+"dz_vs_mompmcs",10,0,2.5,100,-0.5,0.5)
    dux_vs_mompmcs = ROOT.TH2F("h"+tc+"dux_vs_mompmcs","h"+tc+"dux_vs_mompmcs",10,0,2.5,100,-0.1,0.1)
    duy_vs_mompmcs = ROOT.TH2F("h"+tc+"duy_vs_mompmcs","h"+tc+"duy_vs_mompmcs",10,0,2.5,100,-0.1,0.1)
    duz_vs_mompmcs = ROOT.TH2F("h"+tc+"duz_vs_mompmcs","h"+tc+"duz_vs_mompmcs",10,0,2.5,100,-0.1,0.1)
    dr_vs_mompmcs = ROOT.TH2F("h"+tc+"dr_vs_mompmcs","h"+tc+"dr_vs_mompmcs",10,0,2.5,100,-0.5,0.5)
    dur_vs_mompmcs = ROOT.TH2F("h"+tc+"dur_vs_mompmcs","h"+tc+"dur_vs_mompmcs",10,0,2.5,100,-0.1,0.1)
    dphi_vs_mompmcs = ROOT.TH2F("h"+tc+"dphi_vs_mompmcs","h"+tc+"dphi_vs_mompmcs",10,0,2.5,100,-0.3,0.3)
    dtheta_vs_mompmcs = ROOT.TH2F("h"+tc+"dtheta_vs_mompmcs","h"+tc+"dtheta_vs_mompmcs",10,0,2.5,100,-0.2,0.2)
    dp_vs_mompmcs = [dx_vs_mompmcs, dy_vs_mompmcs, dz_vs_mompmcs, dux_vs_mompmcs, duy_vs_mompmcs, duz_vs_mompmcs, dr_vs_mompmcs, dur_vs_mompmcs, dphi_vs_mompmcs, dtheta_vs_mompmcs]
    h_dp_vs_mompmcs.append(dp_vs_mompmcs)
    # track resolution vs momprng
    dx_vs_momprng = ROOT.TH2F("h"+tc+"dx_vs_momprng","h"+tc+"dx_vs_momprng",10,0,2.5,100,-0.5,0.5)
    dy_vs_momprng = ROOT.TH2F("h"+tc+"dy_vs_momprng","h"+tc+"dy_vs_momprng",10,0,2.5,100,-0.5,0.5)
    dz_vs_momprng = ROOT.TH2F("h"+tc+"dz_vs_momprng","h"+tc+"dz_vs_momprng",10,0,2.5,100,-0.5,0.5)
    dux_vs_momprng = ROOT.TH2F("h"+tc+"dux_vs_momprng","h"+tc+"dux_vs_momprng",10,0,2.5,100,-0.1,0.1)
    duy_vs_momprng = ROOT.TH2F("h"+tc+"duy_vs_momprng","h"+tc+"duy_vs_momprng",10,0,2.5,100,-0.1,0.1)
    duz_vs_momprng = ROOT.TH2F("h"+tc+"duz_vs_momprng","h"+tc+"duz_vs_momprng",10,0,2.5,100,-0.1,0.1)
    dr_vs_momprng = ROOT.TH2F("h"+tc+"dr_vs_momprng","h"+tc+"dr_vs_momprng",10,0,2.5,100,-0.5,0.5)
    dur_vs_momprng = ROOT.TH2F("h"+tc+"dur_vs_momprng","h"+tc+"dur_vs_momprng",10,0,2.5,100,-0.1,0.1)
    dphi_vs_momprng = ROOT.TH2F("h"+tc+"dphi_vs_momprng","h"+tc+"dphi_vs_momprng",10,0,2.5,100,-0.3,0.3)
    dtheta_vs_momprng = ROOT.TH2F("h"+tc+"dtheta_vs_momprng","h"+tc+"dtheta_vs_momprng",10,0,2.5,100,-0.2,0.2)
    dp_vs_momprng = [dx_vs_momprng, dy_vs_momprng, dz_vs_momprng, dux_vs_momprng, duy_vs_momprng, duz_vs_momprng, dr_vs_momprng, dur_vs_momprng, dphi_vs_momprng, dtheta_vs_momprng]
    h_dp_vs_momprng.append(dp_vs_momprng)
    
for t in tree: 
    passSelII = tree.passSelII
    tk_nhits = tree.tk_nhits
    tk1_nhits = tree.tk1_nhits
    tk2_nhits = tree.tk2_nhits
    tk_mommurng = tree.tk_mommurng
    tk_mommumcs = tree.tk_mommumcs
    tk_mommumcserr = tree.tk_mommumcserr
    tk_momprng = tree.tk_momprng
    tk_mompmcs = tree.tk_mompmcs
    tk_mompmcserr = tree.tk_mompmcserr
    tk_id = tree.tk_id
    tk1_id = tree.tk1_id
    tk2_id = tree.tk2_id
    mc_mom = tree.mc_mom
    mc_pid = tree.mc_pid
    vx_ntks = tree.vx_ntks
    tk_length = tree.tk_length
    tk1_length = tree.tk1_length
    tk2_length = tree.tk2_length
    tk_contain = tree.tk_contain
    tk_drvtx = tree.tk_drvtx
    #
    # Make sure all vectors are filled
    #
    if len(tk_gpar)<6: continue
    if len(tk1_gpar_vtx)<6: continue
    if len(tk1_gpar_end)<6: continue
    if len(tk2_gpar_vtx)<6: continue
    if len(tk2_gpar_end)<6: continue
    #
    if len(tk2_tk1_gpar_vtx)<6: continue
    #
    if len(tk1_tk_gpar_mid)<6: continue
    if len(tk2_tk_gpar_mid)<6: continue
    #
    if not isdata and len(tk1_mc_gpar_mid)<6: continue
    if not isdata and len(tk2_mc_gpar_mid)<6: continue
    if not isdata and len(tk_mc_gpar_vtx)<6: continue
    if not isdata and len(tk_mc_gpar_end)<6: continue
    #
    # Cuts
    #
    # minimum number of hits
    if tk_nhits<30 or tk1_nhits<10 or tk2_nhits<10: continue
    # sanity cuts on lengths
    if tk_length<(tk1_length+tk2_length): continue
    if tk_length>(50+tk1_length+tk2_length): continue
    #
    if tk_contain == 1: continue
    if tk_gpar[2]<100: continue
    if tk_gpar[4]<-0.4 or tk_gpar[4]>0.4: continue
    #if tk_gpar[5]<0.: continue
    if vx_ntks<3: continue
    #
    # temporary cut to test mc with and without cosmics!
    #if mc_mom<0: continue
    #if mc_mom>0: continue
    #if tk_length<100: continue
    #if abs(mc_pid)!=13: continue
    #if abs(mc_pid)!=2212: continue
    #
    #if tk_contain == 0: continue
    #if not passSelII: continue
    #if vx_ntks<2 or vx_ntks>4: continue
    #if tk_gpar[0]>175: continue
    #if tk_gpar[1]<-75 or tk_gpar[1]>75: continue
    #if tk_gpar[4]<-0.8 or tk_gpar[4]>0.8: continue
    #
    #if tk_gpar[4]<-0.6 or tk_gpar[4]>0.6: continue
    #if vx_ntks<3: continue
    # test hi mom
    #if tk_mommumcs<1.0: continue
    #
    #if not passSelII: continue
    #tmpif vx_ntks<2 or vx_ntks>4: continue
    #tmpif tk_gpar[1]<-75 or tk_gpar[1]>75: continue
    #tmpif tk_gpar[0]>175: continue
    #tmpif tk_gpar[5]<0.2: continue
    #if tk_drvtx>5.: continue
    #
    #if tk_mommumcs<1.0: continue
    #if vtxpos[0]<50 or vtxpos[0]>200: continue
    #if vtxpos[1]<-70 or vtxpos[1]>70: continue
    #if vtxpos[2]<100 or vtxpos[2]>1000: continue
    #if tk_mommumcs>1.0: continue
    #if tk_gpar[4]<-0.4 or tk_gpar[4]>0.4: continue
    #if tk_gpar[4]<-0.6 or tk_gpar[4]>0.6: continue
    #if tk_gpar[4]>-0.6 and tk_gpar[4]<0.6: continue
    #if tk_mommumcs>1.0: 
    #    if tk_gpar[4]<-0.4 or tk_gpar[4]>0.4: continue
    #
    # make sure the two split tracks have the same direction (minor effect)
    dot_12 = tk1_gpar_vtx[3]*tk2_gpar_vtx[3]+tk1_gpar_vtx[4]*tk2_gpar_vtx[4]+tk1_gpar_vtx[5]*tk2_gpar_vtx[5]
    dot_12 = dot_12/abs(dot_12)
    if dot_12<0: continue
    #
    # make sure the split tracks have the same direction as the original track
    #dot_t1 = tk_gpar[3]*tk1_gpar_vtx[3]+tk_gpar[4]*tk1_gpar_vtx[4]+tk_gpar[5]*tk1_gpar_vtx[5]
    #dot_t1 = dot_t1/abs(dot_t1)
    #if dot_t1<0: continue
    #
    dot_mc1 = 1
    dot_mc2 = 1
    dot_mct = 1
    if not isdata: 
        dot_mc1 = mc_gpar_vtx[3]*tk1_gpar_vtx[3]+mc_gpar_vtx[4]*tk1_gpar_vtx[4]+mc_gpar_vtx[5]*tk1_gpar_vtx[5]
        dot_mc2 = mc_gpar_vtx[3]*tk2_gpar_vtx[3]+mc_gpar_vtx[4]*tk2_gpar_vtx[4]+mc_gpar_vtx[5]*tk2_gpar_vtx[5]
        dot_mct = mc_gpar_vtx[3]*tk_gpar[3]+mc_gpar_vtx[4]*tk_gpar[4]+mc_gpar_vtx[5]*tk_gpar[5]
    #
    # now, let's figure out what are the parameters to use
    p1 = []
    p2 = []
    d1 = []
    d2 = []
    p1m = []
    pm1 = []
    d1m = []
    dm1 = []
    #
    # get parameters at midpoint for splitting at half
    #
    p1 = [tk1_tk_gpar_mid[0], tk1_tk_gpar_mid[1], tk1_tk_gpar_mid[2]]
    p2 = [tk2_tk_gpar_mid[0], tk2_tk_gpar_mid[1], tk2_tk_gpar_mid[2]]
    d1 = [tk1_tk_gpar_mid[3], tk1_tk_gpar_mid[4], tk1_tk_gpar_mid[5]]
    d2 = [tk2_tk_gpar_mid[3], tk2_tk_gpar_mid[4], tk2_tk_gpar_mid[5]]        

    if not isdata:
        p1m = [tk1_mc_gpar_mid[0], tk1_mc_gpar_mid[1], tk1_mc_gpar_mid[2]]
        pm1 = [mc_tk_gpar_mid[0], mc_tk_gpar_mid[1], mc_tk_gpar_mid[2]]
        dm1 = [mc_tk_gpar_mid[3], mc_tk_gpar_mid[4], mc_tk_gpar_mid[5]]
        if dot_mc1>0:
            d1m = [tk1_mc_gpar_mid[3], tk1_mc_gpar_mid[4], tk1_mc_gpar_mid[5]] # *** no change for end of track 1 !!!!!
        else:
            d1m = [-tk1_mc_gpar_mid[3], -tk1_mc_gpar_mid[4], -tk1_mc_gpar_mid[5]] # *** reverse for vtx of track 1 !!!!!

        # make sure tk1_mc_gpar_mid is not in the middle of track1, nor track2 (i.e. projection to a point in between the two tracks, just like the split case)
        midvtx1 = [tk1_mc_gpar_mid[0]-tk1_gpar_vtx[0], tk1_mc_gpar_mid[1]-tk1_gpar_vtx[1], tk1_mc_gpar_mid[2]-tk1_gpar_vtx[2]]
        midend1 = [tk1_mc_gpar_mid[0]-tk1_gpar_end[0], tk1_mc_gpar_mid[1]-tk1_gpar_end[1], tk1_mc_gpar_mid[2]-tk1_gpar_end[2]]
        fotmidvtxend1 = midvtx1[0]*midend1[0] + midvtx1[1]*midend1[1] + midvtx1[2]*midend1[2]
        midvtx2 = [tk1_mc_gpar_mid[0]-tk2_gpar_vtx[0], tk1_mc_gpar_mid[1]-tk1_gpar_vtx[1], tk1_mc_gpar_mid[2]-tk2_gpar_vtx[2]]
        midend2 = [tk1_mc_gpar_mid[0]-tk2_gpar_end[0], tk1_mc_gpar_mid[1]-tk1_gpar_end[1], tk1_mc_gpar_mid[2]-tk2_gpar_end[2]]
        fotmidvtxend2 = midvtx2[0]*midend2[0] + midvtx2[1]*midend2[1] + midvtx2[2]*midend2[2]
        if fotmidvtxend1<=0 or fotmidvtxend2<=0: continue
   
    #print 'fotmidvtxend1', fotmidvtxend1
    #print 'fotmidvtxend2', fotmidvtxend2
    #print 'tk1_gpar_vtx: ', tk1_gpar_vtx[0], tk1_gpar_vtx[1], tk1_gpar_vtx[2], tk1_gpar_vtx[3], tk1_gpar_vtx[4], tk1_gpar_vtx[5]
    #print 'tk1_gpar_end: ', tk1_gpar_end[0], tk1_gpar_end[1], tk1_gpar_end[2], tk1_gpar_end[3], tk1_gpar_end[4], tk1_gpar_end[5]

    #print 'tk2_gpar_vtx: ', tk2_gpar_vtx[0], tk2_gpar_vtx[1], tk2_gpar_vtx[2], tk2_gpar_vtx[3], tk2_gpar_vtx[4], tk2_gpar_vtx[5]
    #print 'tk2_gpar_end: ', tk2_gpar_end[0], tk2_gpar_end[1], tk2_gpar_end[2], tk2_gpar_end[3], tk2_gpar_end[4], tk2_gpar_end[5]
    #print ''
    ##print 'tk2_gpar_vtx: ', tk2_gpar_vtx[0], tk2_gpar_vtx[1], tk2_gpar_vtx[2], tk2_gpar_vtx[3], tk2_gpar_vtx[4], tk2_gpar_vtx[5]
    ##print 'tk2_gpar_end: ', tk2_gpar_end[0], tk2_gpar_end[1], tk2_gpar_end[2], tk2_gpar_end[3], tk2_gpar_end[4], tk2_gpar_end[5]

    ptm = []
    pmt = []
    dtm = []
    dmt = []
    if not isdata:
        if dot_mct>0:
            ptm = [tk_mc_gpar_vtx[0], tk_mc_gpar_vtx[1], tk_mc_gpar_vtx[2]]
            pmt = [mc_tk_gpar_vtx[0], mc_tk_gpar_vtx[1], mc_tk_gpar_vtx[2]]
            dtm = [tk_mc_gpar_vtx[3], tk_mc_gpar_vtx[4], tk_mc_gpar_vtx[5]]
            dmt = [mc_tk_gpar_vtx[3], mc_tk_gpar_vtx[4], mc_tk_gpar_vtx[5]]
        else:
            ptm = [tk_mc_gpar_end[0], tk_mc_gpar_end[1], tk_mc_gpar_end[2]]
            pmt = [mc_tk_gpar_end[0], mc_tk_gpar_end[1], mc_tk_gpar_end[2]]
            dtm = [-tk_mc_gpar_end[3], -tk_mc_gpar_end[4], -tk_mc_gpar_end[5]]
            dmt = [mc_tk_gpar_end[3], mc_tk_gpar_end[4], mc_tk_gpar_end[5]]

    # normalize momentum to direction
    ap1 = math.sqrt(d1[0]*d1[0]+d1[1]*d1[1]+d1[2]*d1[2])
    ap2 = math.sqrt(d2[0]*d2[0]+d2[1]*d2[1]+d2[2]*d2[2])
    ap1m = 0.
    apm1 = 0.
    aptm = 0.
    apmt = 0.
    if not isdata:
        ap1m = math.sqrt(d1m[0]*d1m[0]+d1m[1]*d1m[1]+d1m[2]*d1m[2])
        apm1 = math.sqrt(dm1[0]*dm1[0]+dm1[1]*dm1[1]+dm1[2]*dm1[2])
        aptm = math.sqrt(dtm[0]*dtm[0]+dtm[1]*dtm[1]+dtm[2]*dtm[2])
        apmt = math.sqrt(dmt[0]*dmt[0]+dmt[1]*dmt[1]+dmt[2]*dmt[2])
    # sanity check
    checkok = True
    for i in range(0,3):
        if math.fabs(p1[i])<1e-09 or math.isinf(p1[i]) or math.isnan(p1[i]): checkok = False
        if math.fabs(p2[i])<1e-09 or math.isinf(p2[i]) or math.isnan(p2[i]): checkok = False
        if math.fabs(d1[i])<1e-09 or math.isinf(d1[i]) or math.isnan(d1[i]): checkok = False
        if math.fabs(d2[i])<1e-09 or math.isinf(d2[i]) or math.isnan(d2[i]): checkok = False
    if not checkok: continue
    #
    for i in range(0,3):
        d1[i] = d1[i]/ap1
        d2[i] = d2[i]/ap2
        if not isdata:
            d1m[i] = d1m[i]/ap1m
            dm1[i] = dm1[i]/apm1
            dtm[i] = dtm[i]/aptm
            dmt[i] = dmt[i]/apmt
    #
    #print p1
    #print p2
    #print d1
    #print d2
    #
    #print testcase, d1m, dm1
    #print tree.run, tree.subrun, tree.eventid
    for i in range(0,3):
        if math.fabs(d1[i]+d2[i])/2. > 0.9: continue
        h12dp[i].Fill((p1[i]-p2[i])/math.sqrt(2.))
        d_12 = (d1[i]-d2[i])
        h12dp[i+3].Fill(d_12/math.sqrt(2.))
        d_12 = d_12/(0.5*(d1[i]+d2[i]))
        h12dp[i+6].Fill(d_12/math.sqrt(2.))
        #
        dpos = [(p1[i]-p2[i])/math.sqrt(2.)]
        ddir = [(d1[i]-d2[i])/math.sqrt(2.)]
        dr  = [ (math.sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2])-math.sqrt(p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]))/math.sqrt(2.) ]
        dur = [ (math.sqrt(d1[0]*d1[0]+d1[1]*d1[1]+d1[2]*d1[2])-math.sqrt(d2[0]*d2[0]+d2[1]*d2[1]+d2[2]*d2[2]))/math.sqrt(2.) ]
        dp = ( math.atan2(d1[1],d1[0]) - math.atan2(d2[1],d2[0]) )/math.sqrt(2.)
        if dp<(-2*math.pi): dp = dp+2*math.pi
        if dp>(2*math.pi): dp = dp-2*math.pi
        dphi = [ dp ]
        dt = ( math.atan2(math.sqrt(d1[0]*d1[0]+d1[1]*d1[1]),d1[2]) - math.atan2(math.sqrt(d2[0]*d2[0]+d2[1]*d2[1]),d2[2]) )/math.sqrt(2.)
        if dt<(-2*math.pi): dt = dt+2*math.pi
        if dt>(2*math.pi): dt = dt-2*math.pi
        dtheta = [ dt ]
        if not isdata: 
            dpos = [(p1[i]-p2[i])/math.sqrt(2.), p1m[i]-pm1[i], ptm[i]-pmt[i]]
            ddir = [(d1[i]-d2[i])/math.sqrt(2.), d1m[i]-dm1[i], dtm[i]-dmt[i]]
            dr.append( math.sqrt(p1m[0]*p1m[0]+p1m[1]*p1m[1]+p1m[2]*p1m[2])-math.sqrt(pm1[0]*pm1[0]+pm1[1]*pm1[1]+pm1[2]*pm1[2]) )
            dr.append( math.sqrt(ptm[0]*ptm[0]+ptm[1]*ptm[1]+ptm[2]*ptm[2])-math.sqrt(pmt[0]*pmt[0]+pmt[1]*pmt[1]+pmt[2]*pmt[2]) )
            dur.append( math.sqrt(d1m[0]*d1m[0]+d1m[1]*d1m[1]+d1m[2]*d1m[2])-math.sqrt(dm1[0]*dm1[0]+dm1[1]*dm1[1]+dm1[2]*dm1[2]) )
            dur.append( math.sqrt(dtm[0]*dtm[0]+dtm[1]*dtm[1]+dtm[2]*dtm[2])-math.sqrt(dmt[0]*dmt[0]+dmt[1]*dmt[1]+dmt[2]*dmt[2]) )
            dp = math.atan2(d1m[1],d1m[0]) - math.atan2(dm1[1],dm1[0])
            if dp<(-2*math.pi): dp = dp+2*math.pi
            if dp>(2*math.pi): dp = dp-2*math.pi
            dphi.append( dp )
            dp = math.atan2(dtm[1],dtm[0]) - math.atan2(dmt[1],dmt[0])
            if dp<(-2*math.pi): dp = dp+2*math.pi
            if dp>(2*math.pi): dp = dp-2*math.pi
            dphi.append( dp )
            dt = math.atan2(math.sqrt(d1m[0]*d1m[0]+d1m[1]*d1m[1]),d1m[2]) - math.atan2(math.sqrt(dm1[0]*dm1[0]+dm1[1]*dm1[1]),dm1[2])
            if dt<(-2*math.pi): dt = dt+2*math.pi
            if dt>(2*math.pi): dt = dt-2*math.pi
            dtheta.append( dt )
            dt = math.atan2(math.sqrt(dtm[0]*dtm[0]+dtm[1]*dtm[1]),dtm[2]) - math.atan2(math.sqrt(dmt[0]*dmt[0]+dmt[1]*dmt[1]),dmt[2])
            if dt<(-2*math.pi): dt = dt+2*math.pi
            if dt>(2*math.pi): dt = dt-2*math.pi
            dtheta.append( dt )
        for tci in range(0,len(testcase)):
            length = tk1_length
            nhits = tk1_nhits
            #vpos = [tk_gpar[0], tk_gpar[1], tk_gpar[2]] #p1 #FIXME!!!
            vpos = p1
            if testcase[tci] is "tm":
                length = tk_length
                nhits = tk_nhits
                vpos = ptm
            h_dp_vs_len[tci][i].Fill(length,dpos[tci])
            h_dp_vs_len[tci][i+3].Fill(length,ddir[tci])
            #
            h_dp_vs_nh[tci][i].Fill(nhits,dpos[tci])
            h_dp_vs_nh[tci][i+3].Fill(nhits,ddir[tci])
            #
            h_dp_vs_x[tci][i].Fill(vpos[0],dpos[tci])
            h_dp_vs_x[tci][i+3].Fill(vpos[0],ddir[tci])
            #
            h_dp_vs_y[tci][i].Fill(vpos[1],dpos[tci])
            h_dp_vs_y[tci][i+3].Fill(vpos[1],ddir[tci])
            #
            h_dp_vs_z[tci][i].Fill(vpos[2],dpos[tci])
            h_dp_vs_z[tci][i+3].Fill(vpos[2],ddir[tci])
            #
            h_dp_vs_ux[tci][i].Fill(tk_gpar[3],dpos[tci])
            h_dp_vs_ux[tci][i+3].Fill(tk_gpar[3],ddir[tci])
            #
            h_dp_vs_uy[tci][i].Fill(tk_gpar[4],dpos[tci])
            h_dp_vs_uy[tci][i+3].Fill(tk_gpar[4],ddir[tci])
            #
            h_dp_vs_uz[tci][i].Fill(tk_gpar[5],dpos[tci])
            h_dp_vs_uz[tci][i+3].Fill(tk_gpar[5],ddir[tci])
            #
            h_dp_vs_mcmom[tci][i].Fill(mc_mom,dpos[tci])
            h_dp_vs_mcmom[tci][i+3].Fill(mc_mom,ddir[tci])
            #
            coszx = tk_gpar[5]/math.sqrt(tk_gpar[3]*tk_gpar[3]+tk_gpar[5]*tk_gpar[5])
            h_dp_vs_coszx[tci][i].Fill(coszx,dpos[tci])
            h_dp_vs_coszx[tci][i+3].Fill(coszx,ddir[tci])
            #
            coszy = tk_gpar[5]/math.sqrt(tk_gpar[4]*tk_gpar[4]+tk_gpar[5]*tk_gpar[5])
            h_dp_vs_coszy[tci][i].Fill(coszy,dpos[tci])
            h_dp_vs_coszy[tci][i+3].Fill(coszy,ddir[tci])
            #
            if mc_pid!=-999:
                ismu = 1 if abs(mc_pid)==13 else 0
                h_dp_vs_ismu[tci][i].Fill(ismu,dpos[tci])
                h_dp_vs_ismu[tci][i+3].Fill(ismu,ddir[tci])
            #
            ntkvtx = 1 if vx_ntks==-999 else vx_ntks
            h_dp_vs_ntkvtx[tci][i].Fill(ntkvtx,dpos[tci])
            h_dp_vs_ntkvtx[tci][i+3].Fill(ntkvtx,ddir[tci])
            #
            h_dp_vs_contain[tci][i].Fill(tk_contain,dpos[tci])
            h_dp_vs_contain[tci][i+3].Fill(tk_contain,ddir[tci])
            #
            h_dp_vs_mommumcs[tci][i].Fill(tk_mommumcs,dpos[tci])
            h_dp_vs_mommumcs[tci][i+3].Fill(tk_mommumcs,ddir[tci])
            #
            h_dp_vs_mommurng[tci][i].Fill(tk_mommurng,dpos[tci])
            h_dp_vs_mommurng[tci][i+3].Fill(tk_mommurng,ddir[tci])
            #
            h_dp_vs_mompmcs[tci][i].Fill(tk_mompmcs,dpos[tci])
            h_dp_vs_mompmcs[tci][i+3].Fill(tk_mompmcs,ddir[tci])
            #
            h_dp_vs_momprng[tci][i].Fill(tk_momprng,dpos[tci])
            h_dp_vs_momprng[tci][i+3].Fill(tk_momprng,ddir[tci])
            if i==0:
                h_dp_vs_len[tci][6].Fill(length,dr[tci])
                h_dp_vs_len[tci][7].Fill(length,dur[tci])
                h_dp_vs_len[tci][8].Fill(length,dphi[tci])
                h_dp_vs_len[tci][9].Fill(length,dtheta[tci])
                #
                h_dp_vs_nh[tci][6].Fill(nhits,dr[tci])
                h_dp_vs_nh[tci][7].Fill(nhits,dur[tci])
                h_dp_vs_nh[tci][8].Fill(nhits,dphi[tci])
                h_dp_vs_nh[tci][9].Fill(nhits,dtheta[tci])
                #
                h_dp_vs_x[tci][6].Fill(vpos[0],dr[tci])
                h_dp_vs_x[tci][7].Fill(vpos[0],dur[tci])
                h_dp_vs_x[tci][8].Fill(vpos[0],dphi[tci])
                h_dp_vs_x[tci][9].Fill(vpos[0],dtheta[tci])
                #
                h_dp_vs_y[tci][6].Fill(vpos[1],dr[tci])
                h_dp_vs_y[tci][7].Fill(vpos[1],dur[tci])
                h_dp_vs_y[tci][8].Fill(vpos[1],dphi[tci])
                h_dp_vs_y[tci][9].Fill(vpos[1],dtheta[tci])
                #
                h_dp_vs_z[tci][6].Fill(vpos[2],dr[tci])
                h_dp_vs_z[tci][7].Fill(vpos[2],dur[tci])
                h_dp_vs_z[tci][8].Fill(vpos[2],dphi[tci])
                h_dp_vs_z[tci][9].Fill(vpos[2],dtheta[tci])
                #
                h_dp_vs_ux[tci][6].Fill(tk_gpar[3],dr[tci])
                h_dp_vs_ux[tci][7].Fill(tk_gpar[3],dur[tci])
                h_dp_vs_ux[tci][8].Fill(tk_gpar[3],dphi[tci])
                h_dp_vs_ux[tci][9].Fill(tk_gpar[3],dtheta[tci])
                #
                h_dp_vs_uy[tci][6].Fill(tk_gpar[4],dr[tci])
                h_dp_vs_uy[tci][7].Fill(tk_gpar[4],dur[tci])
                h_dp_vs_uy[tci][8].Fill(tk_gpar[4],dphi[tci])
                h_dp_vs_uy[tci][9].Fill(tk_gpar[4],dtheta[tci])
                #
                h_dp_vs_uz[tci][6].Fill(tk_gpar[5],dr[tci])
                h_dp_vs_uz[tci][7].Fill(tk_gpar[5],dur[tci])
                h_dp_vs_uz[tci][8].Fill(tk_gpar[5],dphi[tci])
                h_dp_vs_uz[tci][9].Fill(tk_gpar[5],dtheta[tci])
                #
                h_dp_vs_mcmom[tci][6].Fill(mc_mom,dr[tci])
                h_dp_vs_mcmom[tci][7].Fill(mc_mom,dur[tci])
                h_dp_vs_mcmom[tci][8].Fill(mc_mom,dphi[tci])
                h_dp_vs_mcmom[tci][9].Fill(mc_mom,dtheta[tci])
                #
                h_dp_vs_coszx[tci][6].Fill(coszx,dr[tci])
                h_dp_vs_coszx[tci][7].Fill(coszx,dur[tci])
                h_dp_vs_coszx[tci][8].Fill(coszx,dphi[tci])
                h_dp_vs_coszx[tci][9].Fill(coszx,dtheta[tci])
                #
                h_dp_vs_coszy[tci][6].Fill(coszy,dr[tci])
                h_dp_vs_coszy[tci][7].Fill(coszy,dur[tci])
                h_dp_vs_coszy[tci][8].Fill(coszy,dphi[tci])
                h_dp_vs_coszy[tci][9].Fill(coszy,dtheta[tci])
                #
                if mc_pid!=-999:
                    ismu = 1 if abs(mc_pid)==13 else 0
                    h_dp_vs_ismu[tci][6].Fill(ismu,dr[tci])
                    h_dp_vs_ismu[tci][7].Fill(ismu,dur[tci])
                    h_dp_vs_ismu[tci][8].Fill(ismu,dphi[tci])
                    h_dp_vs_ismu[tci][9].Fill(ismu,dtheta[tci])
                #
                h_dp_vs_ntkvtx[tci][6].Fill(ntkvtx,dr[tci])
                h_dp_vs_ntkvtx[tci][7].Fill(ntkvtx,dur[tci])
                h_dp_vs_ntkvtx[tci][8].Fill(ntkvtx,dphi[tci])
                h_dp_vs_ntkvtx[tci][9].Fill(ntkvtx,dtheta[tci])
                #
                h_dp_vs_contain[tci][6].Fill(tk_contain,dr[tci])
                h_dp_vs_contain[tci][7].Fill(tk_contain,dur[tci])
                h_dp_vs_contain[tci][8].Fill(tk_contain,dphi[tci])
                h_dp_vs_contain[tci][9].Fill(tk_contain,dtheta[tci])
                #
                h_dp_vs_mommumcs[tci][6].Fill(tk_mommumcs,dr[tci])
                h_dp_vs_mommumcs[tci][7].Fill(tk_mommumcs,dur[tci])
                h_dp_vs_mommumcs[tci][8].Fill(tk_mommumcs,dphi[tci])
                h_dp_vs_mommumcs[tci][9].Fill(tk_mommumcs,dtheta[tci])
                #
                h_dp_vs_mommurng[tci][6].Fill(tk_mommurng,dr[tci])
                h_dp_vs_mommurng[tci][7].Fill(tk_mommurng,dur[tci])
                h_dp_vs_mommurng[tci][8].Fill(tk_mommurng,dphi[tci])
                h_dp_vs_mommurng[tci][9].Fill(tk_mommurng,dtheta[tci])
                #
                h_dp_vs_mompmcs[tci][6].Fill(tk_mompmcs,dr[tci])
                h_dp_vs_mompmcs[tci][7].Fill(tk_mompmcs,dur[tci])
                h_dp_vs_mompmcs[tci][8].Fill(tk_mompmcs,dphi[tci])
                h_dp_vs_mompmcs[tci][9].Fill(tk_mompmcs,dtheta[tci])
                #
                h_dp_vs_momprng[tci][6].Fill(tk_momprng,dr[tci])
                h_dp_vs_momprng[tci][7].Fill(tk_momprng,dur[tci])
                h_dp_vs_momprng[tci][8].Fill(tk_momprng,dphi[tci])
                h_dp_vs_momprng[tci][9].Fill(tk_momprng,dtheta[tci])

#for i in range(0,9):
#    xlab = "dx"
#    if i==1: xlab = "dy"
#    if i==2: xlab = "dz"
#    if i==3: xlab = "dux"
#    if i==4: xlab = "duy"
#    if i==5: xlab = "duz"
#    if i==6: xlab = "duxr"
#    if i==7: xlab = "duyr"
#    if i==8: xlab = "duzr"
#    #
#    c = ROOT.TCanvas("c","c",800,800)
#    hdd = h12dp[i]
#    hdd.SetTitle("")
#    hdd.GetXaxis().SetTitle(xlab)
#    hdd.GetYaxis().SetTitle("fraction of tracks")
#    hdd.GetYaxis().SetTitleOffset(1.4)
#    hdd.SetLineColor(2)
#    #hdd.Scale(1./hdd.Integral())
#    hdd.Draw()
#    leg = ROOT.TLegend(0.1,0.9,0.9,0.99)
#    leg.SetNColumns(3)
#    leg.AddEntry(hdd,"split 1,2","L")
#    leg.Draw()
#    c.SaveAs(fnme+"_"+xlab+tag+".png")
    
fout.Write()
fout.Close()
