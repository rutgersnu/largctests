import ROOT as R

f = R.TFile.Open("/uboone/app/users/cerati/valid-mcc9-v0/srcs/tbt.root")

t = f.Get("testbacktrack/tree")

h = R.TH1F("nmcp","nmcp",10,0,10)

for hit in t:
    nmcp = 0
    for imcp in range(0,t.nmcp): 
        if t.mcp_ideFraction[imcp]>0.1: nmcp = nmcp+1
    h.Fill(nmcp)

h.DrawNormalized()
