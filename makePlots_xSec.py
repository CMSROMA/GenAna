import ROOT
import os


ROOT.gROOT.SetBatch(True)

def plotFromDifferentTrees(files):
    # Loop over the files
    gr = {}
    legend = ROOT.TLegend(0.75, 0.75, 0.95, 0.95)
    for i, file_name in enumerate(files):
        # Open the file
        f = ROOT.TFile(file_name)
	base_name = os.path.splitext(os.path.basename(file_name))[0]
	split = base_name.split('_')
	for part in split:
		if 'M' in part:
			mass = float(part.split('M')[1])
		if 'Lambda' in part:
			l = (part.split('Lambda')[1])
			print(l)
       	if not mass:
		mass = 0
		print("Warning file name without mass check inputs")
	if not l:
		l=0
		print("Warning file name without couplings check inputs")
	if l not in gr:
		gr['%s' %l] = ROOT.TGraphErrors()
   	# Create the TGraphErrors
   	#gr = ROOT.TGraphErrors()
       # Get the tree
        tree = f.Get("genAnalyzer/xSec")

        # Loop over the events
        for j, ev in enumerate(tree):
            # Fill the TGraphErrors
	    print(mass)
            gr['%s'%l].SetPoint(i*tree.GetEntries()+j, mass, ev.xSec_)
            gr['%s'%l].SetPointError(i*tree.GetEntries()+j, 0, ev.xSec_err_)
    # Draw the TGraphErrors
    c = ROOT.TCanvas()
    for i,key in enumerate(gr):
	if i == 0:
    		gr['%s'%key].SetTitle("Cross Section")
		gr['%s'%key].SetMarkerStyle(20)
		gr['%s'%key].SetMarkerColor(1)
		c.SetLogy(1)
		gr['%s'%key].GetHistogram().SetAxisRange(550., 5500.,"X");
		gr['%s'%key].GetHistogram().SetAxisRange(0.00001, 0.5,"Y");
		#gr['%s'%key].GetYaxis().SetRangeUser(0.00001, 0.5)
    		gr['%s'%key].GetXaxis().SetTitle("Mass")
    		gr['%s'%key].GetYaxis().SetTitle("xSec")	
    		gr['%s'%key].Draw("AP")
    		c.Update()
	else:
                gr['%s'%key].SetMarkerStyle(i+20)
                gr['%s'%key].SetMarkerColor(i+1)
		gr['%s'%key].Draw("PSAME")
		c.Update()
	legend.AddEntry(gr['%s'%key], "cross section %s" %key, "p")
	legend.Draw("SAME")
	c.Update()
    for ext in ['png','pdf']:
    	c.SaveAs("cross_section."+ext)

# Usage example
files = ["GenAnalq_M1000_Lambda0p1.root", "GenAnalq_M2000_Lambda0p1.root","GenAnalq_M1000_Lambda1p0.root"]
plotFromDifferentTrees(files)
