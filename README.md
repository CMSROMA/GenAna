# Generator level analysis on LQ events with CMSSW

```
mkdir GenAnalysis
cd GenAnalysis

scram project -n CMSSW_10_6_28_LQAna CMSSW CMSSW_10_6_28
cd CMSSW_10_6_28_LQAna/src
cmsenv

git clone git@github.com:CMSROMA/GenAna.git CMSROMA/GenAna
cd CMSROMA/GenAna/

scram b
```

Edit config file "python/GenAnalq_cfg.py" as needed (i.e. input file).

Run analysis:
```
cmsRun python/GenAnalq_cfg.py
```

Outputs:
```
GenAnalq.root (containing both histograms and tree)
```

## To run over several GEN files:

Create a alist of GEN files (i.e. list_singleLQ_GEN.txt)
```
root://eoscms///eos/cms/store/group/phys_exotica/lq-LQ-lq/test_GEN_1/singleLQ_13TeV_Pow_Herwig7_M2000_Lambda1p0_GEN.root
root://eoscms///eos/cms/store/group/phys_exotica/lq-LQ-lq/test_GEN_1/singleLQ_13TeV_Pow_Herwig7_M3000_Lambda1p0_GEN.root
...

```

Run analysis on all samples:
```
python makeANAfromGEN.py -c python/GenAnalq_cfg.py -i list_singleLQ_GEN.txt -t /tmp/santanas/ --outputDir `pwd`/TestOutput
```

Make plots from histograms for all samples:

Edit makePlots_singleLQ.C with input root files, specify output folder (www area on afs) and then run:
```
root -l -b makePlots_singleLQ.C
```


