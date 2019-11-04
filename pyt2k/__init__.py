import os
#os.environ["LD_LIBRARY_PATH"] = ":".join((os.path.abspath(os.path.dirname(__file__)), os.environ["LD_LIBRARY_PATH"]))
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True # Prevent ROOT from hi-jacking --help
ROOT.gSystem.AddDynamicPath(os.path.abspath(os.path.dirname(__file__)))

def _loadlib(fname):
    ret = ROOT.gSystem.Load(fname)
    return

#_loadlib("libT2KDataFrameMakeProject")
_loadlib("libT2KDataSource")

class T2K:
    T2KDataSource = ROOT.T2K.T2KDataSource

# always store weights in histograms
#ROOT.TH1.SetDefaultSumw2(True)
