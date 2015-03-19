#
# usage: %prog [opts] --cfg cmssw.py dataset 

# GJets
./submitBatchSinglePho.py --cfg singlePhoAnaBATCH.py GJets_HT-100to200
./submitBatchSinglePho.py --cfg singlePhoAnaBATCH.py GJets_HT-200to400 
./submitBatchSinglePho.py --cfg singlePhoAnaBATCH.py GJets_HT-400to600 
./submitBatchSinglePho.py --cfg singlePhoAnaBATCH.py GJets_HT-600toInf 

# RS Graviton
./submitBatchSinglePho.py --cfg singlePhoAnaBATCH.py RSGravToGG_kMpl-01_M-1500  
./submitBatchSinglePho.py --cfg singlePhoAnaBATCH.py RSGravToGG_kMpl-02_M-3000
./submitBatchSinglePho.py --cfg singlePhoAnaBATCH.py RSGravToGG_kMpl-001_M-1500  
./submitBatchSinglePho.py --cfg singlePhoAnaBATCH.py RSGravToGG_kMpl-01_M-3000
./submitBatchSinglePho.py --cfg singlePhoAnaBATCH.py RSGravToGG_kMpl-001_M-750   
./submitBatchSinglePho.py --cfg singlePhoAnaBATCH.py RSGravToGG_kMpl-02_M-1500
