from CRABClient.UserUtilities import config
config = config()

# config.General.requestName = 'GEM-STA-SliceTest-2018D-Zmu'
# config.General.requestName = 'GEM-STA-SliceTest-FixEMap-VFatQual-2018D-SingleMuon-Lat2-v6'
config.General.requestName = 'GEM-STA-SliceTest-2018D-SingleMuon-GebFlags-v2-Lat1'
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'RAW2STA_SliceTest.py'

# config.Data.inputDataset = '/SingleMuon/Run2018D-ZMu-PromptReco-v2/RAW-RECO'
#config.Data.inputDataset = '/SingleMuon/Run2018D-v1/RAW'
config.Data.inputDataset = '/SingleMuon/Run2018D-v1/RAW'

config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10 # 20 units should be max 4 hours (2018D Lat)
# config.Data.unitsPerJob = 20 # 20 units should be max 4 hours, ZMu
config.Data.publication = False
# This string is used to construct the output dataset name
#config.Data.outputDatasetTag = 'GEM-STA-SliceTest-FixEMap-VFatQual-2018D-SingleMuon-Lat1-v6'
config.Data.outputDatasetTag = 'GEM-STA-SliceTest-2018D-SingleMuon-GebFlags-v2-Lat1'

# These values only make sense for processing data
#    Select input data based on a lumi mask
#config.Data.lumiMask = '2018D_321069.txt' #  'Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'
#    Select input data based on run-ranges
# config.Data.runRange = '321067-321069'  # v1
# config.Data.runRange = '323388-323399'  # Lat3
# config.Data.runRange = '322356-322492'  # Lat2 (and surrounding)
config.Data.runRange = '321908-321909'  # Lat1

#config.Data.runRange = '325057'

# config.Data.runRange = '323778'
# config.Data.runRange = '324021,323997,323980,323978,323976,323940,323857,323841,323794,323790,323778,323775,323755,323727,323726,323725,323702,323701,323700,323697,323696,323693,323526,323525,323524,323523,323491,323488,323487,323475,323474,323473,323472,323471,323470,323422,323421,323420,323419,323418,323417,323416,323415,323414,323413,323399,323398,323397,323396,323395,323394,323393,323391,323388,323367,323365,323364,323363,322633,322625,322617,322616,322605,322603,322602,322599,322510,322492,322487,322485,322484,322483,322480,322430,322407,322381,322356,322355,322348,322332,322324,322323,322322,322319,322317,322252,322222,322204,322201,322022,321990,321988,321975,321973,321961,321960,321917,321909,321908,321887,321880,321879,321834,321833,321832,321831,321820,321819,321818,321817,321816,321815,321813,321796,321795,321794,321781,321780,321778,321777,321776,321775,321774,321773,321760,321759,321758,321755,321735,321732,321731,321730,321729,321712,321710,321709,321475,321461,321457,321434,321433,321432,321431,321415,321414,321396,321393,321313,321312,321311,321310,321305,321296,321295,321294,321283,321262,321261,321233,321232,321231,321230,321221,321219,321218,321178,321177,321167,321166,321165,321164,321162,321149,321140,321138,321134,321126,321123,321122,321121,321119,321069,321068,321067,321055,321051,320996,320995,319348,319347'

# Where the output files will be transmitted to
config.Site.storageSite = 'T3_KR_KISTI'
