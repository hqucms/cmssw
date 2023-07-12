filein=$1
[[ -z ${filein} ]] && filein="file:/eos/cms/store/cmst3/group/hgcal/CMG_studies/psilva/TB2023/output_numEvent100.root"

fileout=$2
[[ -z ${fileout} ]] && fileout="hgcalnano_testcmsdriver.root"

cmsDriver.py NANO \
    -s USER:DPGAnalysis/HGCalNanoAOD/hgcRecHits_cff.hgctbRecHitsTask \
    --datatier NANOAOD \
    --eventcontent NANOAOD \
    --filein ${filein} \
    --fileout ${fileout} \
    -n -1 \
    --nThreads 4 \
    --conditions auto:phase2_realistic_T21 \
    --geometry Extended2026D95 \
    --era Phase2C17I13M9 \
    --python_filename testcmsdriver_cfg.py \
    --no_exec

cmsRun -j testcmsdriver_cfg.py
