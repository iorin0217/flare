'''
RSG & CUR step : job_RSGCUR.sh
    RSG via pyxtal and CUR selection by SOAP kernel via Dscribe
    input : RSGCUR.json
        RSGCUR.json : {CUR hyps, SOAP hyps}
    output : [{composition}_0.pickle], dft_targets.txt
        dft_targets.txt : path to {composition}_0.pickle s

DFT step : job_DFT.sh
    DFT via ASE array job
    input : structures_{init/add}.pickle
    output : datas.pickle
        datas.pickle : [data.pickle : FLARE structure]

MLE step : job_MLE.sh
    update gp db and MLE via FLARE
    input : gp{num}.pickle, datas.pickle, log.txt
        gp{num}.pickle : FLARE gp model
        log.txt : 
    output : gp{num}.pickle, log.txt

MD & GPR step : job_MDGPR.sh
    sample new structure by MD via ASE and GPR via FLARE array job
    input : datas.pickle, gp{num}.pickle
    output : datas.pickle

OTF step : job_OTF.sh
    decide wheter to perform DFT
    input : datas.pickle
    output : structure_add.pickle, log.txt

MHM step : job_MHM.sh
    RSG via pyxtal and MHM via ASE array by MGP via LAMMPS

TI & TPT step : job_TITPT.sh
    TI and TPT by MGP via LAMMPS
'''
