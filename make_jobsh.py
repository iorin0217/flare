'''
RSG & CUR step : job_RSGCUR.sh
    RSG via pyxtal and CUR selection by SOAP kernel via Dscribe
    input : RSGCUR.json
        RSGCUR.json : {Gibbs hyps, CUR hyps, SOAP hyps}
    output : structures_init.npy
        structures_init.npy : [structure_init : ASE structure]

DFT step : job_DFT.sh
    DFT via ASE array job
    input : structures_{init/add}.npy
    output : datas.npy
        datas.npy : [data.npy : FLARE structure]

MLE step : job_MLE.sh
    update gp db and MLE via FLARE
    input : gp{num}.pickle, datas.npy, log.txt
        gp{num}.pickle : FLARE gp model
        log.txt : 
    output : gp{num}.pickle, log.txt

MD & GPR step : job_MDGPR.sh
    sample new structure by MD via ASE and GPR via FLARE array job
    input : datas.npy, gp{num}.pickle
    output : datas.npy

OTF step : job_OTF.sh
    decide wheter to perform DFT
    input : datas.npy
    output : structure_add.npy, log.txt

MHM step : job_MHM.sh
    RSG via pyxtal and MHM via ASE array by MGP via LAMMPS

TI & TPT step : job_TITPT.sh
    TI and TPT by MGP via LAMMPS
'''
