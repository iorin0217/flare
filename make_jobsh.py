'''
RSG & CUR step : job_RSGCUR.sh
    RSG via pyxtal and CUR selection by SOAP kernel via Dscribe
    input : RSGCUR.json
        RSGCUR.json : {CUR hyps, SOAP hyps}
    output : [{composition}_0.pickle], dft_targets_0.txt
        dft_targets_{num}.txt : path to {composition}_0.pickle s \n
                            add atom targets number

DFT step : job_DFT.sh
    DFT via ASE array job
    input : dft_targets_{num}.txt
    output : [{composition}_0_flare.pickle : FLARE structure], log.txt
        log.txt : write "DFT {num}"

MLE step : job_MLE.sh
    update gp update db and MLE via FLARE
    input : gp_{num}.pickle, log.txt
        gp{num}.pickle : FLARE gp model
        log.txt : read "DFT {num}" or "pass {num}"
    output : gp{num+1}.pickle, log.txt
        write "update {num}"

MD & GPR step : job_MDGPR.sh
    sample new structure by MD via ASE and GPR via FLARE array job
    input : md_targets_{num}.txt, gp_{num}.pickle
        md_targets_{num}.txt : path to {composition}_0_flare.pickle s
    output : md_targets_{num+1}.txt

OTF step : job_OTF.sh
    decide wheter to perform DFT
    input : md_targets_{num}.txt
    output : dft_targets_{num}.txt, log.txt


MHM step : job_MHM.sh
    RSG via pyxtal and MHM via ASE array by MGP via LAMMPS

TI & TPT step : job_TITPT.sh
    TI and TPT by MGP via LAMMPS
'''
