'''
RSG & CUR step : job_RSGCUR.sh
    RSG via pyxtal and CUR selection by SOAP kernel via Dscribe
    input : RSGCUR.json
        RSGCUR.json : {CUR hyps, SOAP hyps}
    output : [{composition}_0.pickle], dft_targets_0.txt

DFT step : job_DFT.sh
    DFT via ASE array job
    input : dft_targets_{num}.txt
    output : [{composition}_dft_{num}.pickle], log.txt

MLE step : job_MLE.sh
    update gp update db and MLE via FLARE
    input : gp_{num-1}.pickle, log.txt
        gp{num}.pickle : FLARE gp model
    output : gp{num}.pickle, log.txt

MD & GPR step : job_MDGPR.sh
    sample new structure by MD via ASE and GPR via FLARE array job
    input : md_targets_{num}.txt, gp_{num}.pickle
    output : [{composition}_gpr_{num+1}.pickle], [{composition}_{num+1}.pickle], [md_log_{num}.txt]

OTF step : job_OTF.sh
    collect results and decide wheter to perform DFT
    input : md_targets_{num}.txt, log.txt
    output : md_targets_{num+1}.txt, dft_targets_{num+1}.txt, log.txt

log.txt : OTF 0 start_time YMDHMS
          DFT 0 start_time YMDHMS
          MLE 0 start_time YMDHMS start_hyps 1 2 3 4 5 end_hyps 1 2 3 4 5
          MDGPR 0 start_time YMDHMS
          OTF 1 start_time YMDHMS True
          DFT 1 ...

dft_targets_{num}.txt : path to {composition}_0.pickle s \n
                            add atom targets number
md_targets_{num}.txt : path to {composition}_{num}.pickle & path to {composition}_{method}_{num}.pickle s


MHM step : job_MHM.sh
    RSG via pyxtal and MHM via ASE array by MGP via LAMMPS

TI & TPT step : job_TITPT.sh
    TI and TPT by MGP via LAMMPS
'''
