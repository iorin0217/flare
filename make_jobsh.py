'''
RSG & CUR step : job_RSGCUR.sh
    RSG via pyxtal and CUR selection by SOAP kernel via Dscribe
    input : RSGCUR.json
        RSGCUR.json : {CUR hyps, SOAP hyps}
    output : [{composition}_0.pickle], dft_targets_0.txt

DFT step : job_DFT.sh
    DFT via ASE array job
    input : dft_targets_{num}.txt
    output : [{composition}_dft_{num}.pickle]

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

log.txt : OTF 0 end_time YMD HMS True
          DFT 0 end_time YMD HMS
          MLE 0 end_time YMD HMS start_hyps [1, 2, 3, 4, 5] end_hyps [1, 2, 4, 6, 8]
          MDGPR 0 end_time YMD HMS
          OTF 1 end_time YMD HMS False
          DFT 1 ...

dft_targets_{num}.txt : expdir/{composition0}/{composition0}_{num}/{composition0}_{num}.pickle
                        [0, 2, 4, 5]
                        expdir/{composition1}/{composition1}_{num}/{composition1}_{num}.pickle
                        [3, 8, 12]
                        ...

md_targets_{num}.txt : expdir/{composition0}/{composition0}_{num}/{composition0}_{num}.pickle
                       expdir/{composition0}/{composition0}_{num}/{composition0}_{method}_{num}.pickle
                       expdir/{composition1}/{composition1}_{num}/{composition1}_{num}.pickle
                       expdir/{composition1}/{composition1}_{num}/{composition1}_{method}_{num}.pickle
                       ...


MHM step : job_MHM.sh
    RSG via pyxtal and MHM via ASE array by MGP via LAMMPS

TI & TPT step : job_TITPT.sh
    TI and TPT by MGP via LAMMPS
'''
