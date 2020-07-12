'''
after PRE step, DFT(0) -> MLE(0) -> ...
imaginary, the first line of log.txt is "OTF 0"

log.txt : OTF 0 end_time YMD HMS False
          DFT 0 end_time YMD HMS
          MLE 0 end_time YMD HMS start_hyps [1, 2, 3, 4, 5] end_hyps [1, 2, 4, 6, 8]
          MDGPR 0 end_time YMD HMS
          OTF 1 end_time YMD HMS True
          DFT 1 ...

DFT step : job_DFT_{num}.sh
    DFT via ASE array job
    input : dft_targets_{num}.txt
    output : [{composition}_dft_{num}.pickle]

MLE step : job_MLE_{num}.sh
    update gp update db and MLE via FLARE
    input : gp_{num-1}.pickle, log.txt
        gp{num}.pickle : FLARE gp model
    output : gp{num}.pickle, log.txt
    submit : MDGPR_{num} array -> OTF_{num+1}

MD & GPR step : job_MDGPR.sh
    sample new structure by MD via ASE and GPR via FLARE array job
    input : md_targets_{num}.txt, gp_{num}.pickle
    output : [{composition}_gpr_{num+1}.pickle], [{composition}_{num+1}.pickle], [md_log_{num}.txt]

OTF step : job_OTF.sh
    collect results and decide wheter to perform DFT
    input : md_targets_{num}.txt, log.txt
    output : md_targets_{num+1}.txt, dft_targets_{num+1}.txt, log.txt
    submit : DFT_{num} array -> MLE_{num}

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
'''


if __name__ == "__main__":
    # python make_jobsh.py token user num
    # @ want to save dir
    import sys
    token = sys.argv[1]
    user = sys.argv[2]
    num = sys.argv[3]
    header_common = ["#PBS -q h-regular", "#PBS -l select=1:ncpus=36:mpiprocs=32:ompthreads=1",
                     f"#PBS -W group_list={token}", "#PBS -j oe"]
    header_DFT = ["#PBS -l walltime=1:30:00", "#PBS -ry"]
    header_MLE = ["#PBS -l walltime=1:00:00"]
    header_MDGPR = ["#PBS -l walltime=0:30:00", "#PBS -ry"]
    header_OTF = ["#PBS -l walltime=0:30:00"]
    config = ["cd $PBS_O_WORKDIR",
              f"source /lustre/{token}/{user}/.bash_profile"]
    for i in range(int(num)):
        DFT = [
            f'target=`sed -n $((2*PBS_ARRAY_INDEX-1))"P" dft_targets_{i}.txt`', f"python3 /lustre/{token}/{user}/flare/DFT.py $target"]
        DFT_sh = header_common + header_DFT + \
            [f"#PBS -N DFT_{i}"] + config + DFT
        print(*DFT_sh, sep="\n", end="\n", file=open(f"job_DFT_{i}.sh", "w"))
        MLE = [
            f"python3 /lustre/{token}/{user}/flare/MLE.py $PBS_O_WORKDIR/gp_{i-1}.pickle $PBS_O_WORKDIR/log.txt"]
        MLE_sh = header_common + header_MLE + \
            [f"#PBS -N MLE_{i}"] + config + MLE
        print(*MLE_sh, sep="\n", end="\n", file=open(f"job_MLE_{i}.sh", "w"))
        MDGPR = [f'target=`sed -n $(2*$PBS_ARRAY_INDEX-1)"P" md_targets_{i}.txt`', f'engine=`sed -n $(2*$PBS_ARRAY_INDEX)"P" md_targets_{i}.txt`',
                 f"python3 /lustre/{token}/{user}/flare/MDGPR.py $target $engine $PBS_O_WORKDIR/gp_{i}.pickle"]
        MDGPR_sh = header_common + header_MDGPR + \
            [f"#PBS -N MDGPR_{i}"] + config + MDGPR
        print(*MDGPR_sh, sep="\n", end="\n",
              file=open(f"job_MDGPR_{i}.sh", "w"))
        OTF = [f'target=`sed -n $PBS_ARRAY_INDEX"P" dft_targets_{i}.txt`',
               f"python3 /lustre/{token}/{user}/flare/OTF.py $PBS_O_WORKDIR/md_targets_{i}.txt $PBS_O_WORKDIR/log.txt"]
        OTF_sh = header_common + header_OTF + \
            [f"#PBS -N OTF_{i}"] + config + OTF
        print(*OTF_sh, sep="\n", end="\n", file=open(f"job_OTF_{i}.sh", "w"))
