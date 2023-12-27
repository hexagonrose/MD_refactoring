import shutil as shu
import subprocess as sub
import os
from module_vasp import *
import numpy as np

def conv_check(EFS_ref, EFS_check, NION, criteria):
    if abs(EFS_ref[0]-EFS_check[0])/NION <= criteria["E_tol"] and \
       np.sqrt(np.max(np.sum(np.power(EFS_ref[1]-EFS_check[1],2),axis=1))) <= criteria["F_tol"] and \
       abs(EFS_ref[2]-EFS_check[2])<=criteria["S_tol"] :
        return 1
    else:
        return 0

def do_conv_test(mat_dir, config, vasp_config, log_dir):
    is_kptest=config['kp_use']
    is_cutoff=config['cutoff_use']

    if is_kptest.upper()=='AUTO' or is_cutoff.upper()=='AUTO':
        os.makedirs('reference', exist_ok=True)
        os.chdir('reference')

        # Calculation of reference
        copy_inputs(from_dir=mat_dir+"/premelt", to_dir='./', copy_list=["INCAR","KPOINTS","POTCAR","CONTCAR"])
        shu.move("CONTCAR", "POSCAR")
        edit_INCAR(target="./INCAR", to_change={"ENCUT":config["reference_cutoff"], "PREC":"normal", \
                "NSW":"0", "POTIM":"#", "TEBEG":"#", "TEEND":"#", "SMASS":"#", "MDALGO":"#"})
        if vasp_config['npar']//2>=1 and vasp_config['npar']%2==0:
            edit_INCAR(target="./INCAR", to_change={"KPAR":"2", "NPAR":str(vasp_config['npar']//2)})
        edit_KP(target="./KPOINTS", num=config["reference_KP"])
        run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config['std'])
        with open("./POSCAR",'r') as O:
            ions = O.readlines()[6].split()
        NION = sum([int(n_tmp) for n_tmp in ions])
        EFS_ref = grep_EFS("./OUTCAR", NION)
        os.chdir('../')

        if is_kptest.upper()=='AUTO':
            os.makedirs('kptest', exist_ok=True)
            os.chdir('kptest')
            for current_KP in ['G', 'B']:
                os.makedirs(current_KP, exist_ok=True)
                os.chdir(current_KP)
                copy_inputs(from_dir="../../reference", to_dir="./", copy_list=["INCAR","KPOINTS","POTCAR","POSCAR"])
                edit_INCAR(target="./INCAR", to_change={"KPAR":"#", "NPAR":str(vasp_config['npar'])})
                edit_KP(target='./KPOINTS', special_kp=current_KP)
                if current_KP=="G":
                    run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config['gam'])
                else:
                    run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config['std'])

                EFS_check = grep_EFS("./OUTCAR", NION)
                converged = conv_check(EFS_ref, EFS_check, NION, config)
                if converged:
                    shu.copy("KPOINTS", "../../../Inputs/")
                    shu.copy("INCAR", "../../../Inputs/")
                    os.chdir('../')
                    with open(log_dir, 'a') as s:
                        s.write('Converged KP: %s\n'%current_KP)
                    break
                else:
                    os.chdir('../')
            else:
                with open(log_dir, 'a') as s:
                    s.write('KP is not converged\n')
                return 1 #parameter is not converged and exit the function.

            os.chdir("../")

        if is_cutoff.upper()=='AUTO':
            os.makedirs('cutoff', exist_ok=True)
            os.chdir('cutoff')
            #TODO get ENMAX from POTCAR
            enmax = 250
            e_range =  range(500, enmax-1, -config["cutoff_stepsize"])
            cutoff_conv = -1
            for current_cutoff in e_range:
                os.makedirs(str(current_cutoff), exist_ok=True)
                os.chdir(str(current_cutoff))
                copy_inputs(from_dir="../../reference", to_dir="./", copy_list=["POTCAR","POSCAR"])
                copy_inputs(from_dir="../../../Inputs", to_dir="./", copy_list=["INCAR","KPOINTS"])
                edit_INCAR("./INCAR", {"ENCUT":str(current_cutoff)})
                if current_KP == 'G':
                    run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config['gam'])
                else:
                    run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config['std'])
                EFS_check = grep_EFS("./OUTCAR", NION)
                converged = conv_check(EFS_ref, EFS_check, NION, config)
                if converged:
                    cutoff_conv += 1
                    shu.copy("./INCAR", "../../../Inputs/INCAR")
                    os.chdir('../')
                else:
                    os.chdir('../')
                    with open(log_dir, 'a') as s:
                        s.write('Converged ENCUT: %s\n'%(current_cutoff+config["cutoff_stepsize"]))
                    break
            else:
                with open(log_dir, 'a') as s:
                    s.write('ENCUT is not converged\n')
                # 20231212 hk ENCUT 복구
                edit_INCAR(target=mat_dir+"/Inputs/INCAR", to_change={"ENCUT":"550"})
                print("new code worked")
            # if cutoff_conv!=1:
            #     shu.copy(str(int(config["reference_cutoff"])-config["cutoff_stepsize"]*cutoff_conv)+"/INCAR","./INCAR_converged")
            # else:
            #     shu.copy("..//INCAR","./INCAR_converged")
            os.chdir("../")
