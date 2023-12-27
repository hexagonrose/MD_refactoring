import shutil as shu
import os
from module_vasp import *
import numpy as np

# def get_MSD(xdat_dir,Atoms):
#     Nspecies=len(Atoms)
#     NION=sum(Atoms)
#     vec=np.zeros((3,3))
#     with open(xdat_dir,'r') as O:
#         xdat=O.readlines()
#     for i in range(3):
#         vec[i,:]=np.array([float(j) for j in xdat[i+2].split()])
#     tot_len=(len(xdat)-7)//(NION+1)
    
#     initial=np.zeros((NION,3))
#     tot_coord=np.zeros(((tot_len-1)*NION,3))
#     MSD=np.zeros((tot_len,Nspecies+2))
#     MSD[:,0]=np.arange(0,tot_len)

#     for i in range(8,8+NION):
#         initial[i-8,:]=[float(coord) for coord in xdat[i].split()]
#     ##initial=initial*vec # Cartesian coord
#     for step in range(0,tot_len-1):
#         for atom_coord in range(0,NION):
#             tot_coord[step*NION+atom_coord,:]= \
#             [float(coord) for coord in xdat[(step+1)*(NION+1)+8+atom_coord].split()]
#     ##tot_coord=tot_coord*vec
#     for i in range(0,tot_len-1):
#         tot_coord[i*NION:(i+1)*NION]=np.abs(tot_coord[i*NION:(i+1)*NION]-initial)
#     tot_coord[np.nonzero(tot_coord>=0.5)]-=1

#     for i in range(0,tot_len-1):
#         start_count=0
#         for e_at,at_num in enumerate(Atoms):
#             MSD[i+1,e_at+1]=np.sum(np.sum(np.power(np.matmul(tot_coord[i*NION+start_count:i*NION+start_count+at_num],vec),2),axis=1))/at_num
#             start_count+=at_num
#         for e_at,at_num in enumerate(Atoms):
#             MSD[i+1,-1]+=MSD[i+1,e_at+1]*at_num
#         MSD[i+1,-1]/=NION

#     # for i in range(0,tot_len-1):
#     #     start_count=0
#     #     for e_at,at_num in enumerate(Atoms):
#     #         MSD[i+1,e_at+1]=np.sum(np.sum(np.power(np.matmul(tot_coord[i*NION+start_count:i*NION+start_count+at_num]-initial[start_count:start_count+at_num],vec),2),axis=1))
#     #         start_count+=at_num
#     #     MSD[i+1,-1]=np.sum(MSD[i+1,1:-1])
#     #     for e_at,at_num in enumerate(Atoms):
#     #         MSD[i+1,e_at+1]/=at_num
#     #     MSD[i+1,-1]/=NION
#     return MSD

def save_xdat_unwrap(xdat_dir):
    vec=np.zeros((3,3))
    with open(xdat_dir,'r') as O:
        xdat=O.readlines()
    Atoms=[int(i) for i in xdat[6].split()]
    NION=sum(Atoms)
    for i in range(3):
        vec[i,:]=np.array([float(j) for j in xdat[i+2].split()])
    tot_len=(len(xdat)-7)//(NION+1)

    tot_coord=np.zeros((tot_len,NION,3))
    rel_to_move=np.zeros((tot_len,NION,3))

    for step in range(0,tot_len):
        for atom_coord in range(0,NION):
            tot_coord[step,atom_coord,:]= \
            [float(coord) for coord in xdat[step*(NION+1)+8+atom_coord].split()]
    for i in range(1,tot_len-1):
        rel_to_move[i:,:,:]+=np.round(tot_coord[i,:,:]-tot_coord[i-1,:,:])
    tot_coord+=rel_to_move
    tot_coord=np.matmul(tot_coord.reshape(-1,3),vec).reshape(tot_len,-1,3)
    with open(xdat_dir+'_unwrapped.xyz','w') as O:
        for i in range(7):
            O.write(xdat[i])
        for i in range(tot_len):
            O.write(f"Cartesian configuration= {i+1}\n")
            for j in range(NION):
                O.write(f"{tot_coord[i,j,0]} {tot_coord[i,j,1]} {tot_coord[i,j,2]}\n")

    return vec,Atoms,tot_coord


def read_xdat_unwrap(xdat_dir):
    vec=np.zeros((3,3))
    with open(xdat_dir,'r') as O:
        xdat=O.readlines()
    Atoms=[int(i) for i in xdat[6].split()]
    NION=sum(Atoms)
    for i in range(3):
        vec[i,:]=np.array([float(j) for j in xdat[i+2].split()])
    tot_len=(len(xdat)-7)//(NION+1)

    tot_coord=np.zeros((tot_len,NION,3))

    for step in range(0,tot_len):
        for atom_coord in range(0,NION):
            tot_coord[step,atom_coord,:]= \
            [float(coord) for coord in xdat[step*(NION+1)+8+atom_coord].split()]

    return vec,Atoms,tot_coord

def calc_MSD(sim_info):
    vec=sim_info[0]
    Atoms=sim_info[1]
    Nspecies=len(Atoms)
    NION=sum(Atoms)
    tot_coord=sim_info[2]
    tot_len=tot_coord.shape[0]
    ### TODO: sim_time should be input & msd_len should be configurable
    sim_time=2e-3
    msd_len=1500
    msd_interval=5

    if tot_len>=msd_len+msd_interval:
        MSD=np.zeros((msd_len,Nspecies+2))
        MSD[:,0]=np.arange(0,msd_len)*sim_time

        ensemble_N=(tot_len-msd_len)//msd_interval
        for ensemble in range(ensemble_N):
            initial=tot_coord[ensemble*msd_interval,:,:].copy()
            for step in range(1,msd_len):
                start_count=0
                for e_at,at_num in enumerate(Atoms):
                    MSD[step,e_at+1]+=np.sum(np.power(tot_coord[ensemble*msd_interval+step,start_count:start_count+at_num,:]-initial[start_count:start_count+at_num,:],2))/at_num
                    start_count+=at_num
        for e_at,at_num in enumerate(Atoms):
            MSD[:,-1]+=MSD[:,e_at+1]*at_num
        MSD[:,-1]/=NION
        MSD[:,1:]/=ensemble_N

    else:
        MSD=np.zeros((tot_len,Nspecies+2))
        MSD[:,0]=np.arange(0,tot_len)*sim_time

        initial=tot_coord[0,:,:].copy()
        for step in range(1,tot_len):
            start_count=0
            for e_at,at_num in enumerate(Atoms):
                MSD[step,e_at+1]+=np.sum(np.power(tot_coord[step,start_count:start_count+at_num,:]-initial[start_count:start_count+at_num,:],2))/at_num
                start_count+=at_num
                for e_at,at_num in enumerate(Atoms):
                    MSD[step,-1]+=MSD[step,e_at+1]*at_num
                MSD[step,-1]/=NION
    writer=""
    for j in MSD:
        for k in j:
            writer+=f"{k} "
        writer+='\n'
    with open("XDATCAR.msd",'w') as O:
        O.write("#Time ")
        for i in range(Nspecies):
            O.write(f"Atom{i+1} ")
        O.write("Average\n")
        O.write(writer)
    return MSD

def pearson(MSD):
    msd=MSD[:,-1]

    xy_dat=np.zeros((len(msd)+1,2),dtype=float)
    xy_dat[:,0]=np.arange(0,len(msd)+1).reshape(1,-1)
    for e,i in enumerate(msd):
        xy_dat[e+1,1]=i
    #print(xy_dat)
    return  np.sum((xy_dat[:,0]-np.mean(xy_dat[:,0]))*(xy_dat[:,1]-np.mean(xy_dat[:,1]))) /np.sqrt(np.sum(  np.power(xy_dat[:,0]-np.mean(xy_dat[:,0]),2.)  )) /np.sqrt(np.sum(  np.power(xy_dat[:,1]-np.mean(xy_dat[:,1]),2.)  ))


def find_Tm(continue_dir, input_Tm, vasp_config, vasp_version, log_dir):
    Tm = 0
    done = False
    find_done = False
    iter = 0
    with open(log_dir, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if 'K melting is done' in line:
                Tm = int(line.split()[0])
            elif 'Predicted Tm:' in line:
                Tm = int(line.split()[2])
                find_done = True
            elif 'iter final relaxation is done' in line:
                iter = int(line.split()[0])
            elif 'MD setting is done' in line:
                done = True
    iter += 2

    if done:
        return 1
    
    if Tm != 0:
        continue_dir = os.getcwd()+f"/Step1_T_{Tm}"
        start_temp = Tm - input_Tm['StepSize_T']
    else:
        start_temp = input_Tm['Start_T']

    ###find MSD at each T
    if find_done == False:
        for i in range(start_temp, input_Tm['End_T']-input_Tm['StepSize_T'], -input_Tm['StepSize_T']):
            os.makedirs(f"Step1_T_{i}", exist_ok=True)
            os.chdir(f"Step1_T_{i}")
            copy_inputs(continue_dir, './', ['INCAR','KPOINTS','POTCAR','CONTCAR'])
            shu.move("./CONTCAR","./POSCAR")
            scale_velo("./POSCAR", i)
            edit_INCAR('./INCAR',{"NSW":str(input_Tm['duration']),
                        "IBRION":"0","MDALGO":"2","SMASS":"0","TEBEG":str(i),
                        "TEEND":str(i),"ISIF":"2","POTIM":"2"})
            run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config[vasp_version])
            if os.path.isfile("./XDATCAR_unwrapped.xyz"):
                sim_info=read_xdat_unwrap("./XDATCAR_unwrapped.xyz")
            else:
                sim_info=save_xdat_unwrap("./XDATCAR")
            msd=calc_MSD(sim_info)
            if msd[-1,-1]/msd[-1,0]/6/2*1e-8 >= 4e-9:   # reference 필요.
                continue_dir=os.getcwd()
                shu.copy("./INCAR","../../Inputs/INCAR")
                Tm=i
                velocity_from_poscar='/'.join([os.getcwd(),"CONTCAR"])
                os.chdir('../')
                with open(log_dir, 'a') as s:
                    s.write('%s K melting is done\n'%i)
            else:
                os.chdir("../")
                with open(log_dir, 'a') as s:
                    s.write('%s K melting is done\n'%i)
                    s.write('Predicted Tm: %s\n'%Tm)
                break
        else:
            with open(log_dir, 'a') as s:
                s.write('Tm is less than 1000K\n'%i)
                s.write('Predicted Tm: 1000\n')

    if find_done == True and Tm == 0:
        os.makedirs("Step1_T_4000", exist_ok=True)
        os.chdir("Step1_T_4000")
        copy_inputs(continue_dir, './', ['INCAR','KPOINTS','POTCAR','CONTCAR'])
        shu.move("./CONTCAR","./POSCAR")
        scale_velo("./POSCAR", 4000)
        edit_INCAR('./INCAR',{"NSW":str(input_Tm['duration']),
                    "IBRION":"0","MDALGO":"2","SMASS":"0","TEBEG":"4000",
                    "TEEND":"4000","ISIF":"2","POTIM":"2"})
        run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config[vasp_version])
        continue_dir=os.getcwd()
        shu.copy("./INCAR","../../Inputs/INCAR")
        Tm = 4000
        os.chdir("../")
        with open(log_dir, 'a') as s:
            s.write('%s K melting is done\n'%4000)
            s.write('Predicted Tm: %s\n'%Tm)


    # Final volume relax
    continue_dir = os.getcwd()+f"/Step1_T_{Tm}"
    #pressures = np.array([float(p_line) for p_line in grep_nth_item("external", continue_dir+"/OUTCAR",3)[250:]])
    pressures = np.array([float(p_line) for p_line in grep_nth_item("external", continue_dir+"/OUTCAR",3)[:]])
    if  np.abs(np.mean(pressures))>=30 :
        ### post-relax after finding Tm if needed
        if iter != 2:
            continue_dir = os.getcwd()+f"/Step{iter}_relax"
            iter += 1
        velocity_from_poscar='/'.join([continue_dir,"CONTCAR"])
            
        for i in range(iter, 2+input_Tm['maximum_post_rlx']):
            os.makedirs(f'Step{i}_relax', exist_ok=True)
            os.chdir(f'Step{i}_relax')
            copy_inputs(continue_dir,'./',['INCAR','KPOINTS','POTCAR','CONTCAR'])
            shu.move('./CONTCAR','./POSCAR')
            
            edit_INCAR("./INCAR",{"IBRION":"2","NSW":"10","ISIF":"7","SMASS":"#","MDALGO":"#","POTIM":"#","TEBEG":"#","TEEND":"#"})
            run_vasp(vasp_config['mpicommand'],vasp_config['num_tasks'],vasp_config[vasp_version])
            shu.move("./OUTCAR",'./OUTCAR_rlx')
            
            shu.move('./CONTCAR','./POSCAR')
            get_velo(velocity_from_poscar, "./POSCAR")
            edit_INCAR("./INCAR",{"IBRION":"0","NSW":"500","ISIF":"2","MDALGO":"2","SMASS":"0","TEBEG":f"{Tm}","TEEND":f"{Tm}","POTIM":"2"})
            run_vasp(vasp_config['mpicommand'],vasp_config['num_tasks'],vasp_config[vasp_version])
            continue_dir=os.getcwd()
            pressures=np.array([float(p_line) for p_line in grep_nth_item("external","./OUTCAR",3)[250:]])
            if  np.abs(np.mean(pressures))<=30 :
                shu.copy("./CONTCAR",'../../Inputs/POSCAR_to_melt')
                os.chdir("../")
                with open(log_dir, 'a') as s:
                    s.write('%s iter final relaxation is done\n'%(i-2))
                    s.write('Pressure is equilibrated\n')
                break
            else:
                velocity_from_poscar='/'.join([os.getcwd(),"CONTCAR"])
                os.chdir("../")
                with open(log_dir, 'a') as s:
                    s.write('%s iter final relaxation is done\n'%(i-2))
            #TODO if converged pressure, break
    else:
        shu.copy('Step1_T_'+str(Tm)+"/CONTCAR",'../Inputs/POSCAR_to_melt')
        #os.chdir("../")
        with open(log_dir, 'a') as s:
            s.write('Pressure is equilibrated\n')

def initial_rlx(initial_config, vasp_config, vasp_version, log_dir):
    # check LOG and find last relaxation index
    iter = 0
    done = False
    with open(log_dir, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if 'Pressure is equilibrated' in line:
                done = True
            elif 'iter initial relaxation is done' in line:
                iter = int(line.split()[0])
    if done:
        return 1
        
    os.makedirs('Step0_relax', exist_ok=True)
    os.chdir('Step0_relax')
    if iter != 0:
        os.chdir(str(iter))
        initial_config["INCAR_dir"]=os.getcwd()+"/INCAR"
        initial_config["KP_dir"]=os.getcwd()+"/KPOINTS"
        initial_config["POT_dir"]=os.getcwd()+"/POTCAR"
        initial_config["POS_dir"]=os.getcwd()+"/CONTCAR"
        os.chdir('../')
        iter += 1

    # Main relaxation loop
    velocity_from_poscar = initial_config["POS_dir"]
    for i in range(iter, 5):
        os.makedirs(str(i), exist_ok=True)
        os.chdir(str(i))
        shu.copy(initial_config["INCAR_dir"],'./INCAR')
        shu.copy(initial_config["KP_dir"],'./KPOINTS')
        shu.copy(initial_config["POT_dir"],'./POTCAR')
        shu.copy(initial_config["POS_dir"],'./POSCAR')
        edit_INCAR("./INCAR",{"IBRION":"2","NSW":"10","ISIF":"7","SMASS":"#","MDALGO":"#","POTIM":"#","TEBEG":"#","TEEND":"#","PSTRESS":"30"})

        run_vasp(vasp_config['mpicommand'],vasp_config['num_tasks'],vasp_config[vasp_version])
        shu.move("./CONTCAR","./POSCAR")

        get_velo(velocity_from_poscar,"./POSCAR")
        shu.move('./OUTCAR','./OUTCAR_rlx')
        edit_INCAR("./INCAR",{"IBRION":"0","NSW":"500","ISIF":"2","SMASS":"0","MDALGO":"2","POTIM":"2","TEBEG":"4500","TEEND":"4500","PSTRESS":""})
        run_vasp(vasp_config['mpicommand'],vasp_config['num_tasks'],vasp_config[vasp_version])
        
        ##get pressure after 0.5 ps to equilibriate
        pressures=np.array([float(p_line) for p_line in grep_nth_item("external","./OUTCAR",3)[250:]]) 
        if np.average(pressures)<=50 and np.average(pressures)>=-10:
            os.chdir("../")
            with open(log_dir, 'a') as s:
                s.write('%s iter initial relaxation is done\n'%i)
                s.write('Pressure is equilibrated\n')
            break
        else:
            velocity_from_poscar='/'.join([os.getcwd(),"CONTCAR"])
            initial_config["INCAR_dir"]=os.getcwd()+"/INCAR"
            initial_config["KP_dir"]=os.getcwd()+"/KPOINTS"
            initial_config["POT_dir"]=os.getcwd()+"/POTCAR"
            initial_config["POS_dir"]=os.getcwd()+"/CONTCAR"
            os.chdir("../")
            with open(log_dir, 'a') as s:
                s.write('%s iter initial relaxation is done\n'%i)
            
            ### TODO: if not converged?
    os.chdir('../')



    
