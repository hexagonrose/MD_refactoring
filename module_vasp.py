import shutil as shu
import os
import subprocess as sub
import numpy as np


def copy_inputs(from_dir,to_dir,copy_list): 
    for i in copy_list:
        shu.copy(from_dir+"/"+i,to_dir)

def edit_INCAR(target,to_change):
    with open(target,'r') as O:
        lines= O.readlines()    
    
    for key_change in to_change.keys():
        for e,ll in enumerate(lines):
            if " "+key_change+" " in ll:
                if to_change[key_change]=="#":
                    if "#" != ll.split()[0]:
                        lines[e]="#"+lines[e]
                elif to_change[key_change]=="":
                    lines[e]=""
                else:
                    lines[e]="   "+key_change+" = "+str(to_change[key_change])+"\n"
                break
        else:
            if to_change[key_change]!="" and to_change[key_change]!="#":
                lines.append("   "+key_change+" = "+str(to_change[key_change])+"\n")
    

    with open(target,'w') as O:
        for i in lines:
            O.write(i)

def edit_KP(target,pack_scheme='MP',num='1 1 1',special_kp=''):
    if special_kp=='B':
        with open(target,'w') as O:
            O.write("KP\n1\nR\n0.25 0.25 0.25 1")
    elif special_kp=='G':
        with open(target,'w') as O:
            O.write("KP\n0\nMP\n1 1 1\n0 0 0")

    else:
        with open(target,'w') as O:
            O.write("KP\n0\n"+pack_scheme+"\n"+num+"\n0 0 0")

def run_vasp(mpi_command,num_tasks,vasp_bin):
    f_stdout=open("stdout.x",'w')
    call_result=sub.call([mpi_command,'-np',num_tasks,vasp_bin],stdout=f_stdout,stderr=sub.STDOUT)
    ##call_result=sub.call(["echo","vasp",'may','run'],stdout=f_stdout,stderr=sub.STDOUT)
    ##shu.copy("./POSCAR","./CONTCAR")
    f_stdout.close()

def get_velo(from_poscar,to_poscar):
    
    with open(from_poscar,'r') as O:
        from_p=O.readlines()
    atoms=sum([int(i) for i in from_p[6].split()])
    velocity=from_p[10+atoms:10+2*atoms]

    with open(to_poscar,'r') as O2:
        to_p=O2.readlines()[:9+atoms]
    with open(to_poscar,'w') as O2:
        for i in to_p:
            O2.write(i)
        O2.write('\n')
        for i in velocity:
            O2.write(i)

def scale_velo(poscar_to_edit,T):
    
    with open(poscar_to_edit,'r') as O:
        from_p=O.readlines()
    each_atom_N=[int(i) for i in from_p[6].split()]
    N_tot=sum(each_atom_N)
    velocity_tmp=from_p[10+N_tot:10+2*N_tot]
    velocity=np.zeros((N_tot,3))
    for i in range(len(velocity_tmp)):
        v_atom=velocity_tmp[i].split()
        for j in range(3):
            velocity[i,j]=float(v_atom[j])
    
    pomass=np.array([float(i.split(";")[0]) for i in grep_nth_item("POMASS","./POTCAR",2)]) #kg
    Ek=0
    start_num=0
    for species,N_at in enumerate(each_atom_N):
       Ek+=0.5*pomass[species]*np.sum(np.power(velocity[start_num:start_num+N_at,:],2))
       start_num+=N_at
    Ek*=1e10*1.6605402E-27*6.241509e18
    T_current=Ek/8.617333e-5/(N_tot-1)/3*2
    T_scale=T/T_current
    velocity*=np.sqrt(T_scale)

    with open(poscar_to_edit,'w') as O2:
        for i in from_p[:9+N_tot]:
            O2.write(i)
        O2.write('\n')
        for i in velocity:
            O2.write(f"    {i[0]} {i[1]} {i[2]}\n")

def grep_EFS(out_dir,NION):
    with open(out_dir,'r') as O:
        outcar=O.readlines()
    EFS=[]
    loop_stop=0
    Energy=0.0
    Force=[]
    Stress=0.0
    for i in range(len(outcar)-1,0,-1):
        if "external" in outcar[i]:
            Stress=float(outcar[i].split()[3])
            loop_stop+=1
        elif "free  " in outcar[i]:
            Energy=float(outcar[i].split()[4])
            loop_stop+=1
        elif "TOTAL-FORCE" in outcar[i]:
            Force=np.zeros((NION,3),dtype=float)
            for force_ll in range(i+2,i+2+NION):
                Force[force_ll-i-2,:]=np.array([float(force_component) for force_component in outcar[force_ll].split()[3:]])
            loop_stop+=1

        if loop_stop==3:
           return [Energy, Force, Stress]

def grep_nth_item(what_to_grep, read_file, item_to_get):
    with open(read_file,'r') as O:
        contents=O.readlines()
    
    to_return=[]
    for i in range(len(contents)):
        if what_to_grep in contents[i]:
            to_return.append(contents[i].split()[item_to_get])
    
    return to_return
