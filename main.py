import yaml, os, sys, re
import shutil as shu
import module_input
import rand_structure
from premelt import do_premelt
from conv_test import *
from module_Tm import *
from module_vasp import *
from melt_quench import do_mq

# Open yaml file and update
src_dir = os.path.dirname(os.path.abspath(__file__))
inp_yaml = module_input.inp_load(src_dir, sys.argv[1])
inp_yaml['vasp_config']['num_tasks'] = str(sys.argv[2])
print(inp_yaml)
# Setting directories
os.makedirs('Output', exist_ok=True)
os.chdir('Output')
PWD=os.getcwd()
os.makedirs(inp_yaml['composition'], exist_ok=True)
os.chdir(inp_yaml['composition'])
module_input.write_inp(inp_yaml)
mat_dir=os.getcwd()
os.makedirs("Inputs", exist_ok=True)

log = module_input.check_log(inp_yaml)

# Premelt
if inp_yaml['Actions']['premelt']:
    os.makedirs('premelt', exist_ok=True)
    os.chdir('premelt')

    # Random spraying
    rand_structure.create_rand_structure(composition=inp_yaml['composition'], pot_dir=inp_yaml['pot_dir'],\
            target_num=inp_yaml['target_num']) #TO DO check min distance, unary or quaternary support
    shu.copy("POTCAR", "../Inputs/")
    do_premelt(src_dir=src_dir, vasp_config=inp_yaml['vasp_config'],\
            premelt_config=inp_yaml['premelt_config'], mat_name=inp_yaml['composition'])
    os.chdir('../')
    with open(log, 'a') as s:
        s.write('Premelting is done\n')

if inp_yaml['Actions']['convergence_test']:
    os.makedirs('conv_test', exist_ok=True)
    os.chdir('conv_test')
    do_conv_test(mat_dir=mat_dir, config=inp_yaml['conv_test_config'], vasp_config=inp_yaml['vasp_config'], log_dir=log)
    os.chdir("../")
    with open(log, 'a') as s:
        s.write('Convergence test is done\n')

# Check if KP is gamma or not
vasp_ver = 'std'
with open(log, 'r') as f:
    lines = f.readlines()
    for line in lines:
        if 'Converged KP:' in line:
            if line.strip().split()[2] == 'G':
                vasp_ver = 'gam'

if inp_yaml['Actions']['find_Tm']:
    os.makedirs('find_Tm', exist_ok=True)
    os.chdir('find_Tm')
    initial_config={"POS_dir":mat_dir+'/premelt/CONTCAR', "INCAR_dir":mat_dir+'/Inputs/INCAR',\
            "POT_dir":mat_dir+'/premelt/POTCAR', "KP_dir":mat_dir+'/Inputs/KPOINTS'}

    initial_rlx(initial_config=initial_config, vasp_config=inp_yaml['vasp_config'], vasp_version=vasp_ver, log_dir=log)
    
    continue_dir='/'.join( [os.getcwd(), 'Step0_relax', str(len(os.listdir('Step0_relax'))-1)])
    # initial_config={"POS_dir":continue_dir+'/CONTCAR',"INCAR_dir":continue_dir+'/INCAR',
    # "POT_dir":continue_dir+'/POTCAR',"KP_dir":continue_dir+'/KPOINTS'}
    find_Tm(continue_dir=continue_dir, input_Tm=inp_yaml['find_Tm_config'], vasp_config=inp_yaml['vasp_config'], vasp_version=vasp_ver, log_dir=log)
    os.chdir("../")
    with open(log, 'a') as s:
        s.write('MD setting is done\n')

p = re.compile('\d+')
def check_last():
    step = 0
    last_idx = 0
    for n in os.listdir('.'):
        if 'OUTCAR' in n:
            if n != 'OUTCAR':
                idx = int(p.findall(n)[0])
                if idx > last_idx:
                    last_idx = idx
            with open(n, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if 'free  ' in line:
                        step += 1
    if 'OUTCAR' in os.listdir('.'):
        last_idx += 1
        os.rename('OUTCAR', 'OUTCAR%s'%last_idx)
        os.rename('XDATCAR', 'XDATCAR%s'%last_idx)
        os.rename('POSCAR', 'POSCAR%s'%last_idx)
        os.rename('CONTCAR', 'POSCAR')

    return step, last_idx

if inp_yaml['Actions']['melt']:
    os.makedirs("melt", exist_ok=True)
    os.chdir("melt")

    step, last_idx = check_last()
    if step == 0:
        copy_inputs("../Inputs", './', ["POSCAR_to_melt","POTCAR",'KPOINTS',"INCAR"])
        os.rename("POSCAR_to_melt", "POSCAR")
        edit_INCAR("./INCAR", {"NSW": inp_yaml['melt_config']['steps']})
    else:
        new_nsw = int(inp_yaml['melt_config']['steps']) - step
        edit_INCAR("./INCAR", {"NSW": f'{new_nsw}'})

    run_vasp(inp_yaml['vasp_config']['mpicommand'], inp_yaml['vasp_config']['num_tasks'], inp_yaml['vasp_config'][vasp_ver])
    os.chdir("../")
    with open(log, 'a') as s:
        s.write('Melting is done\n')

if inp_yaml['Actions']['quench']:
    os.makedirs("quench", exist_ok=True)
    os.chdir("quench")
    step, last_idx = check_last()
    if step == 0:
        copy_inputs("../melt",'./', ["CONTCAR","POTCAR",'KPOINTS',"INCAR"])
        os.rename("CONTCAR", "POSCAR")
        nsw=int((float(grep_nth_item('TEBEG','./INCAR',2)[0]) - inp_yaml['quench_config']['Temp_end'])/inp_yaml['quench_config']['quenching_rate']*1000/2)
        edit_INCAR("./INCAR",{"NSW":f"{nsw}", "TEEND":inp_yaml['quench_config']['Temp_end']})
    else:
        nsw=int((float(grep_nth_item('TEBEG','./INCAR',2)[0]) - inp_yaml['quench_config']['Temp_end'])/inp_yaml['quench_config']['quenching_rate']*1000/2)
        new_nsw = nsw-step
        new_tebeg = 300 + new_nsw/nsw*int((float(grep_nth_item('TEBEG','../melt/INCAR',2)[0]) - inp_yaml['quench_config']['Temp_end']))
        edit_INCAR("./INCAR",{"NSW":f"{new_nsw}", "TEBEG": f'{new_tebeg}'})

    run_vasp(inp_yaml['vasp_config']['mpicommand'], inp_yaml['vasp_config']['num_tasks'], inp_yaml['vasp_config'][vasp_ver])
    os.chdir("../")
    with open(log, 'a') as s:
        s.write('Quenching is done\n')

if inp_yaml['Actions']['anneal']:
    os.makedirs("anneal", exist_ok=True)
    os.chdir("anneal")
    step, last_idx = check_last()
    if step == 0:
        copy_inputs("../quench",'./', ["CONTCAR","POTCAR",'KPOINTS',"INCAR"])
        os.rename("CONTCAR", "POSCAR")
        edit_INCAR("./INCAR", {"NSW":inp_yaml['annealing_config']['steps'], "TEBEG":inp_yaml['annealing_config']['Temp_start'],\
                "TEEND":inp_yaml['annealing_config']['Temp_end']})
    else:
        new_nsw = int(inp_yaml['annealing_config']['steps']) - step
        edit_INCAR("./INCAR", {"NSW": f'{new_nsw}'})

    run_vasp(inp_yaml['vasp_config']['mpicommand'], inp_yaml['vasp_config']['num_tasks'], inp_yaml['vasp_config'][vasp_ver])
    os.chdir("../")
    with open(log, 'a') as s:
        s.write('Annealing is done\n')

with open(log, 'a') as s:
    s.write('Calculation done\n')
os.chdir(PWD)
