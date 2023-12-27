import subprocess as sub
import os
import shutil as shu
import multiprocessing
from module_vasp import *

def do_premelt(src_dir, vasp_config, premelt_config, mat_name):
    #TODO copy INCAR & POSCAR from user's input
    shu.copy(src_dir+"/INCAR_premelt", "./INCAR")
    edit_INCAR("./INCAR", {"SYSTEM":mat_name})
    shu.copy(src_dir+"/KPOINTS", "./KPOINTS")
    #TODO if potcar is not here, check and copy POTCAR
    run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config['gam'])
    #TODO check if MD is good or not

