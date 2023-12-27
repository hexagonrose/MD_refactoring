import yaml, collections, os

def inp_update(inp0, inp):    
    for key in list(inp.keys()):
        if isinstance(inp0, collections.Mapping):
            if isinstance(inp[key], collections.Mapping) and inp[key]:
                returned = inp_update(inp0.get(key, {}), inp[key])
                inp0[key] = returned
            else:
                inp0[key] = inp[key]
        else:
            inp0 = {key: inp[key]}

    return inp0

def convert(inp0, key):
    check = inp0[key].upper()
    if check[0] == '.':
        if check[1] == 'T':
            check = True
        elif check[1] == 'F':
            check = False
    elif check[0] == 'T':
        check = True
    elif check[0] == 'F':
        check = False
    
    if isinstance(check, bool):
        inp0[key] = check

def inp_load(inp0_dir, inp_dir):
    # Open default yaml and update yaml
    with open(inp0_dir+'/configure_default.yaml', 'r') as O:
        inp0=yaml.safe_load(O)
    with open(inp_dir, 'r') as O:
        inp=yaml.safe_load(O)
    inp_update(inp0, inp)

    # Convert boolean strings to boolean (ex. 'T' -> True, '.F.' -> False)
    #for key in inp0.keys():
        #if not isinstance(inp0[key], bool) and isinstance(inp0[key], str):
            #convert(inp0, key)

    return inp0

def write_inp(inp_yaml):
    with open('./final_configure.yaml', 'w') as O:
        yaml.safe_dump(inp_yaml, O, default_flow_style=False)

def check_log(inp_yaml):
    if 'LOG' in os.listdir('./'):
        with open('LOG', 'r') as f:
            lines = f.readlines()
            for line in lines:
                if 'Premelting is done' in line:
                    inp_yaml['Actions']['premelt'] = False
                elif 'Convergence test is done' in line:
                    inp_yaml['Actions']['convergence_test'] = False
                elif 'MD setting is done' in line:
                    inp_yaml['Actions']['find_Tm'] = False
                elif 'Melting is done' in line:
                    inp_yaml['Actions']['melt'] = False
                elif 'Quenching is done' in line:
                    inp_yaml['Actions']['quench'] = False
                elif 'Annealing is done' in line:
                    inp_yaml['Actions']['anneal'] = False
    else:
        with open('LOG', 'a') as s:
            s.write('Calculation start\n')

    log = os.getcwd()+'/LOG' 

    return log