#%% Librarys
import numpy as np
from pathlib import Path
from string import ascii_letters
#%% Funções
PATH_WRITE = Path().absolute().parent / 'observers'
if not PATH_WRITE.exists():
    PATH_WRITE.mkdir()

def probes(number_of_probes:int = 30, lim:tuple = (2 , 102)):
    p = np.linspace(lim[0], lim[1], number_of_probes)
    with open(PATH_WRITE/'probes.txt', 'w') as file:
        file.write('probeLocations\n\t(\n')
        for i in p:
            file.write(f'\t\t({round(i,3)} 0 0)\n')
        file.write('\t);')

def microphones(number_of_observer:int = 30, lim:tuple = (2,102), 
                lenght:float=1):
    m = np.linspace(lim[0], lim[1], number_of_observer)
    letter = list(ascii_letters)
    with open(PATH_WRITE/'microphones.txt','w') as file:
        file.write('observers\n{\n')
        for k,i in enumerate(m):
            file.write(f'\tR-{letter[k]}\n' + '\t{\n')
            file.write(f'\t\tposition\t({round(i,3)} 0 {lenght/2});\n')
            file.write('\t\tpRef\t2.0e-5;\n')
            file.write('\t\tfftFreq\t1024;\n\t}\n')
        file.write('}')

#%% Run
probes()
microphones()
# %%
