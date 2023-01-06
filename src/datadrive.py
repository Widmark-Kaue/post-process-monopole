#%%
from pathlib import Path

PATH_DATA   = Path('data')
PATH_WRITE  = Path('observers')
PATH_DATA.mkdir(exist_ok=True)
PATH_WRITE.mkdir(exist_ok=True)

arq = PATH_DATA.joinpath('teste.txt').absolute()
print(arq.with_suffix('.dat'))
arq.replace(arq.with_suffix('.dat'))
print(arq)
# for child in PATH_DATA.glob("**/p*"):
    
#     if child.is_file():
#         child.rename('p.dat')
#         print(child)
 #%%