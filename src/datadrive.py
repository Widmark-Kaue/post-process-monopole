#%%
from pathlib import Path

PATH_DATA   = Path().cwd().with_name('data')
PATH_TEST  = Path().cwd().with_name('test')
PATH_WRITE  = Path().cwd().with_name('observers')
PATH_DATA.mkdir(exist_ok=True)
PATH_TEST.mkdir(exist_ok=True)
PATH_WRITE.mkdir(exist_ok=True)
#%%
arq =PATH_TEST.joinpath('teste.txt').absolute()
print(arq.with_suffix('.dat'))
arq.replace(arq.with_suffix('.dat'))
print(arq)
#%%
for child in PATH_DATA.glob("**/p*"):
    
    if child.is_file():
        print(child)
        child.replace(child.with_suffix('.dat'))
        print(child)
 #%%