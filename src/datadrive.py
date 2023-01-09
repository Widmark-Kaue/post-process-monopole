#%%
from pathlib import Path

PATH_DATA   = Path().cwd().joinpath('data')
PATH_TEST   = Path().cwd().joinpath('test')
PATH_WRITE  = Path().cwd().joinpath('observers')

PATH_DATA   .mkdir(exist_ok=True)
PATH_TEST   .mkdir(exist_ok=True)
PATH_WRITE  .mkdir(exist_ok=True)
