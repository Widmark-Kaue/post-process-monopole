from pathlib import Path

# from google.colab import drive

PATH_DATA = Path('data')

PATH_WRITE = Path('observers')
if not PATH_WRITE.exists():
    PATH_WRITE.mkdir()
