import os
import sys

def makedir(path):
    try:
        if not os.path.exists(path):
            os.makedirs(path)
    except OSError:
        print('Have no permission to build a directory')
        sys.exit()

