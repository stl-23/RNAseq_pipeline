import os
import sys
import subprocess
import logging
import configparser

logging.basicConfig(
    format='[%(asctime)s - %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)

def makedir(path):
    try:
        if not os.path.exists(path):
            os.makedirs(path)
    except OSError:
        print('Have no permission to build a directory')
        sys.exit()

def getConfig(section, key):
    config = configparser.ConfigParser()
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'softwares.config')
    config.read(path)
    return config.get(section, key)

def run_shell_cmd(cmd):
    run = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    if run.returncode == 0:
        print("STDOUT:", run.stdout)
    else:
        print("STDOUT:", run.stdout)
        print("STDERR:", run.stderr)


