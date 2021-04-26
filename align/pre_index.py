import os

def build(tool,ref,prefix=None,index=False):
    index_script = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'index.sh')
    statistics_script = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'statistics.sh')
    if prefix and index:
        cmd = (f"sh {index_script} {ref} {tool} {prefix} "
               f"sh {statistics_script} {ref} {prefix}".format(**locals()))
    else:
        cmd = ''
    return cmd