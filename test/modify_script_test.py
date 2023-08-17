import os
import sys
import numpy as np

def write_to_file(file_path, text):
    """
    Write to file.
    """
    with open(file_path, 'w') as file:
        file.write(text)
        file.close()
    return

def read_file(file_path):
    """
    Read file.
    """
    with open(file_path, 'r') as file:
        filedata = file.read()
        file.close()
    return filedata

def prepend_text_to_file(file_path, text):
    """
    Prepend text to a file.
    """
    with open(file_path, 'r+') as file:
        content = file.read()
        file.seek(0, 0)
        file.write(text.rstrip('\r\n') + '\n' + content)
        file.close()
    return

def search_and_modify_file(file_path, search_string, modify_string):
    """
    Search for a string in a file and replace it with another string.
    """
    with open(file_path, 'r') as file:
        filedata = file.read()
    filedata = filedata.replace(search_string, modify_string)
    with open(file_path, 'w') as file:
        file.write(filedata)
        file.close()
    return

def main():
    cellpose_python_filepath = 'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\python.exe'

    # STEP1: Generate path name
    cp_dir = os.path.dirname(cellpose_python_filepath)
    dir = os.path.join(cp_dir, 'Lib', 'site-packages', 'cp_copy')

    model_parent_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
    model_path = os.path.join(model_parent_dir, 'models', 'models.py')
    model_dir = os.path.dirname(model_path)
    paths = [model_path, os.path.join(dir, 'dynamics.py'), os.path.join(dir, 'cli.py')]
    txt_paths = [os.path.join(model_dir, 'params.txt'), os.path.join(dir,'params.txt'),os.path.join(dir,'params.txt'),os.path.join(dir,'params.txt')]
    dirs = [model_dir, dir, dir, dir]
    search_strings = ["'--flow_threshold', default=0.4", "'--cellprob_threshold', default=0", \
                      'flow_threshold=0.4', 'cellprob_threshold=0.0',  \
                      'flow_threshold= 0.4', 'flow_threshold = 0.4', \
                        'cellprob_threshold= 0.0', 'cellprob_threshold = 0.0' ]
    replace_strings = ["'--flow_threshold', default=flow_param", "'--cellprob_threshold', default=cellprop_param", \
                          'flow_threshold=flow_param', 'cellprob_threshold=cellprop_param', \
                            'flow_threshold= flow_param', 'flow_threshold = flow_param', \
                            'cellprob_threshold= cellprop_param', 'cellprob_threshold = cellprop_param' ]


    flow_param, cellprop_param = 0.3, -0.5
    # STEP 0: gen file to hold target values (both places)
    


    #TRY ON JUST ONE FILE
    # STEP 1: Preprend text to file
    msg = f"\nimport numpy as np\nimport os\ndir = os.path.dirname(os.path.abspath(__file__))\nflow_param, cellprop_param = np.loadtxt(os.path.join(dir, 'params.txt')) \n\n"  
    #prepend_text_to_file(model_path, msg)

    # STEP 2: Search and modify file
    # First do it for model_path
    if not os.path.isfile(os.path.join(model_dir,'params.txt')):
        np.savetxt(os.path.join(model_dir,'params.txt'), np.array([flow_param, cellprop_param]))
        prepend_text_to_file(model_path, msg)
        for search_string, replace_string in zip(search_strings, replace_strings):
            search_and_modify_file(model_path, search_string, replace_string)

    # Then do it for the rest

    if not os.path.isfile(os.path.join(dir,'params.txt')):
        np.savetxt(os.path.join(dir,'params.txt'), np.array([flow_param, cellprop_param]))
        for path in paths:
            prepend_text_to_file(path, msg)
            for search_string, replace_string in zip(search_strings, replace_strings):
                search_and_modify_file(path, search_string, replace_string)



if __name__ == '__main__':
    main()
