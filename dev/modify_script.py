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
    path = "C:\\Users\\Simon Andersen\\Projects\\Projects\\CellSegmentationTracker\\test\\modify_script_dummy.py"
    txt_path = "C:\\Users\\Simon Andersen\\Projects\\Projects\\CellSegmentationTracker\\test\\params.txt"

    search_strings = ['target_height = 648', 'target_width = 698', 'target_height=648', 'target_width=698',  \
                      'target_width, target_height = 698, 648']
    replace_strings = ['target_height = target_height', 'target_width = target_width', 'target_height=target_height', 'target_width=target_width', \
                       'target_height, target_width = target_height, target_width' ]
    
    # STEP 0: gen file to hold target values
    np.savetxt(txt_path, np.array([410, 311]), fmt='%i')

    # STEP 0.5: Read target values from file and convert to integers
    a,b = np.loadtxt(txt_path, dtype=int)
    print(a,b)

    dir = os.path.dirname(os.path.abspath(__file__))
    txt_path = os.path.join(dir, 'params.txt')

    # STEP 1: Preprend text to file
    msg = f"\nimport numpy as np\nimport os\ndir = os.path.dirname(os.path.abspath(__file__))\ntarget_height, target_width = np.loadtxt(os.path.join(dir, 'params.txt')) \n"  
    prepend_text_to_file(path, msg)

    # STEP 2: Search and modify file
    for search_string, replace_string in zip(search_strings, replace_strings):
        search_and_modify_file(path, search_string, replace_string)


if __name__ == '__main__':
    main()
