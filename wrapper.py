import os, sys
import platform

## run jython_trial from this script, while keeping the username path variable

def get_trackmate_path():
    if platform.system() == 'Windows':
        bracket = '\\'
    else:
        bracket = '/'
        
    user_path = os.environ['USERPROFILE']
    class_path_arg = os.path.join(user_path, 'Fiji.app' + bracket + 'jars' + bracket + '*')
    return class_path_arg


def main():
    script_name = sys.argv[-1]

    trackmate_path = get_trackmate_path()

    # run jython script
    os.system(f'jython -J-cp "{trackmate_path}" {script_name}')

if __name__ == '__main__':
    main()
