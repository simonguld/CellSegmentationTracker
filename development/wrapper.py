import os, sys
import platform

## run jython_trial from this script, while keeping the username path variable

def get_trackmate_path():
    if 0:
        if platform.system() == 'Windows':
            bracket = '\\'
        else:
            bracket = '/'
    bracket = '\\'
    user_path = os.environ['USERPROFILE']
    class_path_arg = os.path.join(user_path, 'Fiji.app' + bracket + 'jars' + bracket + '*')
    return class_path_arg


def main():
    script_name = sys.argv[-1]
    bracket = '\\'
    trackmate_path = get_trackmate_path()
    user_path = os.environ['USERPROFILE']

    path_from_fiji = os.path.join(user_path, 'Projects' + bracket + 'Projects' + bracket + 'CellSegmentationTracking\\jython_trial.py')

    #path_from_fiji = 'C:\\Users\\Simon Andersen\\Projects\\Projects\\CellSegmentationTracking\\jython_trial.py'
    #path_from_fiji = 'jython_trial.py'
    # run jython script
    #os.system(f'jython -J-cp "{trackmate_path}" {script_name}')
    #trackmate_path = 'C:\\Users\\Simon\n Andersen\\Fiji.app\\jars\\*'
    path_final = 'jython -J-cp ' + '"' +trackmate_path + '"'+ " " + '"' +path_from_fiji +'"'
    #print(path_final)
    os.system(path_final)
    print(trackmate_path, "\n")
    print(path_from_fiji)

    #os.system(f'jython -J-cp "{trackmate_path}" "{path_from_fiji}"')

if __name__ == '__main__':
    main()
