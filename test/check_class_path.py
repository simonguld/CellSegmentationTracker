import os, sys

class_path = os.path.join(os.environ['USERPROFILE'], 'Fiji.app' + os.sep + 'jars')

def main():
    print(class_path)

    # Extract all folders from class_path
    folders = os.listdir(class_path)

    #print(folders)

    # Load all jars from class_path.txt separated by space:
    jars = [jar for jar in open('class_path.txt').read().split(' ')]
    
    # Find the jars that are not in class path
    missing_jars = [jar for jar in jars if jar not in folders]

    print(missing_jars)

if __name__ == '__main__':
    main()