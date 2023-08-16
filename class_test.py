class CellSegmentationTracker:
    def __init__(self, path1, path2):
        self.path1 = path1
        self._path2 = path2
        
    def print_path1(self):
        print("Path 1:", self.path1)
        
    def print_path2(self):
        print("Path 2:", self._path2)



def main():
    cst = CellSegmentationTracker('x', "y")
    from subprocess import Popen, PIPE
    import os, sys

    print(cst.path1)
    print(cst._path2)
    print(os.getcwd())

    p1 = Popen(['python','cwd_test.py'], stdout=PIPE)
    p1.wait()
    text = p1.communicate()[0]
    print(text)
    
  #  os.system('python ./cwd_test.py')
    cst.print_path2()
    

if __name__ == "__main__":  
    main()

