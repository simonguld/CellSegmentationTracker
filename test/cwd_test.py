import time
import sys, os


path1 = os.getcwd()
path2 = os.path.dirname(os.path.abspath(__file__))
print(os.path.abspath(os.path.join(path2, os.pardir)))
path3 = "C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\python.exe"
print(path1)
print(path2)
print(print(os.path.abspath(os.path.join(path3, os.pardir))))
print(os.path.relpath(path3,os.pardir))
executable = list(os.path.split(path3))[-1]

path4 = "C:\\Users\\Simon Andersen\\Fiji.app\\ImageJ-win64.exe"
suffix2 = list(os.path.split(path4))[-1]
print(suffix2)
path5  = "C:\\Users\\Simon Andersen\\Fiji.app\\ImageJ-win64.tif"



print("hej med dig")
try:
    time.sleep(8)
except: 
    KeyboardInterrupt
    sys.exit(0)

print("jeg hedder kaj")