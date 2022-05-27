import subprocess
import os

path  = "C:\\Users\\KChan\\source\\repos\\MPM-Snow-FORK\\renders\\frames\\"

all_files = []

for f in os.listdir():
    if (f[-4:] == ".xml"):
        all_files.append(f)

print(all_files)

for file in all_files:
    subprocess.call(["mitsuba", file])

    
        

    