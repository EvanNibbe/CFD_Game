# CFD_Game
This pygame has a lot of lag due to being written in Python, and uses my own variation on Navier-Stokes
To run the game:

Install Python3 from https://www.python.org/ftp/python/3.10.4/python-3.10.4-amd64.exe

In the cmd prompt:

python3-pip install pygame #if that doesn't work, just do "pip install pygame", if that doesn't work, download MSYS2, download scoop from scoop.sh and install it the way it is written on Windows Powershell, then install Python with scoop install python3, then pip install pygame

python3 fluid_sim.py



#####
To make the game faster, use the "fast" or "good" folders with the following to compile the C code for your system (please upload the compiled code in a pull request so I can include it here):
cc -fPIC -shared -o fluid_step.so fluid_step.c
Then use:
python3 fast_fluid_sim.py
