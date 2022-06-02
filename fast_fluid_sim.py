SIM_SIZE=256

import pygame
from ctypes import c_void_p, c_double, c_int, cdll
import numpy as np
from numpy.ctypeslib import ndpointer
from math import sqrt
from os import system
system("pwd > full_path.txt")
fp=open("full_path.txt", "r")
p=fp.read().replace("\n", "").replace("\r", "").replace(" ", "")+"/"
fp.close()
print("Trying to access ", p+"fluid_step.so")
lib=cdll.LoadLibrary(p+"fluid_step.so")
fluid_step=lib.fluid_step #fluid_step is the name of our C function
fluid_step.restype=ndpointer(dtype=c_double, ndim=3, shape=(SIM_SIZE,SIM_SIZE,5))
get_p=lib.get_p
get_p.restype=ndpointer(dtype=c_double, ndim=2, shape=(SIM_SIZE,))
get_meanT=lib.get_meanT
get_meanT.restype=c_double
get_meanD=lib.get_meanD
get_meanD.restype=c_double
get_sd_T=lib.get_sd_T
get_sd_T.restype=c_double
get_sd_D=lib.get_sd_D
get_sd_D.restype=c_double
magnetic=lib.magnetic
#magnetic.restype= #void


fluid_sim=[[[0,0,0, 0, 0] for i in range(SIM_SIZE)] for j in range(SIM_SIZE)] #each position in 2D has velocity vector, density, temperature  
fluid_sim=np.array(fluid_sim, dtype=c_double)
#fluid_sim.shape() #should be (SIM_SIZE, SIM_SIZE, 5)
fluid_sim[25][10][3]=3
fluid_sim[25][10][4]=10
t: c_double=0.01
k: c_double=0.08
fluid_sim=fluid_step(c_double(t), c_void_p(fluid_sim.ctypes.data), c_double(k), c_void_p(get_p().ctypes.data))



#used the simple pygame example code
pygame.init()
screen = pygame.display.set_mode([SIM_SIZE*2, SIM_SIZE*2])
running = True


#the fluid sim variable is fluid_sim
#the fluid sim function is fluid_step(t: float, array: list, k: float, p:list=[])
#fluid_sim is passed as argument 2, p is maintained for later (replaced by the second output of the function
# Fill the background with grey
screen.fill((55, 55, 55))
# Draw a solid blue circle in the center
#pygame.draw.circle(screen, (0, 0, 255), (250, 250), 75) 
# Flip the display (sends stuff to screen)

fluid_sim[25][20][3]=5
fluid_sim[25][20][4]=10
fluid_sim[200][100][3]=5
fluid_sim[200][100][4]=10
def coloring(pygame, fluid_sim):
	meanD: c_double=get_meanD()
	sd_d: c_double=get_sd_D(c_double(meanD))
	#sd_d=sqrt(sd_d) #already done in the C code
	meanT: c_double=get_meanT()
	sd: c_double=get_sd_T(c_double(meanT))
	#sd=sqrt(sd)
	for i in range(len(fluid_sim)):
		for j in range(len(fluid_sim)):
			green=int(min(255, max(0, fluid_sim[i][j][3]-meanD)/sd_d*10))
			if fluid_sim[i][j][4]>meanT:
				pygame.draw.circle(screen, (int(min(255, (fluid_sim[i][j][4]-meanT)/sd*10)), green, 0), (j*2, i*2), 1)
			elif fluid_sim[i][j][4]<meanT:
				pygame.draw.circle(screen, (0, green, (int(min(255, (meanT-fluid_sim[i][j][4])/sd*10)))), (j*2, i*2), 1)
	
prev_mouse=pygame.mouse.get_pos()
while running:
	# Did the user click the window close button?
	for event in pygame.event.get():
		if event.type == pygame.QUIT:
			running = False
	coloring(pygame, fluid_sim)
	new_mouse=pygame.mouse.get_pos()
	if (max(max(new_mouse), max(prev_mouse))<2*SIM_SIZE) and (min(min(new_mouse), min(prev_mouse))>0):
		xV=-prev_mouse[0]+new_mouse[0]
		yV=-prev_mouse[1]+new_mouse[1]
		fluid_sim[prev_mouse[1]//2][prev_mouse[0]//2][0]=xV
		fluid_sim[prev_mouse[1]//2][prev_mouse[0]//2][1]=yV
		fluid_sim[new_mouse[1]//2][new_mouse[0]//2][0]=xV
		fluid_sim[new_mouse[1]//2][new_mouse[0]//2][1]=yV
		fluid_sim[new_mouse[1]//2][new_mouse[0]//2][3]+=sqrt(xV**2+yV**2)
		fluid_sim[new_mouse[1]//2][new_mouse[0]//2][4]+=1
		#the magnetic function has the problem of wiping everything off the map
		#magnetic(c_double(prev_mouse[0]//2), c_double(prev_mouse[1]//2), c_double(0), c_double(new_mouse[0]//2), c_double(new_mouse[1]//2), c_double(0))
	prev_mouse=new_mouse
	fluid_sim=fluid_step(c_double(t), c_void_p(fluid_sim.ctypes.data), c_double(k), c_void_p(get_p().ctypes.data))
	pygame.display.flip()
# Done! Time to quit.
pygame.quit()

