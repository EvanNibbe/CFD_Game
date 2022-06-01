import pygame
SIM_SIZE=256
fluid_sim=[[[{'x':0, 'y':0, 'z':0}, 0, 0] for i in range(SIM_SIZE)] for j in range(SIM_SIZE)] #each position in 2D has velocity vector, density, temperature  
from math import sqrt, sin, cos, tan, atan, atan2, log, exp 
from secrets import randbelow as rand
def fluid_step(t: float, array: list, k: float, p:list): #k is the positive inverse viscosity  
	#first, need to create the new array from diffusion 
	#and advection from array 
	#second, I need to make the new array go through the 
	#process of finding curl-free field of that new array 
	#third, I need to subtract the curl-free field from new array 
	#finally, return new array 
	def move(arg: list, dest: list, weight: float): 
		#print("Passed to move: arg", arg, "\t dest:", dest, "\t weight", weight)
		#orig_v=sqrt(dest[0]["x"]**2+dest[0]["y"]**2+dest[0]["z"]**2)
		#new_v=sqrt(arg[0]["x"]**2+arg[0]["y"]**2+arg[0]["z"]**2)
		
		orig_m=dest[1] 
		dest[1]+=arg[1]*weight 
		orig_t=dest[2] 
		if abs(orig_m)<.00001: 
			dest[0]["x"]+=arg[0]["x"]*weight 
			dest[0]["y"]+=arg[0]["y"]*weight 
			dest[0]["z"]+=arg[0]["z"]*weight 
			dest[2]+=arg[2]*weight 
		elif abs(orig_m+arg[1])>0.001: 
			dest[0]["x"]=(orig_m*dest[0]["x"]+arg[1]*arg[0]["x"]*weight)/(orig_m+arg[1]*weight) 
			dest[0]["y"]=(orig_m*dest[0]["y"]+arg[1]*arg[0]["y"]*weight)/(orig_m+arg[1]*weight)
			dest[0]["z"]=(orig_m*dest[0]["z"]+arg[1]*arg[0]["z"]*weight)/(orig_m+arg[1]*weight)
			dest[2]=(orig_m*orig_t+arg[1]*arg[2]*weight)/(orig_m+arg[1]*weight) 
		else: 
			dest[2]=(orig_m*orig_t+arg[1]*arg[2]*weight)/max(max(.001, max(orig_m, arg[1]), orig_m+arg[1]*weight), min(.01, weight)) 
		return dest 
	def diffusion(t: float, array: list, k: float): 
	#note: one problem with the naive approach in move (the function above) is that it doesn't 
	#consider the effect of thermal mass, that the more stuff you have, the more its temperature 
	#affects the stuff around it (as compared to the less massive stuff) 
	#e.g. a closed system with a candle placed on an iceberg results in the candle losing its 
	#high temperature, and the iceberg not noticably changing temperature. 
		n_arr=[[[{'x':0, 'y':0, 'z':0}, 0, 0] for i in range(len(array))] for j in range(len(array[0]))] 
		for i in range(len(array)): 
			for j in range(len(array[0])): 
				avg=[{'x':0, 'y':0, 'z':0}, 0, 0] 
				div=4.0 
				if i>0 and i<len(array)-1: 
					avg=move(array[i+1][j], avg, 1) 
					avg=move(array[i-1][j], avg, 1) 
				elif i==0: 
					div-=1 
					avg=move(array[i+1][j], avg, 1) 
				else: 
					div-=1 
					avg=move(array[i-1][j], avg, 1) 
				if j>0 and j<len(array[0])-1: 
					avg=move(array[i][j+1], avg, 1) 
					avg=move(array[i][j-1], avg, 1) 
				elif j==0: 
					div-=1 
					avg=move(array[i][j+1], avg, 1) 
				else: 
					div-=1 
					avg=move(array[i][j-1], avg, 1) 
				avg[0]["x"]/=div 
				avg[0]["y"]/=div 
				avg[0]["z"]/=div 
				avg[1]/=div 
				avg[2]/=div 
				n_arr[i][j]=move(avg, n_arr[i][j], 1/(1+k)) 
		return n_arr 
	def advection(t: float, array: list, k: float): 
		n_arr=[[[{'x':0, 'y':0, 'z':0}, 0, 0] for i in range(len(array))] for j in range(len(array[0]))] 
		for i in range(len(array)): 
			for j in range(len(array[0])): 
				avg=[{'x':0, 'y':0, 'z':0}, 0, 0] 
				div=4.0 
				frm=(j-array[i][j][0]["x"]*t, i-array[i][j][0]["y"]*t) 
				dxA=frm[0]-int(frm[0]) 
				dyA=frm[1]-int(frm[1]) 
				dA=sqrt(dxA**2+dyA**2) 
				dxB=int(frm[0]+1)-frm[0] 
				dyB=int(frm[1])-frm[1] 
				dB=sqrt(dxB**2+dyB**2) 
				dxC=int(frm[0])-frm[0] 
				dyC=int(frm[1]+1)-frm[1] 
				dC=sqrt(dxC**2+dyC**2) 
				dxD=int(frm[0]+1)-frm[0] 
				dyD=int(frm[1]+1)-frm[1] 
				dD=sqrt(dxD**2+dyD**2) 
				tot=dA+dB+dC+dD 
				if frm[0]>0 and frm[0]<len(array[0])-1: 
					if frm[1]>0 and frm[1]<len(array)-1: 
						avg=move(array[int(frm[1])][int(frm[0])], avg, (dB+dC+dD)/(3.0*tot)) 
						avg=move(array[int(frm[1])+1][int(frm[0])], avg, (dA+dC+dD)/(3.0*tot)) 
						avg=move(array[int(frm[1])][int(frm[0])+1], avg, (dA+dB+dD)/(3.0*tot)) 
						avg=move(array[int(frm[1])+1][int(frm[0])+1], avg, (dA+dB+dC)/(3.0*tot)) 
				n_arr[i][j]=move(avg, n_arr[i][j], 1) 
		return n_arr 
	def curl_free(array: list, p:list): #p is preferrably the p scalar matrix from the previous iteration 
		if len(p)<len(array): 
			p=[[0 for i in range(len(array))] for j in range(len(array[0]))] 
		n_arr=[[[{'x':0, 'y':0, 'z':0}, 0, 0] for i in range(len(array))] for j in range(len(array[0]))] 
		curl_v=[[0 for i in range(len(array))] for j in range(len(array[0]))] 
		for i in range(len(array)): 
			for j in range(len(array)): 
				div=2.0 
				avg=0.0 
				if i>0 and i<len(array)-1 and j>0 and j<len(array)-1: 
					avg+=array[i+1][j][0]["y"]-array[i-1][j][0]["y"] 
					avg+=array[i][j+1][0]["x"]-array[i-1][j][0]["x"] 
				avg/=div 
				curl_v[i][j]=avg 
		for k in range(3): 
			prev_p=[[i for i in j] for j in p] 
			for i in range(len(array)): 
				for j in range(len(array)): 
					avg=0 
					div=4.0 
					if i>0 and i<len(array)-1 and j>0 and j<len(array)-1: 
						avg+=prev_p[i+1][j] 
						avg+=prev_p[i-1][j] 
						avg+=prev_p[i][j-1]
						avg+=prev_p[i][j+1] 
					avg-=curl_v[i][j] 
					avg/=div 
					p[i][j]=avg 
					if k==2 and i>0 and i<len(array)-1 and j>0 and j<len(array)-1: 
						div_p=[(p[i][j+1]-p[i][j-1])/2, (p[i+1][j]-p[i-1][j])/2] 
						n_arr[i][j][0]["x"]=array[i][j][0]["x"]-div_p[0] 
						n_arr[i][j][0]["y"]=array[i][j][0]["y"]-div_p[1] 
		return (n_arr, p) 
	def temperate(array: list, k: float): #use temperature and pressure vs inverse viscosity k 
		#then use density*temperature compared to adjacent spot to determine outward flow, temperature may be negative 
		#if there is outward flow, then the adjacent spots affected get their velocity changed by a small amount according to  
		#the gradient, and the positive temperature randomizes the resulting velocity's angle by that many degrees. 
		#a small k acts as a moderating force. 
		n_arr=[[move(i, [{"x":0, "y":0, "z":0}, 0, 0], 1) for i in j] for j in array] 
		delta_v=[[[0, 0] for i in j] for j in array] 
		for i in range(len(array)): 
			for j in range(len(array[0])): 
				if i>0 and i<len(array)-1 and j>0 and j<len(array)-1: 
					#we have to alter the velocities at the 8 adjacent spots 
					dM=array[i][j][1] 
					tM=array[i][j][2] 
					Sur=[array[i-1][j-1], array[i-1][j], array[i-1][j+1], array[i][j-1], array[i][j+1], array[i+1][j-1], array[i+1][j], array[i+1][j+1]] 
					if (dM*tM>Sur[0][1]*Sur[0][2]): 
						coef=min(2*k, k*(dM*tM/max(.1, Sur[0][1]*Sur[0][2])-1))/sqrt(2) #sqrt 2 is distance between corner to middle 
						delta_v[i-1][j-1][0]+=coef*(-1) #-1 is the direction along x-axis 
						delta_v[i-1][j-1][1]+=coef*(-1) 
					if (dM*tM>Sur[1][1]*Sur[1][2]): 
						coef=min(2*k, k*(dM*tM/max(.1, Sur[1][1]*Sur[1][2])-1)) 
						delta_v[i-1][j][0]+=0 
						delta_v[i-1][j][1]+=coef*(-1) 
					if (dM*tM>Sur[2][1]*Sur[2][2]): 
						coef=min(2*k, k*(dM*tM/max(.1, Sur[2][1]*Sur[2][2])-1))/sqrt(2) 
						delta_v[i-1][j+1][0]+=coef*(1) 
						delta_v[i-1][j+1][1]+=coef*(-1) 
					if (dM*tM>Sur[3][1]*Sur[3][2]): 
						coef=min(2*k, k*(dM*tM/max(.1, Sur[3][1]*Sur[3][2])-1) )
						delta_v[i][j-1][0]+=coef*(-1) 
						delta_v[i][j-1][1]+=0 
					if (dM*tM>Sur[4][1]*Sur[4][2]): 
						coef=min(2*k, k*(dM*tM/max(.1, Sur[4][1]*Sur[4][2])-1) )
						delta_v[i][j+1][0]+=coef*(1) 
						delta_v[i][j+1][1]+=0 
					if (dM*tM>Sur[5][1]*Sur[5][2]): 
						coef=min(2*k, k*(dM*tM/max(.1, Sur[5][1]*Sur[5][2])-1))/sqrt(2) 
						delta_v[i+1][j-1][0]+=coef*(-1) 
						delta_v[i+1][j-1][1]+=coef*(1) 
					if (dM*tM>Sur[6][1]*Sur[6][2]): 
						coef=min(2*k, k*(dM*tM/max(.1, Sur[6][1]*Sur[6][2])-1) )
						delta_v[i+1][j][0]+=0 
						delta_v[i+1][j][1]+=coef*(1) 
					if (dM*tM>Sur[7][1]*Sur[7][2]): 
						coef=min(2*k, k*(dM*tM/max(.1, Sur[7][1]*Sur[7][2])-1))/sqrt(2) 
						delta_v[i+1][j+1][0]+=coef*(1) 
						delta_v[i+1][j+1][1]+=coef*(1) 
		#now we slightly randomize velocities in delta_v according to temperature 
		#we are changing the direction, not the magnitude. 
		for i in range(len(array)): 
			for j in range(len(array[0])): 
				mag=sqrt(delta_v[i][j][0]**2+delta_v[i][j][1]**2) 
				if mag>0.001: 
					delta_v[i][j][0]+=(2.0-rand(1000000)/1000000.0)*k*array[i][j][2] #multiply by temperature 
					delta_v[i][j][1]+=(2.0-rand(1000000)/1000000.0)*k*array[i][j][2] 
					n_dist=sqrt(delta_v[i][j][0]**2+delta_v[i][j][1]**2) 
					if n_dist<.001: 
						delta_v[i][j][0]=mag 
					else: 
						delta_v[i][j]=[delta_v[i][j][0]*mag/n_dist, delta_v[i][j][1]*mag/n_dist] 
		p=[[0 for i in range(len(array))] for j in range(len(array[0]))] 
		curl_v=[[0 for i in range(len(array))] for j in range(len(array[0]))] 
		for i in range(len(array)): 
			for j in range(len(array)): 
				div=2.0 
				avg=0.0 
				if i>0 and i<len(array)-1 and j>0 and j<len(array)-1: 
					avg+=delta_v[i+1][j][1]-delta_v[i-1][j][1] 
					avg+=delta_v[i][j+1][0]-delta_v[i-1][j][0] 
				avg/=div 
				curl_v[i][j]=avg 
		new_p=[[[i for i in j] for j in p] for k in range(5)] 
		for k in range(4): 
			for i in range(len(array)): 
				for j in range(len(array)): 
					avg=0 
					div=4.0 
					if i>0 and i<len(array)-1 and j>0 and j<len(array)-1: 
						avg+=new_p[k][i+1][j] 
						avg+=new_p[k][i-1][j] 
						avg+=new_p[k][i][j-1] 
						avg+=new_p[k][i][j+1] 
					avg-=curl_v[i][j] 
					avg/=div 
					new_p[k+1][i][j]=avg 
			if k==3: 
				p=[[(new_p[2][i][j]+2*new_p[3][i][j]+4*new_p[4][i][j])/7.0 for j in range(len(array[0]))] for i in range(len(array))] 
		#now we can use those divergence/convergence scalars to increase and decrease density 
		#note: positive values in p mean that density is decreasing 
		#negative values in p mean that density is increasing. 
		for i in range(len(array)): 
			for j in range(len(array[0])): 
				#avg=0 
				#div=0.0 
				#note, you can only get density from adjacent points 
				#due to my avoidance of convolutions to maintain 
				#conservation of mass, there does exist the possibility 
				#of creating matter from nothing with this method, 
				#but it should be quickly drawn back by the negative  
				#densities. 
				if i>0 and i<len(array)-1 and j>0 and j<len(array)-1: 
					partial=-p[i][j]*k*array[i-1][j][1]/4 
					if partial>0: 
						n_arr[i][j][1]+=partial
						n_arr[i-1][j][1]-=partial 
					partial=-p[i][j]*k*array[i+1][j][1]/4 
					if partial>0: 
						n_arr[i][j][1]+=partial 
						n_arr[i+1][j][1]-=partial 
					partial=-p[i][j]*k*array[i][j-1][1]/4 
					if partial>0: 
						n_arr[i][j][1]+=partial 
						n_arr[i][j-1][1]-=partial 
					partial=-p[i][j]*k*array[i][j+1][1]/4 
					if partial>0: 
						n_arr[i][j][1]+=partial 
						n_arr[i][j+1][1]-=partial 
		#Now we proceed with using the change in density to cause a change  
		#in temperature 
		for i in range(len(array)): 
			for j in range(len(array[0])): 
				if abs(n_arr[i][j][1]*array[i][j][1])>.001: 
					if n_arr[i][j][1]<0: 
						n_arr[i][j][2]+=min(1, max(.0001, n_arr[i][j][1]-array[i][j][1])) #increase temperature directly with more material 
						#this is not a case that should happen, but this is the least bad (in terms of gameplay) way to go about 
						#fixing the issue (remember that a positive temperature*negative density will mean a negative 
						#pressure that pulls in outer material, which will then, with this method, accelerate until 
						#the current position has positive density (at which point it may explode a bit with the high 
						#temperature) 
					else: 
						n_arr[i][j][2]+=min(1, max(-1, k*(n_arr[i][j][1]-array[i][j][1])/abs(array[i][j][1]))) 
						#this is trying to be consistent with the ideal gas law,  
						#since squeezing more stuff into an area necessarily increases temperature 
						#but I don't want temperature to explode 
		#Now we have to make sure that temperature does not become negative, by subtracting the minimum value from all values 
		#bringing the minimum to 0 always.
		minT=min([i[2] for thing in n_arr for i in thing])
		for i in range(len(array)):
			for j in range(len(array[i])):
				n_arr[i][j][2]-=minT
		return n_arr 	
	def add_tensor(arg1: list, arg2: list): 
		n_arr=[[move(i, [{"x":0, "y":0, "z":0}, 0, 0], .5) for i in j] for j in arg1] 
		for i in range(min(len(arg1), len(arg2))): 
			for j in range(min(len(arg1[i]), len(arg2[i]))): 
				move(arg2[i][j], n_arr[i][j], .5) 
		return n_arr 
	def sub_tensor(arg1: list, arg2: list): 
		n_arr=[[move(i, [{"x":0, "y":0, "z":0}, 0, 0], 1) for i in j] for j in arg1] 
		for i in range(min(len(arg1), len(arg2))): 
			for j in range(min(len(arg1[i]), len(arg2[i]))): 
				in_sp=sqrt(n_arr[i][j][0]["x"]**2+n_arr[i][j][0]["y"]**2+n_arr[i][j][0]["z"]**2)
				#move(arg2[i][j], n_arr[i][j], -1) 
				n_arr[i][j][0]["x"]-=arg2[i][j][0]["x"] 
				n_arr[i][j][0]["y"]-=arg2[i][j][0]["y"] 
				n_arr[i][j][0]["z"]-=arg2[i][j][0]["z"] 
				n_sp=sqrt(n_arr[i][j][0]["x"]**2+n_arr[i][j][0]["y"]**2+n_arr[i][j][0]["z"]**2)
				if (in_sp*n_sp>.000001):
					n_arr[i][j][0]["x"]*=in_sp/n_sp
					n_arr[i][j][0]["y"]*=in_sp/n_sp
					n_arr[i][j][0]["z"]*=in_sp/n_sp
		return n_arr 
	#Now we proceed with the expected computation 
	diff=diffusion(t, array, k) 
	divergence=curl_free(diff, p)
	p=divergence[1]
	divergence=divergence[0]
	curl=sub_tensor(diff, divergence) 
	mix=temperate(curl, k/10)  
	adv=advection(t, mix, k) 
	#adv=add_tensor(diff, adv) 
	#p=[] #put into fluid_step args 
	divergence=curl_free(adv, p)
	p=divergence[1]
	divergence=divergence[0] 
	curl=sub_tensor(adv, divergence) 
	#mix=temperate(curl, k/10)  
	#diffusion(t: float, array: list, k: float): 
	#advection(t: float, array: list, k: float): 
	#def curl_free(array: list, p=[]) 
	#def temperate(array: list, k: float) 
	#def fluid_step(t: float, array: list, k: float, p:list=[]) 
	return (curl, p)

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
# Flip the display
p=[]
fluid_sim[25][20][1]=5
fluid_sim[25][20][2]=10
fluid_sim[200][100][1]=5
fluid_sim[200][100][2]=10
def coloring(pygame, fluid_sim):
	meanD=sum([sum([j[1] for j in i]) for i in fluid_sim])/(len(fluid_sim)**2)
	sd_d=sum([(item[1]-meanD)**2 for sublist in fluid_sim for item in sublist])/(len(fluid_sim)**2-1)
	sd_d=sqrt(sd_d)
	meanT=sum([sum([j[2] for j in i]) for i in fluid_sim])/(len(fluid_sim)**2)
	sd=sum([(item[2]-meanT)**2 for sublist in fluid_sim for item in sublist])/(len(fluid_sim)**2-1)
	sd=sqrt(sd)
	for i in range(len(fluid_sim)):
		for j in range(len(fluid_sim)):
			green=int(min(255, max(0, fluid_sim[i][j][1]-meanD)/sd*10))
			if fluid_sim[i][j][2]>meanT:
				pygame.draw.circle(screen, (int(min(255, (fluid_sim[i][j][2]-meanT)/sd*10)), green, 0), (j*2, i*2), 1)
			elif fluid_sim[i][j][2]<meanT:
				pygame.draw.circle(screen, (0, green, (int(min(255, (meanT-fluid_sim[i][j][2])/sd*10)))), (j*2, i*2), 1)
	
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
		fluid_sim[prev_mouse[1]//2][prev_mouse[0]//2][0]["x"]=xV
		fluid_sim[prev_mouse[1]//2][prev_mouse[0]//2][0]["y"]=yV
		fluid_sim[new_mouse[1]//2][new_mouse[0]//2][0]["x"]=xV
		fluid_sim[new_mouse[1]//2][new_mouse[0]//2][0]["y"]=yV
		fluid_sim[new_mouse[1]//2][new_mouse[0]//2][1]+=sqrt(xV**2+yV**2)
		fluid_sim[new_mouse[1]//2][new_mouse[0]//2][2]=1
	prev_mouse=new_mouse
	stuff=fluid_step(.01, fluid_sim, .08, p)
	fluid_sim=stuff[0]
	p=stuff[1]
	pygame.display.flip()
# Done! Time to quit.
pygame.quit()

