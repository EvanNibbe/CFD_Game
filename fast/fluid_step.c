#define MAX(x, y) (x>y?x:y)
#define MIN(x, y) (x<y?x:y)
#include <string.h> //to get memset, a fast way to 0 arrays.
#include <stdlib.h> //to get the rand function.
#define SIM_SIZE 256
#define SIM_TOT SIM_SIZE*SIM_SIZE*5+SIM_SIZE*SIM_SIZE
#define and &&
#include <math.h>
void move(double *arg, double *dest, double weight) {
	double orig_m=dest[3];
	dest[3]+=arg[3]*weight;
	double orig_t=dest[4]; //temp
	if (fabs(orig_m)<.00001) {
		dest[0]+=arg[0]*weight;
		dest[1]+=arg[1]*weight;
		dest[2]+=arg[2]*weight;
		dest[4]+=arg[4]*weight;
	} else if (fabs(orig_m+dest[3]>.001)) {
		dest[0]=(orig_m*dest[0]+arg[3]*arg[0]*weight)/(orig_m+arg[3]*weight);
		dest[1]=(orig_m*dest[1]+arg[3]*arg[1]*weight)/(orig_m+arg[3]*weight);
		dest[2]=(orig_m*dest[2]+arg[3]*arg[2]*weight)/(orig_m+arg[3]*weight);
		dest[4]=(orig_m*dest[4]+arg[3]*arg[4]*weight)/(orig_m+arg[3]*weight);
	} else {
		dest[4]=(orig_m*dest[4]+arg[3]*arg[4]*weight)/MAX(MAX(.001, MAX(orig_m, arg[3])), MAX(orig_m+arg[3]*weight, MIN(.01, weight)));
	}
}
/*
This code needs to replicate the functionality of
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
*/

double fluid_step_mem[SIM_TOT]; //This is where results are stored.
double diffusion_mem[SIM_SIZE*SIM_SIZE*5]; //diffusion results stored
double *diffusion(double t, double *array, double k) {
	double *n_arr=diffusion_mem;
	memset(n_arr, 0,  SIM_SIZE*SIM_SIZE*5* sizeof (double));
	for (int i=0; i<SIM_SIZE; i++) {
		for (int j=0; j<SIM_SIZE; j++) {
			double avg[5]={0, 0, 0, 0, 0};
			double div=4;
			if (i>0 and i<SIM_SIZE-1) {
				move(array+((i+1)*SIM_SIZE+j)*5, avg, 1);
				move(array+((i-1)*SIM_SIZE+j)*5, avg, 1);
			} else if (i==0) {
				div-=1;
				move(array+((i+1)*SIM_SIZE+j)*5, avg, 1);
			} else {
				div-=1;
				move(array+((i-1)*SIM_SIZE+j)*5, avg, 1);
			}
			if (j>0 and j<SIM_SIZE-1) {
				move(array+((i)*SIM_SIZE+j+1)*5, avg, 1);
				move(array+((i)*SIM_SIZE+j-1)*5, avg, 1);
			} else if (j==0) {
				div-=1;
				move(array+((i)*SIM_SIZE+j+1)*5, avg, 1);
			} else {
				div-=1;
				move(array+((i)*SIM_SIZE+j-1)*5, avg, 1);
			}
			avg[0]/=div;
			avg[1]/=div;
			avg[2]/=div;
			avg[3]/=div;
			avg[4]/=div;
			move(array +(i*SIM_SIZE+j)*5, n_arr+(i*SIM_SIZE+j)*5, 1-4*k/(1+k));
			move(avg, n_arr+(i*SIM_SIZE+j)*5, k/(1+k));
		}
	}
	return n_arr;
}
/*
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
*/
double advection_mem[SIM_SIZE*SIM_SIZE*5];
double* advection(double t, double *array, double k) {
	double *n_arr=advection_mem;
	memset(n_arr, 0,  SIM_SIZE*SIM_SIZE*5* sizeof (double));
	for (int i=0; i<SIM_SIZE; i++) {
		for (int j=0; j<SIM_SIZE; j++) {
			double avg[5]={0, 0, 0, 0, 0};
			double div=4;
			double frm[]={j-array[(i*SIM_SIZE+j)*5]*t, i-array[(i*SIM_SIZE+j)*5+1]*t};
			if (frm[0]>0 and frm[0]<SIM_SIZE-1 and frm[1]>0 and frm[1]<SIM_SIZE-1) {
				double dxA=frm[0]-(int)frm[0];
				double dyA=frm[1]-(int)frm[1];
				double dA=sqrt(dxA*dxA+dyA*dyA);
				double dxB=(int)(frm[0]+1)-frm[0];
				double dyB=dyB-(int)frm[1];
				double dB=sqrt(dxB*dxB+dyB*dyB);
				double dxC=dxA;
				double dyC=(int)(frm[1]+1)-frm[1];
				double dC=sqrt(dxC*dxC+dyC*dyC);
				double dxD=dxB;
				double dyD=dyC;
				double dD=sqrt(dxD*dxD+dyD*dyD);
				double tot=dA+dB+dC+dD;
				move(array+(((int)(frm[1]))*SIM_SIZE+(int)(frm[0]))*5, avg, (dB+dC+dD)/(3.0*tot));
				move(array+(((int)(frm[1])+1)*SIM_SIZE+(int)(frm[0]))*5, avg, (dB+dC+dD)/(3.0*tot));
				move(array+(((int)(frm[1]))*SIM_SIZE+(int)(frm[0]+1))*5, avg, (dB+dC+dD)/(3.0*tot));
				move(array+(((int)(frm[1]+1))*SIM_SIZE+(int)(frm[0]+1))*5, avg, (dB+dC+dD)/(3.0*tot));
				move(avg, n_arr+5*(i*SIM_SIZE+j), 1);
			}
			
			
		}
	}
	return n_arr;
}
/*
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
*/
double curl_free_mem[SIM_TOT+SIM_SIZE*SIM_SIZE];
double *curl_free(double *array, double *p) {
	double *n_arr=curl_free_mem;
	double *prev_p=curl_free_mem+SIM_SIZE*SIM_SIZE*5;
	double *curl_v=prev_p+SIM_SIZE*SIM_SIZE;
	memset(n_arr, 0,  SIM_TOT+SIM_SIZE*SIM_SIZE* sizeof (double));
	for (int i=0; i<SIM_SIZE; i++) {
		for (int j=0; j<SIM_SIZE; j++) {
			double div=2.0;
			double avg=2.0;
			if (i>0 and i<SIM_SIZE-1 and j>0 and j<SIM_SIZE-1) {
				avg+=array[5*((i+1)*SIM_SIZE+j)+1]-array[5*((i-1)*SIM_SIZE+j)+1];
				avg+=array[5*((i+1)*SIM_SIZE+j)]-array[5*((i-1)*SIM_SIZE+j)];
			}
			avg/=div;
			curl_v[i*SIM_SIZE+j]=avg;
		}
	}
	for (int k=0; k<5; k++) {
		memcpy(prev_p, p, sizeof(double)*SIM_SIZE*SIM_SIZE);
		for (int i=0; i<SIM_SIZE; i++) {
			for (int j=0; j<SIM_SIZE; j++) {
				double avg=0;
				double div=4;
				if (i>0 and i<SIM_SIZE-1 and j>0 and j<SIM_SIZE-1) {
					avg+=prev_p[(i+1)*SIM_SIZE+j];
					avg+=prev_p[(i-1)*SIM_SIZE+j];
					avg+=prev_p[(i)*SIM_SIZE+j+1];
					avg+=prev_p[(i)*SIM_SIZE+j-1];
				}
				avg-=curl_v[i*SIM_SIZE+j];
				avg/=div;
				p[i*SIM_SIZE+j]=avg;
				if (k==4 and i>0 and i<SIM_SIZE-1 and j>0 and j<SIM_SIZE-1) {
					double div_p[]={(p[i*SIM_SIZE+j+1]-p[i*SIM_SIZE+j-1])/2, (p[(i+1)*SIM_SIZE+j]-p[(i-1)*SIM_SIZE+j])/2};
					n_arr[(i*SIM_SIZE+j)*5]=array[(i*SIM_SIZE+j)*5]-div_p[0];
					n_arr[(i*SIM_SIZE+j)*5+1]=array[(i*SIM_SIZE+j)*5+1]-div_p[1];
				}
			}
		}
	}
	return curl_free_mem; //the second part of this array (after SIM_SIZE*SIM_SIZE*5) is p
}
/*
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
*/
					//for n_arr			//for delta_v		 //for p			//curl_v		//new_p
double temperate_mem[SIM_SIZE*SIM_SIZE*5+2*SIM_SIZE*SIM_SIZE+SIM_SIZE*SIM_SIZE+SIM_SIZE*SIM_SIZE+SIM_SIZE*SIM_SIZE*5];
double *temperate(double *array, double k) {
	double *n_arr=temperate_mem;
	double *delta_v=n_arr+SIM_SIZE*SIM_SIZE*5;
	double *p=delta_v+2*SIM_SIZE*SIM_SIZE;
	double *curl_v=p+SIM_SIZE*SIM_SIZE;
	double *new_p=curl_v+SIM_SIZE*SIM_SIZE;
	memcpy(n_arr, array, sizeof(double)*SIM_SIZE*SIM_SIZE*5);
	for (int i=0; i<SIM_SIZE; i++) {
		for (int j=0; j<SIM_SIZE; j++) {
			if (i>0 and i<SIM_SIZE-1 and j>0 and j<SIM_SIZE-1) {
				double dM=array[(i*SIM_SIZE+j)*5+3];
				double tM=array[(i*SIM_SIZE+j)*5+4];
				double *Sur[]={array+((i-1)*SIM_SIZE+j-1)*5, array+((i-1)*SIM_SIZE+j)*5, array+((i-1)*SIM_SIZE+j+1)*5, array+((i)*SIM_SIZE+j-1)*5, array+((i)*SIM_SIZE+j+1)*5, array+((i+1)*SIM_SIZE+j-1)*5, array+((i+1)*SIM_SIZE+j)*5, array+((i+1)*SIM_SIZE+j+1)*5};
				if (dM*tM>Sur[0][3]*Sur[0][4]) {
					double coef=MIN(2*k, k*(dM*tM/MAX(.1, Sur[0][3]*Sur[0][4])-1))/sqrt(2); //sqrt(2) is distance
					delta_v[((i-1)*SIM_SIZE+j-1)*2]+=coef*(-1); //-1 is dir along x-axis
					delta_v[((i-1)*SIM_SIZE+j-1)*2+1]+=coef*(-1);
				}
				if (dM*tM>Sur[1][3]*Sur[1][4]) {
					double coef=MIN(2*k, k*(dM*tM/MAX(.1, Sur[1][3]*Sur[1][4])-1));
					delta_v[((i-1)*SIM_SIZE+j)*2]+=0;
					delta_v[((i-1)*SIM_SIZE+j)*2+1]+=coef*(-1);
				}
				if (dM*tM>Sur[2][3]*Sur[2][4]) {
					double coef=MIN(2*k, k*(dM*tM/MAX(.1, Sur[2][3]*Sur[2][4])-1))/sqrt(2); //sqrt(2) is distance
					delta_v[((i-1)*SIM_SIZE+j+1)*2]+=coef*(1); //1 is dir along x-axis
					delta_v[((i-1)*SIM_SIZE+j+1)*2+1]+=coef*(-1);
				}
				if (dM*tM>Sur[3][3]*Sur[3][4]) {
					double coef=MIN(2*k, k*(dM*tM/MAX(.1, Sur[3][3]*Sur[3][4])-1));
					delta_v[((i)*SIM_SIZE+j-1)*2]+=coef*(-1); //-1 is dir along x-axis
					delta_v[((i)*SIM_SIZE+j-1)*2+1]+=0;
				}
				if (dM*tM>Sur[4][3]*Sur[4][4]) {
					double coef=MIN(2*k, k*(dM*tM/MAX(.1, Sur[4][3]*Sur[4][4])-1));
					delta_v[((i)*SIM_SIZE+j+1)*2]+=coef*(1); //1 is dir along x-axis
					delta_v[((i)*SIM_SIZE+j+1)*2+1]+=0;
				}
				if (dM*tM>Sur[5][3]*Sur[5][4]) {
					double coef=MIN(2*k, k*(dM*tM/MAX(.1, Sur[5][3]*Sur[5][4])-1))/sqrt(2); //sqrt(2) is distance
					delta_v[((i+1)*SIM_SIZE+j-1)*2]+=coef*(-1); //-1 is dir along x-axis
					delta_v[((i+1)*SIM_SIZE+j-1)*2+1]+=coef*(1);
				}
				if (dM*tM>Sur[6][3]*Sur[6][4]) {
					double coef=MIN(2*k, k*(dM*tM/MAX(.1, Sur[6][3]*Sur[6][4])-1));
					delta_v[((i+1)*SIM_SIZE+j)*2]+=0;
					delta_v[((i+1)*SIM_SIZE+j)*2+1]+=coef*(1);
				}
				if (dM*tM>Sur[7][3]*Sur[7][4]) {
					double coef=MIN(2*k, k*(dM*tM/MAX(.1, Sur[7][3]*Sur[7][4])-1))/sqrt(2); //sqrt(2) is distance
					delta_v[((i+1)*SIM_SIZE+j+1)*2]+=coef*(1); //1 is dir along x-axis
					delta_v[((i+1)*SIM_SIZE+j+1)*2+1]+=coef*(1);
				}
			}
		}
	}
	//now we slightly randomize velocities according to temperature
	//change only the direction, not the magnitude.
	for (int i=0; i<SIM_SIZE; i++) {
		for (int j=0; j<SIM_SIZE; j++) {
			double mag=sqrt(delta_v[(i*SIM_SIZE+j)*2]*delta_v[(i*SIM_SIZE+j)*2]+delta_v[(i*SIM_SIZE+j)*2+1]*delta_v[(i*SIM_SIZE+j)*2+1]);
			if (mag>0.001) {
				delta_v[(i*SIM_SIZE+j)*2]+=(2.0-(rand()+.000001)/(RAND_MAX+0.00001))*k*array[(i*SIM_SIZE+j)*5+4]; //multiply by temperature
				delta_v[(i*SIM_SIZE+j)*2+1]+=(2.0-(rand()+.000001)/(RAND_MAX+0.00001))*k*array[(i*SIM_SIZE+j)*5+4];
				double n_dist=sqrt(delta_v[(i*SIM_SIZE+j)*2]*delta_v[(i*SIM_SIZE+j)*2]+delta_v[(i*SIM_SIZE+j)*2+1]*delta_v[(i*SIM_SIZE+j)*2+1]);
				if (n_dist<.00001) {
					delta_v[(i*SIM_SIZE+j)*2]=mag;
				} else {
					delta_v[(i*SIM_SIZE+j)*2]=delta_v[(i*SIM_SIZE+j)*2]*mag/n_dist;
					delta_v[(i*SIM_SIZE+j)*2+1]=delta_v[(i*SIM_SIZE+j)*2+1]*mag/n_dist;
				}
			}
		}
	}
	for (int i=0; i<SIM_SIZE; i++) {
		for (int j=0; j<SIM_SIZE; j++) {
			double div=2;
			double avg=0;
			if (i>0 and i<SIM_SIZE-1 and j>0 and j<SIM_SIZE-1) {
				avg+=delta_v[((i+1)*SIM_SIZE+j)*2+1]-delta_v[((i-1)*SIM_SIZE+j)*2+1];
				avg+=delta_v[((i)*SIM_SIZE+j+1)*2]-delta_v[((i)*SIM_SIZE+j-1)*2];
			}
			avg/=div;
			curl_v[i*SIM_SIZE+j]=avg;
		}
	}
	for (int k=0; k<4; k++) {
		for (int i=0; i<SIM_SIZE; i++) {
			for (int j=0; j<SIM_SIZE; j++) {
				double div=4;
				double avg=0;
				if (i>0 and i<SIM_SIZE-1 and j>0 and j<SIM_SIZE-1) {
					avg+=new_p[k*SIM_SIZE*SIM_SIZE+(i+1)*SIM_SIZE+(j)];
					avg+=new_p[k*SIM_SIZE*SIM_SIZE+(i-1)*SIM_SIZE+(j)];
					avg+=new_p[k*SIM_SIZE*SIM_SIZE+(i)*SIM_SIZE+(j-1)];
					avg+=new_p[k*SIM_SIZE*SIM_SIZE+(i)*SIM_SIZE+(j+1)];
				}
				avg-=curl_v[i*SIM_SIZE+j];
				avg/=div;
				new_p[(k+1)*SIM_SIZE*SIM_SIZE+(i)*SIM_SIZE+(j)]=avg;
				if (k==3) {
					p[SIM_SIZE*i+j]=new_p[2*SIM_SIZE*SIM_SIZE+i*SIM_SIZE+j]+2*new_p[3*SIM_SIZE*SIM_SIZE+i*SIM_SIZE+j]+4*new_p[4*SIM_SIZE*SIM_SIZE+i*SIM_SIZE+j];
					p[SIM_SIZE*i+j]/=7.0;
				}
			}
		}
	}
	//we can now use those convergence/divergence scalars to increase and decrease density
	//Note: positive values in p mean that density is decreasing
	//negative values in p mean that density is increasing
	for (int i=0; i<SIM_SIZE; i++) {
		for (int j=0; j<SIM_SIZE; j++) {
			if (i>0 and i<SIM_SIZE-1 and j>0 and j<SIM_SIZE-1) {
				double partial=-p[i*SIM_SIZE+j]*k*array[((i-1)*SIM_SIZE+j)*5+3];
				if (partial>0) {
					n_arr[(i*SIM_SIZE+j)*5+3]+=partial;
					n_arr[((i-1)*SIM_SIZE+j)*5+3]-=partial;
				}
				partial=-p[i*SIM_SIZE+j]*k*array[((i+1)*SIM_SIZE+j)*5+3];
				if (partial>0) {
					n_arr[(i*SIM_SIZE+j)*5+3]+=partial;
					n_arr[((i+1)*SIM_SIZE+j)*5+3]-=partial;
				}
				partial=-p[i*SIM_SIZE+j]*k*array[((i)*SIM_SIZE+j-1)*5+3];
				if (partial>0) {
					n_arr[(i*SIM_SIZE+j)*5+3]+=partial;
					n_arr[((i)*SIM_SIZE+j-1)*5+3]-=partial;
				}
				partial=-p[i*SIM_SIZE+j]*k*array[((i)*SIM_SIZE+j+1)*5+3];
				if (partial>0) {
					n_arr[(i*SIM_SIZE+j)*5+3]+=partial;
					n_arr[((i)*SIM_SIZE+j+1)*5+3]-=partial;
				}
			}
		}
	}
	//Now we proceed with using the change in density to
	//change the temperature
	for (int i=0; i<SIM_SIZE; i++) {
		for (int j=0; j<SIM_SIZE; j++) {
			if (fabs(n_arr[(i*SIM_SIZE+j)*5+3]*array[(i*SIM_SIZE+j)*5+3])>.001) {
				if (n_arr[(i*SIM_SIZE+j)*5+3]<0) {
					n_arr[(i*SIM_SIZE+j)*5+4]+=MIN(1, MAX(.0001, n_arr[(i*SIM_SIZE+j)*5+3]-array[(i*SIM_SIZE+j)*5+3]));
					//increase temperature directly with more material
					//this is not a case that should happen, but with negative density and
					//positive temperature, material should be drawn to this point.
				} else {
					n_arr[(i*SIM_SIZE+j)*5+4]+=MIN(1, MAX(-1, k*(n_arr[(i*SIM_SIZE+j)*5+3]-array[(i*SIM_SIZE+j)*5+3])/(fabs(array[(i*SIM_SIZE+j)*5+3])+.0000001)));
					//this is trying to be consistent with the ideal gas law, but I don't want temperature values
					//to blow up due to 0 densities (or close to 0 densities)
				}
			}
		}
	}
	//Now we have to make sure that temperature does not become negative,
	//by subtracting the minimum value from all values
	//bringing the minimum to 0 always.
	double minT=n_arr[0*SIM_SIZE+0+4];
	for (int i=0; i<SIM_SIZE*SIM_SIZE; i++) {
		if (n_arr[i*5+4]<minT) {
			minT=n_arr[i*5+4];
		}
	}
	for (int i=0; i<SIM_SIZE*SIM_SIZE; i++) {
		n_arr[i*5+4]-=minT;
	}
	return n_arr;
}
/*
	def temperate(array: list, k: float): #use temperature and pressure vs inverse viscosity k 
		#then use density*temperature compared to adjacent spot to determine outward flow, temperature may be negative 
		#if there is outward flow, then the adjacent spots affected get their velocity changed by a small amount according to  
		#the gradient, and the positive temperature randomizes the resulting velocity's angle by that many degrees. 
		#a small k acts as a moderating force. 
		n_arr=[[move(i, [{"x":0, "y":0, "z":0}, 0, 0], 1) for i in j] for j in array] 
		delta_v=[[[0, 0] for i in j] for j in array] 
		<removing original Python to make space to see the rest of temperate>
		
		minT=min([i[2] for thing in n_arr for i in thing])
		for i in range(len(array)):
			for j in range(len(array[i])):
				n_arr[i][j][2]-=minT
		return n_arr 	
*/
double sub_tensor_mem[SIM_SIZE*SIM_SIZE*5];
double *sub_tensor(double *arg1, double *arg2) {
	double *n_arr=sub_tensor_mem;
	memcpy(n_arr, arg1, SIM_SIZE*SIM_SIZE*5*sizeof(double));
	for (int i=0; i<SIM_SIZE*SIM_SIZE; i++) {
		double in_sp=sqrt(n_arr[i*5]*n_arr[i*5]+n_arr[i*5+1]*n_arr[i*5+1]+n_arr[i*5+2]*n_arr[i*5+2]);
		n_arr[i*5]-=arg2[i*5];
		n_arr[i*5+1]-=arg2[i*5+1];
		n_arr[i*5+2]-=arg2[i*5+1];
		double n_sp=sqrt(n_arr[i*5]*n_arr[i*5]+n_arr[i*5+1]*n_arr[i*5+1]+n_arr[i*5+2]*n_arr[i*5+2]);
		if (in_sp*n_sp>.000001) {
			n_arr[i*5]*=in_sp/n_sp;
			n_arr[i*5+1]*=in_sp/n_sp;
			n_arr[i*5+2]*=in_sp/n_sp;
		}
	}
	return n_arr;
}
/*
	def sub_tensor(arg1: list, arg2: list): 
		n_arr=[[move(i, [{"x":0, "y":0, "z":0}, 0, 0], 1) for i in j] for j in arg1] 
		for i in range(min(len(arg1), len(arg2))): 
			for j in range(min(len(arg1[i]), len(arg2[i]))): 
				in_sp=sqrt(n_arr[i][j][0]["x"]**2+n_arr[i][j][0]["y"]**2+n_arr[i][j][0]["z"]**2)
				#move(arg2[i][j], n_arr[i][j], -1) 
				n_arr[i][j][0]["x"]-=arg2[i][j][0]["x"] 
				n_arr[i][j][0]["y"]-=arg2[i][j][(array+(i_sp
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
	divergence=curl_free(adv, p)
	p=divergence[1]
	divergence=divergence[0] 
	curl=sub_tensor(adv, divergence)  
	#diffusion(t: float, array: list, k: float): 
	#advection(t: float, array: list, k: float): 
	#def curl_free(array: list, p=[]) 
	#def temperate(array: list, k: float) 
	#def fluid_step(t: float, array: list, k: float, p:list=[]) 
	return (curl, p)
	*/
/*
Evan Nibbe
fluid_step.c
May 31, 2022
This code needs to be compiled with 
cc -fPIC -shared -o fluid_step.so fluid_step.c
in order to be used with Python
the C function returns the pointer to an array, so we use ndpointer as return restype.
#in Python:
from ctypes import c_void_p, c_double, c_int, cdll
import numpy as np
from numpy.ctypeslib import ndpointer
lib=cdll.LoadLibrary("fluid_step.so")
fluid_step=lib.fluid_step #fluid_step is the name of our C function
fluid_step.restype=ndpointer(dtype=c_double, shape=(n,))
#result=fluid_step(c_void_p(matrix.ctypes.data),c_int(n),c_int(m))
#line 326 is wrong in our case because the C function declaration is different from the medium post by Matias Aravena Gamboa Feb 26, 2020
#However, we find it useful, particularly matrix.ctypes.data, 
#Matrix is the 2D array he is using, but because it is an np.array
#it has an attribute called ctypes, and of the ctypes attribute is
#an attribute called data, which has the pointer that C uses
*/

/*
Implementation notes:
The array is 256 by 256, with 5 values at each location: x velocity, y velocity, z velocity, density, temperature
327680 doubles
We then have a 256 by 256 double array p, which just has a single double at the locations.
What we will need to do on the Python side is call an additional function to get back the pointer 
to the start of p, since Python can't do pointer arithmetic.

Danger: Do not use this with multithreading of different fluid steps at once, since all rely on the same static array for results
to avoid memory leaks.
*/

double *get_p() {
	return fluid_step_mem+SIM_SIZE*SIM_SIZE*5;
}
double get_meanT() { //mean temp
	double res=0;
	for (int i=1; i<=SIM_SIZE*SIM_SIZE; i++) {
		res+=fluid_step_mem[i*5-1];
	}
	res/=SIM_SIZE*SIM_SIZE;
	return res;
}
double get_meanD() { //mean density
	double res=0;
	for (int i=1; i<=SIM_SIZE*SIM_SIZE; i++) {
		res+=fluid_step_mem[i*5-2];
	}
	res/=SIM_SIZE*SIM_SIZE;
	return res;
}
double get_sd_T(double meanT) {
	double res=0;
	for (int i=1; i<=SIM_SIZE*SIM_SIZE; i++) {
		double diff=(fluid_step_mem[i*5-1]-meanT);
		res+=diff*diff;
	}
	res/=SIM_SIZE*SIM_SIZE-1;
	res=sqrt(res);
	return res;
}
double get_sd_D(double meanD) {
	double res=0;
	for (int i=1; i<=SIM_SIZE*SIM_SIZE; i++) {
		double diff=(fluid_step_mem[i*5-2]-meanD);
		res+=diff*diff;
	}
	res/=SIM_SIZE*SIM_SIZE-1;
	res=sqrt(res);
	return res;
}
double cross_prod_mem[3];
double* cross_prod(double v1[3], double *v2) {
	double *n_arr=cross_prod_mem;
	/*

    cx = aybz − azby
    cy = azbx − axbz
    cz = axby − aybx

	*/
	n_arr[0]=v1[1]*v2[2]-v1[2]*v2[1];
	n_arr[1]=v1[2]*v2[0]-v1[0]*v2[2];
	n_arr[2]=v1[0]*v2[1]-v1[1]*v2[0];
	return n_arr;
}
//this magnetic function seems to be causing a problem of wiping everything off the map.
void magnetic(double x1, double y1, double z1, double x2, double y2,double z2) {
	double *n_arr=fluid_step_mem; //changing these velocities
	double v1[]={x2-x1, y2-y1, z2-z1};
	//double *force=cross_prod(v1, n_arr+(i*SIM_SIZE+j)*5)
	for (int i=0; i<SIM_SIZE; i++) {
		for (int j=0; j<SIM_SIZE; j++) {
			double dist=fabs((j-x1)*(j-x1)*(j-x2)*(j-x2)+(i-y1)*(i-y1)*(i-y2)*(i-y2));
			double v[]={v1[0]/dist, v1[1]/dist, v1[2]/dist};
			double *force=cross_prod(v, n_arr+(i*SIM_SIZE+j)*5);
			n_arr[(i*SIM_SIZE+j)*5]+=force[0];
			n_arr[(i*SIM_SIZE+j)*5+1]+=force[1];
			//maybe this is what is causing falls to 0: 
			//actually not, what was causing falls to 0 was the way that move was placed outside the if statement in advection.
			n_arr[(i*SIM_SIZE+j)*5+2]+=force[2]*(2-rand()/(.00001+RAND_MAX))/1000;
		}
	}
}
double *fluid_step(double t, double *array, double k, double *p) {
	/*
	#Now we proceed with the expected computation 
	diff=diffusion(t, array, k) 
	divergence=curl_free(diff, p)
	p=divergence[1]
	divergence=divergence[0]
	curl=sub_tensor(diff, divergence) 
	mix=temperate(curl, k/10)  
	adv=advection(t, mix, k)  
	divergence=curl_free(adv, p)
	p=divergence[1]
	divergence=divergence[0] 
	curl=sub_tensor(adv, divergence)  
	#diffusion(t: float, array: list, k: float): 
	#advection(t: float, array: list, k: float): 
	#def curl_free(array: list, p=[]) 
	#def temperate(array: list, k: float) 
	#def fluid_step(t: float, array: list, k: float, p:list=[]) 
	return (curl, p)
	*/
	double *n_arr=fluid_step_mem;
	double *new_p=n_arr+SIM_SIZE*SIM_SIZE*5;
	double *diff=diffusion(t, array, k);
	double *divergence=curl_free(diff, p);
	double *div_p=divergence+SIM_SIZE*SIM_SIZE*5;
	memcpy(new_p, div_p, SIM_SIZE*SIM_SIZE*sizeof(double));
	double *curl=sub_tensor(diff, divergence);
	double *mix=temperate(curl, k/10);
	double *adv=advection(t, mix, k);
	divergence=curl_free(adv, new_p);
	curl=sub_tensor(adv, divergence);
	memcpy(n_arr, curl, SIM_SIZE*SIM_SIZE*5*sizeof(double));
	memcpy(new_p, div_p, SIM_SIZE*SIM_SIZE*sizeof(double));
	memcpy(p, div_p, SIM_SIZE*SIM_SIZE*sizeof(double));
	return n_arr;

}
