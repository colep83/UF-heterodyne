# this code will apply the three PSAs, and then average both the intensity plots and the various phase plots.
# GW log post 13590 on PSA is good reference

# ** CURRENTLY: **
# 
# NOVAK adds all phase plots to an array (NOVAK_PHA_PLT) 
# CARRE displays first phase plot
# FOUR POINT displays all phase plots
# all are saved into the Phase_Camera_Images folder

import os
import PySpin
import numpy
import matplotlib.pyplot as plt
import time
from PIL import Image
from numpy import array, empty, ravel, where, ones, reshape, arctan2
from matplotlib.pyplot import plot, draw, show, ion

#read in some number of images 
File_Name = input("Enter the name of the file you wish to retrieve the images from: ") #File_Name in save_aquire
NUM_IMAGES = int(input("Please enter the number of images in the file: "))

picList = []
NOVAK_PHA_PLT =[]	

#averaging method 
def Avg(picMatrix, NUM): #ahh this wont work for phase maps as they are just saved as plots
	result = 2
	avgMat = numpy.empty(388800) #this is the number of pixels in the camera
	
		#do we want to threshold anything or take the average of all the pixels?
	for i in range(388800):
		for j in range(NUM):
			avgMat[i] += picMatrix[j][i] #the value of the jth image array at pixel i 
			
		avgMat[i] = avgmat[i]/NUM
			
	avgMat = numpy.reshape(phase,(540,720))
		
	plt.show()
	#maybe do this outside of the method
	name = 'C:/Users/localadmin/Desktop/Phase_Camera_Images/MyFirstMEGAPHASEPLOT'
	plt.savefig(name)
		
	
	return result #so this is probably going to be really cluncky. (reshape, threshold could help)
	

#standard deviation 


#NOVAK, overlap of 1 is good for 90 degree phase step 
def NOVAK(picList, NUM):
	print('NOVAK\n')
	c = 0 
	for i in range(NUM):
		if i % 4 == 1 and i > 1:
			c += 1
			arr1 = numpy.ravel(numpy.array(picList[i-4],dtype='int'))
			arr2 = numpy.ravel(numpy.array(picList[i-3],dtype='int')) #converts to numpy arrays for faster operation
			arr3 = numpy.ravel(numpy.array(piclist[i-2],dtype='int'))
			arr4 = numpy.ravel(numpy.array(picList[i-1],dtype='int'))
			arr5 = numpy.ravel(numpy.array(picList[i],dtype='int'))
			
			# how much threshold would be good? how about none for now
			phase = numpy.empty(388800)
			
			den = 2*arr3-arr1-arr5
			A = arr2-arr4
			B = arr1-arr5
			num = numpy.sqrt(abs(4*A**2-B**2))
			pm = numpy.sign(A)
			phase = numpy.arctan2(pm*num,den)
			
			NOVAK_PHA_PLT.append(phase)
			
			phase = numpy.reshape(phase,(540,720))
			
			#save externally 
			pltName = 'C:/Users/localadmin/Desktop/Phase_Camera_Images/' + File_Name + 'NOVAK' + '-%d.jpg' % c
			plt.savefig(pltName)
	return c
			
	
#CARRE
def CARRE(picList, NUM):
	print('CARRE\n')
	c = 0 
	for i in range(NUM):
		if i % 4 == 0 and i > 0:
			c += 1
			arr1 = numpy.ravel(numpy.array(picList[i-3],dtype='int'))
			arr2 = numpy.ravel(numpy.array(picList[i-2],dtype='int')) #converts to numpy arrays for faster operation
			arr3 = numpy.ravel(numpy.array(piclist[i-1],dtype='int'))
			arr4 = numpy.ravel(numpy.array(picList[i],dtype='int'))
			
			# how much threshold would be good? how about none for now
			phase = numpy.empty(388800)
			
			den = arr2 + arr3 - arr1 - arr4
			A = arr2 - arr3
			B = arr1 - arr4 
			num = (A+B)*(3*A-B)
			num = numpy.sqrt(abs(num))
			pm = numpy.sign(A)
			phase = numpy.arctan2(pm*num,den)
			
			phase = numpy.reshape(phase,(540,720))
			
			#save externally 
			pltName = 'C:/Users/localadmin/Desktop/Phase_Camera_Images/' + File_Name + 'CARRE' + '-%d.jpg' % c
			plt.savefig(pltName)
			
			#ehh not every PSA plot needs to be shown. Just copy and paste if you want to see these
			if c == 1:
				#print(phase.dtype)
				#print(numpy.shape(phase))
				plt.ion()						
				plt.imshow(phase, cmap = 'jet')
				#cbar = plt.colorbar()
				#plt.clim(vmin=-numpy.pi,vmax=numpy.pi)
				#cbar.set_label("Phase Shift (rad)")
				plt.xlabel("Pixels(x)")
				plt.ylabel("Pixels(y)")
				#plt.xlim([300,500])
				#plt.ylim([200,400])
				plt.pause(0.00001)                            
				plt.show()
				plt.clf()
	return c
				

#Four-Point
def FourPoint(picList, NUM):
	print('FourPoint\n')
	c = 0 
	for i in range(NUM):
		if i % 4 == 0 and i > 0:
			c += 1
			arr1 = numpy.ravel(numpy.array(picList[i-3],dtype='int'))
			arr2 = numpy.ravel(numpy.array(picList[i-2],dtype='int')) #converts to numpy arrays for faster operation
			arr3 = numpy.ravel(numpy.array(piclist[i-1],dtype='int'))
			arr4 = numpy.ravel(numpy.array(picList[i],dtype='int'))
			
			# how much threshold would be good? how about none for now
			phase = numpy.empty(388800)
			
			num = arr4 - arr2
			den = arr3 - arr1
			phase = numpy.arctan2(num,den)
			
			phase = numpy.reshape(phase,(540,720))
		
			#save externally, move this after all that other plt stuff? 
			pltName = 'C:/Users/localadmin/Desktop/Phase_Camera_Images/' + File_Name + 'FourPoint' + '-%d.jpg' % c
			plt.savefig(pltName)
			
			#show plot ( maybe add a method for this?)
			#print(phase.dtype)
			#print(numpy.shape(phase))
			plt.ion()						
			plt.imshow(phase, cmap = 'jet')
			#cbar = plt.colorbar()
			#plt.clim(vmin=-numpy.pi,vmax=numpy.pi)
			#cbar.set_label("Phase Shift (rad)")
			plt.xlabel("Pixels(x)")
			plt.ylabel("Pixels(y)")
			#plt.xlim([300,500])
			#plt.ylim([200,400])
			plt.pause(0.00001)                            
			plt.show()
			plt.clf()
	return c
			

#unwrapping (is there a numpy method?????)

#plots (3D?)

#do I need a main or can I just code and call methods as needed 
def main():
	result = True
	for i in range(NUM_IMAGES):
		result = True
		try:
			#opening each image and assigning to an array
			imgName = 'C:/Users/localadmin/Desktop/Phase_Camera_Images/' + File_Name + '-%d.jpg' % i
			IMG = Image.open(imgName)
			IMG.show() # not a single image is being shown so hmm
			imgArray = IMG.GetNDArray()
			picList.append(imgArray)
			if i == 123:
				print (imgArray)
			
		except PySpin.SpinnakerException as ex: # so PySpin shouldn't be needed in this but it is here, and I do not know how to make a TRY/EXCEPT otherwise, someone smarter should clean this up :) 
			print('Error: %s' % ex)
			return False
			
			
	FourPoint(picList, NUM)
	CARRE(picList, NUM)
	C = NOVAK(picList, NUM)
	Avg(NOVAK_PHA_PLT, C)
	return results 