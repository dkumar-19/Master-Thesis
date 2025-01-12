#!/usr/bin/python3.8

import os
import sys
import re
import pprint
#~ import matplotlib
#~ print (matplotlib.__version__)
#~ print (matplotlib.__file__)
import pandas as pd
import datetime

from OpenADB.MVA import *
from OpenADB.Time import *
from OpenADB.Retracking.ImprovedThreshold import ImprovedThreshold
import matplotlib.pylab as plt


# Please input path number, latitude range, threshold for range calculation as well path to in-situ file

passnr = "193"
#cycle = "102"
min_lat = 41.45
max_lat = 42.7
#range_threshold = 85
data=pd.read_csv("/home/deepankar/test/In-situ/Lake Erie/height.csv")



mva_path = "/DGFI8/A/MVA/"

def outlier_rejection(height,jday,stdheight,points):
	median_height=np.median(height)
	crange=np.array(height-median_height)
	q75,q25=np.percentile(crange,[75,25])
	intr_qr=q75-q25
	
	max_threshold=q75+(1.5*intr_qr)
	min_threshold=q25-(1.5*intr_qr)
	
	for x in crange:
		if abs(x)>max_threshold or abs(x)<min_threshold:
			del jday[height.index(x+median_height)]
			del stdheight[height.index(x+median_height)]
			del points[height.index(x+median_height)]
			height.remove(x+median_height)
			continue
	return height,jday,stdheight,points

def height_function(mission,passnr,min_lat,max_lat,range_threshold):
	cycles = sorted(os.listdir(mva_path+'/'+mission))
	points=[]
	stdheight=[]
	mean_heightmatrix= []
	jday_matrix = []
	points_retracked=[]
	stdheight_retracked=[]
	jday_matrix_retracked=[]
	mean_heightmatrix_retracked=[]

	for c in cycles: 	
		if re.match('^[0-9]{3}$',c) == None:
			continue
#extracting the specific files from mva directory for orbit height, instrument error, geoid height, wet and dry troposphere delay, ionospheric delay, tide delay and time (to plot the time sOntarios) for the specific path of the jason 3 mission. 
		path_orbit = mva_path+'/'+mission+'/'+c+'/'+c+'_'+passnr+'orbit.00'
		if not os.path.isfile(path_orbit):
			print (path_orbit,'not found!')
			continue 
		path_instr = mva_path+'/'+mission+'/'+c+'/'+c+'_'+passnr+'instr.00'
		if not os.path.isfile(path_instr):
			print (path_instr,'not found!')
			continue 
		path_geoh = mva_path+'/'+mission+'/'+c+'/'+c+'_'+passnr+'geoh.10'
		if not os.path.isfile(path_geoh):
			print (path_geoh,'not found!')
			continue 
		path_tropw = mva_path+'/'+mission+'/'+c+'/'+c+'_'+passnr+'tropw.06'
		if not os.path.isfile(path_tropw):
			print (path_tropw,'not found!')
			continue 
		path_tropd = mva_path+'/'+mission+'/'+c+'/'+c+'_'+passnr+'tropd.06'
		if not os.path.isfile(path_tropd):
			print (path_tropd,'not found!')
			continue 
		path_ionos = mva_path+'/'+mission+'/'+c+'/'+c+'_'+passnr+'ionos.09'
		if not os.path.isfile(path_ionos):
			print (path_ionos,'not found!')
			continue 
		path_tidep = mva_path+'/'+mission+'/'+c+'/'+c+'_'+passnr+'tidep.06'
		if not os.path.isfile(path_tidep):
			print (path_tidep,'not found!')
			continue 
		path_tides = mva_path+'/'+mission+'/'+c+'/'+c+'_'+passnr+'tides.06'
		if not os.path.isfile(path_tides):
			print (path_tides,'not found!')
			continue 
		path_time = mva_path+'/'+mission+'/'+c+'/'+c+'_'+passnr+'time.00'
		if not os.path.isfile(path_time):
			print(path_time,'not found')
			continue
		path_uralt = mva_path+'/'+mission+'/'+c+'/'+c+'_'+passnr+'uralt.00'
		if not os.path.isfile(path_uralt):
			print(path_uralt,'not found')
			continue
		path_waveform = None
		if mission == 'jason3_f_hf':
			path_waveform = mva_path+'/'+mission+'/'+c+'/'+c+'_'+passnr+'waveform104.00'
			if not os.path.isfile(path_waveform):
				print(path_waveform,'not found')
				continue
		if mission == 'sentinel6a_LR_NTC_F_v2_hf':
			path_waveform = mva_path+'/'+mission+'/'+c+'/'+c+'_'+passnr+'waveform256.00'
			if not os.path.isfile(path_waveform):
				print(path_waveform,'not found')
				continue
		if mission == 'sentinel6a_HR_NTC_F_v2_hf':
			path_waveform = mva_path+'/'+mission+'/'+c+'/'+c+'_'+passnr+'waveform512.00'
			if not os.path.isfile(path_waveform):
				print(path_waveform,'not found')
				continue



#specifying the min and max latitude of observational areas.
		options = {}
		options['glat'] = ('+',min_lat,max_lat)
		records = mask_MVA(path_orbit,options)
	
#extracting the values specific for the latitudes mentioned
		options = {}
		options['records'] = records
		obj_orbit = MVA(path_orbit,options)
		obj_instr = MVA(path_instr,options)
		obj_geoh = MVA(path_geoh,options)
		obj_tropw = MVA(path_tropw,options)
		obj_tropd = MVA(path_tropd,options)
		obj_tidep = MVA(path_tidep,options)
		obj_tides = MVA(path_tides,options)
		obj_ionos = MVA(path_ionos,options)
		obj_time = MVA(path_time,options)
		obj_uralt = MVA(path_uralt,options)
		obj_waveform = MVA(path_waveform,options)

		jday = np.array(obj_time['jday'])
		if len(jday)==0:
			del jday
			continue

		jday_average=np.mean(jday)

		rralt = []
	
		for j in range(0,len(obj_uralt['uralt'])):
			uralt = obj_uralt['uralt'][j]
			waveform = obj_waveform['waveform'][j]
		
			input = {}
			input['waveform'] = waveform
			input['uralt'] = uralt
			if mission == 'jason3_f_hf':
				input['refbin'] = 31.0
				input['binsize'] = 3.125e-9
		
			elif mission == 'sentinel6a_LR_NTC_F_v2_hf':
				input['refbin'] = 128.0
				input['binsize'] = 1.0/395000000.0
			
			elif mission == 'sentinel6a_HR_NTC_F_v2_hf':
				input['refbin'] = 256.0
				input['binsize'] = (1.0/395000000.0)/2.0

			input['threshold'] = range_threshold
			input['option'] = 'best'
		
			retracker = ImprovedThreshold(input)
		#retracker.plot()
			rralt.append(retracker.range)
	
		rralt = np.array(rralt)
	
	
	

		glon = np.array(obj_orbit['glon'])
		glat = np.array(obj_orbit['glat'])
		hsat = np.array(obj_orbit['hsat'])


		ralt = np.array(obj_instr['ralt'])
	

		geoh = np.array(obj_geoh['geoh'])


		wtrop = np.array(obj_tropw['wtrop'])


		dtrop = np.array(obj_tropd['dtrop'])   


		ptide = np.array(obj_tidep['ptide'])
	

		etide = np.array(obj_tides['etide'])


		ionos = np.array(obj_ionos['ionos'])
	
		height = hsat - ralt - geoh - wtrop - dtrop - ptide - etide - ionos
	
		height_retracked = hsat - rralt - geoh - wtrop - dtrop - ptide - etide - ionos
		jday_average_retracked=jday_average
	
		height=height[np.logical_not(np.isnan(height))]
		if len(height)==0:
			del height
			del jday_average
			continue	 
	
		median_height=np.median(height)
		crange=np.array(height-median_height)
		q75,q25 = np.percentile(crange,[75,25])
		intr_qr=q75-q25

		max_threshold=q75+(1.5*intr_qr)
		min_threshold=q25-(1.5*intr_qr)
	
		lst_height=height.tolist()
		lst_range=crange.tolist()
		for x in lst_range:
			if abs(x)>max_threshold or abs(x)<min_threshold:
				if not (x+median_height) in lst_height:
					print(x+median_height,' not in height')
					continue
				else:
					lst_height.remove(x+median_height)
			
				continue
				
		height=np.array(lst_height)
		std_height=np.std(height)
		stdheight.append(std_height)
		points.append(len(height))
		mean_height=np.mean(height)
		mean_heightmatrix.append(mean_height)	
		jday_matrix.append(jday_average)
	
	
		height_retracked=height_retracked[np.logical_not(np.isnan(height_retracked))]
		if len(height_retracked)==0:
			del height_retracked
			del jday_average_retracked
			continue	 
	
		median_height_retracked=np.median(height_retracked)
		crange=np.array(height_retracked-median_height_retracked)
		q75,q25 = np.percentile(crange,[75,25])
		intr_qr=q75-q25

		max_threshold=q75+(1.5*intr_qr)
		min_threshold=q25-(1.5*intr_qr)
		
		lst_height_retracked=height_retracked.tolist()
		lst_crange=crange.tolist()
		for x in lst_crange:
			if abs(x)>max_threshold or abs(x)<min_threshold:
				if not (x+median_height_retracked) in lst_height_retracked:
					print(x+median_height_retracked,' not in height')
					continue
				else:
					lst_height_retracked.remove(x+median_height_retracked)
				continue
				
		height_retracked=np.array(lst_height_retracked)
		std_height=np.std(height_retracked)
		stdheight_retracked.append(std_height)
		points_retracked.append(len(height_retracked))
		mean_height=np.mean(height_retracked)
		mean_heightmatrix_retracked.append(mean_height)	
		jday_matrix_retracked.append(jday_average_retracked)
	
	
	mean_heightmatrix,jday_matrix,stdheight,points=outlier_rejection(mean_heightmatrix,jday_matrix,stdheight,points)
	mean_heightmatrix_retracked,jday_matrix_retracked,stdheight_retracked,points_retracked=outlier_rejection(mean_heightmatrix_retracked,jday_matrix_retracked,stdheight_retracked,points_retracked)
	
	


	return mean_heightmatrix,jday_matrix,stdheight,points,mean_heightmatrix_retracked,jday_matrix_retracked,stdheight_retracked,points_retracked


rmse_valJ3=[]
rmse_valS6HR=[]
rmse_valS6LR=[]
r2_valJ3=[]
r2_valS6HR=[]
r2_valS6LR=[]
offset_valJ3=[]
offset_valS6HR=[]
offset_valS6LR=[]
offset_valJ3_S6HR=[]
offset_valJ3_S6LR=[]

	
for range_threshold in range(20,65,5):


	heightJ3,jdayJ3,stdheightJ3,pointsJ3,heightJ3_retracked,jdayJ3_retracked,stdheightJ3_retracked,pointsJ3_retracked=height_function("jason3_f_hf",passnr ,min_lat,max_lat,range_threshold)
	heightS6HR,jdayS6HR,stdheightS6HR,pointsS6HR,heightS6HR_retracked,jdayS6HR_retracked,stdheightS6HR_retracked,pointsS6HR_retracked=height_function("sentinel6a_HR_NTC_F_v2_hf",passnr ,min_lat,max_lat,range_threshold)
	heightS6LR,jdayS6LR,stdheightS6LR,pointsS6LR,heightS6LR_retracked,jdayS6LR_retracked,stdheightS6LR_retracked,pointsS6LR_retracked=height_function("sentinel6a_LR_NTC_F_v2_hf",passnr ,min_lat,max_lat,range_threshold)



	time_insitu=[]
	for index,row in data.iterrows():
		cdate=datetime.datetime.strptime(row['Date'],"%Y/%m/%d").date()
		time_insitu.append(cdate)



	def overlap_period(jday):
	
		year=[]
		month=[]
		day=[]
	
		for i in jday:
			t=julianDayDate(i)
			year.append(t['year'])
			month.append(t['month'])
			day.append(t['day'])
		
		ctime=[]
		for i in range(len(year)):
			cdate=datetime.date(year[i],month[i],day[i])
			ctime.append(cdate)
		
		return ctime
	

	timeJ3=overlap_period(jdayJ3)
	timeJ3_retracked=overlap_period(jdayJ3_retracked)
	timeS6HR=overlap_period(jdayS6HR)
	timeS6HR_retracked=overlap_period(jdayS6HR_retracked)
	timeS6LR=overlap_period(jdayS6LR)
	timeS6LR_retracked=overlap_period(jdayS6LR_retracked)


	overlaptime=[]


	s1=set(time_insitu)
	s2=set(timeJ3)
	s3=set(timeS6HR)
	s4=set(timeS6LR)
	s5=set(timeJ3_retracked)
	s6=set(timeS6HR_retracked)
	s7=set(timeS6LR_retracked)

	set1=s1.intersection(s2)
	set2=set1.intersection(s3)
	set3=set2.intersection(s4)
	set4=set3.intersection(s5)
	set5=set4.intersection(s6)
	finalset=set5.intersection(s7)



	overlaptime=list(finalset)
	overlaptime=sorted(overlaptime)

	#print(overlaptime)
	#print(len(overlaptime))


	heighttrue=[]
	heightJ3_overlap=[]
	heightS6HR_overlap=[]
	heightS6LR_overlap=[]
	heightJ3_retracked_overlap=[]
	heightS6HR_retracked_overlap=[]
	heightS6LR_retracked_overlap=[]
	pointsJ3_overlap=[]
	stdheightJ3_overlap=[]
	pointsS6HR_overlap=[]
	stdheightS6HR_overlap=[]
	pointsS6LR_overlap=[]
	stdheightS6LR_overlap=[]
	heightJ3_retracked_overlap=[]
	pointsJ3_retracked_overlap=[]
	stdheightJ3_retracked_overlap=[]
	pointsS6HR_retracked_overlap=[]
	stdheightS6HR_retracked_overlap=[]
	pointsS6LR_retracked_overlap=[]
	stdheightS6LR_retracked_overlap=[]

	for i in overlaptime:
		index=time_insitu.index(i)
		height=data.height[index]
		heighttrue.append(height)
    
		indexJ3=timeJ3.index(i)
		height=heightJ3[indexJ3]
		point=pointsJ3[indexJ3]
		stdev=stdheightJ3[indexJ3]
		heightJ3_overlap.append(height)
		pointsJ3_overlap.append(point)
		stdheightJ3_overlap.append(stdev)
    
		indexS6HR=timeS6HR.index(i)
		height=heightS6HR[indexS6HR]
		point=pointsS6HR[indexS6HR]
		stdev=stdheightS6HR[indexS6HR]
		heightS6HR_overlap.append(height)
		pointsS6HR_overlap.append(point)
		stdheightS6HR_overlap.append(stdev)
    
		indexS6LR=timeS6LR.index(i)
		height=heightS6LR[indexS6LR]
		point=pointsS6LR[indexS6LR]
		stdev=stdheightS6LR[indexS6LR]
		heightS6LR_overlap.append(height)
		pointsS6LR_overlap.append(point)
		stdheightS6LR_overlap.append(stdev)
    
		indexJ3_retracked=timeJ3_retracked.index(i)
		height=heightJ3_retracked[indexJ3_retracked]
		point=pointsJ3_retracked[indexJ3_retracked]
		stdev=stdheightJ3_retracked[indexJ3_retracked]
		heightJ3_retracked_overlap.append(height)
		pointsJ3_retracked_overlap.append(point)
		stdheightJ3_retracked_overlap.append(stdev)
    
		indexS6HR_retracked=timeS6HR_retracked.index(i)
		height=heightS6HR_retracked[indexS6HR_retracked]
		point=pointsS6HR_retracked[indexS6HR_retracked]
		stdev=stdheightS6HR_retracked[indexS6HR_retracked]
		heightS6HR_retracked_overlap.append(height)
		pointsS6HR_retracked_overlap.append(point)
		stdheightS6HR_retracked_overlap.append(stdev)
    
		indexS6LR_retracked=timeS6LR_retracked.index(i)
		height=heightS6LR_retracked[indexS6LR_retracked]
		point=pointsS6LR_retracked[indexS6LR_retracked]
		stdev=stdheightS6LR_retracked[indexS6LR_retracked]
		heightS6LR_retracked_overlap.append(height)
		pointsS6LR_retracked_overlap.append(point)
		stdheightS6LR_retracked_overlap.append(stdev)

	diff=np.median(heightJ3_retracked_overlap)-np.median(heighttrue)
	offset_valJ3.append(diff)
	corr_coeff=np.corrcoef(heightJ3_retracked_overlap,heighttrue)
	r2=corr_coeff[0,1]**2
	r2_valJ3.append(r2)
	rmse=np.sqrt(np.square(np.subtract((heighttrue+diff),heightJ3_retracked_overlap)).mean())
	rmse_valJ3.append(rmse)
	
	diff=np.median(heightS6HR_retracked_overlap)-np.median(heighttrue)
	offset_valS6HR.append(diff)
	corr_coeff=np.corrcoef(heightS6HR_retracked_overlap,heighttrue)
	r2=corr_coeff[0,1]**2
	r2_valS6HR.append(r2)
	rmse=np.sqrt(np.square(np.subtract((heighttrue+diff),heightS6HR_retracked_overlap)).mean())
	rmse_valS6HR.append(rmse)
	
	diff=np.median(heightS6LR_retracked_overlap)-np.median(heighttrue)
	offset_valS6LR.append(diff)
	corr_coeff=np.corrcoef(heightS6LR_retracked_overlap,heighttrue)
	r2=corr_coeff[0,1]**2
	r2_valS6LR.append(r2)
	rmse=np.sqrt(np.square(np.subtract((heighttrue+diff),heightS6LR_retracked_overlap)).mean())
	rmse_valS6LR.append(rmse)
	
	diffJ3=np.median(heightS6HR_retracked_overlap)-np.median(heightJ3_retracked_overlap)
	offset_valJ3_S6HR.append(diffJ3)

	diffJ3=np.median(heightS6LR_retracked_overlap)-np.median(heightJ3_retracked_overlap)
	offset_valJ3_S6LR.append(diffJ3)

diffHR_LR=np.array(offset_valJ3_S6HR)-np.array(offset_valJ3_S6LR)
print(diffHR_LR)
print('mean offset: ',np.mean(diffHR_LR))
print('median offset: ',np.median(diffHR_LR))

plt.subplot(3,1,1)
plt.title('Lake Erie')
plt.plot(range(20,65,5),rmse_valJ3,label='jason 3')
plt.plot(range(20,65,5),rmse_valS6HR,label='sentinel s6hr')
plt.plot(range(20,65,5),rmse_valS6LR,label='sentinel s6Lr')
plt.ylabel('RMSE')
plt.xlabel('threshold')
plt.legend()

plt.subplot(3,1,2)
plt.plot(range(20,65,5),r2_valJ3,label='jason 3')
plt.plot(range(20,65,5),r2_valS6HR,label='sentinel s6hr')
plt.plot(range(20,65,5),r2_valS6LR,label='sentinel s6Lr')
plt.ylabel('R2')
plt.xlabel('threshold')
plt.legend()

plt.subplot(3,1,3)
plt.plot(range(20,65,5),offset_valJ3,label='jason 3')
plt.plot(range(20,65,5),offset_valS6HR,label='sentinel s6hr')
plt.plot(range(20,65,5),offset_valS6LR,label='sentinel s6Lr')
plt.ylabel('Offset')
plt.xlabel('threshold')
plt.legend()
plt.show()



sys.exit(0)

diff=np.median(heightJ3_retracked_overlap)-np.median(heighttrue)
print('Correlation coefficient of Jason-3 retracked with insitu:\n', np.corrcoef(heightJ3_retracked_overlap,heighttrue))
print('offset:',diff)
print('RMSE: ', np.sqrt(np.square(np.subtract((heighttrue+diff),heightJ3_retracked_overlap)).mean()))
print('R2: ', (np.corrcoef(heightJ3_retracked_overlap,heighttrue)[0,1]**2))

diff=np.median(heightJ3_overlap)-np.median(heighttrue)
print('Correlation coefficient of Jason-3 with insitu:\n', np.corrcoef(heightJ3_overlap,heighttrue))
print('offset:',diff)
print('RMSE: ', np.sqrt(np.square(np.subtract((heighttrue+diff),heightJ3_overlap)).mean()))
print('R2: ', (np.corrcoef(heightJ3_overlap,heighttrue)[0,1]**2))

diffJ3=np.median(heightS6HR_overlap)-np.median(heightJ3_overlap)      
diff=np.median(heightS6HR_overlap)-np.median(heighttrue)
print('Correlation coefficient of Sentinel-6aHR with insitu:\n', np.corrcoef(heightS6HR_overlap,heighttrue))
print('offset:',diff)
print('offset from J3: ',diffJ3)
print('RMSE: ', np.sqrt(np.square(np.subtract((heighttrue+diff),heightS6HR_overlap)).mean()))
print('R2: ', (np.corrcoef(heightS6HR_overlap,heighttrue)[0,1]**2))

diffJ3=np.median(heightS6LR_overlap)-np.median(heightJ3_overlap)
diff=np.median(heightS6LR_overlap)-np.median(heighttrue)
print('Correlation coefficient of Sentinel-6a LR with insitu:\n', np.corrcoef(heightS6LR_overlap,heighttrue))
print('offset:',diff)
print('offset from J3: ',diffJ3)
print('RMSE: ', np.sqrt(np.square(np.subtract((heighttrue+diff),heightS6LR_overlap)).mean()))
print('R2: ', (np.corrcoef(heightS6LR_overlap,heighttrue)[0,1]**2))

diffJ3_retracked=np.median(heightS6HR_retracked_overlap)-np.median(heightJ3_retracked_overlap)
diff=np.median(heightS6HR_retracked_overlap)-np.median(heighttrue)
print('Correlation coefficient of Sentinel-6aHR retracked with insitu:\n', np.corrcoef(heightS6HR_retracked_overlap,heighttrue))
print('offset:',diff)
print('offset from J3 retracked: ',diffJ3_retracked)
print('RMSE: ', np.sqrt(np.square(np.subtract((heighttrue+diff),heightS6HR_retracked_overlap)).mean()))
print('R2: ', (np.corrcoef(heightS6HR_retracked_overlap,heighttrue)[0,1]**2))

diffJ3_retracked=np.median(heightS6LR_retracked_overlap)-np.median(heightJ3_retracked_overlap)
diff=np.median(heightS6LR_retracked_overlap)-np.median(heighttrue)
print('Correlation coefficient of Sentinel-6a LR retracked with insitu:\n', np.corrcoef(heightS6LR_retracked_overlap,heighttrue))
print('offset:',diff)
print('offset from J3 retracked: ',diffJ3_retracked)
print('RMSE: ', np.sqrt(np.square(np.subtract((heighttrue+diff),heightS6LR_retracked_overlap)).mean()))
print('R2: ', (np.corrcoef(heightS6LR_retracked_overlap,heighttrue)[0,1]**2))


plt.plot(overlaptime,heighttrue,marker=".",color='r',label='in situ',linewidth=0.75)
plt.plot(overlaptime,heightJ3_overlap,marker=".",color='g',label='J3',linewidth=0.75)
plt.plot(overlaptime,heightS6HR_overlap,marker=".",color='c',label='S6HR',linewidth=0.75)
plt.plot(overlaptime,heightS6LR_overlap,marker=".",label='S6LR',color='y',linewidth=0.75)
plt.plot(overlaptime,heightJ3_retracked_overlap,marker=".",color='g',linestyle='dashed',label='J3_retracked',linewidth=0.75)
plt.plot(overlaptime,heightS6HR_retracked_overlap,marker=".",color='c',linestyle='dashed',label='S6HR_retracked',linewidth=0.75)
plt.plot(overlaptime,heightS6LR_retracked_overlap,marker=".",color='y',linestyle='dashed',label='S6LR_retracked',linewidth=0.75)
plt.xticks(fontsize='7')
plt.xlabel('time (yyyy-mm)',fontsize='15')
plt.ylabel('height (m)',fontsize='15')
plt.title('Lake Michigan, threshold=85')
plt.legend(fontsize=8)
plt.show()


#Jason-3
plt.subplot(2,1,1)
plt.title('Lake Michigan standard deviation, point density plot - Jason-3')
plt.errorbar(overlaptime,heightJ3_overlap,stdheightJ3_overlap,fmt='',capsize=3,elinewidth=0.5)
plt.xticks(fontsize='7')
plt.xlabel('time (yyyy-mm)',fontsize='15')
plt.ylabel('standard deviation in height (m)',fontsize='15')

plt.subplot(2,1,2)
plt.bar(overlaptime,pointsJ3_overlap,width=0.5,align='center')
plt.xticks(fontsize='7')
plt.xlabel('time (yyyy-mm)',fontsize='15')
plt.ylabel('point density',fontsize='15')
plt.show()

#Sentinel-6a HR
plt.subplot(2,1,1)
plt.title('Lake Michigan standard deviation, point density plot - Sentinel-6aHR')
plt.errorbar(overlaptime,heightS6HR_overlap,stdheightS6HR_overlap,fmt='',capsize=3,elinewidth=0.5)
plt.xticks(fontsize='7')
plt.xlabel('time (yyyy-mm)',fontsize='15')
plt.ylabel('standard deviation in height (m)',fontsize='15')

plt.subplot(2,1,2)
plt.bar(overlaptime,pointsS6HR_overlap,width=0.5,align='center')
plt.xticks(fontsize='7')
plt.xlabel('time (yyyy-mm)',fontsize='15')
plt.ylabel('point density',fontsize='15')
plt.show()

#Sentinel-6a LR
plt.subplot(2,1,1)
plt.title('Lake Michigan standard deviation, point density plot - Sentinel-6a LR')
plt.errorbar(overlaptime,heightS6LR_overlap,stdheightS6LR_overlap,fmt='',capsize=3,elinewidth=0.5)
plt.xticks(fontsize='7')
plt.xlabel('time (yyyy-mm)',fontsize='15')
plt.ylabel('standard deviation in height (m)',fontsize='15')

plt.subplot(2,1,2)
plt.bar(overlaptime,pointsS6LR_overlap,width=0.5,align='center')
plt.xticks(fontsize='7')
plt.ylabel('point density',fontsize='15')
plt.xlabel('time (yyyy-mm)',fontsize='15')
plt.show()
