from numpy import loadtxt
import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import os, sys
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, FuncFormatter, NullFormatter)
import math
import re


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]



src_dir       = '1665'
parent_dir    = os.path.basename(os.getcwd())



for file in sorted(os.listdir(src_dir)):
    src_path = os.path.join(src_dir, file)
    src_name = file[:-10]
    print(src_path)
    print(src_name)

    mylines = []
    #for data in open(src_path + '/recon_1665.py', 'rt').readlines():
    #	if data.__contains__('x1_lim, x2_lim '):
    #		data_1 = data
    #		# with open("data_1.py", "w") as py_file_1:
    #		with open(src_path + "/data_1.txt", "w") as py_file_1:
    #			py_file_1.write(data_1)
    #			py_file_1.close()
    #			sys.path.insert(1, '/home/chikaedu/Desktop/Fitting_mainlines/polarisation_spectra/1665/'+src_path)
 	   	#	print(data_1)
    #	if data.__contains__('tmj_tk, tmn_tk'):
    #		data_2 = data
    		# with open("data_2.py", "w") as py_file_2:
    #		with open(src_path + "/data_2.txt", "w") as py_file_2:
    #			py_file_2.write(data_2)
    #			py_file_2.close()
    #			sys.path.insert(1, '/home/chikaedu/Desktop/Fitting_mainlines/polarisation_spectra/1665/'+src_path)


		data1 = loadtxt('/home/chikaedu/Desktop/Fitting_mainlines/polarisation_spectra/1665/'+src_path + '/i.dat')
		data2 = loadtxt('/home/chikaedu/Desktop/Fitting_mainlines/polarisation_spectra/1665/'+src_path + '/q.dat')
		data3 = loadtxt('/home/chikaedu/Desktop/Fitting_mainlines/polarisation_spectra/1665/'+src_path + '/u.dat')
		data4 = loadtxt('/home/chikaedu/Desktop/Fitting_mainlines/polarisation_spectra/1665/'+src_path + '/v.dat')

		ai, bi = np.array(data1[:, 0]), np.array(data1[:, 1])
		aq, bq = np.array(data2[:, 0]), np.array(data2[:, 1])
		au, bu = np.array(data3[:, 0]), np.array(data3[:, 1])
		av, bv = np.array(data4[:, 0]), np.array(data4[:, 1])

		a_file = open('/home/chikaedu/Desktop/Fitting_mainlines/polarisation_spectra/1665/'+src_path + '/data_1.txt', "r")
		list_of_lists = [(line.strip()).split() for line in a_file]
		a_file.close()
		digit_1 = list(re.findall(r"[-+]?\d*\.\d+|\d+", str(list_of_lists)))

		b_file = open('/home/chikaedu/Desktop/Fitting_mainlines/polarisation_spectra/1665/'+src_path + '/data_2.txt', "r")
		list_of_lists = [(line.strip()).split() for line in b_file]
		a_file.close()
		digit_2 = list(re.findall(r"[-+]?\d*\.\d+|\d+", str(list_of_lists)))


		print(digit_1[2])
		print(digit_1[3])
		print(digit_2[0])
		print(digit_2[1])


		x1_lim, x2_lim   =  float(digit_1[2]), float(digit_1[3])


		ul = find_nearest(ai, x2_lim)
		ll = find_nearest(ai, x1_lim)

		upper = int(np.where(ai == ul)[0])
		lower = int(np.where(ai == ll)[0])
		print('upper', upper)
		print('lower', lower)
		bi_idx    = bi[lower:upper]
		bv_idx    = bv[lower:upper]

		print('max_I', bi_idx.max())
		print('max_V', bv_idx.max())



		###############################################################################################
		
		if 0 < bi_idx.max() < 0.5:
			ty1_lim, ty2_lim =  bv_idx.min()-0.1, bi_idx.max()+0.25		# y-axis limits for top panel
			tmj_tk, tmn_tk   =  0.5, 0.25
		elif 0.5 < bi_idx.max() < 1:
			ty1_lim, ty2_lim =  bv_idx.min()-0.1, bi_idx.max()+0.25		# y-axis limits for top panel
			tmj_tk, tmn_tk   =  0.5, 0.25			
		elif 1 < bi_idx.max() < 5:
			ty1_lim, ty2_lim =  bv_idx.min()-1, bi_idx.max()+1		# y-axis limits for top panel
			tmj_tk, tmn_tk   =  1.0, 0.5
		elif 5 < bi_idx.max() < 10:
			ty1_lim, ty2_lim =  bv_idx.min()-2.5, bi_idx.max()+2.5		# y-axis limits for top panel
			tmj_tk, tmn_tk   =  2.0, 1.0
		elif 10 < bi_idx.max() < 20:
			ty1_lim, ty2_lim =  bv_idx.min()-5, bi_idx.max()+5		# y-axis limits for top panel
			tmj_tk, tmn_tk   =  5.0, 2.5
		elif 20 < bi_idx.max() < 50:
			ty1_lim, ty2_lim =  bv_idx.min()-10, bi_idx.max()+10		# y-axis limits for top panel
			tmj_tk, tmn_tk   =  10.0, 5.0
		elif 50 < bi_idx.max() < 100:
			ty1_lim, ty2_lim =  bv_idx.min()-20, bi_idx.max()+20		# y-axis limits for top panel
			tmj_tk, tmn_tk   =  20.0, 10.0
		elif 100 < bi_idx.max() < 200:
			ty1_lim, ty2_lim =  bv_idx.min()-20, bi_idx.max()+20		# y-axis limits for top panel
			tmj_tk, tmn_tk   =  50.0, 25.0
		else:
			ty1_lim, ty2_lim =  bv_idx.min()-50, bi_idx.max()+50		# y-axis limits for top panel
			tmj_tk, tmn_tk   =  100.0, 50.0

		
		# tmj_tk, tmn_tk   =  float(digit_2[0]), float(digit_2[1])
		###############################################################################################

		fts  = 15
		lbsz = 14
		lgds = 16

		plt.rc('font', weight='bold')
		plt.rc('text', usetex=True)
		plt.rc('xtick', labelsize=12)
		plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
		mpl.rcParams['axes.linewidth'] = 2.
		legend_properties = {'weight':'bold'}

		f, axarr = plt.subplots(figsize=(6,3))

		axarr.plot(ai, bi,   label=r'\textbf{Stokes-I}', color='k', linestyle='-')
		axarr.plot(aq, bq, label=r'\textbf{Stokes-Q}', color='B', linestyle='-')
		axarr.plot(au, bu, label=r'\textbf{Stokes-U}', color='G', linestyle='-')
		axarr.plot(av, bv, label=r'\textbf{Stokes-V}', color='R', linestyle='-')
		plt.xlabel(r'\textbf{LSR Velocity (km $s^{-1}$)}',fontsize=14)
		plt.ylabel(r'\textbf{Flux density (Jy)}', fontsize=14)
		axarr.xaxis.set_major_locator(MultipleLocator(5))
		axarr.xaxis.set_minor_locator(MultipleLocator(2.5))
		axarr.yaxis.set_major_locator(MultipleLocator(tmj_tk))
		axarr.yaxis.set_minor_locator(MultipleLocator(tmn_tk))
		axarr.tick_params(axis='x', labelsize=lbsz, pad=8)
		axarr.tick_params(axis='y', labelsize=lbsz)
		axarr.tick_params(which='both', width=1)
		axarr.tick_params(which='major', length=8, direction='in')
		axarr.tick_params(which='minor', length=4, direction='in')
		axarr.yaxis.set_major_formatter(FormatStrFormatter(r'$\mathbf{%.1f}$'))
		axarr.xaxis.set_major_formatter(FormatStrFormatter(r'$\mathbf{%.1f}$'))
		axarr.set_ylim(ty1_lim, ty2_lim)
		axarr.set_ylim(ty1_lim, ty2_lim)
		axarr.set_xlim(x1_lim, x2_lim)
		axarr.set_title(r'$\bf{'+src_name+'}$' + ' ' + r'$\bf{' + '(' + src_dir + ' ' + '~MHz'  + ')}$', size=14)
		plt.legend(loc='upper right', fontsize=12.)
		src_name_new 	= src_name.replace('.', '_')
		plt.savefig('IQUV_images/'+src_name_new + '_IQUV_' +  src_dir  + '.pdf', format = 'pdf', dpi=plt.gcf().dpi, bbox_inches = 'tight')


sys.exit()


