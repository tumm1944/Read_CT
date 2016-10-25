#!/usr/bin/env python
# Filename: Read_qchem_CT.py
#
import os, sys
import numpy as np
import re
import glob
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.stats as scs
from scipy.stats import kurtosis
import statsmodels.sandbox.distributions.extras as extras
from math import sqrt
from dihedrals import *

def call_read_CT(argv):
	FileCounter = 0
	myPath = os.path.dirname(os.path.realpath(__file__))
#	print myPath
	if len(argv) < 3:
        	print "\nUsage: Read_qchem_CT.py <filename extension> #(atoms in molecule1)  #ROOTS  <#state of interest>" 
        	sys.exit(1)
	else:
        	inFile = sys.argv[1]
		AM1 = sys.argv[2]
		ROOTS = sys.argv[3]
		State_I= sys.argv[4]
		SearchFile = "*." + inFile
		
		FileCounter = len(glob.glob1(myPath,SearchFile))
		if FileCounter < 1:
			print ' No files of this extension '
			sys.exit(1)
		else:
#	if not os.path.isfile('%s.custom2' %inFile):
#            print ' %s.custom2 DOES NOT EXIST !!!' %inFile
#	if not os.path.isfile('%s.custom3' %inFile):
#            print ' %s.custom3 DOES NOT EXIST !!!' %inFile
			dplot(inFile, AM1, ROOTS, State_I)

def dplot(extname, Atoms1, NRT, CT_state):
	SString4 = "Excited state"
        SString5 = "Excited State"
        SString3 = "CIS_N_Roots"
        SString1 = "$molecule"
        SString2 = "$end"


	myPath = os.path.dirname(os.path.abspath(__file__))
#	print myPath

        SearchFile = myPath + "/*." + extname
#	print SearchFile
	Filenames = glob.glob(SearchFile)
#	print Filenames
#	Filenames = glob.glob('/scratch/rajesh/Anthracene_CIF_30fs/CIF_AllNeigh_0x/*.dat')
	Fcounter  = len(Filenames) 
	print "# of files", Fcounter
	NRT = int(NRT)
	CT_state = int(CT_state)
	MolWt = {'C': 12.011, 'S': 32.011, 'H': 1.008}

	CT1 = [[0 for x in range(NRT)] for y in range(Fcounter)] # CT state 
        CTE = [[0 for x in range(NRT)] for y in range(Fcounter)] # CT state energy
        ESE = [[0 for x in range(NRT)] for y in range(Fcounter)] # Excited state energy
        TCT = [[0 for x in range(NRT)] for y in range(Fcounter)] # total charge transfer	
	GCT = [0 for x in range(Fcounter)]			 # Ground charge transfer
	FC1 = [[0 for x in range(NRT)] for y in range(Fcounter)] # Frag1 charge transfer
	FC2 = [[0 for x in range(NRT)] for y in range(Fcounter)] # Frag2 charge transfer
	FC3 = [[0 for x in range(NRT)] for y in range(Fcounter)] # Frag3 charge transfer
	FC4 = [[0 for x in range(NRT)] for y in range(Fcounter)] # Frag4 charge transfer
	FC5 = [[0 for x in range(NRT)] for y in range(Fcounter)] # Frag5 charge transfer
	DF1 = [[0 for x in range(5)] for y in range(Fcounter)] # Distances for fragements
	BD1 = [[0 for x in range(5)] for y in range(Fcounter)] # Distances for C-C bonds
	DH1 = [[0 for x in range(5)] for y in range(Fcounter)] # Dihedral between fragments
	Arr_CTnum = []
	Arr_Ener = []
	Arr_SCF = []
	Arr_TCT = []
	index_dist_min = []
	




#	print Fcounter
#	print Filenames

# 0= Volume, 1=faces, 2=facearea, 3= effective dia, 4=aphericity, 5=Q (36*Pi*V^3/S^2)
#rf=open('%s.custom2' %systemname, 'r')
#qf=open('%s.custom3' %systemname, 'r')
#   while True:
	for j in range(Fcounter):
#		print Filenames[j]
		with open(Filenames[j], 'r') as inF:
		    XC = []
		    YC = []
		    ZC = []
		    MC = []
		    count = 0
		    end_count=0
		    ex_count=0
		    ct_count=0
		    while True:
			line = inF.readline()
			if not line: 
				break
			count=count+1
        		if SString1 in line:
#				print line
				inF.readline()
                                for k in range(int(Atoms1)+94):
                                        line1 = inF.readline()
                                        line1 = line1.strip().split()
                                        XC.append(float(line1[-3])*MolWt[line1[0]] )
					YC.append(float(line1[-2])*MolWt[line1[0]] )
					ZC.append(float(line1[-1])*MolWt[line1[0]] )
					MC.append(MolWt[line1[0]] )
#				line = line.strip().split()
			if SString2 in line:
#				if end_count == 0:
#					atom_count=count-atom_count1-3
				end_count=end_count+1
			if SString3 in line:
				line = line.strip().split()
				rootsc = line[-1:]
			if SString4 in line:
				ex_count=ex_count+1
				line = line.strip(': ').split()
				CT1[j][ex_count-1] = ex_count
#				print Filenames[j], ex_count, j
				CTE[j][ex_count-1] = float(line[-1])
				line1 = inF.readline()
				line1 = line1.strip().split()
				ESE[j][ex_count-1] = float(line1[-1])
			if SString5 in line:	
				line = line.strip().split()
				ct_count=ct_count+1
				for k in range(3):
					inF.readline()
				for k in range(int(Atoms1)):
					line1 = inF.readline()
					line1 = line1.strip().split()
#		                        TCT[j][ct_count-1] = float(line1[-1]) + TCT[j][ct_count-1]
					if k < 10 or k==50 :
						FC1[j][ct_count-1] = float(line1[-1]) + FC1[j][ct_count-1]
                                        if k > 9 and k < 20:
                                           	FC2[j][ct_count-1] = float(line1[-1]) + FC2[j][ct_count-1]    
                                        if k > 19 and k < 30:
	                                        FC3[j][ct_count-1] = float(line1[-1]) + FC3[j][ct_count-1]
                                        if k > 29 and k < 40:
        					FC4[j][ct_count-1] = float(line1[-1]) + FC4[j][ct_count-1]                                     
                                        if k > 39 and k < 50: 
                                            	FC5[j][ct_count-1] = float(line1[-1]) + FC5[j][ct_count-1]
					if k == 51:
						FC5[j][ct_count-1] = float(line1[-1]) + FC5[j][ct_count-1]
				TCT[j][ct_count-1] = FC1[j][ct_count-1] + FC2[j][ct_count-1] + FC3[j][ct_count-1] + FC4[j][ct_count-1] + FC5[j][ct_count-1]
					
                        if "Ground-State" in line:
                                line = line.strip().split()
                                for k in range(3):
                                        inF.readline()
                                for k in range(int(Atoms1)):
                                        line1 = inF.readline()
                                        line1 = line1.strip().split()
                                        GCT[j] = float(line1[-1]) + GCT[j]
					
		XC = np.array(XC)
		YC = np.array(YC)
		ZC = np.array(ZC)
		MC = np.array(MC)
		C60X = sum(XC[52:111])/sum(MC[52:111])
		C60Y = sum(YC[52:111])/sum(MC[52:111])
		C60Z = sum(ZC[52:111])/sum(MC[52:111]) 
		DF1[j][0]= sqrt( (sum(XC[0:9])/sum(MC[0:9]) - C60X)**2 + (sum(YC[0:9])/sum(MC[0:9]) - C60Y)**2 + (sum(ZC[0:9])/sum(MC[0:9]) - C60Z)**2 )
		DF1[j][1]= sqrt( (sum(XC[10:19])/sum(MC[10:19]) - C60X)**2 + (sum(YC[10:19])/sum(MC[10:19]) - C60Y)**2 + (sum(ZC[10:19])/sum(MC[10:19]) - C60Z)**2 )
		DF1[j][2]= sqrt( (sum(XC[20:29])/sum(MC[20:29]) - C60X)**2 + (sum(YC[20:29])/sum(MC[20:29]) - C60Y)**2 + (sum(ZC[20:29])/sum(MC[20:29]) - C60Z)**2 )
		DF1[j][3]= sqrt( (sum(XC[30:39])/sum(MC[30:39]) - C60X)**2 + (sum(YC[30:39])/sum(MC[30:39]) - C60Y)**2 + (sum(ZC[30:39])/sum(MC[30:39]) - C60Z)**2 )
		DF1[j][4]= sqrt( (sum(XC[40:49])/sum(MC[40:49]) - C60X)**2 + (sum(YC[40:49])/sum(MC[40:49]) - C60Y)**2 + (sum(ZC[40:49])/sum(MC[40:49]) - C60Z)**2 )
		index_dist_min.append( np.argmin(DF1[j][:]))
		BD1[j][0]= distance_2p( np.array([XC[1],YC[1],ZC[1]]/MC[1]),np.array([XC[10],YC[10],ZC[10]])/MC[10])
		BD1[j][4]= distance_2p( np.array([XC[2],YC[2],ZC[2]]/MC[2]),np.array([XC[11],YC[11],ZC[11]])/MC[11])
		if  BD1[j][0] <  BD1[j][4]:
			q1 =  np.array([XC[2],YC[2],ZC[2]]/MC[2])
			q2 =  np.array([XC[1],YC[1],ZC[1]]/MC[1])
			q3 = np.array([XC[10],YC[10],ZC[10]]/MC[10])
			q4 = np.array([XC[12],YC[12],ZC[12]]/MC[12])
			DH1[j][0] = new_dihedral(q1,q2,q3,q4) 

                        q1 =  np.array([XC[12],YC[12],ZC[12]]/MC[12])
                        q2 =  np.array([XC[11],YC[11],ZC[11]]/MC[11])
                        q3 =  np.array([XC[20],YC[20],ZC[20]]/MC[20])
                        q4 =  np.array([XC[22],YC[22],ZC[22]]/MC[22])
                        DH1[j][1] = new_dihedral(q1,q2,q3,q4)

                        q1 = np.array([XC[22],YC[22],ZC[22]]/MC[22])
                        q2 = np.array([XC[21],YC[21],ZC[21]]/MC[21])
                        q3 = np.array([XC[30],YC[30],ZC[30]]/MC[30])
                        q4 = np.array([XC[32],YC[32],ZC[32]]/MC[32])
                        DH1[j][2] = new_dihedral(q1,q2,q3,q4)

                        q1 = np.array([XC[32],YC[32],ZC[32]]/MC[32])
                        q2 = np.array([XC[31],YC[31],ZC[31]]/MC[31])
                        q3 = np.array([XC[40],YC[40],ZC[40]]/MC[40])
                        q4 = np.array([XC[42],YC[42],ZC[42]]/MC[42])
                        DH1[j][3] = new_dihedral(q1,q2,q3,q4)
		
			BD1[j][1]= distance_2p( np.array([XC[11],YC[11],ZC[11]]/MC[11]),np.array([XC[20],YC[20],ZC[20]]) /MC[20])
			BD1[j][2]= distance_2p( np.array([XC[21],YC[21],ZC[21]]/MC[21]),np.array([XC[30],YC[30],ZC[30]]) /MC[30])
			BD1[j][3]= distance_2p( np.array([XC[31],YC[31],ZC[31]]/MC[31]),np.array([XC[40],YC[40],ZC[40]]) /MC[40])
		else :
                        q1 =  np.array([XC[3],YC[3],ZC[2]]/MC[3])
                        q2 =  np.array([XC[2],YC[2],ZC[1]]/MC[2])
                        q3 = np.array([XC[11],YC[11],ZC[10]]/MC[11])
                        q4 = np.array([XC[13],YC[13],ZC[12]]/MC[13])
                        DH1[j][0] = new_dihedral(q1,q2,q3,q4)

                        q1 =  np.array([XC[13],YC[13],ZC[13]]/MC[13])
                        q2 =  np.array([XC[12],YC[12],ZC[12]]/MC[12])
                        q3 =  np.array([XC[21],YC[21],ZC[21]]/MC[21])
                        q4 =  np.array([XC[23],YC[23],ZC[23]]/MC[23])
                        DH1[j][1] = new_dihedral(q1,q2,q3,q4)

                        q1 = np.array([XC[23],YC[23],ZC[23]]/MC[23])
                        q2 = np.array([XC[22],YC[22],ZC[22]]/MC[22])
                        q3 = np.array([XC[31],YC[31],ZC[31]]/MC[31])
                        q4 = np.array([XC[33],YC[33],ZC[33]]/MC[33])
                        DH1[j][2] = new_dihedral(q1,q2,q3,q4)

                        q1 = np.array([XC[33],YC[33],ZC[33]]/MC[33])
                        q2 = np.array([XC[32],YC[32],ZC[32]]/MC[32])
                        q3 = np.array([XC[41],YC[41],ZC[41]]/MC[41])
                        q4 = np.array([XC[43],YC[43],ZC[43]]/MC[43])
                        DH1[j][3] = new_dihedral(q1,q2,q3,q4)

		
			BD1[j][0] =  BD1[j][4]
			BD1[j][1]= distance_2p( np.array([XC[12],YC[12],ZC[12]]/MC[12]),np.array([XC[21],YC[21],ZC[21]]) /MC[21])
                        BD1[j][2]= distance_2p( np.array([XC[22],YC[22],ZC[22]]/MC[22]),np.array([XC[31],YC[31],ZC[31]]) /MC[31])
                        BD1[j][3]= distance_2p( np.array([XC[32],YC[32],ZC[32]]/MC[32]),np.array([XC[41],YC[41],ZC[41]]) /MC[41])



	for i in range(Fcounter):
		print_true=0
		for j in range(NRT):
			if TCT[i][j] > 0.5:					# comment this for CT1
				 print_true = 1  + print_true			# comment this for CT1
			if  TCT[i][j] > 0.5 and print_true == CT_state and DF1[i][index_dist_min[i]] < 10 :
					name = Filenames[i].split('/')
#					if CTE[i][j] < 0.6 :
#					print "CT State =", CT1[i][j], "CT Energy (eV)= ", CTE[i][j], "Total Energy (Hartrees) = ", ESE[i][j], "Total Charge Transfer (e) =", TCT[i][j], name[-1]
					loc_rat = 0
					loc_rat = ((FC1[i][j]/TCT[i][j])/0.2)**2
					loc_rat = loc_rat + ((FC2[i][j]/TCT[i][j])/0.2)**2
                                        loc_rat = loc_rat + ((FC3[i][j]/TCT[i][j])/0.2)**2
                                        loc_rat = loc_rat + ((FC4[i][j]/TCT[i][j])/0.2)**2
                                        loc_rat = loc_rat + ((FC5[i][j]/TCT[i][j])/0.2)**2
#					loc_rat = (loc_rat - 5.0)/20.0
#					print CT1[i][j], "CTE (eV)= ", CTE[i][j],FC1[i][j],FC2[i][j],FC3[i][j],FC4[i][j],FC5[i][j],TCT[i][j], loc_rat, name[-1], DF1[i][index_dist_min[i]], index_dist_min[i]
#					FC_all=np.array([FC1[i][j],FC2[i][j],FC3[i][j],FC4[i][j],FC5[i][j]])
#					print  CT1[i][j],CTE[i][j],FC1[i][j],FC2[i][j],FC3[i][j],FC4[i][j],FC5[i][j],index_dist_min[i],BD1[i][0],BD1[i][1],BD1[i][2], BD1[i][3],  name[-1]
					print  CT1[i][j],CTE[i][j],DH1[i][0],DH1[i][1],DH1[i][2], DH1[i][3], index_dist_min[i],FC1[i][j],FC2[i][j],FC3[i][j],FC4[i][j],FC5[i][j],BD1[i][0],BD1[i][1],BD1[i][2], BD1[i][3],  name[-1] 
					Arr_CTnum.append(CT1[i][j])
					Arr_Ener.append(CTE[i][j])
					Arr_SCF.append(ESE[i][j])
					Arr_TCT.append(TCT[i][j])



	Arr_CTnum = np.array(Arr_CTnum)
	Arr_Ener = np.array(Arr_Ener)
	Arr_SCF = np.array(Arr_SCF)
	Arr_TCT = np.array(Arr_TCT)
	GCT = np.array(GCT)
	index_dist_min=np.array(index_dist_min)
	
	print "       ", "CT State #", "CT_State_Energy (eV)", "Total Charge Transferred", "Ground-State Charge" , "(charge trasnferrred are for P3HT)"
	print "Mean", np.mean(Arr_CTnum), np.mean(Arr_Ener), np.mean(Arr_TCT), np.mean(GCT)
	print "Std. Dev.", np.std(Arr_CTnum), np.std(Arr_Ener), np.std(Arr_TCT), np.std(GCT)


if __name__ == '__main__':
    call_read_CT(sys.argv[1:])

#####  End  #####
