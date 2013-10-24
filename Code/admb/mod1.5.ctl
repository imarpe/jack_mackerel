mod1.dat  
Model_1.5
# Selectivity sharing vector (number_fisheries + number_surveys)  
#Fsh 1 Fsh_2 Fsh_3 Fsh_4 Srv_1 Srv_2 Srv_3 Srv_4 Srv_5 Srv_6 Srv_7 Srv_8  
#N_Chile_Fshry CS_Chile_Fishery Peruvian International Chile_AcousCS Chile_AcousN Chilean_CPUE DEPM Acoustic_Peru Peru_CPUE Chinese_CPUE EU_CPUE #  
# 2 3 4 1 2 3 4 5 6 7 8 9  
1 1 1 1 2 1 1 2 1 1 1 1 1  
1 2 3 4 1 1 2 4 3 3 4 4 4  
#Sr_type 
2  
#AgeError 
0  
#Retro 
0  
#Steepness 
0.8 300 -6  
#SigmaR 
1.0 15 -4  
#yrs_sr 
1970 2011 
#Linf
74.4	0.1	-4																																									
#K
0.16	0.1	-4																																									
#Lo_Len
18	0.1	-4																																									
#Sigma_len
0.09 0.1	-4																																									
#Natural_Mortality 
0.23 0.05 -4   
# NEW npars_mage
0
# NEW Mage_in

# phase_Mage
-5
#Phase_Random_walk_M 
-4
#Nyrs_Random_walk_M 
0
#Random_walk_M_yrs blank if nyrs==0

#Random_walk_M_sigmas blank if nyrs==0

# 76.464 # 70.8  
#catchability 
0.7632  0.0255  5.e-5  0.4370  0.059  0.0105  0.0002  0.0470  0.0081  
1.2  1.2  1.2  1.2  1.2  1.2  1.2  1.2  1.2  
3  5  3  3  3  4  4  4  4  
#q_power                    
1  1  1  1  1  1  1  1  1  
1.2  1.2  1.2  1.2  1.2  1.2  1.2  1.2  1.2  
-1  -1  -1  -1  -1  -1  -1  -1  -1  
#Random_walk_q_phases                    
1  -1  1  -1  1  -1  -1  -1  -1  
#Nyrs_Random_walk_q
1  0  1  0  2  0  0  0  0  
#Random_walk_q_yrs blank if nyrs==0
2002 
2011
1994 1997
#Random_walk_q_sigmas blank if nyrs==0
2.0
2.0
2.0 2.0
#q_agemin                    
2  2  2  2  2  2  2  2  2  
#q_agemax                    
10  10  10  10  10  10  10  10  10  
#junk                    
0.05  
#n_proj_yrs                    
10  
#---------------------------------------------------------
# Fishery 1 N Chile  
1  #selectivity type
9  #n_sel_ages
2  #phase sel
1  #curvature penalty
1  #Dome-shape penalty
# Years of selectivity change Fishery 1 N Chile  
30																														
	1984	1985	1986	1987	1988	1989	1990	1991	1992	1993	1994	1995	1996	1997	1998	1999	2000	2001	2002	2003	2004	2005	2006	2007	2008	2009	2010	2011	2012	2013
	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5
# Initial values for coefficitients at each change (one for every change plus 1)  
# 2 3 4 5 6 7 8 9 10 11 12  
0.2 0.7 1 1 1 1 1 1 1 1 1 1  
#---------------------------------------------------------
# Fishery 2, Central South Chile  
1																														
10																														
3																														
1																														
25																														
#	Years	of	selectivity	change	Fishery	2,	Central	South	Chile																					
30																														
	1984	1985	1986	1987	1988	1989	1990	1991	1992	1993	1994	1995	1996	1997	1998	1999	2000	2001	2002	2003	2004	2005	2006	2007	2008	2009	2010	2011	2012	2013
	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5
#	Initial	values	for	coefficitients	at	each	change	(one	for	every	change	plus	1)																	
#	2	3	4	5	6	7	8	9	10	11	12																			
0.2	0.7	1	1	1	1	1	1	1	1	1	1																			
#---------------------------------------------------------
# Fishery 3 Peru  
1  
7  
4  
1  
12.5  
# Years of selectivity change Fishery 3 Peru  
32 
1981	1982	1983	1984	1985	1986	1987	1988	1989	1990	1991	1992	1993	1994	1995	1996	1997	1998	1999	2000	2001	2002	2003	2004	2005	2006	2007	2008	2009	2010	2011	2012
0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 
# Initial values for coefficitients at each change (one for every change plus 1)  
# 2 3 4 5 6 7 8 9 10 11 12  
0.2 0.7 1 1 1 1 1 1 1 1 1 1  
#---------------------------------------------------------
# Fishery 4 International  
1  
9  
3  
1  
12.5  
# Years of selectivity change Fishery 4 International  
23
1980	1981	1982	1983	1984	1985	1986	1987	1988	1989	1990	1991	2000	2001	2002	2003	2006	2007	2008	2009	2010	2011	2012
0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 
# Initial values for coefficitients at each change (one for every change plus 1)  
# 2 3 4 5 6 7 8 9 10 11 12  
0.2 0.7 1 1 1 1 1 1 1 1 1 1
#---------------------------------------------------------  
# Index number 1 AcousCS  
1  
10  
2  
0.25  
100  
1
2005
0.7
#0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.7 0 0 0 0 0 0
0.3 1 1 1 1 1 1 1 1 1 1 1
#---------------------------------------------------------  
# Index number 2 Acous_N  
1  
10  
-5  
0.25  
100  
0 #0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0.3 1 1 1 1 1 1 1 1 1 1 1 
#--------------------------------------------------------- 
# Index number 3 Chile_CPUE  
1  
10  
-5  
0.25  
100  
0 #0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0.8 1 1 1 1 1 1 1 1 1 1 1  
#---------------------------------------------------------
# Index number 4 DEPM  
1  
10  
3  
0.25  
100 # 1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010
1
2003
0.7
#0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.7 0 0 0 0 0 0 0 0
0.04 0.43 0.93 1 1 1 1 1 1 1 1 1  
#---------------------------------------------------------
# Index number 5 Acoustic_Peru  
1  
10  
-5  
0.25  
100  
0 #0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0.04 0.43 0.93 1 1 1 1 1 1 1 1 1  
#---------------------------------------------------------
# Index number 6 Peru_CPUE  
1  
10  
-5  
0.25  
100  
0 #0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0.04 0.43 0.93 1 1 1 1 1 1 1 1 1  
#---------------------------------------------------------
# Index number 7 EU_CPUE  
1  
10  
-5  
0.25  
100  
0 #0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0.04 0.43 0.93 1 1 1 1 1 1 1 1 1  
#---------------------------------------------------------
# Index number 8 EU_CPUE  
1  
10  
-5  
0.25  
100  
0 #0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0.04 0.43 0.93 1 1 1 1 1 1 1 1 1  
#---------------------------------------------------------
# Index number 8 EU_CPUE  
1  
10  
-5  
0.25  
100  
0 #0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0.04 0.43 0.93 1 1 1 1 1 1 1 1 1  
#---------------------------------------------------------
#Test  
123456789  

