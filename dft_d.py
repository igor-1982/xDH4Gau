#!/bin/env python3
# Module	:: dft_d
# Authors	:: Igor Ying Zhang, and Xin Xu
# Purpose	:: To handle dispersion contribution for molecular systems
__Version__ = 'V2.1(20091101)'
# Revise    :: 2009-11-01
#
# History   :: V1.0(20090907) :: 1) build new class "DFTD" to handle DFT+D scheme, cooperating
#                                   with Gaussian package.
#                                2) add some five functions for dispersion calculation :
#                                        get_Dist,    get_Damp,    get_DampDeri, 
#                                        get_ForcSix, get_HessSix.
#              V2.0(20091021) :: 1) seperate the class "DFTD" into two parts:
#                                   1.1) The first part is the interface to Gaussian package
#                                        which is moved to the modules of "gaussian_manage"
#                                   1.2) The second part is the dispersion calculation which
#                                        is grouped into new class named "DispGrim" in this
#                                        module, since the six-order coefficients adopted here
#                                        is suggested by Grimme et. al..
#                                2) those functions builded in previous version for dispersion
#                                   calculation are grouped in to the new class "DispGrim"
#              V2.1(20091101) :: 1) Function of evaluating general no-damped six-order term has
#                                   been added into the class "DispGrim" utilizing Grimme's six-
#                                   order coefficients.
#                                2) Utilizing the relationship between six-, and twelve-order 
#                                   coefficients in Lennard-Jones potential, extended 
#                                   twelve-order coefficients are derived based on Grimme's
#                                   six-order coefficients and Grimme's vdW radios.
#                                   Consequently, function of evaluating damped and no-damped
#                                   twelve-order terms has been added into the class "DispGrim"
#                                   based on developed twelve-order coefficients.
#              V2.2(20100717) :: 1) Add four input arguments for D parameters fitting:
#                                   1: self.R0Para; 2: self.C6Para; 3: self.C12Para; 4: self.DPara
#                                   


class DispGrim:
    '''\
    This class calculate the dispersion contribution proposed by Grimme et. al.\n\
    Six-order parameters come from the papers followed :\n\
    1) Grimme... Vol.27 , 1787, J. Comp. Chem.\n\
    2) Martin... Vol.113, 8434, J. Phys. Chem. A\n\
    Twelve-order parameters are generated the correspond six-order one and R0\n\
    E_well  = C6Para / (4 * pow(R0,6))\n\
    C12Para = 4 * E_well * pow(R0,12)\
    '''
    Br2Ang  = 0.5291772083 
    Ang2Nm  = 0.1000000000
    Kc2Kj   = 4.1840000000
    Au2Kc   = 627.50950000
    Eau6Con = 10.0**(-3)/(Kc2Kj *Au2Kc) * Ang2Nm**(-6)               # j2Kj/Au2Kj*Ang2Nm**(-6)
    Eau12Con= 10.0**(-3)/(Kc2Kj *Au2Kc) * Ang2Nm**(-12)              # j2Kj/Au2Kj*Ang2Nm**(-12)
    Fau6Con = Eau6Con * Br2Ang
    Hau6Con = Fau6Con * Br2Ang

    AtDict  = {\
    'x':0  ,
    'h':1  , 'he':2  ,
   'li':3  , 'be':4  ,  'b':5  ,  'c':6  ,  'n':7  ,  'o':8  ,  'f':9  , 'ne':10 ,
   'na':11 , 'mg':12 , 'al':13 , 'si':14 ,  'p':15 ,  's':16 , 'cl':17 , 'ar':18 , 
    'k':19 , 'ca':20 , 'sc':21 , 'ti':22 ,  'v':23 , 'cr':24 , 'mn':25 , 'fe':26 ,
   'co':27 , 'ni':28 , 'cu':29 , 'zn':30 , 'ga':31 , 'ge':32 , 'as':33 , 'se':34 , 'br':35 , 'kr': 36 ,
   'rb':37 , 'sr':38 ,  'y':39 , 'zr':40 , 'nb':41 , 'mo':42 , 'tc':43 , 'ru':44 ,
   'rh':45 , 'pd':46 , 'ag':47 , 'cd':48 , 'in':49 , 'sn':50 , 'sb':51 , 'te':52 ,  'i':53 , 'xe': 54\
                  }
    C6Para  = [\
     0.00 ,
     0.14 ,    0.08 ,
     1.61 ,    1.61 ,    3.13 ,    1.75 ,    1.23 ,    0.70 ,    0.75 ,    0.63 ,
     5.71 ,    5.71 ,   10.79 ,    9.23 ,    7.84 ,    5.57 ,    5.07 ,    4.61 ,
    10.08 ,   10.08 ,   10.08 ,   10.08 ,   10.08 ,   10.08 ,   10.08 ,   10.08 ,
    10.08 ,   10.08 ,   10.08 ,   10.08 ,   16.99 ,   17.10 ,   16.37 ,   12.64 ,   12.47 ,  12.01 ,
    24.67 ,   24.67 ,   24.67 ,   24.67 ,   24.67 ,   24.67 ,   24.67 ,   24.67 ,
    24.67 ,   24.67 ,   24.67 ,   24.67 ,   37.32 ,   38.71 ,   38.44 ,   31.74 ,   31.50 ,  29.99\
                  ]                                                  # C6Para (in J/mol*nm^6)
    #   C12Para = [\
    #    0.00 ,
    # 1.41E-07, 8.59E-08,
    # 5.08E-07, 1.25E-05, 3.36E-05, 1.64E-05, 9.14E-06, 4.09E-06, 3.41E-06, 2.32E-06,
    # 1.28E-05, 3.68E-05, 2.09E-04, 2.36E-04, 1.93E-04, 1.27E-04, 9.83E-05, 7.59E-05,
    # 1.08E-04, 1.03E-04\
    #               ]                                                 # C12Para (in J/mol*nm^12)
                                                                     # 4e[(s/r)^-12-(s/r)^-6]
    C12Para = [\
     0.00 ,
  0.704E-07, 4.30E-08,
  2.54E-07, 0.627E-05, 1.68E-05, 0.82E-05, 4.57E-06, 2.044E-06, 1.704E-06, 1.16E-06,
  0.64E-05, 1.84E-05, 1.046E-04, 1.18E-04, 0.963E-04, 0.633E-04, 4.914E-05, 3.795E-05,
  0.54E-04, 0.517E-04\
                  ]                                                  # C12Para (in J/mol*nm^12)
                                                                     # e[(s/r)^-12-2(s/r)^-6]
    R0Para  = [\
     1.000,
     1.001,    1.012,
     0.825,    1.408,    1.485,    1.452,    1.397,    1.342,    1.287,    1.243,
     1.144,    1.364,    1.639,    1.716,    1.705,    1.683,    1.639,    1.595, 
     1.485,    1.474,    1.562,    1.562,    1.562,    1.562,    1.562,    1.562,
     1.562,    1.562,    1.562,    1.562,    1.650,    1.727,    1.760,    1.771,    1.749,    1.727,
     1.628,    1.606,    1.639,    1.639,    1.639,    1.639,    1.639,    1.639,
     1.639,    1.639,    1.639,    1.639,    1.672,    1.804,    1.881,    1.892,    1.892,    1.881\
                  ]                                                  # R0Para (in Anstrom)

    DPara   = 20.0                                                   # Damp parameter

    def __init__(self, iout, clist, ian, iop=0, bugctrl=0,\
        c6para=None, c12para=None, r0para=None, dpara=None):
        '''\
        Initialize parameters in this class\
        '''
        from my_io import print_Error
        self.IOut   = iout
        self.CList  = clist
        self.IAn    = ian
        self.IPrint = bugctrl
        self.ICtrl  = iop                                            # ICtrl = 0 : Grimme Disp.
                                                                     #         1 : 6-order term
                                                                     #         2 : 12-order term
                                                                     #         3 : "1" + "2"
                                                                     #         4 : damped 12-order
        self.EngyReal = 0.0
        self.ForcList = []
        self.HessList = []

        if r0para==None:                                             # Now loading dispersion parameters
            self.R0Para = DispGrim.R0Para[:]
        else:
            self.R0Para = r0para[:]
        if c6para==None:
            self.C6Para = DispGrim.C6Para[:]
        else:
            self.C6Para = c6para[:]
        if c12para==None:
            self.C12Para = DispGrim.C12Para[:]
        else:
            self.C12Para = c12para[:]
        if dpara==None:
            self.DPara = DispGrim.DPara
        else:
            try:
                self.DPara = dpara
            except TypeError:
                print_Error(self.IOut,'dpara should be "float" but not "List"')
        return
    def get_Dist(self,Atom1,Atom2):
        '''\
        Calculate distance between Atom1 and Atom2\
        '''
        from math import sqrt
        dist=0.0
        for k in range(3): 
            dist=dist+pow((Atom1[k]-Atom2[k]),2.0)
        dis=sqrt(dist)
        return dis
    
    def get_Damp(self,d,Rij,Rr):
        '''\
        Calculate the dumping in Grimme\'s dispersion\
        '''
        from math import exp
        damp=pow(1.0+exp(-d*(Rij/Rr-1.0)),-1.0)
        return damp
    
    def get_DampDeri(self,Atom1,Atom2,d,Rij,Rr):
        '''\
        Calculate the damping derivate in Atom1 contributed by Atom2 in Grimme\'s dispersion\
        '''
        from math import pow
        from math import exp
        Factor  = -d / (Rij*Rr) * pow(1.0+exp(-d*(Rij/Rr-1.0)),-2.0)
        XTerm   = Factor * (Atom1[0]-Atom2[0])
        YTerm   = Factor * (Atom1[1]-Atom2[1])
        ZTerm   = Factor * (Atom1[2]-Atom2[2])
        return (XTerm, YTerm, ZTerm)
    def get_EngyReal(self):
        '''\
        Calculate the dispersion energy\
        '''
        from my_io import print_String
        from my_io import print_List
        from my_io import print_List_free
        from math  import sqrt
        if len(self.IAn)==1:
            self.EngyReal   = 0.0
            print_String(self.IOut,'E(%s)%s= %16.8f A.U. '
                %('Disp',' '*5,self.EngyReal),2)
            return self.EngyReal
        if self.ICtrl==0:                                            # for six-order damped disp.
            self.EngyReal   = 0.0
            ResuList    =[\
                ['Atom1','Atom2','Rij','Rr','Cij','Damp','Disp_G']\
                ]
            Scal    = -1.0 * DispGrim.Eau6Con
            for i in range(len(self.IAn)):
                for j in range(i):
                    Rij	= self.get_Dist(self.CList[i],self.CList[j])
                    Rr  = self.R0Para[self.IAn[i]]+\
                        self.R0Para[self.IAn[j]]
                    Cij = sqrt(self.C6Para[self.IAn[i]]*\
                        self.C6Para[self.IAn[j]])
                    Dmp = self.get_Damp(self.DPara,Rij,Rr)
                    Eny = Scal * Dmp * Cij * pow(Rij,-6)
                    self.EngyReal = self.EngyReal + Eny
                    ResuList.append(\
                        [self.IAn[i],i+1,self.IAn[j],j+1,\
                        Rij,Rr,Cij,Dmp,Eny])
            self.EngyReal	= self.EngyReal
            if self.IPrint>=2:
                FormList    = [\
                        '%4d%4dth%4d%4dth%8.4f%8.4f%8.4f%8.4f%8.4f',
                        '%10s%10s%8s%8s%8s%8s%8s']
                print_List_free(self.IOut,ResuList,2,FormList,
                    'Grimme\'s six-order damped dispersion')
            print_String(self.IOut,'E(%s)%s= %16.8f A.U. '
                %('Disp_G',' '*3,self.EngyReal),2)
        elif self.ICtrl==1:                                          # for pure six-order disp.
            self.EngyReal   = 0.0
            ResuList    =[\
                ['Atom1','Atom2','Rij','Cij','Disp_6']\
                ]
            Scal    = -1.0 * DispGrim.Eau6Con
            for i in range(len(self.IAn)):
                for j in range(i):
                    Rij	= self.get_Dist(self.CList[i],self.CList[j])
                    Cij = sqrt(self.C6Para[self.IAn[i]]*\
                        self.C6Para[self.IAn[j]])
                    Eny = Scal * Cij * pow(Rij,-6)
                    self.EngyReal = self.EngyReal + Eny
                    ResuList.append(\
                        [self.IAn[i],i+1,self.IAn[j],j+1,\
                        Rij,Cij,Eny])
            self.EngyReal	= self.EngyReal
            if self.IPrint>=2:
                FormList    = ['%4d%4dth%4d%4dth%16.8f%16.8f%16.8f',
                        '%10s%10s%16s%16s%16s']
                print_List_free(self.IOut,ResuList,2,FormList,
                    'Pure six-order dispersion')
            print_String(self.IOut,'E(%s)%s= %16.8f A.U. '
                %('Disp_6',' '*3,self.EngyReal),2)
        elif self.ICtrl==2:                                          # for pure twelve-order disp.
            self.EngyReal   = 0.0
            ResuList    =[\
                ['Atom1','Atom2','Rij','C12ij','Disp_12']\
                ]
            Scal    = self.Eau12Con
            for i in range(len(self.IAn)):
                for j in range(i):
                    Rij	= self.get_Dist(self.CList[i],self.CList[j])
                    Cij = sqrt(self.C12Para[self.IAn[i]]*\
                        self.C12Para[self.IAn[j]])
                    Eny = Scal * Cij * pow(Rij,-12)
                    self.EngyReal = self.EngyReal + Eny
                    ResuList.append(\
                        [self.IAn[i],i+1,self.IAn[j],j+1,\
                        Rij,Cij,Eny])
            self.EngyReal	= self.EngyReal
            if self.IPrint>=2:
                print_String(self.IOut,'Scal is %16.8f' %Scal ,2)
                FormList    = ['%4d%4dth%4d%4dth%16.8f%16.8f%16.8f',
                        '%10s%10s%16s%16s%16s']
                print_List_free(self.IOut,ResuList,2,FormList,
                    'Pure twelve-order dispersion')
            print_String(self.IOut,'E(%s)%s= %16.8f A.U. '
                %('Disp_12',' '*2,self.EngyReal),2)
        elif self.ICtrl==3:                                          # for 6 + 12 disp.
            self.EngyReal   = 0.0
            ResuList    =[\
                ['Atom1','Atom2','Rij','Cij','Disp.']\
                ]
            ScalSix     = -1.0 * DispGrim.Eau6Con
            ScalTwel    = DispGrim.Eau12Con
            for i in range(len(self.IAn)):
                for j in range(i):
                    Rij	= self.get_Dist(self.CList[i],self.CList[j])
                    C6ij    = sqrt(self.C6Para[self.IAn[i]]*\
                        self.C6Para[self.IAn[j]])
                    C12ij   = sqrt(self.C12Para[self.IAn[i]]*\
                        self.C12Para[self.IAn[j]])
                    Eny = ScalSix * C6ij * pow(Rij,-6) + \
                        ScalTwel * C12ij * pow(Rij,-12)
                    self.EngyReal = self.EngyReal + Eny
                    ResuList.append(\
                        [self.IAn[i],i+1,self.IAn[j],j+1,\
                        Rij,C6ij,C12ij,Eny])
            if self.IPrint>=2:
                FormList    = ['%4d%4dth%4d%4dth%8.4f%8.4f%8.4f%8.4f',
                        '%10s%10s%8s%8s%8s%8s']
                print_List_free(self.IOut,ResuList,2,FormList,
                    '6+12-order dispersion')
            print_String(self.IOut,'E(%s)%s= %16.8f A.U. '
                %('Disp',' '*5,self.EngyReal),2)
        elif self.ICtrl==4:                                          # for damped 12 Disp.
            self.EngyReal   = 0.0
            ResuList    =[\
                ['Atom1','Atom2','Rij','Rr','Cij','Damp','Disp_12s']\
                ]
            Scal    = -1.0 * DispGrim.Eau12Con
            for i in range(len(self.IAn)):
                for j in range(i):
                    Rij	= self.get_Dist(self.CList[i],self.CList[j])
                    Rr  = self.R0Para[self.IAn[i]]+\
                        self.R0Para[self.IAn[j]]
                    Cij = sqrt(self.C12Para[self.IAn[i]]*\
                        self.C12Para[self.IAn[j]])
                    Dmp = self.get_Damp(self.DPara/2,Rij,Rr)
                    Eny = Dmp * Cij * pow(Rij,-12.0)
                    self.EngyReal = self.EngyReal + Eny
                    ResuList.append(\
                        [self.IAn[i],i+1,self.IAn[j],j+1,\
                        Rij,Rr,Cij,Dmp,Eny])
            self.EngyReal	= Scal * self.EngyReal
            print_String(self.IOut,'E(%s)%s= %16.8f A.U. '
                %('Disp_12s',' '*1,self.EngyReal),2)
        return self.EngyReal
# Mark by Igor for uncomplete get_ForcList
    def get_ForcList(self):
        '''\
        Calculate the corresponding force for DFT+D methods\n\
        Contribution of atom j in atom i is :\n\
        Fij_x	= - Scal * (6.0 * (X_i-X_j) ) /|r_i-r_j|^(-8)\n\
        Fij_y	= - Scal * (6.0 * (Y_i-Y_j) ) /|r_i-r_j|^(-8)\n\
        Fij_z	= - Scal * (6.0 * (Z_i-Z_j) ) /|r_i-r_j|^(-8)\
        '''
        from math import sqrt

        from my_io import print_List
        from my_io import print_String
        Scal		= -1.0 * DFTD.FauCon * self.SixPara
        self.ForcList	= []
        if self.IPrint>=2:
            print_String(self.IOut, 
            '%10s%10s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s'
            %('Atom1','Atom2','Rij','Rr','Cij','Damp',
            'Fx_Disp','Fx_FDamp','Fy_Disp','Fy_FDamp',
            'Fz_Disp','Fz_FDamp'),2)
        if self.GauIO.NAtom==1:
            NTT             = self.GauIO.NAtom*3
            self.ForcList   = [0.0]*NTT
            if self.IPrint>=2:
                print_List(self.IOut,self.ForcList,4,
                    'CM Dispersion Cartesian Gradient:')
            return
        for i in range(self.GauIO.NAtom):
            TmpForceX	= 0.0
            TmpForceY	= 0.0
            TmpForceZ	= 0.0
            for j in range(self.GauIO.NAtom):
                if j!=i:
                    Rij	= get_Dist(self.CList[i],self.CList[j])
                    Rr	= DFTD.R0Para[self.IAn[i]]+\
                          DFTD.R0Para[self.IAn[j]]
                    Cijg= DFTD.C6Para[self.IAn[i]]*\
                          DFTD.C6Para[self.IAn[j]]
                    Cij	= sqrt(Cijg)
                    Damp= get_Damp(self.d,Rij,Rr)
                    XTerm1, YTerm1, ZTerm1	=\
                        get_Force(self.CList[i],self.CList[j],Rij)
                    XTerm	= XTerm1 * Damp
                    YTerm	= YTerm1 * Damp
                    ZTerm	= ZTerm1 * Damp
                    XTerm2, YTerm2, ZTerm2	=\
                        get_DampD(self.CList[i],self.CList[j],\
                        self.d,Rij,Rr)
                    XTermD	= XTerm2 * pow(Rij,-6)
                    YTermD	= YTerm2 * pow(Rij,-6)
                    ZTermD	= ZTerm2 * pow(Rij,-6)

                    TmpForceX		+= XTerm * Cij 
                    TmpForceY 		+= YTerm * Cij 
                    TmpForceZ 		+= ZTerm * Cij 

                    TmpForceX		+= XTermD * Cij 
                    TmpForceY 		+= YTermD * Cij 
                    TmpForceZ 		+= ZTermD * Cij 
                    if self.IPrint>=2:
                        self.IOut.write(\
 '%7s%3dth%5s%3dth%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n' \
 %(self.AtLabel[i],i+1,self.AtLabel[j],j+1,\
 Rij,Rr,Cij,Damp,XTerm,XTermD,YTerm,YTermD,ZTerm,ZTermD))

            TmpForceX   = Scal * TmpForceX
            TmpForceY   = Scal * TmpForceY
            TmpForceZ   = Scal * TmpForceZ
            self.ForcList.append(TmpForceX)
            self.ForcList.append(TmpForceY)
            self.ForcList.append(TmpForceZ)
        if self.IPrint>=2:
            print_List(self.IOut,self.ForcList,4,
                'CM Dispersion Cartesian Gradient:')
        return
# Mark by Igor for uncomplete get_HessList
    def get_HessList(self):

        return
# Mark by Igor for uncomplete get_DampHess
    def get_DampHess(self,Atom1,Atom2,d,Rij,Rr):
        '''\
        Calculate the damping Hessian in Atom1 contributed by Atom2 in Grimme\'s dispersion\
        '''
        return (XX, YY, ZZ, XY, XZ, YZ)
# Mark by Igor for uncomplete get_HessSix
    def get_HessSix(self,Atom1,Atom2):
        '''\
        Calculate six-order dispersion Hessian of Atom1 contributed by Atom2\
        '''
        return (XX, YY, ZZ, XY, XZ, YZ)
# Mark by Igor for uncomplete get_ForcSix
    def get_ForcSix(self,Atom1,Atom2,Rij):
        '''\
        Calculate six-order dispersion force of Atom1 contributed by Atom2\
        '''
        from math import pow
        Factor	= -6.0 / pow(Rij,8.0)
        XTerm	= Factor * (Atom1[0]-Atom2[0])
        YTerm	= Factor * (Atom1[1]-Atom2[1])
        ZTerm	= Factor * (Atom1[2]-Atom2[2])
        r         eturn (XTerm, YTerm, ZTerm)

class DFTD:
    '''\
    This class handles DFT+D calculation cooperating with Gaussian package and "dft_d" module.\n\
    NOTE :: "DFTD" should be loaded after modules of  "GauIO" and "OptHandle"\n\
    Scale parameter for six-order terms "s6" following are originated from\n\
        1) Grimme... Vol.27 , 1787, J. Comp. Chem.\n\
        2) Martin... Vol.113, 8434, J. Phys. Chem. A\n\
    \n\
    disp_g  is the pure dispersion of Grimme\n\
    disp_6  is the unscaled Grimme's dispersion\n\
    disp_12 is the unscaled tweleve-order term extended by Grimme's s6 and LG form\n\
    disp_xy\n\
    disp\
    '''
    GrimDict    = {\
  'b3lyp+d' : [1.05 ,  'b3lyp'] ,      'b971+d' : [0.65 ,   'b971'] , 
   'b972+d' : [1.05 ,   'b972'] ,       'bmk+d' : [0.65 ,    'bmk'] , 
   'blyp+d' : [1.20 ,   'blyp'] ,     
  'm062x+d' : [0.06 ,  'm062x'] ,     
    'pbe+d' : [0.75 , 'pbepbe'] ,    'pbepbe+d' : [0.75 , 'pbepbe'] , 
   'pbe0+d' : [0.70 ,   'pbe0'] ,   
   'tpss+d' : [1.00 ,   'tpss'] ,   
  'x3lyp+d' : [0.85 ,  'x3lyp']\
                }
    DispDict    = {
   'disp_g' : [1.00 , 'disp_g'] , 'dispersion_g': [1.00,  'disp_g'] ,
   'disp_6' : [1.00 , 'disp_6'] , 'dispersion_6': [1.00,  'disp_6'] ,
  'disp_12' : [1.00 ,'disp_12'] ,'dispersion_12': [1.00, 'disp_12'] ,
 'disp_12s' : [1.00 ,'disp_12s'] ,'dispersion_12s': [1.00, 'disp_12s'] ,
     'disp' : [1.00 ,   'disp'] ,   'dispersion': [1.00,    'disp']\
                  }

    def __init__(self, iout, GauIO, OptClass, bugctrl=0,\
        c6para=None, c12para=None, r0para=None, dpara=None):
        '''\
        Initialize parameters\
        '''
        self.IOut       = iout                                       # Set the logout file string
        self.GauIO      = GauIO                                      # Loading GauIO class
        self.OptClass   = OptClass                                   # Loading OptHandle class
        self.DispClass  = None                                       # DispClass
 #
        self.IPrint     = bugctrl                                    # INTEGER, control debugging
        self.AtLabel    = self.GauIO.AtLabel                         # List of AtLabel
        self.CList      = self.GauIO.CList                           # List of Geom. Info.
        self.IAn        = self.GauIO.IAn                             # List of IAn
        self.Method     = ''                                         # STRING, name of DFT+D
        self.SixPara    = 1.0                                        # REAL, Parameter of disp.
 #
        self.EngyReal   = 0.0                                        # REAL, dispersion energy
        self.ForcList   = []                                         # List, dispersion force
        self.HessList   = []                                         # List, dispersion hessian
        self.TurnOn     = False                                      # LOGIC, turn on or off DFT+D
        self.PureDisp   = False                                      # LOGIC, pure disp. calc.
 #
        for option in self.GauIO.OptionList:                         # Filter DFT+D method 
            for key in sorted(DFTD.GrimDict.keys()):
                TmpList=option.strip().split('/')
                if TmpList[0].lower()==key:
                    self.TurnOn   = True
                    self.Method   = key
                    self.SixPara=DFTD.GrimDict[key][0]
                    if len(TmpList)==1:
                        self.GauIO.OptionList.append(\
                            DFTD.GrimDict[key][1])
                    elif len(TmpList)==2:
                        self.GauIO.OptionList.append(\
                            '/'.join([DFTD.GrimDict[key][1],\
                            TmpList[1]]))
                    else:
                        print_Error(self.IOut,
                            'Error happens in DFT+D method '+\
                            'determination \"DFTD.__init__\"')
                    self.GauIO.OptionList.remove(option)
                    break
        if not self.TurnOn:                                          # Filter pure disp. calc.
            for option in self.GauIO.OptionList:
                for key in sorted(DFTD.DispDict.keys()):
                    TmpList = option.strip().split('/')
                    if TmpList[0].lower()==key:
                        self.TurnOn = True
                        self.Method = key
                        self.SixPara    = DFTD.DispDict[key][0]
                        if len(TmpList)==1:
                            self.GauIO.OptionList.append(\
                                DFTD.DispDict[key][1])
                        elif len(TmpList)==2:
                            self.GauIO.OptionList.append(\
                                '/'.join([DFTD.DispDict[key][1],\
                                TmpList[1]]))
                        else:
                            print_Error(self.IOut,
                                'Error happens in Disp. method '+\
                                'determination \"DFTD.__init__\"')
                        self.GauIO.OptionList.remove(option)

        self.OptClass.OptionList    = self.GauIO.OptionList[:]       # Sync OptionList between 
                                                                     #  "OptClass" and "GauIO"
        if self.TurnOn:
            if (self.Method == 'disp_g') or\
            (self.Method == 'dispersion_g'):
                self.PureDisp   = True
                self.DispClass  =\
                    DispGrim(self.IOut, GauIO.CList, GauIO.IAn, 0,\
                    self.IPrint,\
                    c6para, c12para, r0para, dpara)
                if self.IPrint>=1:
                    print_String(self.IOut, 'Disp_G is required',1)
                return
            if (self.Method == 'disp_6') or\
            (self.Method == 'dispersion_6'):
                self.PureDisp   = True
                self.DispClass  =\
                    DispGrim(self.IOut, GauIO.CList, GauIO.IAn, 1,\
                    self.IPrint,\
                    c6para, c12para, r0para, dpara)
                if self.IPrint>=1:
                    print_String(self.IOut, 'Disp_6 is required',1)
                return
            if (self.Method == 'disp_12') or\
            (self.Method == 'dispersion_12'):
                self.PureDisp   = True
                self.DispClass  =\
                    DispGrim(self.IOut, GauIO.CList, GauIO.IAn, 2,\
                    self.IPrint,\
                    c6para, c12para, r0para, dpara)
                if self.IPrint>=1:
                    print_String(self.IOut, 'Disp_12 is required',1)
                return
            if (self.Method == 'disp') or\
            (self.Method == 'dispersion'):
                self.PureDisp   = True
                self.DispClass       =\
                    DispGrim(self.IOut, GauIO.CList, GauIO.IAn, 3,\
                    self.IPrint,\
                    c6para, c12para, r0para, dpara)
                if self.IPrint>=1:
                    print_String(self.IOut,'Disp. is required',1)
                return
            if (self.Method == 'disp_12s') or\
            (self.Method == 'dispersion_12s'):
                self.PureDisp   = True
                self.DispClass       =\
                    DispGrim(self.IOut, GauIO.CList, GauIO.IAn, 4,\
                    self.IPrint,\
                    c6para, c12para, r0para, dpara)
                if self.IPrint>=1:
                    print_String(self.IOut,'Disp_12s is required',1)
                return

            for option in self.OptClass.OptionList:
                TmpList=option.strip().split('/')
                if TmpList[0].lower()==self.Method:
                    if len(TmpList)==1:
                        self.OptClass.OptionList.append(\
                            DFTD.GrimDict[self.Method][1])
                    elif len(TmpList)==2:
                        self.OptClass.OptionList.append(\
                            DFTD.GrimDict[self.Method][1])
                        self.OptClass.OptionList.append(\
                            TmpList[1])
                    self.OptClass.OptionList.remove(option)
                    break
            self.DispClass       =\
                DispGrim(self.IOut, GauIO.CList, GauIO.IAn, 0,\
                self.IPrint)
            if self.IPrint>=1:
                print_String(self.IOut,
                    '%s is employed with s_6 = %3.2f'
                    % (self.Method.upper(),self.SixPara),2)
        else:
            self.DispClass       =\
                DispGrim(self.IOut, GauIO.CList, GauIO.IAn, 0,\
                self.IPrint)
            if self.IPrint>=2:
                print_String(self.IOut,
                    'None DFT+D scheme is required',1)
        return
    def __del__(self):

        return
    def get_EngyReal(self,ICtrl=0):
        '''\
        Get DFT+D or Dispersion energy\n\
        NOTE :: ICtrl   = 0 : default, the same as 1\n\
                        = 1 : to run_GauJob() before QM energy collecting\n\
                          2 : bypass run_GauJob(), to get QM energy from given Chkfile\n\
        '''
        if self.PureDisp:
            self.EngyReal   = self.DispClass.get_EngyReal()
        elif ICtrl==0 or ICtrl==1:                                   # run_GauJob() for QM results
            self.EngyReal   = self.DispClass.get_EngyReal() *\
                self.SixPara
            self.GauIO.form_Inp()
            self.GauIO.run_GauJob()
            Result  = ChkHandle(self.IOut,self.GauIO,self.IPrint)
            self.EngyReal   = my_plus(self.EngyReal,
                Result.collect_EngyReal())
            del Result
        elif ICtrl==2:                                               # bypass run_GauJob()
            self.EngyReal   = self.DispClass.get_EngyReal() *\
                self.SixPara
            Result  = ChkHandle(self.IOut,self.GauIO,self.IPrint)
            self.EngyReal   = my_plus(self.EngyReal,
                Result.collect_EngyReal())
            del Result
        return self.EngyReal
    def get_ForcList(self,ICtrl=0):
        '''\
        Get DFT+D or Dispersion force\n\
        NOTE :: ICtrl   = 0 : default, the same as 1\n\
                        = 1 : to run_GauJob() before QM force collecting\n\
                          2 : bypass run_GauJob(), to get QM force from given Chkfile\n\
        '''
        if self.PureDisp:
            self.ForcList   =\
                [x for x in self.DispClass.get_ForcList()]
        elif ICtrl==0 or ICtrl==1:
            self.ForcList   =\
                [x *self.SixPara\
                for x in self.DispClass.get_ForcList()]
            self.GauIO.form_Inp()
            self.GauIO.run_GauJob()
            Result  = ChkHandle(self.IOut,self.GauIO,self.IPrint)
            for i in range(len(self.ForcList)):
                self.ForcList[i] = my_plus(self.ForcList[i],
                    Result[i])
            del Result
        elif ICtrl==2:
            self.ForcList   =\
                [x *self.SixPara\
                for x in self.DispClass.get_ForcList()]
            Result  = ChkHandle(self.IOut,self.GauIO,self.IPrint)
            for i in range(len(self.ForcList)):
                self.ForcList[i] = my_plus(self.ForcList[i],
                    Result[i])
            del Result
        return self.ForcList
            
