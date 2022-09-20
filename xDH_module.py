#!/bin/env python3
try:
    from  gaussian_manage  import print_Error
    from  gaussian_manage  import print_List
    from  gaussian_manage  import print_String
except:
    from os import getenv
    from os.path import isfile
    import sys
    HomeDir    = getenv('HOME')                                         # STRING, Home DIR
    if isfile('%s/.xdh_modules_path' %HomeDir):                       # Load Private Modules DIR
        with open('%s/.xdh_modules_path'\
                %HomeDir,'r') as tmpf:
            ModuDir=tmpf.readline().strip()                              # STRING, PATH of my modules
            sys.path.append(ModuDir)                                     # Append it into "sys.path"
    else:
        print(('Error in loading \"$HOME/.xdh_modules_path\" \n'+\
            'which contains the absolute path for the relevant py modules'))
        sys.exit(1)
    from  gaussian_manage  import print_Error
    from  gaussian_manage  import print_List
    from  gaussian_manage  import print_String

class xDH:
    '''\
    Handel the xDH@DFA calculations\
    '''
    IFLAG       = 'routine'                                          # Flag : routine or develop
    #List the available R5DFTs and corresponding parameters
    xDHList   = [\
       'xDH',     'DFA','DFT_part','Ex_HF','Ex_LDA','Ex_GGA','Ec_LDA','Ec_GGA','Ec_osPT2','Ec_ssPT2'\
                  ]
    xDHDict   = {\
      'XYG3': [ 'B3LYP','xDH@B3LYP', 0.8033,-0.0140, 0.2107, 0.0000, 0.6789, 0.3211, 0.3211], \
   'XYGJ-OS': [ 'B3LYP','xDH@B3LYP', 0.7731, 0.2269, 0.0000, 0.2309, 0.2754, 0.4364, 0.0000], \
   'revXYG3': [ 'B3LYP','xDH@B3LYP', 0.9196,-0.0222, 0.1026, 0.0000, 0.6059, 0.3941, 0.3941], \
      'XYG5': [ 'B3LYP','xDH@B3LYP', 0.9150, 0.0612, 0.0238, 0.0000, 0.4957, 0.4548, 0.2764], \
      'XYG6': [ 'B3LYP','xDH@B3LYP', 0.9105, 0.1576,-0.0681, 0.1800, 0.2244, 0.4695, 0.2426], \
      'XYG7': [ 'B3LYP','xDH@B3LYP', 0.8971, 0.2055,-0.1408, 0.4056, 0.1159, 0.4052, 0.2589], \
   'XDHPBE0': [ 'PBE1PBE','xDH@PBE0',0.8335, 0.1665, 0.1665, 0.5292, 0.5292, 0.5428, 0.0000]  \
              }
    xDHFamily = {\
            'xDH@B3LYP': ['XYG3','XYGJ-OS','revXYG3','XYG5','XYG6','XYG7'],\
            'xDH@PBE0' : ['XDHPBE0']
            }

    xDHComp   = {\
      'xDH@B3LYP': ['100','205','402'], \
      'xDH@PBE0':  ['100','1000','9'] \
              }

                                                                     # For SP calc.
    SP_ExOvLay  = {\
            'AE':['8/7=1,10=90/1;','9/16=-1/6;','6//8;'],\
            'FC':['8/7=1,10=4/1;','9/16=-1/6;','6//8;']\
                  }
    SP_OptList  = ['IOP(5/33=1)','NoSymm']

    def __init__(self, IOut, GauIO, OptClass, bugctrl=1, gauversion=16, syncinterval=12):
        '''\
        Open the current filename, and initialize some variable belonged to current object\
        ''' 
        from re     import compile
        from os.path    import isfile
        from os         import getcwd

        import gaussian_manage as gaum

        self.GauIO      = GauIO
        self.IOut       = IOut
        self.IPrint     = bugctrl
        self.OptClass   = OptClass
        self.GauVersion = gauversion
        self.SyncInterval = syncinterval
        self.WorkDir    = getcwd().strip()                           # STRING, current DIR 

        self.xDHPara    = []
        FC_Method   = 'FC'                                           # Frozen-core xDH is the default choice

        self.TurnOn     = False                                      # Turn on xDH scheme or Not
        tmpList         = []
        for option in GauIO.OptionList:
            for key in list(xDH.xDHDict.keys()):
                p1 = compile(key.lower())
                tmpOption=option.split('/')
                p1p = p1.match(tmpOption[0].strip().lower())
                if p1p:
                    if len(tmpOption) == 2:
                        tmpList.append(tmpOption[1].strip())
                    elif len(tmpOption) == 1:
                        pass
                    else:
                        print_Error(self.IOut,
                                'Error in determing xDH method: %s' %option)
                    self.TurnOn = True
                    self.xDH    = key
                    self.xDHPara= xDH.xDHDict[key][:]
                    print_String(self.IOut,
                        '"%s" is choosen for the question'
                        % self.xDH,1)
                    if self.IPrint>=2:
                        print_List(self.IOut,
                            self.xDHPara,3,
                            '%s\'s xDHPara(%d) :'
                            % (self.xDH,len(self.xDHPara)))

                    # check if the frozen core algorithm turns on
                    # with the keywords of xyg3(fc) or xyg3(full)
                    p2 = compile('\((?P<fc>\S*)\)')
                    p2p = p2.search(tmpOption[0].strip())
                    #print(p2p.group('fc'))
                    if p2p:
                        if p2p.group('fc').lower()=='full':
                            FC_Method = 'AE'
                    break
            else:
                # check if the frozen core algorithm turns on
                # with the keywords of fc or xyg3
                if option.strip().lower()=='fc':
                    pass
                elif option.strip().lower()=='full':
                    FC_Method = 'AE'
                else:
                    tmpList.append(option)

        GauIO.OptionList = tmpList[:]

        if self.TurnOn:
            if not self.OptClass.Opt:                                # For sp calc.
                GauIO.OptionList.append(self.xDHPara[0])
                for iterm in xDH.SP_OptList:
                    GauIO.OptionList.append(iterm)
                GauIO.OptionList.append('ExtraOverlay')
                GauIO.MoreOptionDict['extraoverlay']    = 1
                tmpExOvLay = xDH.SP_ExOvLay[FC_Method][:]
                GauIO.ExOvList.extend(tmpExOvLay)

                # Modify for "gen" basis set input (V1.5)
                if len(GauIO.RestList)==0:
                    for x in xDH.xDHComp[self.xDHPara[1]]:
                        GauIO.RestList.append('%s' % x)
                else:
                    # Fix for "gen" basis set file input (V1.6)
                    p1      = compile('@')
                    if p1.match(GauIO.RestList[-1]):
                        if self.IPrint >= 2:
                            print_String(self.IOut,
                                'Now check for basis set file',1)
                        if isfile(GauIO.RestList[-1][1:]):
                            tmpf = open(GauIO.RestList[-1][1:],'r')
                            tmpList = tmpf.readlines()
                            tmpf.close()
                            if tmpList[-1].strip()=='':
                                for x in xDH.xDHComp[self.xDHPara[1]]:
                                    GauIO.RestList.append('%s' % x)
                            else:
                                tmpDict=xDH.xDHComp[self.xDHPara[1]][:]
                                tmpDict[0]='\n%s' %tmpDict[0]
                                for x in tmpDict:
                                    GauIO.RestList.append('%s' % x)
                        else:
                            print_Error(self.IOut,'The basis set file %s'
                                    % GauIO.RestList[-1][1:] + ' does not exist')
                    else:
                        tmpDict=xDH.xDHComp[self.xDHPara[1]][:]
                        tmpDict[0]='\n%s' %tmpDict[0]
                        for x in tmpDict:
                            GauIO.RestList.append('%s' % x)
                return
            else:                                                    # For geom. opt. calc.
                pass
        return
    def __del__(self):
        '''\
        To del several class variables\
        '''
        if self.OptClass.Opt and self.TurnOn:
            del self.One
            del self.Two
            del self.MP2
            del self.HF
        return
    def collect_EngyReal(self,EngyPos=None):
        '''\
        Calculate xDH energy\
        '''
        if EngyPos!=None: self.IOut.seek(EngyPos)

        if not self.OptClass.Opt:                                    # For SG Calc::GauIO.EngyReal
            self.collect_xDH_component()
            # now print the energies of the scf procedure and the selected xDH
            tmpEnergy = self.EnoXC
            for x,y in zip(self.ExcList,xDH.xDHDict[self.xDH][2:]):
                tmpEnergy = tmpEnergy + x*y
            print_String(self.IOut,
                'E(%s)%s= %16.8f A.U.    E(%s)%s= %16.8f A.U.'  
                %(self.Method,' '*(9-len(self.Method)),\
                self.EngySCF,\
                self.xDH,' '*(9-len(self.xDH)),\
                tmpEnergy),2)
            if len(xDH.xDHFamily[self.xDHPara[1]])>1:
                print_String(self.IOut,'%s belongs to the family of %s'
                        % (self.xDH, self.xDHPara[1]),1)
            for xMethod in xDH.xDHFamily[self.xDHPara[1]]:
                if xMethod.lower() == self.xDH.lower(): continue
                tmpEnergy = self.EnoXC
                #print(xDH.xDHDict[xMethod][2:])
                for x,y in zip(self.ExcList,xDH.xDHDict[xMethod][2:]):
                    tmpEnergy = tmpEnergy + x*y
                print_String(self.IOut,
                    'E(%s)%s= %16.8f A.U.'  
                    %(xMethod,' '*(9-len(xMethod)),\
                    tmpEnergy),1)
        return self.GauIO.EngyReal
    def cut_Log(self,pos,currdir):
        '''get log information in synchronism,\n
        return current end position'''
        with open('%s/Job_%s.log' \
                %(currdir.value,self.GauIO.JobName), 'r') as tf:
            tf.seek(pos.value)
            tmpString   = tf.read()
            if pos.value != tf.tell():
                try:
                    tmpString   = self.filter_Log(tmpString)
                    self.IOut.write(tmpString)
                    self.IOut.flush()
                except:
                    print('Error happens in cut_Log()')
                pos.value   = tf.tell()
            else:
                pass
            return
        print(('Error happens in cut_Log() to open %s/Job_%s.log' \
                % (currdir.value, self.GauIO.JobName)))
        return
    def run_Job(self,sync=True):
        '''\
        Interface to call the Gaussian package\
        '''
        from os     import remove
        from os     import removedirs
        from os     import listdir
        from os.path import isfile
        from multiprocessing    import Process, Manager
        from time               import sleep                         # sync log file by threads
        if sync==True:
            if not self.OptClass.Opt:
                self.GauIO.form_Inp()
                Pos     = Manager().Value('i',0)
                SyncValue = Manager().Value('c','.script')
                MainJob = Process(target = self.GauIO.run_GauJob,
                        args=(self.GauVersion,SyncValue))
                MainJob.start()
                while MainJob.is_alive():
                    if not isfile('%s/Job_%s.log' \
                            %(SyncValue.value,self.GauIO.JobName)):
                            sleep(self.SyncInterval)
                            continue
                    #print(Pos.value, SyncValue.value)
                    SyncJob = Process(target = self.cut_Log,
                        args=(Pos,SyncValue))
                    SyncJob.start()
                    SyncJob.join()
                    #print('sync log during Gaussian procedure')
                    sleep(self.SyncInterval)
                else:
                    #SyncValue.value   = './'
                    SyncJob = Process(target = self.cut_Log,
                        args=(Pos,SyncValue))
                    SyncJob.start()
                    SyncJob.join()
                    #print('sync log after Gaussian procedure')
                if self.IPrint<2:
                    for tmpFile in listdir('%s' % SyncValue.value):
                        remove('%s/%s' % (SyncValue.value,tmpFile))
                    removedirs('%s' %SyncValue.value)
            else:
                pass
        else:                                                        # Do not sync log file
            if not self.OptClass.Opt:
                self.GauIO.form_Inp()
                self.GauIO.run_GauJob(self.GauVersion)
            else:
                pass
        return
    def collect_xDH_component(self):
        '''Collect each energy terms of DFT method'''
        # LogFile of single point calculation is valid but not geometry optimization job
        from re        import compile
        with open('Job_%s.log' %self.GauIO.JobName,'r') as self.rf:
            self.rf.seek(0)
            lsf = self.rf.read()
        tmpP = 'SCF Done:  E\([RU](?P<scf>\S*)\)\s*=\s*(?P<engy>-?\d+.\d+)\s*A.U. after\s*(?P<snum>\d+) cycles'
        p1   = compile(tmpP)
        p1p  = p1.search(lsf)
        if p1p:
            self.Method = p1p.group('scf')
            self.CycNum = int(p1p.group('snum'))
            self.EngySCF = float(p1p.group('engy'))
        else:
            print_Error(self.IOut,\
            'Error happens in collecting the SCF energies from %s/Job_%s.log' \
                    % (self.WorkDir,self.GauIO.JobName))
        if self.Method == 'B+HF-LYP': self.Method = 'B3LYP'

        #print('debug: %s%s%s' %(self.Method,self.CycNum, self.EngyReal))

        tmpP = 'ENTVJ=\s*(?P<EnoSCF>-?\d*.\d*)\s*Ex=\s*(?P<Ex>-?\d*.\d*)\s*Ec'
        tmpP = tmpP+'=\s*(?P<Ec>-?\d*.\d*)\s*ETotM2e=\s*(?P<ETot>-?\d*.\d*)'
        p2   = compile(tmpP)
        p2p  = p2.findall(lsf)
        self.xDHComp = []
        if p2p:
           for xList in p2p:
               #print(xList)
               self.xDHComp.append([float(x) for x in xList])
        else:
            print_Error(self.IOut,\
            'Error happens in collecting the DFT components from %s/Job_%s.log' \
                    % (self.GauIO.JobName))

        self.xDHPT2 = []
        tmpP = 'alpha-beta\s*T2 =\s*(?P<t2>-?\d*.\d*D[-|+]\d\d)\s*E2=\s*(?P<os>-?\d*.\d*D[-|+]\d\d)'
        p3   = compile(tmpP)
        p3p  = p3.search(lsf)
        if p3p:
            self.xDHPT2.append(float(p3p.group('os').replace('D','E')))
        else:
            print_Error(self.IOut,\
            'Error happens in collecting the alpha-beta osPT2 from %s/Job_%s.log' \
                    % (self.GauIO.JobName))
        tmpP = 'alpha-alpha\s*T2 =\s*(?P<t2>-?\d*.\d*D[-|+]\d\d)\s*E2=\s*(?P<ss>-?\d*.\d*D[-|+]\d\d)'
        p4   = compile(tmpP)
        p4p  = p4.search(lsf)
        if p4p:
            self.xDHPT2.append(float(p4p.group('ss').replace('D','E')))
        else:
            print_Error(self.IOut,\
            'Error happens in collecting the alpha-alpha ssPT2 from %s/Job_%s.log' \
                    % (self.GauIO.JobName))
        tmpP = 'beta-beta\s*T2 =\s*(?P<t2>-?\d*.\d*D[-|+]\d\d)\s*E2=\s*(?P<ss>-?\d*.\d*D[-|+]\d\d)'
        p5   = compile(tmpP)
        p5p  = p5.search(lsf)
        if p5p:
            self.xDHPT2[1]=self.xDHPT2[1]+float(p5p.group('ss').replace('D','E'))
        else:
            print_Error(self.IOut,\
            'Error happens in collecting the beta-beta ssPT2 from %s/Job_%s.log' \
                    % (self.GauIO.JobName))
        #Now re-order the xDH components:
        self.EnoXC = self.xDHComp[0][0]

        if self.GauIO.MoreOptionDict['scrf']!=0:
            tmpP = 'Erf\(P\)=\s*(?P<scrf>-?\d*.\d*)'
            p6   = compile(tmpP)
            p6p  = p6.findall(lsf)
            if p6p:
                self.SCRF = float(p6p[-1])
                print_String(self.IOut,
                    'Solvation energy is considered by %s'
                    % self.GauIO.MoreOptionDict['scrf'],1)
                self.EnoXC =  self.EnoXC + self.SCRF
            else:
                print_Error(self.IOut,\
                'Error happens in collecting the solvation energy from %s/Job_%s.log' \
                        % (self.GauIO.JobName))
                       
        self.ExcList = [self.xDHComp[0][1],       #Ex[HF]
                        self.xDHComp[1][1],       #Ex[LDA]
                        self.xDHComp[2][1],       #Ex[GGA]
                        self.xDHComp[1][2],       #Ec[LDA]
                        self.xDHComp[2][2],       #Ec[GGA]
                        self.xDHPT2[0],           #Ec[osPT2]
                        self.xDHPT2[1],           #Ec[ssPT2]
                       ]
        return
    def filter_Log(self,cons):
        '''filter log information for routine use'''
        import re
        p6 = re.compile(r'(?P<cons1>[^$]*)'                  # Input
                r'\sAdditional\soverlay\scards:\n\s*-*\n\s*8/7=1,10=[0-9]{1,2}/1;\s9/16=-[0-9]/6;\s6//8;\n\s*-*\n'                                                          # seperate 1
                r'(?P<cons2>[^$]*)'                          # SCF+MP2
                r'\s\(Enter\s[^\n]*l608.exe\)\n'             #           # seperate 2
                r'(?P<cons3>[^$]*)'                          # Detail DFT portions
                r'\sLeave\sLink\s*608\sat\s[a-z]{1,3}\s[a-z]{1,3}\s*[0-9]{1,2}\s*[0-9]{1,2}:[0-9]{1,2}:[0-9]{1,2}\s[0-9]{1,4},sMaxMem=\s*[0-9]*\scpu:\s*[0-9]*.[0-9]\n'     # seperate 3
                r'(?P<cons4>[^$]*$)',                        # End
                re.VERBOSE|re.IGNORECASE)

        p5 = re.compile(r'(?P<cons1>[^$]*)'                  # Input
                r'\sAdditional\soverlay\scards:\n\s*-*\n\s*8/7=1,10=[0-9]{1,2}/1;\s9/16=-[0-9]/6;\s6//8;\n\s*-*\n'                                                          # seperate 1
                r'(?P<cons2>[^$]*)'                          # SCF+MP2
                r'\s\(Enter\s[^\n]*l608.exe\)\n'             #           # seperate 2
                r'(?P<cons3>[^$]*$)',                        # Detail DFT portions
                re.VERBOSE|re.IGNORECASE)

        p4 = re.compile(r'(?P<cons2>[^$]*)'                  # SCF+MP2
                r'\s\(Enter\s[^\n]*l608.exe\)\n'             #           # seperate 2
                r'(?P<cons3>[^$]*)'                          # Detail DFT portions
                r'\sLeave\sLink\s*608\sat\s[a-z]{1,3}\s[a-z]{1,3}\s*[0-9]{1,2}\s*[0-9]{1,2}:[0-9]{1,2}:[0-9]{1,2}\s[0-9]{1,4},sMaxMem=\s*[0-9]*\scpu:\s*[0-9]*.[0-9]\n'     # seperate 3
                r'(?P<cons4>[^$]*$)',                        # End
                re.VERBOSE|re.IGNORECASE)

        p3 = re.compile(r'(?P<cons1>[^$]*)'                  # Input
                r'\sAdditional\soverlay\scards:\n\s*-*\n\s*8/7=1,10=[0-9]{1,2}/1;\s9/16=-[0-9]/6;\s6//8;\n\s*-*\n'                                                          # seperate 1
                r'(?P<cons2>[^$]*$)',                        # SCF+MP2
                re.VERBOSE|re.IGNORECASE)

        p2 = re.compile(r'(?P<cons2>[^$]*)'                  # SCF+MP2
                r'\s\(Enter\s[^\n]*l608.exe\)\n'             #           # seperate 2
                r'(?P<cons3>[^$]*$)',                        # Detail DFT portions
                re.VERBOSE|re.IGNORECASE)

        p1 = re.compile(r'(?P<cons3>[^$]*)'                  # Detail DFT portions
                r'\sLeave\sLink\s*608\sat\s[a-z]{1,3}\s[a-z]{1,3}\s*[0-9]{1,2}\s*[0-9]{1,2}:[0-9]{1,2}:[0-9]{1,2}\s[0-9]{1,4},sMaxMem=\s*[0-9]*\scpu:\s*[0-9]*.[0-9]\n'     # seperate 3
                r'(?P<cons4>[^$]*$)',                        # End
                re.VERBOSE|re.IGNORECASE)

        if p1.match(cons):                                           # End (- "cons3")
            tmpGroup    = p1.match(cons)
            cons        = tmpGroup.group('cons4')
        elif p2.match(cons):                                         # (SCF+MP2) (- "cons3")
            tmpGroup    = p2.match(cons)
            cons        = tmpGroup.group('cons2')
        elif p3.match(cons):                                         # Input+(SCF+MP2)
            tmpGroup    = p3.match(cons)
            cons        =\
                    tmpGroup.group('cons1') + tmpGroup.group('cons2')
        elif p4.match(cons):                                         # (SCF+MP2)+End (- "cons3")
            tmpGroup    = p4.match(cons)
            cons        = \
                    tmpGroup.group('cons2') + tmpGroup.group('cons4')
        elif p5.match(cons):                                         # Input+(SCF+MP2) (- "cons3")
            tmpGroup    = p5.match(cons)
            cons        =\
                    tmpGroup.group('cons1') + tmpGroup.group('cons2')
        elif p6.match(cons):                                         # Input+(SCF+MP2)+End
            tmpGroup    = p6.match(cons)
            cons        = tmpGroup.group('cons1') + \
                    tmpGroup.group('cons2') + tmpGroup.group('cons4')
        return cons

