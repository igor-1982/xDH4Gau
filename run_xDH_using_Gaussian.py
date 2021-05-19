#!/bin/env python3
#Filename:  run_xDH_using_Gaussian.py
#Author  :  Igor Ying Zhang, Xin Xu 
#Purpose :  Run Rung 5th DFT method calculation, based on Gaussian package
#           Utilize the gaussian input file as it's input file
#
#Usage   :  run_xDH_using_Gaussian.py GauInpFile
#
import os
import sys

__version__    = '1.0(2021)'
iprint         = 1
__gaussian__   = 16
syncinterval  = 6

WorkDir    = os.getcwd().strip()                                 # STRING, current DIR 
HomeDir    = os.getenv('HOME')                                   # STRING, Home DIR
if os.path.isfile('%s/.xdh_modules_path' %HomeDir):            # Load Private Modules DIR
    tmpf    = open('%s/.xdh_modules_path'\
        %HomeDir,'r')
    ModuDir=tmpf.readline().strip()                              # STRING, PATH of my modules
    sys.path.append(ModuDir)                                     # Append it into "sys.path"
    tmpf.close()
else:
    print(('Error for loading \"$HOME/.xdh_modules_path\" \n'+\
       'which contains absolute path of personal python modules'))
    sys.exit(1)
if os.path.isfile('%s/version.txt' %ModuDir):                   # Load Private Modules DIR
    tmpf    = open('%s/version.txt'\
        %ModuDir,'r')
    __version__=tmpf.readline().strip()                      
else:
    __version__='no version info.'
#

def prepare_Info(version):
    info    = ['#Filename:  xDH',
        '#Authors :  Igor Ying Zhang, Xin Xu',
        '#Version :  %s' % version,
        '#Purpose :  1) Perform XYG3-type doubly hybrid (xDH) calculations using the',
        '               Gaussian package (available for Gaussian 03, 09, and 16).',
        '            2) Share the same input/output format of the Gaussian package',
        '#Refs :     If you publish the xDH results using the xDH script, please make',
        '            sure to cite the references properly',
        '            1) For the original XYG3 method, please cite:',
        '               Zhang, Y.; Xu, X.; Goddard, W. A. Doubly Hybrid Density',
        '               Functional for Accurate Descriptions of Nonbond Interactions,',
        '               Thermochemistry, and Thermochemical Kinetics.',
        '               Proc. Natl. Acad. Sci. USA 2009, 106 (13), 4963-4968.',
        '            2) For the lower-scaling XYGJ-OS method, please cite:',
        '               Zhang, I. Y.; Xu, X.; Jung, Y.; Goddard, W. A. A Fast Doubly',
        '               Hybrid Density Functional Method Close to Chemical Accuracy',
        '               Using a Local Opposite Spin Ansatz.',
        '               Proc. Natl. Acad. Sci. USA 2011, 108 (50), 19896-19900.',
        '            3) For the family of xDH@B3LYP methods, please cite:',
        '               Zhang, I. Y.; Xu, X. Exploring the Limits of the XYG3-Type',
        '               Doubly Hybrid Approximations for the Main-Group Chemistry:',
        '               The XDH@B3LYP Model.',
        '               J. Phys. Chem. Lett. 2021, 12, 2638-2644.'\
        ]
    return info

def parse_input(argv,version):
    global iprint
    global __gaussian__
    tmpargv = [x.lower() for x in argv]
    if '-h' in tmpargv:
        print('Usage   :  run_xDH_using_Gaussian.py [options] file')
        print("Try 'run_xDH_using_Gaussian.py --help' for more information")
        sys.exit(1)
    if '-v' in tmpargv or '--version' in tmpargv:
        tmpS ='run_xDH_using_Gaussian.py (powered by Python 3.x) V%s\n' %(version)
        tmpS = tmpS + 'Copyright (C) 2021 the XDFT@Fudan research group.\n'
        tmpS = tmpS + 'License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.\n'
        tmpS = tmpS + 'This is free software: you are free to change and redistribute it.\n'
        tmpS = tmpS + 'There is NO WARRANTY, to the extent permitted by law.\n'
        tmpS = tmpS + '\n'
        tmpS = tmpS + 'Report bugs to Igor Ying Zhang@Fudan University\n'
        tmpS = tmpS + '(Email: igor_zhangying@fudan.edu.cn)'
        print(tmpS)
        sys.exit(1)
    if '--help' in tmpargv:
        with open('%s/Usage' %ModuDir) as rf:
            rfs = rf.read()
            print(rfs)
        sys.exit(1)
    if '-p' in tmpargv:
        i = tmpargv.index('-p') + 1
        try: 
            iprint = int(tmpargv[i])
        except:
            print('Error in specifying the print option "-p"\n')
            print('please use the option of "--help" for more message')
    if '-g' in tmpargv:
        i = tmpargv.index('-g') + 1
        try: 
            __gaussian__ = int(tmpargv[i])
        except:
            print('Error in specifying the version of the selected Gaussian package "-g"\n')
            print('please use the option of "--help" for more message')
    if '-s' in tmpargv:
        i = tmpargv.index('-f') + 1
        try: 
            syncinterval = int(tmpargv[i])
        except:
            print('Error in specifying the syncronous interval of the Gaussian output to the xDH output file "-f"\n')
            print('please use the option of "--help" for more message')
    for xkey in tmpargv:
        if xkey.find('--gaussian-version=')!=-1:
            try:
                __gaussian__ = int(xkey.strip().split('=')[1])
            except:
                print('Error in specifying the version of the selected Gaussian package "--gaussian-version"\n')
                print('please use the option of "--help" for more message')
        if xkey.find('--print-level=')!=-1:
            try:
                iprint = int(xkey.strip().split('=')[1])
            except:
                print('Error in specifying the print option "--print-level"\n')
                print('please use the option of "--help" for more message')
        if xkey.find('--sync-interval=')!=-1:
            try:
                iprint = int(xkey.strip().split('=')[1])
            except:
                print('Error in specifying the sync interval of the Gaussian output to the xDH output file "--sync-interval"\n')
                print('please use the option of "--help" for more message')
    return

def run_xDH(argv=None):
    from os      import remove
    from os.path import isfile
    import sys
    import gaussian_manage as gaum                               # Import private modules
    import xDH_module as xDH                                     # Import private modules
    from  gaussian_manage  import print_Error
    from  gaussian_manage  import print_List
    from  gaussian_manage  import print_String

    if argv is None: argv=sys.argv
         
    FileName    = argv[-1]
    Path, FileName  = os.path.split(os.path.abspath(FileName))
    Name, Extension = os.path.splitext(FileName)

    __info__ = prepare_Info(__version__)
        
    if isfile('%s.xDH' %Name):                                     # Open the output file
        iout    = open('%s.xDH' % Name,'a')
    else:
        iout    = open('%s.xDH' % Name,'w')
        print_List(iout,__info__,2,'%s' % '-'*76+'==')               # Writing the package info.

    print_String(iout,
			'Start the job of "%s" using the Gaussian %02i package'
			% (FileName,__gaussian__),2)
    MainIO    = gaum.GauIO(iout,'%s%s' %(Name,extension),iprint)

    MainIO.KickOptionList    = ['extraoverlay','oniom','opt']         # Disable options for xDH
    MainIO.get_MachAndOpt()
    MainIO.ctrl_Option()
    MainIO.get_TCSGR()


    OptClass    = gaum.OptHandle(iout,MainIO,iprint)
    R5Class    = xDH.xDH(iout,MainIO,OptClass,iprint,__gaussian__,syncinterval)
    if not MainIO.CartesianFlag:
        MainIO.collect_Geom()
    #DFTDClass    = gaum.DFTD(iout,MainIO,OptClass,iprint)
    #if DFTDClass.PureDisp:
    #    DFTDClass.DispClass.get_EngyReal()                           # Get classical disp. engy.
    #    return
    iout.flush()                                                     # Flush the output
    if R5Class.TurnOn:
        #R5Class.gen_OptionList()
        EngyPos = iout.tell()
        iout.write('%s\n' % (' '*640))
        print_String(iout,
            'The following is the output for preparing KS orbitals and density ::',2)
        iout.flush()                                                 # Flush the output
        R5Class.run_Job(sync=True)
        R5Class.collect_EngyReal(EngyPos)                            # Bring energy print to front
        if not OptClass.Opt:
            if iprint>=1:
                print_String(iout,
                    'Job Type :: Single-Point Calculation',1)
            else:
                os.remove('Job_%s.log' % MainIO.JobName)
            iout.write('='*80+'\n')
            iout.write('**%s**\n' % (' '*76))
            iout.write('** THE JOB OF \"%s\" IS DONE%s**\n'
                % (Name,(' '*(54-len(Name)))))
            iout.write('**%s**\n' % (' '*76))
            iout.write('='*80+'\n')
            #del DFTDClass 
            del R5Class 
            del OptClass 
            del MainIO
            return
        else:
            print_Error(iout,'Geom. optimization is not supported'\
                    +' in this version')
    else:
        print_Error(iout, 
            'Normal Gaussian-Job doesn\'t need the xDH package') 

    return

#//////////////////////////////#
#                              #
#                              #
#NOW START TO THE MAIN PROGRAM #
#                              #
#                              #
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
import sys
import re
import os
import os.path

parse_input(sys.argv,__version__)

FileName    = sys.argv[-1]
WorkDir     = os.getcwd().strip()
path, filename  = os.path.split(os.path.abspath(FileName))
name, extension = os.path.splitext(filename)
os.chdir(path)
if os.path.isfile('%s.xDH' %name):
    os.remove('%s.xDH' %name)
f=open('%s%s' %(name,extension),'r')
if f.read().lower().find('--link1--')!=-1:  #For InputFile containing "--link1--"
    f.seek(0)
    SperateInputList=re.split('--[Ll][Ii][Nn][Kk]1--',f.read())
    tmpI=0
    for input in SperateInputList:
        tmpf=open('Link%d.com' %tmpI,'w')
        tmpf.write(input.strip())
        tmpf.write(' \n')
        tmpf.write(' \n')
        tmpf.write(' \n')
        tmpf.close()
        run_xDH(['xDH.py','Link%d.com' %tmpI])
        os.remove('Link%d.com' %tmpI)
        tmpI+=1
else:    #for InputFile without "--Link1--"
   run_xDH(['xDH.py','%s%s' %(name,extension)]) 
os.chdir(WorkDir)
