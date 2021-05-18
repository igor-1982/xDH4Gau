#!/usr/bin/env python3
#Usage          :: ChkReplace.py [Input] [OrgChk] [FinalChk]
#Purpose        :: To replace the string [OrgChk] by [FinalChk] in [Input] file
#Authors        :: Igor Ying Zhang, and Xin Xu
#Version        :: 0.2(20100614)
#History        :: 0.1) Basic functional
#                  0.2) Focus replacement on the first line which contains '%chk'

import sys

ff = open(sys.argv[1],'r')
tmpList     = ff.readlines()
ff.close()
ff = open(sys.argv[1],'w')
for tmpString in tmpList:
    if tmpString.lower().find('%chk')!=-1:
        tmpString = tmpString.replace(sys.argv[2],sys.argv[3])
    ff.write(tmpString)
ff.close()

