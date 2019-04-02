#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 14:01:26 2018
Script to clean all the temperaturay files
@author: walkandthinker
"""
import os
import shutil

currentdir=os.getcwd()

print('*** Start to clean the workspace ...')
print('*** We are in folder:%s\n'%(currentdir))

exe=0;o=0;vtu=0;jit=0;tmp=0;cmake=0;peacock=0;check=0
for subdir,dirs,files in os.walk(currentdir):
    #>>> clean files
    for file in files:
        if 'asfem' in file or 'ASFEM' in file:
            try:
                exe+=1
                removepath=subdir+'/'+file
                os.remove(removepath)
            except:
                print('%s is not here'%(file))
        elif '.vtu' in file:
            try:
                vtu+=1
                removepath=subdir+'/'+file
                os.remove(removepath)
            except:
                print('%s is not here'%(file))
        elif '.o' in file or '.pos' in file:
            try:
                o+=1
                removepath=subdir+'/'+file
                os.remove(removepath)
            except:
                print('%s is not here'%(file))
        elif 'tmp' in file:
            try:
                tmp+=1
                removepath=subdir+'/'+file
                os.remove(removepath)
            except:
                print('%s is not here'%(file))
        elif 'cmake' in file or 'CMakeCache' in file:
            try:
                cmake+=1
                removepath=subdir+'/'+file
                os.remove(removepath)
            except:
                print('%s is not here'%(file))
        elif 'Makefile' in file:
            try:
                cmake+=1
                removepath=subdir+'/'+file
                os.remove(removepath)
            except:
                print('%s is not here'%(file))
    #>> clean folder
    for dir in dirs:
        if '.idea' in dir:
            try:
                removepath=subdir+'/'+dir
                print('remove folder: ',dir)
                shutil.rmtree(removepath)
                IdeaRemove=True
            except:
                if(not IdeaRemove):
                    print('%s is not here'%(dir))
        elif '.jit' in dir:
            try:
                jit+=1
                removepath=subdir+'/'+dir
                print('remove folder: ',dir)
                shutil.rmtree(removepath)
                IdeaRemove=True
            except:
                if(not IdeaRemove):
                    print('%s is not here'%(dir))
        elif ('my_checkpoint' in dir) or ('checkpoint' in dir):
            try:
                check+=1
                removepath=subdir+'/'+dir
                print('remove folder: ',dir)
                shutil.rmtree(removepath)
                IdeaRemove=True
            except:
                if(not IdeaRemove):
                    print('%s is not here'%(dir))
        elif ('cmake-build-debug' in dir) or ('CMakeFiles' in dir):
            try:
                cmake+=1
                removepath=subdir+'/'+dir
                print('remove folder: ',dir)
                shutil.rmtree(removepath)
                IdeaRemove=True
            except:
                if(not IdeaRemove):
                    print('%s is not here'%(dir))

print('Remove %4d exe files!'%(exe))
print('Remove %4d tmp files!'%(tmp))
print('Remove %4d .vtu files!'%(vtu))
print('Remove %4d .o files!'%(o))
print('Remove %4d .jit folder!'%(jit))
print('Remove %4d cmake folder!'%(cmake))
print('Remove %4d check folder!'%(check))
print('***** Clean task finished !!!   *****')
