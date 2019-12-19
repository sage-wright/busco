#!/usr/bin/env python3
#
# This file fills the config.ini file with info matching your environment
#
# python3 busco_configurator.py config.ini.default yourconfig.ini
import sys
import shutil
paths = {}
try:
    sys.argv[1]
    sys.argv[2]
except IndexError:
    print('\nUsage: python3 busco_configurator.py config.ini.default yourconfig.ini\n')
    exit()
for line in open(sys.argv[1]):
    if line.startswith('['):
        name = line.strip().replace('[','').replace(']','')
        if name not in ['busco_run','sepp']:
            paths.update({name:shutil.which(name)})
        elif name == 'sepp':
            paths.update({name:shutil.which('run_sepp.py')})
outp = open(sys.argv[2],'w')
name = ''
for line in open(sys.argv[1]):
    if line.startswith('['):
        name = line.strip().replace('[','').replace(']','')
    if line.startswith('path ='):
        try:
            outp.write('path = %s/\n' % '/'.join(paths[name].split('/')[0:-1]))
        except AttributeError:
            raise SystemExit('Cannot find the path for the command `%s`, add it in your $PATH and rerun this script' % name)
        continue
    elif line.startswith('command ='):
        try:
            outp.write('command = %s\n' % paths[name].split('/')[-1])
        except AttributeError:
            raise SystemExit('Cannot find the path for the command `%s`, add it in your $PATH and rerun this script' % name)
        continue
    outp.write(line)
