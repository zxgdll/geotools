#!/usr/bin/env python

# This script reads the .travis.yml file
# and runs the commands for testing the build
# process.

import yaml
import subprocess as sub

with open('.travis.yml', 'r') as f:
	cmd = yaml.load(f)

stages = ['install', 'before_script', 'script', 'after_script']

with open('travis.sh', 'w') as f:
	f.write('#!/bin/sh -e\n')
	for k in stages:
		s = cmd.get(k, False)
		if s:
			for c in s:
				f.write(c + '\n')

sub.Popen('chmod +x travis.sh', shell=True).wait()
sub.Popen('./travis.sh', shell=True).wait()
sub.Popen('rm travis.sh', shell=True).wait()


