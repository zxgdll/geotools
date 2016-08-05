#!/usr/bin/env python

import sys
import os
import subprocess

sys.path.append(os.path.realpath('../scripts'))
os.environ['PATH'] += os.pathsep + os.path.realpath('../makefiles')
os.environ['LAS2CSRS_DATA'] = os.pathsep + os.path.realpath('../share')
"""
Usage: las2csrs <options> <las files>

This program converts coordinates from LAS files from any reference frame to NAD83(CSRS),
between any two epochs.

If orthometric heights are used, be sure to provide SRSes with geoid parameters, or that the source
files contain such information. SRSes are entered in the form, 'epsg:<horizontal code>+<vertical code>'.
The + is only required if orthometric heights are desired.

 -o     Overwrite existing files. Defaults to false.
 -v     Verbose output.
 -d 	Destination folder. Required.
 -fs    The source SRS. If left out, will use the first las file's SRS. If the file has no SRS, will fail.
 -ts    The destination SRS. Required.
 -fe    The source epoch. Required.
 -te    The destination epoch. Required.
 -f     The source reference frame. Required.
"""
def test():
    
    tests = [
            ('itrf88', 1986., 2011., 26912, 2956, 470000.000, 6520000.000, 200.000, 470000.800, 6519999.074, 200.431), 
            ('itrf89', 1992., 2002., 26910, 3157, 211704.236, 5617008.921,  53.000, 211705.329, 5617008.261,  53.248),
            ('itrf90', 1990., 2001., 32611, 2955, 467473.356, 6430442.112, 981.230, 467474.288, 6430441.270, 981.517),

            ('itrf91', 1986., 2011., 26912, 2956, 470000.000, 6520000.000, 200.000, 470000.805, 6519999.017, 200.416), 
            ('itrf92', 1992., 2002., 26910, 3157, 211704.236, 5617008.921,  53.000, 211705.320, 5617008.252,  53.211),
            ('itrf93', 1990., 2001., 32611, 2955, 467473.356, 6430442.112, 981.230, 467474.308, 6430441.257, 981.494),

            ('itrf94', 1986., 2011., 26912, 2956, 470000.000, 6520000.000, 200.000, 470000.817, 6519999.028, 200.396), 
            ('itrf96', 1992., 2002., 26910, 3157, 211705.325, 5617008.251,  53.000, 211706.414, 5617007.581,  53.196),
            ('itrf97', 1990., 2001., 32611, 2955, 467473.356, 6430442.112, 981.230, 467474.295, 6430441.268, 981.488),

            ('itrf2000',  1986., 2011., 26912, 2956, 470000.000, 6520000.000, 200.000, 470000.809, 6519999.031, 200.395), 
            ('nad83csrs', 1992., 2002., 26910, 3157, 211704.236, 5617008.921,  53.000, 211704.274, 5617008.971,  53.005)
    ]
    
    def rnd(n):
        return int(round(n * 1000.)) / 1000.

    for test in tests:
        
        ffrom, efrom, eto, fsrid, tsrid, x0, y0, z0, x1, y1, z1 = test
        params = map(str, test[:8])
        p = sub.Popen(['las2csrs', '-o', '-v', '-d', 'output', '-fs', test[3], '-ts', test[4], '-fe', test[1], '-te', test[2], '-f', test[0], 'data/las2csrs.las'], stdout=sub.PIPE, stderr=sub.PIPE)
        stdout, stderr = p.communicate()
        x2, y2, z2 = map(float, stdout.split())
        print 'input: %f %f %f; output: %f %f %f; expected: %f %f %f' % (x0, y0, z0, x2, y2, z2, x1, y1, z1)
        if rnd(x1) != rnd(x2) or rnd(y1) != rnd(y2) or rnd(z1) != rnd(z2):
            print 'Failed', test
        else:
            print 'Passed'
 
if __name__ == '__main__':
    
    test()

    
