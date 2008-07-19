#!/usr/bin/env python

# Using CapWords for function, method, and class names
# Using underscored_names for variable names
# Using firstWordLower for module names;
# Using CapWords for package names;
# Using ALL_CAPS_WITH_UNDERSCORES for file handles

#standard imports 
import sys
import os
import re

def PadNone(y):
    for item in y:
        yield item
    while True:
        yield None

def IsOption(s):
    """does string s start with a - ? like -Hap"""
    return s[0] == '-'

def MakeMyPanel():
    mypanel = "mypanel <- function(x, y, ...) {\n"
    mypanel += "x <- x[is.finite(y)]\n"
    mypanel += "y <- y[is.finite(y)]\n"
    mypanel += "panel.xyplot(x, y, ...)\n"
    mypanel += "m <- max(x)\n"
    mypanel += "n <- min(x)\n"
    mypanel += "if (m != n) panel.abline(lm(y~x), ...)\n"
    mypanel += "}\n"
    return mypanel


def MakeRCommands3dPlot(rdata_file, title):
    scan_cmd = 'dqdata <- scan(\"' + rdata_file + '\", list(0, 0, 0))\n'
    lib_cmd  = 'library(lattice)\n'
    x_cmd = 'site1_to_site2_distance <- dqdata[[1]]\n'
    y_cmd = 'site2_to_site3_distance <- dqdata[[2]]\n'
    z_cmd = 'disequilibrium <- dqdata[[3]]\n'
    xyplot1_cmd = 'xyplot( disequilibrium ~ site1_to_site2_distance | site2_to_site3_distance, panel=mypanel, par.strip.text=list(cex=0.7), main=\'' + title + '\')\n'
    xyplot2_cmd = 'xyplot( disequilibrium ~ site2_to_site3_distance | site1_to_site2_distance, panel=mypanel, par.strip.text=list(cex=0.7), main=\'' + title + '\')\n'
    r_cmds = scan_cmd + lib_cmd + x_cmd + y_cmd + z_cmd + xyplot1_cmd + xyplot2_cmd 
    print r_cmds
    return r_cmds

def MakeRCommands2dPlot(rdata_file, xlab, ylab, title):
    scan_cmd = 'dqdata <- scan(\"' + rdata_file + '\", list(0, 0))\n'
    plot_cmd = 'plot(dqdata[[1]], dqdata[[2]], xlog = TRUE, xlab=\'' + xlab + '\', ylab=\'' + ylab + '\', main=\'' + title + '\')\n'
    lm_cmd = 'lmfit <- lm(dqdata[[2]] ~ dqdata[[1]])\n'
    abline_cmd = 'abline(lmfit)\n'
    r_cmds = scan_cmd + plot_cmd + lm_cmd + abline_cmd
    print r_cmds
    return r_cmds



list_of_rdatafiles = sys.argv[1]
r_script = sys.argv[2]

gather_plots = True
pdf_file = 'DQplot.pdf'

if len(sys.argv) > 3:
    assert sys.argv[3] == '-gather'
    if sys.argv[4] == 'no':
        gather_plots = False
    else:
        pdf_file = sys.argv[4]
        

R_SCRIPT = open(r_script, "w")
mypanel = MakeMyPanel()
R_SCRIPT.write(mypanel)
if gather_plots == True:
    R_SCRIPT.write('pdf(\"' + pdf_file + '\")\n')

RDATA_LIST = open(list_of_rdatafiles, "r")

for line in RDATA_LIST.readlines():

    #default values for optional parameters
    level = 2
    xlab = ''
    ylab = ''
    zlab = ''
    title = ''

    #ignore comments. comments work the same as python - indicated by #
    #everything after the first # is ignored.
    real_line = line.split('#', 1)[0]
    if real_line == '':
        continue
    print real_line.split()
    tokens = PadNone(real_line.split())

    # rdata_file will contain the output from DQ in the form that can be
    # 'scan'-ned by R. This option *must* be
    # provided
    rdata_file = tokens.next()
    assert not(rdata_file == None)
    assert not(IsOption(rdata_file))

    if gather_plots == False:
        pdf_file = rdata_file.rstrip('rdata') + 'pdf'

    arg = tokens.next()
    while not(arg == None):
        print arg
        assert IsOption(arg), 'Unknown parameter '+str(arg)

        if (arg == '-level'):
            level = int(tokens.next())
            arg = tokens.next()

        if (arg == '-xlab'):
            arg = tokens.next()
            while not(arg == None or IsOption(arg)):
                xlab += arg+' '
                arg = tokens.next()
            continue

        if (arg == '-ylab'):
            arg = tokens.next()
            while not(arg == None or IsOption(arg)):
                ylab += arg+' '
                arg = tokens.next()
            continue

        if (arg == '-zlab'):
            arg = tokens.next()
            while not(arg == None or IsOption(arg)):
                zlab += arg+' '
                arg = tokens.next()
            continue

        if (arg == '-main'):
            arg = tokens.next()
            while not(arg == None or IsOption(arg)):
                title += arg+' '
                arg = tokens.next()
            continue
    print pdf_file, xlab, ylab, zlab, title    
    #rformat_cmd = 'python formatforR.py ' + rdata_file + ' ' + rdata_file
    #print rformat_cmd
    #os.system(rformat_cmd)
    if (level == 2):
        r_cmds = MakeRCommands2dPlot(rdata_file, xlab, ylab, title)
    else:
        assert level == 3
        r_cmds = MakeRCommands3dPlot(rdata_file, title)

    if gather_plots == False:
        R_SCRIPT.write('pdf(\"' + pdf_file + '\")\n')
    R_SCRIPT.write(r_cmds)
    if gather_plots == False:
        R_SCRIPT.write('dev.off()\n')

if gather_plots == True:
    R_SCRIPT.write('dev.off()\n')
R_SCRIPT.close()    

