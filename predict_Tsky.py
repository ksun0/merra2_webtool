#! /usr/bin/env python

from optparse import OptionParser
import dateutil.parser as dparser
import datetime
import pickle
import os
import cgi, cgitb
import matplotlib
matplotlib.use('Agg')
from pylab import *
import merra2Player as m2p

cgitb.enable()

# TODO: figure out how to run this on holybicep01
# with the pyenv activated from within python
# http://stackoverflow.com/questions/28411960/execute-a-command-on-remote-machine-in-python

def sendData():
    HTML_HEADER = 'Content-type: text/html\n\n'
    print HTML_HEADER
    form = cgi.FieldStorage()
    site = form.getvalue('SiteList')
    if site == "Custom":
        siteKeys = ['siteName','lat','lon','alt']
        siteInfo = [form.getvalue(i) for i in siteKeys]
        site = ','.join(siteInfo)
    tstart = form.getvalue('startdate') ##You're gonna have to change this to customize to different date formats, and potentially the use of time.
    tend = form.getvalue('enddate')
    predictTsky(site = site, tstart = tstart, tend = tend)
    print "<html> \n <body> \n %s %s %s \n </body> </html>"%(site,tstart,tend)
    print '<a href = "merra2_output/tsky-%s-%s-%s.png"><img src = "merra2_output/tsky-%s-%s-%s.png" width = 80%> </a>'%(site,tstart,tend,site,tstart,tend)
    print '</body> \n </html>'

def predictTsky(site = "SouthPole", groundData = "2", tstart = None, tend = None, band = 'tipper850', output = 'tsky'):
    if tstart == None:
        print " ### Error: You must specify a start datetime with the -s option.\n"
        exit()
    if site == None:
        print " ### Error: You must specify a location with the -l option.\n"
        exit()
    if tend == None:
        tend = tstart
        single = True
    else:
        single = False

    os.chdir("/n/home10/anwang16/merra2_work")
    m = m2p.merra2Player()
    m.defineSite(sitestr = site)
    Tsky_merra = m.runMERRA(tstart, tend, gndData = groundData) #t, pwv, BK100, BK150, BK220,BK270 tipper850 all saved as lists


    years = [x.year for x in Tsky_merra['t']]
    months = [x.month for x in Tsky_merra['t']]
    days = [x.day for x in Tsky_merra['t']]
    hours = [x.hour for x in Tsky_merra['t']]

    print "Saving results with pickle or other format"
    print Tsky_merra.keys()

    os.chdir("/cgi-bin")
    if output == "pickle":
        with open("merra2_output/%s_%s_%s_gndData%s.p"%(site,tstart,tend,groundData),"wb") as f:
            pickle.dump(Tsky_merra, f)

    elif output == "tsky":
        data = np.array([years,months,days,hours,Tsky_merra[band]])
        data = np.transpose(data)
        with open("merra2_output/Tsky_%s_%s_%s_gndData%s.csv"%(site,tstart,tend,groundData),"w") as f:
            np.savetxt(f ,data, fmt = "%.4i %.2i %.2i %.2i %3.3f", header = "Year, Month, Day, Hour (UTC), Band-Integrated Brightness Temperature (K_rj)")

    elif output == "pwv":
        data = np.array([years,months,days,hours,Tsky_merra['pwv']])
        data = np.transpose(data)
        with open("merra2_output/PWV_%s_%s_%s_gndData%s.csv"%(site,tstart,tend,groundData),"w") as f:
             np.savetxt(f,data, fmt = "%.4i %.2i %.2i %.2i %3.3f", header = "Year, Month, Day, Hour (UTC), PWV (um)")
    elif output == "spectrum":
        pass

    figure(figsize = (20,15))
    plot(Tsky_merra['t'],Tsky_merra[band],'o-')
    ylabel("Tsky (K_rj)")
    xlabel("Date")
    title("Zenith Brightness Temperature (Tsky): %s %s-%s"%(site, tstart,tend))
    savefig("merra2_output/tsky-%s-%s-%s.png"%(site,tstart,tend))
    return Tsky_merra[band][0]

try:
    sendData()
except:
    cgi.print_exception()
