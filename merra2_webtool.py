from flask import Flask, render_template, flash, request
from wtforms import Form, TextField, TextAreaField, validators, StringField, SubmitField
import datetime
import pickle
import os
from pylab import *
import merra2Player as m2p
app = Flask(__name__,static_folder='merra2_output')
app.config.from_object(__name__)
app.config['SECRET_KEY'] = 'scikit4life'

class ReusableForm(Form):
    site = TextField('site:', validators=[validators.required()])
    tstart = TextField('tstart:', validators=[validators.required()])
    tend = TextField('tend:', validators=[validators.required()])

@app.route("/", methods=['GET', 'POST'])
def hello():
    form = ReusableForm(request.form)
    tskyFileName = ""
    imageFile = "SP_timeseries_averaged_2017_groundData2_vapor.png"
    print(form.errors)
    if request.method == 'POST':
        site=request.form['site']
        groundData=2
        tstart=request.form['tstart']
        tend=request.form['tend']
        band="tipper850"
        output="tsky"
        tskyFileName = "Tsky_%s_%s_%s_gndData%s.csv"%(site,tstart,tend,groundData)
        imageFile = "tsky-%s-%s-%s.png"%(site,tstart,tend)
        print tstart,tend,site,groundData,band,output,tskyFileName,imageFile
        

        if form.validate():
            predictTsky(site, groundData, tstart, tend, band, output)
            flash('Successful')
        else:
            flash('Error categorizing article.')
    return render_template('index.html', form=form, tskyFileName = tskyFileName, imageFile = imageFile)

def predictTsky(site = "SouthPole", groundData = 2, tstart = None, tend = None, band = 'tipper850', output = 'tsky'):
    if tstart == None:
        print(" ### Error: You must specify a start datetime with the -s option.\n")
        exit()
    if site == None:
        print(" ### Error: You must specify a location with the -l option.\n")
        exit()
    if tend == None:
        tend = tstart
        single = True
    else:
        single = False

    m = m2p.merra2Player()
    m.defineSite(sitestr = site)
    Tsky_merra = m.runMERRA(tstart, tend, gndData = groundData) #t, pwv, BK100, BK150, BK220,BK270 tipper850 all saved as lists


    years = [x.year for x in Tsky_merra['t']]
    months = [x.month for x in Tsky_merra['t']]
    days = [x.day for x in Tsky_merra['t']]
    hours = [x.hour for x in Tsky_merra['t']]

    print("Saving results with pickle or other format")
    print(Tsky_merra.keys())

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

if __name__ == "__main__":
    app.run(host='0.0.0.0')
