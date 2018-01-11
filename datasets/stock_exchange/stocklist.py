#coding: utf8

import pandas as pd
import numpy as np
from scipy.stats import levy_stable, norm, levy
import datetime
import time

def calculate():
    
    # load market data
    filename = 'sp500-list'
    io = pd.read_csv(filename+'.csv', sep="\t", header=None, skiprows=1, engine='python')

    io = io[io[6]<='1997-01-01']
    io = io.reset_index()
    print io
    print 'Stocks:', len(io)
    
    # write stocklist
    '''
    io[[0]].to_csv('stocklist.csv')
    '''

    # download data from stooq
    # by visiting urls of type 'https://stooq.pl/q/d/l/?s=abt.us&i=d'
    '''
    import webbrowser
    
    for stockcode in io[0]:
        print stockcode
        url = 'https://stooq.pl/q/d/l/?s='+stockcode.lower()+'.us&i=d'
        webbrowser.open(url, new=0, autoraise=True)
        time.sleep(5)
    '''

    # compile data to one dataset
    
    ALLDATA = pd.DataFrame()
    for stockcode in io[0]:
        sdata = pd.read_csv('raw_data/' +stockcode.lower()+ '_us_d'+'.csv', sep=",", names=['date', 'o', 'h', 'l', stockcode, 'v'], parse_dates=True, index_col = 'date', header=None, skiprows=1, engine='python')

        # drop if earliest date is too late
        if sdata.index[0] > datetime.datetime.strptime('1997-01-01','%Y-%m-%d'):
            continue
        ALLDATA = ALLDATA.join(sdata[stockcode], how='outer')

    ALLDATA = ALLDATA.dropna()
    #print ALLDATA
    
    #ALLDATA.to_csv('alldata.csv')


    # detrend data and write to output file

    ALLDETREND = pd.DataFrame(index=ALLDATA.index)
    for column in ALLDATA:
        print column
        ALLDATA[column+'shift'] = ALLDATA[column].shift(1)
        ALLDATA.dropna()
        
        # detrending as log returns
        ALLDETREND[column] = np.log(ALLDATA[column]/ALLDATA[column+'shift'])
        
        # alternative: detrending as pct changes
        #ALLDETREND[column] = ALLDATA[column].pct_change(1)

    ALLDETREND = ALLDETREND.dropna()
    #print ALLDETREND

    #ALLDETREND.to_csv('detrended-log.csv')


calculate()