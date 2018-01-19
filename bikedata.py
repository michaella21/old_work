import numpy as np
import pandas as pd
import os
import datetime

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

bikes=pd.read_csv(r'big_data.csv')
bikes.columns=['bikeid', 'end_station_id','start_station_id', 'starttime','stoptime']
#change the data format from strings to datetime
bikes['starttime']=pd.to_datetime(bikes['starttime'])
bikes['stoptime']=pd.to_datetime(bikes['stoptime'])

#compute the travel time of each case, save in travelT column
bikes['travelT']=(bikes['stoptime']-bikes['starttime'])

print "A1:Average travel time of each bike ride is %r and the number of total travel is %i " \
    %(bikes['travelT'].mean(),bikes['travelT'].count())


#Find the number of missing data; need to count when the last end station id is not matching with the start station id

b_sorted=bikes.sort(['bikeid','starttime'])
b_sorted['missing']=(b_sorted.bikeid ==b_sorted.bikeid.shift(1)) &(b_sorted.start_station_id != b_sorted.end_station_id.shift(1))
missing=np.sum(b_sorted['missing'])
print "There are %d missing cases and it is %f percentage of the whole data" % (missing,100*float(missing)/len(bikes))
b_sorted['incorrect']=(b_sorted.bikeid ==b_sorted.bikeid.shift(1)) &(b_sorted.starttime < b_sorted.stoptime.shift(1))
incorrect=np.sum(b_sorted['incorrect'])
print "There are %d incorrect cases and it is %f percentage of the whole data" % (incorrect,100*float(incorrect)/len(bikes))

b_sorted['staying_at_s']=(b_sorted.starttime-b_sorted.stoptime.shift(1))*(1-b_sorted['missing'])*(1-b_sorted['incorrect'])
b_sorted['staying_at_e']=b_sorted['staying_at_s'].shift(-1)
correct=int(len(bikes)-missing-incorrect)
total_stay=np.sum(b_sorted['staying_at_s'])
print "Excluding the missing or/and the incorrect data, a bike satys at each station on average %r" %(total_stay/(correct))

#Find the information about the station '8f0f64' and '4a4b61': the number of bikes and average time stayed at each station

b_4a=b_sorted[((b_sorted['start_station_id']=='4a4b61')|(b_sorted['end_station_id']=='4a4b61'))].sort('bikeid')
b_8f=b_sorted[((b_sorted['start_station_id']=='8f0f64')|(b_sorted['end_station_id']=='8f0f64'))].sort('bikeid')

s=datetime.datetime(2013,10,30,0,0,0)
e=datetime.datetime(2013,10,31,0,0,0)

b_4a.drop(b_4a.columns[[5,6]],axis=1, inplace=True)
b_4a=b_4a.reset_index(drop=True)

b_8f.drop(b_8f.columns[[5,6]],axis=1, inplace=True)
b_8f=b_8f.reset_index(drop=True)

b_4a_s=b_4a[(b_4a['start_station_id']=='4a4b61') & (b_4a['starttime']-e<= b_4a['staying_at_s']) & (b_4a['starttime']>= s)]
b_4a_e=b_4a[(b_4a['end_station_id']=='4a4b61') & (b_4a['stoptime']+b_4a['staying_at_e']> s) & (b_4a['stoptime']<e)]

b_8f_s=b_8f[(b_8f['start_station_id']=='8f0f64') & (b_8f['starttime']-e<= b_8f['staying_at_s']) & (b_8f['starttime']>= s)]
b_8f_e=b_8f[(b_8f['end_station_id']=='8f0f64') & (b_8f['stoptime']+b_8f['staying_at_e']> s) & (b_8f['stoptime']<e)]

b_4a_c=pd.DataFrame(columns=['bikeid','stay_from','stay_till','how_long'])
b_8f_c=pd.DataFrame(columns=['bikeid','stay_from','stay_till','how_long'])

for i, row in b_4a_s.iterrows():
    if (b_4a_s.loc[i,'stoptime']-b_4a_s.loc[i,'staying_at_s'])<= s:
        b_4a_c.loc[i,'bikeid']=b_4a_s.loc[i,'bikeid']
        b_4a_c.loc[i,'stay_from']=0
    else: 
        b_4a_c.loc[i,'stay_from']=(b_4a_s.loc[i,'starttime']-b_4a_s.loc[i,'staying_at_s']).to_datetime().hour
        b_4a_c.loc[i,'bikeid']=b_4a_s.loc[i,'bikeid']
    if b_4a_s.loc[i,'starttime']>=e:
        b_4a_c.loc[i,'stay_till']=24
    else:
        b_4a_c.loc[i,'stay_till']=(b_4a_s.loc[i,'starttime'].to_datetime().hour)+1
    b_4a_c.loc[i,'how_long']=b_4a_s.loc[i,'staying_at_s']

for i, row in b_4a_e.iterrows():
    if (b_4a_e.loc[i,'stoptime']+b_4a_e.loc[i,'staying_at_e']) >= e:
        b_4a_c.loc[i+len(b_4a_s),'bikeid']=b_4a_e.loc[i,'bikeid']
        b_4a_c.loc[i+len(b_4a_s),'stay_till']=24
    else: 
        b_4a_c.loc[i+len(b_4a_s),'stay_till']=((b_4a_e.loc[i,'stoptime']+b_4a_e.loc[i,'staying_at_e']).to_datetime().hour)+1
        b_4a_c.loc[i+len(b_4a_s),'bikeid']=b_4a_e.loc[i,'bikeid']
    if b_4a_e.loc[i,'stoptime']<=s:
        b_4a_c.loc[i+len(b_4a_s),'stay_from']=0
    else:
        b_4a_c.loc[i+len(b_4a_s),'stay_from']=(b_4a_e.loc[i,'stoptime'].to_datetime().hour)
    b_4a_c.loc[i+len(b_4a_s),'how_long']=b_4a_e.loc[i,'staying_at_e']
    
b_4a_c=b_4a_c.drop_duplicates()
b_4a_c=b_4a_c.sort(['bikeid']).reset_index(drop=True)

for i, row in b_8f_s.iterrows():
    if (b_8f_s.loc[i,'stoptime']-b_8f_s.loc[i,'staying_at_s'])<= s:
        b_8f_c.loc[i,'bikeid']=b_8f_s.loc[i,'bikeid']
        b_8f_c.loc[i,'stay_from']=0
    else: 
        b_8f_c.loc[i,'stay_from']=(b_8f_s.loc[i,'starttime']-b_8f_s.loc[i,'staying_at_s']).to_datetime().hour
        b_8f_c.loc[i,'bikeid']=b_8f_s.loc[i,'bikeid']
    if b_8f_s.loc[i,'starttime']>=e:
        b_8f_c.loc[i,'stay_till']=24
    else:
        b_8f_c.loc[i,'stay_till']=(b_8f_s.loc[i,'starttime'].to_datetime().hour)+1
    b_8f_c.loc[i,'how_long']=b_8f_s.loc[i,'staying_at_s']

for i, row in b_8f_e.iterrows():
    if (b_8f_e.loc[i,'stoptime']+b_8f_e.loc[i,'staying_at_e']) >= e:
        b_8f_c.loc[i+len(b_8f_s),'bikeid']=b_8f_e.loc[i,'bikeid']
        b_8f_c.loc[i+len(b_8f_s),'stay_till']=24
    else: 
        b_8f_c.loc[i+len(b_8f_s),'stay_till']=((b_8f_e.loc[i,'stoptime']+b_8f_e.loc[i,'staying_at_e']).to_datetime().hour)+1
        b_8f_c.loc[i+len(b_8f_s),'bikeid']=b_8f_e.loc[i,'bikeid']
    if b_8f_e.loc[i,'stoptime']<=s:
        b_8f_c.loc[i+len(b_8f_s),'stay_from']=0
    else:
        b_8f_c.loc[i+len(b_8f_s),'stay_from']=(b_8f_e.loc[i,'stoptime'].to_datetime().hour)
    b_8f_c.loc[i+len(b_8f_s),'how_long']=b_8f_e.loc[i,'staying_at_e']
        
b_8f_c=b_8f_c.drop_duplicates()
b_8f_c=b_8f_c.sort(['bikeid']).reset_index(drop=True)

b_4a_c['how_long']=pd.to_datetime(b_4a_c['how_long'])
b_8f_c['how_long']=pd.to_datetime(b_8f_c['how_long'])
oct30=pd.DataFrame(0,index=range(24), columns=['4a4b61','8f0f64'])

for i, row in b_4a_c.iterrows():
    for j,row in oct30.iterrows():
        if b_4a_c.loc[i,'stay_from']<= j and j< b_4a_c.loc[i,'stay_till']:
            oct30.loc[j,'4a4b61'] +=1        

for i, row in b_8f_c.iterrows():
    for j,row in oct30.iterrows():
        if b_8f_c.loc[i,'stay_from']<= j and j< b_8f_c.loc[i,'stay_till']:
            oct30.loc[j,'8f0f64'] +=1

print "There are %d bikes visited the station 4a4b61 on Oct/30 and each of them stays at the station %r on average" %(len(b_4a_c)+1, b_4a_c['how_long'].mean())
print "There are %d bikes visited the station 8f0f4 on Oct/30 and each of them stays at the station %r on average" %(len(b_8f_c)+1, b_8f_c['how_long'].mean())


