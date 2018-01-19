
#!/Users/michelle/Library/Enthought/Canopy_64bit/User/bin/python
# -*- coding: utf-8 -*-

import MySQLdb
import numpy as np
import pandas as pd

db = MySQLdb.connect(host="dataincubator.celondxt5ghj.us-east-1.rds.amazonaws.com", # your host, usually localhost
user="michaella21", # your username
passwd="dataincubator", # your password
db="mydb") # name of the data base

# you must create a Cursor object. It will let
# you execute all the queries you need
cur = db.cursor() 
cur.execute(#"ALTER TABLE `call_center` ADD INDEX `agency` (`agency`) \
"SELECT count(`agency`)\
FROM call_center")

total_cases=cur.fetchone()
total_cases=float(total_cases[0])

cur.execute("SELECT count(agency) AS cnt, agency\
FROM call_center\
GROUP BY agency\
ORDER BY cnt DESC")
comp_by_agency=cur.fetchall()
list_complaint_agency=[]

for i in range(len(comp_by_agency)):
list_complaint_agency.append(list(comp_by_agency[i]))
list_complaint_agency[i].insert(3,round(list_complaint_agency[i][0]/total_cases,10))

"""
Q1. What fraction of complaints are associated with the 2nd most popular agency? 0.171329548 
"""
print list_complaint_agency[1]


"""
Q2. What is the most 'surprising' complaint type when conditioned on a borough? That is, what is the largest ratio of the conditional probability of a complaint type given a specified borough divided by the unconditioned probability of that complaint type?
18.2795184415
"""

#cur.execute("ALTER TABLE `call_center` ADD INDEX `complaint_type` (`complaint_type`)")
#cur.execute("ALTER TABLE `call_center` ADD INDEX `complaint_borough` (`complaint_type`,`borough`)")

cur.execute("CREATE TABLE `borough_count` (`borough` varchar(255) DEFAULT NULL,\
`borough_count` int(11) DEFAULT NULL")

cur.execute("CREATE TABLE `complaint` (`complaint_type` varchar(255) DEFAULT NULL,\
`count_complaint_type` int(11) DEFAULT NULL") 
cur.execute("CREATE TABLE `complaint_borough` (`complaint_type` varchar(255) DEFAULT NULL,\
`borough` varchar(255) DEFAULT NULL,`count_complaint_type_borough` int(11) DEFAULT NULL")

cur.execute("INSERT INTO complaint SELECT `complaint_type`, count(`complaint_type`) \
FROM call_center GROUP BY `complaint_type`")

cur.execute("INSERT INTO complaint_borough SELECT `complaint_type`, `borough`, count(`id`)\
FROM call_center GROUP BY `complaint_type`, `borough`")
cur.execute("INSERT INTO borough_count SELECT `borough`, count(`id`) FROM call_center GROUP BY `borough`")

cur.execute("SELECT cb.*, bc.*, c.*, (cb.count_complaint_type_borough/bc.borough_count)/(c.count_complaint_type/10138724) as ratio\
FROM complaint_borough AS cb\
INNER JOIN borough_count as bc ON bc.borough=cb.borough\
INNER JOIN complaint as c ON c.complaint_type=cb.complaint_type\
ORDER BY ratio DESC") #10138724 =total_cases

ans=cur.fetchone() 
cond=float(ans[2])/float(ans[4])
uncond=float(ans[6])/total_cases
ratio=round(cond/uncond,10)


"""
Q3. What is the distance (in degrees) between the 90% and 10% percentiles of degrees latitude? 40.85148752696949-40.62212465681701=0.2293628702
"""

cur.execute("SELECT latitude FROM call_center WHERE latitude <> "" ORDER BY latitude LIMIT 1013871,2 ")
per10=cur.fetchall()#both numbers are same as 40.62212465681701

#cur.execute("SELECT `latitude` FROM call_center WHERE `latitude`<>"" ORDER BY latitude DESC LIMIT 1013871, 2")

per90=cur.fetchall()#both numbers are same as 40.85148752696949

ans=round(per90[0]-per10[0])

"""
Q4. Let's estimate the area that 311 supports. Suppose calls are 2D normally distributed on the surface of the earth with mean and standard deviation given by those of the latitude and longitude. How many square kilometers is the single-standard-deviation ellipse?

As noted in the question, suppose that the (population) mean and (population) standard deviation of 2D normal are mean/stdv of given data
"""

cur.execute("SELECT STDDEV_POP(`latitude`),STDDEV_POP(`longitude`) FROM call_center \
WHERE `latitude` IS NOT NULL and `longitude` IS NOT NULL")

stdev=cur.fetchall()

cur.execute("SELECT ( SUM( `latitude` * `longitude` ) - SUM( `latitude` ) * SUM( `longitude` ) / COUNT( `latitude` ) ) / COUNT( `latitude` )\
WHERE `latitude` IS NOT NULL and `longitude` IS NOT NULL")

cov=cur.fetchall()

