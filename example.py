from compoundDB import inputtools as it
import datetime

conn = it.openconnection(host='gea', dbname='compounds', user='postgres', password='DBAdmin')

dbID = it.addSource(conn, 'CosmosDB', version= '1', description= 'Contains repeated \
                 dose toxicity data for cosmetic ingredients, as well as data from in silico \
                 models to predict repeated dose toxicity.', \
                 link= 'http://www.cosmostox.eu/what/COSMOSdb/', \
                 date= datetime.date(2013, 12, 1))

fname = 'cosmos_db_v1_2016_04_02_full.sdf'
it.addSubstanceSDFile(conn, dbID, fname, extIDfield= 'SYSTEM_ID', synonymsFields=('PREFERRED_NAME', 'INCI_NAME'))
