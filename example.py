from inputtools import *

conn = openconnection(host='gea', dbname='compounds', user='postgres', password='DBAdmin')

dbID = addSource(conn, 'CosmosDB', version= '1', description= 'Contains repeated \
                 dose toxicity data for cosmetic ingredients, as well as data from in silico \
                 models to predict repeated dose toxicity.', \
                 link= 'http://www.cosmostox.eu/what/COSMOSdb/', \
                 date= datetime.date(2013, 12, 1))

fname = 'C:\\Users\\bet\\Google Drive\\Inditex\\Compound inventory\\Annotated compounds\\\
CosmosDB\\v1\\db\\cosmos_db_v1_2016_04_02_full.sdf'
addSubstanceFile(conn, dbID, fname, 'sdf', molID= 'SYSTEM_ID', synonyms=('PREFERRED_NAME', 'INCI_NAME'))