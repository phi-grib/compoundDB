from compoundDB import inputtools as it
import sys

passwd = None   # Input the assowrd here
if not passwd:
    print ('You need to set a password for the postgresDB!')
    sys.exit()

conn = it.openconnection(host='gea', dbname='compounds', \
                            user='postgres', password=passwd)
dbID = it.addSource(conn, 'CosmosDB', version= '1', description= 'Contains repeated \
                    dose toxicity data for cosmetic ingredients, as well as data \
                    from in silico models to predict repeated dose toxicity.', \
                    link= 'http://www.cosmostox.eu/what/COSMOSdb/', \
                    date= datetime.date(2013, 12, 1))
fname = 'cosmos_db_v1_2016_04_02_full.sdf'
it.addSubstanceFile(conn, dbID, fname, 'sdf', molID= 'SYSTEM_ID', synonyms=('PREFERRED_NAME', 'INCI_NAME'))