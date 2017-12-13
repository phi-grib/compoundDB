import sys, datetime

import psycopg2
from psycopg2 import extras
from rdkit import Chem
import pandas as pd

def getSubsID(conn, sourceID, extID):
    """
    Return the ID for a substance given the source and external ID.
    Arguments:
          - conn: psycopg2 connection to the database.
          - sourceID: id for the source of origin from the 'source' table.
          - extID: id for the subsance in the source of origin.
    """
    curs = conn.cursor()
    cmd = "SELECT id FROM substance WHERE sourceid = %s AND externalid = %s;"
    curs.execute(cmd, (sourceID, extID))
    subsID = curs.fetchone()
    conn.commit()

    if subsID is not None:
        subsID = subsID[0]

    return subsID

def getSourceID(conn, sourceName, version= None):
    """
    Return the ID for a source given the source name and optionally the version. 
    If no version is provided, the last version will be used.
    Arguments:
          - conn: psycopg2 connection to the database.
          - sourceName: Source name.
          - version: Optional. Source's version (default: None). If None, 
          the last version will be used.
    """
    curs = conn.cursor()
    if version is None:
        # Get the last version        
        cmd = "SELECT id FROM source WHERE name = %s ORDER BY version DESC LIMIT 1;"
        curs.execute(cmd, (sourceName,))
    else:
        # Get a specific version
        cmd = "SELECT id FROM source WHERE name = %s AND version = %s;"
        curs.execute(cmd, (sourceName, version))

    sourceID = curs.fetchone()[0]
    conn.commit()

    return sourceID

def getStructureFromSyn(conn, syn):
    """
    Return the smiles for the given synonym.
    Arguments:
          - conn: psycopg2 connection to the database.
          - syn: Synonym.
    """
    curs = conn.cursor()
    cmd = "SELECT smiles \
            FROM substance AS subs \
            INNER JOIN synonym AS syn ON subs.id = syn.subsid\
            WHERE name = %s;"
    curs.execute(cmd, (syn,))
    smiles = curs.fetchone()
    conn.commit()

    if smiles is not None:
        smiles = smiles[0]

    return smiles

def getAllSourcesWithInfo(conn):
    """
    Return a pandas dataframe with all sources and the number of substances and compounds
    by version.
    """
    cmd = "SELECT name, version, \
                  COUNT(DISTINCT subs.externalid) AS entries, \
                  COUNT(DISTINCT subs.externalid) FILTER (WHERE subs.smiles IS NOT NULL) AS substances, \
                  COUNT(DISTINCT cmpd.inchikey) AS compounds \
	       FROM public.source as s \
           INNER JOIN public.substance AS subs ON s.id = subs.sourceid\
	       LEFT OUTER JOIN public.subs_cmpd AS sc ON sc.subsid = subs.id \
	       LEFT OUTER JOIN public.compound AS cmpd ON sc.cmpdid = cmpd.id \
	       GROUP BY name, version"
    df = pd.read_sql(cmd, con=conn)

    return df

def getAllSourcesWithInfoLimitSubstances(conn, ids):
    """
    Return a pandas dataframe with all sources and the number of substances and compounds
    by version.
    Arguments:
          - conn: psycopg2 connection to the database.
          - ids: Tuple with external IDs of compounds to consider.
    """
    cmd = "SELECT name, version, \
                  COUNT(DISTINCT subs.externalid) AS entries, \
                  COUNT(DISTINCT subs.externalid) FILTER (WHERE subs.smiles IS NOT NULL) AS substances, \
                  COUNT(DISTINCT cmpd.inchikey) AS compounds \
	       FROM public.source as s \
           INNER JOIN public.substance AS subs ON s.id = subs.sourceid\
	       LEFT OUTER JOIN public.subs_cmpd AS sc ON sc.subsid = subs.id \
	       LEFT OUTER JOIN public.compound AS cmpd ON sc.cmpdid = cmpd.id \
           WHERE subs.externalid in {} \
	       GROUP BY name, version".format(str(ids))
    df = pd.read_sql(cmd, con=conn)

    return df

def getOneSourceWithInfo(conn, sourceID= None, sourceName= None, version= None):
    """
    Return a pandas dataframe with all sources and the number of substances and compounds
    by version.
    Arguments:
        - conn: psycopg2 connection to the database.
        - sourceID: Optional. Internal source ID if avaialble. Otherwise, it will be retrieved    from the source name and version.
        - sourceName: Optional. Source name.
        - version: Optional. Source's version (default: None). If None, 
          the last version will be used.
    """
    if sourceID is None:
        sourceID = getSourceID(conn, sourceName, version)

    cmd = "SELECT name, version, \
                  COUNT(DISTINCT subs.externalid) AS entries, \
                  COUNT(DISTINCT subs.externalid) FILTER (WHERE subs.smiles IS NOT NULL) AS substances, \
                  COUNT(DISTINCT cmpd.inchikey) AS compounds \
	       FROM public.source as s \
           INNER JOIN public.substance AS subs ON s.id = subs.sourceid\
	       LEFT OUTER JOIN public.subs_cmpd AS sc ON sc.subsid = subs.id \
	       LEFT OUTER JOIN public.compound AS cmpd ON sc.cmpdid = cmpd.id \
           WHERE s.id = %s \
	       GROUP BY name, version"
    df = pd.read_sql(cmd, con=conn, params=(sourceID,))

    return df

# def getMulipleSourcesWithInfo(conn, sourceID= None, sourceName= None, version= None):
#     """
#     Return a pandas dataframe with all sources and the number of substances and compounds
#     by version.
#     Arguments:
#         - conn: psycopg2 connection to the database.
#         - sourceID: Optional. List of internal source IDs if avaialble. Otherwise, they will be retrieved from the source name and version.
#         - sourceName: Optional. List of source names.
#         - version: Optional. List of source versions (default: None). If None, 
#           the last version will be used.
#     """
#     if sourceID is None:
#         sourceID = getSourceID(conn, sourceName, version)

#     cmd = "SELECT name, version, \
#                   COUNT(DISTINCT subs.externalid) AS entries, \
#                   COUNT(DISTINCT subs.externalid) FILTER (WHERE subs.smiles IS NOT NULL) AS substances, \
#                   COUNT(DISTINCT cmpd.inchikey) AS compounds \
# 	       FROM public.source as s \
#            INNER JOIN public.substance AS subs ON s.id = subs.sourceid\
# 	       LEFT OUTER JOIN public.subs_cmpd AS sc ON sc.subsid = subs.id \
# 	       LEFT OUTER JOIN public.compound AS cmpd ON sc.cmpdid = cmpd.id \
# 	       GROUP BY name, version"
#     df = pd.read_sql(cmd, con=conn)

#     return df

def getSubstancesFromSource(conn, sourceID= None, sourceName= None, version= None, includeempty= False):
    """
    Retruns a pandas dataframe with all substances for the given source.
    The dataframe contains:
        - exernalid
        - smiles
        - inchi
        - inchikey
        - link
    Arguments:
        - conn: psycopg2 connection to the database.
        - sourceID: Optional. Internal source ID if avaialble. Otherwise, it will be retrieved    from the source name and version.
        - sourceName: Optional. Source name.
        - version: Optional. Source's version (default: None). If None, 
          the last version will be used.
        - includeempty: Optional. If True, returns all substances, including those
        without a chemical structure.
    """
    if sourceID is None:
        sourceID = getSourceID(conn, sourceName, version)

    if includeempty:
        cmd = "SELECT externalid, smiles, inchi, inchikey, subs.link \
                FROM public.substance AS subs \
                INNER JOIN public.source AS s ON s.id = subs.sourceid \
                WHERE s.id = %s"
    else:
        cmd = "SELECT externalid, smiles, inchi, inchikey, subs.link \
                FROM public.substance AS subs \
                INNER JOIN public.source AS s ON s.id = subs.sourceid \
                WHERE s.id = %s AND smiles IS NOT NULL"

    df = pd.read_sql(cmd, con=conn, params=(sourceID,))

    return df
    
def getSubstancesFromSynonyms(conn, synList):
    """
    Retruns a pandas dataframe with all substances for the given source.
    The dataframe contains:
        - synonym
        - exernalid
        - smiles
        - inchi
        - inchikey
        - link 
    Arguments:
        - conn: psycopg2 connection to the database.
        - synList: List or tuple with synonyms.
    """
    cmd = "SELECT DISTINCT ON (synonym, externalid) \
            name AS synonym, externalid, smiles, inchi, inchikey, link \
            FROM public.synonym AS syn \
            LEFT JOIN public.substance AS subs ON syn.subsid = subs.id \
            WHERE name in {}".format(tuple(synList))

    df = pd.read_sql(cmd, con=conn)

    return df


def getCompoundsFromSource(conn, sourceID= None, sourceName= None, version= None, includeempty= False):
    """
    Retruns a pandas dataframe with all compounds for the given source.
    Arguments:
        - conn: psycopg2 connection to the database.
        - sourceID: Optional. Internal source ID if avaialble. Otherwise, it will be retrieved    from he source name and version.
        - sourceName: Source name.
        - version: Optional. Source's version (default: None). If None, 
          the last version will be used.
        - includeempty: Optional. If True, returns all substances, including those
        without a chemical structure.
    """
    if sourceID is None:
        sourceID = getSourceID(conn, sourceName, version)

    if includeempty:
        cmd = "SELECT externalid, smiles, inchi, inchikey, subs.link \
                FROM public.substance AS subs \
                INNER JOIN public.source AS s ON s.id = subs.sourceid \
                WHERE s.id = %s"
    else:
        cmd = "SELECT externalid, smiles, inchi, inchikey, subs.link \
                FROM public.substance AS subs \
                INNER JOIN public.source AS s ON s.id = subs.sourceid \
                WHERE s.id = %s AND smiles IS NOT NULL"

    df = pd.read_sql(cmd, con=conn)

    return df

def deleteSource(conn, sourceID):
    """
    Deletes and source and all its associated substances.
    """
    curs = conn.cursor()
    cmd = "DELETE FROM public.source WHERE id=%s"
    curs.execute(cmd, (sourceID,))
    conn.commit()