import sys, datetime, os, pickle
import pandas as pd
import psycopg2
import pubchempy as pcp
import time

from phitools import moleculeHelper as mh
from standardiser import process_smiles as ps
from compoundDB import querytools as qt

from typing import *

from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.Crippen import MolLogP

annClassification_file = "annD.pkl"
with open(os.path.join(os.path.dirname(__file__), "data",  annClassification_file), 'rb') as annC_fh:
    annD = pickle.load(annC_fh)

def openconnection(host='gea', dbname='compounds', user='postgres', password=''):
    """
        Opens a psycopg2 connection to the PostgreSQL DB and returns the connection object.
        Arguments:
          - host
          - dbname
          - user
          - password
    """
    conn = psycopg2.connect(host=host, dbname=dbname, user=user, password=password)
    curs = conn.cursor()
    curs.execute('create extension if not exists rdkit;')

    return conn

def standardiseAnnotation(ann):
    """
    From a given toxicity annotation return a dictionary with:
        - Corrected annotation: Standardised annotation.
        - Category: Toxicity group (CMR, PBT, Physical, etc.)
        - Type: Annotation type (Confirmed, Suspected, Negative)
        - General annotation: For certain annotations, a more general classification of the toxicity. i.e. for Carc. 1A the corresponding General Annotation would be 'Carcinogenic'. This allows to group certain types of toxicity with more detail than the category (in this case CMR).
    """
    if ann in annD:
        result = annD[ann]
    else:
        result = {'Category': '',
                    'Corrected annotation': ann,
                    'GeneralAnnotation': 'Mutagen',
                    'Type': 'Confirmed'}

    return annD[ann]

def addAnnotation(conn, subsID, ann, annType=None, annCategory=None, generalAnn=None, sourceID=None):
    """
        Add the annotation provided for a given substance to the corresponding table.
        Arguments:
          - conn: psycopg2 connection to the database.
          - subsID: id for the partent subsance from the 'substance' table.
          - ann: Toxicity annotation. If it cannot be standadized no category or general annotation will be assigned.
          - annType: Optional. Provide it to overwrite the annotation type from the annD.
          - annCategory: Optional. Provide it to overwrite the annotation category from the annD.
          - generalAnn: Optional. Provide it to overwrite the general annotation from the annD.
    """
    ann = ann.strip().strip('\"').strip()
    try:
        d = standardiseAnnotation(ann)
    except:
        std_ann = ann
        if not annType:
            annType = 'Confirmed'
        annCategory = ''
        generalAnn = ''
    else:
        std_ann = d['Corrected annotation']
        if not annType:
            annType = d['Type']
        if not annCategory:
            annCategory = d['Category']
        if not generalAnn:
            generalAnn = d['GeneralAnnotation']
    curs = conn.cursor()

    cmd = "SELECT id FROM annotation WHERE annotation = %s;"
    curs.execute(cmd, (std_ann,))
    annID = curs.fetchone()
    conn.commit()
    
    if annID is None:      
        cmd = "INSERT INTO annotation (annotation, general, category)\
               VALUES (%s, %s, %s)"
        curs.execute(cmd, (std_ann, generalAnn, annCategory))
        conn.commit()
        cmd = "SELECT currval('annotation_id_seq');"
        curs.execute(cmd)
        annID = curs.fetchone()[0]
        conn.commit()
    else:
        annID = annID[0]

    if sourceID is None:
        cmd = "SELECT sourceid FROM substance \
                WHERE id = %s;"
        curs.execute(cmd, (subsID,))
        sourceID = curs.fetchone()[0]
        conn.commit()

    cmd = "INSERT INTO subs_ann (subsid, annid, sourceid, original_annotation, type)\
           SELECT %s, %s, %s, %s, %s \
           WHERE NOT EXISTS (SELECT subsid, annid, sourceid, original_annotation, type \
               FROM subs_ann \
               WHERE subsid= %s AND annid= %s)"
    curs.execute(cmd, (subsID, annID, sourceID, ann, annType, subsID, annID))
    conn.commit()

def addSynonyms(conn, subsID, synD):
    """
        Add all synonyms provided for a given substance.
        Arguments:
          - conn: psycopg2 connection to the database.
          - subsID: id for the partent subsance from the 'substance' table.
          - synD: Dictionary with the substance synonyms' and their types {type: synonym}.
    """
    curs = conn.cursor()
    cmd = "INSERT INTO public.synonym (subsid, type, name)\
           SELECT %s, %s, %s\
           WHERE NOT EXISTS (SELECT subsid, type, name\
               FROM public.synonym \
               WHERE subsid= %s AND type= %s \
               AND name= %s)"
    curs.executemany(cmd, [(subsID,syntype,syn,subsID,syntype,syn) for syntype in synD for syn in synD[syntype]])
    conn.commit()

def addSource(conn, sourceName, version= None, description= None, link= None, \
              date= datetime.datetime.now().date()):
    """
        Insert a source to the DB.
        Arguments:
          - conn: psycopg2 connection to the database.
          - sourceName: Source name.
          - version: Optional. Source's version (default: None). If None, if another source with the same name \
            exists the version will be increased by one, else the version will be 1.
          - description: Optional. Verbose description of the source (default: None).
          - link: Optional. Link to the source's home page (default: None).
          - date: Optional. Date of source's creation (default: Today).

        Returns source id from the 'source' table.
    """
    curs = conn.cursor()
    # Check if version is provided, otherwise generate it
    if not version:
        cmd = "SELECT id FROM source WHERE name = %s ORDER BY version DESC;"
        curs.execute(cmd, (sourceName,))
        oldVersion = curs.fetchone()
        conn.commit()
        if oldVersion is None:
            version = '1'
        else:
            version = int(oldVersion[0]+1)
    version = str(version)
    
    cmd = "SELECT id FROM source WHERE name = %s AND version = %s;"
    curs.execute(cmd, (sourceName, version))
    sourceID = curs.fetchone()
    conn.commit()
    
    if sourceID is None:
        cmd = "INSERT INTO source (name, version, description, link, added)\
               VALUES (%s, %s, %s, %s, %s)"
        curs.execute(cmd, (sourceName, version, description, link, date))
        conn.commit()
        cmd = "SELECT currval('source_id_seq');"
        curs.execute(cmd)
        sourceID = curs.fetchone()[0]
        conn.commit()
    else:
        sourceID = sourceID[0]
    
    ### Check latest. Set to False if the version is not the last one. Otherwise, True

    cmd= "UPDATE source\
            SET latest= CASE\
                WHEN version = (SELECT max(version) FROM source WHERE name='{}')\
                THEN True\
                ELSE False\
                END\
            where name='{}';"
    curs.execute(cmd.format(sourceName,sourceName))
    conn.commit()

    return sourceID

def addEmptyCompound(conn, subsID, smiles):
    """
        Insert an empty compound to the DB, in case a compound's smiles is available but it can't be converted into an RDkit mol object.
        Also inserts an instance to the 'subs_cmpd' table with the parent's substance id and the compound id.
        Arguments:
          - conn: psycopg2 connection to the database.
          - subsID: id for the partent subsance from the 'substance' table.
          - smiles: Smiles string for the compound.

        Returns compound id from the 'compound' table.
    """
    curs = conn.cursor()
    cmd = "SELECT id FROM compound WHERE smiles = %s;"
    curs.execute(cmd, (smiles,))
    cmpdID = curs.fetchone()
    conn.commit()
        
    if cmpdID is None:  
        # /// Where not exists
        cmd = "INSERT INTO compound (smiles)\
           VALUES (%s)"
        curs.execute(cmd, (smiles,))
        conn.commit()
        cmd = "SELECT currval('compound_id_seq');"
        curs.execute(cmd)
        cmpdID = curs.fetchone()[0]
        conn.commit()
    else:
        cmpdID = cmpdID[0]
            
    cmd = "INSERT INTO subs_cmpd (subsid, cmpdid)\
       VALUES (%s, %s)"
    curs.execute(cmd, (subsID, cmpdID))
    conn.commit()
    
    return cmpdID

def addCompound(conn, subsID, smiles= None, mol= None, ismetal=False):
    """
        Insert one standardised compound to the DB. Must provide either a smiles string or an RDKit mol object.
        Also inserts an instance to the 'subs_cmpd' table with the parent's substance id and the compound id.
        Arguments:
          - conn: psycopg2 connection to the database.
          - subsID: id for the partent subsance from the 'substance' table.
          - smiles: Optional. Smiles string for the compound (default: None).
          - mol: Optional. RDkit mol object (default: None).
          - ismetal: Optional. Boolean flag indicating whether the compound is a metal ion (default: False).

        Returns compound id from the 'compound' table.
    """
    curs = conn.cursor()
    if smiles is None and mol is None:
        cmpdID = addEmptyCompound(conn, subsID, smiles)
    elif mol is None:
        try:
            mol = Chem.MolFromSmiles(smiles)
        except:
            print (smiles)
            cmpdID = addEmptyCompound(conn, subsID, smiles)
            
    if mol is None:
        # Molecule object has been created, but only as an empty instnce
        cmpdID = addEmptyCompound(conn, subsID, smiles)            
    else:
        inchi = Chem.MolToInchi(mol)
        inchikey = Chem.InchiToInchiKey(inchi)

        cmd = "SELECT id FROM compound WHERE inchikey = %s;"
        curs.execute(cmd, (inchikey,))
        cmpdID = curs.fetchone()
        conn.commit()

        if cmpdID is None:
            if smiles is None:
                smiles = Chem.MolToSmiles(mol)
            mw = MolWt(mol)
            logp = MolLogP(mol)

            cmd = "INSERT INTO compound (smiles, inchi, inchikey, mol, mw, logp, ecfp4, \
                                         fcfp4, ismetalion)\
                   VALUES (%s, %s, %s, mol_from_smiles(%s), %s, %s, \
                           morgan_fp(mol_from_smiles(%s), 4), \
                           featmorgan_fp(mol_from_smiles(%s), 4), %s)"
            curs.execute(cmd, (smiles, inchi, inchikey, smiles, mw, logp, \
                              smiles, smiles, ismetal))
            conn.commit()
            cmd = "SELECT currval('compound_id_seq');"
            curs.execute(cmd)
            cmpdID = curs.fetchone()[0]
            conn.commit()
        else:
            cmpdID = cmpdID[0]
    
    cmd = "SELECT subsid, cmpdid FROM subs_cmpd WHERE subsid = %s AND cmpdid = %s;"
    curs.execute(cmd, (subsID, cmpdID))
    subs_cmpid = curs.fetchone()
    conn.commit()
    
    if subs_cmpid is None:
        cmd = "INSERT INTO subs_cmpd (subsid, cmpdid)\
        VALUES (%s, %s)"
        curs.execute(cmd, (subsID, cmpdID))
        conn.commit()
    
    return cmpdID

def addEmptySubstance(conn, sourceID, extID, link= None):
    """
        Insert an empty substance to the DB, in case a substance is present in the source but it doesn't have a structure or it can't be converted into an RDkit mol object.
        Arguments:
          - conn: psycopg2 connection to the database.
          - sourceID: id for the source of origin from the 'source' table.
          - extID: id for the subsance in the source of origin.
          - link: Optional. Link to information on the substance (default: None).

        Returns a tupple with the substance id from the 'substance' table and None, since no RDKit mol object is available.
    """
    curs = conn.cursor()
    cmd = "SELECT id FROM substance WHERE sourceid = %s \
            AND externalid = %s;"
    
    curs.execute(cmd, (sourceID, extID))
    subsID = curs.fetchone()
    conn.commit()

    if not subsID:
        curs = conn.cursor()
        cmd = "INSERT INTO substance (sourceid, externalid, link)\
        VALUES (%s, %s, %s)"
        curs.execute(cmd, (sourceID, extID, link))
        conn.commit()
        cmd = "SELECT currval('substance_id_seq');"
        curs.execute(cmd)
        subsID = curs.fetchone()[0]
        conn.commit()
    else:
        subsID = subsID[0]
    
    return (subsID, None)

def addSubstance(conn, sourceID, extID, smiles= None, mol= None, link= None):
    """
        Standardise and insert a substance to the DB. Must provide either a smiles string or an RDKit mol object.
        Arguments:
          - conn: psycopg2 connection to the database.
          - sourceID: id for the source of origin from the 'source' table.
          - extID: id for the subsance in the source of origin.
          - smiles: Optional. Smiles string for the substance (default: None).
          - mol: Optional. RDkit mol object (default: None).
          - link: Optional. Link to information on the substance (default: None).

        Returns a tupple with the substance id from the 'substance' table and the RDKit mol object of the standardised substance if available.
    """
    curs = conn.cursor()    
            
    cmd = "SELECT id FROM substance WHERE sourceid = %s \
            AND externalid = %s;"
    curs.execute(cmd, (sourceID, extID))
    subsID = curs.fetchone()
    conn.commit()

    if smiles is None and mol is None:
        (subsID, mol) = addEmptySubstance(conn, sourceID, extID, link)
    elif mol is None:
        try:
            mol = Chem.MolFromSmiles(smiles)
        except:
            (subsID, mol) = addEmptySubstance(conn, sourceID, extID, link)

    if not subsID:
        if mol is None:
            # Molecule object has been created, but only as an empty instnce
            (subsID, mol) = addEmptySubstance(conn, sourceID, extID, link)
        else:
            inchi = Chem.MolToInchi(mol)
            inchikey = Chem.InchiToInchiKey(inchi)
            if smiles is None:
                smiles = Chem.MolToSmiles(mol)

            cmd = "INSERT INTO substance (sourceid, externalid, smiles, mol, inchi, inchikey, link)\
            VALUES (%s, %s, %s, mol_from_smiles(%s), %s, %s, %s)"
            curs.execute(cmd, (sourceID, extID, smiles, smiles, inchi, inchikey, link))
            conn.commit()
            cmd = "SELECT currval('substance_id_seq');"
            curs.execute(cmd)
            subsID = curs.fetchone()[0]
            conn.commit()

            try:
                stdD = ps.std(mol, returnMetals=True)
                for cmpd_smiles in stdD:
                    (cmpd_mol, ismetal, passed, errmessage) = stdD[cmpd_smiles]
                    if ismetal or passed or errmessage == 'Multiple non-salt/solvate components':
                        cmpdID = addCompound(conn, subsID, smiles= cmpd_smiles, ismetal=ismetal)
            except:
                cmpdID = addCompound(conn, subsID)

    else:
        subsID = subsID[0]
            
    return (subsID, mol)

def addSubstanceFromSmilesFile(conn, sourceID, fname, 
                            extIDindex= None, extIDfield= None, 
                            smilesIndex= 1, smilesField= 'smiles', 
                            ann= None,
                            annIndex= None, annField= None, 
                            annType = None, annTypeField= None, 
                            annTypeIndex= None,
                            linkIndex= None, linkField= None, 
                            synonymsIndices= None, synonymsFields= None, header= True):
    """
        Process a text file with smiles strings of substances from a given source.
        Arguments:
          - conn: psycopg2 connection to the database.
          - sourceID: id for the source of origin from the 'source' table.
          - fname: Input file name.
          - extIDindex: Optional. Index of the column containing the id of the substance in the source of origin (default: None). If None, an id will be generated with a substance counter.
          - extIDfield: Optional. Name of the header of the column containing the substance id (default: None). If None, an id will be generated with a substance counter.
          - smilesIndex: Optional. Index of the column containing the substance's smiles string (default: 1). 
          - smilesField: Optional. Name of the header of the column containing the substance's smiles string (default: 'smiles').
          - annIndex: Optional. Index of the column containing the substance's annotation (default: 1). 
          - annField: Optional. Name of the header of the column containing the substance's annotation (default: 'Annotation').
          - annType: Optional. Type of all the annotations of the substances for this source (default: None). It should be one of the following: 'Confirmed', 'Suspected', 'Negative'.
          - linkIndex: Optional. Index of the column containing a link to the substance information page (default: None). 
          - linkField: Optional. Name of the header of the column containing a link to the substance information page (default: None).
          - synonymsIndices: Optional. List of indices of the column(s) containing synonyms of the substance (default: None). Synonym type will be 'Name'.
          - synonymsFields: Optional. List of name(s) of the header of the column(s) containing synonyms of the substance (default: None). 
          - header: Boolean indicating if the file has a header (default: False).
    """         
    with open(fname) as f:

        if header:
            header = f.readline().rstrip().split('\t')
            if extIDfield:
                extIDindex = header.index(extIDfield)

            if smilesField:
                smilesIndex = header.index(smilesField)

            if annField:
                annIndex = header.index(annField)

            if annTypeField:
                annTypeIndex = header.index(annTypeField)

            if linkField:
                linkIndex = header.index(linkField)

            if synonymsFields:
                synTypes = synonymsFields
                synonymsIndices = []
                for t in synTypes:
                    synonymsIndices.append(header.index(t))

            elif synonymsIndices:
                synTypes = []
                for i in synonymsIndices:
                    synTypes.append(header[i].strip())

        else:
            # If the file has no header and some columns contain synonyms, 
            # set each synonym type to 'Name'
            if synonymsIndices:
                synonymsIndices = synonymsIndices
                synTypes = []
                for i in synonymsIndices:
                    synTypes.append('Name')

        molcount = 0
        for line in f:
            molcount += 1
            fields = line.rstrip().split('\t')

            try:
                smi = fields[smilesIndex]
            except:
                smi = None

            if extIDindex is None:
                # No field with the ID of the substance in the source of origin 
                # has been provided so one will be generated.
                extID = 'mol%0.8d'%molcount
            else:
                try:
                    extID = fields[extIDindex]
                except:
                    extID = 'mol%0.8d'%molcount

            if not linkIndex: link = None
            else: link= fields[linkIndex]

            # Add the subsance
            try:
                mol = Chem.MolFromSmiles(smi)
            except:
                (subsID, mol) = addEmptySubstance(conn, sourceID, extID, link)
            else:
                (subsID, mol) = addSubstance(conn, sourceID, extID= extID, smiles= smi, 
                                                mol= mol, link= link)
                
            # Add synonyms
            synD = {'ExternalID': set([extID])}
            if synonymsIndices:
                for i in range(len(synonymsIndices)):
                    sindex = synonymsIndices[i]
                    stype = synTypes[i]
                    if len(fields) <= sindex: continue
                    syn = fields[sindex].strip()
                    if syn == 'N/A' or syn == '': continue
                    if stype not in synD:
                        synD[stype] = set([syn])
                    else:
                        synD.add(syn)
            addSynonyms(conn, subsID, synD)

            # Add annotations
            annotation = None
            if not ann:
                try:
                    annotation = fields[annIndex]
                except:
                    annotation = None
                else:
                    if annotation == '****':
                        annotation = None
            else:
                annotation = ann
            
            if annotation:
                if annTypeIndex is not None:
                    try:
                        ann_type = fields[annTypeIndex]
                    except:
                        ann_type = None
                else:
                    ann_type = annType
                addAnnotation(conn, subsID, annotation, ann_type)

def processCASfromPandas(cas_:str) -> Union[str,list]:
    """
        Process CAS field from the dataframe, either getting a clean string of one CAS or
        a list of CAS

        Arguments:
            - cas_: content from CAS field in the dataframe
        
        Returns:
            - cas_proc: is either a string or a list containing the processed CAS related to
            the substance
    """

    preproc_cas = cas_.split(',')

    if len(preproc_cas) > 1:
        cas_proc = [cas_element.replace("'","").strip() for cas_element in preproc_cas]
        if '-' in cas_proc:
            cas_proc.remove('-')
    else:
        cas_proc = str(preproc_cas).strip('[').strip(']').replace("'","")

    return cas_proc

def getresultsfromPubChem(casrn: str) -> list:
    """
        Get results from PubChem using the API. Checks for HTTPErrors and puts to sleep the
        code in order to avoid timeout error or bad gateaway error.

        Arguments:
            - casrn: CAS number
        
        Returns:
            - results: list with names from PubChem
        
        TODO: make counter up to 5 requests then time.sleep for 2 seconds. Catch exception.
    """
    try:
        results = pcp.get_compounds(casrn, 'name')
    except:
        results = None
    
    return results

def rerieveSMILESfromMHorPubChem(casrn: str, conn: psycopg2.extensions.connection) -> Optional[str]:
    """
        If molecule helper can't get SMILES, it uses PubChem API to check for SMILES
        string from CAS number

        Arguments:
            - casrn: CAS number
            - conn: connector to Database
        
        Return:
            - smi: canonical SMILES from CAS in PubChem
    """
    smi = mh.resolveCAS(cas= casrn, conn= conn)
    
    if smi:
        return smi
    else:
        struc_list = []
        results = getresultsfromPubChem(casrn)
        if results:
            for compound in results:
                smiles = compound.canonical_smiles
                struc_list.append(smiles)
            if len(struc_list) > 1:
                struc_list = list(set(struc_list))
            smi = struc_list[0]
        else:
            smi = None
        return smi

def getSMILESfromCAS(row_df: pd.DataFrame, casrn:str, updatesmi: Optional[str], extID: str, sourceID: int,
                    conn: psycopg2.extensions.connection,
                    curs: psycopg2.extensions.cursor, link: str, smilesField: Optional[str] = None) -> str:
    """
        Retrieves SMILES from CAS

        Arguments:
            - casrn: CAS number
            - updatesmi:
            - extID:
            - sourceID:
            - conn:
            - curs:
            - link:
        
        Return:
            - smi: SMILES string from CAS
            - subsID: substance id from the db
    """

    if not updatesmi:
        cmd = "SELECT id \
                FROM substance \
                WHERE externalid = %s AND sourceid = %s;"
        curs.execute(cmd, (extID, sourceID))
        subsID = curs.fetchone()
        conn.commit()

        if subsID is None:
            subsfound = False
        else:
            subsID = subsID[0]
            subsfound = True
    else:
        subsfound = None

    if not subsfound:
        if smilesField:
            smi = row_df[smilesField]                
        elif casrn:
            smi = rerieveSMILESfromMHorPubChem(casrn, conn)
        else:
            smi = None
        
        # Add the subsance
        try:
            mol = Chem.MolFromSmiles(smi)
        except:
            (subsID, mol) = addEmptySubstance(conn, sourceID, extID, link)
        else:
            (subsID, mol) = addSubstance(conn, sourceID, extID= extID, smiles= smi, 
                                        mol= mol, link= link)

    return subsID

def addSynonymsFromCASpandasDf(row_df: pd.Series, synonym_fields: list, extID: str, 
                                subsID: int, conn: psycopg2.extensions.connection):
    """
        Adds synonyms to the database

        Arguments:
            - row_df: Dataframe row with substance information
            - synonym_fields: column names where Synonyms are
            - extID: CAS number
            - subsID: substance ID given in the database
    """

    synD = {'ExternalID': set([extID])}
    if synonym_fields:
        for synonym in synonym_fields:
            syn = row_df[synonym].strip()
            if syn == 'N/A' or syn == '': 
                continue
            if synonym not in synD:
                synD[synonym] = set([syn])
            else:
                synD[synonym].add(syn)
    
    addSynonyms(conn, subsID, synD)

def addAnnotationsFromPandasDf(row_df: pd.Series, annotation_field: str, subsid: int, 
                            annotation_type: str, conn:psycopg2.extensions.connection):
    """
        Add annotations from pandas dataframe

        Arguments:
            - row_df: dataframe row
            - annotation_field: column in the dataframe with annotation information
            - subsid: substance id from the database
            - annotation_type: type of annotation to check
            - conn: connection to the database
    """

    try:
        annotation = row_df[annotation_field]
    except:
        annotation = None
    else:
        if annotation == '****':
            annotation = None
    
    if annotation:
        if annotation_type is not None:
            try:
                ann_type = row_df[annotation_type]
            except:
                ann_type = None
        else:
            ann_type = annotation_type
        
        addAnnotation(conn, subsid, annotation, ann_type)

def addSubstanceFromPandasDf(conn, sourceID, dataframe, 
                        extIDindex= None, extIDfield= None, 
                        CASindex= None, CASfield= None, 
                        ann= None, smilesField= None,
                        annIndex= None, annField= None, 
                        annType = None, 
                        annTypeField= None, annTypeIndex= None, 
                        linkIndex= None, linkField= None, 
                        synonymsIndices= None, synonymsFields= None,
                        updatesmi= False):
    """
        Extracts the necessary fields from a Pandas Dataframe that comes from
        processing the source file.

        Arguments:
            - conn: psycopg2 connection to the database.
            - sourceID: id for the source of origin from the 'source' table.
            - dataframe: input dataframe
            - extIDindex: Optional. Index of the column containing the id of the substance in the source of origin (default: None). If None, the CAS number will be used as external ID.
            - extIDfield: Optional. Name of the header of the column containing the substance id (default: None). If None, the CAS number will be used as external ID.
            - CASindex: Optional. Index of the column containing the substance's CAS number (default: None). 
            - CASfield: Optional. Name of the header of the column containing the substance's CAS number (default: None).
            - smilesField: Optional. Name of the header of the column containing the substance's smiles string (default: None).
            - annIndex: Optional. Index of the column containing the substance's annotation (default: None). 
            - annField: Optional. Name of the header of the column containing the substance's annotation (default: None).
            - annType: Optional. Type of all the annotations of the substances for this source (default: None). It should be one of the following: 'Confirmed', 'Suspected', 'Negative'.
            - annTypeField: Optional. Name of the header of the column containing all the annotations of the substances for this source (default: None). It should be one of the following: 'Confirmed', 'Suspected', 'Negative'.
            - annTypeIndex: Optional. Index of the column containing all the annotations of the substances for this source (default: None). It should be one of the following: 'Confirmed', 'Suspected', 'Negative'.
            - linkIndex: Optional. Index of the column containing a link to the substance information page (default: None). 
            - linkField: Optional. Name of the header of the column containing a link to the substance information page (default: None).
            - synonymsIndices: Optional. List of indices of the column(s) containing synonyms of the substance (default: None). Synonym type will be 'Name'.
            - synonymsFields: Optional. List of name(s) of the header of the column(s) containing synonyms of the substance (default: None). 
            - updatesmi: Even if the substance is already in the DB, force lookup of its structure.
    """

    curs = conn.cursor()

    if extIDfield is None:
        extIDfield = CASfield

    molcount = 0
    for i, rows in dataframe.iterrows():
        molcount += 1
        try:
            CAS = processCASfromPandas(rows[CASfield])
        except:
            CAS = None
        else:
            if CAS == '-':
                CAS = None

        if extIDfield is None:
            # No field with the ID of the substance in the source of origin 
            # has been provided so one will be generated.
            extID = CAS
        else:
            try:
                extID = rows[extIDfield]
            except:
                extID = CAS

        if extIDfield is None:
            # No field with the ID of the substance in the source of origin 
            # has been provided so one will be generated.
            extID = 'mol%0.8d'%molcount

        if not linkField: 
            link = None
        else: 
            link = dataframe[linkField]
        
        # Get smiles from CAS
        if isinstance(CAS, list):
            subid_list = []
            for casrn,extid in zip(CAS,extID):
                subsID = getSMILESfromCAS(rows, casrn, updatesmi, extid, sourceID, conn, 
                                        curs, link, smilesField)
                subid_list.append((extid,subsID))
        elif isinstance(CAS, str):
            subsID = getSMILESfromCAS(rows, CAS, updatesmi, extID, sourceID, conn, curs, 
                                    link, smilesField)

        # Add synonyms and annotations
        
        if isinstance(extID, list):
            for extid, subsid in subid_list:
                addSynonymsFromCASpandasDf(rows, synonymsFields, extid, subsid, conn)
                addAnnotationsFromPandasDf(rows, annField, subsid, annType, conn)
        elif isinstance(extID, str):
            addSynonymsFromCASpandasDf(rows, synonymsFields, extID, subsID, conn)
            addAnnotationsFromPandasDf(rows, annField, subsID, annType, conn)

def addSubstanceFromCASFile(conn, sourceID, fname, 
                        extIDindex= None, extIDfield= None, 
                        CASindex= None, CASfield= None, 
                        ann= None,
                        annIndex= None, annField= None, 
                        annType = None, 
                        annTypeField= None, annTypeIndex= None, 
                        linkIndex= None, linkField= None, 
                        synonymsIndices= None, synonymsFields= None,
                        updatesmi= False, header= True):
    """
        Process a text file with smiles strings of substances from a given source.
        Arguments:
          - conn: psycopg2 connection to the database.
          - sourceID: id for the source of origin from the 'source' table.
          - fname: Input file name.
          - extIDindex: Optional. Index of the column containing the id of the substance in the source of origin (default: None). If None, the CAS number will be used as external ID.
          - extIDfield: Optional. Name of the header of the column containing the substance id (default: None). If None, the CAS number will be used as external ID.
          - CASindex: Optional. Index of the column containing the substance's CAS number (default: 1). 
          - CASfield: Optional. Name of the header of the column containing the substance's CAS number (default: 'smiles').
          - annIndex: Optional. Index of the column containing the substance's annotation (default: 1). 
          - annField: Optional. Name of the header of the column containing the substance's annotation (default: 'Annotation').
          - annType: Optional. Type of all the annotations of the substances for this source (default: None). It should be one of the following: 'Confirmed', 'Suspected', 'Negative'.
          - annTypeField: Optional. Name of the header of the column containing all the annotations of the substances for this source (default: None). It should be one of the following: 'Confirmed', 'Suspected', 'Negative'.
          - annTypeIndex: Optional. Index of the column containing all the annotations of the substances for this source (default: None). It should be one of the following: 'Confirmed', 'Suspected', 'Negative'.
          - linkIndex: Optional. Index of the column containing a link to the substance information page (default: None). 
          - linkField: Optional. Name of the header of the column containing a link to the substance information page (default: None).
          - synonymsIndices: Optional. List of indices of the column(s) containing synonyms of the substance (default: None). Synonym type will be 'Name'.
          - synonymsFields: Optional. List of name(s) of the header of the column(s) containing synonyms of the substance (default: None). 
          - updatesmi: Even if the substance is already in the DB, force lookup of its structure.
          - header: Boolean indicating if the file has a header (default: False).
    """
    curs = conn.cursor()
    with open(fname) as f:
        if header: 
            header = f.readline().rstrip().split('\t')
            if extIDfield:
                extIDindex = header.index(extIDfield)

            if CASfield:
                CASindex = header.index(CASfield)

            if annField:
                annIndex = header.index(annField)

            if annTypeField:
                annTypeIndex = header.index(annTypeField)

            if linkField:
                linkIndex = header.index(linkField)

            if synonymsFields:
                synTypes = synonymsFields
                synonymsIndices = []
                for t in synTypes:
                    synonymsIndices.append(header.index(t))

            elif synonymsIndices:
                synTypes = []
                for i in synonymsIndices:
                    synTypes.append(header[i].strip())

        else:
            # If the file has no header and some columns contain synonyms, 
            # set each synonym type to 'Name'
            if synonymsIndices:
                synonymsIndices = synonymsIndices
                synTypes = []
                for i in synonymsIndices:
                    synTypes.append('Name')

        if extIDindex is None:
            extIDindex = CASindex

        molcount = 0
        for line in f:
            molcount += 1
            fields = line.rstrip().split('\t')
            try:
                CAS = fields[CASindex]
            except:
                CAS = None
            else:
                if CAS == '-':
                    CAS = None

            if extIDindex is None:
                # No field with the ID of the substance in the source of origin 
                # has been provided so one will be generated.
                extID = CAS
            else:
                try:
                    extID = fields[extIDindex]
                except:
                    extID = CAS
            if extIDindex is None:
                # No field with the ID of the substance in the source of origin 
                # has been provided so one will be generated.
                extID = 'mol%0.8d'%molcount

            if not linkIndex: link = None
            else: link= fields[linkIndex]

            # Get smiles from CAS
            if not updatesmi:
                cmd = "SELECT id \
                        FROM substance \
                        WHERE externalid = %s AND sourceid = %s;"
                curs.execute(cmd, (extID, sourceID))
                subsID = curs.fetchone()
                conn.commit()

                if subsID is None:
                    subsfound = False
                else:
                    subsID = subsID[0]
                    subsfound = True
            else:
                subsfound = None

            if not subsfound:                
                if CAS:
                    smi = mh.resolveCAS(cas= CAS, conn= conn)
                else:
                    smi = None

                # Add the subsance
                try:
                    mol = Chem.MolFromSmiles(smi)
                except:
                    (subsID, mol) = addEmptySubstance(conn, sourceID, extID, link)
                else:
                    (subsID, mol) = addSubstance(conn, sourceID, extID= extID, smiles= smi, mol= mol, link= link)
                
            # Add synonyms
            synD = {'ExternalID': set([extID])}
            if synonymsIndices:
                for i in range(len(synonymsIndices)):
                    sindex = synonymsIndices[i]
                    stype = synTypes[i]
                    if len(fields) <= sindex: continue
                    syn = fields[sindex].strip()
                    if syn == 'N/A' or syn == '': continue
                    if stype not in synD:
                        synD[stype] = set([syn])
                    else:
                        synD.add(syn)
            addSynonyms(conn, subsID, synD)

            # Add annotations
            annotation = None
            if not ann:
                try:
                    annotation = fields[annIndex]
                except:
                    annotation = None
                else:
                    if annotation == '****':
                        annotation = None
            else:
                annotation = ann

            if annotation:
                if annTypeIndex is not None:
                    try:
                        ann_type = fields[annTypeIndex]
                    except:
                        ann_type = None
                else:
                    ann_type = annType
                addAnnotation(conn, subsID, annotation, ann_type)

def addSubstanceSDFile(conn, sourceID, fname, extIDfield= None, 
                    linkField= None, ann= None,
                    annField= None, annType = None, annTypeField= None, synonymsFields= None):
    """
        Process an SD file of substances from a given source.
        Arguments:
          - conn: psycopg2 connection to the database.
          - sourceID: id for the source of origin from the 'source' table.
          - fname: Input file name.
          - extIDfield: Optional. Name of the field containing the substance id (default: None). If None, an id will be generated with a substance counter.
          - synonymsFields: Optional. List of name(s) of the field(s) containing synonyms of the substance (default: None). 
          - annField: Optional. Name of the field with the substance's annotation (default: 'Annotation').
          - annType: Optional. Type of all the annotations of the substances for this source (default: None). It should be one of the following: 'Confirmed', 'Suspected', 'Negative'.
          - annTypeField: Optional. Name of the field with the type of the substance's annotation (default: None). 
    """
    suppl = Chem.SDMolSupplier(fname)
    molcount = 0
    for mol in suppl:
        molcount += 1
        if mol is None:
            extID = mh.getNameFromEmpty(suppl, molcount-1, extIDfield)
            (subsID, mol) = addEmptySubstance(conn, sourceID, extID, link)
            
        else:                
            extID = mh.getName(mol, molcount, extIDfield)                    
            link = None
            if linkField:
                try:
                    link = mol.GetProp(linkField)
                except:
                    pass
            (subsID, mol) = addSubstance(conn, sourceID, extID= extID, mol= mol, link= link)
        
        # Add synonyms
        synD = {'ExternalID': set([extID])}
        if synonymsFields:
            for synType in synonymsFields:
                try:
                    syn = mol.GetProp(synType)
                    if syn != 'N/A':
                        if synType not in synD:
                            synD[synType] = set([syn])
                        synD[synType].add(syn)
                except:
                    pass
        addSynonyms(conn, subsID, synD)
        
        # Add annotations
        annotation = None
        if not ann:
            if annField:
                try:
                    ann = mol.GetProp(annField)
                except:
                    annotation = None
        else:
            annotation = ann

        if annotation:
            if annTypeIndex is not None:
                try:
                    ann_type = fields[annTypeIndex]
                except:
                    ann_type = None
            else:
                ann_type = annType
                
            ann_type = annType
            if annTypeField:
                try:
                    ann_type = mol.GetProp(annField)
                except:
                    pass
            addAnnotation(conn, subsID, ann, ann_type)

def addSubstanceFromQuery(conn, sourceID, cmd, host='gea', dbname='chembl_23', user='postgres', \
                          password='', extIDf= 'ID', smilesF= 'smiles', \
                          linkF= None):
    """
        Add substances from a given source from a PostgeSQL DB. 
        Arguments:
          - conn: psycopg2 connection to the database where the substances will be input.
          - sourceID: id for the source of origin from the 'source' table.
          - cmd: Postgres statement to extract the substances.
          Connection details for the DB containing the substances:
            - host: 
            - dbname
            - user
            - password
          - fname: Input file name.
          - ftype: Input file type (sdf / smi).
          - extIDf: Name of the DB field containing the id of the substance in the source of origin (default: ID). 
          - smilesF: Name of the DB field containing the smiles string (default: smiles). 
          - linkF: Optional. Name of the DB field containing the link to the substance information page (default: None). 
    """
    
    # Open connection to the query DB
    qconn = psycopg2.connect(host=host, dbname=dbname, user=user, password=password)
    df = pd.read_sql(cmd, con=qconn)
    for index, row in df.iterrows():
        extID = row[extIDf]
        smi = row[smilesF]
        if not linkF: link = None
        else: link= row[linkF]
        try:
            mol = Chem.MolFromSmiles(smi)
        except:
            addEmptySubstance(conn, sourceID, extID, link)
            continue
        addSubstance(conn, sourceID, extID= extID, smiles= smi, \
                     mol= mol, link= link)

def addSynonymsFromFile(conn, fname, sourceID= None, sourceName= None, version= None,extIDindex= None, extIDfield= None, synonymsIndices= None, synonymsFields= None, header= False):
    """
    Add synonyms for substances already in the DB.
    Arguments:
        - conn: psycopg2 connection to the database.
        - fname: Input file name.
        - sourceID: Optional. Internal source ID if avaialble. Otherwise, it will be retrieved from he source name and version.
        - sourceName: Source name.
        - version: Optional. Source's version (default: None). If None, 
        the last version will be used.
        - extIDindex: Optional. Index of the column containing the id of the substance in the source of origin (default: None). If None, an id will be generated with a substance counter.
        - extIDfield: Optional. Name of the header of the column containing the substance id (default: None). If None, an id will be generated with a substance counter.
        - synonymsIndices: Optional. List of indices of the column(s) containing synonyms of the substance (default: None). Synonym type will be 'Name'.
        - synonymsFields: Optional. List of name(s) of the header of the column(s) containing synonyms of the substance (default: None). 
        - header: Boolean indicating if the file has a header (default: False).
    """
    if sourceID is None:
        sourceID = getSourceID(conn, sourceName, version)
    with open(fname) as f:
        if header: 
            header = f.readline().rstrip().split('\t')
            if extIDfield:
                extIDindex = header.index(extIDfield)

            if synonymsFields:
                synTypes = synonymsFields
                synonymsIndices = []
                for t in synTypes:
                    synonymsIndices.append(header.index(t))

            elif synonymsIndices:
                synTypes = []
                for i in synonymsIndices:
                    synTypes.append(header[i].strip())

        else:
            # If the file has no header and some columns contain synonyms, 
            # set each synonym type to 'Name'
            if synonymsIndices:
                synonymsIndices = synonymsIndices
                synTypes = []
                for i in synonymsIndices:
                    synTypes.append('Name')

        molcount = 0
        for line in f:
            molcount += 1
            fields = line.rstrip().split('\t')
            if extIDindex is None:
                # No field with the ID of the substance in the source of origin 
                # has been provided so one will be generated.
                extID = 'mol%0.8d'%molcount
            else:
                try:
                    extID = fields[extIDindex]
                except:
                    if smi is None:
                        continue
                    else:
                        extID = 'mol%0.8d'%molcount

            # Get substance ID
            subsID = qt.getSubsID(conn, sourceID, extID)
            if not subsID:
                # The substance is not in the DB
                continue
                
            # Add synonyms
            synD = {'ExternalID': set([extID])}
            if synonymsIndices:
                for i in range(len(synonymsIndices)):
                    sindex = synonymsIndices[i]
                    stype = synTypes[i]
                    if len(fields) <= sindex: continue
                    syn = fields[sindex].strip()
                    if syn == 'N/A' or syn == '': continue
                    if stype not in synD:
                        synD[stype] = set([syn])
                    else:
                        synD.add(syn)
            addSynonyms(conn, subsID, synD)

def addAnnotationsFromFile(conn, sourceID, fname, 
                            extIDindex= None, extIDfield= None, 
                            ann= None,
                            annIndex= None, annField= None, 
                            annType = None, annTypeField= None, 
                            annTypeIndex= None,
                            header= True):
    """
        Process a text file with smiles strings of substances from a given source.
        Arguments:
          - conn: psycopg2 connection to the database.
          - sourceID: id for the source of origin from the 'source' table.
          - fname: Input file name.
          - extIDindex: Optional. Index of the column containing the id of the substance in the source of origin (default: None). If None, an id will be generated with a substance counter.
          - extIDfield: Optional. Name of the header of the column containing the substance id (default: None). If None, an id will be generated with a substance counter.
          - smilesIndex: Optional. Index of the column containing the substance's smiles string (default: 1). 
          - smilesField: Optional. Name of the header of the column containing the substance's smiles string (default: 'smiles').
          - annIndex: Optional. Index of the column containing the substance's annotation (default: 1). 
          - annField: Optional. Name of the header of the column containing the substance's annotation (default: 'Annotation').
          - annType: Optional. Type of all the annotations of the substances for this source (default: None). It should be one of the following: 'Confirmed', 'Suspected', 'Negative'.
          - linkIndex: Optional. Index of the column containing a link to the substance information page (default: None). 
          - linkField: Optional. Name of the header of the column containing a link to the substance information page (default: None).
          - synonymsIndices: Optional. List of indices of the column(s) containing synonyms of the substance (default: None). Synonym type will be 'Name'.
          - synonymsFields: Optional. List of name(s) of the header of the column(s) containing synonyms of the substance (default: None). 
          - header: Boolean indicating if the file has a header (default: False).
    """         
    with open(fname) as f:
        if header: 
            header_names = f.readline().rstrip().split('\t')
            if extIDfield:
                extIDindex = header_names.index(extIDfield)

            if annField:
                annIndex = header_names.index(annField)

            if annTypeField:
                annTypeIndex = header_names.index(annTypeField)

        molcount = 0
        for line in f:
            molcount += 1
            fields = line.rstrip().split('\t')

            if extIDindex is None:
                # No field with the ID of the substance in the source of origin 
                # has been provided so one will be generated.
                extID = 'mol%0.8d'%molcount
            else:
                try:
                    extID = fields[extIDindex]
                except:
                    extID = 'mol%0.8d'%molcount

            # Get the subsance ID
            try:
                subsID = qt.getSubsID(conn, sourceID, extID)
            except:
                continue

            # Add annotations
            annotation = None
            if not ann:
                try:
                    annotation = fields[annIndex]
                except:
                    annotation = None
                else:
                    if annotation == '****':
                        annotation = None
            else:
                annotation = ann
            
            if annotation:
                if annTypeIndex is not None:
                    try:
                        ann_type = fields[annTypeIndex]
                    except:
                        ann_type = None
                else:
                    ann_type = annType
                addAnnotation(conn, subsID, annotation, ann_type)