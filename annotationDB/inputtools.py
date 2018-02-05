import sys, os, datetime, pickle
from annotationDB import querytools as qt

import psycopg2

annClassification_file = "annD.pkl"
with open(os.path.join(os.path.dirname(__file__), "data",  annClassification_file), 'rb') as annC_fh:
    annD = pickle.load(annC_fh)

def openconnection(host='gea', dbname='annotations', user='postgres', password=''):
    """
        Opens a psycopg2 connection to the PostgreSQL DB and returns the connection object.
        Arguments:
          - host
          - dbname
          - user
          - password
    """
    conn = psycopg2.connect(host=host, dbname=dbname, user=user, password=password)

    return conn

def standardiseAnnotation(ann):
    """
    From a given toxicity annotation return a dictionary with:
        - Corrected annotation: Standardised annotation.
        - Category: Toxicity group (CMR, PBT, Physical, etc.)
        - Type: Annotation type (Confirmed, Suspected, Negative)
        - General annotation: For certain annotations, a more general classification of the toxicity. i.e. for Carc. 1A the corresponding General Annotation would be 'Carcinogenic'. This allows to group certain types of toxicity with more detail than the category (in this case CMR).
    """
    return annD[ann]

def addAnnotation(conn, subsID, ann, annType=None, annCategory=None, generalAnn=None):
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
            annType = 'Confirmed'
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

    cmd = "INSERT INTO subs_ann (subsid, annid, original_annotation, type)\
           SELECT %s, %s, %s, %s \
           WHERE NOT EXISTS (SELECT subsid, annid, original_annotation, type \
               FROM subs_ann \
               WHERE subsid= %s AND annid= %s)"
    curs.execute(cmd, (subsID, annID, ann, annType, subsID, annID))
    conn.commit()

def addSynonyms(conn, subsID, synD):
    """
        Add all synonyms provided for a given substance to the corresponding table.
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
          - version: Optional. Source's version (default: None). If None, if another source with the same name exists the version will be increased by one, else the version will be 1.
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
            version += str(int(oldVersion[0]))
            
    cmd = "SELECT id FROM source WHERE name = %s AND version = %s;"
    curs.execute(cmd, (sourceName, str(version)))
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
    
    return sourceID

def addSubstance(conn, sourceID, extID, link= None):
    """
        Insert a substance to the DB. 
        Arguments:
          - conn: psycopg2 connection to the database.
          - sourceID: id for the source of origin from the 'source' table.
          - extID: id for the subsance in the source of origin.
          - link: Optional. Link to information on the substance (default: None).

        Returns the substance id from the 'substance' table.
    """
    curs = conn.cursor()
            
    cmd = "SELECT id FROM substance WHERE sourceid = %s \
            AND externalid = %s;"
    curs.execute(cmd, (sourceID, extID))
    subsID = curs.fetchone()
    conn.commit()

    if not subsID:
        cmd = "INSERT INTO substance (sourceid, externalid, link) \
                VALUES (%s, %s, %s)"
        curs.execute(cmd, (sourceID, extID, link))
        conn.commit()
        cmd = "SELECT currval('substance_id_seq');"
        curs.execute(cmd)
        subsID = curs.fetchone()[0]
        conn.commit()
        addSynonyms(conn, subsID, {'ExternalID': extID})
    else:
        subsID = subsID[0]
            
    return subsID

def addSubstancesFromFile(conn, sourceID, fname, extIDindex= None, extIDfield= None, 
                        annIndex= 1, annField= 'smiles', annType = None, typeIndex= None, typeField= None, linkIndex= None, linkField= None, synonymsIndices= None, synonymsFields= None, header= True):
    """
        Process a text file with substances and annoations from a given source.
        Arguments:
          - conn: psycopg2 connection to the database.
          - sourceID: id for the source of origin from the 'source' table.
          - fname: Input file name.
          - extIDindex: Optional. Index of the column containing the id of the substance in the source of origin (default: None). If None, an id will be generated with a substance counter.
          - extIDfield: Optional. Name of the header of the column containing the substance id (default: None). If None, an id will be generated with a substance counter.
          - annIndex: Optional. Index of the column containing the substance's annotation (default: 1). 
          - annField: Optional. Name of the header of the column containing the substance's annotation (default: 'Annotation').
          - annType: Optional. Type of all the annotations of the substances for this source (default: None). It should be one of the following: 'Confirmed', 'Suspected', 'Negative'.
          - typeIndex: Optional. Index of the column containing the type of the annotation in the source of origin (default: None). It should be one of the following: 'Confirmed', 'Suspected', 'Negative'.
          - typeField: Optional. Name of the header of the column containing the type of the annotation (default: None). It should be one of the following: 'Confirmed', 'Suspected', 'Negative'.
          - linkIndex: Optional. Index of the column containing a link to the substance information page (default: None). 
          - linkField: Optional. Name of the header of the column containing a link to the substance information page (default: None).
          - synonymsIndices: Optional. List of indices of the column(s) containing synonyms of the substance (default: None). Synonym type will be 'Name'.
          - synonymsFields: Optional. List of name(s) of the header of the column(s) containing synonyms of the substance (default: None). 
          - header: Boolean indicating if the file has a header (default: False).
    """
    curs = conn.cursor()          
    with open(fname) as f:
        if header: 
            header = f.readline().rstrip().split('\t')
            if extIDfield:
                extIDindex = header.index(extIDfield)

            if annField:
                annIndex = header.index(annField)

            if typeField:
                typeIndex = header.index(typeField)

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
                ann = fields[annIndex]
            except:
                continue
            else:
                if ann == '****':
                    continue

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
            subsID = addSubstance(conn, sourceID, extID= extID, link= link)

            # Add the annotation
            if typeIndex is not None:
                try:
                    ann_type = fields[typeIndex]
                except:
                    ann_type = None
            else:
                ann_type = annType
            addAnnotation(conn, subsID, ann, ann_type)
                
            # Add synonyms
            synD = {}
            if synonymsIndices:
                for i in range(len(synonymsIndices)):
                    sindex = synonymsIndices[i]
                    stype = synTypes[i]
                    if len(fields) <= sindex: continue
                    syn = fields[sindex].strip()
                    if syn in ('N/A', '-', '', '--', '---'): continue
                    if stype not in synD:
                        synD[stype] = set([syn])
                    else:
                        synD.add(syn)
            addSynonyms(conn, subsID, synD)

def addSubstanceFromQuery(conn, sourceID, cmd, host='gea', dbname='compounds', user='postgres',
                          password='', extIDf= 'ID', linkF= None):
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
          - linkF: Optional. Name of the DB field containing the link to the substance information page (default: None). 
    """
    
    # Open connection to the query DB
    qconn = psycopg2.connect(host=host, dbname=dbname, user=user, password=password)
    df = pd.read_sql(cmd, con=qconn)
    for index, row in df.iterrows():
        extID = row[extIDf]
        if not linkF: link = None
        else: link= row[linkF]
        try:
            mol = Chem.MolFromSmiles(smi)
        except:
            addEmptySubstance(conn, sourceID, extID, link)
            continue
        addSubstance(conn, sourceID, extID= extID, mol= mol, link= link)

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
    curs = conn.cursor()
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
