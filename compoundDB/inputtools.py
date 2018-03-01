import sys, datetime, os, pickle
import psycopg2

from phitools import moleculeHelper as mh
from standardiser import process_smiles as ps
from compoundDB import querytools as qt

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
    if smiles is None and mol is None:
        (subsID, mol) = addEmptySubstance(conn, sourceID, extID, link)
    elif mol is None:
        try:
            mol = Chem.MolFromSmiles(smiles)
        except:
            (subsID, mol) = addEmptySubstance(conn, sourceID, extID, link)
            
    cmd = "SELECT id FROM substance WHERE sourceid = %s \
            AND externalid = %s;"
    curs.execute(cmd, (sourceID, extID))
    subsID = curs.fetchone()
    conn.commit()

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
        
            stdD = ps.std(mol)
            for smiles in stdD:
                (cmpd, ismetal, passed, errmessage) = stdD[smiles]
                cmpdID = addCompound(conn, subsID, smiles= smiles, mol= cmpd, ismetal=ismetal)
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

            if not ann:
                try:
                    ann = fields[annIndex]
                except:
                    ann = None
                else:
                    if ann == '****':
                        ann = None

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
            if ann:
                if annTypeIndex is not None:
                    try:
                        ann_type = fields[annTypeIndex]
                    except:
                        ann_type = None
                else:
                    ann_type = annType
                addAnnotation(conn, subsID, ann, ann_type)

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

            if not ann:
                try:
                    ann = fields[annIndex]
                except:
                    ann = None
                else:
                    if ann == '****':
                        ann = None

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
            if ann:
                if annTypeIndex is not None:
                    try:
                        ann_type = fields[annTypeIndex]
                    except:
                        ann_type = None
                else:
                    ann_type = annType
                addAnnotation(conn, subsID, ann, ann_type)

def addSubstanceSDFile(conn, sourceID, fname, extIDfield= None, linkField= None, 
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
        if annField:
            try:
                ann = mol.GetProp(annField)
            except:
                continue
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
