import sys, datetime
from phitools import moleculeHelper as mh
from standardiser import process_smiles as ps
from compoundDB import querytools as qt

import psycopg2
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.Crippen import MolLogP

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
            version += str(int(oldVersion[0]))
            
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
            
    if mol is None:
        # Molecule object has been created, but only as an empty instnce
        (subsID, mol) = addEmptySubstance(conn, sourceID, extID, link)
    else:
        cmd = "SELECT id FROM substance WHERE sourceid = %s \
                AND externalid = %s;"
        curs.execute(cmd, (sourceID, extID))
        subsID = curs.fetchone()
        conn.commit()

        if not subsID:
            inchi = Chem.MolToInchi(mol)
            inchikey = Chem.InchiToInchiKey(inchi)
            if smiles is None:
                smiles = Chem.MolToSmiles(mol)

            cmd = "INSERT INTO substance (sourceid, externalid, smiles, mol, inchi, inchikey)\
            VALUES (%s, %s, %s, mol_from_smiles(%s), %s, %s)"
            curs.execute(cmd, (sourceID, extID, smiles, smiles, inchi, inchikey))
            conn.commit()
            cmd = "SELECT currval('substance_id_seq');"
            curs.execute(cmd)
            subsID = curs.fetchone()[0]
            conn.commit()
        
            stdD = ps.std(mol)
            for smiles in stdD:
                (cmpd, ismetal) = stdD[smiles]
                cmpdID = addCompound(conn, subsID, smiles= smiles, mol= cmpd, ismetal=ismetal)
        else:
            subsID = subsID[0]
            
    return (subsID, mol)

def addSubstanceFromSmilesFile(conn, sourceID, fname, extIDindex= None, extIDfield= None, smilesIndex= 1, smilesField= 'smiles', linkIndex= None, linkField= None, synonymsIndices= None, synonymsFields= None, header= False):
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

            if smilesField:
                smilesIndex = header.index(smilesField)

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
                    if smi is None:
                        continue
                    else:
                        extID = 'mol%0.8d'%molcount

            if not linkIndex: link = None
            else: link= fields[linkIndex]

            # Add the subsance
            try:
                mol = Chem.MolFromSmiles(smi)
            except:
                (subsID, mol) = addEmptySubstance(conn, sourceID, extID, link)
            else:
                (subsID, mol) = addSubstance(conn, sourceID, extID= extID, smiles= smi, \
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

def addSubstanceFromCASFile(conn, sourceID, fname, extIDindex= None, extIDfield= None, CASindex= None, CASfield= None, linkIndex= None, linkField= None, synonymsIndices= None, synonymsFields= None, header= False):
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

            if CASfield:
                CASindex = header.index(CASfield)

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
                continue

            if extIDindex is None:
                # No field with the ID of the substance in the source of origin 
                # has been provided so one will be generated.
                extID = CAS
            else:
                try:
                    extID = fields[extIDindex]
                except:
                    extID = CAS

            if not linkIndex: link = None
            else: link= fields[linkIndex]

            # Get smiles from CAS
            # First check if it's already in the DB
            smi = qt.getStructureFromSyn(conn, syn= CAS)
            if not smi:
                # Otherwise, try to resolve it through web services
                smi = mh.resolveCAS(cas= CAS)

            # Add the subsance
            try:
                mol = Chem.MolFromSmiles(smi)
            except:
                (subsID, mol) = addEmptySubstance(conn, sourceID, extID, link)
            else:
                (subsID, mol) = addSubstance(conn, sourceID, extID= extID, smiles= smi, \
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

def addSubstanceSDFile(conn, sourceID, fname, extIDfield= None, linkField= None, synonymsFields= None):
    """
        Process an SD file of substances from a given source.
        Arguments:
          - conn: psycopg2 connection to the database.
          - sourceID: id for the source of origin from the 'source' table.
          - fname: Input file name.
          - extIDfield: Optional. Name of the field containing the substance id (default: None). If None, an id will be generated with a substance counter.
          - synonymsFields: Optional. List of name(s) of the field(s) containing synonyms of the substance (default: None). 
    """
    curs = conn.cursor()
    suppl = SDMolSupplier(fname)
    molcount = 0
    for mol in suppl:
        molcount += 1
        if mol is None:
            extID = mh.getNameFromEmpty(suppl, molcount-1, extID)
            (subsID, mol) = addEmptySubstance(conn, sourceID, extID, link)
            
        else:                
            extID = mh.getName(mol, molcount, extID)                    
            link = None
            if linkI:
                try:
                    link = mol.GetProp(linkI)
                except:
                    pass
            (subsID, mol) = addSubstance(conn, sourceID, extID= extID, mol= mol, link= link)
        
        # Add synonyms
        synD = {'ExternalID': set([extID])}
        if synonyms:
            for synType in synonyms:
                try:
                    syn = mol.GetProp(synType)
                    if syn != 'N/A':
                        if synType not in synD:
                            synD[synType] = set([syn])
                        synD[synType].add(syn)
                except:
                    pass
        addSynonyms(conn, subsID, synD)

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
