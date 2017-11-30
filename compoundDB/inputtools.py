import sys, datetime
from phitools import moleculeHelper as mh
from standardiser import process_smiles as ps

import psycopg2
from psycopg2 import extras
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
    cmd = "INSERT INTO synonym (subsid, type, name) \
            VALUES (%s, %s, %s)"
    curs.executemany(cmd, [(subsID,syntype,syn) for syntype in synD for syn in synD[syntype]])
    conn.commit()

def addSource(conn, name, version= None, description= None, link= None, \
              date= datetime.datetime.now().date()):
    """
        Insert a source to the DB.
        Arguments:
          - conn: psycopg2 connection to the database.
          - name: Source name.
          - version: Optional. Source's version (default: None). If None, if another source with the same name exists the version will be increased by one, else the version will be 1.
          - description: Optional. Verbose description of the source (default: None).
          - link: Optional. Link to the source's home page (default: None).
          - date: Optional. Date of source's creation (default: Today).

        Returns source id from the 'source' table.
    """
    curs = conn.cursor()
    # Check if version is provided, otherwise generate it
    if version is None:
        cmd = "SELECT id FROM source WHERE name = %s ORDER BY version DESC;"
        curs.execute(cmd, (sourcename,))
        oldVersion = curs.fetchone()
        conn.commit()
        if oldVersion is None:
            version = 1
        else:
            version += int(oldVersion[0])
            
    cmd = "SELECT id FROM source WHERE name = %s AND version = %s;"
    curs.execute(cmd, (name, version))
    sourceID = curs.fetchone()
    conn.commit()
    
    if sourceID is None:
        cmd = "INSERT INTO source (name, version, description, link, added)\
               VALUES (%s, %s, %s, %s, %s)"
        curs.execute(cmd, (name, version, description, link, date))
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
            inchi = MolToInchi(mol)
            inchikey = InchiToInchiKey(inchi)
            if smiles is None:
                smiles = MolToSmiles(mol)

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

def addSubstanceFile(conn, sourceID, fname, ftype, extID= None, smilesI= 1, linkI= None, header= False, synonyms= None):
    """
        Process a file with substances from a given source. The file can either be in SD format or a text file with smiles strings.
        Arguments:
          - conn: psycopg2 connection to the database.
          - sourceID: id for the source of origin from the 'source' table.
          - fname: Input file name.
          - ftype: Input file type (sdf / smi).
          - extID: Field containing the id of the substance in the source of origin (default: None). 
            - If the input format is 'sdf', it can be a field name or None. If None, the substance name from the first line will be used and if this is empty an id will be generated with a substance counter. 
            - If the input format is 'smi', it can either be an integer corresponding to the column index or the name of the header of the column containing the substance id. If None, an id will be generated with a substance counter.
          - smilesI: Optional. Field containing substance's smiles string. It can either be an integer corresponding to the column index or the name of the header of the column containing the smiles string (default: 1).
          - linkI: Optional. Field containing a link to the substance information page (default: None). 
            - If the input format is 'sdf', it can be a field name or None.  
            - If the input format is 'smi', it can either be an integer corresponding to the column index or the name of the header of the column containing the link. 
    """
    curs = conn.cursor()
    if ftype == 'smi':                
        with open(fname) as f:
            if header: 
                header = f.readline().rstrip().split('\t')
                if synonyms is not None:
                    if all(isinstance(n, int) for n in synonyms):
                        synIndices = list(synonyms)
                        synTypes = []
                        for i in synIndices:
                            synTypes.append(header[i])
                    else:
                        synTypes = list(synonyms)
                        for t in synTypes:
                            synIndices.append(header.index(t))
                if not isinstance(extID, int):
                    extID = header.index(extID)
                if not isinstance(smilesI, int):
                    smilesI = header.index(smilesI)
            molcount = 0
            for line in f:
                molcount += 1
                fields = line.rstrip().split('\t')
                smi = fields[smilesI]
                if extID is None:
                    extID = 'mol%0.8d'%molcount
                else:
                    extID = fields[extID]
                if linkI is None: link = None
                else: link= fields[linkI]

                # Add the subsance
                try:
                    mol = Chem.MolFromSmiles(smi)
                except:
                    (subsID, mol) = addEmptySubstance(conn, dbID, extID, link)
                else:
                    (subsID, mol) = addSubstance(conn, sourceID, extID= extID, smiles= smi, \
                                                 mol= mol, link= link)
                    
                # Add synonyms
                synD = {'ExternalID': set([extID])}
                if synonyms is None:
                    addSynonyms(conn, subsID, synD)
                else:
                    for i in range(synonyms):
                        syn = fields[synIndices[i]]
                        if syn != 'N/A': synD[synTypes[i]] = set([])
                    addSynonyms(conn, subsID, synD)
                    
    elif ftype == 'sdf':
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
                if linkI is not None:
                    try:
                        link = mol.GetProp(linkI)
                    except:
                        pass
                (subsID, mol) = addSubstance(conn, sourceID, extID= extID, mol= mol, link= link)
            
            # Add synonyms
            synD = {'ExternalID': set([extID])}
            if synonyms is not None:
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
    qcurs = qconn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    qcurs.execute(cmd)
    row = qcurs.fetchone()
    while row:
        extID = row[extIDf]
        smi = row[smilesF]
        if linkF is None: link = None
        else: link= row[linkF]
        try:
            mol = Chem.MolFromSmiles(smi)
        except:
            addEmptySubstance(conn, dbID, extID, link)
            continue
        addSubstance(conn, sourceID, extID= extID, smiles= smi, \
                     mol= mol, link= link)
        row = qcurs.fetchone()
    qcurs.commit()

