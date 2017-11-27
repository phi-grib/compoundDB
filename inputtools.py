import sys, datetime
from moleculeHelper import *
from process_smiles import *

import psycopg2
from psycopg2 import extras
from rdkit.Chem import *
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.Crippen import MolLogP

def openconnection(host='gea', dbname='compounds', user='postgres', password='DBAdmin'):
    conn = psycopg2.connect(host=host, dbname=dbname, user=user, password=password)
    curs = conn.cursor()
    curs.execute('create extension if not exists rdkit;')

    return conn

def addSource(conn, name, version= None, description= None, link= None, \
              date= datetime.datetime.now().date()):
    curs = conn.cursor()
    # Check if version is provided, otherwise generate it
    if version is None:
        cmd = "SELECT id FROM source WHERE name = %s ORDER BY version DESC;"
        curs.execute(cmd, (sourcename,))
        oldVersion = curs.fetchone()[0]
        conn.commit()
        if oldVersion is None:
            version = 1
        else:
            version += int(oldVersion)
            
    cmd = "SELECT id FROM source WHERE name = %s AND version = %s;"
    curs.execute(cmd, (name, version))
    sourceID = curs.fetchone()[0]
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
    
    return sourceID

def addEmptyCompound(conn, subsID, smiles):
    curs = conn.cursor()
    cmd = "SELECT id FROM compound WHERE smiles = %s;"
    curs.execute(cmd, (smiles,))
    cmpdID = curs.fetchone()[0]
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
            
    cmd = "INSERT INTO subs_cmpd (subsid, cmpdid)\
       VALUES (%s, %s)"
    curs.execute(cmd, (subsID, cmpdID))
    conn.commit()
    
    return (cmpdID, None)

def addCompound(conn, sourceID, subsID, smiles= None, mol= None, ismetal=False):
    curs = conn.cursor()
    if smiles is None and mol is None:
        (cmpdID, mol) = addEmptyCompound(conn, subsID, smiles)
    elif mol is None:
        try:
            mol = MolFromSmiles(smiles)
        except:
            print (smiles)
            (cmpdID, mol) = addEmptyCompound(conn, subsID, smiles)
            
    if mol is None:
        # Molecule object has been created, but only as an empty instnce
        (cmpdID, mol) = addEmptyCompound(conn, subsID, smiles)            
    else:
        inchi = MolToInchi(mol)
        inchikey = InchiToInchiKey(inchi)

        cmd = "SELECT id FROM compound WHERE inchikey = %s;"
        curs.execute(cmd, (inchikey,))
        cmpdID = curs.fetchone()[0]
        conn.commit()

        if cmpdID is None:
            if smiles is None:
                smiles = MolToSmiles(mol)
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
            
    cmd = "INSERT INTO subs_cmpd (subsid, cmpdid)\
       VALUES (%s, %s)"
    curs.execute(cmd, (subsID, cmpdID))
    conn.commit()
    
    return (cmpdID, mol)

def addEmptySubstance(conn, sourceID, extID, molID= None, link= None):
    curs = conn.cursor()
    cmd = "SELECT id FROM substance WHERE sourceid = %s \
            AND externalid = %s;"
    curs.execute(cmd, (sourceID, extID))
    subsID = curs.fetchone()[0]
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
    
    return (subsID, None)

def addSubstance(conn, sourceID, extID, smiles= None, mol= None, link= None):
    curs = conn.cursor()
    if smiles is None and mol is None:
        (subsID, mol) = addEmptySubstance(curs, conn, sourceID, extID, link)
    elif mol is None:
        try:
            mol = MolFromSmiles(smiles)
        except:
            (subsID, mol) = addEmptySubstance(curs, conn, sourceID, extID, link)
            
    if mol is None:
        # Molecule object has been created, but only as an empty instnce
        (subsID, mol) = addEmptySubstance(curs, conn, sourceID, extID, link)
    else:
        cmd = "SELECT id FROM substance WHERE sourceid = %s \
                AND externalid = %s;"
        curs.execute(cmd, (sourceID, extID))
        subsID = curs.fetchone()[0]
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
        
            stdD = std(mol)
            for smiles in stdD:
                (cmpd, ismetal) = stdD[smiles]
                (cmpdID, cmpdmol) = addCompound(conn, sourceID, subsID, smiles= smiles, mol= cmpd, ismetal=ismetal)
            
    return (subsID, mol)

def addSubstanceFile(conn, sourceID, fname, ftype, molID= None, smilesI= 1, linkI= None, header= False, synonyms= None):
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
            molcount = 0
            for line in f:
                molcount += 1
                fields = line.rstrip().split('\t')
                smi = fields[smilesI]
                if molID is None:
                    extID = 'mol%0.8d'%molcount
                else:
                    extID = fields[molID]
                if linkI is None: link = None
                else: link= fields[linkI]

                # Add the subsance
                try:
                    mol = MolFromSmiles(smi)
                except:
                    (subsID, mol) = addEmptySubstance(conn, dbID, extID, link)
                else:
                    (subsID, mol) = addSubstance(conn, sourceID, extID= extID, smiles= smiles, \
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
                extID = getNameFromEmpty(suppl, molcount-1, molID)                    
                (subsID, mol) = addEmptySubstance(conn, dbID, extID, link)
                
            else:                
                extID = getName(mol, molcount, molID)                    
                link = None
                if linkI is not None:
                    try:
                        link = GetProp(mol, linkI)
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
                          password='DBAdmin', extIDf= 'ID', smilesF= 'smiles', \
                          linkF= None):
    
    # Open connection to the query DB
    qconn = psycopg2.connect(host=host, dbname=dbname, user=user, password=password)
    qcurs = qconn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    qcurs.execute(cmd)
    row = qcurs.fetchone()[0]
    while row:
        extID = row[extIDf]
        smi = row[smilesF]
        if linkF is None: link = None
        else: link= row[linkF]
        try:
            mol = MolFromSmiles(smi)
        except:
            addEmptySubstance(conn, dbID, extID, link)
            continue
        addSubstance(conn, sourceID, extID= extID, smiles= smi, \
                     mol= mol, link= link)
        row = qcurs.fetchone()[0]
    qcurs.commit()

