# Input output - version 02/11/2017

#IMPORT PACKAGES
import os
import csv
import numpy as np


#IMPORT CSV FILE
def loadcsv(file_name, dict1=True, dict2=True):
    """
    FUNCTION to read a .CSV file into a data structure of choice
    INPUT
    - <file_name>: csv file name, including extension. CSV file is assumed to have a header row
    OPTIONAL INPUT
    - <dict1> determines whether data rows are treated as individual dictionaries (true) or lists (false). 
              If a dictionary is used, header rows function and dictionary keys
    - <dict2> determined wheters columns are treated as dictionaries (true) or lists (false)
              if dict1 = True
                 if dict2 = True  --> dictionary with dictionaries. First data column is used as main dictionary keys
                 if dict2 = False --> single dictionary with column data in dictionary fields
              if dict1 = False
                 if dict2 = True  --> list with a single dictionary per row
                 if dict2 = False --> list of row data lists
    OUTPUT
    - Data structure
    """
    
    #open file
    reader = csv.reader(open(file_name))
    
    #Declare some memory space
    if dict1 == False:
        #load data to list format
        dat = []
    else:
        #load data to dictionary format
        dat = {}
        
    #Loop through rows in input data file
    for row in reader:
        if reader.line_num == 1:
            #first row = header row
            headers = row

            #define space in dictionary with lists, so that data can be appended
            if dict1 == True and dict2 == False:
                for i in headers:
                    dat[i] = []            

        else:
            #not first row
            #convert data to doubles if possible
            for el in range(len(row)):
                try:
                    row[el] = float(row[el])
                except:
                    row[el] = row[el]
        
            if dict1 == False:
                #load data to list format
                if dict2 == True:
                    #Load data to single list with dictionaries
                    dat.append(dict(zip(headers, row)))
                else:
                    #Load data to list with lists
                    dat.append(row)
    
            else:
            #load data to dictionary format, using first column as main dictionary key            
                if dict2 == True:
                    #Load data to dictionary with dictionaries, using first column nas main dictionary key
                    dat[int(row[0])] = dict(zip(headers[1:], row[1:]))
                else:
                    #Load data to single dictionary, containing lists of data
                    for i,j in zip(headers,row):
                        dat[i].append(j)

    #return data
    return(dat)


#LOAD NODAL AND SEGMENT DATA
def loadnodesegment(file_name_node, file_name_segment):
    """
    FUNCTION to load nodal and segment data
    INPUT
    - <file_name_node>: csv file name of nodal data file, including extension. CSV file is assumed to have a header row
    - <file_name_segment>: csv file name of segment data file, including extension. CSV file is assumed to have a header row
    OUTPUT
    - <dnode> dictionary with nodal data
    - <dsegm> dictionary with segment data
    """    
    
    #load nodal data
    dnode = loadcsv(file_name_node, dict1=True, dict2=True)
    for ikey,i in dnode.items():
        #convert boundary condition input to boolean
        i['bound_X']     = bool(i['bound_X'])
        i['bound_Y']     = bool(i['bound_Y'])
        i['bound_Theta'] = bool(i['bound_Theta'])

    #load segment data
    dsegm = loadcsv(file_name_segment, dict1=True, dict2=True)
    for ikey,i in dsegm.items():
        #convert nodeID's to integer
        i['NodeID1']   = int(i['NodeID1'])
        i['NodeID2']   = int(i['NodeID2'])
        #find begin and end coordinates for segment
        i['X1']        = dnode[i['NodeID1']]['X']
        i['X2']        = dnode[i['NodeID2']]['X']
        i['Y1']        = dnode[i['NodeID1']]['Y']
        i['Y2']        = dnode[i['NodeID2']]['Y']
        #additional segment properties (length, angle)
        i['L']     = np.sqrt((i['X2']-i['X1'])**2 + (i['Y2']-i['Y1'])**2)
        if i['X2']-i['X1'] == 0:
            if i['Y2'] >= i['Y1']:
                i['Theta'] = np.pi/2.0
            else:
                i['Theta'] = 3.0/2.0 * np.pi
        elif i['X2']-i['X1'] > 0:
            i['Theta'] = np.arctan((i['Y2']-i['Y1']) / (i['X2']-i['X1']))
        else:
            i['Theta'] = np.arctan((i['Y2']-i['Y1']) / (i['X2']-i['X1'])) + np.pi
        #additional segment properties (area, second moment of intertia)
        i['A'] = np.pi / 4.0  * i['d']**2
        i['I'] = np.pi / 64.0 * i['d']**4

    #find all segments connected to node + side to which they connect (side=0, first node, side=1, second node
    for ikey,i in dnode.items():
        i['SegmentID']   = [j for j,k in zip(dsegm.keys(), [s['NodeID1'] for skey,s in dsegm.items()]) if ikey==k]  +  [j for j,k in zip(dsegm.keys(), [s['NodeID2'] for skey,s in dsegm.items()]) if ikey==k]
        i['SegmentSide'] = [0 for j,k in zip(dsegm.keys(), [s['NodeID1'] for skey,s in dsegm.items()]) if ikey==k]  +  [1 for j,k in zip(dsegm.keys(), [s['NodeID2'] for skey,s in dsegm.items()]) if ikey==k]

    #Determine node type ('end' or 'middle' node) (end node=only 1 segment connected, 'middle'=more than 1 segment connected)
    for ikey,i in dnode.items():
        i['NodeType']   = 'end' if len(i['SegmentID']) == 1 else 'middle'

    #return
    return(dnode, dsegm)


#Function to transform data in a dictionary to numpy array
def dict2matrix(dic, field_names=None, rows=False):
    """
    FUNCTION to transform dictionary data to numpy array
    INPUT
    - <dic>:          dictionary
    OPTIONAL INPUT
    - <field_names>:  list of names of fields taken into account (can also be single string with one field)
    - <rows>:         Boolean indicating whether dictionary data outputted as rows in matrix (default = columns)
    OUTPUT
    - <out>:          numpy array with data
    """
    #get list of field names
    if field_names is None:
        #no input field_names --> take all fields
        field_names = list(dic.keys())
    elif isinstance(field_names, str):
        #single field string --> take only this field
        field_names = [field_names]

    #get data
    if rows is True:
        #dictionary items are rows in output array
        return(np.array([i for ikey,i in dic.items()]))
    else:
        #dictionary items are columns in output array
        return(np.transpose(np.array([i for ikey,i in dic.items()])))


def dictdict2dict(dic, field_names1=None, field_names2=None, key1_header='key1'):
    """
    FUNCTION to create single dictionary out of dictionary of dictionaries
    INPUT
    - <dic>: dictionary to be written to file
    OPTIONAL INPUT
    - <field_names1>: write data for certain keys in outer dictionary only (list of key names)
    - <field_names2>: write data for certain keys in inner dictionaries only (list of key names)
    - <key1_header>: column name for key of outer dictionary
    OUTPUT
    - none, but writes a file
    """
    #get list of field names1
    if field_names1 is None:
        #no input field_names1 --> take all fields
        field_names1 = list(dic.keys())
    elif isinstance(field_names1, str):
        #single field string --> take only this field
        field_names1 = [field_names1]

    #get list of field names2
    if field_names2 is None:
        #no input field_names --> take all fields
        field_names2 = list(dic[field_names1[0]].keys())
    elif isinstance(field_names2, str):
        #single field string --> take only this field
        field_names2 = [field_names2]
        
    #output
    out = dict()
    #loop through fieldnames2
    if isinstance(dic[field_names1[0]][field_names2[0]], np.ndarray) or isinstance(dic[field_names1[0]][field_names2[0]], list):
        for i in field_names2:    
            out[i] = np.concatenate([dic[j][i] for j in field_names1])
        out[key1_header] = [j for j in field_names1 for i in range(0,len(dic[j][field_names2[0]]))]
    else:
        for i in field_names2:    
            out[i] = [dic[j][i] for j in field_names1]
        out[key1_header] = field_names1
    #add key1
    #out[key1_header] = np.concatenate([np.repeat(j, len(dic[j][field_names2[0]])) for j in field_names1])
    
            
    #return
    return(out)


#EXPORT GENERIC DICTIONARY OF LISTS TO CSV FILE
def dict2csv(file_name, dic, directory=None, field_names=None, rows=False):
    """
    FUNCTION to write data in dictionary to a .CSV file. Dictionary keys function as csv headers
    INPUT
    - <file_name>: filename of the output file, including extension
    - <dic>: dictionary to be written to file
    OPTIONAL INPUT
    - <directory>: directory where the file should be saved
    - <field_names>: write data for certain keys only (list of key names)
    OUTPUT
    - none, but writes a file
    """
    #get list of field names1
    if field_names is None:
        #no input field_names1 --> take all fields
        field_names = list(dic.keys())
    elif isinstance(field_names, str):
        #single field string --> take only this field
        field_names = [field_names]

    #create folder if not alreay existing
    if not os.path.exists(directory):
        os.makedirs(directory)
    if directory is not None:
        file_name = directory + '/' + file_name
        
    #open file
    with open(file_name, 'w', encoding='utf8', newline='') as file:
        #create writer object
        writer = csv.writer(file)
        #write headers
        writer.writerow(dic.keys()) 
        #loop through length of lists in dictionary (rows in output file)
        for j in dict2matrix(dic, rows=rows):
            writer.writerow(list(j))