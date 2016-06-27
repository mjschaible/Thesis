# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 10:50:41 2016

Parses and converts lammps data dump files to calculate aggregate values for
user specified data columns. 

@author: ksiller
"""

import pandas as pd
import numpy as np
import sys
import string


def createDataframeFromDump(dumpfile):
    # create dummy headers for 14 columns. 
    # Chis is needed because rows are of unequal length.
    # Missing values will be assigned as NaN
    my_cols=string.ascii_uppercase[0:14]
    dataframe = pd.read_table(dumpfile, delimiter=' ', names=my_cols) #, engine='python')
    # Pull out all row indices that contain 'TIMESTEP' in 2nd column ('B')
    timestep_startrow=dataframe[dataframe.B == 'TIMESTEP'].index
    # Get list of TIMESTEP values
    timesteps = dataframe['A'].loc[timestep_startrow + 1]
    
    # Create Series for TIMESTEPS
    rowsPerTimestep=timestep_startrow[1]-timestep_startrow[0]
    matrix = [ ([i] * rowsPerTimestep) for i in timesteps ]
    timeSeriesData=np.hstack(np.asarray(matrix)) # flatten to 1D array
    # Add TIMESTEP series as new column 
    dataframe['timestep']=pd.Series(timeSeriesData,index=dataframe.index)
    
    # Rename columns
    dataframe.columns=['type','id', 'mol', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'c_pe_atom', 'c_ke_atom', 'L','M','N','timestep']
    
    # Drop all rows that are headers. Headers contain NaN in column 'vx' 
    dataframe=dataframe.dropna(subset=['vx'])
    dataframe=dataframe[dataframe.type != 'ITEM:']
    # Remove additional columns
    dataframe=dataframe.dropna(axis=1)
    # Convert entire dataframe to numeric data type
    dataframe=dataframe.apply(pd.to_numeric)
    return dataframe
    
    
def createDataframeFromConvDump(convDumpfile):    
    dataframe = pd.read_table(convDumpfile, delimiter=' ')
    return dataframe            

                
def main():                
    if len(sys.argv) < 2:
        print ("File name required")
        exit()
    
    # create DataFrame and save as csv --> converted dump file
    dataframe=createDataframeFromDump(sys.argv[1])
    dataframe.to_csv(sys.argv[1]+'_conv',sep=' ', index=False)
    # comment lines 59 & 60 and use the function call line 63 instead if reading 
    # from converted dump file 
    #dataframe=createDataframeFromConvDump(sys.argv[1])
    
    # aggregate some numbers for data grouped by TIMESTEPS and mol
    groups=dataframe.groupby(['mol', 'timestep'])
    sums=groups['c_pe_atom','c_ke_atom'].sum()
    print(sums)

    # example: pull out c_ke_atom sum values for mol=1 for all timesteps 
    mol_groups=sums.groupby(level='mol')
    mol_no=1
    print(mol_groups.get_group(mol_no)['c_ke_atom'])

if __name__ == "__main__":
    main()
 

