%read in an excel spreadsheet 
filename = '/home/flagg/Downloads/Rrs_correction.xls'
sheet = 'Rrs Optimization'
xlsRange = 'D11:D115'
wavelengths = xlsread(filename,sheet,xlsRange);