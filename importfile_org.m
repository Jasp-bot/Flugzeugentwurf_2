function aircraftDataorg = importfile_org(workbookFile, sheetName, dataLines)
%IMPORTFILE Import data from a spreadsheet
%  AIRCRAFTDATAORG = IMPORTFILE(FILE) reads data from the first
%  worksheet in the Microsoft Excel spreadsheet file named FILE.
%  Returns the numeric data.
%
%  AIRCRAFTDATAORG = IMPORTFILE(FILE, SHEET) reads from the specified
%  worksheet.
%
%  AIRCRAFTDATAORG = IMPORTFILE(FILE, SHEET, DATALINES) reads from the
%  specified worksheet for the specified row interval(s). Specify
%  DATALINES as a positive scalar integer or a N-by-2 array of positive
%  scalar integers for dis-contiguous row intervals.
%
%  Example:
%  aircraftDataorg = importfile("C:\Users\jaspe\OneDrive\UNI\Flugzeugentwurf\Aufgabe\aircraftData_org.xlsx", "aircraftData", [2, 183]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 24-Nov-2022 16:58:56

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 2
    dataLines = [2, 183];
end

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 14);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = "A" + dataLines(1, 1) + ":N" + dataLines(1, 2);

% Specify column names and types
opts.VariableNames = ["typeOfAircraft", "entryIntoServicea", "PAXdesignPoint", "wingLoadingNm", "wingAream", "maximumTakeOffMasskg", "maximumLandingMasskg", "maximumZeroFuelMasskg", "operationEmptyMasskg", "cruiseSpeedkmh", "physRangeByMaximumPayloadkm", "maximumPayloadMasskg", "physRangeByMaximumFuelkm", "payloadMassByMaximumFuelkg"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
aircraftDataorg = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = "A" + dataLines(idx, 1) + ":N" + dataLines(idx, 2);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    aircraftDataorg = [aircraftDataorg; tb]; %#ok<AGROW>
end

%% Convert to output type
aircraftDataorg = table2array(aircraftDataorg);
end