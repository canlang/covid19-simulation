% clc;clearvars;
%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: /Users/sjlee/Developer/COVID19/SIR_sim/KCDC_data.xlsx
%    Worksheet: SummaryByDay
%
% Auto-generated by MATLAB on 12-Oct-2020 14:30:45

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 16);

% Specify sheet and range
opts.Sheet = "SummaryByDay";
opts.DataRange = "A1:P230";

% Specify column names and types
opts.VariableNames = ["date", "PressRelease_Time", "PressRelease_DateTime", "ReportDate", "total", "confirmed_total", "confirmed_discharged", "confirmed_isolated", "confirmed_deceased", "piu_total", "piu_beingtested", "piu_testednegative", "VarName13", "VarName14", "VarName15", "VarName16"];
opts.VariableTypes = ["datetime", "categorical", "datetime", "string", "double", "double", "double", "double", "double", "string", "double", "double", "string", "string", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, ["ReportDate", "piu_total", "VarName13", "VarName14"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["PressRelease_Time", "ReportDate", "piu_total", "VarName13", "VarName14"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "date", "InputFormat", "");
opts = setvaropts(opts, "PressRelease_DateTime", "InputFormat", "");

% Import the data
if ~exist('T_KCDC','var')
    T_KCDC = readtable("KCDC_data.xlsx", opts, "UseExcel", false);
end


%% Clear temporary variables
clear opts
%%
% T = table(Date,InfectedCum,Isolated);
T = table(T_KCDC.date,T_KCDC.confirmed_total,T_KCDC.confirmed_discharged,T_KCDC.confirmed_isolated,...
    'VariableNames',{'Date' 'InfectedCum' 'Discharged' 'Isolated'});
