% This script loops through the analyzed_spFR_results and groups
% them by experimental conditions: 
% For example: (these indice should match 'Cat' in the cell_id excel file) 
    % Ctrl- 1
    % APV- 2

% The experimental condition of a cell can be uniquely determined by the date, 
% and the cell_num assigned on that day via a look-up cell_id_index table.  

%% %% Change these accordingly based on how you want to group the data
% in each anaylzed_spFR_results file, dates can be extracted from the last six digits

%save results
save_results = 1;

%experiment name (correpsonding to the sheet name in the cell_id_index
% excel file)
exp_name = 'APV_acute';

%file name
saved_file_name = 'APV_acute.mat';

%where to save grouped files
fp_grouped_data = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/spFR_data_by_groups';

%experimental conditions
exp_con = {'NT_Ctrl','APV_acute'};

%import cell_id_index table 
cd('/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/spFR_data_by_groups')
cell_id_index = readtable('cell_id_index.xlsx','Sheet',exp_name);

%location of analyzed spFR results
fp_analyzed_spFR = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/analyzed_spFR_results/';

%% extract and group data
data_struct_temp = cell(1,numel(exp_con));

cd(fp_analyzed_spFR)

all_files_temp = dir;

for fi = 1:size(all_files_temp,1)
    if strcmp('.DS_Store',all_files_temp(fi).name)
        delete '.DS_Store'
        continue
    end
end

all_files = dir;
file_num = numel(all_files)-2;

for fi = 1:file_num
    curr_name = all_files(fi+2).name;
    if ismember('test',curr_name(6:11))
        continue
    else
        curr_date = curr_name(6:11);
    end
    
    
    cy = strcat('20',curr_date(1:2));
    cM = curr_date(3:4);
    cday = curr_date(5:6);
    
    date = datetime(strcat(cy,'-',cM,'-',cday),'InputFormat','y-M-d',...
        'Format', 'M/d/y');
    if ~(ismember(date,cell_id_index.Date))
        continue
    else
        load(all_files(fi+2).name)
        
        for ci = 1:size(cell_id,2)
            if isempty(cell_id{1,ci})
                continue
            else
                if isempty(find(cell_id_index.Date == date & cell_id_index.Cell_num == ci,1))
                    continue
                else
                
                    row = find(cell_id_index.Date == date & cell_id_index.Cell_num == ci);
                    cond_i = cell_id_index{row,'Cat'};
                    cell_ID = cell_id_index{row,'Cell_ID'};
                    ti_start = cell_id{1,ci}(1,1);
                    trace_ct = size(cell_id{1,ci},1);
                    ti_end = cell_id{1,ci}(trace_ct,1);
                    
                    data_struct_temp{cond_i}.meanFR(cell_ID,1) = meanFR(ci,1);  
                    data_struct_temp{cond_i}.total_recording_time(cell_ID,1) = total_recording_time(ci,1);

                    c_V_th = NaN(1,numel(V_th{ci}));
                    c_AP_peak = NaN(1,numel(V_th{ci}));

                    for ti = 1:numel(V_th{ci})
                        c_V_th(ti) = mean(V_th{ci}{ti},'omitnan');
                        c_AP_peak(ti) = mean(AP_peak{ci}{ti},'omitnan');
                    end

                    data_struct_temp{cond_i}.Vth(cell_ID,1) = mean(c_V_th,'omitnan');
                    data_struct_temp{cond_i}.AP_peak(cell_ID,1) = mean(c_AP_peak,'omitnan');
       
                    %for Rin and Cm
                    if isnan(cell_stats.Rin_stats(ci).rev_ave)                      
                        data_struct_temp{cond_i}.Rin(cell_ID,1) = cell_stats.Rin_stats(ci).ave;
                    else
                        data_struct_temp{cond_i}.Rin(cell_ID,1) = cell_stats.Rin_stats(ci).rev_ave;
                    end
                    
                    if isnan(cell_stats.Cm_stats(ci).rev_ave)                    
                        data_struct_temp{cond_i}.Cm(cell_ID,1) = cell_stats.Cm_stats(ci).ave;
                    else
                        data_struct_temp{cond_i}.Cm(cell_ID,1) = cell_stats.Cm_stats(ci).rev_ave;
                    end 

                    if isnan(cell_stats.Vm_stats(ci).rev_ave)                    
                        data_struct_temp{cond_i}.Vm(cell_ID,1) = cell_stats.Vm_stats(ci).ave;
                    else
                        data_struct_temp{cond_i}.Vm(cell_ID,1) = cell_stats.Vm_stats(ci).rev_ave;
                    end 

                end     
            end
        end
    end
end
 
%%%% CHANGE FOR EACH RUN %%%%

NT_Ctrl = data_struct_temp{1,1};
APV_acute = data_struct_temp{1,2};

%% save grouped data
if save_results == 1
    cd(fp_grouped_data)
    save(saved_file_name,exp_con{1,:},'exp_con','cell_id_index')
end
