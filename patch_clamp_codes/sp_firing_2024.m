%%%% This script loops over spontaneous firing traces from all cells
%%%% Each trace includes a 0.5s hyperpolarizing step of -50mV as the
%%%% seal test (timestamp 0.5s to 1s) and a ~30s recording of spontaneous
%%%% firing under the current clamp

%%%% raw traces are excised to 

%%%% Assess the baseline quality for each trace
%%%% Calculate passive properties from each trace and assess the average
%%%% for a cell

%% Experimental settings

%Name of the data folder
experiment = '241220';

%location of the excised file
excised_data_fp = '/Users/wwneuro/My_Drive/Lab/Data/culture_experiments/excised_data/';

%Filepath where you keep data folders
fp_data = '/Users/wwneuro/My_Drive/Lab/Data/culture_experiments/spontaneousFR/';

%location to save the analyzed results
fp_analyzed_data = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/analyzed_spFR_results/';

%location to save the analyzed passive property data
fp_pp = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/analyzed_seal_test/';

%whether to save results
save_results = 1;

%whether to show figures for each spFR trace; 0 = no, 1 = yes.
figures_on_spFR = 0;

%whether to show figures for passive property analysis; 0 = no, 1 = yes.
figures_on_pp = 0;

%% Parameters

%sampling rate (in Hz)
sample_rate = 10000;

%holding potential (in mV)
v_hold = -50;

%%% seal test settings
%%latency of pulse onset (in s)
step_start_pp = 0.5;

%duration of the seal test pulse (in seconds)
pulse_pp = 0.5;

%pulse amplitude (in A)
I_step = -50*10^-12;

%pp CV cutoff (in % of the mean)
pp_cv_cutoff = 15;

%%% spike identification criteria
%dV/dt (in V/s, must be higher than this value)
dVdt_low = 20;

%peak amplitude cutoff (in mV, must be higher than this value)
spike_peak_low = -30;

%spike threshold cutoff (in mV, must be lower than this value)
spike_thre_high = -20;

%% Data readout 

%file names of the analyzed data were they to be saved
save_name_spFR = strcat('spFR_', experiment, '.mat'); %spFR
save_name_pp = strcat('pp_', experiment, '.mat'); %pp

%import excised regions
excised_file_name = strcat('excised_',experiment,'.mat');
excised_data = load(strcat(excised_data_fp,excised_file_name));
excised_points = excised_data.excised_points;

current_folder = strcat(fp_data, experiment);
[raw_h5_files, data, cell_id, cell_num, filename] = h5_file_readout(current_folder);

%data structure initialization

aDAT = cell(1, max(cell_num));
spike_index = cell(1,max(cell_num));
spike_count = cell(1, max(cell_num));
duration = cell(1, max(cell_num));
FR = cell(1,max(cell_num));
V_th = cell(1, max(cell_num));
AP_peak = cell(1, max(cell_num));
total_recording_time = NaN(max(cell_num),1);
meanFR = NaN(max(cell_num),1);

Rin = cell(1, max(cell_num));
Cm = cell(1, max(cell_num));
Vm = cell(1, max(cell_num));

Rin_stats =struct('ave','NaN','stdev','NaN','CV','NaN','ACCEPT','NaN',...,
    'rev_ave','NaN','rev_stdev','NaN','rev_CV','NaN');
Cm_stats = struct('ave','NaN','stdev','NaN','CV','NaN','ACCEPT','NaN',...,
    'rev_ave','NaN','rev_stdev','NaN','rev_CV','NaN');
Vm_stats = struct('ave','NaN','stdev','NaN','CV','NaN','ACCEPT','NaN',...,
    'rev_ave','NaN','rev_stdev','NaN','rev_CV','NaN');

Rin_outliers = cell(1,max(cell_num));
Cm_outliers = cell(1,max(cell_num));
Vm_outliers = cell(1,max(cell_num));

%% spike detection and analysis

for ci = 1:max(cell_num)
    if isempty(cell_id{ci})
        continue
    else
        trace_start = cell_id{ci}(1,1);
        trace_end = cell_id{ci}(end,1);

        for ti = trace_start:trace_end
            if ~ismember(ti, cell_id{ci})
                disp(strcat('Cell',num2str(ci),' Trace',num2str(ti),' does not exist!'))

                aDAT{ci}{ti} = [];
                spike_index{ci}{ti} = [];
                spike_count{ci}(ti,1) = NaN;
                duration{ci}(ti,1) = NaN;
                V_th{ci}{ti} = [];
                AP_peak{ci}{ti} = [];
                Rin{ci}(ti,1) = NaN;
                Cm{ci}(ti,1) = NaN;
                Vm{ci}(ti,1) = NaN; %the holding membrane potential

                continue
            else

                vals = nonzeros(data{ci}(:,ti));
                aDAT{ci}{ti} = vals;

                %calculate passive properties (data should be smoothed)
                window = 100;
                sm_vals_pp = smoothdata(vals(1:4*sample_rate),'gaussian',window);

                [Rin_curr, Cm_curr, ~, ~, ~] ...,
                    = get_PP_I_clamp(sm_vals_pp, step_start_pp, pulse_pp, I_step, sample_rate, figures_on_pp);

                Rin{ci}(ti,1) = Rin_curr;
                Cm{ci}(ti,1) = Cm_curr;

                %get Vm from the rest of the trace
                start_end = int32(excised_points{1}{ci}{ti});
                if start_end(2) > numel(vals)
                    sm_vals_sp = smoothdata(vals(start_end(1):end),'gaussian',window);
                else
                    sm_vals_sp = smoothdata(vals(start_end(1):start_end(2)),'gaussian',window);
                end
                Vm_curr = mean(sm_vals_sp,'omitnan');
                Vm{ci}(ti,1) = Vm_curr;

                %Use raw data for spike detection
                if start_end(2)>numel(vals)
                    vals_sp = vals(start_end(1):end);
                else
                    vals_sp = vals(start_end(1):start_end(2));
                end

                dV = diff(vals_sp);
                dV_sec = smooth(dV.*(sample_rate/1000)); % in V/s

                sp_ct = 0;
                spike_ind = struct;
                di_last = 1;
               
                for di = 1:numel(vals_sp)-1
                    if dV_sec(di) > dVdt_low && dV_sec(di+1) < dVdt_low && vals_sp(di) > spike_peak_low
                       
                        %peak
                        if di+100 > numel(dV_sec)
                            continue
                        else
                            pi = find(dV_sec(di:di+100) <= 0, 1, 'first')+di-1;
                        end

                        if sp_ct > 1 && pi == spike_index{ci}{ti}(sp_ct).peak
                            continue
                        else
                            sp_ct = sp_ct+1;
                            spike_ind(sp_ct).peak = pi;
                            spike_index{ci}{ti} = spike_ind;
                            AP_peak{ci}{ti}(sp_ct) = vals_sp(pi);
                        end

                        %threshold

                        if isempty(find(dV_sec(di_last:di) <= 5,1,'last'))
                            thi = 2+di_last;
                        else
                            thi = find(dV_sec(di_last:di) <= 5,1,'last')+ di_last;
                        end

                        spike_ind(sp_ct).threshold = thi;
                        V_th{ci}{ti}(sp_ct) = vals_sp(thi);

                        di_last = di;
                    end
                end
            end

            if figures_on_spFR == 1
                figure('position',[119 171 1000 611])
                hold on
                plot(vals_sp)

                if isfield(spike_ind, 'threshold')
                
                    for ppi = 1:numel(spike_ind)
                        thre = spike_ind(ppi).threshold;
                        pk = spike_ind(ppi).peak;
                        scatter(thre,vals_sp(thre),'ro','filled')
                        scatter(pk,vals_sp(pk),'ko','filled')
                    end
                end
                hold off
    
                title(strcat('Cell', num2str(ci), ' Trace', num2str(ti)))
            end
    
            %spike_index{ci}{ti} = spike_ind;
            spike_count{ci}(ti,1) = sp_ct;
            duration{ci}(ti,1) = numel(vals_sp)/sample_rate; 
            FR{ci}(ti,1) = spike_count{ci}(ti,1)/duration{ci}(ti,1);

        end %trace

    end
        
    % Passive Property Evalutaion
    Rin_curr = nonzeros(Rin{ci});
    Cm_curr = nonzeros(Cm{ci});
    Vm_curr = nonzeros(Vm{ci});
    
    [Rin_stats_temp,Rin_outliers_temp] = pp_eval(Rin_curr,pp_cv_cutoff);
    [Cm_stats_temp,Cm_outliers_temp] = pp_eval(Cm_curr,pp_cv_cutoff);
    [Vm_stats_temp,Vm_outliers_temp] = pp_eval(Vm_curr,pp_cv_cutoff);

    Rin_stats(ci) = Rin_stats_temp;
    Cm_stats(ci) = Cm_stats_temp;
    Vm_stats(ci) = Vm_stats_temp;
    
    Rin_outliers{ci} = Rin_outliers_temp;
    Cm_outliers{ci} = Cm_outliers_temp;
    Vm_outliers{ci} = Vm_outliers_temp;

    %assess deviation of Vm from the holding Vm
    %CV < 10%, deviation from V_hold < 5%
    if Vm_stats(ci).ACCEPT == 1
        Vm_dev = abs(Vm_stats(ci).ave-v_hold)/abs(v_hold)*100;
        Vm_CV = Vm_stats(ci).CV;
    else
        Vm_dev = abs(Vm_stats(ci).rev_ave-v_hold)/abs(v_hold)*100;
        Vm_CV = Vm_stats(ci).rev_CV;
    end

    if Vm_dev>5 || Vm_CV > 10
        Vm_stats(ci).ACCEPT = 0;
    else
        Vm_stats(ci).ACCEPT = 1;
    end

    %total recording time (in minutes)
    total_recording_time(ci,1) = sum(duration{ci},'omitnan')/60;

    %mean FR for the cell
    meanFR(ci,1) = mean(nonzeros(FR{ci}),'omitnan');
        
        
end %cell

%save PP stats
outliers.Rin = Rin_outliers;
outliers.Cm = Cm_outliers;
outliers.Vm = Vm_outliers;
   
cell_stats.Rin_stats = Rin_stats;
cell_stats.Cm_stats = Cm_stats;
cell_stats.Vm_stats = Vm_stats;
cell_stats.outliers = outliers;

%% save results
if save_results

    %%% save a copy for passive properties only
    cd(fp_pp)
    save(save_name_pp, 'Vm','Vm_stats','Rin','Cm','Rin_stats',...,
        'Cm_stats', 'outliers')

    %%% save everything else in the analyzed data file
    cd(fp_analyzed_data)
    save(save_name_spFR,'cell_stats','cell_id','aDAT','spike_index','spike_count',...,
        'duration','V_th','AP_peak','total_recording_time','Rin','Cm','Vm','FR','meanFR')
end
