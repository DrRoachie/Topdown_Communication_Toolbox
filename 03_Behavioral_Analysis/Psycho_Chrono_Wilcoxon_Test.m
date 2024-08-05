
% Performs Wilcoxon rank-sum test on psychometric/chronometric curves for
% specified animal and condition (script must be run 4 times to perform
% test on all data: 2 animals, 2 conditions)
% 
% Input: Run Psychometric_Chronometric_Analysis.m first on desired 
% animal/condition to generate psychometric/chronometric data arrays
% (total_proportion_H and total_RTs). 
% 
% Output: Table called T holds all p-values organized by SNR ratio (rows).
% Columns are labeled and correspond to a psychometric or chronometric
% curve. Each animal outputs a table, so T must be saved for each animal
% before running the test on the next animal. 
% 
% Execution: Each part of the 2 if-else statements must be reached to
% generate a complete table of p-values. The if statement in the 
% psychometric and chronometric chunks is reached after running
% Psychometric_Chronometric_Analysis.m with the 'prior' condition. The
% elseif statement is reached after running
% Psychometric_Chronometric_Analysis.m with the 'pretone' condition. 
% 
% For example, to output a complete table for MrCassius: 
% 1. Run Psychometric_Chronometric_Analysis.m in 'prior' condition 
% 2. Run Psycho_Chrono_Wilcoxon_Test.m 
% 3. Run Psychometric_Chronometric_Analysis.m in 'pretone' condition 
% 4. Run Psycho_Chrono_Wilcoxon_Test.m 
% 5. Save table T

%% Rank sum test on psychometric data

if ~exist('total_proportion_H', 'var')
    error("Missing the variable 'total_proportion_H'. Unable to run the Wilcoxon rank sum test on psychometric data.");
end

if length(total_proportion_H(:,1,1)) == 3   % total_proportion_H has 3 rows in the prior only condition
    % run test
    PriorOnly_Psycho_p_high = [];
    PriorOnly_Psycho_p_low  = [];

    for i = 1:length(total_proportion_H(1,:,1)) % for every SNR ratio
        % separate mean proportion H values by high, low, neutral prior
        high_data       = squeeze(total_proportion_H(1,i,:));
        low_data        = squeeze(total_proportion_H(2,i,:));
        neutral_data    = squeeze(total_proportion_H(3,i,:));

        PriorOnly_Psycho_p_high(i)  = ranksum(high_data, neutral_data);
        PriorOnly_Psycho_p_low(i)   = ranksum(low_data, neutral_data);
    end

elseif length(total_proportion_H(:,1,1)) == 2   % total_proportion_H has 2 rows in the pretone only condition
    % run test
    PretoneOnly_Psycho_p = [];

    for i = 1:length(total_proportion_H(1,:,1)) % for every SNR ratio
        % separate mean proportion H values by high, low, neutral prior
        high_data       = squeeze(total_proportion_H(1,i,:));
        low_data        = squeeze(total_proportion_H(2,i,:));

        PretoneOnly_Psycho_p(i)  = ranksum(high_data, low_data);
    end
end


%% Rank sum test on chronometric data

if ~exist('total_RTs', 'var')
    error("Missing the variable 'total_RTs'. Unable to run the Wilcoxon rank sum test on chronometric data.");
end

if length(total_RTs(:,1,1)) == 3   % total_RTs has 3 rows in the prior only condition
    % run test
    PriorOnly_Chrono_p_high = [];
    PriorOnly_Chrono_p_low  = [];

    for i = 1:length(total_RTs(1,:,1)) % for every SNR ratio
        % separate mean proportion H values by high, low, neutral prior
        low_data        = squeeze(total_RTs(1,i,:));
        neutral_data    = squeeze(total_RTs(2,i,:));
        high_data       = squeeze(total_RTs(3,i,:));

        PriorOnly_Chrono_p_high(i)  = ranksum(high_data, neutral_data);
        PriorOnly_Chrono_p_low(i)   = ranksum(low_data, neutral_data);
    end

elseif length(total_RTs(:,1,1)) == 2   % total_RTs has 2 rows in the pretone only condition
    % run test
    PretoneOnly_Chrono_p = [];

    for i = 1:length(total_RTs(1,:,1)) % for every SNR ratio
        % separate mean proportion H values by high, low, neutral prior
        high_data       = squeeze(total_RTs(1,i,:));
        low_data        = squeeze(total_RTs(2,i,:));

        PretoneOnly_Chrono_p(i)  = ranksum(high_data, low_data);
    end
end

%% create table of saved p-values

SNR_names = strrep(cellstr(num2str(SNR_list)),' ','');

T = table(PriorOnly_Psycho_p_high', PriorOnly_Psycho_p_low', PriorOnly_Chrono_p_high', PriorOnly_Chrono_p_low', ...
          PretoneOnly_Psycho_p', PretoneOnly_Chrono_p', 'VariableNames', ...
          {'PriorOnly_Prop_H_high_v_neutral', 'PriorOnly_Prop_H_low_v_neutral', ...
           'PriorOnly_RT_high_v_neutral', 'PriorOnly_RT_low_v_neutral', ...
           'PretoneOnly_Prop_H_high_v_low', 'PretoneOnly_RT_high_v_low'}, ...
           'RowNames', cellfun(@(x) sprintf('SNR_%s', x), SNR_names, 'UniformOutput', false));
