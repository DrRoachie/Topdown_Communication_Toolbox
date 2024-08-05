
% plots summary figure of spectrogram, coherence + statistical test,
% granger + statistical test, and cross correlation. 

close all;

% Define data

rootdir = 'F:\2024_07_24_Analysis';
sessions = dir(fullfile(rootdir, '19*'));

for i = 1:length(sessions)

    RecDate = sessions(i).name;

    animals = {'MrCassius', 'MrM'};

    for j = 1:length(animals)

        Animal = animals{j};
        Epoch  = 'testToneOnset';

        figFiles = dir(fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), '*Coh_Spec.fig'));

        for k = 1:length(figFiles)
            
            tokens = regexp(figFiles(k).name, '(.*?)_', 'tokens');
            if isempty(tokens)
                continue;
            end

            Send_Cort   = tokens{4}{1}; 
            Send_Num    = tokens{5}{1};
            Rec_Cort    = tokens{6}{1};
            Rec_Num     = tokens{7}{1};

            Send_Chan   = [Send_Cort '_' Send_Num];
            Rec_Chan    = [Rec_Cort '_' Rec_Num];

            Freq_Band   = tokens{8}{1};
            Behavior    = tokens{9}{1};    
                      

            datadir     = fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'));

            % plot summary figure
            summary_fig = figure;
            
            % super title
            Send_Chan_lb   = strrep(Send_Chan,'_ch','');
            Rec_Chan_lb    = strrep(Rec_Chan,'_ch','');
            
            sgtitle(sprintf('%s %s %s-%s Connectivity Analysis in %s', Animal, RecDate, Send_Chan_lb, Rec_Chan_lb, Freq_Band));
            
            
            % Prior Spectrogram
            file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Coh_Prior_Spectrogram.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);
            
            h = openfig(fullfile(datadir, file_name),'invisible');
            axesHandles = findall(h, 'type', 'axes');
            ax1 = axesHandles(2);
            ax2 = axesHandles(1);
            
            figure(summary_fig);
            s1 = subplot(4,8,1);
            s2 = subplot(4,8,2);
            
            fig1 = get(ax1,'children');
            fig2 = get(ax2,'children');
            
            copyobj(fig1,s1);
            copyobj(fig2,s2);
            
            % fix x and y limits
            set(s1, 'XLim', get(ax1, 'XLim'));
            set(s1, 'YLim', get(ax1, 'YLim'));
            set(s1, 'CLim', get(ax1, 'CLim'));
            set(s1, 'View', get(ax1, 'View'));
            set(s2, 'XLim', get(ax2, 'XLim'));
            set(s2, 'YLim', get(ax2, 'YLim'));
            set(s2, 'CLim', get(ax2, 'CLim'));
            set(s2, 'View', get(ax2, 'View'));
            
            % add labels and titles
            xlabel(s1, get(get(ax1,'xlabel'),'string'));
            ylabel(s1, get(get(ax1,'ylabel'),'string'));
            title(s1, 'Prior PFC Spec');
            xlabel(s2, get(get(ax1,'xlabel'),'string'));
            ylabel(s2, get(get(ax1,'ylabel'),'string'));
            title(s2, 'Prior AC Spec');
            
            
            % Pretone Spectrogram
            file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Coh_Pretone_Spectrogram.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);
            
            h = openfig(fullfile(datadir, file_name), 'invisible');
            axesHandles = findall(h, 'type', 'axes');
            ax1 = axesHandles(2);
            ax2 = axesHandles(1);
            
            figure(summary_fig);
            s1 = subplot(4,8,3);
            s2 = subplot(4,8,4);
            
            fig1 = get(ax1,'children');
            fig2 = get(ax2,'children');
            
            copyobj(fig1,s1);
            copyobj(fig2,s2);
            
            % fix x and y limits
            set(s1, 'XLim', get(ax1, 'XLim'));
            set(s1, 'YLim', get(ax1, 'YLim'));
            set(s1, 'CLim', get(ax1, 'CLim'));
            set(s1, 'View', get(ax1, 'View'));
            set(s2, 'XLim', get(ax2, 'XLim'));
            set(s2, 'YLim', get(ax2, 'YLim'));
            set(s2, 'CLim', get(ax2, 'CLim'));
            set(s2, 'View', get(ax2, 'View'));
            
            % add labels and titles
            xlabel(s1, get(get(ax1,'xlabel'),'string'));
            ylabel(s1, get(get(ax1,'ylabel'),'string'));
            title(s1, 'Pretone PFC Spec');
            xlabel(s2, get(get(ax1,'xlabel'),'string'));
            ylabel(s2, get(get(ax1,'ylabel'),'string'));
            title(s2, 'Pretone AC Spec');
            
            colormap jet;
            
            
            % Coherence Spectrum
            file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Coh_Spec.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);
            
            openfig(fullfile(datadir, file_name),'invisible');
            ax1 = gca;
            
            figure(summary_fig);
            s1 = subplot(4,6,7);
            
            fig1 = get(ax1,'children');
            
            copyobj(fig1,s1);
            
            % fix x and y limits
            x_tick_label = str2double(get(ax1, 'XTickLabel'));
            set(s1, 'XLim', get(ax1, 'XLim'));
            set(s1, 'XTickLabel', x_tick_label);
            set(s1, 'YLim', get(ax1, 'YLim'));
            set(s1, 'CLim', get(ax1, 'CLim'));
            set(s1, 'View', get(ax1, 'View'));
            
            % add labels and titles
            xlabel(s1, get(get(ax1,'xlabel'),'string'));
            ylabel(s1, get(get(ax1,'ylabel'),'string'));
            title(s1, 'Coherence');
            
            
            % Coherence zthresh
            file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Coh_zthresh.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);
            
            openfig(fullfile(datadir, file_name),'invisible');
            ax1 = gca;
            
            figure(summary_fig);
            s1 = subplot(4,6,8);
            
            fig1 = get(ax1,'children');
            
            copyobj(fig1,s1);
            
            % fix x and y limits
            x_tick_label = str2double(get(ax1, 'XTickLabel'));
            set(s1, 'XLim', get(ax1, 'XLim'));
            set(s1, 'XTickLabel', x_tick_label);
            set(s1, 'YLim', get(ax1, 'YLim'));
            set(s1, 'CLim', get(ax1, 'CLim'));
            set(s1, 'View', get(ax1, 'View'));
            
            % add labels and titles
            xlabel(s1, get(get(ax1,'xlabel'),'string'));
            ylabel(s1, get(get(ax1,'ylabel'),'string'));
            title(s1, 'Coherence z threshold');
            
            
            % Coherence montecarlo
            file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Coh_montecarlo.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);
            
            openfig(fullfile(datadir, file_name), 'invisible');
            ax1 = gca;
            
            figure(summary_fig);
            s1 = subplot(4,6,9);
            
            fig1 = get(ax1,'children');
            
            copyobj(fig1,s1);
            
            % fix x and y limits
            set(s1, 'XLim', get(ax1, 'XLim'));
            set(s1, 'YLim', get(ax1, 'YLim'));
            set(s1, 'CLim', get(ax1, 'CLim'));
            set(s1, 'View', get(ax1, 'View'));
            
            % add labels and titles
            xlabel(s1, get(get(ax1,'xlabel'),'string'));
            ylabel(s1, get(get(ax1,'ylabel'),'string'));
            title(s1, 'Coherence montecarlo');
            
            % add p value
            file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Coh_data.mat', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);
            load(fullfile(datadir, file_name), 'pvalue_c');
            p_value_text = sprintf('p = %.4f', pvalue_c.(Freq_Band));
            annotation('textbox', [0.46, 0.58, 0.1, 0.1], 'String', p_value_text, 'EdgeColor', 'black', 'BackgroundColor', 'white', 'FontSize', 8, 'FitBoxToText', 'on');
            
            
            % Granger Spectrum PFC to AC
            file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Granger_Spec_PFC_to_AC.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);
            
            h = openfig(fullfile(datadir, file_name), 'invisible');
            ax1 = gca;
            
            figure(summary_fig);
            s1 = subplot(4,6,13);
            
            fig1 = get(ax1,'children');
            
            copyobj(fig1,s1);
            
            % fix x and y limits
            set(s1, 'XLim', get(ax1, 'XLim'));
            set(s1, 'YLim', get(ax1, 'YLim'));
            set(s1, 'CLim', get(ax1, 'CLim'));
            set(s1, 'View', get(ax1, 'View'));
            
            % add labels and titles
            xlabel(s1, get(get(ax1,'xlabel'),'string'));
            ylabel(s1, get(get(ax1,'ylabel'),'string'));
            title(s1, 'Granger PFC-to-AC');
            
            
            % Granger zthresh PFC to AC
            file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Granger_zthresh_PFC_to_AC.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);
            
            openfig(fullfile(datadir, file_name), 'invisible');
            ax1 = gca;
            
            figure(summary_fig);
            s1 = subplot(4,6,14);
            
            fig1 = get(ax1,'children');
            
            copyobj(fig1,s1);
            
            % fix x and y limits
            x_tick_label = str2double(get(ax1, 'XTickLabel'));
            set(s1, 'XLim', get(ax1, 'XLim'));
            set(s1, 'XTickLabel', x_tick_label);
            set(s1, 'CLim', get(ax1, 'CLim'));
            set(s1, 'View', get(ax1, 'View'));
            
            % add labels and titles
            xlabel(s1, get(get(ax1,'xlabel'),'string'));
            ylabel(s1, get(get(ax1,'ylabel'),'string'));
            title(s1, 'Granger z thresh PFC-to-AC');
            
            
            % Granger montecarlo PFC to AC
            file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Granger_montecarlo_PFC_to_AC.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);
            
            openfig(fullfile(datadir, file_name), 'invisible');
            ax1 = gca;
            
            figure(summary_fig);
            s1 = subplot(4,6,15);
            
            fig1 = get(ax1,'children');
            
            copyobj(fig1,s1);
            
            % fix x and y limits
            set(s1, 'XLim', get(ax1, 'XLim'));
            set(s1, 'YLim', get(ax1, 'YLim'));
            set(s1, 'CLim', get(ax1, 'CLim'));
            set(s1, 'View', get(ax1, 'View'));
            
            % add labels and titles
            xlabel(s1, get(get(ax1,'xlabel'),'string'));
            ylabel(s1, get(get(ax1,'ylabel'),'string'));
            title(s1, 'Granger montecarlo PFC-to-AC');
            
            % add p value
            file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Granger_data.mat', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);
            load(fullfile(datadir, file_name), 'pvalue_PFC_AC');
            p_value_text = sprintf('p = %.4f', pvalue_PFC_AC.(Freq_Band));
            annotation('textbox', [0.46, 0.37, 0.1, 0.1], 'String', p_value_text, 'EdgeColor', 'black', 'BackgroundColor', 'white', 'FontSize', 8, 'FitBoxToText', 'on');
            
            
            % Granger Spectrum AC to PFC
            file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Granger_Spec_AC_to_PFC.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);
            
            openfig(fullfile(datadir, file_name), 'invisible');
            ax1 = gca;
            
            figure(summary_fig);
            s1 = subplot(4,6,19);
            
            fig1 = get(ax1,'children');
            
            copyobj(fig1,s1);
            
            % fix x and y limits
            set(s1, 'XLim', get(ax1, 'XLim'));
            set(s1, 'YLim', get(ax1, 'YLim'));
            set(s1, 'CLim', get(ax1, 'CLim'));
            set(s1, 'View', get(ax1, 'View'));
            
            % add labels and titles
            xlabel(s1, get(get(ax1,'xlabel'),'string'));
            ylabel(s1, get(get(ax1,'ylabel'),'string'));
            title(s1, 'Granger AC-to-PFC');
            
            
            % Granger zthresh AC to PFC
            file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Granger_zthresh_AC_to_PFC.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);
            
            openfig(fullfile(datadir, file_name), 'invisible');
            ax1 = gca;
            
            figure(summary_fig);
            s1 = subplot(4,6,20);
            
            fig1 = get(ax1,'children');
            
            copyobj(fig1,s1);
            
            % fix x and y limits
            x_tick_label = str2double(get(ax1, 'XTickLabel'));
            set(s1, 'XLim', get(ax1, 'XLim'));
            set(s1, 'XTickLabel', x_tick_label);
            set(s1, 'YLim', get(ax1, 'YLim'));
            set(s1, 'CLim', get(ax1, 'CLim'));
            set(s1, 'View', get(ax1, 'View'));
            
            % add labels and titles
            xlabel(s1, get(get(ax1,'xlabel'),'string'));
            ylabel(s1, get(get(ax1,'ylabel'),'string'));
            title(s1, 'Granger z thresh AC-to-PFC');
            
            
            % Granger montecarlo AC to PFC
            file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Granger_montecarlo_AC_to_PFC.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);
            
            openfig(fullfile(datadir, file_name), 'invisible');
            ax1 = gca;
            
            figure(summary_fig);
            s1 = subplot(4,6,21);
            
            fig1 = get(ax1,'children');
            
            copyobj(fig1,s1);
            
            % fix x and y limits
            set(s1, 'XLim', get(ax1, 'XLim'));
            set(s1, 'YLim', get(ax1, 'YLim'));
            set(s1, 'CLim', get(ax1, 'CLim'));
            set(s1, 'View', get(ax1, 'View'));
            
            % add labels and titles
            xlabel(s1, get(get(ax1,'xlabel'),'string'));
            ylabel(s1, get(get(ax1,'ylabel'),'string'));
            title(s1, 'Granger montecarlo AC-to-PFC');
            
            % add p value
            file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Granger_data.mat', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);
            load(fullfile(datadir, file_name), 'pvalue_AC_PFC');
            p_value_text = sprintf('p = %.4f', pvalue_AC_PFC.(Freq_Band));
            annotation('textbox', [0.46, 0.15, 0.1, 0.1], 'String', p_value_text, 'EdgeColor', 'black', 'BackgroundColor', 'white', 'FontSize', 8, 'FitBoxToText', 'on');
            
            
            % XCorr
            file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_XCorr.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);
            
            openfig(fullfile(datadir, file_name), 'invisible');
            ax1 = gca;
            
            figure(summary_fig);
            s1 = subplot(2,2,2);
            
            fig1 = get(ax1,'children');
            
            copyobj(fig1,s1);
            
            % fix x and y limits
            set(s1, 'XLim', get(ax1, 'XLim'));
            set(s1, 'YLim', get(ax1, 'YLim'));
            set(s1, 'CLim', get(ax1, 'CLim'));
            set(s1, 'View', get(ax1, 'View'));
            
            % add labels and titles
            xlabel(s1, get(get(ax1,'xlabel'),'string'));
            ylabel(s1, get(get(ax1,'ylabel'),'string'));
            title(s1, get(get(ax1,'title'),'string'));

        end
    end
end


