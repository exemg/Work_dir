%% Part 1: Service variables
% NG = selection of the fit type: linear bkg with 1,2,3 gaussians (NG = 1,2,3)
NG = 1;
% valcond = Specifies whether the plots refer to ROI averages values (if=0)  
% or to single-pixel values(if=1)
valcond = 1; 
% wdthmax = Limit value above which a bead width is discarded (in units of the z-stack)
wdthmax = 30;
% rawnetbit = 'raw' or 'net' /  small string bit to indicate whether or not bkg subtraction has been performed prior to analysis
rawnetbit = 'net'; 
% phys = factor containing the physical length of a pixel size
phys = 0.2/0.86;

%% Part 2: Preliminary folder investigation
% From this section, we expect to:
%  1) derive a list of the CSV files we require
%  2) create an analysis folder with the date, time, NG number and analysis
%  mode
if valcond == 0
    if strcmpi(rawnetbit,'raw')
        dirlist_csv = ls('*Zprof*avg_raw*.csv');  % list of .csv files with the "*yyy" particles in the title 
        suffix = strcat('_avg_raw_NG',num2str(NG));
    elseif strcmpi(rawnetbit,'net')
        dirlist_csv = ls('*Zprof*avg_net*.csv');
        suffix = strcat('_avg_net_NG',num2str(NG));
    else
        error('Check the rawnetbit variable, as it has a nonvalid argument');
    end
elseif valcond == 1
    if strcmpi(rawnetbit,'raw')
        dirlist_csv = ls('*Zprof*pks_raw*.csv');  % list of .csv files with the "frm" bit in their title
        suffix = strcat('_pks_raw_NG',num2str(NG));
    elseif strcmpi(rawnetbit,'net')
        dirlist_csv = ls('*Zprof*pks_net_2*.csv');
        suffix = strcat('_pks_net_NG',num2str(NG));
    else
        error('Check the rawnetbit variable, as it has a nonvalid argument');
    end
else
    disp('invalid valcond variable');
    return
end
% Ncsv = size of the list, Nchar = maximum size of a title string
[Ncsv, Nchar] = size(dirlist_csv);

% creates the list of the save folders that will be used to store the
% analysis data
dtst = datestr(now,'yyyymmdd_HH_MM');
sv_pth_gen = strcat('Analysis',suffix,'_',dtst,'\'); % path of the analysis folder ("general" path)
[statusmkdirgen,msggen] = mkdir(sv_pth_gen); % creates the folder (if it fails, the reason is explained in "msggen"
% We need to create as many analysis folders as there ar .csv files ("local folders")
if Ncsv == 1
    pth_loc_id = dirlist_csv(end-19:end-4);
    sv_pth_loc_list = cat(2,sv_pth_gen,pth_loc_id,'\');
else
    pth_loc_id = dirlist_csv(:,end-19:end-4);
    sv_pth_loc_list = cat(2,repmat(sv_pth_gen,Ncsv,1),pth_loc_id,repmat('\',Ncsv,1));
end

% For each ROI CSV table:
%  1) takes the data, X and Y
%  2) makes a single-table analysis report for each csv, together with plots
%  3) saves figures and tables in each of the "sv_pth_loc_list" subfolders
%  4) keeps track of the averages and errors of the overall beads results and saves
%    them on a separate table ("gentable") in the MAIN folder ("sv_pth_gen")
gentable = zeros(Ncsv,3*(3*NG-2));  % sizes are given by: (avg+err+std [3]) X (wdth [NG] + Dz [NG-1] + wdr [NG-1])

%% Part 3: Analysis of the .csv files
for ixT = 1:Ncsv
    %% Part 4: Create local folder
    prefix = sv_pth_loc_list(ixT,:); % "prefix" gives the path where to save
    foldname = pth_loc_id(ixT,:);    % local folder NAME (not path)
    [statusmkdir,msg] = mkdir(sv_pth_gen, foldname); % creates the local folder (failure is explained in "msg")
    
    %% Part5: Data analysis preparation
    Data = csvread(dirlist_csv(ixT,:),1,0);  % gets the datsa from the csv
    
    % Definition of the fit type: linear bkg with 1, 2, or 3 gaussians
    % The S# parameters are exactly the FWHMs thanks to the 1.6651 factor. 
    % Use 0.8326 if you want the HWHM.
    if NG == 3 % three peaks
        fitType = fittype('A + B*x + C1*exp(-(1.6651*(x-X1)/S1)^2) + C2*exp(-(1.6651*(x-X2)/S2)^2) + C3*exp(-(1.6651*(x-X3)/S3)^2)',...
            'dependent',{'y'},'independent',{'x'},...
            'coefficients',{'A','B','C1','S1','X1','C2','S2','X2','C3','S3','X3'});
    elseif NG == 2 % two peaks
        fitType = fittype('A + B*x + C1*exp(-(1.6651*(x-X1)/S1)^2) + C2*exp(-(1.6651*(x-X2)/S2)^2)',...
            'dependent',{'y'},'independent',{'x'},...
            'coefficients',{'A','B','C1','S1','X1','C2','S2','X2'}); 
    elseif NG == 1 % one peak
        fitType = fittype('A + B*x + C1*exp(-(1.6651*(x-X1)/S1)^2)',...
            'dependent',{'y'},'independent',{'x'},...
            'coefficients',{'A','B','C1','S1','X1'}); 
    end
    fo = fitoptions(fitType);          % variable with the fit options
    Ncff = numel(coeffnames(fitType)); % gets the number of coefficients
    N = size(Data,2)-1; % number of columns (that is, beads) in the (Y) data matrix 

    % Arrays of positions and their errors
    pkspos    = zeros(N,NG);
    pksposerr = zeros(N,NG);
    % Arrays of widths and their errors
    pkswdths    = zeros(N,NG);
    pkswdthserr = zeros(N,NG);
    % Arrays of distances and their errors
    Dz     = zeros(N,max(NG-1,1));
    Dz_err = zeros(N,max(NG-1,1));
    
    fitR2 = zeros(N,1); % "diagnostic" array with all the file's beads R^2 factors
    
    % smoothing filter kernel for helping the peak identification
    nfilt = 3;
    filt=(ones(1,nfilt))/nfilt;
    
    % gets the X data and scales it to its physical size
    Xdata = Data(:,1); 
    Xdata = phys*(Xdata-Xdata(1))/(Xdata(2)-Xdata(1)); 
    
    % loop for all the beads
    for ix = 1:N
        ignorecond = 0;       % =0 if the bead is valid (can be fitted with the desired amount of peaks), =1 otherwise
        Ydata = Data(:,1+ix); % define current profile
        
        % Get starting values for the background parameters
        Bkg_0 = Ydata(1);
        Bkg_1 = (Ydata(end)-Bkg_0)/(Xdata(end)-Xdata(1));
        
        D = 1; % this is the initial value for the width of the peaks
        
        % finds the peaks of the profile (convolution helps smoothing the noise))
        Ydata_c = conv(Ydata,filt,'same');
        [pks,locs,w,p] = findpeaks(Ydata_c,Xdata);
        [pks,isort] = sort(pks,'descend'); % sorts the maxima from higher to smaller
        p = p(isort);
        locs = locs(isort);
        w = w(isort);
        
        % if no peaks are found and NG = 1, takes the maximum
        % if the peaks number is not at least equal to NG, bead is ignored
        % otherwise we take the NG highest peaks
        if isempty(pks)&& NG == 1
            [PK,PI] = max(Ydata);
            LOCS = Xdata(PI);
            PROM = 0;
            W = 0;
        elseif numel(pks) < NG
            ignorecond = 1;
        else
            PK = pks(1:NG); 
            LOCS = locs(1:NG);
            PROM = p(1:NG);
            W = w(1:NG);
        end
        
        % if we ignore the bead, the parameters array can be filled with zeros
        % otherwise, we continue the analysis
        if ignorecond == 1
            fitR2(ix) = 0;  
            pkspos(ix,:) =      zeros(1,NG);
            pksposerr(ix,:) =   zeros(1,NG);
            pkswdths(ix,:) =    zeros(1,NG);
            pkswdthserr(ix,:) = zeros(1,NG);
            if NG > 1
                Dz(ix,:) =      zeros(1,NG-1);
                Dz_err(ix,:) =  zeros(1,NG-1);
            end
        else
            if NG == 3
                [LOCS,isort2] = sort(LOCS,'ascend'); % let's sort the peaks in depth order
                PKS = PK(isort2);
                WS = W(isort2);
                PROMS = PROM(isort2);
                P = 0.8*(PKS-Bkg_0); % we need to correct the peak height not to get overshoots
%                 fo.StartPoint = [Bkg_0 Bkg_1 P(1) WS(1) LOCS(1) P(2) WS(2) LOCS(2) P(3) WS(3) LOCS(3)];
                fo.StartPoint = [Bkg_0 Bkg_1 P(1) D 6 P(2) D 15 P(3) D 24];
            elseif NG == 2
                [LOCS,isort2] = sort(LOCS,'ascend'); % let's sort them in depth order
                PKS = PK(isort2);
                WS = W(isort2);
                PROMS = PROM(isort2);
                P = 0.5*(PKS-Bkg_0);
                fo.StartPoint = [Bkg_0 Bkg_1 P(1) D LOCS(1) P(2) D LOCS(2)];
            elseif NG == 1
                P = PK-Bkg_0;
                fo.StartPoint = [Bkg_0 Bkg_1 P D LOCS];
            end
            
            %% Part 6: bead data analysis
            [fitobj,gof] = fit(Xdata,Ydata,fitType,fo);
            % parameter values, confidence intervals and related error
            % (half the confidence interval size), R^2 values
            cval = coeffvalues(fitobj);
            cintval = confint(fitobj);
            cerrval = 0.5*diff(cintval,1);
            fitR2(ix) = gof.rsquare;
            
            if NG == 3
                % positions
                pkspos(ix,:) = cval(end-6:3:end);
                pksposerr(ix,:) = cerrval(end-6:3:end);
                % peak-to-peak (P2P) distances
                Dz(ix,:) = diff(pkspos(ix,:));
                Dz_err(ix,:) = sqrt(pksposerr(ix,1:2).^2 + pksposerr(ix,2:3).^2);
                % peak widths
                pkswdths(ix,:) = cval(end-7:3:end-1);
                pkswdthserr(ix,:) = cerrval(end-7:3:end-1);
            elseif NG == 2
                pkspos(ix,:) = cval(end-3:3:end);
                pksposerr(ix,:) = cerrval(end-3:3:end);
                Dz(ix) = diff(pkspos(ix,:));
                Dz_err(ix) = sqrt(sum(pksposerr(ix,:).^2));
                pkswdths(ix,:) = cval(end-4:3:end-1);
                pkswdthserr(ix,:) = cerrval(end-4:3:end-1);
            elseif NG == 1
                pkspos(ix) = cval(end);
                pksposerr(ix) = cerrval(end);
                Dz(ix,:) = 0;
                Dz_err(ix,:) = 0;
                pkswdths(ix,:) = cval(end-1);
                pkswdthserr(ix,:) = cerrval(end-1);
            end
        end
        
        % counter to track the software progress
        if mod(ix,10)==0
            disp(strcat('Analyzed',{' '},num2str(ix),' beads (',num2str(100*ix/N,'%.1f%%'),')'));
        end
    end
    
    %% Part7: data collection and plot
    % X label for the plots
    Xlab = 'Beads n.';
    
    %% Plot of the fitted widths of all the beads
    febw = figure('Name','beads fitted width (average intensity plots)');
        errorbar(pkswdths,pkswdthserr);
        xlabel(Xlab)
        ylabel('FWHM (\mum)')
        if NG == 2
            legend({'Top FWHM','Bottom FWHM'},'Location','northeast');
        elseif NG == 3
            legend({'Top FWHM','Central FWHM','Bottom FWHM'},'Location','northeast');
        end
        triplesave(febw,strcat(prefix,'bds_wdths',suffix));
    close(febw);
    
    %% Selection of the "best" beads
% ---------- This part can be customized to accomodate almost any      ----
% ---------- selection criteria. It all boils down to get the idx_sel  ----
% ---------- logical array and use the AND(idx_sel,<criteria>) as many ----
% ---------- times as needed                                           ----
    selstdfactor = 3;
    if NG == 1
        idx_sel = (fitR2>0.75);  % we select all the beads whose fit has R^2 higher than a certain amount
%         idx_sel = and(idx_sel,pkswdthserr<errfact*pkswdths);   % We ignore values with too big an uncertainty
        idx_sel = and(idx_sel,pkswdthserr<wdthmax); % We ignore the values that surpass reasonable levels
        idx_sel = and(idx_sel,pkswdths<wdthmax);
        idx_sel = and(idx_sel,abs(mean(pkswdths)-pkswdths) < selstdfactor*std(pkswdths));
    else
        idx_sel = (fitR2>0.75);  % we select all the beads whose fit has R^2 higher than a certain amount
%         idx_sel = ones(N,1);
        for ix = 1:NG
            spw = squeeze(pkswdths(:,ix));
%             idx_sel = and(idx_sel,(pkswdthserr(:,ix)<errfact*pkswdths(:,ix)));
            idx_sel = and(idx_sel,(squeeze(pkswdthserr(:,ix))<wdthmax));
            idx_sel = and(idx_sel,(spw < wdthmax));
            idx_sel = and(idx_sel,abs(spw-mean(spw)) < selstdfactor*std(spw));
        end
    end

    %% Plot of the fitted widths of the selected beads
    pkspos_sel    = pkspos(idx_sel,:);
    pksposerr_sel = pksposerr(idx_sel,:);
    pkswdths_sel    = pkswdths(idx_sel,:);
    pkswdthserr_sel = pkswdthserr(idx_sel,:);
    febwsel = figure('Name','"Selected" beads fitted width (average intensity plots)');
        errorbar(pkswdths_sel,pkswdthserr_sel);
        xlabel(Xlab)
        ylabel('FWHM (\mum)')
        if NG == 2
            legend({'Top FWHM','Bottom FWHM'},'Location','northeast');
        elseif NG == 3
            legend({'Top FWHM','Central FWHM','Bottom FWHM'},'Location','northeast');
        end
        triplesave(febwsel,strcat(prefix,'Sel_bds_wdths',suffix));
    close(febwsel);

    pkswdthavg     = mean(pkswdths_sel,1);
    pkswdthavg_err = avg_err(pkswdthserr_sel,1);
    pkswdthavg_std = std(pkswdths_sel,1);
    
    if NG == 1
        % creation of the values table
        save(strcat(prefix,'pks_main_values',suffix,'.mat'),...
            'pkspos','pksposerr',...
            'pkswdths','pkswdthserr',...         % info on positions and widths
            'pkspos_sel','pksposerr_sel',...
            'pkswdths_sel','pkswdthserr_sel',...
            'pkswdthavg','pkswdthavg_err','pkswdthavg_std');
        % filling of the GENERAL table of values
        gentable(ixT,:) = [pkswdthavg;pkswdthavg_err;pkswdthavg_std];
    else
        %% Plot of the peak distances in the Z-plots of the selected beads
        Dz_sel    = Dz(idx_sel,:);
        Dzerr_sel = Dz_err(idx_sel,:);
        fGd = figure('Name','"Selected" beads peaks distances'); 
            axebd = errorbar(Dz_sel,Dzerr_sel);
            xlabel(Xlab)
            ylabel('Peak-to-peak distance (\mum)')
            if NG == 3
                legend({'Top distance','Bottom distance'},'Location','northeast');
            end
            triplesave(fGd,strcat(prefix,'Sel_bds_dists',suffix));
        close(fGd);
        
        %% Plot of the ratio between peaks widths and peak distances in the Z-plots of the selected beads
        % Sums of the HWHM of adjacent peaks (only for selected beads)
        pksHwdthS_sel    = 0.5*(pkswdths_sel(:,1:end-1)+pkswdths_sel(:,2:end));
        pksHwdthSerr_sel = 0.5*sqrt(pkswdthserr_sel(:,1:end-1).^2 + pkswdthserr_sel(:,2:end).^2);
        % Ratio between HWHM(i)+HWHM(i+1) and d(i,i+1)
        pkswdr_sel    = pksHwdthS_sel./Dz_sel;
        pkswdrerr_sel = pkswdr_sel.*sqrt((pksHwdthSerr_sel./pksHwdthS_sel).^2 + ...
            (Dzerr_sel./Dz_sel).^2);
        fwdr = figure('Name','Width-to-distance ratio'); 
            axebwdr = errorbar(pkswdr_sel,pkswdrerr_sel);
            xlabel(Xlab)
            ylabel('FWHM/distance')
            if NG == 3
                legend({'Top tail','Bottom tail'},'Location','northeast');
            end
            triplesave(fwdr,strcat(prefix,'Sel_bds_wdr',suffix));
        close(fwdr);
        
        Dzavg      = mean(Dz_sel,1);
        Dzavg_err  = avg_err(Dzerr_sel,1);
        Dzavg_std  = std(Dz_sel,1);
        pkswdravg      = mean(pkswdr_sel,1);
        pkswdravg_err  = avg_err(pkswdrerr_sel,1);
        pkswdravg_std  = std(pkswdr_sel,1);
        
        % creation of the values table
        save(strcat(prefix,'pks_main_values',suffix,'.mat'),'pkspos','pksposerr',...
            'pkswdths','pkswdthserr',...      % info on positions and widths...
            'pkswdths_sel','pkswdthserr_sel',...
            'pkspos_sel','pksposerr_sel',...
            'pkswdthavg','pkswdthavg_err','pkswdthavg_std',...
            'Dz_sel','Dzerr_sel',...          % on distances...
            'Dzavg','Dzavg_err','Dzavg_std',...
            'pkswdr_sel','pkswdrerr_sel',...  % ...and their ratio
            'pkswdravg','pkswdravg_err','pkswdravg_std');
        
        % filling of the GENERAL table of values
        if NG == 2
            gentable(ixT,:) = cat(2,...
            pkswdthavg(1),pkswdthavg_err(1),pkswdthavg_std(1),...
            pkswdthavg(2),pkswdthavg_err(2),pkswdthavg_std(2),...
            Dzavg     ,Dzavg_err     ,Dzavg_std     ,...
            pkswdravg ,pkswdravg_err ,pkswdravg_std);
        else
            gentable(ixT,:) = cat(2,...
            pkswdthavg(1),pkswdthavg_err(1),pkswdthavg_std(1),...
            pkswdthavg(2),pkswdthavg_err(2),pkswdthavg_std(2),...
            pkswdthavg(3),pkswdthavg_err(3),pkswdthavg_std(3),...
            Dzavg(1)     ,Dzavg_err(1)     ,Dzavg_std(1)     ,...
            Dzavg(2)     ,Dzavg_err(2)     ,Dzavg_std(2)     ,...
            pkswdravg(1) ,pkswdravg_err(1) ,pkswdravg_std(1) ,...
            pkswdravg(2) ,pkswdravg_err(2) ,pkswdravg_std(2) );
        end
    end
    
    % counter to track the software progress
    disp(strcat('Analysis on',{' '},foldname,{' '},'completed. (',num2str(100*ixT/Ncsv,'%.1f%%'),')'));
end

% translates the table array to an actual table
if NG == 1
    T_gentbl = array2table(gentable,'VariableNames',...
        {'FWHM_um','FWHM_err_um','FWHM_std_um'});
elseif NG == 2
    T_gentbl = array2table(gentable,'VariableNames',...
        {'top_FWHM_um','top_FWHM_err_um','top_FWHM_std_um',...
        'bottom_FWHM_um','bottom_FWHM_err_um','bottom_FWHM_std_um',...
        'P2P_distance_um','P2P_distance_err_um','P2P_distance_std_um',...
        'FWHM_dist_ratio','FWHM_dist_ratio_err','FWHM_dist_ratio_std'});
else
    T_gentbl = array2table(gentable,'VariableNames',...
        {'top_FWHM_um','top_FWHM_err_um','top_FWHM_std_um',...
        'central_FWHM_um','central_FWHM_err_um','central_FWHM_std_um',...
        'bottom_FWHM_um','bottom_FWHM_err_um','bottom_FWHM_std_um',...
        'top_P2P_distance_um','top_P2P_distance_err_um','top_P2P_distance_std_um',...
        'bottom_P2P_distance_um','bottom_P2P_distance_err_um','bottom_P2P_distance_std_um',...
        'top_FWHM_dist_ratio','top_FWHM_dist_ratio_err','top_FWHM_dist_ratio_std',...
        'bottom_FWHM_dist_ratio','bottom_FWHM_dist_ratio_err','bottom_FWHM_dist_ratio_std'});
end

% writes the table to a csv file
writetable(T_gentbl,strcat(sv_pth_gen,'pks_main_values',suffix,'.csv'));

% additional plots required when more than a single .csv file is found
if Ncsv == 1
    return;
else
    %% Plot of the fitted widths of all the beads
    X_range = 1:Ncsv;
    fgw = figure('Name','Beads average FWHM vs depth');
        if NG == 1
            errorbar(X_range,gentable(:,1),gentable(:,3));
        elseif NG == 2
            errorbar(X_range,gentable(:,1:3:4),gentable(:,3:3:6));
        elseif NG == 3
            errorbar(X_range,gentable(:,1:3:7),gentable(:,3:3:9));
        end
        xlabel("Cam position n.") % xlabel("Z (um)")
        ylabel('FWHM (\mum)')
        if NG == 2
            legend({'Top FWHM','Bottom FWHM'},'Location','northeast');
        elseif NG == 3
            legend({'Top FWHM','Central FWHM','Bottom FWHM'},'Location','northeast');
        end
        triplesave(fgw,strcat(sv_pth_gen,'gen_bds_avg_wdths',suffix));
    close(fgw);

    if NG > 1
        %% Plot of the width_to_distance ratio of all the beads
        fgwdr = figure('Name','Beads average width-to-distance ratio');
            if NG == 2
                errorbar(X_range,gentable(:,10),gentable(:,12));
            elseif NG == 3
                errorbar(X_range,gentable(:,10:3:13),gentable(:,12:3:15));
            end
            xlabel("Cam position n.") % xlabel("Z (um)")
            ylabel('FWHM/distance')
                if NG == 3
                    legend({'Top tail','Bottom tail'},'Location','northeast');
                end
            triplesave(fgwdr,strcat(sv_pth_gen,'gen_bds_wdr',suffix));
        close(fgwdr);
    end
end

%% ---------- FUNCTIONS ---------------------------------------------------
function errarrout = avg_err(errarrin,dim)
if ~isnumeric(errarrin)
    error('avg_err accepts only numeric arrays');
elseif ~any(size(errarrin))
    error('nonpositive array dimension found in avg_err array');
elseif prod(dim)==1
    warning('avg_err used on single value array. Verify the code does not add excessive complexity');
end

errarrout = sqrt(sum(errarrin.^2,dim))/size(errarrin,dim);
end

function triplesave(handle,filename)
if isgraphics(handle,'Axis')
    saveas(handle.Parent,strcat(filename,'.fig'));
    saveas(handle.Parent,strcat(filename,'.png'));
    saveas(handle.Parent,strcat(filename,'.jpg'));
elseif isgraphics(handle,'Figure')
    saveas(handle,strcat(filename,'.fig'));
    saveas(handle,strcat(filename,'.png'));
    saveas(handle,strcat(filename,'.jpg'));
else
    error('handle argument must be either a ''Axis'' or a ''Figure'' object');
end 
end