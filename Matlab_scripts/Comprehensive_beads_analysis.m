clear all
% close all
clc;

%% constants used during the script run
ex_name = 'Good_beads_exc.xlsx';
phys_x = 0.162;   % the optical sizes of the voxel
phys_y = 0.162;
phys_z = 0.2;
dil = floor(2);   % "dilatation" factor to enlarge the substack width
numcond = 0;      % in questionloop software uses 1/0/2 (if=1) or y/n/e (if=0) options
fitcond = 1;      % the fit involves the central values (if=0) or the peak ones (if=1)
beadplotcond = 0; % chooses whether to plot (if=1) or not (if=0) the substack plots
condoff = 0;      % chooses whether to change coordinates from PSFj to 
                  % Imagej/Matlab (if=1) or not (if=0)

%% Preliminary folder investigation
dirlist = ls;               % list of files inside the current folder
[Nlist, Nchar] = size(dirlist);
dirlist_tif = ls('*.tif');  % list of .tif files inside the current folder
dirlist_fold_arr = zeros(Nlist,1); % list of folders
Ndir = 0;
for ixd = 1:Nlist
    if isfolder(dirlist(ixd,:))
        dirlist_fold_arr(ixd) = ixd; % list of folders
        Ndir = Ndir + 1;
    end
end
dirlist_fold_arr = dirlist_fold_arr(dirlist_fold_arr~=0);
dirlist_fold = cell(Ndir,1);
for ixdir = 1:Ndir
    dirlist_fold{ixdir} = dirlist(dirlist_fold_arr(ixdir),:); 
end

% asks whether to analyse the files inside the folders (OR NOT)
if (size(dirlist_fold,1))~=0
    if numcond == 1
        foldcond = questionloop(input('Analyse the files inside the subfolders? 1(Yes)/0(No) [press 2 to exit]'),...
            'Invalid character. Use 1 or 0 to continue, 2 to exit.',[1,0,2]);
    else
        foldcondc = questionloopchar(input('Analyse the files inside the subfolders? y/n [press "e" to exit]','s'),...
            'Invalid character. Use "y" or "n" to continue, "e" to exit.',{'y','n','e'});
        foldcond = convertcharnum(foldcondc,{'y','n','e'},[1,0,2]);
    end
end


%% Main "loop"
if foldcond == 2        % If asked to exit, exits
    disp('Closing program.');
    return;
elseif foldcond == 0    % We are analyzing only the files in the current folder
    if size(dirlist_tif,1)==0
        disp('No file found in directory. Closing program');
        return
    elseif size(dirlist_tif,1)> 1
        if numcond == 1
            multifilecond = questionloop(input('Analyse all the files inside the current folder? 1(Yes)/0(No) [press 2 to exit]'),...
                'Invalid character. Use 1 or 0 to continue, 2 to exit.',[1,0,2]);
        else
            multifilecondc = questionloopchar(input('Analyse all the files inside the current folder? y/n [press "e" to exit]','s'),...
                'Invalid character. Use "y" or "n" to continue, "e" to exit.',{'y','n','e'});
            multifilecond = convertcharnum(multifilecondc,{'y','n','e'},[1,0,2]);
        end
        % if asked to exit, exits
        if multifilecond == 2 % strcmpi(multifilecond,'e')
            disp('Closing program.');
            return
        elseif multifilecond == 1 % multifilecond == 0 % strcmpi(multifilecond,'n') % Analyzes all the .tiff files in the directory
            disp('Sorry. Software functionality still in development. Closing program.');
            return
        else % Analyzes only one .tiff file in the directory
            disp('Sorry. Software functionality still in development. Closing program.');
            return;
            tiffname = uigetfile('*.tif','Select the file...'); % Selects the desired file
            Im = getmultitiff(tiffname); % We need to open a multitiff
        end
    else % Analyzes the only .tiff file in the directory
        % We want to fit everything with a linear background value and
        % three gaussians (here the S# are exactly the FWHMs thanks to the 
        % 1.6651 factor). Use 0.8326 if you want the HWHM.
        fitType = fittype('A + B*x + C1*exp(-(1.6651*(x-X1)/S1)^2) + C2*exp(-(1.6651*(x-X2)/S2)^2) + C3*exp(-(1.6651*(x-X3)/S3)^2)',...
            'dependent',{'y'},'independent',{'x'},...
            'coefficients',{'A','B','C1','S1','X1','C2','S2','X2','C3','S3','X3'});
        Ncff = numel(coeffnames(fitType)); % number of coefficients
        fo = fitoptions(fitType);
        % We get the image and a table with ROIs central coordinates and widths
        [Im,row,col,frm] = getmultitiff(dirlist_tif(1,:));
        TIm = readtable(ex_name);
        NamesIm = table2cell(TIm(1,:));
        DataIm = table2array(TIm(2:end,:));
        Nc = size(DataIm,1);
        
        %careful! the order is x,y,z here
        phys_arr = [phys_x,phys_y,phys_z];
        % gets the coordinates in pixel units (with 0 as the left/top sides)
        coordlist = floor(DataIm(:,1:3)./repmat(phys_arr,Nc,1));
        if condoff == 1
            offset = floor(([col,row,0])/2)+1;
            coordlist(:,1)= offset(1)+coordlist(:,1);
            coordlist(:,2)= offset(2)-coordlist(:,2);
        end
        % gets the widths in pixel units, and multiplies them 
        % to get the half-range of the stacks 
        widthlist = dil*floor(DataIm(:,4:6)./repmat(phys_arr,Nc,1));
        
        % Preallocates memory for the (cell) arrays where data will be stored
        %substacks = cell(Nc,1);
        max_Zlist = cell(Nc,1);
        min_Zlist = cell(Nc,1);
        avg_Zlist = cell(Nc,1);
        contr_Zlist = cell(Nc,1);
        Z_plot = cell(Nc,1);
        cval(Nc,Ncff) = 0;
        cerrval(Nc,Ncff) = 0;
        fitR2(Nc) = 0;
        
        if beadplotcond == 1
            % creates the figures and axes where the plots will be shown
            fmax = figure('Name','Peak Intensity');
            axmax = axes(fmax);     hold(axmax,'on');
            favg = figure('Name','Average intensity');
            axavg = axes(favg);     hold(axavg,'on');
            fcontr = figure('Name','Contrast');
            axcontr = axes(fcontr); hold(axcontr,'on');
            fprof = figure('Name','Centre Z-profile');
            axprof = axes(fprof);   hold(axprof,'on');
        end
        
        fitcond = questionloop(fitcond,...
            'Invalid fit condition. Type 0 for fit of central data, 1 for fit of peak data, 2 to exit',...
            [1,0,2]);
        
        if fitcond == 2
            disp('Closing program.');
            return
        end
        
        for ixb = 1:Nc
            % We create the ranges over which we'll span in the three
            % dimensions
            xrng = (max(coordlist(ixb,1)-widthlist(ixb,1),1):...
                min(coordlist(ixb,1)+widthlist(ixb,1),col))';
            yrng = (max(coordlist(ixb,2)-widthlist(ixb,2),1):...
                min(coordlist(ixb,2)+widthlist(ixb,2),row))';
            zrng = (max(coordlist(ixb,3)-widthlist(ixb,3),1):...
                min(coordlist(ixb,3)+widthlist(ixb,3),frm))';
            Czrng = floor(zrng - mean([zrng(1) zrng(end)])); % centred z-range
            
            % We create a substack over these ranges
            %substacks{ixb} = Im(xrng,yrng,zrng);
            substack = Im(yrng,xrng,zrng); %remember? Matlab orders matrices as y,x,z
            
            %max_Zlist = max2d(substacks{ixb});
            max_Zlist{ixb}  = max2d(substack,1,2);  % peak intensity over Z(should have one peak)
            avg_Zlist{ixb}  = avg2d(substack,1,2);  % total energy over Z(should be constant)
            min_Zlist{ixb}  = min2d(substack,1,2);  % min intensity over Z
            contr_Zlist{ixb}= contr(max_Zlist{ixb},min_Zlist{ixb}); % contrast over Z
            % Zplot over the assumed position of the bead
            Z_plot{ixb}     = squeeze(substack(widthlist(ixb,2)+1,widthlist(ixb,1)+1,:));
            if beadplotcond == 1    
                plot(axmax,Czrng,max_Zlist{ixb}); 
                plot(axavg,Czrng,avg_Zlist{ixb});
                plot(axcontr,Czrng,contr_Zlist{ixb});
                plot(axprof,Czrng,Z_plot{ixb});
            end
            
            % Now, we prepare the initial parameters to fit the Z-plot
            if fitcond == 0
                fitYdata = Z_plot{ixb};
            elseif fitcond == 1
                fitYdata = max_Zlist{ixb};
            end
            % we initially assume a constant bkg
            Bkg_0 = min(fitYdata);
            Bkg_1 = 0;
            % we assume similarly scaled tails, with 2/3rds of the maximum height 
            C1 = 0.66*max(fitYdata)-Bkg_0; C2 = C1; C3 = C1;
            % we assume 3 gaussians with FWHM as large as 1/8th of the substack range
            S1 = 0.25*max(Czrng); S2 = S1; S3 = S1; 
            % we assume equispaced gaussians, one at the centre, separated
            % by 1 FWHM (approximation a bit on the eccess)
            X2 = 0; X1 = X2-S1; X3 = X2+S1;
            
            % set starting fit parameters
            fo.StartPoint = [Bkg_0 Bkg_1 C1 S1 X1 C2 S2 X2 C3 S3 X3];
            % fit the data
            [fitobj,gof] = fit(Czrng,fitYdata,fitType,fo);
            
            % get the parameter values and errors, and the R^2 value
            cval(ixb,:) = coeffvalues(fitobj);
            cerrval(ixb,:) = 0.5*diff(confint(fitobj),1);
            fitR2(ixb) = gof.rsquare;
                        
            % small counter to show the progress of the work
            if mod(ixb,10)==0
                disp(strcat('Analyzed',{' '},num2str(ixb),' beads (',num2str(100*ixb/Nc,'%.1f%%'),')'));
            end
        end
        if beadplotcond == 1
            % Stops holding the figures and saves the figures (derived from
            % <axes_name>.Parent) with their onw titles names (derived from 
            % <axes_name>.Parent.Name)
            hold(axmax,'off');   triplesave(axmax.Parent,axmax.Parent.Name);
            hold(axavg,'off');   triplesave(axavg.Parent,axavg.Parent.Name);
            hold(axcontr,'off'); triplesave(axcontr.Parent,axcontr.Parent.Name);
            hold(axprof,'off');  triplesave(axprof.Parent,axprof.Parent.Name);
        end
        
        % get the peaks positions with errors
        peakspos    = cval(:,end-6:3:end);
        peaksposerr = cerrval(:,end-6:3:end);
        peakswdth    = cval(:,end-7:3:end-1);
        peakswdtherr = cerrval(:,end-7:3:end-1);
        
        %indices where the position error in ALL the three gaussian is < 10 frames
        thrE = 2; % median(peaksposerr(:));
        sel_idx = and(and(peaksposerr(:,1)<thrE,peaksposerr(:,2)<thrE),peaksposerr(:,3)<thrE);
        peaksposerr_sel = peaksposerr(sel_idx,:);
        peakspos_sel    = peakspos(sel_idx,:);
        peakswdth_sel    = peakswdth(sel_idx,:);
        peakswdtherr_sel = peakswdtherr(sel_idx,:);
        
        distpks     = diff(peakspos_sel,1,2);
        distpkserr  = sqrt( peaksposerr_sel(:,1:2).^2 + peaksposerr_sel(:,2:3).^2 );
        
        avgdpks     = mean(distpks,1);
        avgdpks_err = avg_err(distpkserr,1);
        avgdpks_std = std(distpks,1);
        avgwpks     = mean(peakswdth_sel,1);
        avgwpks_err = avg_err(peakswdtherr_sel,1);
        avgwpks_std = std(peakswdth_sel,1);
        
        % Here we take the ratio between the widths and the distances of
        % the peaks. To make sense of it, each distance divides not a
        % proper width, but the sum of the HWHMs of the two correspondent peaks
        pkswdr     = (peakswdth_sel(:,1:2) + peakswdth_sel(:,2:3))./(2*distpks);
        pkswdr_err = sqrt((peakswdtherr_sel(:,1:2)./(2*distpks)).^2 +... 
            (peakswdtherr_sel(:,2:3)./(2*distpks)).^2 +...
            (pkswdr.*distpkserr./distpks).^2);
        pkswdr_std = std(pkswdr,1);
        
        Xlab = 'Bead n.';
        fGw = figure('Name','Gaussians widths'); 
        axebw = errorbar(peakswdth_sel*phys_z,peakswdtherr_sel*phys_z);
        xlabel(Xlab)
        ylabel('FWHM (\mum)')
        legend('Bottom peak','Central peak','Top peak','Location','SouthEast');
        triplesave(fGw,fGw.Name);
        fGd = figure('Name','Gaussians distances'); 
        axebd = errorbar(distpks*phys_z,distpkserr*phys_z);
        xlabel(Xlab)
        ylabel('Peak-to-peak distance (\mum)')
        legend('Bottom distance','Top distance','Location','SouthEast');
        triplesave(fGd,fGd.Name);
        fwdr = figure('Name','Width-to-distance ratio'); 
        axebwdr = errorbar(pkswdr,pkswdr_err);
        xlabel(Xlab)
        ylabel('FWHM/distance')
        legend('Bottom half','Top half','Location','southeast');
        triplesave(fwdr,fwdr.Name);
    end
else % We are analyzing the files inside the subfolders
    disp('Sorry. Software functionality still in development. Closing program.');
    return;
end

%% ------------ Functions -------------------------------------------------
function answerout = questionloop(answerin,msg,answerarr)
if ~isnumeric(answerarr)
    error('Error! answerarr variable not a numeric array!');
end
loop = 0;
while loop==0
    % 1. determines matches
    boolarr = (answerin==answerarr);
    % 2. returns the answer if matches are found, otherwise repeats loop
    if any(boolarr)
        loop = 1;
    else
        answerin = input(strcat(msg,{''}));
    end
end
answerout = answerin;
end

function answercharout = questionloopchar(answercharin,msg,answercell)
if ~iscell(answercell)
    error('Error! answercell variable not a cell array!');
end
loop = 0;
while loop==0
    % 1. determines matches
    boolarr = strcmpi(answercharin,answercell);
    % 2. returns the answer if matches are found, otherwise repeats loop
    if any(boolarr)
        loop = 1;
    else
        answercharin = input(strcat(msg,{''}),'s');
    end
end
answercharout = answercharin;
end

function numout = convertcharnum(charin,chararr,numarr)
if ~iscell(chararr)
    error('Error! chararr variable not a cell array!');
elseif ~isnumeric(numarr)
    error('Error! numarr variable not a numeric array!');
elseif numel(chararr)~=numel(numarr)
    error('Error! Choices arrays need to have the same sizes');
end

numout = numarr(strcmpi(charin,chararr));
end

function max_arr = max2d(stack, dim1, dim2) % max
if numel(size(stack))<3
    error('The array has not enough dimensions');
elseif dim1<1||dim2<1
    error('Dimensions need to be positive integer values');
else
    max_arr = squeeze(max(max(stack,[],dim2),[],dim1));
end
end

function min_arr  = min2d(stack, dim1, dim2)  % min
if numel(size(stack))<3
    error('The array has not enough dimensions');
elseif dim1<1||dim2<1
    error('Dimensions need to be positive integer values');
else
    min_arr = squeeze(min(min(stack,[],dim2),[],dim1));
end
end

function avg_arr = avg2d(stack, dim1, dim2)  % average
if numel(size(stack))<3
    error('The array has not enough dimensions');
elseif dim1<1||dim2<1
    error('Dimensions need to be positive integer values');
else
    Navg = (size(stack,dim1)*size(stack,dim2));
    avg_arr = squeeze(sum(sum(stack,dim2),dim1)/Navg);
end
end

function contr_arr = contr(max_arr,min_arr) % contrast
contr_arr = (max_arr - min_arr)./(max_arr+min_arr);
end

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