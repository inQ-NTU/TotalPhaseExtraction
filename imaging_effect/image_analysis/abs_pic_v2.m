classdef abs_pic_v2 < dynamicprops
    %-----------------------------------------------------------------------------
    % File name  : abs_pic.m
    % Created on : 16-Mar-2010
    % Description: object to store an absorption pic and calculate its basic
    %              properties
    %----------------------------------------------------------------------
    % Updated with the following features: (should still be backwards compatible)
    %  - Now supports overlapping ROI and RBC 
    %  - atompic and backpic can be set manually by calling the set_pics()
    %    function. To use manually set pictures, set the variable useLoadedPics
    %    to 0 (zero).
    
    properties (SetAccess = public)
      pic_path = '';
      pic_number = 1;
      
      ROI = []; %coordinates for ROI
      RBC = []; %coordinates for RBC
      
      ROI_log = []; % logical matrix for ROI
      RBC_log = []; % logical matrix for RBC
      
      background_correction = 1; %do background corr. on ROI_data? 0 = false, 1 = true
      gain = 3;

      atompic=[];
      backpic=[];

      useBackPics   = 2; % set to 0 for no back pics, 1 for back pics or to 2 to substract the usual average counts
      useLoadedPics = 1; % if 1, load pics as usual. If 0, set atompic and backpic manually
            
      abs_data = [];
      ROI_data = [];
      RBC_data = [];
      
      shotnoise_data = [];
      
      TOF = 20; %needed for fringe analysis
      
      %imaging system properties. set automatically by choosing the actual
      %imaging system. see below!
      
      imaging_system = 'longitudinalAndor' % or transversalAndor or Pixelfly
      
      pixel2um=0;
      magnification=0;
      pixelsize=0;  %cm
      abs_cross=0;  %cm^2
      backgroundCounts=0;
      alpha=0;
      ff=0;
      epc=0;
      qe=0;
      texp=0;
      kx=0; % standing wave [um^-1] 
      phi=0; % standing wave offset
      tilt=0; % probe beam tilt relative to chip
      inSitu=0; % flag for toggling in-Situ absorption formula
      
      scalefact = 0;
      atomval=0;
      backval=0;
    end
    
    methods 
    
    % constructor
         function obj = abs_pic_v2(props)
                
            % Create a chip_trap object
            obj=obj@dynamicprops;

            % define object properties
            if nargin > 0
                fn = fieldnames(props);
                for l = 1:length(fn)                      
                    if isprop(obj, fn{l}) % only copy field if defined among properties
                        eval(['obj.',fn{l},' = props.',fn{l},';']);
                    end
                end
            end
                
            % format picture path according to OS
            if (~exist(obj.pic_path))             
                %if pic path not existent make path mac/pc compatible
                macroot = '/Volumes/KRbData'; 
                pcroot = 'G:';
                linuxroot = '/mnt/krbdata/';
                if (exist('G:\data_analysis'))
                    pcroot = 'G:'; 
                elseif (exist('X:\data_analysis'))
                    pcroot = 'X:';
                end

                if filesep == '/'   %mac
                    %replace slashes
                    obj.pic_path = strrep(obj.pic_path,'\',filesep);
                    %replace root folder
                    if ismac
                        if length(obj.pic_path)>=length(macroot)
                            if ~strcmp(obj.pic_path(1:length(macroot)),macroot)
                                obj.pic_path = strrep(obj.pic_path,obj.pic_path(1:length(pcroot)),macroot);
                            end
                        end
                    elseif isunix
                        if length(obj.pic_path)>=length(linuxroot)
                            if ~strcmp(obj.pic_path(1:length(linuxroot)),linuxroot)
                                obj.pic_path = strrep(obj.pic_path,obj.pic_path(1:length(pcroot)),linuxroot);
                            end
                        end
                    end
                elseif filesep == '\' %PC
                    %replace slashes
                    obj.pic_path = strrep(obj.pic_path,'/',filesep);
                    %replace root folder
                    if length(obj.pic_path)>=length(pcroot)
                        if ~strcmp(obj.pic_path(1:length(pcroot)),pcroot)
                            obj.pic_path = strrep(obj.pic_path,obj.pic_path(1:length(macroot)),pcroot);
                        end
                    end
                end
            end
            %add path to fitting functions
            %addpath(['.' filesep '@abs_pic' filesep 'fitfunctions'])

            % set camera properties
            if strcmpi(obj.imaging_system,'longitudinalAndor') || strcmpi(obj.imaging_system,'LAndor')
                obj.magnification=5.3;
                obj.pixelsize=13*10^-4;  %cm
                obj.pixel2um=obj.pixelsize*1e4/obj.magnification;
                obj.abs_cross=2.907*1e-9;  %cm^2
                obj.backgroundCounts=398;
                obj.alpha=1;
                obj.ff=1;
                obj.qe=0.2; %note: measured value for quantum-eff/epc
                obj.epc=1; %electrons per count, note: not really 1, value for epc contained in obj.qe
                obj.texp=75e-6;
            elseif strcmpi(obj.imaging_system,'verticalAndor')|| strcmpi(obj.imaging_system,'VAndor')
                obj.magnification=8.22;         
                obj.pixelsize=16*10^-4; %cm     
                obj.pixel2um=obj.pixelsize*1e4/obj.magnification;
                obj.abs_cross=2.907*1e-9; %cm^2´trans 
                obj.backgroundCounts=1001;
                obj.alpha=1;
                obj.ff=1;
                obj.epc=1.84; %electrons per count
                obj.qe=0.8;   %quantum efficiency
                obj.texp=50e-6;
            elseif strcmpi(obj.imaging_system,'transversalAndor')|| strcmpi(obj.imaging_system,'TAndor')
                obj.magnification=12.39;         
                obj.pixelsize=13*10^-4; %cm     
                obj.pixel2um=obj.pixelsize*1e4/obj.magnification;
                obj.abs_cross=2.907*1e-9; %cm^2´trans 
                obj.backgroundCounts=458;
                obj.alpha=1/0.54;
                obj.ff=1;
                obj.qe=0.57; %note: measured value for quantum-eff/epc
                obj.epc=2.87; %electrons per count, note: not really 1, value for epc contained in obj.qe
                obj.texp=75e-6;
            elseif strcmpi(obj.imaging_system,'Pixelfly')
                obj.magnification=1.68;         
                obj.pixelsize=6.45*10^-4; %cm
                obj.pixel2um=obj.pixelsize*1e4/obj.magnification;
                obj.abs_cross=2.907*1e-9; %cm^2´trans
                obj.backgroundCounts=0;
                obj.useBackPics=2;
                obj.alpha=1/0.54;
                obj.ff=0;
            end

         end % end constructor
        
         
    % all other functions
         %----------------------------------------------------------------%
         function do_background_correction(obj,value)
             obj.background_correction = value;
         end
         
         
         %----------------------------------------------------------------%
         function set_gain(obj,value)
             obj.gain = value;
             %obj.calc_ROI_data();
         end
         
         %----------------------------------------------------------------%
        function setLogicalRegions(obj)
             % determine picture size
            if ~isempty(obj.atompic)
                picsize = size(obj.atompic);
            elseif ~isempty(obj.abs_data)
                picsize = size(obj.abs_data);
            else
                disp('Could not determine size of picture!')
                return
            end
            
            function log_mat = createLogical(matsize, region)
                log_mat = false( matsize );
                for i = 1:size(region,3)
                    reg_slice = squeeze(region(:,:,i));
                    log_mat(reg_slice(1):reg_slice(2) , reg_slice(3):reg_slice(4)) = true;
                end
            end
            
  
            if isempty(obj.RBC)
                % if no ROI specified, set whole logical region to true,
                % i.e use full picture as ROI
                obj.ROI_log = true( picsize );
            else
                obj.ROI_log = createLogical(picsize, obj.ROI);
            end
            
            if isempty(obj.RBC)
                % if no RBC specified, set whole logical region to false,
                % i.e. use no RBC
                obj.RBC_log = false( picsize );
            else
                obj.RBC_log = createLogical(picsize, obj.RBC);
            end
            
            obj.RBC_log( obj.RBC_log & obj.ROI_log ) = false; % subtract ROI from RBC
        end
         
         %----------------------------------------------------------------%
         function setROI(obj,ROI_new)            
            % update object properties
            obj.ROI = ROI_new;
            obj.setLogicalRegions();
         end
         
         %----------------------------------------------------------------%
         function setRBC(obj,RBC_new)            
            % update object properties
            obj.RBC = RBC_new;
            obj.setLogicalRegions();
         end
         
         %----------------------------------------------------------------%

         function load_pics(obj,pics2load)
            % loads the pictures from the disc
            % pics2load is optinal and defines the pictures which are loaded. 
            % 0=atom and back pic are loaded
            % 1=only atom pic is loaded
            % 2=only back pic is loaded
             
            if not(exist('pics2load','var'))
                pics2load=0;
            end
            

            try
                 
                if pics2load == 1
                    obj.atompic = double(imread([obj.pic_path ,filesep,num2str(obj.pic_number),'-atomcloud','.tif']));
                elseif pics2load == 2
                    obj.backpic = double(imread([obj.pic_path ,filesep,num2str(obj.pic_number),'-withoutatoms','.tif']));
                else
                    obj.atompic = double(imread([obj.pic_path ,filesep,num2str(obj.pic_number),'-atomcloud','.tif']));
                    obj.backpic = double(imread([obj.pic_path ,filesep,num2str(obj.pic_number),'-withoutatoms','.tif']));
                end  
                
                

            catch
                %rethrow(lasterror);
                disp(['failed to load pic ' ,obj.pic_path ,filesep,num2str(obj.pic_number),'-atomcloud','.tif'])
            end
             
             % update ROI and RBC logicals
             obj.setLogicalRegions();
         end
         
         
         function set_pics(obj, atompic_in, backpic_in)
             % set pictures directly
             
             if ~isempty(atompic_in)
                 obj.atompic = atompic_in;
             else
                 disp('Warning: No atompic set!')
             end
             
             if ~isempty(backpic_in)
                 obj.backpic = backpic_in;
             else
                 disp('Warning: No backpic set!')
             end
             
             % update ROI and RBC logicals
             obj.setLogicalRegions();
         end
         
         
         %----------------------------------------------------------------%
         function output = load_data(obj)
             
             %define some constants 
             h = 6.6261e-34; %Placks constant
             c = 2.99e8; %speed of light
             
             if obj.ff == 1 
                 I_sat = 16.6933; %Rb saturation intensity [W/m2]
             else
                 I_sat = inf; % not taking saturation into account 
             end

             
             if obj.useLoadedPics
                 % try to load pic
                 obj.load_pics();
             else
                 % assume that pictures have been set manually
                 if isempty(obj.backpic)
                     error('No backpic set! Call function set_pics() to directly set pictures.')
                 end
                 if isempty(obj.atompic)
                     error('No atompic set! Call function set_pics() to directly set pictures.')
                 end
             end
             
            % adjust for background if desired
            if obj.useBackPics == 0 %use no backpics
                atom = obj.atompic;
                back = obj.backpic;
            elseif obj.useBackPics == 1 %use backpics from file                
                error('Use backpics from file (useBackPics==1) currently not supported!')
            else %use a fixed number as background count                    
                atomback = obj.backgroundCounts;
                atom = obj.atompic-atomback;
                atom = atom.*(atom>=0);

                backback = obj.backgroundCounts;
                back = obj.backpic-backback;
                back = back.*(back>=0);
             end  
             

             if max(max(atom)) == -1000
                output = 0;
                return
             end
             
             try

                Natom = atom*(obj.epc);
                Nback = back*(obj.epc);
                 
                % convert to intensities
                atom = atom.*(h*c/(780e-9*obj.pixel2um^2*1e-12*obj.texp))*(obj.epc/obj.qe);
                back = back.*(h*c/(780e-9*obj.pixel2um^2*1e-12*obj.texp))*(obj.epc/obj.qe);

                if obj.background_correction == 1 %if wanted do background correction
                    %the old rescaler function, here used locally
                    if isempty(obj.RBC)
                        scalefactor = 1;
                    else
                        atom_value = sum(atom(obj.RBC_log));
                        back_value = sum(back(obj.RBC_log));
                        scalefactor = atom_value./back_value;

                        obj.atomval = atom_value;
                        obj.backval = back_value;
                        obj.scalefact = scalefactor;
                    end

                    back = back*scalefactor;
                    Nback = Nback*scalefactor;
                end
                 
                % calculate absorption data
                obj.abs_data   = -obj.alpha*log(atom./(back+eps)) + (back - atom)/I_sat;  % full formula
                obj.shotnoise_data  = (obj.alpha + back./I_sat).^2 .* ( 1./(Natom+eps) + 1./(Nback+eps) );
                
                obj.calc_ROI_data();
                output = 1;
             catch
                rethrow(lasterror);
                disp(['failed to calculate ' ,obj.pic_path ,filesep,num2str(obj.pic_number),'-atomcloud','.tif'])
                output=0;
             end                 
         end
         
         %----------------------------------------------------------------%
         function calc_ROI_data(obj)
             % apply logical filter
             tmp_data = obj.abs_data;
             tmp_data( ~obj.ROI_log ) = 0;
             
             % remore 0-padding
             tmp_data( ~any(tmp_data,2), : ) = [];  %rows
             tmp_data( :, ~any(tmp_data,1) ) = [];  %columns
             
             obj.ROI_data = tmp_data;
         end
         
         %----------------------------------------------------------------%
         function calc_RBC_data(obj)
             % apply logical filter
             tmp_data = obj.abs_data;
             tmp_data( ~obj.RBC_log ) = 0;
             
             % remore 0-padding
             tmp_data( ~any(tmp_data,2), : ) = [];  %rows
             tmp_data( :, ~any(tmp_data,1) ) = [];  %columns
             
             obj.RBC_data = tmp_data;
         end
         
         %----------------------------------------------------------------%
         function counts = getROIPhotonCounts(obj)
             % apply logical filter
             counts = obj.backpic;
             counts( ~obj.ROI_log ) = 0;
             
             % remore 0-padding
             counts( ~any(counts,2), : ) = [];  %rows
             counts( :, ~any(counts,1) ) = [];  %columns
             
             counts = counts - obj.backgroundCounts;
         end

         %----------------------------------------------------------------%
         function linescan = get_linescan_x(obj) 
             %projects everything onto the x-axis
             linescan = sum(obj.ROI_data,1);
         end
         
         %----------------------------------------------------------------%
         function linescan = get_linescan_y(obj) 
            linescan = sum(obj.ROI_data,2); 
         end
         
         %----------------------------------------------------------------%
         function density = get_density_x(obj) 
            density = obj.pixelsize/obj.magnification/obj.abs_cross * 1e2 * sum(obj.ROI_data,1); % density in 1/m
         end
         
         %----------------------------------------------------------------%
         function density = get_density_y(obj) 
            density = obj.pixelsize/obj.magnification/obj.abs_cross * 1e2 * sum(obj.ROI_data,2); % density in 1/m
         end
         
         %----------------------------------------------------------------%
         
         function density = get_density_xy(obj) 
            density = 1/obj.abs_cross * 1e4 * obj.ROI_data; % density in 1/m
         end
         
         %----------------------------------------------------------------%
         
         function shotnoise = get_shotnoise_x(obj)
             % apply logical filter
             tmp_data = obj.shotnoise_data;
             tmp_data( ~obj.ROI_log ) = 0;
             
             % remore 0-padding
             tmp_data( ~any(tmp_data,2), : ) = [];  %rows
             tmp_data( :, ~any(tmp_data,1) ) = [];  %columns
             
             shotnoise = (obj.pixelsize/obj.magnification)^4/obj.abs_cross^2 * sum(tmp_data,1);
         end
         
         %----------------------------------------------------------------%
         
         function shotnoise = get_shotnoise_y(obj)
             % apply logical filter
             tmp_data = obj.shotnoise_data;
             tmp_data( ~obj.ROI_log ) = 0;
             
             % remore 0-padding
             tmp_data( ~any(tmp_data,2), : ) = [];  %rows
             tmp_data( :, ~any(tmp_data,1) ) = [];  %columns
             
             shotnoise = (obj.pixelsize/obj.magnification)^4/obj.abs_cross^2 * sum(tmp_data,2);
         end
         
         function [sf, av, bv] = get_scalefactor(obj)
            sf = obj.scalefact;
            av = obj.atomval;
            bv = obj.backval;
         end
         
         
         %----------------------------------------------------------------%
         function atomnumber = get_atomnumber(obj,roionly)      
             % NOTE: roionly part of legacy code - no longer has a function
             
             % try to determine atomnumber, might fail if no image was
             % loaded, in this case: set atomnumber to zero
             try
                if isempty(obj.ROI_data)
                    calc_ROI_data(obj);
                end
                atomnumber = obj.pixelsize^2/obj.magnification^2/obj.abs_cross*(sum(sum(obj.ROI_data)));
           
             catch
                disp('Cant determine atomnumber!')
                atomnumber = 0; %set to zero
             end
         end
         
         %----------------------------------------------------------------%
         function data = get_fringe_data(obj,fringeSpacingGuess)
            
            %constants
            muB      = 9.2740154e-28;         % J/G  Bohr magneton
            hbar = 1.054e-34;                 % Js   
            amuKg    = 1.6605387e-27;         % g   mass of a proton
            mRb=87*amuKg;
            
            try
                %do the fit
                yy=obj.get_linescan_x();
                xx=1:length(yy); 
                xx=xx*obj.pixel2um;

                clear GCF
                clear CGCF
            
                [GCF CGCF smoothedCurve]=GaussCosineFitfun(xx,yy,1.7,obj.pixel2um,fringeSpacingGuess);       
                %results 
                data.xinput=xx;
                data.yinput=yy;
                
                data.GCF=GCF;
                data.CGCF=CGCF;

                data.TotalFit=GCF.FitValues;
                data.GaussCenter=GCF.GaussFit.Center;
                data.Contrast=GCF.CosineFit.Contrast;
                data.FringeSpacing=GCF.CosineFit.FringeSpacing;
                data.Phase=GCF.CosineFit.Phase;
                data.Displacement=2*pi*hbar*obj.TOF/data.FringeSpacing/mRb*1e9;   % calculate splitting in µm, TOF[ms]*1e-3, c2.b[µm]*1e-6, displacement[m]*1e6 --> *1e9
                data.ComDisplacement=2*pi*hbar*obj.TOF/data.CGCF.FringeSpacing/mRb*1e9;
                data.smoothedCurve=smoothedCurve;

                %errors
                dummyB=diff(confint(GCF.CosineFit.FitResult,0.682));
                data.ContrastError=dummyB(1)/2;
                
            catch
                disp('Can´t fit');
                
                yy=obj.get_linescan_x();
                xx=1:length(yy); 
                xx=xx*obj.pixel2um;
                data.xinput=xx;
                data.yinput=yy;
                
                data.TotalFit=0;
                data.GaussCenter=0;
                data.Contrast=0;
                data.FringeSpacing=0;
                data.Phase=0;
                data.Displacement=0;  
                data.ContrastError=0;
                data.CGCF.Center=0;
                data.CGCF.Sigma=0;
                data.CGCF.Contrast=0;
                data.CGCF.FringeSpacing=0;
                data.CGCF.Phase=0;
                data.CGCF.GOF.rsquare=0;
                data.CGCF.FitValues=0;
                data.CGCF.GOF.sse=0;
                data.ComGaussCenter=0;
                data.ComContrast=-1;
                data.ComFringeSpacing=0;
                data.ComPhase=0;
                data.ComDisplacement=0;  
                data.ComContrastError=0;
            end
         end
         
         %----------------------------------------------------------------%
         function plot_image(obj)
             image(obj.gain.*256.*obj.ROI_data);
         end
                
         %----------------------------------------------------------------%
         function plot_image_smoothed(obj,spansize)
             if nargin>1
                span = spansize; % Size of the averaging window
             else
                span = 4;
             end
             window = ones(span,span)/span^2; 
             abs_pic2 = convn(obj.ROI_data,window,'same');            
             image(obj.gain.*256.*abs_pic2);
         end
         
         
         %----------------------------------------------------------------%
         function plot_image_3d(obj)
             surf(obj.gain.*256.*obj.ROI_data);
             shading interp;
             axis off;
         end
         
         %----------------------------------------------------------------%
         function plot_image_3d_smooth(obj,spansize)
             if nargin>1
                span = spansize; % Size of the averaging window
             else
                span = 4;
             end
             window = ones(span,span)/span^2; 
             abs_pic2 = convn(obj.ROI_data,window,'same');
             surf(obj.gain.*256.*abs_pic2);
             shading interp;
             axis off;
         end
         
         
         %----------------------------------------------------------------%
         function plot_whole_image(obj)
             image(obj.gain.*256.*obj.abs_data);
             hold on
             visboundaries(obj.RBC_log, 'Color', [0 0 0 ], 'EnhanceVisibility', false)
             visboundaries(obj.ROI_log, 'Color', [1 1 1 ], 'EnhanceVisibility', false, 'LineWidth', 1.5)
         end
         
         %----------------------------------------------------------------%
         function load_ROI_data(obj,newdata)
             obj.ROI_data = newdata;
             %do not change anything in the particular data (with the gain!).
             %obj.gain = 1;
         end
         
         %----------------------------------------------------------------%
         function load_rescaled_data(obj,newdata)
             obj.abs_data = newdata;
             obj.calc_ROI_data;
             %do not change anything in the particular data (with the gain!).
             %obj.gain = 1;
         end
         
         %----------------------------------------------------------------%
         function data = get_ROI_data(obj)
             data = obj.ROI_data;
         end
         
         %------------------------------------experimental------------------%
         function bimodal_fit(obj) 
             if strcmp(obj.imaging_system , 'longitudinalAndor') ...
                 || strcmp(obj.imaging_system , 'Pixelfly')
                 x=obj.ROI(3):obj.ROI(4);
                 y=obj.get_linescan_x();
             elseif strcmp(obj.imaging_system , 'transversalAndor') 
                 x=obj.ROI(1):obj.ROI(2);    
                 y=obj.get_linescan_y();
             end
             
             [lower,upper,GaussFit,ParaFit,sseGauss,ssePara,GaussError,bogoBEC] = GaussParabolaFit(x,y);
         end
         
         %----------------------------------------------------------------%
         function [hor_pos, ver_pos, hor_width, ver_width,fit_row,fit_col,hor_fit_amp,ver_fit_amp] = gauss_fit(obj)
            
            [abs_pic_rows,abs_pic_cols]=size(obj.ROI_data);
            rowAverage=obj.get_linescan_x();
            colAverage=obj.get_linescan_y()';
            %fit gaussian to get width
            %% make the fit

            f2 = fittype('a1*exp( -( (x-b1)^2/(2*c1^2) ) )+d1' );

            %[d,g]=fit(x',log(meanatomnum),f2)
            [~, indcol] = max(colAverage);
            [~, indrow] = max(rowAverage);

            options_col = fitoptions(f2);
            options_col.Lower = [0 0 eps -1000];
            options_col.Upper = [1000 300 150 1000];
            %options_col.StartPoint=[50 70 100 0];
            options_col.StartPoint=[300 indcol 100 0];

            options_row = fitoptions(f2);
            options_row.Lower = [0 0 eps -1000];
            options_row.Upper = [1000 300 150 1000];
            options_row.StartPoint=[300 indrow 100 0];

            %avoid inf values in the fit
            cols=[1:abs_pic_cols];
            rows=[1:abs_pic_rows];
            fit_row=fit(cols(rowAverage~=Inf)',rowAverage(rowAverage~=Inf)',f2,options_row);
            fit_col=fit(rows(colAverage~=Inf)',colAverage(colAverage~=Inf)',f2,options_col);


            hor_pos=fit_row.b1*obj.pixel2um;
            ver_pos=fit_col.b1*obj.pixel2um;
            hor_width=fit_row.c1*obj.pixel2um;
            ver_width=fit_col.c1*obj.pixel2um;
            hor_fit_amp=fit_row.a1;
            ver_fit_amp=fit_col.a1;
         end   
    end
end