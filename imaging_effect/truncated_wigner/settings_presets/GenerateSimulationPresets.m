clear all

%% Presets for double-well
preset_name             = 'dw_LL_thermo_preset';

% physical parameters
settings.atomnumber     = 8000;
settings.no_gridpoints  = 401;
settings.boxlength      = 100e-6;
settings.T_si           = 50e-9; % temperature of condensate
J_Hz_const              = 0;
settings.Jfunc_si       = 0; %@(t) 2*pi*J_Hz_const; % coupling as function of time
settings.periodic_BC    = 0; % periodic boundary conditions for LL
settings.omegaT_si      = 2*pi*1.4e3; % transverse trapping frequency
settings.double_well    = 1; % flag for double-well
settings.psf_DMD        = 1.0e-6;

% settings for saving and plotting
settings.plot_flag   = 1;
settings.save_flag   = 0;
settings.save_folder = '';
settings.save_name   = '';

% save settings as preset file
save(preset_name, '-struct', 'settings')


%% Presets for dressed single-well
preset_name             = 'dsw_LL_thermo_preset';

% physical parameters
settings.atomnumber     = 8000;
settings.no_gridpoints  = 401;
settings.boxlength      = 100e-6;
settings.T_si           = 50e-9; % temperature of condensate
J_Hz_const              = 0;
settings.Jfunc_si       = 0; %@(t) 2*pi*J_Hz_const; % coupling as function of time
settings.periodic_BC    = 0; % periodic boundary conditions for LL
settings.omegaT_si      = 2*pi*1.37e3; % transverse trapping frequency
settings.double_well    = 0; % flag for double-well
settings.psf_DMD        = 1.0e-6;

% settings for saving and plotting
settings.plot_flag   = 1;
settings.save_flag   = 0;
settings.save_folder = '';
settings.save_name   = '';

% save settings as preset file
save(preset_name, '-struct', 'settings')