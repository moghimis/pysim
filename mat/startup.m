% ##############  Saeed Moghimi
my_root='/home/server/pi/homes/moghimi/work/opt/matlab_tools/';
path(path, fullfile(my_root, 'mexcdf/mexnc', ''))
path(path, fullfile(my_root, 'mexcdf/snctools', ''))
path(path, fullfile(my_root, 'geo', ''))
path(path, fullfile(my_root, 'm_map', ''))
path(path, fullfile(my_root, 'seawater', ''))
path(path, fullfile(my_root, 'tidal_ellipse', ''))
path(path, fullfile(my_root, 't_tide', ''))
path(path, fullfile(my_root, 'LP_Bathymetry/Mfiles', ''))
path(path, fullfile(my_root, 'roms_swan_bathy_tools', ''))
path(path, fullfile(my_root, 'utils', ''))
path(path, fullfile(my_root, 'matlab_netCDF_OPeNDAP', ''))
path(path, fullfile(my_root, 'ncx/ncx', ''))


% OSU stuff ###################################################

cil1= '/home/ruby/matlab/CIL';
path(path, fullfile(cil1, 'UTM', ''))
path(path, fullfile(cil1, 'makeRelease/r201105/argusDB', ''))
path(path, fullfile(cil1, 'makeRelease/r201105/argusDB.mysql2', ''))

opengl neverselect

%addpath(genpath('/home/server/student/homes/gwilson/matlab'));


%addpath('/home/ruby/matlab/CIL');
%CILpaths;


CILPATH='/home/ruby/matlab/CIL';
addpath( [ CILPATH filesep 'local' ] );
addpath( [ CILPATH filesep 'CILTools' ] );
addpath( [ CILPATH filesep 'OCR' ] );
addpath( [ CILPATH filesep 'argusDB' ] );
addpath( [ CILPATH filesep 'argusDB.mysql2' ] );
addpath( [ CILPATH filesep 'geometry6/tools' ] );
addpath( [ CILPATH filesep 'geometry6' ] );
addpath( [ CILPATH filesep 'runup' ] );
addpath( [ CILPATH filesep 'cBathy' ] );
addpath( [ CILPATH filesep 'SD97' ] );
addpath( [ CILPATH filesep 'merge2' ] );
addpath( [ CILPATH filesep 'merge' ] );
addpath( [ CILPATH filesep 'pixel' ] );
addpath( [ CILPATH filesep 'pixel/apps' ] );
addpath( [ CILPATH filesep 'LOESS' ] );
addpath( [ CILPATH filesep 'Duck' ] );
addpath( [ CILPATH filesep 'TSA' ] );
addpath( [ CILPATH filesep 'ephemeris' ] );
addpath( [ CILPATH filesep 'CelestialMechMat' ] );
addpath( [ CILPATH filesep 'curveFit' ] );
addpath( [ CILPATH filesep 'dots' ] );
addpath( [ CILPATH filesep 'diwasp' ] );
addpath( [ CILPATH filesep 'mpiv' ] );
addpath( [ CILPATH filesep 'UTM' ] );
addpath( [ CILPATH filesep 'keyword' ] );
addpath( [ CILPATH filesep 'autoGeom' ] );
addpath( [ CILPATH filesep 'reseed' ] );
addpath( [ CILPATH filesep 'homogeneous' ] );
addpath( [ CILPATH filesep 'newCamCal' ] );
addpath( [ CILPATH filesep 'makeTemplates' ] );
addpath( [ CILPATH filesep 'surveyGUI' ] );
addpath( [ CILPATH filesep 'stereo' ] );
addpath( [ CILPATH filesep 'polarization' ] );
%addpath( [ CILPATH filesep 'netcdf/mexnc' ] );
%addpath( [ CILPATH filesep 'netcdf/snctools' ] );
addpath( [ CILPATH filesep 'othercolor' ] );
addpath( [ CILPATH filesep 'GoogleMap' ] );
%addpath( [ CILPATH filesep 'desijnTool/argusDesignTool' ] );
clear CILPATH
DBConnect blossom;


