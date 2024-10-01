%-----------------------------------------------------------------------------%
%				                FUNDACION FULGOR
%-----------------------------------------------------------------------------%

%-----------------------------------------------------------------------------%
%                                   ADD PATH 
%-----------------------------------------------------------------------------%

proj_dir = mfilename('fullpath');
proj_dir = proj_dir(1: end - length(mfilename));

addpath(genpath([proj_dir, 'modules/']))
addpath(genpath([proj_dir, 'modules_tests/']))
addpath(genpath([proj_dir, 'tests/']))
addpath(genpath([proj_dir, 'tools/']))